"""
build_gc_grid.py
================

Example: build APOGEE-range grids for FERRE (grid mode). Parameters come
from FERRE on these grids; abundances afterwards from small on-the-fly
syntheses (apogee_gc_test.py / turbofit.analyse(mode='grid')).

Two grids, as in ASPCAP, because the MARCS Teff steps change at 4000 K:

    M   Teff 3000-4000 K in 100 K steps
    GK  Teff 3500-5500 K in 250 K steps

Atmosphere holes are left out at the synthesis stage; the missing SPECTRA
are filled afterwards with fill_grid_holes.py (RBF interpolation of
neighbouring spectra), so no atmospheres are ever interpolated.

The instrument is optional at synthesis time (INSTRUMENT below):

    None   unconvolved spectra at the synthesis step (default). Apply the
           instrument later with `convolve`, as many times as you like
           (different LSFs or R) without re-running Turbospectrum.
    'lsf'  the APOGEE LSF (GRID_LSF) applied straight away, onto the apStar
           pixels.
    R      a number: a Gaussian of that resolving power, onto the apStar pixels.

Steps (NAME = M or GK)
-----
    python build_gc_grid.py count NAME               # grid size, nothing written
    python build_gc_grid.py setup NAME               # grid table + Slurm array script
    python build_gc_grid.py time NAME 3              # synthesise 3 points here, timed
    sbatch <GRIDS>/NAME/submit_array.slurm           # the syntheses

  if INSTRUMENT is None, add the instrument afterwards (makes NAME_lsf or NAME_R22500):
    python build_gc_grid.py convolve NAME lsf        # or: convolve NAME 22500

  then, on NAME (instrument at synthesis) or NAME_lsf / NAME_R22500:
    python build_gc_grid.py collect NAME_lsf         # -> grid.npz (absolute fluxes)
    python build_gc_grid.py fill NAME_lsf            # holes -> grid_filled.npz (RBF)
    python build_gc_grid.py ferre NAME_lsf           # -> FERRE grid file, normalised
                                                     #    as set by NORMALISE

The grids are kept unnormalised; with NORMALISE = 'data' (default) every
model is normalised with exactly the same procedure and settings as the
observed spectra when the FERRE file is written.

Settings shared with apogee_gc_test.py (line lists, MARCS directory, LSF
files, air/vacuum) are imported from it.
"""

import os
import sys
import glob
import time
import shutil
import pickle
import numpy as np

import turbogrid
import broaden
import apogee_lsf
import apogee_gc_test as G

# =====================================================================
# CONFIG
# =====================================================================

GRIDS = G.BASE + '/gc_grids'                 # control files + collected grids
OUTPUTS = G.BASE + '/gc_grids_outputs'       # one file per grid point
PARTITION = 'nodes'
CPUS_PER_NODE = 64                           # points run in parallel per array task
POINTS_PER_TASK = 640                        # grid points per array task

COMMON = dict(
    mh=(11, -2.0, 0.25),                     # -2.00 .. +0.50
    am=(4, -0.25, 0.25),                     # -0.25 .. +0.50
    cm=(7, -1.0, 0.25),                      # -1.00 .. +0.50
    nm=(6, -0.5, 0.5),                       # -0.50 .. +2.00
    vmicro=(3, 0.0, 0.2),                    # log10: 1.00, 1.58, 2.51 km/s
)
GRID_DEFS = {
    'M':  dict(teff=[float(t) for t in range(3000, 4001, 100)], logg=(6, -0.5, 0.5), **COMMON),
    'GK': dict(teff=[float(t) for t in range(3500, 5501, 250)], logg=(8, 0.0, 0.5), **COMMON),
}
INSTRUMENT = None     # None (unconvolved, add later with `convolve`), 'lsf', or R (e.g. 22500)
VMAC = [4.0]          # km/s; one value = fixed, several (e.g. [2., 4., 8.]) = a FERRE axis
WRANGE = (15080.0, 17010.0)                  # synthesis range (apStar + margins)
DW = 0.05                                    # synthesis step (A)

# Normalisation of the FERRE grid. The collected grids keep ABSOLUTE fluxes
# (flux and continuum); the choice is made only when the FERRE file is written:
#   'data'  every model normalised with exactly the procedure used for the
#           observed spectra: continuum.normalise(wave, flux, chunks=CONT_CHUNKS,
#           deg=CONT_DEG) -- use the same chunks/deg for the data (grid mode:
#           turbofit.analyse(..., cont_chunks=CONT_CHUNKS, cont_deg=CONT_DEG))
#   'cont'  divided by the theoretical continuum
#   'none'  absolute fluxes (let FERRE normalise both with its CONT option)
NORMALISE = 'data'
CONT_CHUNKS = [(15140., 15810.), (15850., 16440.), (16470., 16960.)]   # APOGEE chips (vacuum)
CONT_DEG = 4

# LSF for 'lsf': LCO, average of the calibrations these stars' visits use
GRID_LSF = [('as', '22940020', '22670000'), ('as', '29380037', '27980000')]

# =====================================================================


def _axis(v):
    if isinstance(v, tuple) and len(v) == 3:
        n, s, d = v
        return s + d * np.arange(int(n))
    return np.atleast_1d(np.asarray(v, float))


def count(name):
    g = GRID_DEFS[name]
    n_atm = 1
    for k in ('teff', 'logg', 'mh', 'am', 'cm'):
        n_atm *= len(_axis(g[k]))
    n_syn = n_atm * len(_axis(g['nm'])) * len(_axis(g['vmicro']))
    npix_raw = int((WRANGE[1] - WRANGE[0]) / DW) + 1
    print(f'{name}: {n_atm} atmosphere nodes, {n_syn} syntheses, '
          f'{n_syn * len(VMAC)} spectra in the FERRE grid')
    print(f'   disk: ~{n_syn * npix_raw * 8 / 1e9:,.0f} GB unconvolved (float32), '
          f'~{n_syn * len(VMAC) * 8575 * 8 / 1e9:,.0f} GB on the apStar pixels')
    for k in ('teff', 'logg', 'mh', 'am', 'cm', 'nm', 'vmicro'):
        print(f'   {k:7s} {np.round(_axis(g[k]), 3)}')
    return n_syn


def grid_lsf():
    wv = apogee_lsf.apStarWavegrid()          # LSF on the detector (vacuum) grid
    lsf_sum, x = 0.0, None
    for prefix, lid, wid in GRID_LSF:
        L = apogee_lsf.ApogeeLSF(G.LSF_DIR, lsfid=lid, waveid=wid, prefix=prefix)
        lsf, x = L.lsf_array(wv, fiber='combo')
        lsf_sum = lsf_sum + lsf
    lsf = lsf_sum / len(GRID_LSF)
    return lsf / lsf.sum(axis=1, keepdims=True), x


def instrument_kw(instrument):
    """broaden.convolve() keywords for 'lsf' or a resolving power, onto apStar pixels."""
    if instrument == 'lsf':
        lsf, x = grid_lsf()
        return dict(lsf=lsf, lsf_offsets=x, wave_out=apogee_lsf.apStarWavegrid())
    return dict(R=float(instrument), wave_out=apogee_lsf.apStarWavegrid())


def setup(name):
    g = GRID_DEFS[name]
    outdir = os.path.join(GRIDS, name)
    lookup = turbogrid.make_atmod_lookup(G.MARCS_DIR, flat=True)
    params = turbogrid.make_grid_params(g['teff'], g['logg'], g['mh'], am=g['am'], cm=g['cm'],
                                        nm=g['nm'], vmicro=g['vmicro'], vmicro_log=True,
                                        atmod_lookup=lookup)
    if INSTRUMENT is None:
        kw = dict(float32=True)                           # unconvolved, vmac added later too
    else:
        kw = dict(broaden_kw=instrument_kw(INSTRUMENT), broaden_axis=('vmac', VMAC))
    turbogrid.write_grid_table(params, outdir, wrange=WRANGE, dw=DW, linelists=G.LINELISTS,
                               outputs_dir=os.path.join(OUTPUTS, name),
                               air_to_vac=G.LINELISTS_AIR,     # grid always in vacuum, like apStar
                               **kw)
    n_tasks = -(-len(params) // POINTS_PER_TASK)
    turbogrid.write_slurm_array(outdir, jobname=f'grid{name}', time='48:00:00', mem='64G',
                                partition=PARTITION, array_step=POINTS_PER_TASK,
                                parallel=CPUS_PER_NODE)
    print(f'{len(params)} points with an atmosphere ({n_tasks} array tasks); '
          f'the rest are holes, filled later from the spectra')


def time_points(name, n=3):
    outdir = os.path.join(GRIDS, name)
    rows, _ = turbogrid._read_grid_table(outdir)
    idx = np.linspace(0, len(rows) - 1, n).astype(int)
    t = []
    for i in idx:
        t0 = time.time()
        turbogrid.run_one(outdir, int(i), verbose=False)
        t.append(time.time() - t0)
        print(f'point {i}: {t[-1]:.0f} s')
    tm = np.median(t)
    cpu_h = tm * len(rows) / 3600
    print(f'median {tm:.0f} s per synthesis -> {cpu_h:,.0f} CPU-hours for {len(rows)} points '
          f'= {cpu_h / CPUS_PER_NODE:,.0f} node-hours at {CPUS_PER_NODE} cores per node')


# ---------------------------------------------------------------------
# instrument after synthesis
# ---------------------------------------------------------------------

_KW = None


def _conv_one(args):
    src, dst = args
    if os.path.exists(dst):
        return dst
    d = np.load(src)
    fl, co, w_out = [], [], None
    for v in VMAC:
        w_out, fc = broaden.convolve(d['wave'], np.vstack([d['flux'], d['cont']]),
                                     vmac=v, **_KW)
        fl.append(fc[0]); co.append(fc[1])
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    tmp = dst[:-4] + '.tmp.npz'
    np.savez(tmp, flux=np.array(fl, np.float32), cont=np.array(co, np.float32), wave=w_out,
             **{k: d[k] for k in ('teff', 'logg', 'mh', 'am', 'cm', 'nm', 'vmicro')},
             bname='vmac', bvals=np.array(VMAC, float))
    os.replace(tmp, dst)
    return dst


def convolve(name, instrument, nproc=None):
    """Apply vmac (VMAC) and the instrument to every finished point of the
    unconvolved grid NAME, onto the apStar pixels -> grid NAME_lsf / NAME_R<R>."""
    global _KW
    from multiprocessing import Pool
    tag = 'lsf' if instrument == 'lsf' else f'R{float(instrument):.0f}'
    src_dir, dst_name = os.path.join(GRIDS, name), f'{name}_{tag}'
    dst_dir = os.path.join(GRIDS, dst_name)
    os.makedirs(dst_dir, exist_ok=True)
    with open(os.path.join(src_dir, 'grid_meta.pkl'), 'rb') as f:
        meta = pickle.load(f)
    if meta.get('broaden_kw') or meta.get('broaden_axis'):
        sys.exit(f'{name} already has an instrument applied')
    src_out = meta['outputs_dir']
    dst_out = os.path.join(OUTPUTS, dst_name)
    meta = dict(meta, outputs_dir=dst_out, instrument=tag, broaden_axis=('vmac', list(VMAC)))
    with open(os.path.join(dst_dir, 'grid_meta.pkl'), 'wb') as f:
        pickle.dump(meta, f)
    shutil.copy(os.path.join(src_dir, 'grid_params.txt'), dst_dir)

    _KW = instrument_kw(instrument)
    files = sorted(glob.glob(os.path.join(src_out, '**', 'spec_*.npz'), recursive=True))
    jobs = [(f, os.path.join(dst_out, os.path.relpath(f, src_out))) for f in files]
    nproc = nproc or int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))
    print(f'{len(jobs)} spectra -> {dst_out} ({tag}, vmac = {VMAC}), {nproc} processes', flush=True)
    with Pool(nproc) as pool:
        for i, _ in enumerate(pool.imap_unordered(_conv_one, jobs, chunksize=16)):
            if (i + 1) % 10000 == 0:
                print(f'{i + 1}/{len(jobs)}', flush=True)
    print(f'done; next: collect {dst_name}')


# ---------------------------------------------------------------------
# collect / fill / FERRE
# ---------------------------------------------------------------------

def collect(name):
    turbogrid.collect_grid(os.path.join(GRIDS, name), gridfile='grid.npz')


def fill(name, nproc=16):
    import fill_grid_holes
    outdir = os.path.join(GRIDS, name)
    fill_grid_holes.fill_grid_holes(os.path.join(outdir, 'grid.npz'),
                                    outfile=os.path.join(outdir, 'grid_filled.npz'),
                                    labels=('teff', 'logg', 'mh', 'am', 'cm', 'nm', 'lgvmicro')
                                    + (('lgvmac',) if len(VMAC) > 1 else ()),
                                    nproc=nproc, dtype='float32')


def _norm_like_data(wave, flux):
    import continuum
    return continuum.normalise(wave, flux, chunks=CONT_CHUNKS, deg=CONT_DEG)[1]


def ferre(name):
    outdir = os.path.join(GRIDS, name)
    norm = dict(data=_norm_like_data, cont=True, none=False)[NORMALISE]
    labels = ['lgvmicro', 'cm', 'nm', 'am', 'mh', 'logg', 'teff']
    if len(VMAC) > 1:
        labels = ['lgvmac'] + labels
    turbogrid.write_ferre_grid(os.path.join(outdir, 'grid_filled.npz'),
                               os.path.join(outdir, f'f_gc_{name}.dat'), labels=tuple(labels),
                               normalize=norm, vacuum=1)  # apStar pixels, vacuum
    if NORMALISE == 'data':
        print(f'normalise the data the same way: cont_chunks={CONT_CHUNKS}, cont_deg={CONT_DEG}')


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    what, name = sys.argv[1], sys.argv[2]
    if what == 'count':
        count(name)
    elif what == 'setup':
        setup(name)
    elif what == 'time':
        time_points(name, int(sys.argv[3]) if len(sys.argv) > 3 else 3)
    elif what == 'convolve':
        convolve(name, sys.argv[3] if len(sys.argv) > 3 else 'lsf')
    elif what == 'collect':
        collect(name)
    elif what == 'fill':
        fill(name)
    elif what == 'ferre':
        ferre(name)
    else:
        sys.exit(__doc__)
