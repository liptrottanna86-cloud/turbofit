"""
turbogrid.py
============

Grid-building utilities for your Turbospectrum-based synthesis code
(`synthesize()` / `do_turbospec()`), mirroring the workflow Synple uses
for building spectral grids:

    1) build a table of grid points (Teff, logg, [M/H], ...) + the
       atmosphere file and linelists that go with each point
    2) set up jobs that can be run in parallel -- either locally
       (multiprocessing) or on a cluster (GNU parallel / Slurm array)
    3) run the jobs (each one just calls `synthesize()` once)
    4) collect the finished spectra into a single grid file

Design notes / assumptions
---------------------------
* This assumes your synthesis code lives in a module that exposes
  `synthesize(teff, logg, mh, am, cm, nm, vmicro, wrange, dw, atmod,
  linelists, ...)` returning `(flux, cont, wave)`, exactly like the
  `synthesize()` you already have. Adjust TURBOSPEC_MODULE below (or
  the import line) if your module has a different name.
* Unlike Synple's `collect_kurucz`/`collect_marcs`, this module does
  NOT try to guess your atmosphere-file naming convention -- you pass
  in an `atmod_lookup(teff, logg, mh)` callable that returns the path
  to the right atmosphere file. Wire that up to whatever you already
  use in `utils`/`atmos` for finding/interpolating atmospheres.
* Rather than one job-script per grid point (Synple's `polysyn`
  approach), grid points are rows in a single table (`grid_params.txt`)
  and there's one generic driver (`run_one`) that is told *which row*
  to do by an integer index. This scales much better to large grids
  under Slurm array jobs / GNU parallel, and makes restarts trivial
  (a finished point is just skipped).
* Per-point outputs (outputs/spec_NNNNNN.npz) are written under a
  dedicated `outputs_dir`, which you should point at your HPC's data
  or scratch storage rather than home/project space -- see
  write_grid_table(). They're also SHARDED into subdirectories of
  `shard_size` points each, instead of one flat directory, since a
  directory with order-of-a-million files is a common way to overwhelm
  a shared filesystem's metadata server (slow `ls`/opens, sometimes an
  outright per-directory file-count quota). `outputs_dir` and
  `shard_size` are recorded in grid_meta.pkl, so run_one() (called
  standalone from job scripts) always finds the right place regardless
  of the job's working directory.
"""

import os
import sys
import glob
import json
import pickle
import linecache
import traceback
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

import synthesis  # <-- rename to match your module (the one with synthesize())
import broaden
import atmos       # Holtzman's atmos.py -- filename()/kurucz2turbo()/marcs2turbo()
import utils

# ---------------------------------------------------------------------
# 1) BUILD THE GRID TABLE
# ---------------------------------------------------------------------

def _expand(vals):
    """
    Accept either an explicit list/array of values, or a FERRE-style
    (n, start, step) triplet like Synple's tteff=(5, 5000., 250.), and
    return a 1D numpy array of values.
    """
    if isinstance(vals, (list, np.ndarray)) and not (
        isinstance(vals, tuple) and len(vals) == 3
    ):
        return np.atleast_1d(np.array(vals, dtype=float))
    if isinstance(vals, tuple) and len(vals) == 3:
        n, start, step = vals
        return start + step * np.arange(int(n))
    # single scalar
    return np.atleast_1d(np.array([vals], dtype=float))


def make_atmod_lookup(griddir, model='marcs', vers='x3', cache_dir=None, trim=0,
                       flat=False):
    """
    Build an `atmod_lookup(teff, logg, mh, am, cm) -> path` callable on
    top of your existing atmos.py, using:
        atmos.filename()     -- resolves (teff,logg,[M/H],[C/M],[a/M])
                                 to the raw grid file path (MARCS- or
                                 Kurucz-style directory layout)
        atmos.marcs2turbo() / atmos.kurucz2turbo()
                              -- converts it into the format Turbospectrum
                                 expects, and transparently falls back to
                                 the "<file>.filled" symlink left behind by
                                 atmos.fill_holes() when the exact model
                                 wasn't computed.

    Parameters
    ----------
    griddir : str
        Root of the raw atmosphere grid.
    model : 'marcs' or 'kurucz'
        Atmosphere family (passed through to atmos.filename()).
    vers : str
        MARCS filename version tag (only used for MARCS, default 'x3').
    flat : bool
        atmos.filename() normally returns a path that includes a
        subdirectory (e.g. 'mod_z+0.75/p6500_..._.mod' for MARCS,
        'm...c...o.../am....mod' for Kurucz). Set flat=True if your
        atmosphere files are NOT in that subdirectory layout -- i.e.
        every file sits directly in griddir -- and only the basename
        will be used to look inside griddir.
    cache_dir : str, optional
        Where converted, Turbospectrum-ready files are written. Defaults
        to '<griddir>/_turbo_cache'. A file already present there is
        reused rather than reconverted (also cached in memory for the
        life of the process), so this is cheap to call many times.
    trim : int
        Layers to trim -- passed straight through to the conversion
        routines.

    Returns
    -------
    callable(teff, logg, mh, am=0.0, cm=0.0) -> str or None
        None (with a printed warning) if neither the exact file nor its
        hole-filled stand-in exists.
    """
    if cache_dir is None:
        cache_dir = os.path.join(griddir, '_turbo_cache')
    os.makedirs(cache_dir, exist_ok=True)
    _mem_cache = {}

    def lookup(teff, logg, mh, am=0.0, cm=0.0):
        key = (round(teff, 1), round(logg, 3), round(mh, 3), round(am, 3), round(cm, 3))
        if key in _mem_cache:
            return _mem_cache[key]

        # atmos.filename(teff, logg, z, c, a, ...): z=[M/H], c=[C/M], a=[alpha/M]
        relfile = atmos.filename(teff, logg, mh, cm, am, model=model, vers=vers)
        infile = os.path.join(griddir, os.path.basename(relfile) if flat else relfile)
        outfile = os.path.join(cache_dir, os.path.basename(relfile) + '.turbo')

        if not os.path.exists(outfile):
            try:
                if model.lower() == 'marcs':
                    ret = atmos.marcs2turbo(infile, outfile, trim=trim)
                    if ret == -1:
                        raise FileNotFoundError(infile)
                else:
                    atmos.kurucz2turbo(infile, outfile, trim=trim)
            except Exception as e:
                # atmos.py's own fallback-to-".filled" logic has a couple of
                # bugs (undefined names) on the "nothing at all found" path,
                # so this catches broadly rather than just FileNotFoundError.
                print(f'no atmosphere for teff={teff} logg={logg} mh={mh} '
                      f'am={am} cm={cm} ({infile}): {e}')
                _mem_cache[key] = None
                return None

        _mem_cache[key] = outfile
        return outfile

    return lookup


def _holtzman_vmicro(logg):
    """Holtzman et al. (2015, AJ 150, 148) vmicro(logg) relation, km/s."""
    return 2.478 - 0.325 * logg


def nlte_atmod_lookup(teff, logg, mh, am=0.0, cm=0.0):
    """atmod_lookup for NLTE grids: no file, see make_grid_params()."""
    return 'nlte'


def make_grid_params(teff_grid, logg_grid, mh_grid, am=0.0, cm=0.0, nm=0.0,
                      vmicro=2.0, atmod_lookup=None, vmicro_log=False):
    """
    Build the cartesian-product list of grid points.

    Parameters
    ----------
    teff_grid, logg_grid, mh_grid : array-like or (n, start, step)
        Grid axes. am/cm/nm can each be a single value, an array, or a
        triplet too -- if given as arrays they are also expanded as grid
        axes (multiplying the grid size), otherwise they're held fixed
        for every point. Note am/cm also select *which atmosphere* is
        used (structure depends on them), so varying them re-triggers
        atmod_lookup; nm never affects the atmosphere and only varies
        the synthesis.
    vmicro : float, array-like, (n,start,step), or 'holtzman'
        Either a fixed value / grid axis like the others, or the string
        'holtzman' to derive it from logg via the Holtzman et al. (2015)
        relation vmicro = 2.478 - 0.325*logg instead of treating it as a
        free grid dimension (one value per teff/logg point, not an axis
        that multiplies the grid size).
    vmicro_log : bool
        If True, the vmicro values/triplet are log10(vmicro) (an even axis in
        log vmicro, as ASPCAP's LOG10VDOP); e.g. vmicro=(4, -0.3, 0.2),
        vmicro_log=True gives 0.50, 0.79, 1.26, 2.00 km/s. Use the label
        'lgvmicro' in write_ferre_grid()/fill_grid_holes.py.
    atmod_lookup : callable(teff, logg, mh, am, cm) -> str
        Required. Returns the path to the matching (Turbospectrum-ready)
        atmosphere-model file, or None if there isn't one -- see
        make_atmod_lookup() for a ready-made implementation built on
        atmos.py. Points with no match are dropped with a warning.
        For an NLTE grid (write_grid_table(nlte=..., nlte_elems=...)) use
        nlte_atmod_lookup: the atmosphere is then interpolated from the
        standard MARCS grid at run time, so am/cm do not change it.

    Returns
    -------
    list of dict, one per grid point, with keys:
        teff, logg, mh, am, cm, nm, vmicro, atmod
    """
    if atmod_lookup is None:
        raise ValueError("atmod_lookup(teff, logg, mh, am, cm) -> atmosphere file path is required")

    holtzman = isinstance(vmicro, str) and vmicro.lower() == 'holtzman'

    teffs = _expand(teff_grid)
    loggs = _expand(logg_grid)
    mhs = _expand(mh_grid)

    def axis_or_fixed(x):
        if isinstance(x, (list, np.ndarray, tuple)):
            return _expand(x), True
        return np.array([x]), False

    ams, _ = axis_or_fixed(am)
    cms, _ = axis_or_fixed(cm)
    nms, _ = axis_or_fixed(nm)
    vmicros = None if holtzman else axis_or_fixed(vmicro)[0]
    if vmicro_log and vmicros is not None:
        vmicros = 10 ** vmicros

    points = []
    skipped = 0
    for teff in teffs:
        for logg in loggs:
            vs = [_holtzman_vmicro(logg)] if holtzman else vmicros
            for mh in mhs:
                for a in ams:
                    for c in cms:
                        atmod = atmod_lookup(float(teff), float(logg), float(mh), float(a), float(c))
                        # 'nlte': atmosphere comes from the NLTE interpolator at run time
                        if not atmod or (atmod != 'nlte' and not os.path.exists(atmod)):
                            skipped += 1
                            continue
                        for n in nms:
                            for v in vs:
                                points.append(dict(
                                    teff=float(teff), logg=float(logg), mh=float(mh),
                                    am=float(a), cm=float(c), nm=float(n),
                                    vmicro=float(v), atmod=atmod,
                                ))
    if skipped:
        print(f'Warning: {skipped} (teff,logg,mh,am,cm) combinations had no matching atmosphere and were skipped')
    return points


def write_grid_table(params, outdir, wrange=(15000.0, 17000.0), dw=0.1,
                      linelists=None, solarisotopes=False, dospherical=True,
                      atmos_type=None, outputs_dir=None, shard_size=1000,
                      R_out=None, dw_out=None, wout=None, elems=None,
                      nlte=None, nlte_elems=None, broaden_kw=None, c12c13=None,
                      broaden_axis=None):
    """
    Write the grid-point table and shared settings to disk so that
    `run_one()` can be called standalone (from a job script / array
    task) without needing `params` passed in-process.

    Creates:
        outdir/grid_params.txt   -- one row per grid point
        outdir/grid_meta.pkl      -- shared settings (wrange, dw, linelists,
                                      outputs_dir, shard_size, ...)
        outputs_dir/               -- per-point results land here, sharded
                                       into subdirectories (see below)

    Parameters
    ----------
    outdir : str
        Where the (small) control files grid_params.txt/grid_meta.pkl go.
        This can stay on a regular home/project filesystem.
    outputs_dir : str, optional
        Where the (potentially huge -- one file per grid point) spectra
        actually get written. Defaults to '<outdir>/outputs', but for
        anything beyond a small grid you should point this explicitly at
        your HPC's data/scratch storage instead, e.g.
        outputs_dir='/scratch/you/grid1_outputs'. Stored as an ABSOLUTE
        path in grid_meta.pkl, so run_one() (invoked standalone by job
        scripts, from whatever directory Slurm/GNU parallel happens to
        start in) always writes to the right place.
    nlte, nlte_elems : nlte.NLTEConfig and list of elements, optional
        Build an NLTE grid (standard MARCS atmospheres; use
        atmod_lookup=nlte_atmod_lookup in make_grid_params()).
    broaden_axis : ('vsini' or 'vmac', values), optional
        An extra grid axis in stellar broadening (km/s), at NO extra synthesis
        cost: every synthesised spectrum is broadened with each value (after
        which broaden_kw is applied). collect_grid() then expands the grid and
        adds the columns '<name>' and 'lg<name>' (log10); use e.g. the label
        'lgvsini' (LGVSINI, as ASPCAP) with log-spaced values such as
        10**np.arange(0.0, 1.81, 0.3).
    c12c13 : float, optional
        12C/13C number ratio used for every spectrum of the grid (default:
        the solarisotopes behaviour of synthesize()).
    broaden_kw : dict, optional
        Post-synthesis broadening/resampling with broaden.convolve(), e.g.
        dict(R=22500, dlogw_out=6e-6, wave_out=None) or an LSF:
        dict(lsf=lsf_array, lsf_offsets=np.arange(-7, 8), wave_out=apogee_wave).
        Replaces R_out/dw_out/wout (which still work as before).
    shard_size : int
        Output files are split into subdirectories of at most this many
        points each: outputs_dir/00000/spec_000000.npz,
        outputs_dir/00000/spec_000999.npz, outputs_dir/00001/spec_001000.npz,
        etc., rather than one flat directory. A single directory holding
        order-of-a-million files is a common way to bring a shared HPC
        filesystem's metadata server (Lustre/GPFS) to a crawl -- slow
        `ls`, slow opens, sometimes an outright per-directory file-count
        quota -- so this keeps each directory's file count bounded no
        matter how large the grid is. Default 1000 (so e.g. a ~1e6-point
        grid ends up in ~1000 subdirectories); lower it further if your
        filesystem's per-directory limit is tighter than that.
    """
    os.makedirs(outdir, exist_ok=True)

    if outputs_dir is None:
        outputs_dir = os.path.join(outdir, 'outputs')
    outputs_dir = os.path.abspath(outputs_dir)
    os.makedirs(outputs_dir, exist_ok=True)

    cols = ['teff', 'logg', 'mh', 'am', 'cm', 'nm', 'vmicro', 'atmod']
    with open(os.path.join(outdir, 'grid_params.txt'), 'w') as f:
        f.write('# idx ' + ' '.join(cols) + '\n')
        for i, p in enumerate(params):
            f.write('{:d} {teff:.2f} {logg:.3f} {mh:.3f} {am:.3f} {cm:.3f} {nm:.3f} '
                    '{vmicro:.6f} {atmod}\n'.format(i, **p))   # 6 dp: keeps log-vmicro axes even

    meta = dict(wrange=list(wrange), dw=dw, linelists=linelists,
                solarisotopes=solarisotopes, dospherical=dospherical,
                atmos_type=atmos_type, outputs_dir=outputs_dir,
                shard_size=int(shard_size),
                # optional: convolve each spectrum to R_out and resample onto
                # wout[0]..wout[1] in steps of dw_out right after synthesis
                R_out=R_out, dw_out=dw_out, wout=list(wout) if wout else None,
                # optional [X/M] overrides passed to synthesize(elems=...), e.g.
                # reset Ne..Ti to 0 so the 'am' axis changes only O
                elems=elems,
                nlte=nlte, nlte_elems=list(nlte_elems) if nlte_elems else None,
                broaden_kw=broaden_kw, c12c13=c12c13,
                broaden_axis=(None if broaden_axis is None else
                              (str(broaden_axis[0]), [float(v) for v in broaden_axis[1]])))
    with open(os.path.join(outdir, 'grid_meta.pkl'), 'wb') as f:
        pickle.dump(meta, f)

    n_shards = (len(params) + shard_size - 1) // max(shard_size, 1)
    print(f'{len(params)} grid points written to {outdir}/grid_params.txt')
    print(f'outputs will be written to {outputs_dir} ({n_shards} shard subdirectories)')
    return os.path.join(outdir, 'grid_params.txt')


def _read_grid_meta(outdir):
    """Just the (small) shared-settings pickle -- see _read_grid_table for
    the full per-point table plus this."""
    with open(os.path.join(outdir, 'grid_meta.pkl'), 'rb') as f:
        meta = pickle.load(f)
    # backward-compatible defaults for grid_meta.pkl files written before
    # sharded/relocatable outputs were added
    meta.setdefault('outputs_dir', os.path.abspath(os.path.join(outdir, 'outputs')))
    meta.setdefault('shard_size', 1000)
    return meta


def _read_grid_table(outdir):
    rows = []
    path = os.path.join(outdir, 'grid_params.txt')
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            idx = int(parts[0])
            teff, logg, mh, am, cm, nm, vmicro = (float(x) for x in parts[1:8])
            atmod = parts[8]
            rows.append(dict(idx=idx, logg=logg, teff=teff, mh=mh, am=am, cm=cm,
                              nm=nm, vmicro=vmicro, atmod=atmod))
    meta = _read_grid_meta(outdir)
    return rows, meta


def _read_grid_row(outdir, idx):
    """
    Read a single row (by idx) of grid_params.txt without parsing the
    whole table. Matters once the table has ~1e6 rows: each Slurm array
    task handles many idx values (see array_step in write_slurm_array),
    and re-parsing the full multi-million-line file for every single
    point would itself become the kind of needless HPC filesystem load
    this module is otherwise trying to avoid.
    """
    path = os.path.join(outdir, 'grid_params.txt')
    # line 1 is the '# idx ...' header, and rows were written in order,
    # so row i sits at physical line i+2 -- read it directly rather than
    # scanning from the top.
    line = linecache.getline(path, idx + 2)
    if line:
        parts = line.split()
        if parts and parts[0].isdigit() and int(parts[0]) == idx:
            teff, logg, mh, am, cm, nm, vmicro = (float(x) for x in parts[1:8])
            atmod = parts[8]
            return dict(idx=idx, teff=teff, logg=logg, mh=mh, am=am, cm=cm,
                        nm=nm, vmicro=vmicro, atmod=atmod)
    # table wasn't written in strict idx-per-line order for some reason --
    # fall back to a full scan rather than silently returning the wrong point
    rows, _ = _read_grid_table(outdir)
    row = next((r for r in rows if r['idx'] == idx), None)
    if row is None:
        raise IndexError(f'no grid point with idx={idx} in {path}')
    return row


# ---------------------------------------------------------------------
# 2) RUN ONE GRID POINT (this is what gets parallelised)
# ---------------------------------------------------------------------

def _shard_subdir(idx, shard_size):
    return f'{idx // shard_size:05d}'


def _spec_paths(outputs_dir, idx, shard_size):
    """
    Output/fail-file paths for grid point `idx`, sharded into
    subdirectories of `shard_size` points each -- see write_grid_table's
    `shard_size` docs for why.
    """
    shard_dir = os.path.join(outputs_dir, _shard_subdir(idx, shard_size))
    outfile = os.path.join(shard_dir, f'spec_{idx:06d}.npz')
    failfile = os.path.join(shard_dir, f'spec_{idx:06d}.failed')
    return shard_dir, outfile, failfile


def run_one(outdir, idx, overwrite=False, verbose=True):
    """
    Compute a single grid point (row `idx` of grid_params.txt) and
    save the result to <outputs_dir>/<shard>/spec_{idx:06d}.npz, where
    outputs_dir and the shard layout come from grid_meta.pkl (see
    write_grid_table).

    Safe to call many times / from many workers at once: a point that
    already has a finished output file is skipped unless
    overwrite=True, so this is restart-friendly for both GNU parallel
    and Slurm array jobs.

    verbose : bool
        If True (default), prints exactly which outputs_dir/outfile this
        call resolved to and confirms the shard directory exists, BEFORE
        calling synthesize(). synthesize() itself can run for minutes with
        no output of its own, so without this there is no way to tell,
        from the log alone, where (or whether) a point is actually being
        written -- this print is what lets you `grep` a job's log and
        immediately see the real resolved path instead of guessing
        between old/new outputs_dir values. Set False to silence it (e.g.
        very high-throughput local runs where stdout itself becomes a
        bottleneck).
    """
    meta = _read_grid_meta(outdir)
    outputs_dir, shard_size = meta['outputs_dir'], meta['shard_size']
    shard_dir, outfile, failfile = _spec_paths(outputs_dir, idx, shard_size)

    if verbose:
        print(f'[run_one idx={idx}] outputs_dir={outputs_dir} '
              f'shard_size={shard_size} -> outfile={outfile}', flush=True)

    if os.path.exists(outfile) and not overwrite:
        if verbose:
            print(f'[run_one idx={idx}] already exists, skipping', flush=True)
        return outfile

    row = _read_grid_row(outdir, idx)
    os.makedirs(shard_dir, exist_ok=True)
    if verbose:
        print(f'[run_one idx={idx}] shard dir ready: '
              f'{shard_dir} (exists={os.path.isdir(shard_dir)})', flush=True)

    try:
        flux, cont, wave = synthesis.synthesize(
            row['teff'], row['logg'], mh=row['mh'], am=row['am'], cm=row['cm'],
            nm=row['nm'], vmicro=row['vmicro'], atmod=row['atmod'],
            atmos_type=meta['atmos_type'], wrange=meta['wrange'], dw=meta['dw'],
            linelists=meta['linelists'], solarisotopes=meta['solarisotopes'],
            dospherical=meta['dospherical'], elems=meta.get('elems'),
            nlte=meta.get('nlte'), nlte_elems=meta.get('nlte_elems'),
            c12c13=meta.get('c12c13'),
        )
        if np.isscalar(flux) and flux == 0.:
            raise RuntimeError('synthesize() returned a failure sentinel (0.,0.,0.)')
        extra = {}
        if meta.get('broaden_axis'):
            # one output spectrum per broadening value (no extra synthesis)
            bname, bvals = meta['broaden_axis']
            if bname not in ('vsini', 'vmac'):
                raise ValueError("broaden_axis name must be 'vsini' or 'vmac'")
            fl, co, w_out = [], [], None
            for v in bvals:
                w_out, fc = broaden.convolve(wave, np.vstack([flux, cont]),
                                             **dict(meta.get('broaden_kw') or {}, **{bname: v}))
                fl.append(fc[0]); co.append(fc[1])
            flux, cont, wave = np.array(fl), np.array(co), w_out
            extra = dict(bname=bname, bvals=np.array(bvals))
        elif meta.get('broaden_kw'):
            wave, fc = broaden.convolve(wave, np.vstack([flux, cont]), **meta['broaden_kw'])
            flux, cont = fc[0], fc[1]
        elif meta.get('R_out'):
            # convolve on the fine grid, then trim/resample to the output grid
            # (padding in wrange beyond wout absorbs the convolution edges)
            w_out = np.arange(meta['wout'][0], meta['wout'][1] + 0.5 * meta['dw_out'],
                              meta['dw_out'])
            # working ln-grid must be at least as fine as the input sampling,
            # otherwise narrow lines are point-sampled (aliased), not resolved
            osamp = max(5.0, 1.2 * wave[-1] / (meta['R_out'] * np.min(np.diff(wave))))
            wln, fc = convolve_to_R(wave, np.vstack([flux, cont]), meta['R_out'],
                                    dw_out=None, oversample=osamp)
            flux = np.interp(w_out, wln, fc[0])
            cont = np.interp(w_out, wln, fc[1])
            wave = w_out
        np.savez(outfile, flux=flux, cont=cont, wave=wave,
                 teff=row['teff'], logg=row['logg'], mh=row['mh'], am=row['am'],
                 cm=row['cm'], nm=row['nm'], vmicro=row['vmicro'], **extra)
        if verbose:
            print(f'[run_one idx={idx}] wrote {outfile}', flush=True)
        if os.path.exists(failfile):
            os.remove(failfile)
        return outfile
    except Exception:
        with open(failfile, 'w') as f:
            f.write(f'idx={idx}\n')
            f.write(traceback.format_exc())
        raise


# ---------------------------------------------------------------------
# 3a) LOCAL PARALLEL RUN (single multi-core machine)
# ---------------------------------------------------------------------

def build_grid_local(outdir, nproc=None, overwrite=False):
    """
    Run every grid point using a local process pool. Good for a
    workstation / single node. For a cluster, use
    write_slurm_array() or write_parallel_joblist() instead.

    Returns (n_ok, n_failed).
    """
    rows, meta = _read_grid_table(outdir)
    outputs_dir, shard_size = meta['outputs_dir'], meta['shard_size']
    todo = [r['idx'] for r in rows]
    if not overwrite:
        todo = [i for i in todo
                if not os.path.exists(_spec_paths(outputs_dir, i, shard_size)[1])]

    n_ok, n_failed = 0, 0
    with ProcessPoolExecutor(max_workers=nproc) as ex:
        futures = {ex.submit(run_one, outdir, i, overwrite): i for i in todo}
        for fut in as_completed(futures):
            i = futures[fut]
            try:
                fut.result()
                n_ok += 1
            except Exception as e:
                n_failed += 1
                print(f'grid point {i} failed: {e}')
    print(f'done: {n_ok} ok, {n_failed} failed, {len(rows) - len(todo)} already present')
    return n_ok, n_failed


# ---------------------------------------------------------------------
# 3b) CLUSTER PARALLEL RUN (GNU parallel or Slurm array)
# ---------------------------------------------------------------------

def write_parallel_joblist(outdir, joblist='joblist.txt', python_exe=sys.executable):
    """
    Write a job list you can feed to GNU parallel:

        parallel -j 32 :::: <outdir>/joblist.txt
    """
    rows, _ = _read_grid_table(outdir)
    this_file = os.path.abspath(__file__)
    path = os.path.join(outdir, joblist)
    with open(path, 'w') as f:
        for r in rows:
            f.write(f"{python_exe} {this_file} run_one {os.path.abspath(outdir)} {r['idx']}\n")
    print(f'wrote {len(rows)} commands to {path}')
    print(f'run with:  parallel -j <N> :::: {path}')
    return path


def write_slurm_array(outdir, jobname='tsgrid', time='04:00:00', mem='4G',
                       partition=None, cpus_per_task=1, array_step=None,
                       max_concurrent=None, max_array_size=None,
                       python_exe=sys.executable, script_name='submit_array.slurm',
                       parallel=1):
    """
    Write a Slurm array-job script. Submits one task per grid point
    (or per batch of array_step points, if given, to reduce the
    number of array tasks for very large grids).

    Parameters
    ----------
    array_step : int, optional
        Grid points handled per array task (a bash loop inside one task).
        Required for large grids -- see max_array_size below.
    max_concurrent : int, optional
        If given, caps how many array tasks run AT THE SAME TIME via
        Slurm's `--array=0-N%max_concurrent` syntax -- good shared-cluster
        etiquette so a big grid doesn't monopolize the partition, even
        though squeue already collapses the array into one line
        regardless of size.
    max_array_size : int, optional
        If given, checked against the number of array tasks implied by
        array_step and raises a clear error (rather than a confusing
        Slurm rejection later) if you'd exceed it -- e.g. pass 999 if
        your cluster's MaxArraySize is 999.

    Submit with:  sbatch <outdir>/submit_array.slurm
    """
    rows, _ = _read_grid_table(outdir)
    n = len(rows)
    this_file = os.path.abspath(__file__)
    outdir_abs = os.path.abspath(outdir)
    path = os.path.join(outdir, script_name)

    if parallel > 1:
        # each array task = one node: run `parallel` grid points at a time
        if array_step is None or array_step < parallel:
            raise ValueError('with parallel>1, array_step must be >= parallel')
        cpus_per_task = parallel
        n_batches = (n + array_step - 1) // array_step
        body = (
            f'export OMP_NUM_THREADS=1\n'
            f'START=$(( SLURM_ARRAY_TASK_ID * {array_step} ))\n'
            f'END=$(( START + {array_step} - 1 ))\n'
            f'if [ $END -ge {n} ]; then END=$(( {n} - 1 )); fi\n'
            f'seq $START $END | xargs -P {parallel} -I{{}} '
            f'{python_exe} {this_file} run_one {outdir_abs} {{}}\n'
            f'echo "task $SLURM_ARRAY_TASK_ID done: points $START-$END"\n'
        )
    elif array_step is None or array_step <= 1:
        n_batches = n
        body = f'{python_exe} {this_file} run_one {outdir_abs} $SLURM_ARRAY_TASK_ID\n'
    else:
        n_batches = (n + array_step - 1) // array_step
        body = (
            f'START=$(( SLURM_ARRAY_TASK_ID * {array_step} ))\n'
            f'END=$(( START + {array_step} - 1 ))\n'
            f'if [ $END -ge {n} ]; then END=$(( {n} - 1 )); fi\n'
            f'for idx in $(seq $START $END); do\n'
            f'    {python_exe} {this_file} run_one {outdir_abs} $idx\n'
            f'done\n'
        )

    if max_array_size is not None and n_batches > max_array_size:
        raise ValueError(
            f'{n_batches} array tasks (array_step={array_step or 1}) exceeds '
            f'max_array_size={max_array_size} -- raise array_step. '
            f'Minimum array_step to fit: {-(-n // max_array_size)}'
        )

    array_range = f'0-{n_batches - 1}'
    if max_concurrent:
        array_range += f'%{max_concurrent}'

    lines = [
        '#!/bin/bash',
        f'#SBATCH --job-name={jobname}',
        f'#SBATCH --array={array_range}',
        f'#SBATCH --time={time}',
        f'#SBATCH --mem={mem}',
        f'#SBATCH --cpus-per-task={cpus_per_task}',
        f'#SBATCH --output={outdir_abs}/slurm_%A_%a.out',
    ]
    if partition:
        lines.append(f'#SBATCH --partition={partition}')
    lines.append('')
    lines.append(body)

    with open(path, 'w') as f:
        f.write('\n'.join(lines))
    os.chmod(path, 0o755)
    print(f'wrote {path} ({n} grid points, {n_batches} array tasks, array={array_range})')
    print(f'submit with:  sbatch {path}')
    return path


# ---------------------------------------------------------------------
# 4) COLLECT RESULTS INTO A SINGLE GRID FILE
# ---------------------------------------------------------------------

def collect_grid(outdir, gridfile='grid.npz', require_common_wave=True):
    """
    Gather all finished <outputs_dir>/<shard>/spec_*.npz files into one
    grid file (written to outdir/gridfile) containing:
        params : structured array (teff, logg, mh, am, cm, nm, vmicro,
                 lgvmicro [, vsini/vmac, lgvsini/lgvmac if broaden_axis])
        wave   : 1D wavelength array (shared across all points)
        flux   : 2D array (npoints, npix)
        cont   : 2D array (npoints, npix)

    Also reports any grid points that are missing or marked .failed.
    """
    rows, meta = _read_grid_table(outdir)
    outputs_dir = meta['outputs_dir']
    files = sorted(glob.glob(os.path.join(outputs_dir, '**', 'spec_*.npz'), recursive=True))

    if not files:
        raise RuntimeError(f'no finished outputs found under {outputs_dir}')

    wave_ref = None
    params, fluxes, conts, used_idx = [], [], [], []
    for fpath in files:
        idx = int(os.path.basename(fpath).split('_')[1].split('.')[0])
        d = np.load(fpath)
        if wave_ref is None:
            wave_ref = d['wave']
        elif require_common_wave and not np.array_equal(wave_ref, d['wave']):
            raise ValueError(
                f'{fpath} has a different wavelength grid than the rest -- '
                'make sure all points used the same wrange/dw, or set '
                'require_common_wave=False and interpolate afterwards.'
            )
        base = (float(d['teff']), float(d['logg']), float(d['mh']), float(d['am']),
                float(d['cm']), float(d['nm']), float(d['vmicro']),
                float(np.log10(d['vmicro'])))
        if 'bvals' in d.files:
            bname = str(d['bname'])
            for k, v in enumerate(d['bvals']):
                params.append(base + (float(v), float(np.log10(v)) if v > 0 else -np.inf))
                fluxes.append(d['flux'][k]); conts.append(d['cont'][k])
        else:
            bname = None
            params.append(base)
            fluxes.append(d['flux']); conts.append(d['cont'])
        used_idx.append(idx)

    dtype = [('teff', 'f8'), ('logg', 'f8'), ('mh', 'f8'), ('am', 'f8'),
             ('cm', 'f8'), ('nm', 'f8'), ('vmicro', 'f8'), ('lgvmicro', 'f8')]
    if bname is not None:
        dtype += [(bname, 'f8'), ('lg' + bname, 'f8')]
    params_arr = np.array(params, dtype=dtype)
    flux_arr = np.vstack(fluxes)
    cont_arr = np.vstack(conts)

    outpath = os.path.join(outdir, gridfile)
    np.savez(outpath, params=params_arr, wave=wave_ref, flux=flux_arr, cont=cont_arr)

    all_idx = {r['idx'] for r in rows}
    missing = sorted(all_idx - set(used_idx))
    failed = sorted(glob.glob(os.path.join(outputs_dir, '**', '*.failed'), recursive=True))
    print(f'collected {len(used_idx)}/{len(rows)} grid points into {outpath}')
    if missing:
        print(f'{len(missing)} points not yet computed, e.g. idx={missing[:10]}')
    if failed:
        print(f'{len(failed)} points failed, see *.failed files, e.g. {failed[:5]}')
    return outpath


def read_grid(gridfile):
    """Load a grid file written by collect_grid()."""
    d = np.load(gridfile, allow_pickle=True)
    return d['params'], d['wave'], d['flux'], d['cont']


# ---------------------------------------------------------------------
# 4b) FILL HOLES IN THE COLLECTED GRID (spectra stage, not the atmod stage)
# ---------------------------------------------------------------------
#
# Lives in its own file, fill_grid_holes.py, deliberately with NO
# dependency on turbogrid.py/synthesis.py/atmos.py/utils.py -- it only
# needs numpy/scipy and the collected .npz, so it can run standalone on
# any node, and doesn't require the synthesis toolchain to be importable
# just to post-process spectra you already have. See that file's
# docstring for what it does and why (RBF interpolation of missing grid
# nodes from their neighbors, edge-vs-interior hole classification, etc).
# Imported here only so the tg.fill_grid_holes(...) pipeline call still
# works if fill_grid_holes.py sits next to turbogrid.py; run it directly
# (`python fill_grid_holes.py grid1/grid.npz ...`) if it doesn't.
try:
    from fill_grid_holes import fill_grid_holes  # noqa: F401
except ImportError:
    def fill_grid_holes(*args, **kwargs):
        raise ImportError(
            'fill_grid_holes.py not found next to turbogrid.py (or not on '
            'PYTHONPATH). It is a standalone script with no turbogrid '
            'dependency -- either place it alongside turbogrid.py, or run '
            'it directly: python fill_grid_holes.py <gridfile> --labels ...'
        )


# ---------------------------------------------------------------------
# 5) OPTIONAL: CONVOLVE TO A GIVEN RESOLVING POWER R
# ---------------------------------------------------------------------

def convolve_to_R(wave, flux, R, dw_out=None, oversample=5.0):
    """
    Convolve spectra synthesized at high resolution (e.g. dw=0.01) down to
    a constant resolving power R = lambda/FWHM, and put the result on a
    linear wavelength grid with step dw_out.

    Constant R means the line-spread function has constant FWHM in
    ln(wavelength) -- equivalently, a constant velocity width FWHM_v = c/R
    -- regardless of the input wavelength grid's spacing. So this:
      1) resamples onto a uniform ln(wave) grid fine enough to resolve R
         (oversample pixels per FWHM),
      2) convolves with a Gaussian of that constant width,
      3) resamples the result onto a final linear grid with step dw_out
         (or leaves it on the ln(wave) working grid if dw_out is None).

    Parameters
    ----------
    wave : 1D array
        Input wavelength grid, as produced by synthesize() (linear).
    flux : 1D or 2D array
        One spectrum, or (nspec, npix) stacked spectra, on `wave`.
    R : float
        Target resolving power (lambda / FWHM).
    dw_out : float, optional
        Output linear wavelength step, e.g. 0.1. If None, the convolved
        spectrum is left on its (finer, non-linear-in-wave) working grid
        -- normally you want to set this.
    oversample : float
        ln(wave) pixels per resolution element used for the convolution
        itself (not the final output sampling). 5 is usually plenty;
        raise it if you see ringing near sharp features.

    Returns
    -------
    wave_out : 1D array
    flux_out : 1D or 2D array (matches the dimensionality of `flux`)
    """
    from scipy.ndimage import gaussian_filter1d

    wave = np.asarray(wave, dtype=float)
    single = (np.ndim(flux) == 1)
    flux2d = np.atleast_2d(np.asarray(flux, dtype=float))

    fwhm_lnw = 1.0 / R  # constant FWHM in ln(wave) = FWHM in velocity/c
    d_lnw = fwhm_lnw / oversample
    lnw0, lnw1 = np.log(wave[0]), np.log(wave[-1])
    n_ln = int(np.ceil((lnw1 - lnw0) / d_lnw)) + 1
    lnw_grid = lnw0 + np.arange(n_ln) * d_lnw
    wave_ln = np.exp(lnw_grid)

    flux_ln = np.empty((flux2d.shape[0], n_ln))
    for i in range(flux2d.shape[0]):
        flux_ln[i] = np.interp(wave_ln, wave, flux2d[i])

    sigma_pix = (fwhm_lnw / (2.0 * np.sqrt(2.0 * np.log(2.0)))) / d_lnw
    flux_conv = gaussian_filter1d(flux_ln, sigma_pix, axis=1, mode='nearest')

    if dw_out is None:
        wave_out, flux_out = wave_ln, flux_conv
    else:
        n_out = int(round((wave[-1] - wave[0]) / dw_out)) + 1
        wave_out = wave[0] + np.arange(n_out) * dw_out
        flux_out = np.empty((flux2d.shape[0], n_out))
        for i in range(flux2d.shape[0]):
            flux_out[i] = np.interp(wave_out, wave_ln, flux_conv[i])

    return wave_out, (flux_out[0] if single else flux_out)


def convolve_grid_to_R(gridfile, outfile, R, dw_out=None, oversample=5.0):
    """
    Convenience wrapper: load a grid .npz from collect_grid(), convolve
    both flux and cont to resolving power R (dw_out sets the output
    sampling), and save a new .npz with the same layout (params unchanged,
    wave/flux/cont replaced) -- so it can be fed straight into
    write_ferre_grid() afterwards. Purely optional/one-off; the main
    build_grid_local -> collect_grid pipeline doesn't need this at all.
    """
    params, wave, flux, cont = read_grid(gridfile)
    wave_out, flux_out = convolve_to_R(wave, flux, R, dw_out=dw_out, oversample=oversample)
    _, cont_out = convolve_to_R(wave, cont, R, dw_out=dw_out, oversample=oversample)
    np.savez(outfile, params=params, wave=wave_out, flux=flux_out, cont=cont_out)
    print(f'convolved {gridfile} to R={R:g}'
          + (f', resampled to dw={dw_out}' if dw_out else '')
          + f' -> {outfile}  ({len(wave)} -> {len(wave_out)} pix)')
    return outfile


# ---------------------------------------------------------------------
# 5) WRITE A FERRE-FORMAT GRID
# ---------------------------------------------------------------------

def write_ferre_grid(gridfile, ferre_root,
                      labels=('vmicro', 'cm', 'nm', 'am', 'mh', 'logg', 'teff'),
                      normalize=True, resolution=None, vacuum=0, comment=None,
                      rtol=1e-4):
    """
    Write a FERRE-format ASCII "synth" grid file (header + flux data) from
    the .npz produced by collect_grid() -- teff/logg/vmicro/[C/M]/[alpha/M]
    matched to the atmosphere, [N/M] (and any broadening axis) varied only at
    synthesis.

    The header follows the FERRE SYNTH namelist (N_OF_DIM, N_P, STEPS,
    LLIMITS, LABEL, NPIX, WAVE, LOGW, ...), and the flux block is stored with
    the LAST label varying fastest, as FERRE expects.

    Parameters
    ----------
    gridfile : str
        Path to the .npz file from collect_grid().
    ferre_root : str
        Output file path. FERRE ascii synth files are a single file
        (header followed by the flux data), so this is written as one
        plain-text file at exactly this path.
    labels : tuple of str
        Which grid axes to expose to FERRE and in what order -- this
        becomes LABEL(1), LABEL(2), ... in the header. The LAST label varies
        fastest in the flux block -- verified with FERRE itself (pure
        interpolation at every node of a test grid returns that node). Any axis here that doesn't
        actually vary in the collected grid is dropped automatically
        (e.g. if you only ever used one vmicro).
    normalize : bool
        Divide flux by continuum before writing -- FERRE grids are almost
        always continuum-normalized. Default True.
    resolution : float, optional
        If given, written as the header's RESOLUTION keyword.
    vacuum : int
        Header VACUUM flag (0 = air wavelengths, matching Turbospectrum's
        output for wave > 2000 A). Default 0.
    comment : str, optional
        Free-text ID string for the header (defaults to the grid file name).
    rtol : float
        Relative tolerance used when checking that each axis is evenly
        spaced (required -- FERRE regular grids are defined by LLIMITS +
        STEPS*index) and that the requested grid is complete.

    Returns
    -------
    str : path to the written grid file (== ferre_root).
    """
    d = np.load(gridfile, allow_pickle=True)
    params, wave, flux, cont = d['params'], d['wave'], d['flux'], d['cont']
    npix = len(wave)

    dwave = np.diff(wave)
    dlog = np.diff(np.log10(wave))
    if np.allclose(dwave, dwave[0], rtol=rtol):
        logw, wave0, dw = 0, float(wave[0]), float(dwave[0])
    elif np.allclose(dlog, dlog[0], rtol=rtol):
        # log10-uniform (e.g. broaden_kw dlogw_out / APOGEE-style pixels)
        logw, wave0, dw = 1, float(np.log10(wave[0])), float(dlog[0])
    else:
        raise ValueError('wavelength grid is neither linear nor log10-linear -- '
                          'FERRE synth files need one of the two')

    # drop any requested axis that doesn't actually vary in this grid
    axes = []
    for lab in labels:
        vals = np.unique(np.round(params[lab], 6))
        if len(vals) > 1:
            axes.append((lab, vals))
    if not axes:
        raise ValueError('none of the requested labels vary in this grid')

    n_p, steps, llimits = [], [], []
    for lab, vals in axes:
        d_v = np.diff(vals)
        if not np.allclose(d_v, d_v[0], rtol=rtol):
            raise ValueError(f"axis '{lab}' is not evenly spaced ({vals}); "
                              'FERRE regular grids require uniform steps per axis')
        n_p.append(len(vals))
        steps.append(float(d_v[0]))
        llimits.append(float(vals[0]))

    ntotal = int(np.prod(n_p))
    if ntotal != len(params):
        raise ValueError(
            f'grid is not complete/regular for FERRE: expected {ntotal} points '
            f'({" x ".join(str(n) for n in n_p)}) but collected {len(params)}. '
            'FERRE needs every combination present -- either re-run the missing '
            'indices (see collect_grid()\'s report), drop an axis, or call '
            'fill_grid_holes() first to RBF-interpolate stand-ins for the gaps.'
        )

    # strides for "labels[0] fastest" storage order
    strides = [1] * len(axes)
    for i in range(len(axes) - 2, -1, -1):
        strides[i] = strides[i + 1] * n_p[i + 1]

    flux_out = np.empty((ntotal, npix), dtype=np.float32)
    filled = np.zeros(ntotal, dtype=bool)
    for row in range(len(params)):
        flat = 0
        for (lab, vals), stride in zip(axes, strides):
            idx = int(np.argmin(np.abs(vals - params[lab][row])))
            flat += idx * stride
        spec = flux[row] / cont[row] if normalize else flux[row]
        flux_out[flat] = spec
        filled[flat] = True
    if not filled.all():
        raise ValueError(f'{(~filled).sum()} grid slots ended up unfilled -- '
                          'this should not happen if the completeness check above passed')

    ferre_names = dict(mh='METALS', am='ALPHA', cm='C', nm='N', teff='TEFF', logg='LOGG',
                       vmicro='VMICRO', lgvmicro='LOG10VDOP', vsini='VSINI',
                       lgvsini='LGVSINI', vmac='VMACRO', lgvmac='LGVMACRO')
    label_names = [ferre_names.get(lab, lab.upper()) for lab, _ in axes]

    lines = ['&SYNTH']
    lines.append(f" ID = '{comment or os.path.basename(gridfile)}'")
    lines.append(f' N_OF_DIM = {len(axes)}')
    lines.append(' N_P = ' + ' '.join(str(n) for n in n_p))
    lines.append(' STEPS = ' + ' '.join(f'{s:.6f}' for s in steps))
    lines.append(' LLIMITS = ' + ' '.join(f'{l:.6f}' for l in llimits))
    for i, name in enumerate(label_names):
        lines.append(f" LABEL({i + 1}) = '{name}'")
    lines.append(f' NPIX = {npix}')
    if logw:
        lines.append(f' WAVE = {wave0:.10f} {dw:.10e}')
    else:
        lines.append(f' WAVE = {wave0:.6f} {dw:.6f}')
    lines.append(f' LOGW = {logw}')
    lines.append(f' VACUUM = {vacuum}')
    if resolution is not None:
        lines.append(f' RESOLUTION = {resolution}')
    lines.append(' /')

    with open(ferre_root, 'w') as f:
        f.write('\n'.join(lines) + '\n')
        for i in range(ntotal):
            row = flux_out[i]
            for j in range(0, npix, 8):
                f.write(' '.join(f'{v:14.6E}' for v in row[j:j + 8]) + '\n')

    print(f'wrote FERRE grid ({" x ".join(str(n) for n in n_p)} = {ntotal} points, '
          f'{npix} pix) to {ferre_root}')
    print(f'sanity-check with:  ascii2bin {ferre_root}   (from the FERRE distribution)')
    return ferre_root


# ---------------------------------------------------------------------
# CLI: lets job scripts call `python turbogrid.py run_one <outdir> <idx>`
# ---------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) >= 2 and sys.argv[1] == 'run_one':
        _outdir, _idx = sys.argv[2], int(sys.argv[3])
        run_one(_outdir, _idx)
    else:
        print(__doc__)