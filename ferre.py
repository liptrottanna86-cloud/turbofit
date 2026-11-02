"""
ferre.py
========

Thin driver for FERRE (Allende Prieto et al. 2006) on grids written by
turbogrid.write_ferre_grid(): resample observed spectra onto the grid's
pixels (rest frame), write the FERRE inputs, run it, read the results.

    import ferre
    hdr = ferre.read_header('grid1/grid.ferre')
    out = ferre.fit([wave], [flux], [err], rv=[-121.4], grid='grid1/grid.ferre',
                    ids=['star1'], workdir='ferre_run', nthreads=8,
                    ferre_exe='ferre.x')
    out['star1']['params']    # {'TEFF': ..., 'LOGG': ..., 'METALS': ..., 'C': ..., ...}
    out['star1']['errors'], out['star1']['chi2'], out['star1']['model']

File conventions follow the FERRE manual (and the APOGEE/ASPCAP driver):
.ipf (id + starting parameters in grid-label order), .frd/.err (one
spectrum per line), &LISTA namelist; the output parameter file holds
[id, params (n), errors (n), frac, S/N, log10 reduced chi2, covariance
(COVPRINT=1)]. Tested against FERRE built from
https://github.com/callendeprieto/ferre on a turbogrid-written grid.
Check the namelist against your own FERRE setup the first time you use it:
anything in `extra` overrides or adds keywords, e.g.
extra={'INTER': 3, 'ALGOR': 3, 'CONT': 1, 'NCONT': 0}.
"""

import os
import re
import subprocess
import numpy as np

import broaden


def read_header(path):
    """Parse the &SYNTH header of a FERRE ASCII grid."""
    hdr, txt = {}, []
    with open(path) as f:
        for line in f:
            txt.append(line)
            if line.strip() == '/':
                break
    body = ''.join(txt)
    for key, val in re.findall(r"^\s*([A-Z_]+(?:\(\d+\))?)\s*=\s*(.+?)\s*$", body, flags=re.M):
        val = val.strip().strip("'").strip()
        try:
            nums = [float(v) for v in val.split()]
            hdr[key] = nums if len(nums) > 1 else nums[0]
        except ValueError:
            hdr[key] = val
    n = int(hdr['N_OF_DIM'])
    hdr['labels'] = [hdr[f'LABEL({i + 1})'] for i in range(n)]
    hdr['N_P'] = [int(v) for v in np.atleast_1d(hdr['N_P'])]
    hdr['LLIMITS'] = list(np.atleast_1d(hdr['LLIMITS']))
    hdr['STEPS'] = list(np.atleast_1d(hdr['STEPS']))
    hdr['path'] = os.path.abspath(path)
    return hdr


def grid_wave(hdr):
    w0, dw = hdr['WAVE']
    npix = int(hdr['NPIX'])
    if int(hdr.get('LOGW', 0)) == 1:
        return 10 ** (w0 + dw * np.arange(npix))
    return w0 + dw * np.arange(npix)


def to_grid_pixels(wave, flux, err, rv, hdr, bad_err=1e10):
    """Shift an observed spectrum to the rest frame and resample it onto the
    grid's pixels; uncovered/bad pixels get err = bad_err."""
    gw = grid_wave(hdr)
    wr = np.asarray(wave, float) / broaden.doppler(1.0, rv)
    good = np.isfinite(flux) & np.isfinite(err) & (np.asarray(err) > 0)
    f = np.interp(gw, wr[good], np.asarray(flux)[good], left=np.nan, right=np.nan)
    e = np.interp(gw, wr[good], np.asarray(err)[good], left=np.nan, right=np.nan)
    # pixels whose nearest observed pixel was bad
    near = np.interp(gw, wr, good.astype(float), left=0, right=0) < 0.999
    miss = ~np.isfinite(f) | ~np.isfinite(e) | near
    f[miss] = 1.0
    e[miss] = bad_err
    return f, e


def write_inputs(workdir, ids, flux, err, hdr, init=None, nthreads=1, name='ferre',
                 extra=None):
    """Write <name>.ipf/.frd/.err/input.nml in workdir. init: list of dicts
    {label: value} (grid-label names) or None (grid centre)."""
    os.makedirs(workdir, exist_ok=True)
    labs, ll, st, npt = hdr['labels'], hdr['LLIMITS'], hdr['STEPS'], hdr['N_P']
    lo = np.array(ll); hi = lo + np.array(st) * (np.array(npt) - 1)
    with open(os.path.join(workdir, name + '.ipf'), 'w') as f:
        for k, i in enumerate(ids):
            p = [(init[k].get(l, 0.5 * (a + b)) if init else 0.5 * (a + b))
                 for l, a, b in zip(labs, lo, hi)]
            p = np.clip(p, lo + 0.001, hi - 0.001)
            f.write('{:40s}'.format(str(i)) + ''.join('{:12.3f}'.format(v) for v in p) + '\n')
    for ext, arr in (('.frd', flux), ('.err', err)):
        with open(os.path.join(workdir, name + ext), 'w') as f:
            for row in arr:
                f.write(' '.join('{:14.6e}'.format(v) for v in row) + '\n')
    n = len(labs)
    nml = dict(NDIM=n, NOV=n, INDV=' '.join(str(i + 1) for i in range(n)),
               SYNTHFILE_1_="'{}'".format(hdr['path']),
               PFILE=f"'{name}.ipf'", FFILE=f"'{name}.frd'", ERFILE=f"'{name}.err'",
               OPFILE=f"'{name}.spm'", OFFILE=f"'{name}.mdl'",
               ERRBAR=1, COVPRINT=1, INTER=int(min(3, min(npt) - 1)), ALGOR=3, INIT=1,
               INDINI=' '.join('1' for _ in range(n)), NRUNS=1,
               F_FORMAT=0, F_ACCESS=0, NTHREADS=int(nthreads))
    for k, v in (extra or {}).items():
        nml[k.replace('(', '_').replace(')', '_')] = v
    with open(os.path.join(workdir, 'input.nml'), 'w') as f:
        f.write('&LISTA\n')
        for k, v in nml.items():
            key = re.sub(r'_(\d+)_$', r'(\1)', k)
            f.write(f' {key} = {v}\n')
        f.write(' /\n')
    return os.path.join(workdir, 'input.nml')


def read_output(workdir, ids, hdr, name='ferre'):
    """Parse <name>.spm (+ .mdl/.frd if present) into {id: dict}."""
    labs = hdr['labels']
    n = len(labs)
    out = {}
    spm = os.path.join(workdir, name + '.spm')
    rows = {}
    with open(spm) as f:
        for line in f:
            t = line.split()
            if t:
                rows[t[0]] = np.array([float(x) for x in t[1:]])
    mdl = os.path.join(workdir, name + '.mdl')
    models = np.loadtxt(mdl, ndmin=2) if os.path.exists(mdl) else None
    lo = np.array(hdr['LLIMITS']); st = np.array(hdr['STEPS'])
    hi = lo + st * (np.array(hdr['N_P']) - 1)
    for k, i in enumerate(ids):
        v = rows.get(str(i))
        if v is None:
            out[i] = None
            continue
        p, e = v[:n], v[n:2 * n]
        d = dict(params=dict(zip(labs, p)), errors=dict(zip(labs, e)),
                 chi2=10 ** v[2 * n + 2] if len(v) > 2 * n + 2 else np.nan,
                 raw=v, wave=grid_wave(hdr))
        # COVPRINT=1: the manual describes the upper triangle, but the current
        # FERRE (github.com/callendeprieto/ferre) writes the full n x n matrix;
        # accept either
        rest = v[2 * n + 3:]
        if len(rest) >= n * n:
            d['cov'] = rest[:n * n].reshape(n, n)
        elif len(rest) >= n * (n + 1) // 2:
            cov = np.zeros((n, n))
            cov[np.triu_indices(n)] = rest[:n * (n + 1) // 2]
            d['cov'] = cov + np.triu(cov, 1).T
        d['edge'] = [l for l, a, b, s, x in zip(labs, lo, hi, st, p)
                     if (x - a) < s / 8 or (b - x) < s / 8]
        d['failed'] = bool(np.any(p < -999) or np.any(e < -0.01))
        if models is not None and k < len(models):
            d['model'] = models[k]
        out[i] = d
    return out


def fit(waves, fluxes, errs, rv, grid, ids=None, init=None, workdir='ferre_run',
        ferre_exe='ferre.x', nthreads=1, extra=None, name='ferre'):
    """Resample, write, run FERRE, read. Lists of spectra in, {id: result} out."""
    hdr = read_header(grid)
    ids = ids or [f'obj{i:05d}' for i in range(len(fluxes))]
    rv = np.broadcast_to(np.atleast_1d(rv), (len(fluxes),))
    F, E = [], []
    for w, f, e, v in zip(waves, fluxes, errs, rv):
        ff, ee = to_grid_pixels(w, f, e, v, hdr)
        F.append(ff); E.append(ee)
    write_inputs(workdir, ids, F, E, hdr, init=init, nthreads=nthreads, name=name, extra=extra)
    res = subprocess.run([ferre_exe], cwd=workdir, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    with open(os.path.join(workdir, name + '.log'), 'wb') as f:
        f.write(res.stdout)
    # NB: FERRE exits with 'STOP 1' (non-zero status) even after a normal run,
    # so success is judged by the output file, not the return code
    spm = os.path.join(workdir, name + '.spm')
    if not os.path.exists(spm) or os.path.getsize(spm) == 0:
        raise RuntimeError(f'FERRE failed, see {workdir}/{name}.log')
    return read_output(workdir, ids, hdr, name=name)


# FERRE label -> turbofit star-dict key
LABEL2KEY = {'TEFF': 'teff', 'LOGG': 'logg', 'METALS': 'mh', 'ALPHA': 'am', 'C': 'cm',
             'N': 'nm', 'VMICRO': 'vmicro', 'LOG10VDOP': 'vmicro'}


def to_star(result, rv=0.0, vmicro=None):
    """Convert one FERRE result into the star/star_err dicts linefit wants."""
    star, err = {}, {}
    for lab, v in result['params'].items():
        k = LABEL2KEY.get(lab.upper())
        if k is None:
            continue
        if lab.upper() == 'LOG10VDOP':
            v = 10 ** v
        star[k] = float(v)
        err[k] = float(result['errors'][lab])
    star['rv'] = rv
    if 'vmicro' not in star:
        star['vmicro'] = vmicro if vmicro is not None else 2.478 - 0.325 * star.get('logg', 2.5)
    return star, err
