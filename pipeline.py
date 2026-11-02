"""
pipeline.py
===========

The full chain for one star, each step usable on its own:

    RV (CCF + chi^2, against a grid spectrum)
      -> FERRE on a turbogrid CNO grid (Teff, logg, [M/H], [a/M], [C/M], [N/M])
      -> automatic line masks for the requested elements
      -> line-by-line abundances (CNO-LTE, or NLTE where available) with
         per-line statistical and parameter-propagated systematic errors

    import pipeline, linefit
    lf = linefit.LineFitter(...)             # see linefit.py
    out = pipeline.analyse(wave, flux, err, grid='grid.ferre', lf=lf,
                           elements=['C', 'N', 'O', 'Mg', 'Al', 'Fe'],
                           ferre_exe='ferre.x', workdir='star1')
    out['star'], out['star_err'], out['lines'].summary
"""

import os
import numpy as np

import ferre
import rv as rvmod


def analyse(wave, flux, err, grid, lf, elements, star_id='star', ferre_exe='ferre.x',
            workdir='turbofit_run', template_params=None, mask=None, nthreads=1,
            ferre_extra=None, vmicro=None, rv_range=(-600, 600), fit_labels=None,
            ferre_init=None, params_only=False, **fit_kw):
    """
    Parameters
    ----------
    grid : FERRE grid file (from turbogrid.write_ferre_grid), or a LIST of
        grid files (e.g. separate cool and warm grids, as ASPCAP). With a list,
        the star is fitted on every grid and the best one is kept: the lowest
        chi^2 among the fits not pinned to a Teff/logg edge of their grid
        (or simply the lowest chi^2 if all are). The RV template is
        FERRE-interpolated from each grid.
    lf : linefit.LineFitter
    template_params : dict of grid labels for the RV template (default: grid centre).
    vmicro : microturbulence for the line fits (default: from the grid if it
        has a vmicro axis, else Holtzman et al. 2015 vmicro(logg)).
    fit_kw : passed to lf.fit (e.g. vmac='global', iterate=1, sys_vmicro=True).
    fit_labels / ferre_init : FERRE labels to fit and starting/fixed values
        (e.g. fit_labels=[...all but 'LOGG'], ferre_init={'LOGG': 2.45}).
    params_only : stop after FERRE (no masks / line fits; lf may be None).
    """
    os.makedirs(workdir, exist_ok=True)
    if isinstance(grid, (list, tuple)):
        return _analyse_multi(wave, flux, err, list(grid), lf, elements, star_id=star_id,
                              ferre_exe=ferre_exe, workdir=workdir,
                              template_params=template_params, mask=mask, nthreads=nthreads,
                              ferre_extra=ferre_extra, vmicro=vmicro, rv_range=rv_range,
                              fit_labels=fit_labels, ferre_init=ferre_init,
                              params_only=params_only, **fit_kw)
    hdr = ferre.read_header(grid)
    gw = ferre.grid_wave(hdr)

    # 1) RV template: FERRE pure interpolation at template_params
    tdir = os.path.join(workdir, 'template')
    pts = [template_params] if template_params else None
    dummy = np.ones((1, int(hdr['NPIX'])))
    ferre.write_inputs(tdir, ['template'], dummy, dummy, hdr, init=pts, extra={'NOV': 0})
    nml = open(os.path.join(tdir, 'input.nml')).read()
    nml = '\n'.join(l for l in nml.splitlines() if not l.strip().startswith('INDV'))
    open(os.path.join(tdir, 'input.nml'), 'w').write(nml + '\n')
    import subprocess
    subprocess.run([ferre_exe], cwd=tdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    tflux = np.loadtxt(os.path.join(tdir, 'ferre.mdl'), ndmin=2)[0]

    v, ve, _ = rvmod.ccf_rv(wave, flux, gw, tflux, err=err, vmin=rv_range[0], vmax=rv_range[1])
    v, ve, _ = rvmod.chi2_rv(wave, flux, err, gw, tflux, v)

    # 2) FERRE
    res = ferre.fit([wave], [flux], [err], rv=v, grid=grid, ids=[star_id],
                    workdir=os.path.join(workdir, 'ferre'), ferre_exe=ferre_exe,
                    nthreads=nthreads, extra=ferre_extra, fit_labels=fit_labels,
                    init=[ferre_init] if ferre_init else None)[star_id]
    star, star_err = ferre.to_star(res, rv=v, vmicro=vmicro)
    star_err['rv'] = ve
    if params_only:
        return dict(rv=v, rv_err=ve, ferre=res, star=star, star_err=star_err)

    # 3) masks + 4) line-by-line
    if mask is None:
        mask = lf.automask(elements, star, wave)
    lines = lf.fit(wave, flux, err, star, elements=elements, mask=mask,
                   star_err=star_err, **fit_kw)
    return dict(rv=v, rv_err=ve, ferre=res, star=star, star_err=star_err,
                mask=mask, lines=lines)


def _edge_labels(res):
    """Teff / logg edges hit by a FERRE result (other edges are allowed)."""
    return [l for l in (res.get('edge') or []) if l.upper() in ('TEFF', 'LOGG')]


def _analyse_multi(wave, flux, err, grids, lf, elements, star_id='star', workdir='turbofit_run',
                   params_only=False, mask=None, **kw):
    """Fit on every grid, keep the best, then (optionally) the line stage once."""
    fit_kw = {k: kw.pop(k) for k in list(kw) if k not in (
        'ferre_exe', 'template_params', 'nthreads', 'ferre_extra', 'vmicro', 'rv_range',
        'fit_labels', 'ferre_init')}
    tries = []
    for i, g in enumerate(grids):
        try:
            o = analyse(wave, flux, err, g, None, [], star_id=star_id,
                        workdir=os.path.join(workdir, f'grid{i}'), params_only=True, **kw)
        except Exception as e:                      # e.g. star outside this grid's range
            print(f'grid {g}: {e!r}')
            continue
        r = o['ferre']
        if r is None or r.get('failed'):
            continue
        o['grid'] = g
        o['edge_hit'] = _edge_labels(r)
        tries.append(o)
    if not tries:
        raise RuntimeError(f'FERRE failed on every grid for {star_id}')
    inside = [o for o in tries if not o['edge_hit']]
    best = min(inside or tries, key=lambda o: o['ferre']['chi2'])
    best['grids'] = {o['grid']: dict(chi2=o['ferre']['chi2'], edge=o['edge_hit']) for o in tries}
    if params_only:
        return best
    star, star_err = best['star'], best['star_err']
    if mask is None:
        mask = lf.automask(elements, star, wave)
    best['mask'] = mask
    best['lines'] = lf.fit(wave, flux, err, star, elements=elements, mask=mask,
                           star_err=star_err, **fit_kw)
    return best
