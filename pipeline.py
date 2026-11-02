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
            ferre_extra=None, vmicro=None, rv_range=(-600, 600), **fit_kw):
    """
    Parameters
    ----------
    grid : FERRE grid file (from turbogrid.write_ferre_grid); its .npz twin
        ('<grid without extension>.npz' or given via template_params) is not
        needed -- the RV template is FERRE-interpolated from the grid.
    lf : linefit.LineFitter
    template_params : dict of grid labels for the RV template (default: grid centre).
    vmicro : microturbulence for the line fits (default: from the grid if it
        has a vmicro axis, else Holtzman et al. 2015 vmicro(logg)).
    fit_kw : passed to lf.fit (e.g. vmac='global', iterate=1, sys_vmicro=True).
    """
    os.makedirs(workdir, exist_ok=True)
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
                    nthreads=nthreads, extra=ferre_extra)[star_id]
    star, star_err = ferre.to_star(res, rv=v, vmicro=vmicro)

    # 3) masks + 4) line-by-line
    if mask is None:
        mask = lf.automask(elements, star, wave)
    lines = lf.fit(wave, flux, err, star, elements=elements, mask=mask,
                   star_err=star_err, **fit_kw)
    return dict(rv=v, rv_err=ve, ferre=res, star=star, star_err=star_err,
                mask=mask, lines=lines)
