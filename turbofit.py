"""
turbofit.py
===========

One call for one spectrum: stellar parameters -> line masks -> abundances.

    import turbofit, turbogrid

    out = turbofit.analyse(
        wave, flux, err,                       # observed frame; flux need not be normalised
        linelists=[...],
        inst=dict(R=22500.),                   # or any LSF (broaden.py / apogee_lsf.py)
        elements=['C', 'N', 'O', 'Mg', 'Al', 'Fe', 'C12C13'],
        mode='local',                          # or 'grid'
        nodes=dict(teff=(15, 3000., 250.), logg=(12, -0.5, 0.5), mh=(14, -2.5, 0.25),
                   am=(7, -1.0, 0.25), cm=(11, -1.5, 0.25)),
        atmod_lookup=turbogrid.make_atmod_lookup('/data/marcs_apogee', flat=True),
        guess=dict(teff=4600, logg=2.2, mh=-0.3),
        n_jobs=16, cache_dir='/scratch/me/turbofit_cache')

    out['star'], out['star_err']        # parameters (incl. vmicro, vmac/vsini, rv)
    out['lines'].summary                # abundances (and 12C/13C)

mode='local'  starfit.StarFitter: syntheses made on demand around this star,
              no grid needed (best for one star or modest samples).
mode='grid'   FERRE on a turbogrid grid (grid=...), then the same line-by-line
              stage (best for large samples that share one grid). grid can be
              a list of grids (e.g. cool and warm): each star is fitted on all
              of them and the best fit is kept automatically.

Nothing is calibrated empirically and there are no instrument presets: you
describe your instrument with inst=.
"""

import os
import numpy as np

import continuum
import linefit


def analyse(wave, flux, err, linelists, inst, elements, mode='local',
            nodes=None, atmod_lookup=None, atmos_type='marcs',
            guess=None, fixed=None, broadening='vmac',
            grid=None, ferre_exe='ferre.x', fit_labels=None, ferre_extra=None,
            nlte=None, nlte_elems=(), mask=None, automask_kw=None,
            cont_deg=3, cont_chunks=None, c12c13=None,
            n_jobs=1, cache_dir=None, workdir='turbofit_run',
            mc=0, seed=None, starfit_kw=None, linefit_kw=None, fit_kw=None,
            verbose=True):
    """
    Parameters (main ones)
    ----------------------
    wave, flux, err : one spectrum, observed frame (flux need not be normalised).
    linelists, inst : line lists and instrument (see broaden.py).
    elements : elements for the line-by-line stage (and 'C12C13' for the
        carbon isotope ratio); [] to stop after the stellar parameters.
    mode : 'local' (starfit) or 'grid' (FERRE on grid=..., one file or a list).
    nodes, atmod_lookup, atmos_type : the CNO atmosphere grid (needed for
        local mode and for CNO-LTE line fitting).
    guess, fixed : starting / fixed parameters (local mode; in grid mode,
        `fixed` holds grid labels via FERRE, e.g. fixed={'LOGG': 2.45}).
    broadening : 'vmac', 'vsini' (or 'both' in local mode).
    nlte, nlte_elems : NLTE for the abundance stage (see nlte.py).
    mask : your own line windows (else automask; automask_kw passed on).
    mc : Monte Carlo realisations for the parameter and abundance errors.
    starfit_kw / linefit_kw / fit_kw : extra options for StarFitter(),
        LineFitter() and LineFitter.fit().

    Returns
    -------
    dict(star, star_err, params (StarFitResult or FERRE result), mask, lines,
         norm_flux, norm_err, mode)
    """
    wave = np.asarray(wave, float); flux = np.asarray(flux, float); err = np.asarray(err, float)
    os.makedirs(workdir, exist_ok=True)

    if mode == 'local':
        import starfit
        if nodes is None or atmod_lookup is None:
            raise ValueError("mode='local' needs nodes= and atmod_lookup=")
        sf = starfit.StarFitter(linelists, inst, nodes, atmod_lookup, atmos_type=atmos_type,
                                broadening=broadening, cont_deg=cont_deg, cont_chunks=cont_chunks,
                                c12c13=c12c13, cache_dir=cache_dir, n_jobs=n_jobs,
                                verbose=verbose, **(starfit_kw or {}))
        res = sf.fit(wave, flux, err, guess=guess, fixed=fixed, mc=mc, seed=seed)
        star, star_err = res.star, res.star_err
        # normalise the whole spectrum with the continuum fitted alongside the model
        cont = np.full_like(flux, np.nan)
        for a, b in continuum.find_chunks(res.wave):
            m = (wave >= a) & (wave <= b)
            mm = (res.wave >= a) & (res.wave <= b)
            cont[m] = np.interp(wave[m], res.wave[mm], res.continuum[mm])
        nflux, nerr = flux / cont, err / cont
        params = res
        if verbose:
            print(res.summary())
    elif mode == 'grid':
        import ferre
        import rv as rvmod
        import pipeline
        if grid is None:
            raise ValueError("mode='grid' needs grid= (a FERRE grid from turbogrid)")
        _, nflux, nerr = continuum.normalise(wave, flux, err, chunks=cont_chunks, deg=cont_deg)
        out = pipeline.analyse(wave, nflux, nerr, grid=grid, lf=None, elements=[],
                               ferre_exe=ferre_exe, workdir=workdir, nthreads=n_jobs,
                               ferre_extra=ferre_extra, fit_labels=fit_labels,
                               ferre_init=fixed, params_only=True)
        params = out['ferre']
        if 'grids' in out:                   # several grids: record which one won
            params = dict(params, grid=out['grid'], grids=out['grids'])
        star, star_err = out['star'], out['star_err']
        if c12c13 is not None:
            star['c12c13'] = c12c13
    else:
        raise ValueError("mode must be 'local' or 'grid'")
    if c12c13 is not None:
        star['c12c13'] = c12c13

    lines = None
    if elements:
        lkw = dict(nodes=nodes, atmod_lookup=atmod_lookup, atmos_type=atmos_type,
                   nlte=nlte, nlte_elems=nlte_elems, n_jobs=n_jobs,
                   broadening='vsini' if broadening == 'vsini' else 'vmac',
                   workdir=os.path.join(workdir, 'linefit'), verbose=verbose)
        lkw.update(linefit_kw or {})
        lf = linefit.LineFitter(linelists, inst, **lkw)
        if mask is None:
            mask = lf.automask(elements, star, wave, **(automask_kw or {}))
        kind = 'vsini' if broadening == 'vsini' else 'vmac'
        # use the parameter-stage broadening if it was measured, else fit it globally
        fkw = dict(star_err=star_err, mc=mc, seed=seed, iterate=1,
                   vmac='star' if kind in star else 'global')
        fkw.update(fit_kw or {})
        lines = lf.fit(wave, nflux, nerr, star, elements=elements, mask=mask, **fkw)
    return dict(star=star, star_err=star_err, params=params, mask=mask, lines=lines,
                norm_flux=nflux, norm_err=nerr, mode=mode)
