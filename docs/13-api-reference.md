# 13. API reference

Arguments are listed with their defaults. Full descriptions are in each function's docstring.

## synthesis

**`synthesize(teff, logg, mh=0.0, am=0.0, cm=0.0, nm=0.0, vmicro=2.0, elems=None, wrange=[15000, 17000], dw=0.1, atmod=None, atmos_type=None, dospherical=True, linelists=None, solarisotopes=None, workdir=None, save=False, verbose=False, nlte=None, nlte_elems=None, segments=None, c12c13=None)`**
Returns `(flux, cont, wave)`.

**`make_abundances(mh=0.0, am=0.0, cm=0.0, nm=0.0, elems=None)`**
Returns the A(X) array for Z = 1…99 (index Z−1).

**`do_turbospec(root, atmod, linelists, mh, am, abundances, wrange, dw, save=False, solarisotopes=False, babsma=None, atmos_type='marcs', spherical=True, vmicro=2.0, tfactor=1.0, verbose=False, nlte_info=None, segments=None, c12c13=None)`**
The low-level babsma + bsyn run in the current directory.

**`isotope_block(solarisotopes=False, c12c13=None)`**
The Turbospectrum `ISOTOPES` lines for carbon.

## abund

**`synthesize_sets(teff, logg, mh=0.0, am=0.0, cm=0.0, nm=0.0, vmicro=2.0, abund_sets=([],), elems=None, wrange=(15000, 17000), dw=0.1, atmod=None, atmos_type=None, dospherical=True, linelists=None, solarisotopes=None, segments=None, nlte=None, nlte_elems=None, rerun_babsma=True, workdir=None, save=False, c12c13=None)`**
Returns `(wave, flux[nset, npix], cont[nset, npix])`.

**`ladder(teff, logg, elem, offsets=(-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6), mh=0.0, am=0.0, cm=0.0, nm=0.0, elems=None, **kw)`**
Returns `(xm, wave, flux, cont)`.

**`isotope_ladder(teff, logg, ratios=(4, 6, 8, 10, 15, 20, 30, 50, 90), **kw)`**
Returns `(ratios, wave, flux, cont)`.

**`base_xm(elem, am=0.0, cm=0.0, nm=0.0, elems=None)`**
The starting [X/M] of an element.

## nlte

**`NLTEConfig(marcs_dir, nlte_dir, config_file=None, model_atom_dir=None, interpol_dir=None, comp='st')`**
Methods: `.available_elements()`, `.grid`, `.aux_length(el)`.

**`bracket_models(cfg, teff, logg, mh)`**
Returns `(list of 8 model paths, family)`.

**`prepare(cfg, teff, logg, mh, abundances, nlte_elems, workdir='.', verbose=False)`**
Returns `(atmod, nlte_info_file)`.

## broaden

**`convolve(wave, flux, R=None, lsf=None, lsf_offsets=None, lsf_halfwidth=None, wave_out=None, dw_out=None, dlogw_out=None, rv=0.0, vmac=0.0, vsini=0.0, eps=0.6, pixel_integrate=False)`**
Returns `(wave_out, flux_out)`.

**`lsf_matrix(wave_in, wave_out, R=None, lsf=None, lsf_offsets=None, lsf_halfwidth=None, rv=0.0, nsig=5.0, pixel_integrate=False, nsub=5)`**
Returns a sparse matrix of shape `(n_out, n_in)`.

**`stellar(wave, flux, vmac=0.0, vsini=0.0, eps=0.6)`**
Stellar broadening on the input grid.

**`doppler(wave, rv)`**, **`output_grid(wave, dw_out=None, dlogw_out=None, wrange=None)`**, **`rotation_kernel(dlnw, vsini, eps=0.6)`**

## apogee_lsf

**`ApogeeLSF(datadir=None, lsfid='05440020', waveid='02420038')`**
`datadir` defaults to `$APOGEE_LSF_DIR`.

* `.inst(wave_out, fiber='combo', offsets=None)` returns `dict(lsf=..., lsf_offsets=...)`.
* `.lsf_array(wave_out, fiber='combo', offsets=None)` returns `(lsf[n_out, n_k], offsets)`.
* `.wave2pix(wave, chip, fiber=300)`

**`apStarWavegrid()`**
The 8575-pixel apStar log-λ grid.

**`raw(x, xcenter, params)`**, **`unpack_lsf_params(lsfarr)`**
Low-level LSF evaluation.

## continuum

**`find_chunks(wave, flux=None, gap=None)`**
Returns a list of `(wmin, wmax)`.

**`normalise(wave, flux, err=None, chunks=None, deg=3, niter=10, lo=1.0, hi=3.0)`**
Returns `(cont, flux/cont, err/cont)`.

**`running_percentile(wave, flux, width, q=90.0)`**

## starfit

**`StarFitter(linelists, inst, nodes, atmod_lookup, atmos_type='marcs', dw=None, nm_step=0.5, lgvmicro_step=0.15, broadening='vmac', cont_deg=3, cont_chunks=None, c12c13=None, cache_dir=None, n_jobs=1, verbose=True)`**

* `.fit(wave, flux, err, guess=None, fixed=None, fit_rv=True, max_iter=6, mc=0, seed=None, wrange=None, segments=None, mask=None)` returns a `StarFitResult`.

**`StarFitResult`**

* Attributes: `params`, `errors`, `errors_mc`, `cov`, `corr`, `free`, `fixed`, `chi2r`, `rv`, `rv_err`, `flags`, `wave`, `flux`, `err`, `model`, `continuum`, `norm_flux`, `norm_model`, `n_synth`.
* Properties: `.star`, `.star_err`.
* Method: `.summary()`.

**`plot_fit(res, outfile=None, chunk_width=None, colors=('#3A1450', '#D6458F'))`**

## turbofit

**`analyse(wave, flux, err, linelists, inst, elements, mode='local', nodes=None, atmod_lookup=None, atmos_type='marcs', guess=None, fixed=None, broadening='vmac', grid=None, ferre_exe='ferre.x', fit_labels=None, ferre_extra=None, nlte=None, nlte_elems=(), mask=None, automask_kw=None, cont_deg=3, cont_chunks=None, c12c13=None, n_jobs=1, cache_dir=None, workdir='turbofit_run', mc=0, seed=None, starfit_kw=None, linefit_kw=None, fit_kw=None, verbose=True)`**
Returns a dict with keys `star, star_err, params, mask, lines, norm_flux, norm_err, mode`.

## rv

**`ccf_rv(wave, flux, twave, tflux, err=None, vmin=-600.0, vmax=600.0, dv=None, mask=None)`**
Returns `(rv, rv_err, (v, ccf))`.

**`chi2_rv(wave, flux, err, twave, tflux, v0, half=None, nstep=41, mask=None)`**
Returns `(rv, rv_err, (v, chi2))`.

## ferre

**`fit(waves, fluxes, errs, rv, grid, ids=None, init=None, workdir='ferre_run', ferre_exe='ferre.x', nthreads=1, extra=None, name='ferre', fit_labels=None, resample='cubic')`**
Returns `{id: result}`.

**`to_star(result, rv=0.0, vmicro=None)`**
Returns `(star, star_err)`.

**Lower-level functions:**

* `read_header(path)`
* `grid_wave(hdr)`
* `to_grid_pixels(wave, flux, err, rv, hdr, bad_err=1e10, kind='cubic')`
* `write_inputs(workdir, ids, flux, err, hdr, init=None, nthreads=1, name='ferre', extra=None, fit_labels=None)`
* `read_output(workdir, ids, hdr, name='ferre')`

## linefit

**`LineFitter(linelists, inst, nodes=None, atmod_lookup=None, atmos_type='marcs', nlte=None, nlte_elems=(), offsets=(-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6), dw=0.02, n_jobs=1, drv_max=5.0, vmac=3.0, vsini=0.0, cont_deg=1, max_extend=2, workdir=None, verbose=True)`**

* `.automask(elements, star, wave_obs, dx=0.3, min_response=0.02, max_blend=0.3, cut=0.25, exclude=None, wrange=None)` returns a list of window dictionaries.
* `.fit(wave, flux, err, star, elements=None, mask=None, star_err=None, vmac='global', iterate=1, sys_vmicro=False, sys_nlte=False, min_pix=4, profile=True, mc=0, seed=None)` returns a `LineFitResult`.
* `.mode(elem)` returns `'CNO-LTE'` or `'NLTE'`.

Use the pseudo-element `'C12C13'` (`linefit.ISO`) in `elements` or in masks to fit log₁₀(¹²C/¹³C), and `star['c12c13']` to set the ratio. `linefit.default_c12c13(star)` and `linefit.base_value(el, star)` give the starting values.

**`LineFitResult`**

* Attributes: `.lines`, `.summary`, `.vmac`, `.vsini`.
* Methods: `.table()` (astropy Table) and `.write(path)`.

**`read_mask(path, element=None)`**, **`write_mask(mask, path)`**

**`plot_lines(res, outfile=None, ncols=4, element=None, colors=('k', 'r'))`**

**`bracket(star, nodes)`**, **`corner_weights(br, point)`**
The node-cell geometry used by CNO-LTE mode.

## pipeline

**`analyse(wave, flux, err, grid, lf, elements, star_id='star', ferre_exe='ferre.x', workdir='turbofit_run', template_params=None, mask=None, nthreads=1, ferre_extra=None, vmicro=None, rv_range=(-600, 600), fit_labels=None, ferre_init=None, params_only=False, **fit_kw)`**
Returns a dict with keys `rv, rv_err, ferre, star, star_err, mask, lines`.

## turbogrid (new or changed arguments only)

**`make_grid_params(..., vmicro_log=False)`**

**`write_grid_table(..., nlte=None, nlte_elems=None, broaden_kw=None, c12c13=None, broaden_axis=None)`**

**`collect_grid(...)`**
Also stores `lgvmicro`, and expands `broaden_axis` grids with `vsini`/`lgvsini` (or `vmac`/`lgvmac`) columns.

**`nlte_atmod_lookup(teff, logg, mh, am=0.0, cm=0.0)`**
Returns `'nlte'`.

**`write_ferre_grid(...)`**
Writes `LOGW = 1` automatically for log₁₀-uniform wavelength grids.

For the rest of turbogrid (`make_atmod_lookup`, `make_grid_params`, `run_one`, `build_grid_local`, `write_slurm_array`, `write_parallel_joblist`, `collect_grid`, `convolve_grid_to_R`, `fill_grid_holes`), see `turbogrid_user_guide.pdf`.
