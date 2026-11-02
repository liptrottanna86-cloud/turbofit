# 8. Stellar parameters in local mode (no grid)

`starfit.py` fits a spectrum ASPCAP-style, without a precomputed grid. It suits single stars and modest samples, with any instrument. The full spectrum (or the windows you choose) is fitted for:

* Teff, log g, [M/H], [α/M], [C/M], [N/M];
* **vmicro**;
* **vmac and/or vsini**;
* the **RV** and a smooth **continuum**.

The model is built from Turbospectrum syntheses made on demand around the star, on the APOGEE MARCS atmospheres with C and α nodes. No empirical calibrations are applied.

## 8.1 Running it

```python
import turbogrid, starfit
sf = starfit.StarFitter(
    linelists=['atomic.list', 'molecular.list', 'Hlinedata'],
    inst=dict(R=22500.),                       # any instrument from chapter 4
    nodes=dict(teff=(15, 3000., 250.), logg=(12, -0.5, 0.5), mh=(14, -2.5, 0.25),
               am=(7, -1.0, 0.25), cm=(11, -1.5, 0.25)),
    atmod_lookup=turbogrid.make_atmod_lookup('/data/marcs_apogee', flat=True),
    broadening='vmac',                          # or 'vsini', or 'both'
    cache_dir='/scratch/me/turbofit_cache', n_jobs=16)

res = sf.fit(wave, flux, err,                   # observed frame, flux NOT normalised
             guess=dict(teff=4600, logg=2.2, mh=-0.3),
             fixed=dict(logg=2.45),             # optional, e.g. from asteroseismology
             mask=good_pixels,                  # optional: tellurics, bad regions
             mc=50)                             # optional Monte Carlo errors
print(res.summary())
res.star, res.star_err                          # ready for linefit
starfit.plot_fit(res, 'star1_fit.pdf')
```

Starting values that you don't give default to the middle of the node range for Teff and log g, 0 for the abundances, 1.5 km/s for vmicro and 5 km/s for vmac. The RV comes from a cross-correlation if you don't give one.

## 8.2 How the model is built

* **Atmospheres:** the 2⁵ nodes (Teff, log g, [M/H], [α/M], [C/M]) of the cell of the atmosphere grid that contains the current estimate. Each node is synthesised with its own [M/H], [α/M] and [C/M], exactly as a grid point would be.
* **[N/M] and vmicro** only affect the synthesis. Each node is computed at a 3-point ladder in [N/M] (default step 0.5 dex) and in log₁₀ vmicro (default step 0.15 dex).
* **Interpolation is in spectra, never atmospheres:** multilinear across the cell, and quadratic (Lagrange) in [N/M] and log vmicro.
* **Broadening** (vmac and/or vsini) is applied analytically to the interpolated spectrum, so it needs no syntheses.
* **The instrument** (R, R(λ), a per-pixel LSF array, or the APOGEE LSF) is applied exactly at your observed pixels, at the current RV. There is no resampling of the data.
* **The continuum** is a polynomial per chunk (degree `cont_deg`, default 3; chunks as in chapter 5 or `cont_chunks=`). It multiplies the model and is solved linearly at every step.
* **Fitting:** `scipy.optimize.least_squares` within the current cell and ladders. Between iterations, the RV is refined from χ²(RV) of the full model, with the instrument applied exactly.
* **Moving the cell:** if the solution sits on a face of the cell, or away from the middle of a ladder, the cell moves one node and the ladders re-centre, and the fit is repeated (up to `max_iter`, default 6).

**Caching.** Every synthesis is cached on disk, keyed by node, [N/M], vmicro, wavelength range, step, line lists and ¹²C/¹³C. Later iterations, later stars and other instruments with the same wavelength coverage reuse them.

## 8.3 Cost

* One cell is 32 nodes × 3 × 3 ladder points = **288 syntheses**. Fewer are needed when some axes are fixed or have a single node.
* Each cell move adds the new face's nodes. Moving in [N/M] or vmicro only adds new ladder points.
* By default only the wavelength ranges covered by your data are synthesised (one segment per data chunk, with a margin for the RV), so gaps between chips, gratings or setups cost nothing. Pass `segments=` to restrict it further to the windows you want to fit.

The first star in a new part of parameter space takes longest; neighbours that share nodes, and re-runs with different options, reuse the cache.

## 8.4 Results

`res` is a `StarFitResult`:

| attribute | meaning |
|---|---|
| `params`, `errors` | best values and 1σ errors (covariance, scaled by √χ²_red when that exceeds 1); vmicro in km/s |
| `errors_mc` | Monte Carlo errors, if `mc` > 0 |
| `cov`, `corr`, `free` | covariance and correlation matrices, and the order of the free parameters |
| `rv`, `rv_err` | radial velocity and its error |
| `chi2r` | reduced χ² |
| `flags` | e.g. `teff_at_grid_edge`, `not_converged` |
| `wave`, `flux`, `err`, `model`, `continuum`, `norm_flux`, `norm_model` | the fitted pixels and the fit |
| `star`, `star_err` | dictionaries for `linefit` and `synthesis` |
| `n_synth` | syntheses run in this session |

## 8.5 Choices worth thinking about

* **`broadening`:**
  * `'vmac'` (default) suits most giants;
  * `'vsini'` suits rotating stars;
  * `'both'` only makes sense at high resolution and S/N. At R ≈ 20 000 the two are nearly degenerate below about 10 km/s, and the correlation matrix will show it.
* **Fixed parameters:** for example, log g from asteroseismology, or Teff from photometry. Fixed values need not lie on a node; the cell still brackets them.
* **What to fit:** for very cool or crowded spectra, pass `segments=` (rest-frame windows to synthesise) and/or `mask=` (pixels to use), so that regions with known line-list problems or tellurics are left out.
