# 7. Stellar parameters with FERRE

`ferre.py` is a thin driver for [FERRE](https://github.com/callendeprieto/ferre) (Allende Prieto et al. 2006) on grids written by `turbogrid.write_ferre_grid()`. It:

1. resamples spectra onto the grid's pixels in the rest frame;
2. writes the FERRE inputs;
3. runs FERRE;
4. reads the results back.

## 7.1 Fitting

```python
import ferre
out = ferre.fit([wave1, wave2], [flux1, flux2], [err1, err2], rv=[-121.4, 33.0],
                grid='gridm31/grid.ferre', ids=['star1', 'star2'],
                workdir='ferre_run', ferre_exe='/path/to/ferre.x', nthreads=8)
r = out['star1']
r['params']      # {'TEFF': 4521., 'LOGG': 1.93, 'METALS': -0.41, 'C': ..., 'ALPHA': ..., 'N': ...}
r['errors']      # same keys, 1-sigma
r['chi2']        # reduced chi^2
r['cov']         # covariance matrix (grid-label order)
r['model']       # best-fitting model on the grid pixels (r['wave'])
r['edge']        # labels within STEPS/8 of a grid edge
r['failed']      # True if FERRE flagged the fit as failed
```

Spectra must be continuum-normalised in the same way as the grid (the grid is normalised by default). Pixels outside your data, or marked bad, get an error of 10¹⁰, so FERRE ignores them.

## 7.2 Controlling FERRE

Default namelist:

| keyword | value | meaning |
|---|---|---|
| `NDIM`, `NOV`, `INDV` | all grid dimensions | fit every parameter |
| `INTER` | `min(3, min(N_P) − 1)` | cubic Bézier where the grid allows |
| `ALGOR` | 3 | Nelder–Mead |
| `ERRBAR`, `COVPRINT` | 1, 1 | errors and covariance |
| `INIT`, `INDINI`, `NRUNS` | 1, 1…, 1 | start from the `.ipf` values (default: grid centre) |
| `F_FORMAT`, `F_ACCESS` | 0, 0 | ASCII grid loaded into memory |

Anything else goes through `extra=`, which overrides or adds keywords:

```python
ferre.fit(..., extra={'CONT': 1, 'NCONT': 0,       # FERRE's own continuum renormalisation
                      'INDINI': '2 2 2 2 2 1', 'NRUNS': 32})
```

Starting values can be given per star as `init=[{'TEFF': 4500., 'LOGG': 2.0}, ...]`, using grid-label names. They are nudged 0.001 inside the grid edges, because FERRE rejects starting points that sit exactly on an edge.

## 7.3 From FERRE to line-by-line fitting

```python
star, star_err = ferre.to_star(r, rv=v)
# star     = {'teff':..., 'logg':..., 'mh':..., 'am':..., 'cm':..., 'nm':..., 'vmicro':..., 'rv':...}
# star_err = {'teff':..., 'logg':..., 'mh':..., ...}
```

If the grid has no vmicro axis, `vmicro` comes from the Holtzman et al. (2015) relation, 2.478 − 0.325 log g, unless you pass `vmicro=`. `star` and `star_err` go straight into `LineFitter.fit()`.

## 7.4 The lower-level pieces

* `read_header(grid)` parses the `&SYNTH` header.
* `grid_wave(hdr)` gives the pixel wavelengths.
* `to_grid_pixels(wave, flux, err, rv, hdr)` resamples one spectrum.
* `write_inputs(workdir, ids, flux, err, hdr, init, nthreads, extra=...)` writes the input files.
* `read_output(workdir, ids, hdr)` reads the results.

## 7.5 Three FERRE behaviours to know about

* FERRE ends with `STOP 1`, a **non-zero exit status, even after a normal run**. turbofit judges success by the output file instead of the exit code.
* With `COVPRINT=1`, the current FERRE writes the **full n×n covariance**, not the upper triangle described in the manual. turbofit reads either.
* The output parameter file holds the id, n parameters, n errors, then three more columns (fraction, S/N, log₁₀ reduced χ²), then the covariance.
