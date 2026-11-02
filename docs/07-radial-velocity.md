# 7. Radial velocities

`rv.py` measures the velocity of a normalised spectrum against a synthetic template. The template should be at the instrument's resolution, for example a grid spectrum or `broaden.convolve()` output.

## 7.1 Cross-correlation, then χ²

```python
import rv
v0, e0, (vgrid, ccf) = rv.ccf_rv(wave, flux, twave, tflux, err=err,
                                 vmin=-600., vmax=600.)
v, e, (vv, chi2) = rv.chi2_rv(wave, flux, err, twave, tflux, v0)
```

**`ccf_rv`** works as follows:

* Both spectra are put on a common log-λ grid (step 1/3 of an observed pixel), the weighted, normalised cross-correlation is computed over `vmin..vmax`, and a parabola refines the peak.
* The uncertainty follows Zucker (2003), with N set to the number of *observed* pixels used.
* A `mask` (True = use the pixel) excludes tellurics or bad regions.
* A peak at the edge of the search range gives `rv_err = nan`. Widen the range if that happens.

**`chi2_rv`** scans around the CCF velocity, fits a parabola to χ²(v), and returns the velocity where χ² is lowest. Its error is where χ² rises by 1, scaled by the reduced χ² when that is above 1. It is the recommended final step.

## 7.2 How well calibrated are the errors?

On 30 simulated spectra (R = 22 500, S/N = 100, random velocities within ±300 km/s), the pulls (v − v_true)/σ had:

* a mean of 0.02 and a standard deviation of 0.88 for `ccf_rv`;
* a mean of 0.02 and a standard deviation of 0.76 for `chi2_rv`.

Both error estimates are therefore realistic to slightly conservative.

## 7.3 Per-line RV offsets

The line-by-line fitter also fits a small velocity offset per line, within ±`drv_max` km/s (default 5). This absorbs wavelength-calibration and line-list errors, and the offsets are reported as `drv` in the results table.
