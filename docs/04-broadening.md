# 4. Broadening, LSFs and resampling

Everything between a raw Turbospectrum spectrum and something comparable with data goes through `broaden.py`. The same code is used by the grids, the RV code and the line fits, so they always agree.

## 4.1 The one call you need: `broaden.convolve()`

```python
import broaden
wave_out, flux_out = broaden.convolve(wave, flux,
        R=22500.,                 # instrument (see 4.2)
        rv=-45.3,                 # km/s, applied to the model
        vmac=4.0, vsini=0.0,      # stellar broadening (km/s)
        dlogw_out=6e-6)           # output grid (see 4.3)
```

`flux` can be one spectrum `(npix,)` or many `(nspec, npix)`. The same matrix is then applied to every row, so convolving a whole grid needs only one matrix build.

The steps, in order:

1. **Stellar broadening** (`broaden.stellar`):
    * Gaussian macroturbulence with profile exp(−(v/vmac)²), so σ = vmac/√2;
    * rotation with Gray's kernel (limb darkening `eps = 0.6`).
2. **Instrumental LSF**, the Doppler shift and resampling, all in **one sparse matrix** (`broaden.lsf_matrix`). The LSF is evaluated at each *output* pixel's own wavelength, which is what makes a varying R or a per-pixel LSF exact.

## 4.2 Describing your instrument

| your instrument | arguments |
|---|---|
| constant resolving power | `R=22500.` (Gaussian, FWHM = λ/R) |
| resolving power varying with wavelength | `R=(wave_R, R_values)` (interpolated) or `R=lambda lam: ...` |
| your own kernel shape | `lsf=lambda dlam, lam: weights, lsf_halfwidth=2.0` (Å) |
| a per-pixel LSF array (APOGEE, NIRSpec, …) | `lsf=array` of shape `(n_output_pixels, n_k)` and `lsf_offsets` of length `n_k`, in *output pixels* |
| detector pixel not negligible vs LSF | add `pixel_integrate=True` (averages the kernel over each output pixel) |

Kernels are normalised automatically.

> [!TIP]
> Sample LSF arrays finely, ≤ 0.2 pixel (e.g. `lsf_offsets=np.arange(-7, 7.01, 0.1)`), because the kernel is linearly interpolated between offsets. At 1-pixel sampling, a σ ≈ 1-pixel LSF is widened by about 8 %. A fitted vmac absorbs this, and abundances are unaffected because equivalent widths are conserved, but vmac is then biased low.

## 4.3 Choosing the output grid

Give exactly one of:

* `wave_out=array`: e.g. the observed pixels. This is what the fitting code uses.
* `dw_out=0.1`: a linear grid in Å.
* `dlogw_out=6e-6`: log₁₀-uniform, like APOGEE/ASPCAP pixels.

With none of them, the output stays on the input grid (shifted if `rv` is set). A useful rule of thumb is at least 2–3 output pixels per resolution element: λ/(R·Δλ) ≥ 2.

## 4.4 Doppler shifts

`broaden.doppler(wave, rv)` uses the relativistic formula λ_obs = λ·√((1+β)/(1−β)). With `rv` given, the model is shifted to the observed frame *before* the LSF is applied, so the LSF stays defined in the detector frame.

## 4.5 The APOGEE LSF

`apogee_lsf.py` evaluates the real APOGEE line-spread function at any output wavelength grid, for one fiber or the "combo" average. It has a Gauss–Hermite core plus wings, varies along each detector, and is integrated over the detector pixel. The result plugs straight into `broaden`, grids and line fits:

```python
import apogee_lsf, broaden
alsf = apogee_lsf.ApogeeLSF('/data/apogee_lsf')     # or set $APOGEE_LSF_DIR
wave = apogee_lsf.apStarWavegrid()                   # 8575 log-λ pixels, vacuum

inst = alsf.inst(wave, fiber='combo')                # dict(lsf=..., lsf_offsets=...)
w, f = broaden.convolve(wsyn, fsyn, wave_out=wave, vmac=4.0, **inst)
```

`ApogeeLSF` needs the six APOGEE calibration files, from the APOGEE data products or the `apogee` package:

* `apLSF-a/b/c-05440020.fits`;
* `apWave-a/b/c-02420038.fits`.

Other IDs can be passed as `lsfid=` and `waveid=`.

Details:

* `fiber` is `'combo'` (fibers 50, 100, …, 300 averaged), a single fiber 1–300, or a list of fibers to average.
* The LSF is sampled at 0.1-pixel offsets from −7 to +7 output pixels by default (`offsets=`).
* Rows that fall in chip gaps or off the detectors take the nearest valid row, so the convolution is defined everywhere. Mask those pixels in the data as usual.
* The profile already includes the detector pixel, so do **not** add `pixel_integrate=True`.
* turbofit's vmac is the Gaussian 1/e half-width (σ = vmac/√2). In Bovy's `apogee.spec.lsf`, `vmacro` is a FWHM; a FWHM of F there corresponds to vmac ≈ F/1.665 here.

The LSF evaluation is adapted from `apogee.spec.lsf` (J. Bovy, BSD-3), which is based on the APOGEE reduction code by J. Holtzman and D. Nidever. The adapted code was checked against analytic Gaussian LSF parameters: the second moment of the evaluated rows matches the expected σ² + 1/12 pixel² (pixel integration), converted to output pixels, to 0.1 %, on all three chips. With real APOGEE files it has not been checked against Bovy's original output.

## 4.6 How accurate is it?

The following were checked numerically:

* **Line width:** the measured σ of a convolved line is 0.29027 Å against 0.29028 Å expected.
* **Equivalent width** is conserved to a few parts in 10⁹.
* **R(λ):** a constant R(λ) matches a scalar R exactly.
* **LSF arrays:** a finely sampled Gaussian LSF array matches R to 8 × 10⁻⁵ in flux.
* **RV shift:** a 50 km/s shift is recovered to 0.004 km/s.
