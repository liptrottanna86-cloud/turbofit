"""
rv.py
=====

Radial velocity by cross-correlation against a synthetic template.

    import rv
    v, verr, (vgrid, ccf) = rv.ccf_rv(wave, flux, twave, tflux, err=err)

The template is typically a turbofit synthesis (or a FERRE grid spectrum) at
roughly the right parameters, already broadened to the instrument with
broaden.convolve(). Both spectra should be continuum-normalised. The
velocity is measured on a common log-lambda grid, the CCF peak is refined
with a parabola, and the uncertainty follows Zucker (2003, MNRAS 342, 1291).
"""

import numpy as np

C_KMS = 299792.458


def _loggrid(w0, w1, dv):
    dln = dv / C_KMS
    n = int(np.floor(np.log(w1 / w0) / dln)) + 1
    return w0 * np.exp(dln * np.arange(n))


def ccf_rv(wave, flux, twave, tflux, err=None, vmin=-600.0, vmax=600.0, dv=None,
           mask=None):
    """
    Parameters
    ----------
    wave, flux : observed spectrum (normalised); err optional (1-sigma).
    twave, tflux : template at rest (normalised).
    vmin, vmax : velocity search range (km/s).
    dv : velocity step of the working grid; default 1/3 of the observed pixel.
    mask : boolean array, True = use pixel (e.g. exclude tellurics).

    Returns
    -------
    rv, rv_err, (v, ccf)
    """
    wave = np.asarray(wave, float); flux = np.asarray(flux, float)
    good = np.isfinite(flux)
    if err is not None:
        good &= np.isfinite(err) & (np.asarray(err) > 0)
    if mask is not None:
        good &= np.asarray(mask, bool)
    if dv is None:
        dv = np.median(np.diff(wave[good]) / wave[good][1:]) * C_KMS / 3.0
    # overlap of data with the template shifted over the whole velocity range
    w0 = max(wave[good][0], twave[0] * (1 + vmax / C_KMS))
    w1 = min(wave[good][-1], twave[-1] * (1 + vmin / C_KMS))
    if w1 <= w0:
        raise ValueError('template does not cover the data over the requested velocity range')
    lw = _loggrid(w0, w1, dv)
    fo = np.interp(lw, wave[good], flux[good])
    wt = np.ones_like(lw)
    if err is not None:
        e = np.interp(lw, wave[good], np.asarray(err, float)[good])
        wt = 1.0 / np.maximum(e, 1e-6)
    a = (1.0 - fo)
    a = (a - np.average(a, weights=wt)) * wt
    lags = np.arange(int(np.floor(vmin / dv)), int(np.ceil(vmax / dv)) + 1)
    ccf = np.empty(len(lags))
    na = np.sqrt(np.sum(a ** 2))
    for k, L in enumerate(lags):
        # template shifted by v = L*dv: observed lambda = rest * exp(L*dv/c)
        ft = np.interp(lw * np.exp(-L * dv / C_KMS), twave, tflux)
        b = (1.0 - ft)
        b = (b - np.average(b, weights=wt)) * wt
        nb = np.sqrt(np.sum(b ** 2))
        ccf[k] = np.sum(a * b) / (na * nb) if nb > 0 else 0.0
    v = lags * dv
    i = int(np.argmax(ccf))
    if 0 < i < len(ccf) - 1:
        y0, y1, y2 = ccf[i - 1], ccf[i], ccf[i + 1]
        den = y0 - 2 * y1 + y2
        off = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        vbest = (lags[i] + off) * dv
        cpk = y1 - 0.25 * (y0 - y2) * off
        c2 = den / dv ** 2                        # second derivative (km/s)^-2
    else:
        vbest, cpk, c2 = v[i], ccf[i], np.nan
    # Zucker (2003) eq. 5, with N = number of independent (observed) pixels
    N = int(np.sum(good & (wave >= w0) & (wave <= w1)))
    if np.isfinite(c2) and c2 < 0 and 0 < cpk < 1:
        verr = np.sqrt(-1.0 / (N * c2 / cpk * cpk ** 2 / (1 - cpk ** 2)))
    else:
        verr = np.nan
    if i in (0, len(ccf) - 1):
        verr = np.nan     # peak at the edge of the search range: not trustworthy
    # Doppler: log grid shift -> relativistic-consistent velocity for small v
    vbest = C_KMS * (np.exp(vbest / C_KMS) ** 2 - 1) / (np.exp(vbest / C_KMS) ** 2 + 1)
    return vbest, verr, (v, ccf)


def chi2_rv(wave, flux, err, twave, tflux, v0, half=None, nstep=41, mask=None):
    """
    Refine a velocity by chi^2 against the template (already at instrument
    resolution), scanning v0 +- half km/s and fitting a parabola to chi^2.
    The 1-sigma error is where chi^2 rises by 1 (scaled by the reduced chi^2
    when that exceeds 1). Recommended after ccf_rv().

    Returns rv, rv_err, (v, chi2)
    """
    wave = np.asarray(wave, float); flux = np.asarray(flux, float); err = np.asarray(err, float)
    good = np.isfinite(flux) & np.isfinite(err) & (err > 0)
    if mask is not None:
        good &= np.asarray(mask, bool)
    if half is None:
        half = 2.0 * np.median(np.diff(wave) / wave[1:]) * C_KMS
    vv = v0 + np.linspace(-half, half, nstep)
    chi = np.empty(nstep)
    for k, v in enumerate(vv):
        b = v / C_KMS
        ft = np.interp(wave[good] / np.sqrt((1 + b) / (1 - b)), twave, tflux)
        chi[k] = np.sum(((flux[good] - ft) / err[good]) ** 2)
    i = int(np.argmin(chi))
    sel = slice(max(0, i - 5), min(nstep, i + 6))
    c2, c1, c0 = np.polyfit(vv[sel] - vv[i], chi[sel], 2)
    if c2 <= 0:
        return vv[i], np.nan, (vv, chi)
    vbest = vv[i] - c1 / (2 * c2)
    chimin = c0 - c1 ** 2 / (4 * c2)
    red = max(1.0, chimin / max(good.sum() - 1, 1))
    return vbest, np.sqrt(red / c2), (vv, chi)
