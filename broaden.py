"""
broaden.py
==========

Everything between a raw Turbospectrum spectrum and something you can
compare to data, used identically by turbogrid (grids) and linefit (fits):

    1) stellar broadening      vmac (Gaussian macroturbulence) and/or vsini
    2) instrumental broadening one of
         R=22500                    constant resolving power (Gaussian LSF)
         R=callable / (wave, R)     resolving power that varies with wavelength
         lsf=callable               your own kernel: lsf(dlam, lam) -> weights
         lsf=array + lsf_offsets    per-output-pixel kernels sampled at pixel
                                    offsets (APOGEE / NIRSpec style LSF arrays)
    3) Doppler shift            rv (km/s)
    4) resampling               onto wave_out (e.g. observed pixels), a linear
                                step dw_out, or a log10 step dlogw_out

Steps 2-4 are done in ONE sparse matrix (convolve + resample), so the LSF is
evaluated at each OUTPUT pixel's own wavelength -- which is what makes a
wavelength-dependent R or a per-pixel LSF exact, rather than an
approximation.

    import broaden
    w, f = broaden.convolve(wave, flux, R=22500, dlogw_out=6e-6)
    w, f = broaden.convolve(wave, flux, R=(wave_R, R_of_wave), wave_out=wobs, rv=-45.)
    w, f = broaden.convolve(wave, flux, lsf=lsf_arr, lsf_offsets=np.arange(-7, 8), wave_out=wobs)

flux may be 1D (npix,) or 2D (nspec, npix): the same matrix is applied to
every spectrum, so convolving a whole grid costs one matrix build.
"""

import numpy as np
from scipy import sparse
from scipy.ndimage import gaussian_filter1d

C_KMS = 299792.458
FWHM2SIG = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


# ---------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------

def doppler(wave, rv):
    """Observed wavelength for rest wavelength `wave` and radial velocity rv
    (km/s, positive = receding), relativistic."""
    b = rv / C_KMS
    return wave * np.sqrt((1.0 + b) / (1.0 - b))


def output_grid(wave, dw_out=None, dlogw_out=None, wrange=None):
    """Linear (dw_out, A) or log10-linear (dlogw_out) grid spanning wrange
    (default: the input range)."""
    w0, w1 = (wave[0], wave[-1]) if wrange is None else wrange
    if dw_out is not None:
        n = int(np.floor((w1 - w0) / dw_out + 1e-9)) + 1
        return w0 + dw_out * np.arange(n)
    if dlogw_out is not None:
        l0, l1 = np.log10(w0), np.log10(w1)
        n = int(np.floor((l1 - l0) / dlogw_out + 1e-9)) + 1
        return 10 ** (l0 + dlogw_out * np.arange(n))
    raise ValueError('give wave_out, dw_out or dlogw_out')


def _to_lngrid(wave, flux2d, dlnw=None):
    """Resample onto a uniform ln(wave) grid at least as fine as the input."""
    if dlnw is None:
        dlnw = np.min(np.diff(np.log(wave)))
    lnw = np.arange(np.log(wave[0]), np.log(wave[-1]), dlnw)
    f = np.array([np.interp(lnw, np.log(wave), fl) for fl in flux2d])
    return lnw, f, dlnw


# ---------------------------------------------------------------------
# 1) stellar broadening
# ---------------------------------------------------------------------

def rotation_kernel(dlnw, vsini, eps=0.6):
    """Gray (2005) rotational broadening profile on a ln(wave) grid of step dlnw."""
    dv = dlnw * C_KMS
    n = int(np.ceil(vsini / dv))
    x = np.arange(-n, n + 1) * dv / vsini
    x = np.clip(x, -1, 1)
    k = (2 * (1 - eps) * np.sqrt(1 - x ** 2) + 0.5 * np.pi * eps * (1 - x ** 2))
    return k / k.sum()


def stellar(wave, flux, vmac=0.0, vsini=0.0, eps=0.6):
    """
    Macroturbulent (Gaussian, profile exp(-(v/vmac)^2), i.e. sigma =
    vmac/sqrt(2)) and rotational broadening, in velocity space. Returns flux
    on the input wavelength grid.
    """
    if vmac <= 0 and vsini <= 0:
        return flux
    single = np.ndim(flux) == 1
    f2 = np.atleast_2d(flux).astype(float)
    lnw, fl, dlnw = _to_lngrid(wave, f2)
    if vmac > 0:
        sig_pix = (vmac / np.sqrt(2.0)) / C_KMS / dlnw
        fl = gaussian_filter1d(fl, sig_pix, axis=1, mode='nearest')
    if vsini > 0:
        k = rotation_kernel(dlnw, vsini, eps)
        pad = len(k) // 2
        fp = np.pad(fl, ((0, 0), (pad, pad)), mode='edge')
        fl = np.array([np.convolve(r, k, mode='valid') for r in fp])
    out = np.array([np.interp(np.log(wave), lnw, r) for r in fl])
    return out[0] if single else out


# ---------------------------------------------------------------------
# 2-4) instrumental LSF + Doppler shift + resampling, as one sparse matrix
# ---------------------------------------------------------------------

def _resolution_func(R):
    """Accept a scalar, callable R(lam), or (wave, R) arrays; return R(lam)."""
    if callable(R):
        return R
    if np.ndim(R) == 0:
        return lambda lam: np.full_like(np.asarray(lam, float), float(R))
    wR, RR = R
    wR, RR = np.asarray(wR, float), np.asarray(RR, float)
    return lambda lam: np.interp(lam, wR, RR)


def lsf_matrix(wave_in, wave_out, R=None, lsf=None, lsf_offsets=None,
               lsf_halfwidth=None, rv=0.0, nsig=5.0, pixel_integrate=False,
               nsub=5):
    """
    Sparse (n_out, n_in) matrix M so that M @ flux_in = flux observed at
    wave_out, convolved with the instrumental LSF and Doppler shifted by rv.

    The LSF is defined in the OBSERVED frame at each output pixel, so a
    varying R or a per-pixel LSF is handled exactly.

    Parameters
    ----------
    R : float, callable R(lam), or (wave, R) arrays
        Gaussian LSF with FWHM = lam/R.
    lsf : callable or 2D array
        callable lsf(dlam, lam) -> kernel weights at wavelength offsets dlam
        (A) from output wavelength lam; needs lsf_halfwidth (A).
        2D array of shape (n_out, n_k): kernel for each output pixel, sampled
        at lsf_offsets (n_k, in units of OUTPUT pixels). The kernel is linearly
        interpolated between offsets, so sample it finely (<= 0.2 pixel, e.g.
        np.arange(-7, 7.01, 0.1)): at 1-pixel sampling the interpolation
        widens a ~1-pixel-sigma LSF by ~8%, which a fitted vmac absorbs.
    pixel_integrate : bool
        Also average the kernel over each output pixel's width (use when the
        LSF is given per unit wavelength but the detector pixel is not
        negligible; APOGEE-style LSF arrays already include the pixel).
    """
    wave_in = np.asarray(wave_in, float)
    wave_out = np.asarray(wave_out, float)
    # rest-frame -> observed frame for the model pixels
    win = doppler(wave_in, rv) if rv else wave_in
    dwin = np.gradient(win)                       # trapezoid-ish weights
    dwout = np.gradient(wave_out)

    if lsf is None and R is None:
        raise ValueError('give R or lsf')

    if lsf is None:
        Rf = _resolution_func(R)
        sig = wave_out / Rf(wave_out) * FWHM2SIG
        half = nsig * sig
        kern = lambda d, j: np.exp(-0.5 * (d / sig[j][:, None]) ** 2)
    elif callable(lsf):
        if lsf_halfwidth is None:
            raise ValueError('callable lsf needs lsf_halfwidth (Angstrom)')
        half = np.full(len(wave_out), float(lsf_halfwidth))
        kern = lambda d, j: lsf(d, wave_out[j][:, None])
    else:
        arr = np.asarray(lsf, float)
        if arr.shape[0] != len(wave_out) or lsf_offsets is None:
            raise ValueError('lsf array must be (n_out, n_k) with lsf_offsets of length n_k')
        offs = np.asarray(lsf_offsets, float)
        half = (np.max(np.abs(offs)) + 0.5) * np.abs(dwout)

        def kern(d, j):
            x = d / np.abs(dwout[j])[:, None]         # offset in output pixels
            out = np.empty_like(x)
            for r in range(x.shape[0]):
                out[r] = np.interp(x[r], offs, arr[j[r]], left=0.0, right=0.0)
            return out

    lo = np.searchsorted(win, wave_out - half, side='left')
    hi = np.searchsorted(win, wave_out + half, side='right')
    kmax = int(np.max(hi - lo)) if len(lo) else 0
    if kmax == 0:
        raise ValueError('LSF narrower than the model sampling -- synthesise on a finer grid')

    rows, cols, vals = [], [], []
    chunk = max(1, int(2e6 // max(kmax, 1)))
    for s in range(0, len(wave_out), chunk):
        j = np.arange(s, min(s + chunk, len(wave_out)))
        idx = lo[j][:, None] + np.arange(kmax)[None, :]
        valid = idx < hi[j][:, None]
        idx = np.minimum(idx, len(win) - 1)
        d = win[idx] - wave_out[j][:, None]
        if pixel_integrate:
            subs = (np.arange(nsub) + 0.5) / nsub - 0.5
            w = np.mean([kern(d - u * dwout[j][:, None], j) for u in subs], axis=0)
        else:
            w = kern(d, j)
        w = np.where(valid, w * dwin[idx], 0.0)
        norm = w.sum(axis=1)
        ok = norm > 0
        w[ok] /= norm[ok][:, None]
        rr, cc = np.nonzero(valid & (w != 0))
        rows.append(j[rr]); cols.append(idx[rr, cc]); vals.append(w[rr, cc])
    rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
    return sparse.csr_matrix((vals, (rows, cols)), shape=(len(wave_out), len(wave_in)))


def convolve(wave, flux, R=None, lsf=None, lsf_offsets=None, lsf_halfwidth=None,
             wave_out=None, dw_out=None, dlogw_out=None, rv=0.0,
             vmac=0.0, vsini=0.0, eps=0.6, pixel_integrate=False):
    """
    One-stop: stellar broadening -> LSF -> Doppler shift -> resample.

    With neither R nor lsf, only stellar broadening, shift and (linear)
    resampling are applied.

    Returns (wave_out, flux_out); flux_out has the same dimensionality as flux.
    """
    wave = np.asarray(wave, float)
    single = np.ndim(flux) == 1
    f2 = np.atleast_2d(np.asarray(flux, float))
    f2 = np.atleast_2d(stellar(wave, f2, vmac=vmac, vsini=vsini, eps=eps))

    if wave_out is None:
        if dw_out is None and dlogw_out is None:
            wave_out = doppler(wave, rv) if rv else wave.copy()
            # keep the output inside the shifted model range
        else:
            wr = (doppler(wave[0], rv), doppler(wave[-1], rv)) if rv else None
            wave_out = output_grid(wave, dw_out=dw_out, dlogw_out=dlogw_out, wrange=wr)
    wave_out = np.asarray(wave_out, float)

    if R is None and lsf is None:
        win = doppler(wave, rv) if rv else wave
        out = np.array([np.interp(wave_out, win, r) for r in f2])
    else:
        M = lsf_matrix(wave, wave_out, R=R, lsf=lsf, lsf_offsets=lsf_offsets,
                       lsf_halfwidth=lsf_halfwidth, rv=rv,
                       pixel_integrate=pixel_integrate)
        out = (M @ f2.T).T
    return wave_out, (out[0] if single else out)
