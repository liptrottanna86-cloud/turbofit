"""
continuum.py
============

Generic pseudo-continuum normalisation for any spectrograph (no instrument
presets). Used to get a first normalised spectrum (RV, first guesses, line
masks). The parameter fit (starfit.py) and the line fits (linefit.py) do NOT
rely on it: they fit a smooth continuum multiplied into the model at the
same time as the stellar parameters, which avoids the usual bias of
normalising data and models differently.

    import continuum
    cont, nflux, nerr = continuum.normalise(wave, flux, err)              # automatic chunks
    cont, nflux, nerr = continuum.normalise(wave, flux, err, deg=4,
                                            chunks=[(15140, 15810), (15850, 16440), (16470, 16960)])
    chunks = continuum.find_chunks(wave)     # the automatic split (at gaps in the data)
"""

import numpy as np


def find_chunks(wave, flux=None, gap=None):
    """
    Split a spectrum into contiguous pieces at gaps (chip gaps, missing
    orders, masked regions). A gap is a jump larger than `gap` (default: 20
    times the median pixel size) or a run of non-finite flux.

    Returns a list of (wmin, wmax).
    """
    wave = np.asarray(wave, float)
    good = np.isfinite(wave)
    if flux is not None:
        good &= np.isfinite(flux)
    w = wave[good]
    if len(w) < 2:
        raise ValueError('not enough valid pixels')
    d = np.diff(w)
    g = gap if gap is not None else 20 * np.median(d)
    cut = np.nonzero(d > g)[0]
    starts = np.concatenate([[0], cut + 1])
    ends = np.concatenate([cut, [len(w) - 1]])
    return [(w[s], w[e]) for s, e in zip(starts, ends)]


def _fit_chunk(w, f, e, deg, niter, lo, hi):
    x = (w - w.mean()) / max(np.ptp(w), 1e-9) * 2
    keep = np.isfinite(f) & np.isfinite(e) & (e > 0)
    c = np.full_like(f, np.nan)
    for _ in range(niter):
        if keep.sum() <= deg + 1:
            break
        p = np.polynomial.polynomial.polyfit(x[keep], f[keep], deg, w=1 / e[keep])
        c = np.polynomial.polynomial.polyval(x, p)
        r = (f - c) / e
        s = 1.4826 * np.nanmedian(np.abs(r[keep]))
        new = np.isfinite(r) & (r > -lo * s) & (r < hi * s)
        if np.array_equal(new, keep):
            break
        keep = new
    return c


def normalise(wave, flux, err=None, chunks=None, deg=3, niter=10, lo=1.0, hi=3.0):
    """
    Iterative polynomial pseudo-continuum, separately for each chunk.

    Asymmetric sigma-clipping (default: reject points 1 sigma below / 3 sigma
    above) pushes the fit up onto the continuum and away from absorption
    lines. Choose deg according to chunk length and line density: low for
    short chunks or crowded spectra (cool stars, H band), higher for long,
    clean chunks.

    Returns (cont, flux/cont, err/cont).
    """
    wave = np.asarray(wave, float)
    flux = np.asarray(flux, float)
    err = np.ones_like(flux) if err is None else np.asarray(err, float)
    if chunks is None:
        chunks = find_chunks(wave, flux)
    cont = np.full_like(flux, np.nan)
    for w0, w1 in chunks:
        m = (wave >= w0) & (wave <= w1)
        if m.sum() > deg + 2:
            cont[m] = _fit_chunk(wave[m], flux[m], err[m], deg, niter, lo, hi)
    return cont, flux / cont, err / cont


def running_percentile(wave, flux, width, q=90.0):
    """Continuum as the q-th percentile of the flux in a running window of
    `width` Angstrom (robust, no polynomial; good for crowded spectra)."""
    wave = np.asarray(wave, float); flux = np.asarray(flux, float)
    out = np.full_like(flux, np.nan)
    lo = np.searchsorted(wave, wave - width / 2)
    hi = np.searchsorted(wave, wave + width / 2)
    for i in range(len(wave)):
        seg = flux[lo[i]:hi[i]]
        seg = seg[np.isfinite(seg)]
        if len(seg):
            out[i] = np.percentile(seg, q)
    return out
