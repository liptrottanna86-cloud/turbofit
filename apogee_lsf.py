"""
apogee_lsf.py
=============

The APOGEE line-spread function (Gauss-Hermite core + wings, varying along
each detector, per fiber), packaged for turbofit's broaden.py. Use it
anywhere turbofit takes an instrument: grids, line fits, RV templates.

    import apogee_lsf, broaden, linefit
    alsf = apogee_lsf.ApogeeLSF('/data/apogee_lsf')       # apLSF-*.fits, apWave-*.fits
    wave = apogee_lsf.apStarWavegrid()                     # 8575 log-lambda pixels

    inst = alsf.inst(wave, fiber='combo')                  # dict(lsf=..., lsf_offsets=...)
    w, f = broaden.convolve(wsyn, fsyn, wave_out=wave, **inst)

    lf = linefit.LineFitter(..., inst=alsf.inst(wave_obs, fiber=150))
    tg.write_grid_table(..., broaden_kw=dict(wave_out=wave, **alsf.inst(wave)))

The LSF is evaluated directly at each OUTPUT pixel (any wavelength grid,
not only apStar), at sub-pixel offsets, and already includes the detector
pixel (the Gauss-Hermite profile is integrated over each pixel), so do not
also set pixel_integrate=True.

Needs, in one directory (from the APOGEE data products / apogee package):
    apLSF-a-05440020.fits  apLSF-b-05440020.fits  apLSF-c-05440020.fits
    apWave-a-02420038.fits apWave-b-02420038.fits apWave-c-02420038.fits
(other LSF/wavelength-solution IDs can be passed in).

The LSF evaluation (raw(), unpack_lsf_params(), wave2pix() and the binned
Gauss-Hermite/wing integrals) is adapted from apogee.spec.lsf by Jo Bovy
(https://github.com/jobovy/apogee, BSD-3), itself based on code by Jon
Holtzman and David Nidever in the APOGEE reduction pipeline.
"""

import os
import math
import numpy
from scipy import special, interpolate

_SQRTTWO = numpy.sqrt(2.)
_LOG10LAMBDA0 = 4.179
_DLOG10LAMBDA = 6.e-6
_NLAMBDA = 8575


def apStarWavegrid():
    """The apStar log-lambda wavelength grid (vacuum, 8575 pixels)."""
    return 10. ** (_LOG10LAMBDA0 + _DLOG10LAMBDA * numpy.arange(_NLAMBDA))


class ApogeeLSF:
    """
    Parameters
    ----------
    datadir : directory with the apLSF-[abc]-<lsfid>.fits and
        apWave-[abc]-<waveid>.fits files (default: $APOGEE_LSF_DIR).
    lsfid, waveid : file IDs (defaults as in apogee.spec.lsf).
    """

    def __init__(self, datadir=None, lsfid='05440020', waveid='02420038'):
        from astropy.io import fits
        datadir = datadir or os.environ.get('APOGEE_LSF_DIR')
        if not datadir:
            raise ValueError('give datadir= or set APOGEE_LSF_DIR')
        self.lsfpars, self.wavepix = {}, {}
        for chip in 'abc':
            with fits.open(os.path.join(datadir, f'apLSF-{chip}-{lsfid}.fits')) as h:
                self.lsfpars[chip] = numpy.array(h[0].data, dtype=float)
            with fits.open(os.path.join(datadir, f'apWave-{chip}-{waveid}.fits')) as h:
                self.wavepix[chip] = numpy.array(h[2].data, dtype=float)

    # -- detector geometry -------------------------------------------
    def wave2pix(self, wave, chip, fiber=300):
        """Wavelength -> detector pixel for one chip/fiber (NaN off the chip)."""
        wave0 = self.wavepix[chip][300 - fiber]
        pix0 = numpy.arange(len(wave0))
        s = numpy.argsort(wave0)
        wave0, pix0 = wave0[s], pix0[s]
        base = numpy.polynomial.Polynomial.fit(wave0, pix0, 1)
        ip = interpolate.InterpolatedUnivariateSpline(wave0, pix0 / base(wave0), k=3)
        wave = numpy.atleast_1d(numpy.asarray(wave, float))
        out = base(wave) * ip(wave)
        out[(wave > wave0[-1]) | (wave < wave0[0])] = numpy.nan
        return out

    # -- the LSF ------------------------------------------------------
    def lsf_array(self, wave_out, fiber='combo', offsets=None):
        """
        LSF at every pixel of wave_out, sampled at `offsets` (in OUTPUT pixels).

        fiber : 'combo' (average of fibers 50,100,...,300, as in apogee.spec.lsf),
            an int 1..300, or a list of fibers to average.
        offsets : default numpy.arange(-7, 7.001, 0.1).

        Returns (lsf, offsets): lsf has shape (len(wave_out), len(offsets)),
        each row normalised to unit sum. Rows off the detectors (chip gaps,
        beyond the ends) are filled with the nearest valid row.
        """
        if isinstance(fiber, str) and fiber.lower() == 'combo':
            fibers = [50, 100, 150, 200, 250, 300]
        elif isinstance(fiber, (int, numpy.integer)):
            fibers = [int(fiber)]
        else:
            fibers = [int(f) for f in fiber]
        x = numpy.arange(-7., 7.001, 0.1) if offsets is None else numpy.asarray(offsets, float)
        wave_out = numpy.asarray(wave_out, float)
        nout = len(wave_out)
        # wavelength one output pixel further on (for the pixel-scale conversion)
        dl = numpy.empty(nout)
        dl[:-1] = wave_out[1:]
        dl[-1] = wave_out[-1] * wave_out[-1] / wave_out[-2]
        out = numpy.zeros((nout, len(x)))
        for chip in 'abc':
            for fib in fibers:
                pix = self.wave2pix(wave_out, chip, fib)
                dpix = self.wave2pix(dl, chip, fib) - pix      # detector px per output px
                gd = numpy.isfinite(pix) & numpy.isfinite(dpix)
                if not gd.any():
                    continue
                xs = x[None, :] * dpix[gd][:, None]             # offsets in detector px
                out[gd] += raw(xs, pix[gd], self.lsfpars[chip][:, 300 - fib])
        out[out < 0] = 0.
        s = out.sum(axis=1)
        good = s > 0
        if not good.any():
            raise ValueError('wave_out does not fall on the APOGEE detectors')
        out[good] /= s[good][:, None]
        if (~good).any():
            ig = numpy.nonzero(good)[0]
            for i in numpy.nonzero(~good)[0]:
                out[i] = out[ig[numpy.argmin(numpy.abs(ig - i))]]
        return out, x

    def inst(self, wave_out, fiber='combo', offsets=None):
        """dict(lsf=..., lsf_offsets=...) for broaden / LineFitter(inst=...) /
        turbogrid broaden_kw (add wave_out=wave_out there)."""
        lsf, x = self.lsf_array(wave_out, fiber=fiber, offsets=offsets)
        return dict(lsf=lsf, lsf_offsets=x)


# =====================================================================
# raw LSF evaluation (adapted from apogee.spec.lsf, J. Bovy, BSD-3)
# =====================================================================

def raw(x, xcenter, params):
    """Raw APOGEE LSF on the native detector-pixel scale at offsets x
    ([ncenter, nx] or 1D) around detector positions xcenter."""
    if len(x.shape) == 1:
        x = numpy.tile(x, (len(xcenter), 1))
    params = unpack_lsf_params(params)
    wingparams = numpy.empty((params['nWpar'], len(xcenter)))
    for ii in range(params['nWpar']):
        poly = numpy.polynomial.Polynomial(params['Wcoefs'][ii])
        wingparams[ii] = poly(xcenter + params['Xoffset'])
    ghparams = numpy.empty((params['Horder'] + 2, len(xcenter)))
    for ii in range(params['Horder'] + 2):
        if ii == 1:
            ghparams[ii] = 1. - wingparams[0]            # fixed, corrected for the wing
        else:
            poly = numpy.polynomial.Polynomial(params['GHcoefs'][ii - (ii > 1)])
            ghparams[ii] = poly(xcenter + params['Xoffset'])
        if ii > 0:
            ghparams[ii] /= numpy.sqrt(2. * numpy.pi * math.factorial(ii - 1))
    out = _gausshermitebin(x, ghparams, params['binsize'])
    out += _wingsbin(x, wingparams, params['binsize'], params['Wproftype'])
    return out


def _gausshermitebin(x, params, binsize):
    """Gauss-Hermite profile integrated over each pixel."""
    ncenter = params.shape[1]
    out = numpy.empty((ncenter, x.shape[1]))
    integ = numpy.empty((params.shape[0] - 1, x.shape[1]))
    for ii in range(ncenter):
        poly = numpy.polynomial.HermiteE(params[1:, ii])
        poly = poly.convert(kind=numpy.polynomial.Polynomial)
        coef = numpy.zeros(params.shape[0] - 1)
        coef[:len(poly.coef)] = poly.coef[:len(coef)]          # convert() may trim zeros
        w1 = (x[ii] - 0.5 * binsize) / params[0, ii]
        w2 = (x[ii] + 0.5 * binsize) / params[0, ii]
        eexp1 = numpy.exp(-0.5 * w1 ** 2.)
        eexp2 = numpy.exp(-0.5 * w2 ** 2.)
        integ[0] = numpy.sqrt(numpy.pi / 2.) * (special.erf(w2 / _SQRTTWO) - special.erf(w1 / _SQRTTWO))
        out[ii] = coef[0] * integ[0]
        if params.shape[0] > 2:
            integ[1] = -eexp2 + eexp1
            out[ii] += coef[1] * integ[1]
        for jj in range(2, params.shape[0] - 1):
            integ[jj] = (-w2 ** (jj - 1) * eexp2 + w1 ** (jj - 1) * eexp1) + (jj - 1) * integ[jj - 2]
            out[ii] += coef[jj] * integ[jj]
    return out


def _wingsbin(x, params, binsize, Wproftype):
    """LSF wings integrated over each pixel (Gaussian wings, Wproftype 1)."""
    ncenter = params.shape[1]
    out = numpy.zeros((ncenter, x.shape[1]))
    if Wproftype == 1:
        for ii in range(ncenter):
            w1 = (x[ii] - 0.5 * binsize) / params[1, ii]
            w2 = (x[ii] + 0.5 * binsize) / params[1, ii]
            out[ii] = params[0, ii] / 2. * (special.erf(w2 / _SQRTTWO) - special.erf(w1 / _SQRTTWO))
    return out


def unpack_lsf_params(lsfarr):
    """Unpack an APOGEE LSF parameter array (see apogee.spec.lsf)."""
    out = {}
    out['binsize'] = lsfarr[0]
    out['Xoffset'] = lsfarr[1]
    out['Horder'] = int(lsfarr[2])
    out['Porder'] = lsfarr[3:out['Horder'] + 4].astype('int')
    nGHcoefs = numpy.sum(out['Porder'] + 1)
    maxPorder = numpy.amax(out['Porder'])
    GHcoefs = numpy.zeros((out['Horder'] + 1, maxPorder + 1))
    GHpar = lsfarr[out['Horder'] + 4:out['Horder'] + 4 + nGHcoefs]
    CoeffStart = numpy.hstack((0, numpy.cumsum(out['Porder'] + 1)))
    for ii in range(out['Horder'] + 1):
        GHcoefs[ii, :out['Porder'][ii] + 1] = GHpar[CoeffStart[ii]:CoeffStart[ii] + out['Porder'][ii] + 1]
    out['GHcoefs'] = GHcoefs
    wingarr = lsfarr[3 + out['Horder'] + 1 + nGHcoefs:]
    out['Wproftype'] = int(wingarr[0])
    out['nWpar'] = int(wingarr[1])
    out['WPorder'] = wingarr[2:2 + out['nWpar']].astype('int')
    maxWPorder = numpy.amax(out['WPorder'])
    Wcoefs = numpy.zeros((out['nWpar'], maxWPorder + 1))
    Wpar = wingarr[out['nWpar'] + 2:]
    WingcoeffStart = numpy.hstack((0, numpy.cumsum(out['WPorder'] + 1)))
    for ii in range(out['nWpar']):
        Wcoefs[ii, :out['WPorder'][ii] + 1] = Wpar[WingcoeffStart[ii]:WingcoeffStart[ii] + out['WPorder'][ii] + 1]
    out['Wcoefs'] = Wcoefs
    return out
