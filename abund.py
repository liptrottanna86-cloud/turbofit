"""
abund.py
========

Many spectra of ONE star (one atmosphere) with different abundances, for
the cost of one babsma run plus one bsyn run per abundance set.

    import abund
    # arbitrary sets of [X/M] overrides
    wave, flux, cont = abund.synthesize_sets(
        4500., 2.0, mh=-0.5, am=0.2, vmicro=1.5,
        abund_sets=[[], [['Mg', 0.4]], [['Mg', 0.4], ['C', -0.3]]],
        atmod=atm, linelists=ll, wrange=[15100, 15200], dw=0.01)

    # a ladder in one element around its base value
    xm, wave, flux, cont = abund.ladder(4500., 2.0, elem='Mg', mh=-0.5, am=0.2,
                                        offsets=[-0.4, -0.2, 0, 0.2, 0.4], ...)

By default babsma is rerun for every abundance set (rerun_babsma=True).
Reusing the first set's continuum-opacity file is ~25% faster, but bsyn also
takes the electron pressure from that file, so the ionisation balance is
frozen: in tests (4500 K giant, +-0.3 dex) reuse gave normalised-flux errors
of 0.09 for Fe (a major electron donor) and ~0.001-0.003 for Ti and Ni.
Only set rerun_babsma=False for trace elements and small offsets.

NLTE elements: if the varied element is in nlte_elems, its departure
coefficients are re-interpolated to each A(X) (cheap), as in TSFitPy.
"""

import os
import shutil
import tempfile
import numpy as np

import atomic
import synthesis

ALPHA = ('O', 'Ne', 'Mg', 'Si', 'S', 'Ar', 'Ca', 'Ti')


def base_xm(elem, am=0.0, cm=0.0, nm=0.0, elems=None):
    """[X/M] that synthesize() would use for `elem` given am/cm/nm/elems."""
    if elems:
        for e, v in elems:
            if e == elem:
                return float(v)
    if elem == 'C':
        return cm
    if elem == 'N':
        return nm
    if elem in ALPHA:
        return am
    return 0.0


def synthesize_sets(teff, logg, mh=0.0, am=0.0, cm=0.0, nm=0.0, vmicro=2.0,
                    abund_sets=([],), elems=None, wrange=(15000.0, 17000.0), dw=0.1,
                    atmod=None, atmos_type=None, dospherical=True, linelists=None,
                    solarisotopes=None, segments=None, nlte=None, nlte_elems=None,
                    rerun_babsma=True, workdir=None, save=False, c12c13=None):
    """
    Synthesize one spectrum per entry of abund_sets (each a list of
    [element, [X/M]] overrides applied on top of `elems`).

    c12c13 : None, a number (same 12C/13C for every set), or a list with one
        ratio per set (e.g. to vary only the isotope ratio, use
        abund_sets=[[]]*n with c12c13=[...]; see isotope_ladder()).

    Returns
    -------
    wave : (npix,), flux : (nset, npix), cont : (nset, npix)
    """
    if workdir is None:
        workdir = tempfile.mkdtemp(prefix='turboab')
    os.makedirs(workdir, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(workdir)
    try:
        if type(linelists) is str:
            linelists = [linelists]
        for l in linelists:
            if not os.path.exists(l):
                raise FileNotFoundError(l)

        if c12c13 is None or np.ndim(c12c13) == 0:
            ratios = [c12c13] * len(abund_sets)
        else:
            ratios = list(c12c13)
            if len(ratios) != len(abund_sets):
                raise ValueError('c12c13 list must have one entry per abundance set')
        sets_abund = []
        for s in abund_sets:
            el = list(elems or [])
            names = [e for e, _ in s]
            el = [e for e in el if e[0] not in names] + [list(e) for e in s]
            sets_abund.append(synthesis.make_abundances(mh, am, cm, nm, el))

        # atmosphere (and NLTE departures for the first set)
        nlte_info = None
        if nlte_elems:
            import nlte as nlte_mod
            atmod, nlte_info = nlte_mod.prepare(nlte, teff, logg, mh, sets_abund[0],
                                                nlte_elems, workdir='.')
            atmos_type = 'marcs_interpol'
        if atmod is None or not os.path.exists(atmod):
            raise FileNotFoundError(atmod)
        if atmos_type is None:
            import utils
            atmos_type = utils.identify_atmos(atmod)
        if solarisotopes is None:
            base = os.path.basename(atmod).lower()
            solarisotopes = base.startswith('p') if base[:1] in 'ps' else logg > 3.001
        spherical = dospherical and ('marcs' in atmos_type) and logg <= 3.001

        # which NLTE elements change between sets (need new departure files)?
        varied = set()
        for s in abund_sets:
            varied |= {e for e, _ in s}
        renlte = sorted(varied & set(nlte_elems or []))

        flux, cont, wave = [], [], None
        opac = None
        for i, ab in enumerate(sets_abund):
            if i > 0 and renlte:
                atmod, nlte_info = nlte_mod.prepare(nlte, teff, logg, mh, ab,
                                                    nlte_elems, workdir='.')
            root = 'ab{:04d}'.format(i)
            f, c, w = synthesis.do_turbospec(
                root, atmod, linelists, mh, am, ab, wrange, dw,
                solarisotopes=solarisotopes,
                babsma=None if (opac is None or rerun_babsma) else opac,
                atmos_type=atmos_type, spherical=spherical, vmicro=vmicro,
                nlte_info=nlte_info, segments=segments, c12c13=ratios[i])
            if np.isscalar(f):
                raise RuntimeError(f'Turbospectrum failed for abundance set {i}: {abund_sets[i]}')
            if opac is None:
                opac = root + 'opac'
            if wave is None:
                wave = w
            elif len(w) != len(wave) or not np.allclose(w, wave):
                f, c = np.interp(wave, w, f), np.interp(wave, w, c)
            flux.append(f)
            cont.append(c)
    finally:
        os.chdir(cwd)
        if not save:
            shutil.rmtree(workdir, ignore_errors=True)
    return wave, np.array(flux), np.array(cont)


def ladder(teff, logg, elem, offsets=(-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6),
           mh=0.0, am=0.0, cm=0.0, nm=0.0, elems=None, **kw):
    """
    Spectra with [elem/M] = base + offsets, where base is what am/cm/nm/elems
    imply for that element. Extra keywords go to synthesize_sets().

    Returns xm (the [X/M] values), wave, flux, cont.
    """
    b = base_xm(elem, am, cm, nm, elems)
    xm = b + np.asarray(offsets, float)
    sets = [[[elem, float(x)]] for x in xm]
    wave, flux, cont = synthesize_sets(teff, logg, mh=mh, am=am, cm=cm, nm=nm,
                                       elems=elems, abund_sets=sets, **kw)
    return xm, wave, flux, cont


def isotope_ladder(teff, logg, ratios=(4, 6, 8, 10, 15, 20, 30, 50, 90), **kw):
    """
    Spectra of one star with different 12C/13C ratios (everything else
    fixed). Extra keywords go to synthesize_sets() (mh, am, cm, nm, vmicro,
    elems, wrange, dw, atmod, linelists, segments, ...).

    Returns ratios, wave, flux (nratio, npix), cont.
    """
    ratios = np.asarray(ratios, float)
    wave, flux, cont = synthesize_sets(teff, logg, abund_sets=[[] for _ in ratios],
                                       c12c13=list(ratios), **kw)
    return ratios, wave, flux, cont
