"""
linefit.py
==========

Line-by-line abundances (plus a per-line RV offset and broadening) for one
star whose stellar parameters are already known -- typically from FERRE.

Two ways each element can be synthesised, chosen per element:

  CNO-LTE  (default)  LTE on the APOGEE MARCS grid with [C/M] and [alpha/M]
                      atmosphere nodes. The spectrum at the star's parameters
                      is interpolated from the 2^5 bracketing atmosphere
                      NODES (teff, logg, [M/H], [alpha/M], [C/M]) -- spectra
                      are interpolated, atmospheres never are. All corners are
                      synthesised with the STAR's abundances, so only the
                      atmospheric structure is being interpolated.
  NLTE                for elements in nlte_elems that have departure grids:
                      standard MARCS + departure coefficients, interpolated to
                      the exact parameters by Turbospectrum's interpolator
                      (one synthesis per abundance step, no corners).

For each element and atmosphere, a short abundance ladder is synthesised
over just the line windows (Turbospectrum segments), then every fit
iteration is pure interpolation + broadening -- no Turbospectrum calls.

Line masks
----------
* your own:      read_mask('mg_lines.txt')  (also reads TSFitPy linemasks)
* automatic:     lf.automask(['Mg','C','N','O'], star, wave_obs)
  synthesises the star (base, [X/M]+dx, and X removed), broadens to the
  instrument and keeps windows where X has a strong response and other
  species contribute little. Works for molecular features (C from CO/CN/CH,
  N from CN, O from OH) exactly as for atomic lines. Each window records its
  response and blend fraction, and write_mask() saves it for inspection.

Systematics
-----------
With star_err=dict(teff=..., logg=..., mh=..., am=..., cm=...), the CNO-LTE
mode re-weights the corner spectra it already has to teff+-sigma etc and
refits, so every line gets d[X/H]/dp and a systematic error at no extra
synthesis cost. vmicro (and, in NLTE mode, all parameters) need extra
ladders and are only done if sys_vmicro / sys_nlte are set.

Typical use
-----------
    import turbogrid, linefit, nlte
    lookup = turbogrid.make_atmod_lookup('/data/marcs_apogee', flat=True)
    lf = linefit.LineFitter(
        linelists=[...], inst=dict(R=22500),
        nodes=dict(teff=(15, 3000., 250.), logg=(12, -0.5, 0.5), mh=(14, -2.5, 0.25),
                   am=(7, -1.0, 0.25), cm=(11, -1.5, 0.25)),
        atmod_lookup=lookup,
        nlte=nlte.NLTEConfig('/data/marcs_std', '/data/nlte'), nlte_elems=['Mg', 'Ca'],
        n_jobs=16)
    star = dict(teff=4520., logg=1.9, mh=-0.42, am=0.18, cm=-0.2, nm=0.3, vmicro=1.6, rv=-121.4)
    mask = lf.automask(['Mg', 'Al', 'C', 'N', 'O'], star, wave)
    res = lf.fit(wave, flux, err, star, mask=mask,
                 star_err=dict(teff=50, logg=0.1, mh=0.05, am=0.03, cm=0.05))
    print(res.summary)
    res.write('star_lines.csv')
    linefit.plot_lines(res, 'star_lines.pdf')
"""

import os
import tempfile
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares
from scipy.signal import find_peaks

import abund
import broaden

C_KMS = 299792.458
ATM_AXES = ('teff', 'logg', 'mh', 'am', 'cm')

# pseudo-element for the carbon isotope ratio: fitted quantity is log10(12C/13C)
ISO = 'C12C13'
SOLAR_C12C13 = 0.989 / 0.011     # Turbospectrum's default (solar) mix


def default_c12c13(star):
    """12C/13C used for `star` when star['c12c13'] is not given: solar for
    logg > 3, 15 ('giant') otherwise -- the synthesize() default."""
    r = star.get('c12c13')
    if r is not None:
        return float(r)
    return SOLAR_C12C13 if star.get('logg', 3.0) > 3.001 else 15.0


def base_value(el, star):
    """Starting value of the fitted quantity: [X/M], or log10(12C/13C)."""
    if el == ISO:
        return float(np.log10(default_c12c13(star)))
    return abund.base_xm(el, star.get('am', 0.0), star.get('cm', 0.0), star.get('nm', 0.0),
                         star.get('elems'))


# =====================================================================
# line masks
# =====================================================================

def read_mask(path, element=None):
    """
    Read a line mask. Either turbofit format (from write_mask)

        # element  center  wmin  wmax  [response  blend]
        Mg  15740.716  15739.9  15741.5  0.081  0.05

    or a TSFitPy-style mask (';' comments, columns center wmin wmax), in
    which case pass element='Mg'. Wavelengths are REST-frame, in the same
    convention (air/vacuum) as your linelists.
    """
    out = []
    with open(path) as f:
        for line in f:
            s = line.split('#')[0].split(';')[0].strip()
            if not s:
                continue
            t = s.split()
            if t[0][0].isalpha():
                el, nums = t[0], [float(x) for x in t[1:]]
            else:
                if element is None:
                    raise ValueError(f'{path} has no element column: pass element=')
                el, nums = element, [float(x) for x in t]
            d = dict(element=el, center=nums[0], wmin=nums[1], wmax=nums[2],
                     response=nums[3] if len(nums) > 3 else np.nan,
                     blend=nums[4] if len(nums) > 4 else np.nan)
            out.append(d)
    return out


def write_mask(mask, path):
    with open(path, 'w') as f:
        f.write('# element  center  wmin  wmax  response  blend\n')
        for m in mask:
            f.write('{:3s} {:11.4f} {:11.4f} {:11.4f} {:8.4f} {:6.3f}\n'.format(
                m['element'], m['center'], m['wmin'], m['wmax'],
                m.get('response', np.nan), m.get('blend', np.nan)))
    return path


# =====================================================================
# atmosphere-node geometry (CNO-LTE mode)
# =====================================================================

def _axis_values(v):
    if isinstance(v, tuple) and len(v) == 3:
        n, s, d = v
        return s + d * np.arange(int(n))
    return np.atleast_1d(np.asarray(v, float))


def bracket(star, nodes):
    """For each atmosphere axis: (name, lo, hi). Exactly-on-node values still
    get a real cell [node, next node] so derivatives (systematics) exist."""
    out = []
    for ax in ATM_AXES:
        vals = np.sort(_axis_values(nodes[ax])) if ax in nodes else np.array([star.get(ax, 0.0)])
        x = float(star.get(ax, 0.0))
        if len(vals) == 1:
            out.append((ax, vals[0], vals[0]))
            continue
        if x < vals[0] - 1e-9 or x > vals[-1] + 1e-9:
            raise ValueError(f'{ax}={x} is outside the atmosphere nodes {vals[0]}..{vals[-1]}')
        i = int(np.clip(np.searchsorted(vals, x, side='right') - 1, 0, len(vals) - 2))
        out.append((ax, vals[i], vals[i + 1]))
    return out


def corner_weights(br, point):
    """Multilinear weights of the 2^n corners of cell `br` at `point` (dict);
    point is clipped into the cell. Returns {corner_tuple: weight}."""
    fr = []
    for ax, lo, hi in br:
        x = float(np.clip(point.get(ax, lo), min(lo, hi), max(lo, hi)))
        fr.append(0.0 if hi == lo else (x - lo) / (hi - lo))
    w = {}
    n = len(br)
    for v in range(2 ** n):
        key, wt = [], 1.0
        for k, (ax, lo, hi) in enumerate(br):
            bit = (v >> k) & 1
            key.append(hi if bit else lo)
            wt *= fr[k] if bit else (1.0 - fr[k])
        key = tuple(key)
        w[key] = w.get(key, 0.0) + wt
    return w


# =====================================================================
# synthesis worker (module level so it can be pickled)
# =====================================================================

def _ladder_worker(a):
    """a: dict with everything abund.synthesize_sets needs. Returns (wave, normflux)."""
    if a['elem'] == ISO:
        sets = [[] for _ in a['xm']]
        ratios = [10 ** float(x) for x in a['xm']]
    else:
        sets = [[[a['elem'], float(x)]] for x in a['xm']]
        ratios = a.get('c12c13')
    w, f, c = abund.synthesize_sets(
        a['teff'], a['logg'], mh=a['mh'], am=a['am'], cm=a['cm'], nm=a['nm'],
        vmicro=a['vmicro'], abund_sets=sets, elems=a.get('elems'),
        wrange=a['wrange'], dw=a['dw'], atmod=a['atmod'], atmos_type=a['atmos_type'],
        linelists=a['linelists'], solarisotopes=a['solarisotopes'],
        segments=a['segments'], nlte=a['nlte'], nlte_elems=a['nlte_elems'],
        c12c13=ratios, workdir=tempfile.mkdtemp(prefix='tfl', dir=a['tmpdir']))
    return w, f / c


def _sets_worker(a):
    w, f, c = abund.synthesize_sets(
        a['teff'], a['logg'], mh=a['mh'], am=a['am'], cm=a['cm'], nm=a['nm'],
        vmicro=a['vmicro'], abund_sets=a['sets'], elems=a.get('elems'),
        wrange=a['wrange'], dw=a['dw'], atmod=a['atmod'], atmos_type=a['atmos_type'],
        linelists=a['linelists'], solarisotopes=a['solarisotopes'], segments=None,
        nlte=a['nlte'], nlte_elems=a['nlte_elems'], c12c13=a.get('c12c13'),
        workdir=tempfile.mkdtemp(prefix='tfm', dir=a['tmpdir']))
    return w, f / c


def _merge_segments(wins, pad):
    segs = sorted((w0 - pad, w1 + pad) for w0, w1 in wins)
    out = []
    for s0, s1 in segs:
        if out and s0 <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], s1))
        else:
            out.append((s0, s1))
    return out


# =====================================================================
# results
# =====================================================================

class LineFitResult:
    """Per-line results (.lines, list of dicts) and per-element .summary."""

    def __init__(self, lines, star, star_err):
        self.lines = lines
        self.star = star
        self.star_err = star_err or {}
        self.summary = self._summarise()

    def _summarise(self):
        out = {}
        for el in sorted({l['element'] for l in self.lines}):
            L = [l for l in self.lines if l['element'] == el and l['ok']]
            if not L:
                out[el] = dict(n=0)
                continue
            x = np.array([l['xm'] for l in L])
            e = np.array([max(l['err_stat'], 0.01) for l in L])
            w = 1 / e ** 2
            mean = np.sum(w * x) / np.sum(w)
            # systematic on the mean: lines respond coherently to a parameter
            sys2 = 0.0
            for p, sig in self.star_err.items():
                d = [l['dxh_dp'].get(p, np.nan) for l in L]
                if np.all(np.isfinite(d)):
                    dm = np.sum(w * np.array(d)) / np.sum(w)
                    if p == 'mh' and el != ISO:
                        dm -= 1.0          # [X/M] = [X/H] - [M/H]
                    sys2 += (dm * sig) ** 2
            err_mean = np.sqrt(1 / np.sum(w))
            std = np.std(x, ddof=1) if len(L) > 1 else np.nan
            err_scatter = std / np.sqrt(len(L)) if len(L) > 1 else np.nan
            # measurement error of the mean: the larger of the propagated
            # per-line errors and the observed line-to-line scatter
            err_meas = np.nanmax([err_mean, err_scatter])
            err_sys = np.sqrt(sys2) if self.star_err else np.nan
            if el == ISO:
                # fitted quantity is log10(12C/13C); also give the ratio itself
                extra = dict(c12c13=10 ** mean, c12c13_lo=10 ** (mean - err_meas),
                             c12c13_hi=10 ** (mean + err_meas))
            else:
                extra = {}
            out[el] = dict(extra, n=len(L), xm=mean,
                           xh=np.nan if el == ISO else mean + self.star.get('mh', 0.0),
                           err_mean=err_mean, std=std, err_scatter=err_scatter,
                           err_meas=err_meas, err_sys=err_sys,
                           err_total=np.sqrt(err_meas ** 2 + (err_sys if np.isfinite(err_sys) else 0.0) ** 2),
                           mode=L[0]['mode'])
        return out

    def table(self):
        from astropy.table import Table
        cols = ['element', 'center', 'wmin', 'wmax', 'mode', 'xm', 'xh', 'err_stat',
                'err_lo', 'err_hi', 'err_mc', 'err_sys', 'drv', 'vmac', 'chi2r', 'npix', 'snr_resp', 'response',
                'blend', 'flags', 'ok']
        t = Table({c: [l[c] for l in self.lines] for c in cols})
        if any(l['element'] == ISO for l in self.lines):
            t['c12c13'] = [10 ** l['xm'] if l['element'] == ISO else np.nan for l in self.lines]
        for p in self.star_err:
            t['dxh_d' + p] = [l['dxh_dp'].get(p, np.nan) for l in self.lines]
        return t

    def write(self, path):
        self.table().write(path, overwrite=True)
        return path


# =====================================================================
# the fitter
# =====================================================================

class LineFitter:
    """
    Parameters
    ----------
    linelists : list of str
        Turbospectrum line lists (atomic, molecular, H).
    inst : dict
        Instrument for broaden.lsf_matrix: R=..., or lsf=callable +
        lsf_halfwidth=..., or lsf=(n_obs_pix, n_k) array + lsf_offsets (rows
        must match the observed pixels you will fit), pixel_integrate=...
    nodes : dict, optional
        Atmosphere-node axes for CNO-LTE mode, e.g. teff=(15,3000.,250.) or
        arrays, for 'teff','logg','mh','am','cm' (as in turbogrid).
    atmod_lookup : callable(teff, logg, mh, am, cm) -> path or None
        e.g. turbogrid.make_atmod_lookup(...).
    atmos_type : str
        Atmosphere file type for the node models (default 'marcs').
    nlte : nlte.NLTEConfig, optional ; nlte_elems : elements to do in NLTE.
    offsets : abundance-ladder offsets (dex) around the starting [X/M].
    dw : synthesis step for the full-range runs in automask (A); line windows
        use Turbospectrum's segment sampling (R=500000).
    n_jobs : parallel Turbospectrum processes.
    drv_max : per-line RV-offset bound (km/s).
    vmac, vsini : stellar broadening (km/s); vmac is the starting value / the
        value used by automask (see fit(vmac=...) for how it is fitted).
    cont_deg : per-window continuum polynomial degree (0, 1) or None (none).
    """

    def __init__(self, linelists, inst, nodes=None, atmod_lookup=None, atmos_type='marcs',
                 nlte=None, nlte_elems=(), offsets=(-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6),
                 dw=0.02, n_jobs=1, drv_max=5.0, vmac=3.0, vsini=0.0,
                 cont_deg=1, max_extend=2, workdir=None, verbose=True):
        self.linelists = list(linelists)
        self.inst = dict(inst)
        self.nodes = nodes
        self.lookup = atmod_lookup
        self.atmos_type = atmos_type
        self.nlte = nlte
        self.nlte_elems = list(nlte_elems or [])
        self.offsets = np.asarray(offsets, float)
        self.dw = dw
        self.n_jobs = n_jobs
        self.drv_max = drv_max
        self.vmac = vmac
        self.vsini = vsini
        self.cont_deg = cont_deg
        self.max_extend = max_extend
        self.tmpdir = workdir or tempfile.mkdtemp(prefix='turbofit_')
        os.makedirs(self.tmpdir, exist_ok=True)
        self.verbose = verbose
        self._cache = {}
        self._nlte_ok = set(nlte.available_elements()) if nlte is not None else set()
        if self.verbose and self.nlte_elems:
            miss = [e for e in self.nlte_elems if e not in self._nlte_ok]
            if miss:
                print(f'linefit: no NLTE data for {miss} -- they will be fitted in CNO-LTE mode')

    # ------------------------------------------------------------------
    def mode(self, elem):
        if elem == ISO:
            # molecular features: always LTE on the CNO atmospheres
            if self.nodes is not None and self.lookup is not None:
                return 'CNO-LTE'
            raise ValueError('fitting 12C/13C needs atmosphere nodes and atmod_lookup')
        if elem in self.nlte_elems and elem in self._nlte_ok:
            return 'NLTE'
        if self.nodes is not None and self.lookup is not None:
            return 'CNO-LTE'
        raise ValueError(f'{elem}: no NLTE data and no atmosphere nodes/lookup configured')

    def _lsf_half(self, lam):
        if self.inst.get('R') is not None:
            R = broaden._resolution_func(self.inst['R'])(np.atleast_1d(lam))
            return float(np.max(5 * lam / R * broaden.FWHM2SIG))
        if self.inst.get('lsf_halfwidth') is not None:
            return float(self.inst['lsf_halfwidth'])
        return 2.0

    def _pad(self, lam):
        return (self._lsf_half(lam) + lam * (self.drv_max + 3 * max(self.vmac, 10.0)
                                             + self.vsini) / C_KMS + 0.3)

    def _base_kw(self, star):
        return dict(mh=star.get('mh', 0.0), am=star.get('am', 0.0), cm=star.get('cm', 0.0),
                    nm=star.get('nm', 0.0), vmicro=star.get('vmicro', 1.5),
                    elems=star.get('elems'), linelists=self.linelists,
                    solarisotopes=star.get('logg', 3.0) > 3.001,
                    c12c13=star.get('c12c13'), nlte=self.nlte, tmpdir=self.tmpdir)

    def _atmospheres(self, star, elem):
        """[(key, atmod, atmos_type, teff_for_geometry, logg_for_geometry)], bracket."""
        if self.mode(elem) == 'NLTE':
            return [('nlte', None, None, star['teff'], star['logg'])], None
        br = bracket(star, self.nodes)
        out, missing = [], []
        for key in corner_weights(br, star):
            p = dict(zip(ATM_AXES, key))
            atm = self.lookup(p['teff'], p['logg'], p['mh'], p['am'], p['cm'])
            if not atm or not os.path.exists(atm):
                missing.append(key)
                continue
            out.append((key, atm, self.atmos_type, p['teff'], p['logg']))
        if missing:
            raise FileNotFoundError(f'no atmosphere for corner(s) {missing}; fill the grid holes '
                                    '(atmos.fill_holes) or move the star')
        return out, br

    @staticmethod
    def _ckey(star, elem, key, segs, vmicro):
        # the synthesis depends on the star's abundances (all corners use them)
        # and, for NLTE, on the exact parameters too
        ab = tuple(round(float(star.get(k, 0.0)), 4) for k in ('mh', 'am', 'cm', 'nm'))
        extra = (round(star['teff'], 2), round(star['logg'], 4)) if key == 'nlte' else ()
        return (elem, key, tuple(segs), round(vmicro, 4), ab, repr(star.get('elems')),
                star.get('logg', 3.0) > 3.001, star.get('c12c13')) + extra

    def _ladders(self, star, elem, segs, xm, atms, vmicro=None):
        """Synthesise (cached) the ladder `xm` for every atmosphere; returns
        wave, {key: (nx, npix) normalised flux}."""
        kw = self._base_kw(star)
        if vmicro is not None:
            kw['vmicro'] = vmicro
        wr = (segs[0][0] - 2.0, segs[-1][1] + 2.0)
        jobs, keys = [], []
        for key, atm, at, tg, lg in atms:
            ck = self._ckey(star, elem, key, segs, kw['vmicro'])
            have = self._cache.get(ck, {})
            todo = [x for x in xm if round(x, 4) not in have]
            if todo:
                a = dict(kw, teff=tg, logg=lg, elem=elem, xm=todo, atmod=atm, atmos_type=at,
                         wrange=wr, dw=self.dw, segments=segs,
                         nlte_elems=[elem] if key == 'nlte' else None)
                if key == 'nlte':
                    a.update(teff=star['teff'], logg=star['logg'])
                jobs.append(a); keys.append((ck, todo))
        if jobs:
            if self.verbose:
                print(f'linefit: {elem}: {len(jobs)} atmosphere(s) x {len(keys[0][1])} '
                      f'abundance steps over {len(segs)} segment(s)')
            if self.n_jobs > 1 and len(jobs) > 1:
                with ProcessPoolExecutor(max_workers=self.n_jobs) as ex:
                    res = list(ex.map(_ladder_worker, jobs))
            else:
                res = [_ladder_worker(a) for a in jobs]
            for (ck, todo), (w, f) in zip(keys, res):
                d = self._cache.setdefault(ck, {})
                d['_wave'] = d.get('_wave', w)
                for x, row in zip(todo, f):
                    if len(w) != len(d['_wave']) or not np.allclose(w, d['_wave']):
                        row = np.interp(d['_wave'], w, row)
                    d[round(x, 4)] = row
        wave, cubes = None, {}
        for key, atm, at, tg, lg in atms:
            ck = self._ckey(star, elem, key, segs, kw['vmicro'])
            d = self._cache[ck]
            if wave is None:
                wave = d['_wave']
            arr = np.array([d[round(x, 4)] for x in xm])
            if len(d['_wave']) != len(wave) or not np.allclose(d['_wave'], wave):
                arr = np.array([np.interp(wave, d['_wave'], r) for r in arr])
            cubes[key] = arr
        return wave, cubes

    # ------------------------------------------------------------------
    def automask(self, elements, star, wave_obs, dx=0.3, min_response=0.02, max_blend=0.3,
                 cut=0.25, exclude=None, wrange=None):
        """
        Choose line windows for `elements` tailored to this star and instrument.

        Parameters
        ----------
        wave_obs : observed wavelength array (observed frame; star['rv'] used).
        dx : abundance step used to measure each element's response (dex).
        min_response : minimum |flux change| for +dx dex, after instrumental
            broadening, for a window to be kept.
        max_blend : maximum fraction of the line depth (response-weighted
            over the window) contributed by OTHER species.
        cut : window edges where the response falls below cut x its peak.
        exclude : list of (w0, w1) OBSERVED-frame ranges to avoid (tellurics,
            chip gaps, bad sky residuals).
        """
        rv = star.get('rv', 0.0)
        wave_obs = np.asarray(wave_obs, float)
        if wrange is None:
            wrange = (wave_obs[0] / broaden.doppler(1.0, rv) - 2,
                      wave_obs[-1] / broaden.doppler(1.0, rv) + 2)
        # atmosphere: nearest CNO node if available, else NLTE/standard MARCS
        kw = self._base_kw(star)
        if self.nodes is not None and self.lookup is not None:
            p = {}
            for ax in ATM_AXES:
                vals = _axis_values(self.nodes[ax]) if ax in self.nodes else np.array([star.get(ax, 0.0)])
                p[ax] = float(vals[np.argmin(np.abs(vals - star.get(ax, 0.0)))])
            atm = self.lookup(p['teff'], p['logg'], p['mh'], p['am'], p['cm'])
            atype, tg, lg, nl = self.atmos_type, p['teff'], p['logg'], None
        else:
            atm, atype, tg, lg = None, None, star['teff'], star['logg']
            nl = [e for e in elements if e in self._nlte_ok] or None
            if nl is None:
                raise ValueError('automask needs atmosphere nodes, or NLTE data for at least one element')
        # 'X removed': -9 dex in LTE; in NLTE stay inside the departure grids
        # (which span about -2..+1 dex), so use -1.9 (lines ~80x weaker)
        xoff = -9.0 if nl is None else -1.9
        r0 = default_c12c13(star)
        sets, ratios = [[]], [r0]
        for el in elements:
            if el == ISO:
                # more 13C (ratio lower by dx dex) and essentially no 13C
                sets += [[], []]
                ratios += [r0 / 10 ** dx, 1e6]
                continue
            b = abund.base_xm(el, kw['am'], kw['cm'], kw['nm'], kw['elems'])
            sets += [[[el, b + dx]], [[el, xoff if (nl is None or el not in nl) else b - 1.9]]]
            ratios += [r0, r0]
        kw['c12c13'] = ratios
        # split the range over n_jobs processes
        n = max(1, self.n_jobs)
        edges = np.linspace(wrange[0], wrange[1], n + 1)
        jobs = []
        for i in range(n):
            jobs.append(dict(kw, teff=tg, logg=lg, sets=sets, atmod=atm, atmos_type=atype,
                             wrange=(edges[i] - 1.0, edges[i + 1] + 1.0), dw=self.dw,
                             nlte_elems=nl))
        if self.verbose:
            print(f'linefit.automask: {len(sets)} syntheses over {wrange[0]:.0f}-{wrange[1]:.0f} A '
                  f'in {n} chunk(s)')
        if n > 1:
            with ProcessPoolExecutor(max_workers=n) as ex:
                res = list(ex.map(_sets_worker, jobs))
        else:
            res = [_sets_worker(jobs[0])]
        ws, fs = [], []
        for i, (w, f) in enumerate(res):
            keep = (w >= edges[i]) & ((w < edges[i + 1]) if i < n - 1 else (w <= edges[i + 1]))
            ws.append(w[keep]); fs.append(f[:, keep])
        w = np.concatenate(ws); f = np.concatenate(fs, axis=1)

        inst = {k: v for k, v in self.inst.items()}
        wo, fo = broaden.convolve(w, f, wave_out=wave_obs, rv=rv, vmac=self.vmac,
                                  vsini=self.vsini, **inst)
        base = fo[0]
        bad = np.zeros(len(wave_obs), bool)
        for a, b in (exclude or []):
            bad |= (wave_obs >= a) & (wave_obs <= b)
        # model can't be trusted where it was extrapolated beyond the synthesis range
        bad |= (wave_obs < broaden.doppler(w[0], rv) + 1) | (wave_obs > broaden.doppler(w[-1], rv) - 1)

        dpix = np.median(np.diff(wave_obs))
        fwhm_pix = max(2, int(round(self._lsf_half(np.median(wave_obs)) / 5 * 2.3548 / dpix)))
        mask = []
        for k, el in enumerate(elements):
            fplus, fnox = fo[1 + 2 * k], fo[2 + 2 * k]
            r = np.where(bad, 0.0, base - fplus)
            pk, _ = find_peaks(r, height=min_response, distance=max(2, fwhm_pix))
            wins = []
            for i in pk:
                lo = i
                while lo > 0 and r[lo - 1] > cut * r[i] and not bad[lo - 1]:
                    lo -= 1
                hi = i
                while hi < len(r) - 1 and r[hi + 1] > cut * r[i] and not bad[hi + 1]:
                    hi += 1
                lo, hi = max(lo - 1, 0), min(hi + 1, len(r) - 1)
                tot = np.clip(1 - base[lo:hi + 1], 1e-6, None)
                own = np.clip(fnox[lo:hi + 1] - base[lo:hi + 1], 0, None)
                oth = np.clip((tot - own) / tot, 0, 1)
                rr = np.clip(r[lo:hi + 1], 0, None)
                blend = float(np.sum(rr * oth) / np.sum(rr))
                wins.append([lo, hi, i, float(r[i]), blend])
            # merge overlapping windows
            wins.sort()
            merged = []
            for wv in wins:
                if merged and wv[0] <= merged[-1][1]:
                    m = merged[-1]
                    if wv[3] > m[3]:
                        m[2], m[3] = wv[2], wv[3]
                    m[4] = max(m[4], wv[4])
                    m[1] = max(m[1], wv[1])
                else:
                    merged.append(list(wv))
            for lo, hi, i, resp, blend in merged:
                if blend > max_blend:
                    continue
                f_ = broaden.doppler(1.0, rv)
                mask.append(dict(element=el, center=wave_obs[i] / f_, wmin=wave_obs[lo] / f_,
                                 wmax=wave_obs[hi] / f_, response=resp, blend=blend))
            if self.verbose:
                nk = sum(1 for m in mask if m['element'] == el)
                print(f'linefit.automask: {el}: {len(merged)} candidate windows, {nk} kept '
                      f'(blend <= {max_blend})')
        return mask

    # ------------------------------------------------------------------
    def _window_data(self, wobs, fobs, eobs, win, rv):
        f_ = broaden.doppler(1.0, rv)
        sel = (wobs >= win['wmin'] * f_) & (wobs <= win['wmax'] * f_) & np.isfinite(fobs) \
            & np.isfinite(eobs) & (eobs > 0)
        idx = np.nonzero(sel)[0]
        inst = dict(self.inst)
        if inst.get('lsf') is not None and not callable(inst['lsf']):
            inst['lsf'] = np.asarray(inst['lsf'])[idx]
        return idx, inst

    def _fit_one(self, wm, cube, xm, wo, fo, eo, inst, rv, x0, vm0, fit_vmac, fixed=None,
                 profile=False, mc=0, rng=None):
        """Fit one window. cube: (nx, npix) normalised model on rest wave wm.

        profile : also return asymmetric 1-sigma errors from the chi^2 profile
            (abundance stepped, RV offset [and vmac] re-optimised at each step).
        mc : number of Monte Carlo refits with Gaussian noise (from eo) added
            to the data; returns their standard deviation as err_mc.
        """
        lam = np.median(wo)
        pad = self._pad(lam)
        f_ = broaden.doppler(1.0, rv)
        cols = (wm >= wo[0] / f_ - pad) & (wm <= wo[-1] / f_ + pad)
        ws, cs = wm[cols], cube[:, cols]
        spl = CubicSpline(xm, cs, axis=0)
        t = (wo - lam) / lam
        deg = self.cont_deg

        def model(p, fdat=fo):
            x, drv, vm = p[0], p[1], (p[2] if fit_vmac else vm0)
            s = spl(x)
            if vm > 0 or self.vsini > 0:
                s = broaden.stellar(ws, s, vmac=vm, vsini=self.vsini)
            M = broaden.lsf_matrix(ws, wo, rv=rv + drv, **inst)
            m = M @ s
            if deg is None:
                return m
            A = np.vstack([m * t ** k for k in range(deg + 1)]).T
            c, *_ = np.linalg.lstsq(A / eo[:, None], fdat / eo, rcond=None)
            return A @ c

        def resid(p, fdat=fo):
            return (fdat - model(p, fdat)) / eo

        lo = [xm[0], -self.drv_max] + ([0.0] if fit_vmac else [])
        hi = [xm[-1], self.drv_max] + ([30.0] if fit_vmac else [])
        xs = [0.1, 1.0] + ([1.0] if fit_vmac else [])
        ds = [1e-3, 1e-3] + ([1e-3] if fit_vmac else [])
        p0 = [np.clip(x0, xm[0] + 1e-3, xm[-1] - 1e-3), 0.0] + ([np.clip(vm0, 0.1, 29)] if fit_vmac else [])
        if fixed is not None:
            # refit only the abundance, other parameters frozen (systematics)
            def r1(q):
                return resid([q[0]] + list(fixed[1:]))
            r = least_squares(r1, [fixed[0]], bounds=([xm[0]], [xm[-1]]), x_scale=[0.1])
            return r.x[0]
        r = least_squares(resid, p0, bounds=(lo, hi), x_scale=xs, diff_step=ds)
        chi2 = 2 * r.cost
        dof = max(len(fo) - len(p0) - (0 if deg is None else deg + 1), 1)
        scale = max(1.0, chi2 / dof)
        try:
            cov = np.linalg.inv(r.jac.T @ r.jac) * scale
            ex = float(np.sqrt(cov[0, 0]))
        except np.linalg.LinAlgError:
            ex = np.nan
        # how detectable is 0.1 dex here?
        p_up = list(r.x); p_up[0] = min(r.x[0] + 0.1, xm[-1])
        p_dn = list(r.x); p_dn[0] = max(r.x[0] - 0.1, xm[0])
        dm = (model(p_up) - model(p_dn)) / max(p_up[0] - p_dn[0], 1e-6) * 0.1
        snr_resp = float(np.sqrt(np.sum((dm / eo) ** 2)))
        out = dict(x=r.x, xm=float(r.x[0]), err=ex, chi2r=chi2 / dof, model=model(r.x),
                   snr_resp=snr_resp, err_lo=np.nan, err_hi=np.nan, err_mc=np.nan)

        if profile:
            # chi^2 profile in [X/M]: nuisance parameters re-optimised at each step
            xb = float(r.x[0])
            s0 = ex if np.isfinite(ex) and ex > 0 else 0.05
            grid = np.unique(np.clip(xb + s0 * np.linspace(-4, 4, 17), xm[0], xm[-1]))
            prof = []
            for xg in grid:
                if len(r.x) > 1:
                    rr = least_squares(lambda q: resid([xg] + list(q)), r.x[1:],
                                       bounds=(lo[1:], hi[1:]), x_scale=xs[1:], diff_step=ds[1:])
                    prof.append(2 * rr.cost)
                else:
                    prof.append(float(np.sum(resid([xg]) ** 2)))
            dchi = (np.array(prof) - chi2) / scale
            dchi = np.maximum(dchi, 0.0)

            def crossing(side):
                sel = grid >= xb if side > 0 else grid <= xb
                g, d = grid[sel], dchi[sel]
                o = np.argsort(np.abs(g - xb))
                g, d = g[o], d[o]
                for k in range(1, len(g)):
                    if d[k] >= 1.0:
                        # linear interpolation between the bracketing points
                        f = (1.0 - d[k - 1]) / max(d[k] - d[k - 1], 1e-12)
                        return abs(g[k - 1] + f * (g[k] - g[k - 1]) - xb)
                return np.inf    # never rose by 1 within the ladder / scan
            out['err_lo'], out['err_hi'] = crossing(-1), crossing(+1)

        if mc and mc > 0:
            rng = rng or np.random.default_rng()
            xs_mc = []
            for _ in range(int(mc)):
                fk = fo + rng.normal(0.0, eo)
                rk = least_squares(lambda q: resid(q, fk), r.x, bounds=(lo, hi),
                                   x_scale=xs, diff_step=ds)
                xs_mc.append(rk.x[0])
            out['err_mc'] = float(np.std(xs_mc, ddof=1))
        return out

    # ------------------------------------------------------------------
    def fit(self, wave, flux, err, star, elements=None, mask=None, star_err=None,
            vmac='global', iterate=1, sys_vmicro=False, sys_nlte=False, min_pix=4,
            profile=True, mc=0, seed=None):
        """
        Fit every window in `mask` (default: automask for `elements`).

        star : dict teff, logg, mh, am, cm, nm, vmicro, rv (+ optional
            elems=[[El, [X/M]], ...] for known abundances of blending species)
        star_err : dict of 1-sigma errors for systematics, any of teff, logg,
            mh, am, cm (CNO-LTE: free), vmicro (needs sys_vmicro=True), and
            for NLTE elements teff/logg/mh/vmicro (needs sys_nlte=True).
        vmac : 'global' (default) -- fit vmac per line in a first pass, then
            fix it to the response-weighted median of the good lines; 'line'
            -- free per line (degenerate with abundance for weak lines in
            narrow windows); or a number to fix it.
        profile : asymmetric measurement errors err_lo / err_hi from the
            chi^2 profile of each line (final pass only; no extra synthesis).
        mc : number of Monte Carlo noise realisations per line (final pass
            only) -> err_mc. 0 (default) to skip; ~100 is plenty.
        seed : random seed for the Monte Carlo errors.
        iterate : extra passes in which the abundances just measured for the
            fitted elements are fed back into the synthesis (so blends from
            one fitted element into another's windows are modelled with
            the right abundance). Each pass re-synthesises the ladders.
        """
        star = dict(star)
        elems = list(star.get('elems') or [])
        if mask is None:
            if elements is None:
                raise ValueError('give elements or mask')
            mask = self.automask(elements, star, wave)
        if vmac == 'global':
            res = self._fit_pass(wave, flux, err, star, elements, mask, {}, True, self.vmac,
                                 False, False, min_pix)
            v = np.array([l['vmac'] for l in res.lines if l['ok']])
            wt = np.array([l['snr_resp'] for l in res.lines if l['ok']])
            if len(v):
                o = np.argsort(v); c = np.cumsum(wt[o]) / np.sum(wt)
                vglob = float(v[o][np.searchsorted(c, 0.5)])
            else:
                vglob = self.vmac
            if self.verbose:
                print(f'linefit: global vmac = {vglob:.2f} km/s from {len(v)} lines')
            fit_vmac, vm = False, vglob
        elif vmac == 'line':
            fit_vmac, vm = True, self.vmac
        else:
            fit_vmac, vm = False, float(vmac)
        npass = 1 + int(iterate)
        for ip in range(npass):
            last = ip == npass - 1
            res = self._fit_pass(wave, flux, err, star, elements, mask,
                                 star_err if last else {}, fit_vmac, vm,
                                 sys_vmicro and last, sys_nlte and last, min_pix,
                                 profile=profile and last, mc=mc if last else 0, seed=seed)
            if last:
                break
            # feed the measured abundances back in for the next pass
            done = {el: v['xm'] for el, v in res.summary.items() if v.get('n', 0) > 0}
            if ISO in done:
                star['c12c13'] = float(10 ** done.pop(ISO))
            elems = [e for e in elems if e[0] not in done] + [[k, float(v)] for k, v in done.items()]
            star['elems'] = elems
            if self.verbose:
                print('linefit: pass {} abundances '.format(ip + 1) +
                      ' '.join(f'{k}={v:+.3f}' for k, v in done.items()))
        res.vmac = vm
        return res

    def _fit_pass(self, wave, flux, err, star, elements, mask, star_err, fit_vmac, vmac0,
                  sys_vmicro, sys_nlte, min_pix, profile=False, mc=0, seed=None):
        rng = np.random.default_rng(seed)
        wave = np.asarray(wave, float); flux = np.asarray(flux, float); err = np.asarray(err, float)
        star = dict(star)
        rv = star.get('rv', 0.0)
        if mask is None:
            if elements is None:
                raise ValueError('give elements or mask')
            mask = self.automask(elements, star, wave)
        if elements is None:
            elements = sorted({m['element'] for m in mask})
        star_err = dict(star_err or {})
        lines = []
        for el in elements:
            wins = [m for m in mask if m['element'] == el]
            if not wins:
                continue
            md = self.mode(el)
            lam = np.median([w['center'] for w in wins])
            segs = _merge_segments([(w['wmin'], w['wmax']) for w in wins], self._pad(lam) + 0.5)
            atms, br = self._atmospheres(star, el)
            b0 = base_value(el, star)
            xm = np.round(b0 + self.offsets, 4)
            step = np.median(np.diff(self.offsets))

            for it in range(self.max_extend + 1):
                wm, cubes = self._ladders(star, el, segs, list(xm), atms)
                if md == 'NLTE':
                    W = {'nlte': 1.0}
                else:
                    W = corner_weights(br, star)
                cube = sum(W[k] * cubes[k] for k in W)
                results = []
                need_lo = need_hi = False
                for win in wins:
                    idx, inst = self._window_data(wave, flux, err, win, rv)
                    if len(idx) < min_pix:
                        results.append(None)
                        continue
                    r = self._fit_one(wm, cube, xm, wave[idx], flux[idx], err[idx], inst, rv,
                                      b0, vmac0, fit_vmac, profile=profile, mc=mc, rng=rng)
                    results.append((idx, inst, r))
                    if r['xm'] <= xm[0] + 0.02:
                        need_lo = True
                    if r['xm'] >= xm[-1] - 0.02:
                        need_hi = True
                if not (need_lo or need_hi) or it == self.max_extend:
                    break
                add = []
                if need_lo:
                    add += list(xm[0] - step * np.arange(1, 3))
                if need_hi:
                    add += list(xm[-1] + step * np.arange(1, 3))
                xm = np.round(np.sort(np.concatenate([xm, add])), 4)
                if self.verbose:
                    print(f'linefit: {el}: extending ladder to {xm[0]:+.2f}..{xm[-1]:+.2f}')

            # ---- systematics -------------------------------------------
            sys_cubes = {}       # p -> (cube_plus, cube_minus, h_plus, h_minus)
            if star_err and md == 'CNO-LTE':
                for p, sig in star_err.items():
                    if p not in ATM_AXES:
                        continue
                    ax = [b for b in br if b[0] == p][0]
                    if ax[1] == ax[2]:
                        continue
                    pp, pm = dict(star), dict(star)
                    pp[p] = min(star[p] + sig, max(ax[1], ax[2]))
                    pm[p] = max(star[p] - sig, min(ax[1], ax[2]))
                    if pp[p] - pm[p] <= 0:
                        continue
                    Wp, Wm = corner_weights(br, pp), corner_weights(br, pm)
                    sys_cubes[p] = (sum(Wp[k] * cubes[k] for k in Wp),
                                    sum(Wm[k] * cubes[k] for k in Wm),
                                    pp[p] - star[p], star[p] - pm[p])
            if star_err and md == 'NLTE' and sys_nlte:
                for p in ('teff', 'logg', 'mh'):
                    if p not in star_err:
                        continue
                    pp, pm = dict(star), dict(star)
                    pp[p] += star_err[p]; pm[p] -= star_err[p]
                    _, cp = self._ladders(pp, el, segs, list(xm), [('nlte', None, None, pp['teff'], pp['logg'])])
                    _, cm_ = self._ladders(pm, el, segs, list(xm), [('nlte', None, None, pm['teff'], pm['logg'])])
                    sys_cubes[p] = (cp['nlte'], cm_['nlte'], star_err[p], star_err[p])
            if star_err and 'vmicro' in star_err and ((md == 'CNO-LTE' and sys_vmicro) or
                                                       (md == 'NLTE' and sys_nlte)):
                # vmicro derivative from the nearest atmosphere only (cheap)
                if md == 'NLTE':
                    near = [('nlte', None, None, star['teff'], star['logg'])]
                    k0 = 'nlte'
                else:
                    k0 = max(W, key=W.get)
                    near = [a for a in atms if a[0] == k0]
                v0, sv = star.get('vmicro', 1.5), star_err['vmicro']
                _, cvp = self._ladders(star, el, segs, list(xm), near, vmicro=v0 + sv)
                _, cvm = self._ladders(star, el, segs, list(xm), near, vmicro=max(v0 - sv, 0.0))
                sys_cubes['vmicro'] = (cube + cvp[k0] - cubes[k0], cube + cvm[k0] - cubes[k0],
                                       sv, v0 - max(v0 - sv, 0.0))

            for win, rr in zip(wins, results):
                row = dict(element=el, center=win['center'], wmin=win['wmin'], wmax=win['wmax'],
                           mode=md, response=win.get('response', np.nan),
                           blend=win.get('blend', np.nan), dxh_dp={}, flags='')
                if rr is None:
                    row.update(xm=np.nan, xh=np.nan, err_stat=np.nan, err_sys=np.nan,
                               err_lo=np.nan, err_hi=np.nan, err_mc=np.nan, drv=np.nan,
                               vmac=np.nan, chi2r=np.nan, npix=0, snr_resp=np.nan, ok=False,
                               flags='nodata', wave=None, flux=None, model=None)
                    lines.append(row)
                    continue
                idx, inst, r = rr
                flags = []
                if r['xm'] <= xm[0] + 0.02 or r['xm'] >= xm[-1] - 0.02:
                    flags.append('edge')
                # 'weak': 0.1 dex in [X/M] moves the window by < 3 sigma. For
                # log(12C/13C) the natural step is 0.3 dex (a factor 2 in the
                # ratio), i.e. the same test at 1 sigma per 0.1 dex
                if r['snr_resp'] < (1.0 if el == ISO else 3.0):
                    flags.append('weak')
                if r['chi2r'] > 5:
                    flags.append('poorfit')
                sys2 = 0.0
                for p, (cp, cmn, hp, hm) in sys_cubes.items():
                    xp = self._fit_one(wm, cp, xm, wave[idx], flux[idx], err[idx], inst, rv,
                                       r['xm'], vmac0, fit_vmac, fixed=list(r['x']))
                    xn = self._fit_one(wm, cmn, xm, wave[idx], flux[idx], err[idx], inst, rv,
                                       r['xm'], vmac0, fit_vmac, fixed=list(r['x']))
                    d = (xp - xn) / (hp + hm)          # d[X/M]/dp at fixed [M/H] synthesis
                    row['dxh_dp'][p] = d
                    dd = d - 1.0 if (p == 'mh' and el != ISO) else d    # error on [X/M]
                    sys2 += (dd * star_err[p]) ** 2
                row.update(xm=r['xm'], xh=(np.nan if el == ISO else r['xm'] + star.get('mh', 0.0)),
                           err_stat=r['err'],
                           err_lo=r['err_lo'], err_hi=r['err_hi'], err_mc=r['err_mc'],
                           err_sys=np.sqrt(sys2) if sys_cubes else np.nan,
                           drv=float(r['x'][1]), vmac=float(r['x'][2]) if fit_vmac else vmac0,
                           chi2r=r['chi2r'], npix=len(idx), snr_resp=r['snr_resp'],
                           ok=('edge' not in flags and 'weak' not in flags),
                           flags=','.join(flags), wave=wave[idx], flux=flux[idx],
                           model=r['model'])
                lines.append(row)
        return LineFitResult(lines, star, {k: v for k, v in star_err.items()
                                           if k in ATM_AXES or (k == 'vmicro' and (sys_vmicro or sys_nlte))})


# =====================================================================
# plotting
# =====================================================================

def plot_lines(res, outfile=None, ncols=4, element=None, colors=('k', 'r')):
    """Observed (colors[0], default black) vs best-fit model (colors[1],
    default red) for every fitted window, in the observed frame."""
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    rc = {'font.family': 'serif', 'font.size': 14, 'axes.labelsize': 14,
          'axes.titlesize': 14, 'legend.fontsize': 10,
          'xtick.direction': 'in', 'ytick.direction': 'in',
          'xtick.top': True, 'ytick.right': True,
          'xtick.minor.visible': True, 'ytick.minor.visible': True}
    L = [l for l in res.lines if l['wave'] is not None and (element is None or l['element'] == element)]
    if not L:
        return None
    nrows = int(np.ceil(len(L) / ncols))
    with matplotlib.rc_context(rc):
        fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows), squeeze=False)
        for ax, l in zip(axes.flat, L):
            ax.plot(l['wave'], l['flux'], '.', color=colors[0], ms=4, label='data')
            ax.plot(l['wave'], l['model'], '-', color=colors[1], lw=1.5, label='fit')
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.set_title('{} {:.2f}'.format(l['element'], l['center']))
            if l['element'] == ISO:
                lab = '$^{{12}}$C/$^{{13}}$C={:.1f}'.format(10 ** l['xm'])
            else:
                lab = '[{}/M]={:+.2f}$\\pm${:.2f}'.format(l['element'], l['xm'], l['err_stat'])
            ax.text(0.03, 0.06, lab + (' ' + l['flags'] if l['flags'] else ''),
                transform=ax.transAxes, fontsize=10)
        for ax in list(axes.flat)[len(L):]:
            ax.set_visible(False)
        axes[0, 0].legend(loc='lower right')
        for ax in axes[-1]:
            ax.set_xlabel(r'Observed wavelength ($\AA$)')
        for ax in axes[:, 0]:
            ax.set_ylabel('Normalised flux')
        fig.tight_layout()
        if outfile:
            fig.savefig(outfile)
    return fig
