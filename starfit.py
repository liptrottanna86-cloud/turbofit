"""
starfit.py
==========

LOCAL-MODE stellar parameters, spectrum by spectrum, with no precomputed
grid: an ASPCAP-style full-spectrum fit of

    Teff, log g, [M/H], [alpha/M], [C/M], [N/M], vmicro, broadening (vmac
    and/or vsini) -- plus RV and a smooth continuum,

built from Turbospectrum syntheses made on demand around the star.

How it works
------------
* Atmospheres: the 2^5 APOGEE-MARCS nodes (Teff, logg, [M/H], [a/M], [C/M])
  of the grid cell containing the current estimate. Each node is synthesised
  with its own [M/H], [a/M], [C/M] (exactly like a grid point), and at a
  short ladder in [N/M] and in log10(vmicro) (synthesis-only parameters).
* Spectra are interpolated, never atmospheres: multilinear across the cell,
  quadratic (Lagrange, 3 points) in [N/M] and log vmicro.
* Broadening (vmac / vsini) is applied analytically, so it costs nothing.
* The instrument (R, R(lambda), any LSF array, incl. the APOGEE LSF) is
  applied exactly at the observed pixels, at the current RV.
* The continuum is a polynomial per chunk of the spectrum, multiplied into the
  model and solved at every step, so data do NOT need to be normalised first
  (data and model are treated identically).
* If the best fit sits on a face of the cell, the cell (or ladder) moves one
  step and the fit is repeated. Every synthesis is cached on disk and reused
  by later iterations, later stars and other instruments.

    import turbogrid, starfit
    sf = starfit.StarFitter(
        linelists=[...], inst=dict(R=22500.),          # or any LSF (broaden.py)
        nodes=dict(teff=(15, 3000., 250.), logg=(12, -0.5, 0.5), mh=(14, -2.5, 0.25),
                   am=(7, -1.0, 0.25), cm=(11, -1.5, 0.25)),
        atmod_lookup=turbogrid.make_atmod_lookup('/data/marcs_apogee', flat=True),
        cache_dir='/scratch/me/turbofit_cache', n_jobs=16)
    res = sf.fit(wave, flux, err, guess=dict(teff=4600, logg=2.2, mh=-0.3))
    res.params, res.errors, res.star, res.star_err

Use the result with linefit (res.star / res.star_err) for abundances.
"""

import os
import json
import hashlib
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import least_squares

import synthesis
import broaden
import continuum
import rv as rvmod

C_KMS = 299792.458
ATM = ('teff', 'logg', 'mh', 'am', 'cm')
SCALE = dict(teff=50., logg=0.1, mh=0.05, am=0.05, cm=0.05, nm=0.1, lgvmicro=0.05,
             vmac=1.0, vsini=1.0)


def _axis_values(v):
    if isinstance(v, tuple) and len(v) == 3:
        n, s, d = v
        return s + d * np.arange(int(n))
    return np.sort(np.atleast_1d(np.asarray(v, float)))


def _lagrange3(x, nodes):
    """Quadratic Lagrange weights for 3 nodes (or linear for 2, 1 for 1)."""
    n = len(nodes)
    if n == 1:
        return np.array([1.0])
    if n == 2:
        f = (x - nodes[0]) / (nodes[1] - nodes[0])
        return np.array([1 - f, f])
    w = np.ones(3)
    for i in range(3):
        for j in range(3):
            if i != j:
                w[i] *= (x - nodes[j]) / (nodes[i] - nodes[j])
    return w


def _synth_worker(a):
    """Synthesise one node spectrum (normalised) and save it to the cache."""
    path = a.pop('_cache')
    if os.path.exists(path):
        return path
    f, c, w = synthesis.synthesize(**a)
    if np.isscalar(f):
        raise RuntimeError(f'Turbospectrum failed for {a}')
    tmp = path + '.tmp.npz'
    np.savez(tmp, wave=w, flux=(f / c).astype(np.float32))
    os.replace(tmp, path)
    return path


class StarFitResult:
    def __init__(self, **kw):
        self.__dict__.update(kw)

    @property
    def star(self):
        """dict for linefit / synthesis: teff, logg, mh, am, cm, nm, vmicro, rv, vmac, vsini"""
        p = dict(self.params)
        out = {k: p[k] for k in ('teff', 'logg', 'mh', 'am', 'cm', 'nm') if k in p}
        out['vmicro'] = p['vmicro']
        out['rv'] = self.rv
        for k in ('vmac', 'vsini'):
            if k in p:
                out[k] = p[k]
        return out

    @property
    def star_err(self):
        e = {k: v for k, v in self.errors.items() if k in ('teff', 'logg', 'mh', 'am', 'cm', 'nm')}
        if 'vmicro' in self.errors:
            e['vmicro'] = self.errors['vmicro']
        return e

    def summary(self):
        lines = [f'chi2_red = {self.chi2r:.3f}   RV = {self.rv:.3f} +- {self.rv_err:.3f} km/s']
        for k, v in self.params.items():
            e = self.errors.get(k, np.nan)
            lines.append(f'{k:8s} = {v:10.4f} +- {e:.4f}' + ('  (fixed)' if k in self.fixed else ''))
        if self.flags:
            lines.append('flags: ' + ', '.join(self.flags))
        return '\n'.join(lines)


class StarFitter:
    """
    Parameters
    ----------
    linelists : list of Turbospectrum line lists.
    inst : instrument for broaden.lsf_matrix: R=..., R=(wave, R), lsf=callable
        + lsf_halfwidth, or lsf=array (rows = the observed pixels you fit) +
        lsf_offsets. For an LSF array, pass the SAME observed wavelength
        array to fit().
    nodes : atmosphere-node axes (teff, logg, mh, am, cm), as in turbogrid.
    atmod_lookup : callable(teff, logg, mh, am, cm) -> atmosphere file.
    atmos_type : atmosphere file type (default 'marcs').
    dw : synthesis step in A (default: ~1 km/s, lambda_blue / 300000).
    nm_step, lgvmicro_step : spacing of the [N/M] and log10(vmicro) ladders
        (3 points each, re-centred automatically). Default 0.5 dex and 0.15 dex.
    broadening : 'vmac' (default), 'vsini', or 'both'.
    cont_deg : continuum polynomial degree per chunk (default 3); chunks are
        found automatically at gaps, or pass cont_chunks=[(w0, w1), ...]
        (observed frame).
    c12c13 : 12C/13C used in all syntheses (default: synthesize() default).
    cache_dir : where node spectra are stored (reused across runs). Default
        ~/.turbofit_cache.
    n_jobs : parallel Turbospectrum runs.
    """

    def __init__(self, linelists, inst, nodes, atmod_lookup, atmos_type='marcs', dw=None,
                 nm_step=0.5, lgvmicro_step=0.15, broadening='vmac', cont_deg=3,
                 cont_chunks=None, c12c13=None, cache_dir=None, n_jobs=1, verbose=True):
        self.linelists = list(linelists)
        self.inst = dict(inst)
        self.nodes = {ax: _axis_values(nodes[ax]) for ax in ATM if ax in nodes}
        self.lookup = atmod_lookup
        self.atmos_type = atmos_type
        self.dw = dw
        self.nm_step = nm_step
        self.lgv_step = lgvmicro_step
        if broadening not in ('vmac', 'vsini', 'both'):
            raise ValueError("broadening must be 'vmac', 'vsini' or 'both'")
        self.broadening = broadening
        self.cont_deg = cont_deg
        self.cont_chunks = cont_chunks
        self.c12c13 = c12c13
        self.cache_dir = os.path.expanduser(cache_dir or '~/.turbofit_cache')
        os.makedirs(self.cache_dir, exist_ok=True)
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.n_synth = 0

    # ------------------------------------------------------------------
    def _log(self, *a):
        if self.verbose:
            print('starfit:', *a, flush=True)

    def _cell(self, p, ax):
        vals = self.nodes.get(ax, np.array([p[ax]]))
        if len(vals) == 1:
            return (vals[0], vals[0])
        x = float(np.clip(p[ax], vals[0], vals[-1]))
        i = int(np.clip(np.searchsorted(vals, x, side='right') - 1, 0, len(vals) - 2))
        return (vals[i], vals[i + 1])

    def _key(self, node, nm, vmicro):
        d = dict(node=[round(float(node[a]), 4) for a in ATM], nm=round(float(nm), 4),
                 vmicro=round(float(vmicro), 4), wr=[round(x, 3) for x in self.wrange],
                 segs=[[round(a, 3), round(b, 3)] for a, b in (self.segments or [])],
                 dw=self.dw_used, ll=self.linelists, c13=self.c12c13, at=self.atmos_type)
        h = hashlib.sha1(json.dumps(d, sort_keys=True).encode()).hexdigest()[:20]
        return os.path.join(self.cache_dir, f'tf_{h}.npz')

    def _ensure(self, specs):
        """specs: list of (node dict, nm, vmicro). Synthesise missing ones."""
        jobs = []
        for node, nm, vm in specs:
            path = self._key(node, nm, vm)
            if os.path.exists(path):
                continue
            atm = self.lookup(node['teff'], node['logg'], node['mh'], node['am'], node['cm'])
            if not atm or not os.path.exists(atm):
                raise FileNotFoundError(f'no atmosphere for node {node}; fill grid holes '
                                        '(atmos.fill_holes) or restrict the fit range')
            jobs.append(dict(teff=node['teff'], logg=node['logg'], mh=node['mh'], am=node['am'],
                             cm=node['cm'], nm=nm, vmicro=vm, wrange=list(self.wrange),
                             dw=self.dw_used, atmod=atm, atmos_type=self.atmos_type,
                             linelists=self.linelists, segments=self.segments,
                             c12c13=self.c12c13, _cache=path))
        if not jobs:
            return
        self._log(f'{len(jobs)} new syntheses ({len(specs) - len(jobs)} cached)')
        if self.n_jobs > 1 and len(jobs) > 1:
            with ProcessPoolExecutor(max_workers=self.n_jobs) as ex:
                list(ex.map(_synth_worker, jobs))
        else:
            for j in jobs:
                _synth_worker(j)
        self.n_synth += len(jobs)

    def _load(self, node, nm, vm):
        d = np.load(self._key(node, nm, vm))
        return d['wave'], d['flux'].astype(float)

    # ------------------------------------------------------------------
    def _build(self, p):
        """Cell + ladders around p; returns basis description and stacked spectra."""
        cells = {ax: self._cell(p, ax) for ax in ATM}
        nmv = np.round(p['nm'] + self.nm_step * np.array([-1, 0, 1]), 4)
        lgv = np.round(p['lgvmicro'] + self.lgv_step * np.array([-1, 0, 1]), 4)
        corners = []
        for v in range(2 ** len(ATM)):
            node = {ax: cells[ax][(v >> k) & 1] for k, ax in enumerate(ATM)}
            key = tuple(node[a] for a in ATM)
            if key not in [tuple(c[a] for a in ATM) for c in corners]:
                corners.append(node)
        specs = [(c, n, 10 ** g) for c in corners for n in nmv for g in lgv]
        self._ensure(specs)
        wave, S = None, []
        for c, n, vm in specs:
            w, f = self._load(c, n, vm)
            if wave is None:
                wave = w
            elif len(w) != len(wave) or not np.allclose(w, wave):
                f = np.interp(wave, w, f)
            S.append(f)
        return dict(cells=cells, corners=corners, nmv=nmv, lgv=lgv, wave=wave,
                    S=np.array(S, dtype=np.float32))

    def _weights(self, b, p):
        wc = []
        for c in b['corners']:
            wt = 1.0
            for ax in ATM:
                lo, hi = b['cells'][ax]
                if hi == lo:
                    continue
                f = (np.clip(p[ax], lo, hi) - lo) / (hi - lo)
                wt *= f if c[ax] == hi else (1 - f)
            wc.append(wt)
        wn = _lagrange3(p['nm'], b['nmv'])
        wg = _lagrange3(p['lgvmicro'], b['lgv'])
        return (np.array(wc)[:, None, None] * wn[None, :, None] * wg[None, None, :]).ravel()

    def _rv_scan(self, b, p, wo, fo, eo, inst, v0, model_fn, nstep=7):
        """chi^2(RV) of the full model around v0 (step ~ 1/4 resolution element
        or pixel); parabola through the points near the minimum."""
        pix_kms = np.median(np.diff(wo) / wo[1:]) * C_KMS
        step = max(pix_kms / 3.0, 0.05)
        vv = v0 + step * (np.arange(nstep) - nstep // 2)
        chi = []
        for v in vv:
            M = broaden.lsf_matrix(b['wave'], wo, rv=v, **inst)
            mod, _ = model_fn(b, M, {})
            chi.append(np.sum(((fo - mod) / eo) ** 2))
        chi = np.array(chi)
        i = int(np.argmin(chi))
        if i in (0, nstep - 1):
            # minimum outside the scan: move there and let the next iteration refine
            return float(vv[i]), np.nan
        c2, c1, c0 = np.polyfit(vv - vv[i], chi, 2)
        if c2 <= 0:
            return float(vv[i]), np.nan
        vbest = vv[i] - c1 / (2 * c2)
        dof = max(len(fo) - 1, 1)
        red = max(1.0, (c0 - c1 ** 2 / (4 * c2)) / dof)
        return float(vbest), float(np.sqrt(red / c2))

    # ------------------------------------------------------------------
    def fit(self, wave, flux, err, guess=None, fixed=None, fit_rv=True, max_iter=6,
            mc=0, seed=None, wrange=None, segments=None, mask=None):
        """
        Fit one spectrum.

        wave, flux, err : observed frame; flux need NOT be normalised.
        guess : starting values (any of teff, logg, mh, am, cm, nm, vmicro,
            vmac, vsini, rv). Missing ones: grid centre / 0 / 1.5 km/s /
            5 km/s; rv from a cross-correlation if not given.
        fixed : dict of parameters to hold fixed (e.g. dict(logg=2.45)).
        fit_rv : refine the RV against the best model between iterations.
        mc : Monte Carlo noise realisations for the errors (0 = covariance only).
        wrange / segments : rest-frame range or windows to synthesise. Default:
            one segment per contiguous chunk of the data (gaps skipped), with a
            margin for the RV. Pass segments= to skip regions you do not fit.
        mask : boolean array, True = use pixel (tellurics, bad regions off).
        """
        wave = np.asarray(wave, float); flux = np.asarray(flux, float); err = np.asarray(err, float)
        good = np.isfinite(flux) & np.isfinite(err) & (err > 0)
        if mask is not None:
            good &= np.asarray(mask, bool)
        fixed = dict(fixed or {})
        g = dict(guess or {})
        p = {}
        for ax in ATM:
            vals = self.nodes.get(ax, np.array([0.0]))
            p[ax] = float(g.get(ax, fixed.get(ax, 0.5 * (vals[0] + vals[-1]) if ax in ('teff', 'logg')
                                               else 0.0)))
        p['nm'] = float(g.get('nm', fixed.get('nm', 0.0)))
        p['lgvmicro'] = float(np.log10(g.get('vmicro', fixed.get('vmicro', 1.5))))
        if self.broadening in ('vmac', 'both'):
            p['vmac'] = float(g.get('vmac', fixed.get('vmac', 5.0)))
        if self.broadening in ('vsini', 'both'):
            p['vsini'] = float(g.get('vsini', fixed.get('vsini', 5.0 if self.broadening == 'vsini' else 2.0)))
        if 'vmicro' in fixed:
            fixed['lgvmicro'] = float(np.log10(fixed.pop('vmicro')))
        for k, v in fixed.items():
            p[k] = float(v)

        rvel = g.get('rv', fixed.get('rv'))
        f_ = broaden.doppler(1.0, rvel or 0.0)
        self.wrange = wrange or (wave[good][0] / broaden.doppler(1.0, 600.) - 3,
                                 wave[good][-1] / broaden.doppler(1.0, -600.) + 3)
        if segments is None and wrange is None:
            # synthesise only where there are data (chip gaps, separate gratings
            # or setups are skipped), with a +-600 km/s margin for the RV
            segments = [(a / broaden.doppler(1.0, 600.) - 3, b / broaden.doppler(1.0, -600.) + 3)
                        for a, b in continuum.find_chunks(wave[good])]
            merged = []
            for a, b in sorted(segments):
                if merged and a <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], b))
                else:
                    merged.append((a, b))
            segments = merged
        self.segments = segments
        # default synthesis step: ~1 km/s at the blue end, which samples the
        # intrinsic line profiles whatever the instrument
        self.dw_used = self.dw if self.dw is not None else round(self.wrange[0] / 3.0e5, 4)
        chunks = self.cont_chunks or continuum.find_chunks(wave[good])
        cidx = [np.nonzero(good & (wave >= a) & (wave <= b))[0] for a, b in chunks]
        cidx = [c for c in cidx if len(c) > self.cont_deg + 2]
        idx = np.concatenate(cidx)
        wo, fo, eo = wave[idx], flux[idx], err[idx]
        T = []
        for c in cidx:
            wc = wave[c]
            x = (wc - wc.mean()) / max(np.ptp(wc), 1e-9) * 2
            T.append((len(T), x))
        # design blocks for the continuum (one polynomial per chunk)
        nch = len(cidx)
        blocks = np.concatenate([np.full(len(c), k) for k, c in enumerate(cidx)])
        xpoly = np.concatenate([t[1] for t in T])

        inst = dict(self.inst)
        if inst.get('lsf') is not None and not callable(inst['lsf']):
            inst['lsf'] = np.asarray(inst['lsf'])[idx]

        free = [k for k in p if k not in fixed]
        self._log('free parameters:', ', '.join(free))

        # first RV from a cross-correlation (unbroadened model at the guess) if not given
        b = self._build(p)
        if rvel is None:
            s0 = np.dot(self._weights(b, p), b['S']).astype(float)
            _, cn, en = continuum.normalise(wo, fo, eo, chunks=chunks)
            rvel, _, _ = rvmod.ccf_rv(wo, cn, b['wave'], s0, err=en)
            self._log(f'initial RV from CCF: {rvel:.2f} km/s')
        rv_err = np.nan

        def model_fn(b, M, q, fdat=fo):
            pp = dict(p); pp.update(q)
            s = np.dot(self._weights(b, pp), b['S'])
            vm = pp.get('vmac', 0.0); vs = pp.get('vsini', 0.0)
            if vm > 0 or vs > 0:
                s = broaden.stellar(b['wave'], s, vmac=vm, vsini=vs)
            m = M @ s
            # continuum: polynomial per chunk, solved linearly
            out = np.empty_like(m)
            for k in range(nch):
                sel = blocks == k
                A = np.vstack([m[sel] * xpoly[sel] ** j for j in range(self.cont_deg + 1)]).T
                c, *_ = np.linalg.lstsq(A / eo[sel, None], fdat[sel] / eo[sel], rcond=None)
                out[sel] = A @ c
            return out, m

        flags = []
        for it in range(max_iter):
            b = self._build(p)
            M = broaden.lsf_matrix(b['wave'], wo, rv=rvel, **inst)
            lo, hi = [], []
            for k in free:
                if k in ATM:
                    c = b['cells'][k]
                    lo.append(c[0]); hi.append(c[1] if c[1] > c[0] else c[0] + 1e-6)
                elif k == 'nm':
                    lo.append(b['nmv'][0]); hi.append(b['nmv'][-1])
                elif k == 'lgvmicro':
                    lo.append(b['lgv'][0]); hi.append(b['lgv'][-1])
                else:
                    lo.append(0.0); hi.append(60.0)
            x0 = [float(np.clip(p[k], l + 1e-6 * (h - l), h - 1e-6 * (h - l))) for k, l, h in zip(free, lo, hi)]

            def resid(x, fdat=fo):
                mod, _ = model_fn(b, M, dict(zip(free, x)), fdat)
                return (fdat - mod) / eo

            r = least_squares(resid, x0, bounds=(lo, hi), x_scale=[SCALE[k] for k in free],
                              diff_step=1e-3)
            for k, v in zip(free, r.x):
                p[k] = float(v)
            chi2 = 2 * r.cost
            dof = max(len(fo) - len(free) - nch * (self.cont_deg + 1), 1)
            self._log(f'iteration {it + 1}: chi2_red = {chi2 / dof:.3f}  ' +
                      ' '.join(f'{k}={p[k]:.3f}' for k in free))

            # refine the RV: chi^2 of the full model (instrument applied exactly)
            # on a few velocities around the current one, parabola through them
            if fit_rv and 'rv' not in fixed:
                v_new, v_err = self._rv_scan(b, p, wo, fo, eo, inst, rvel, model_fn)
                moved = abs(v_new - rvel) if np.isfinite(v_new) else 0.0
                if np.isfinite(v_new):
                    rvel, rv_err = v_new, v_err
            else:
                moved = 0.0

            # move the cell / ladders if the solution sits on a face
            moved_cell = False
            for k in free:
                if k in ATM:
                    vals = self.nodes.get(k)
                    if vals is None or len(vals) < 2:
                        continue
                    c = b['cells'][k]
                    tol = 1e-3 * (c[1] - c[0])
                    if p[k] <= c[0] + tol and c[0] > vals[0]:
                        p[k] = c[0] - 0.5 * (c[1] - c[0]); moved_cell = True
                    elif p[k] >= c[1] - tol and c[1] < vals[-1]:
                        p[k] = c[1] + 0.5 * (c[1] - c[0]); moved_cell = True
                elif k in ('nm', 'lgvmicro'):
                    lad = b['nmv'] if k == 'nm' else b['lgv']
                    st = self.nm_step if k == 'nm' else self.lgv_step
                    if abs(p[k] - lad[1]) > 0.5 * st:
                        moved_cell = True        # re-centre the ladder on the new value
            if not moved_cell and moved < 0.1:
                break
            if it == max_iter - 1:
                flags.append('not_converged')

        # final evaluation, errors
        b = self._build(p)
        M = broaden.lsf_matrix(b['wave'], wo, rv=rvel, **inst)
        lo = [b['cells'][k][0] if k in ATM else (b['nmv'][0] if k == 'nm' else (b['lgv'][0] if k == 'lgvmicro' else 0.0)) for k in free]
        hi = [b['cells'][k][1] if k in ATM else (b['nmv'][-1] if k == 'nm' else (b['lgv'][-1] if k == 'lgvmicro' else 60.0)) for k in free]
        hi = [h if h > l else l + 1e-6 for l, h in zip(lo, hi)]

        def resid(x, fdat=fo):
            mod, _ = model_fn(b, M, dict(zip(free, x)), fdat)
            return (fdat - mod) / eo

        x0 = [float(np.clip(p[k], l, h)) for k, l, h in zip(free, lo, hi)]
        r = least_squares(resid, x0, bounds=(lo, hi), x_scale=[SCALE[k] for k in free], diff_step=1e-3)
        for k, v in zip(free, r.x):
            p[k] = float(v)
        chi2 = 2 * r.cost
        dof = max(len(fo) - len(free) - nch * (self.cont_deg + 1), 1)
        chi2r = chi2 / dof
        try:
            cov = np.linalg.inv(r.jac.T @ r.jac) * max(1.0, chi2r)
        except np.linalg.LinAlgError:
            cov = np.full((len(free), len(free)), np.nan)
        errs = {k: float(np.sqrt(cov[i, i])) for i, k in enumerate(free)}
        err_mc = {}
        if mc:
            rng = np.random.default_rng(seed)
            draws = []
            for _ in range(int(mc)):
                fk = fo + rng.normal(0, eo)
                rk = least_squares(lambda x: resid(x, fk), r.x, bounds=(lo, hi),
                                   x_scale=[SCALE[k] for k in free], diff_step=1e-3)
                draws.append(rk.x)
            draws = np.array(draws)
            err_mc = {k: float(np.std(draws[:, i], ddof=1)) for i, k in enumerate(free)}

        # report vmicro in km/s
        params = {k: p[k] for k in ATM}
        params['nm'] = p['nm']
        params['vmicro'] = 10 ** p['lgvmicro']
        for k in ('vmac', 'vsini'):
            if k in p:
                params[k] = p[k]
        errors = {}
        for k, e in errs.items():
            if k == 'lgvmicro':
                errors['vmicro'] = params['vmicro'] * np.log(10) * e
            else:
                errors[k] = e
        if err_mc:
            err_mc = {('vmicro' if k == 'lgvmicro' else k):
                      (params['vmicro'] * np.log(10) * v if k == 'lgvmicro' else v)
                      for k, v in err_mc.items()}
        for ax in ATM:
            vals = self.nodes.get(ax)
            if vals is not None and len(vals) > 1 and ax in free:
                if np.isclose(p[ax], vals[0]) or np.isclose(p[ax], vals[-1]):
                    flags.append(f'{ax}_at_grid_edge')
        corr = cov / np.sqrt(np.outer(np.diag(cov), np.diag(cov)))
        mod, m = model_fn(b, M, {})
        cont = mod / np.where(m != 0, m, 1)
        names = ['vmicro' if k == 'lgvmicro' else k for k in free]
        self._log(f'done: chi2_red={chi2r:.3f}, {self.n_synth} syntheses in total this session')
        return StarFitResult(params=params, errors=errors, errors_mc=err_mc, cov=cov, corr=corr,
                             free=names, fixed={('vmicro' if k == 'lgvmicro' else k):
                                                (10 ** v if k == 'lgvmicro' else v)
                                                for k, v in fixed.items()},
                             chi2r=chi2r, rv=float(rvel), rv_err=float(rv_err), flags=flags,
                             wave=wo, flux=fo, err=eo, model=mod, continuum=cont,
                             norm_flux=fo / cont, norm_model=m, n_synth=self.n_synth)


def plot_fit(res, outfile=None, chunk_width=None, colors=('#3A1450', '#D6458F')):
    """Normalised data and best model, one panel per chunk (or per
    chunk_width A), with residuals."""
    import matplotlib
    import matplotlib.pyplot as plt
    rc = {'font.family': 'serif', 'font.size': 14, 'axes.labelsize': 14,
          'axes.titlesize': 14, 'legend.fontsize': 10,
          'xtick.direction': 'in', 'ytick.direction': 'in',
          'xtick.top': True, 'ytick.right': True,
          'xtick.minor.visible': True, 'ytick.minor.visible': True}
    w = res.wave
    edges = continuum.find_chunks(w)
    if chunk_width:
        edges = [(a, min(a + chunk_width, b)) for (c0, b) in edges
                 for a in np.arange(c0, b, chunk_width)]
    with matplotlib.rc_context(rc):
        fig, axes = plt.subplots(2 * len(edges), 1, figsize=(12, 3.2 * len(edges)),
                                 gridspec_kw=dict(height_ratios=[3, 1] * len(edges)),
                                 squeeze=False)
        for k, (a, b_) in enumerate(edges):
            m = (w >= a) & (w <= b_)
            ax, axr = axes[2 * k, 0], axes[2 * k + 1, 0]
            ax.plot(w[m], res.norm_flux[m], '.', color=colors[0], ms=2, label='data')
            ax.plot(w[m], res.norm_model[m], '-', color=colors[1], lw=1, label='model')
            axr.plot(w[m], res.norm_flux[m] - res.norm_model[m], '.', color=colors[0], ms=2)
            axr.axhline(0, color=colors[1], lw=0.8)
            ax.set_ylabel('Normalised flux'); axr.set_ylabel('Resid.')
            ax.set_xlim(a, b_); axr.set_xlim(a, b_)
        axes[0, 0].legend(loc='lower left')
        axes[-1, 0].set_xlabel(r'Observed wavelength ($\AA$)')
        fig.tight_layout()
        if outfile:
            fig.savefig(outfile)
    return fig
