"""
apogee_gc_test.py
=================

Run turbofit on APOGEE DR17 globular-cluster members (from the GC
value-added catalogue) and compare with ASPCAP.

Steps
-----
1. select members of the chosen clusters from the GC VAC, matched to allStar;
2. read each star's apStar spectrum, mask bad pixels, build its LSF from
   the fibres of its visits (allVisit);
3. run turbofit.analyse() in local mode (one result file per star, so the
   run can be restarted or split over a job array);
4. collect everything into one table next to ASPCAP's uncalibrated values;
5. plots: parameter offsets, abundance offsets, cluster spreads against the
   reported errors, and the C-N / Mg-Al / 12C/13C patterns.

Usage
-----
    python apogee_gc_test.py list                      # members per cluster
    python apogee_gc_test.py download                  # their apStar/asStar files (login node)
    python apogee_gc_test.py run M13                   # all selected M13 stars
    python apogee_gc_test.py run M13 --index 4         # one star (job arrays)
    python apogee_gc_test.py collect                   # table + plots

Edit the CONFIG block first. Column names are checked on start-up; if one
is missing, the available columns are printed.
"""

import os
import sys
import json
import argparse
import numpy as np
from astropy.io import fits
from astropy.table import Table

import turbofit
import turbogrid
import apogee_lsf

# =====================================================================
# CONFIG
# =====================================================================

BASE = '/mnt/aridata1/users/arialipt'
D = BASE + '/apogee_dr17'                     # where get_apogee.sh puts the files
GC_VAC = D + '/GC_members_VAC-v1_1.fits'
ALLSTAR = D + '/allStar-dr17-synspec_rev1.fits'
ALLVISIT = D + '/allVisit-dr17-synspec_rev1.fits'
APSTAR = D + '/stars/{telescope}/{field}/{file}'   # filled by 'download'
SAS_URL = 'https://data.sdss.org/sas/dr17/apogee/spectro/redux/dr17/stars/{telescope}/{field}/{file}'
# {file} is allStar's FILE column: apStar-dr17-<id>.fits (APO) or asStar-dr17-<id>.fits (LCO)

# DR17 LSF and wavelength-solution IDs (cal/lsf, cal/wave). Each visit gets
# the most recent LSF taken on or before its MJD; the star's LSF is the
# average over its visits (as the apStar spectrum combines them).
LSF_DIR = os.environ.get('APOGEE_LSF_DIR', BASE + '/apogee_dr17/lsf')
LSF_IDS = {
    'apo25m': ('ap', ['03430016', '07510018', '11130063', '14600018', '18430026',
                      '22330043', '25560065', '28860055', '32900046']),
    'lco25m': ('as', ['22670019', '22940020', '26990075', '29380037', '32690044',
                      '36210017']),
}
WAVE_IDS = {
    'apo25m': ['02380000', '05680000', '09500000', '13140000', '16680000',
               '20380000', '24040000', '27980000', '31380000'],
    'lco25m': ['22670000', '24040000', '27980000', '31380000', '33880000'],
}
FIBRES = 'visits'          # 'visits' (fibres of the star's visits) or 'combo'

CLUSTERS = ['Liller1', 'NGC6752']   # as named in GC_NAME

SNR_MIN = 50.
MAX_PER_CLUSTER = None       # highest S/N first; None for all

LINELISTS = ['/users/arialipt/turbofit/aplinelists/turbospec.20180901t20.atoms',
             '/users/arialipt/turbofit/aplinelists/turbospec.20180901t20.Hlinedata',
             '/users/arialipt/turbofit/aplinelists/turbospec.20180901t20.molec']
MARCS_DIR = '/mnt/aridata1/users/arialipt/marcs_apogee'                         # APOGEE MARCS (C, alpha, M nodes)
LINELISTS_AIR = False      # True if the line lists are in air (apStar is vacuum): Fe I at 15207.53 = air, 15211.68 = vacuum
TEFF_NODES = [float(t) for t in list(range(2800, 4001, 100)) + list(range(4250, 8001, 250))]
NODES = dict(teff=TEFF_NODES, logg=(12, -0.5, 0.5), mh=(14, -2.5, 0.25),
             am=(7, -1.0, 0.25), cm=(11, -1.5, 0.25))
ELEMENTS = ['C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'S', 'K', 'Ca', 'Ti', 'V',
            'Mn', 'Fe', 'Co', 'Ni', 'Ce', 'C12C13']
BROADENING = 'vmac'
VERBOSE = True             # print the fit's progress (cells, chi2, parameters)
FILL_HOLES = True          # missing MARCS node -> interpolated (or nearest) model, see marcs_fill.py
FILL_DIR = BASE + '/marcs_apogee_filled'   # filled models are written here (grid left untouched)
FIX_VMICRO = False         # True: hold vmicro at ASPCAP's value (diagnostic run)
N_JOBS = 16
CACHE = BASE + '/turbofit_cache'
OUT = BASE + '/gc_test'

# APOGEE_PIXMASK bits treated as bad
BAD_BITS = [0, 1, 2, 3, 4, 5, 6, 9, 12]   # BADPIX CRPIX SATPIX UNFIXABLE BADDARK
                                          # BADFLAT BADERR PERSIST_HIGH SIG_SKYLINE

# ASPCAP columns used for the comparison (uncalibrated where it matters)
ASPCAP_PARAMS = dict(teff='TEFF_SPEC', logg='LOGG_SPEC', mh='M_H', am='ALPHA_M',
                     vmicro='VMICRO', vmac='VMACRO')

# =====================================================================
# tables and selection
# =====================================================================

def _need(tab, cols, name):
    miss = [c for c in cols if c not in tab.colnames]
    if miss:
        print(f'{name}: missing columns {miss}\navailable: {tab.colnames}')
        sys.exit(1)


def _str(col):
    return np.char.strip(np.asarray(col).astype(str))


def load_members():
    vac = Table.read(GC_VAC, hdu=1)
    _need(vac, ['APOGEE_ID', 'GC_NAME'], 'GC VAC')
    ok = np.isin(_str(vac['GC_NAME']), CLUSTERS)
    vac = vac[ok]
    vac['APOGEE_ID'] = _str(vac['APOGEE_ID'])
    vac['GC_NAME'] = _str(vac['GC_NAME'])

    star = Table.read(ALLSTAR, hdu=1)
    _need(star, ['APOGEE_ID', 'TELESCOPE', 'FIELD', 'FILE', 'SNR'] + list(ASPCAP_PARAMS.values()),
          'allStar')
    star['APOGEE_ID'] = _str(star['APOGEE_ID'])
    star = star[np.isin(star['APOGEE_ID'], vac['APOGEE_ID'])]
    star = star[np.asarray(star['SNR'], float) >= SNR_MIN]

    # one row per star: keep the highest-S/N entry (stars observed from both sites)
    star.sort('SNR', reverse=True)
    _, first = np.unique(star['APOGEE_ID'], return_index=True)
    star = star[np.sort(first)]

    gc = dict(zip(vac['APOGEE_ID'], vac['GC_NAME']))
    star['GC_NAME'] = [gc[i] for i in star['APOGEE_ID']]
    sel = []
    for c in CLUSTERS:
        t = star[star['GC_NAME'] == c]
        t.sort('SNR', reverse=True)
        sel.append(t if MAX_PER_CLUSTER is None else t[:MAX_PER_CLUSTER])
    from astropy.table import vstack
    return vstack(sel)


_VISITS = None


def star_visits(apogee_id, telescope):
    """(fibre, mjd) of every visit of this star from this telescope."""
    global _VISITS
    if _VISITS is None:                       # read allVisit once, three columns only
        with fits.open(ALLVISIT, memmap=True) as h:
            d = h[1].data
            _VISITS = Table(dict(APOGEE_ID=_str(d['APOGEE_ID']), TELESCOPE=_str(d['TELESCOPE']),
                                 FIBERID=np.asarray(d['FIBERID']),
                                 MJD=np.asarray(d['MJD'])))
    v = _VISITS
    ok = (v['APOGEE_ID'] == apogee_id) & (v['TELESCOPE'] == telescope)
    return [(int(f), int(m)) for f, m in zip(v['FIBERID'][ok], v['MJD'][ok])]

# =====================================================================
# spectra
# =====================================================================

def read_apstar(row):
    path = APSTAR.format(telescope=row['TELESCOPE'].strip(), field=row['FIELD'].strip(),
                         file=row['FILE'].strip())
    with fits.open(path) as h:
        flux = np.atleast_2d(h[1].data)[0].astype(float)    # row 0: combined spectrum
        err = np.atleast_2d(h[2].data)[0].astype(float)
        mask = np.atleast_2d(h[3].data)[0].astype(np.int64)
    wave = apogee_lsf.apStarWavegrid()
    if LINELISTS_AIR:                      # vacuum -> air (Morton 2000)
        s2 = (1e4 / wave) ** 2
        wave = wave / (1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2))
    bad = np.zeros(len(wave), bool)
    for b in BAD_BITS:
        bad |= (mask & (1 << b)) > 0
    bad |= ~np.isfinite(flux) | (flux <= 0) | ~np.isfinite(err) | (err <= 0) | (err > 1e9)
    flux[bad] = np.nan            # chip gaps become separate chunks; single pixels are skipped
    err[bad] = np.nan
    return wave, flux, err


def download(sel):
    """Fetch the apStar files of the selected stars (skips ones already there)."""
    import urllib.request
    for row in sel:
        kw = dict(telescope=row['TELESCOPE'].strip(), field=row['FIELD'].strip(),
                  file=row['FILE'].strip())
        path = APSTAR.format(**kw)
        if os.path.exists(path) and os.path.getsize(path) > 0:
            continue
        os.makedirs(os.path.dirname(path), exist_ok=True)
        try:
            urllib.request.urlretrieve(SAS_URL.format(**kw), path + '.part')
            os.replace(path + '.part', path)
            print('got', path)
        except Exception as e:
            print('FAILED', row['APOGEE_ID'], SAS_URL.format(**kw), repr(e))


def _id_mjd(i):
    return int(i[:4]) + 55562            # APOGEE IDs: first 4 digits = MJD - 55562


def _pick(ids, mjd):
    """Most recent ID on or before mjd (else the earliest)."""
    before = [i for i in ids if _id_mjd(i) <= mjd]
    return max(before, key=_id_mjd) if before else min(ids, key=_id_mjd)


_LSF_CACHE = {}


def instrument(row, wave):
    tel = row['TELESCOPE'].strip()
    prefix, lsf_ids = LSF_IDS[tel]
    visits = star_visits(row['APOGEE_ID'], tel)
    if not visits:
        raise ValueError(f"no visits for {row['APOGEE_ID']} in allVisit")
    # group visits by the (LSF, wave) calibration they get
    groups = {}
    for fib, mjd in visits:
        key = (_pick(lsf_ids, mjd), _pick(WAVE_IDS[tel], mjd))
        groups.setdefault(key, []).append(fib)
    lsf_sum, x, n = 0.0, None, 0
    for (lid, wid), fibs in groups.items():
        if (lid, wid) not in _LSF_CACHE:
            _LSF_CACHE[(lid, wid)] = apogee_lsf.ApogeeLSF(LSF_DIR, lsfid=lid, waveid=wid,
                                                          prefix=prefix)
        fib = 'combo' if FIBRES == 'combo' else [f for f in fibs if 1 <= f <= 300] or 'combo'
        lsf, x = _LSF_CACHE[(lid, wid)].lsf_array(wave, fiber=fib)
        lsf_sum = lsf_sum + lsf * len(fibs)          # weight by number of visits
        n += len(fibs)
    lsf = lsf_sum / n
    lsf /= lsf.sum(axis=1, keepdims=True)
    print(row['APOGEE_ID'], 'LSF:', {k[0]: len(v) for k, v in groups.items()})
    return dict(lsf=lsf, lsf_offsets=x)

# =====================================================================
# one star
# =====================================================================

def _f(x):
    try:
        x = float(x)
        return x if np.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def run_star(row, lookup):
    out = os.path.join(OUT, 'stars', f"{row['GC_NAME']}_{row['APOGEE_ID']}.json")
    if os.path.exists(out):
        print('done already:', out)
        return
    os.makedirs(os.path.dirname(out), exist_ok=True)
    wave, flux, err = read_apstar(row)
    inst = instrument(row, apogee_lsf.apStarWavegrid())   # LSF on the detector (vacuum) grid

    guess = {k: _f(row[c]) for k, c in ASPCAP_PARAMS.items() if k in ('teff', 'logg', 'mh', 'am')}
    guess = {k: v for k, v in guess.items() if v is not None}
    fixed = None
    if FIX_VMICRO and _f(row['VMICRO']):
        fixed = dict(vmicro=_f(row['VMICRO']))

    res = turbofit.analyse(wave, flux, err, linelists=LINELISTS, inst=inst,
                           elements=ELEMENTS, mode='local', nodes=NODES,
                           atmod_lookup=lookup, guess=guess, fixed=fixed,
                           broadening=BROADENING, n_jobs=N_JOBS, cache_dir=CACHE,
                           workdir=os.path.join(OUT, 'work', row['APOGEE_ID']),
                           verbose=VERBOSE)
    p = res['params']
    rec = dict(apogee_id=row['APOGEE_ID'], gc=row['GC_NAME'], telescope=row['TELESCOPE'].strip(),
               snr=_f(row['SNR']), chi2r=_f(p.chi2r), flags=list(p.flags),
               star={k: _f(v) for k, v in res['star'].items()},
               star_err={k: _f(v) for k, v in res['star_err'].items()},
               lines={el: {k: (_f(v) if not isinstance(v, str) else v) for k, v in s.items()}
                      for el, s in res['lines'].summary.items()})
    res['lines'].write(out.replace('.json', '_lines.csv'))
    with open(out, 'w') as f:
        json.dump(rec, f, indent=1)
    print('wrote', out)

# =====================================================================
# collect + compare
# =====================================================================

def collect():
    import glob
    files = sorted(glob.glob(os.path.join(OUT, 'stars', '*.json')))
    if not files:
        sys.exit('no results yet')
    recs = [json.load(open(f)) for f in files]
    star = Table.read(ALLSTAR, hdu=1)
    star['APOGEE_ID'] = _str(star['APOGEE_ID'])
    star['TELESCOPE'] = _str(star['TELESCOPE'])
    idx = {(a, t): i for i, (a, t) in enumerate(zip(star['APOGEE_ID'], star['TELESCOPE']))}

    rows = []
    for r in recs:
        a = star[idx[(r['apogee_id'], r['telescope'])]]
        d = dict(apogee_id=r['apogee_id'], gc=r['gc'], snr=r['snr'], chi2r=r['chi2r'],
                 flags=','.join(r['flags']))
        for k in ('teff', 'logg', 'mh', 'am', 'cm', 'nm', 'vmicro', 'vmac', 'rv'):
            d[k] = r['star'].get(k)
            d[k + '_err'] = r['star_err'].get(k)
        for k, c in ASPCAP_PARAMS.items():
            d['aspcap_' + k] = _f(a[c])
        fe = r['lines'].get('Fe', {}).get('xh')
        for el in ELEMENTS:
            s = r['lines'].get(el, {})
            if el == 'C12C13':
                d['c12c13'] = s.get('c12c13')
                d['c12c13_lo'], d['c12c13_hi'] = s.get('c12c13_lo'), s.get('c12c13_hi')
                continue
            xh = s.get('xh')
            d[el + '_H'] = xh
            d[el + '_Fe'] = None if el == 'Fe' or xh is None or fe is None else xh - fe
            d[el + '_err'] = s.get('err_total')
            d[el + '_n'] = s.get('n')
            col = 'FE_H_SPEC' if el == 'Fe' else el.upper() + '_FE_SPEC'
            d['aspcap_' + (el + '_H' if el == 'Fe' else el + '_Fe')] = \
                _f(a[col]) if col in star.colnames else None
        rows.append(d)
    t = Table(rows=[[r.get(k) for k in rows[0]] for r in rows], names=list(rows[0]),
              dtype=[object] * len(rows[0]))
    for c in t.colnames:                      # numeric columns -> float (None -> nan)
        if c not in ('apogee_id', 'gc', 'flags'):
            t[c] = np.array([np.nan if v is None else v for v in t[c]], float)
        else:
            t[c] = np.array(t[c], str)
    t.write(os.path.join(OUT, 'gc_turbofit_vs_aspcap.fits'), overwrite=True)
    t.write(os.path.join(OUT, 'gc_turbofit_vs_aspcap.csv'), overwrite=True)
    spreads(t)
    plots(t)
    return t


def spreads(t):
    """Per cluster and element: observed scatter vs median reported error."""
    lines = ['gc,element,n,median,std,median_err,ratio']
    for gc in CLUSTERS:
        s = t[t['gc'] == gc]
        if len(s) < 3:
            continue
        for el in [e for e in ELEMENTS if e != 'C12C13']:
            col = 'Fe_H' if el == 'Fe' else el + '_Fe'
            ok = np.isfinite(s[col]) & np.isfinite(s[el + '_err'])
            if ok.sum() < 3:
                continue
            x, e = s[col][ok], s[el + '_err'][ok]
            sd, me = np.std(x, ddof=1), np.median(e)
            lines.append(f'{gc},{el},{ok.sum()},{np.median(x):.3f},{sd:.3f},{me:.3f},{sd / me:.2f}')
    path = os.path.join(OUT, 'gc_spreads.csv')
    open(path, 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))

# =====================================================================
# plots
# =====================================================================

def _style():
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'serif', 'font.size': 14, 'axes.labelsize': 14,
                         'axes.titlesize': 14, 'legend.fontsize': 10,
                         'xtick.direction': 'in', 'ytick.direction': 'in',
                         'xtick.top': True, 'ytick.right': True,
                         'xtick.minor.visible': True, 'ytick.minor.visible': True})
    return plt


def plots(t):
    plt = _style()
    gcs = [g for g in CLUSTERS if np.any(t['gc'] == g)]
    cols = dict(zip(gcs, plt.cm.plasma(np.linspace(0.05, 0.85, len(gcs)))))

    # 1. parameters: turbofit - ASPCAP against ASPCAP Teff
    keys = [k for k in ('teff', 'logg', 'mh', 'am', 'vmicro', 'vmac') if 'aspcap_' + k in t.colnames]
    lab = dict(teff=r'$T_{\rm eff}$', logg=r'$\log g$', mh='[M/H]', am=r'[$\alpha$/M]',
               vmicro=r'$v_{\rm micro}$', vmac=r'$v_{\rm mac}$')
    fig, ax = plt.subplots(len(keys), 1, figsize=(7, 2.2 * len(keys)), sharex=True)
    for a, k in zip(np.atleast_1d(ax), keys):
        for g in gcs:
            s = t[t['gc'] == g]
            a.errorbar(s['aspcap_teff'], s[k] - s['aspcap_' + k], yerr=s[k + '_err'],
                       fmt='o', ms=4, color=cols[g], label=g, alpha=0.8)
        a.axhline(0, color='k', lw=0.8)
        a.set_ylabel(r'$\Delta$' + lab[k])
    np.atleast_1d(ax)[0].legend(ncol=3)
    np.atleast_1d(ax)[-1].set_xlabel(r'ASPCAP $T_{\rm eff}$ (K)')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'gc_params.pdf'))

    # 2. abundances: median and scatter of turbofit - ASPCAP per element
    els = [e for e in ELEMENTS if e != 'C12C13']
    fig, ax = plt.subplots(figsize=(10, 4.5))
    for j, g in enumerate(gcs):
        s = t[t['gc'] == g]
        for i, el in enumerate(els):
            col = 'Fe_H' if el == 'Fe' else el + '_Fe'
            if 'aspcap_' + col not in t.colnames:
                continue
            d = s[col] - s['aspcap_' + col]
            d = d[np.isfinite(d)]
            if len(d):
                x = i + (j - (len(gcs) - 1) / 2) * 0.12
                ax.errorbar(x, np.median(d), yerr=np.std(d) if len(d) > 1 else 0,
                            fmt='o', ms=5, color=cols[g], label=g if i == 0 else None)
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xticks(range(len(els)))
    ax.set_xticklabels(els)
    ax.tick_params(axis='x', which='minor', bottom=False, top=False)
    ax.set_ylabel(r'$\Delta$[X/Fe] (turbofit $-$ ASPCAP)')
    ax.legend(ncol=3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'gc_abundances.pdf'))

    # 3. multiple-population patterns
    pairs = [('C_Fe', 'N_Fe', '[C/Fe]', '[N/Fe]'), ('Mg_Fe', 'Al_Fe', '[Mg/Fe]', '[Al/Fe]'),
             ('logg', 'c12c13', r'$\log g$', r'$^{12}$C/$^{13}$C')]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    for a, (x, y, xl, yl) in zip(ax, pairs):
        for g in gcs:
            s = t[t['gc'] == g]
            xe = s[x.split('_')[0] + '_err']
            if y == 'c12c13':
                ye = [s['c12c13'] - s['c12c13_lo'], s['c12c13_hi'] - s['c12c13']]
            else:
                ye = s[y.split('_')[0] + '_err']
            a.errorbar(s[x], s[y], xerr=xe, yerr=ye, fmt='o', ms=4, color=cols[g],
                       label=g, alpha=0.8)
        a.set_xlabel(xl)
        a.set_ylabel(yl)
    ax[2].invert_xaxis()
    ax[0].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'gc_patterns.pdf'))
    print('plots in', OUT)

# =====================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('what', choices=['list', 'download', 'run', 'collect'])
    ap.add_argument('cluster', nargs='?')
    ap.add_argument('--index', type=int)
    a = ap.parse_args()
    if a.what == 'collect':
        collect()
        return
    sel = load_members()
    if a.what == 'list':
        for g in CLUSTERS:
            s = sel[sel['GC_NAME'] == g]
            print(f'{g:8s} {len(s):3d} stars, S/N {np.min(s["SNR"]) if len(s) else 0:.0f}+')
        return
    if a.cluster:
        sel = sel[sel['GC_NAME'] == a.cluster]
    if a.index is not None:
        sel = sel[a.index:a.index + 1]
    if a.what == 'download':
        download(sel)
        return
    lookup = turbogrid.make_atmod_lookup(MARCS_DIR, flat=True)
    if FILL_HOLES:
        import marcs_fill
        filler = marcs_fill.Filler(MARCS_DIR, FILL_DIR, NODES['teff'])
        lookup = filler.lookup(lookup, turbogrid.make_atmod_lookup(FILL_DIR, flat=True))
    for row in sel:
        try:
            run_star(row, lookup)
        except Exception as e:                  # keep going; note the failure
            import traceback
            traceback.print_exc()
            print('FAILED', row['APOGEE_ID'], repr(e))


if __name__ == '__main__':
    main()