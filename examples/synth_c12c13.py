"""
synth_c12c13.py
===============

Two synthetic spectra that differ only in 12C/13C, on the APOGEE MARCS grid
(C and alpha nodes), over 0.9-3.2 micron, convolved to R = 1400, and a plot.

    python synth_c12c13.py            # synthesise (parallel chunks) + plot
    python synth_c12c13.py plot       # re-plot from the saved spectra

Run it as a job (the full range is long): see N_JOBS.
"""

import os
import sys
import glob
import numpy as np
from concurrent.futures import ProcessPoolExecutor

import synthesis
import turbogrid
import broaden

# =====================================================================
# CONFIG
# =====================================================================

STAR = dict(teff=3800., logg=0.5, mh=-0.5, am=0.25, cm=0.25, nm=0.0, vmicro=2.0)
C12C13 = [15., 20., 50., 89.]           # one spectrum per ratio; differences are
                                        # plotted relative to the last one

MARCS_DIR = '/mnt/aridata1/users/arialipt/marcs_apogee'
LINELIST_DIR = '/users/arialipt/turbofit/lima_linelists'
LINELISTS = sorted(f for f in glob.glob(os.path.join(LINELIST_DIR, '*'))
                   if os.path.isfile(f) and not os.path.basename(f).startswith('.'))

WRANGE = (9000., 32000.)                # Angstrom (in the line lists' frame)
R_SYNTH = 300000.                       # synthesis sampling, lambda/dlambda
CHUNK = 500.                            # Angstrom per Turbospectrum run
R_OUT = 1400.                           # final resolving power
PIX_PER_RES = 3.                        # output pixels per resolution element
N_JOBS = 32                             # chunks synthesised in parallel

OUT = '/mnt/aridata1/users/arialipt/synth_c12c13'

# =====================================================================


def _chunks():
    edges = np.arange(WRANGE[0], WRANGE[1], CHUNK)
    return [(a, min(a + CHUNK, WRANGE[1])) for a in edges]


def _synth_chunk(args):
    (a, b), ratio, atm = args
    pad = 5.0                                           # overlap, trimmed afterwards
    dw = round(a / R_SYNTH, 4)
    f, c, w = synthesis.synthesize(STAR['teff'], STAR['logg'], mh=STAR['mh'], am=STAR['am'],
                                   cm=STAR['cm'], nm=STAR['nm'], vmicro=STAR['vmicro'],
                                   wrange=[a - pad, b + pad], dw=dw, atmod=atm,
                                   atmos_type='marcs', linelists=LINELISTS, c12c13=ratio)
    if np.isscalar(f):
        raise RuntimeError(f'Turbospectrum failed for {a:.0f}-{b:.0f} A, 12C/13C={ratio}')
    keep = (w >= a) & (w < b)
    return w[keep], f[keep], c[keep]


def synthesise():
    os.makedirs(OUT, exist_ok=True)
    if not LINELISTS:
        sys.exit(f'no line lists found in {LINELIST_DIR}')
    print('line lists:', *LINELISTS, sep='\n   ')
    lookup = turbogrid.make_atmod_lookup(MARCS_DIR, flat=True)
    atm = lookup(STAR['teff'], STAR['logg'], STAR['mh'], STAR['am'], STAR['cm'])
    if atm is None:
        sys.exit('no APOGEE MARCS model for this node (a grid hole)')
    print('atmosphere:', atm)

    # one common log-lambda grid for the R = 1400 output
    dlog = np.log10(1 + 1 / (R_OUT * PIX_PER_RES))
    w_out = 10 ** np.arange(np.log10(WRANGE[0] + 20), np.log10(WRANGE[1] - 20), dlog)
    for ratio in C12C13:
        path = os.path.join(OUT, f'synth_c12c13_{ratio:g}.npz')
        if os.path.exists(path):
            print('have', path)
            continue
        jobs = [(ch, ratio, atm) for ch in _chunks()]
        with ProcessPoolExecutor(max_workers=N_JOBS) as ex:
            parts = list(ex.map(_synth_chunk, jobs))
        w = np.concatenate([p[0] for p in parts])
        f = np.concatenate([p[1] for p in parts])
        c = np.concatenate([p[2] for p in parts])
        # uniform log-lambda grid at the synthesis sampling, then R_OUT
        wl = 10 ** np.arange(np.log10(w[0]), np.log10(w[-1]), np.log10(1 + 1 / R_SYNTH))
        fl, cl = np.interp(wl, w, f), np.interp(wl, w, c)
        _, fc = broaden.convolve(wl, np.vstack([fl, cl]), R=R_OUT, wave_out=w_out)
        np.savez(path, wave=w_out, flux=fc[0], cont=fc[1], wave_hr=w, norm_hr=(f / c).astype(np.float32),
                 c12c13=ratio, **STAR)
        print(f'12C/13C = {ratio:g}: {len(w)} synthesis points -> {len(w_out)} at R={R_OUT:g}; saved {path}')


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'serif', 'font.size': 14, 'axes.labelsize': 14,
                         'axes.titlesize': 14, 'legend.fontsize': 10,
                         'xtick.direction': 'in', 'ytick.direction': 'in',
                         'xtick.top': True, 'ytick.right': True,
                         'xtick.minor.visible': True, 'ytick.minor.visible': True})
    sp = [np.load(os.path.join(OUT, f'synth_c12c13_{r:g}.npz')) for r in C12C13]
    w = sp[0]['wave'] / 1e4                                     # micron
    n = [s['flux'] / s['cont'] for s in sp]
    title = (rf"$T_{{\rm eff}}$={STAR['teff']:.0f} K, $\log g$={STAR['logg']}, "
             rf"[M/H]={STAR['mh']:+.2f}, [$\alpha$/M]={STAR['am']:+.2f}, [C/M]={STAR['cm']:+.2f}, "
             rf"R={R_OUT:g}")

    fig, ax = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                           gridspec_kw=dict(height_ratios=[3, 1], hspace=0))
    for s, f in zip(sp, n):
        ax[0].plot(w, f, lw=0.8, label=rf"$^{{12}}$C/$^{{13}}$C = {float(s['c12c13']):g}")
    for s, f in zip(sp[:-1], n[:-1]):
        ax[1].plot(w, f - n[-1], lw=0.8)
    ax[0].set_ylabel('flux / continuum')
    ax[0].set_title(title)
    ax[0].legend()
    ax[1].set_ylabel(rf"$-$ ({float(sp[-1]['c12c13']):g})")
    ax[1].set_xlabel(r'wavelength ($\mu$m)')
    fig.savefig(os.path.join(OUT, 'synth_c12c13_full.pdf'), bbox_inches='tight')

    # zooms: (wmin, wmax in micron, name, [(bandhead, label), ...])
    zooms = [(1.03, 1.12, 'CN', []),
             (1.54, 1.66, 'CO_2nd', [(1.5582, r'$^{12}$CO 3-0'), (1.5780, r'$^{12}$CO 4-1'),
                                     (1.5982, r'$^{12}$CO 5-2'), (1.6189, r'$^{12}$CO 6-3'),
                                     (1.6397, r'$^{12}$CO 7-4')]),
             (2.28, 2.42, 'CO', [(2.2935, r'$^{12}$CO 2-0'), (2.3227, r'$^{12}$CO 3-1'),
                                 (2.3448, r'$^{13}$CO 2-0'), (2.3739, r'$^{13}$CO 3-1')])]
    for w0, w1, name, heads in zooms:
        fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True,
                               gridspec_kw=dict(height_ratios=[3, 1], hspace=0))
        z = (w > w0) & (w < w1)
        for s, f in zip(sp, n):
            ax[0].plot(w[z], f[z], lw=1, label=rf"$^{{12}}$C/$^{{13}}$C = {float(s['c12c13']):g}")
        for f in n[:-1]:
            ax[1].plot(w[z], (f - n[-1])[z], lw=1)
        for x, lab in heads:
            ax[0].axvline(x, color='gray', ls=':', lw=1)
            ax[0].text(x, 1.02, lab, rotation=90, va='bottom', ha='right', fontsize=10,
                       transform=ax[0].get_xaxis_transform())
        ax[0].set_ylabel('flux / continuum')
        ax[0].legend(loc='lower left')
        ax[1].set_ylabel(rf"$-$ ({float(sp[-1]['c12c13']):g})")
        ax[1].set_xlabel(r'wavelength ($\mu$m)')
        fig.savefig(os.path.join(OUT, f'synth_c12c13_{name}.pdf'), bbox_inches='tight')
    print('plots in', OUT)


if __name__ == '__main__':
    if not (len(sys.argv) > 1 and sys.argv[1] == 'plot'):
        synthesise()
    plot()