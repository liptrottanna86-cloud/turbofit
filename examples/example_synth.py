"""One spectrum, an abundance ladder, and an NLTE spectrum.

Run from anywhere; set TURBOSPEC_DIR and put babsma_lu/bsyn_lu on PATH first:
    export TURBOSPEC_DIR=/path/to/Turbospectrum_NLTE
    export PATH=$TURBOSPEC_DIR/exec-gf:$PATH
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
import synthesis, abund, broaden, turbogrid

linelists = ['/path/to/atomic_vac.list', '/path/to/molecular_vac.list',
             os.path.join(os.environ['TURBOSPEC_DIR'], 'DATA', 'Hlinedata')]

# APOGEE C/alpha/M MARCS grid -> Turbospectrum-ready atmosphere for a NODE
lookup = turbogrid.make_atmod_lookup(griddir='/path/to/marcs_apogee', model='marcs', flat=True)
atm = lookup(4500., 2.0, -0.5, 0.25, -0.25)          # teff, logg, [M/H], [a/M], [C/M]

# 1) one spectrum (LTE, CNO atmosphere)
flux, cont, wave = synthesis.synthesize(
    4500., 2.0, mh=-0.5, am=0.25, cm=-0.25, nm=0.5, vmicro=1.8,
    wrange=[15150., 15250.], dw=0.01, atmod=atm, linelists=linelists)

# ... broadened to APOGEE-like R and resampled on a log-lambda grid
w, f = broaden.convolve(wave, flux / cont, R=22500., dlogw_out=6e-6)

# 2) the same star with several abundance patterns (one babsma+bsyn each)
wave2, fluxes, conts = abund.synthesize_sets(
    4500., 2.0, mh=-0.5, am=0.25, cm=-0.25, nm=0.5, vmicro=1.8,
    abund_sets=[[], [['Mg', 0.5]], [['Mg', 0.5], ['Al', -0.3]]],
    wrange=[15150., 15250.], dw=0.01, atmod=atm, linelists=linelists)

# ... or a ladder in one element
xm, wave3, fl, co = abund.ladder(4500., 2.0, 'Mg', offsets=[-0.4, -0.2, 0, 0.2, 0.4],
                                 mh=-0.5, am=0.25, cm=-0.25, nm=0.5, vmicro=1.8,
                                 wrange=[15150., 15250.], dw=0.01, atmod=atm,
                                 linelists=linelists)

# 3) NLTE Mg and Ca (standard MARCS + departure coefficients); atmod not needed
import nlte
cfg = nlte.NLTEConfig(marcs_dir='/path/to/marcs_standard', nlte_dir='/path/to/nlte_data')
print('NLTE data available for:', cfg.available_elements())
fn, cn, wn = synthesis.synthesize(4500., 2.0, mh=-0.5, am=0.25, vmicro=1.8,
                                  wrange=[5500., 5600.], dw=0.01, linelists=linelists,
                                  nlte=cfg, nlte_elems=['Mg', 'Ca'])
