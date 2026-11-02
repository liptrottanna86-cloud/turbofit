"""Line-by-line abundances for one star with known parameters."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
import turbogrid, linefit, nlte

linelists = ['/path/to/atomic_vac.list', '/path/to/molecular_vac.list', '/path/to/Hlinedata']
wave, flux, err = np.loadtxt('my_star_normalised.txt', unpack=True)   # observed frame

lf = linefit.LineFitter(
    linelists=linelists,
    inst=dict(R=22500.),                         # or lsf=..., see broaden.py
    nodes=dict(teff=(15, 3000., 250.), logg=(12, -0.5, 0.5), mh=(14, -2.5, 0.25),
               am=(7, -1.0, 0.25), cm=(11, -1.5, 0.25)),
    atmod_lookup=turbogrid.make_atmod_lookup('/path/to/marcs_apogee', flat=True),
    nlte=nlte.NLTEConfig('/path/to/marcs_standard', '/path/to/nlte_data'),
    nlte_elems=['Mg', 'Al', 'Ca'],                # NLTE where data exist; rest CNO-LTE
    n_jobs=16)

star = dict(teff=4520., logg=1.9, mh=-0.42, am=0.18, cm=-0.20, nm=0.30, vmicro=1.6, rv=-121.4)
star_err = dict(teff=50., logg=0.1, mh=0.05, am=0.03, cm=0.05)

# automatic, star- and instrument-specific line windows (inspect / edit the file!)
mask = lf.automask(['C', 'N', 'O', 'Mg', 'Al', 'Ca', 'Fe'], star, wave,
                   min_response=0.02, max_blend=0.3,
                   exclude=[(15890., 15960.)])  # e.g. a chip gap, observed frame
linefit.write_mask(mask, 'my_star_mask.txt')
# or: mask = linefit.read_mask('my_star_mask.txt'); TSFitPy masks: read_mask(f, element='Fe')

res = lf.fit(wave, flux, err, star, mask=mask, star_err=star_err,
             vmac='global', iterate=1, sys_vmicro=False)
for el, s in res.summary.items():
    print(el, s)
res.write('my_star_lines.csv')
linefit.plot_lines(res, 'my_star_lines.pdf')
