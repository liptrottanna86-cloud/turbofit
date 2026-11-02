"""One spectrum, any instrument: stellar parameters (local or grid mode) ->
line masks -> abundances (+ 12C/13C). No empirical calibrations, no presets."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
import turbogrid, turbofit

wave, flux, err = np.loadtxt('my_star.txt', unpack=True)        # observed frame, not normalised

nodes = dict(teff=(15, 3000., 250.), logg=(12, -0.5, 0.5), mh=(14, -2.5, 0.25),
             am=(7, -1.0, 0.25), cm=(11, -1.5, 0.25))
lookup = turbogrid.make_atmod_lookup('/path/to/marcs_apogee', flat=True)

# describe YOUR instrument: a resolving power, R(lambda), or an LSF array
inst = dict(R=22500.)
# inst = dict(R=(wave_R, R_of_wave))                                     # e.g. NIRSpec R(lambda)
# import apogee_lsf; inst = apogee_lsf.ApogeeLSF('/path/to/lsf').inst(wave)   # APOGEE

out = turbofit.analyse(
    wave, flux, err, linelists=['/path/to/atomic.list', '/path/to/molecular.list', '/path/to/Hlinedata'],
    inst=inst, elements=['C', 'N', 'O', 'Mg', 'Al', 'Fe', 'C12C13'],
    mode='local',                                   # or mode='grid', grid='/path/to/grid.ferre'
    nodes=nodes, atmod_lookup=lookup,
    guess=dict(teff=4600, logg=2.2, mh=-0.3),
    # fixed=dict(logg=2.45),                        # e.g. from asteroseismology
    broadening='vmac',                              # 'vsini' for rotating stars
    n_jobs=16, cache_dir='/scratch/me/turbofit_cache', workdir='my_star_run')

print(out['star'])
print(out['lines'].summary)
out['lines'].write('my_star_lines.csv')
