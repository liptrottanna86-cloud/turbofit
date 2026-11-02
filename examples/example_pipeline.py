"""RV -> FERRE -> automatic masks -> line-by-line, for one star."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
import turbogrid, linefit, pipeline

wave, flux, err = np.loadtxt('my_star_normalised.txt', unpack=True)
lf = linefit.LineFitter(
    linelists=['/path/to/atomic_vac.list', '/path/to/molecular_vac.list', '/path/to/Hlinedata'],
    inst=dict(R=22500.),
    nodes=dict(teff=(15, 3000., 250.), logg=(12, -0.5, 0.5), mh=(14, -2.5, 0.25),
               am=(7, -1.0, 0.25), cm=(11, -1.5, 0.25)),
    atmod_lookup=turbogrid.make_atmod_lookup('/path/to/marcs_apogee', flat=True),
    n_jobs=16)
out = pipeline.analyse(wave, flux, err, grid='/path/to/gridm31/grid.ferre', lf=lf,
                       elements=['C', 'N', 'O', 'Mg', 'Al', 'Fe'],
                       ferre_exe='ferre.x', workdir='star1', nthreads=8)
print(out['rv'], out['star'])
print(out['lines'].summary)
