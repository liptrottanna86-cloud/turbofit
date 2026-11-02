# 9. The full pipeline

`pipeline.analyse()` runs the whole chain for one star:

1. **RV template:** FERRE interpolates a spectrum from your grid, at the grid centre or at `template_params`.
2. **RV:** cross-correlation, then χ² refinement.
3. **FERRE:** fit on the CNO grid, giving Teff, log g, [M/H], [α/M], [C/M] and [N/M] with errors.
4. **Masks:** `automask` for your elements (or your own `mask=`).
5. **Line by line:** `LineFitter.fit()`, using FERRE's errors as `star_err` for the systematics.

```python
import numpy as np, turbogrid, linefit, pipeline

wave, flux, err = np.loadtxt('star1_normalised.txt', unpack=True)
lf = linefit.LineFitter(linelists=ll, inst=dict(R=22500.),
                        nodes=dict(teff=(15, 3000., 250.), logg=(12, -0.5, 0.5),
                                   mh=(14, -2.5, 0.25), am=(7, -1.0, 0.25),
                                   cm=(11, -1.5, 0.25)),
                        atmod_lookup=turbogrid.make_atmod_lookup('/data/marcs_apogee', flat=True),
                        n_jobs=16)

out = pipeline.analyse(wave, flux, err, grid='gridm31/grid.ferre', lf=lf,
                       elements=['C', 'N', 'O', 'Mg', 'Al', 'Fe'],
                       ferre_exe='/path/to/ferre.x', workdir='star1', nthreads=8,
                       vmac='global', iterate=1)          # extra keywords go to lf.fit

out['rv'], out['rv_err']
out['star'], out['star_err']        # from FERRE
out['ferre']                        # the raw FERRE result (chi2, cov, model, edge flags)
out['mask']                         # the windows used
out['lines'].summary                # abundances per element
```

Everything is written under `workdir`:

* `template/`: the FERRE run that made the RV template;
* `ferre/`: the FERRE fit itself, including the log.

Pass `vmicro=` to override the vmicro taken from the grid or the Holtzman relation, and `rv_range=` to change the velocity search window.

> [!NOTE]
> The grid must contain TEFF and LOGG axes for the line-by-line step, because the star dictionary needs them. Any other grid labels (METALS, ALPHA, C, N, VMICRO/LOG10VDOP) are mapped automatically by `ferre.to_star()`.
