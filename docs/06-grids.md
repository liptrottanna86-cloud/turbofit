# 6. Building grids for FERRE

`turbogrid.py` builds Synple-style grids with Turbospectrum. The whole workflow is described step by step in `turbogrid_user_guide.pdf`:

1. grid table;
2. Slurm or GNU parallel jobs;
3. collect;
4. fill holes;
5. convolve;
6. write the FERRE grid.

This chapter covers what is new in turbofit and gives the short version.

## 6.1 The short version

```python
import numpy as np, turbogrid as tg

lookup = tg.make_atmod_lookup('/data/marcs_apogee', model='marcs', flat=True)
params = tg.make_grid_params(teff_grid=(15, 3000., 250.), logg_grid=(12, -0.5, 0.5),
                             mh_grid=(14, -2.5, 0.25), cm=(11, -1.5, 0.25),
                             am=(7, -1.0, 0.25), nm=(6, -0.5, 0.5),
                             vmicro='holtzman', atmod_lookup=lookup)

tg.write_grid_table(params, outdir='gridm31', wrange=(15050., 17050.), dw=0.01,
                    linelists=ll, atmos_type='marcs',
                    outputs_dir='/scratch/me/m31grid', shard_size=1000,
                    broaden_kw=dict(R=22500., dlogw_out=6e-6))      # new, see 5.2
tg.write_slurm_array('gridm31', jobname='m31grid', time='24:00:00', mem='8G',
                     partition='nodes', array_step=1011, max_concurrent=30,
                     max_array_size=999)
#   sbatch gridm31/submit_array.slurm      ... then:
tg.collect_grid('gridm31', gridfile='/scratch/me/m31grid/grid.npz')
#   python fill_grid_holes.py /scratch/me/m31grid/grid.npz --labels teff logg mh cm am nm ...
tg.write_ferre_grid('/scratch/me/m31grid/grid_filled.npz', 'gridm31/grid.ferre',
                    labels=('teff', 'logg', 'mh', 'cm', 'am', 'nm'),
                    resolution=22500, vacuum=1)
```

## 6.2 Broadening on the fly with any LSF: `broaden_kw`

`write_grid_table(..., broaden_kw=dict(...))` passes the dictionary straight to `broaden.convolve()` right after each spectrum is synthesised. Only the broadened spectrum is stored. Anything from [chapter 4](04-broadening.md) works:

```python
broaden_kw = dict(R=22500., dlogw_out=6e-6)                        # constant R, log-λ pixels
broaden_kw = dict(R=(wave_R, R_of_wave), dw_out=0.1)               # R varying with λ
broaden_kw = dict(lsf=lsf_array, lsf_offsets=np.arange(-7, 7.01, 0.1),
                  wave_out=instrument_wave)                         # per-pixel LSF
```

The older `R_out`, `dw_out` and `wout` options still work. `broaden_kw` takes priority if both are given.

> [!TIP]
> Make `wrange` a few resolution elements wider than the range you want in the output, so convolution edge effects fall outside it.

### The APOGEE LSF on the apStar grid

```python
import apogee_lsf
alsf = apogee_lsf.ApogeeLSF('/data/apogee_lsf')
wave = apogee_lsf.apStarWavegrid()
broaden_kw = dict(wave_out=wave, **alsf.inst(wave, fiber='combo'))
```

The resulting grid is log-λ, so `write_ferre_grid()` writes `LOGW = 1` (6.4).

### A fixed ¹²C/¹³C for the whole grid

`write_grid_table(..., c12c13=8.)` uses that ratio for every spectrum. Without it, the `solarisotopes` default of `synthesize()` applies.

## 6.3 NLTE grids

```python
import nlte
cfg = nlte.NLTEConfig('/data/marcs_standard', '/data/nlte_data')
params = tg.make_grid_params(teff_grid=..., logg_grid=..., mh_grid=..., am=..., nm=...,
                             vmicro='holtzman', atmod_lookup=tg.nlte_atmod_lookup)
tg.write_grid_table(params, outdir='grid_nlte', ..., nlte=cfg, nlte_elems=['Mg', 'Ca', 'Fe'])
```

`nlte_atmod_lookup` tells turbogrid that each atmosphere will be interpolated from the standard MARCS grid when the point runs. `am` and `cm` then only change the synthesis abundances.

## 6.4 vmicro and broadening as grid axes (as ASPCAP)

For a grid mode that fits vmicro and the stellar broadening like ASPCAP:

```python
params = tg.make_grid_params(..., vmicro=(4, -0.3, 0.2), vmicro_log=True,   # 0.5..2.0 km/s
                             atmod_lookup=lookup)
tg.write_grid_table(params, ..., broaden_kw=dict(R=22500., dlogw_out=6e-6),
                    broaden_axis=('vsini', 10 ** np.arange(0.0, 1.81, 0.3)))   # or ('vmac', ...)
...
tg.write_ferre_grid('grid_filled.npz', 'grid.ferre',
                    labels=('teff', 'logg', 'mh', 'cm', 'am', 'nm', 'lgvmicro', 'lgvsini'))
```

* **`vmicro_log=True`** makes the vmicro axis even in log₁₀ vmicro. The label `'lgvmicro'` is written as `LOG10VDOP`, ASPCAP's name.
* **`broaden_axis=(name, values)`** adds a vsini or vmac axis at **no extra synthesis cost**. Every synthesised spectrum is broadened with each value, then `broaden_kw` is applied. `collect_grid()` expands the grid and adds the columns `vsini`/`lgvsini` (or `vmac`/`lgvmac`).
  * Use the log label (`'lgvsini'` → `LGVSINI`, `'lgvmac'` → `LGVMACRO`) with log-spaced values. Values must be > 0 for a log axis.
  * Disk space grows with the number of broadening values; synthesis time does not.
* `ferre.to_star()` converts `LOG10VDOP`, `LGVSINI` and `LGVMACRO` back to km/s, with errors propagated.

## 6.5 FERRE grids on log-λ pixels

`write_ferre_grid()` now detects the wavelength grid:

* linear grids get `LOGW = 0` and `WAVE = λ₀ Δλ`;
* log₁₀-uniform grids (e.g. from `dlogw_out`) get `LOGW = 1` and `WAVE = log₁₀λ₀ Δlog₁₀λ`.

## 6.6 Axis order

In the flux block, **the last label in `labels` varies fastest**, as FERRE expects. Use the same `labels` order in `fill_grid_holes.py` and `write_ferre_grid()`.

> [!NOTE]
> Set `vacuum=1` if your line lists are in vacuum wavelengths (APOGEE lists are). The flag should describe your line lists, not Turbospectrum's default.
