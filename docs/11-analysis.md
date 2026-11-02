# 11. Analysing a star: local or grid mode

`turbofit.analyse()` runs the whole chain for one spectrum:

1. stellar parameters, in either **local** or **grid** mode;
2. line masks;
3. line-by-line abundances (plus ¹²C/¹³C if requested), with measurement and systematic errors.

It works with any spectrograph: you describe the instrument with `inst=`. There are no instrument presets and no empirical calibrations.

```python
import turbofit, turbogrid

out = turbofit.analyse(
    wave, flux, err,                                  # observed frame, flux not normalised
    linelists=['atomic.list', 'molecular.list', 'Hlinedata'],
    inst=dict(R=22500.),                              # any instrument (chapter 4)
    elements=['C', 'N', 'O', 'Mg', 'Al', 'Fe', 'C12C13'],
    mode='local',                                     # or 'grid'
    nodes=dict(teff=(15, 3000., 250.), logg=(12, -0.5, 0.5), mh=(14, -2.5, 0.25),
               am=(7, -1.0, 0.25), cm=(11, -1.5, 0.25)),
    atmod_lookup=turbogrid.make_atmod_lookup('/data/marcs_apogee', flat=True),
    guess=dict(teff=4600, logg=2.2, mh=-0.3),
    nlte=None, nlte_elems=[],                          # optional NLTE for abundances
    n_jobs=16, cache_dir='/scratch/me/turbofit_cache', workdir='star1')

out['star'], out['star_err']        # Teff, logg, [M/H], [a/M], [C/M], [N/M], vmicro, vmac/vsini, RV
out['params']                       # the full parameter-stage result
out['lines'].summary                # abundances per element (and 12C/13C)
out['lines'].write('star1_lines.csv')
```

## 11.1 Choosing the mode

| | **local** (`mode='local'`) | **grid** (`mode='grid'`, `grid='grid.ferre'`) |
|---|---|---|
| parameters from | `starfit` (chapter 8): syntheses on demand around this star | FERRE (chapter 9) on a turbogrid grid (chapter 6) |
| needs | the atmosphere nodes and line lists | a grid built beforehand for your instrument and range |
| best for | one star, a handful, or a few hundred; any instrument | large samples that share one instrument setup |
| data normalisation | not needed (continuum fitted with the model) | done for you with `continuum.normalise` |
| fitted | everything in `starfit`, including vmicro and vmac/vsini | whatever axes the grid has (add `lgvmicro` and `lgvsini`/`lgvmac`, 6.4) |
| fixing parameters | `fixed=dict(logg=2.45)` | `fixed={'LOGG': 2.45}` (FERRE labels), with `fit_labels=` |

**The abundance stage is the same in both modes.** It uses the parameter-stage RV and, by default, its broadening (`fit_kw=dict(vmac='star')`), with the parameter errors feeding the systematics. Two further options apply there:

* **NLTE:** `nlte=` / `nlte_elems=` treat those elements in NLTE where data exist.
* **Line windows:** `mask=` uses your own windows instead of `automask` (`automask_kw=` passes options to `automask`).

## 11.2 The pieces it uses

| stage | module |
|---|---|
| continuum (grid mode, and for the line stage) | `continuum.py` (chapter 5) |
| parameters, local | `starfit.StarFitter` (chapter 8) |
| parameters, grid | `pipeline.analyse(..., params_only=True)` → RV from a FERRE-interpolated template, then `ferre.fit` (chapter 9) |
| masks and abundances | `linefit.LineFitter` (chapter 10) |

Extra options go to each piece through `starfit_kw=`, `linefit_kw=` and `fit_kw=`.

## 11.3 Test of the whole chain

The test used one simulated spectrum:

* two regions, 5400–5500 Å and 8000–8090 Å;
* R = 20 000, S/N 150;
* a curved, non-normalised continuum;
* VALD atoms and an artificial CN test list.

It was run through `turbofit.analyse(mode='local', elements=['Fe', 'Ti', 'C12C13'])`, with Teff, log g, [α/M] and [C/M] fixed (only two atmosphere nodes were available here):

| quantity | injected | recovered |
|---|---|---|
| [M/H] | 0.00 | −0.005 |
| [N/M] | 0.25 | 0.247 |
| vmicro | 1.70 km/s | 1.718 |
| vmac | 4.5 km/s | 4.51 |
| RV | 37.3 km/s | 37.33 |
| [Fe/M] (8 lines) | 0.00 | 0.010 ± 0.014 |
| [Ti/M] (2 lines) | 0.00 | −0.013 ± 0.019 |
| ¹²C/¹³C (3 windows) | 15 | 15.4 (13.3–17.7) |

The whole run took 6 minutes on 2 cores, starting from an empty cache.

## 11.4 Many stars

There is nothing survey-specific. Loop over your spectra:

* In local mode, give every star the same `cache_dir`, so stars with similar parameters reuse each other's syntheses.
* In grid mode, give each star its own `workdir`.

`pipeline.analyse()` remains available for the grid route when you already have a `LineFitter` set up.
