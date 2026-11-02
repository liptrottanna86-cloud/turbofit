# turbofit user guide

*Anna Liptrott, Astrophysics Research Institute, Liverpool John Moores University*

turbofit is a Python toolkit around **Turbospectrum v20**, an ASPCAP-style analysis that anyone can run on one spectrum or many, from any spectrograph. It covers:

* synthesising spectra, one at a time or many abundance patterns per star;
* stellar parameters, either in **local mode** (on-demand syntheses around each star, no grid needed) or in **grid mode** (FERRE on a grid you build);
* measuring radial velocities and normalising continua;
* fitting abundances (and ¹²C/¹³C) line by line, with measurement and systematic errors.

There are no empirical calibrations and no instrument presets: you describe your instrument (R, R(λ) or an LSF).

For each element you choose how it is synthesised:

* **CNO-LTE**: APOGEE MARCS atmospheres with [C/M] and [α/M] nodes, in LTE.
* **NLTE**: standard MARCS atmospheres plus departure coefficients, for elements that have NLTE data.

![turbofit workflow](img/workflow.svg)

## Contents

| | chapter | what it covers |
|---|---|---|
| 1 | [Installation](01-installation.md) | Turbospectrum, csh, atmospheres, NLTE data, APOGEE LSF files, FERRE |
| 2 | [Synthesising spectra](02-synthesis.md) | `synthesize()`, abundances, ¹²C/¹³C, segments, abundance sets and ladders |
| 3 | [NLTE](03-nlte.md) | `NLTEConfig`, what happens under the hood, choosing elements |
| 4 | [Broadening, LSFs and resampling](04-broadening.md) | vmac/vsini, R, R(λ), LSF arrays, the APOGEE LSF, Doppler shift, output grids |
| 5 | [Continuum normalisation](05-continuum.md) | automatic chunks, polynomial and running-percentile continua |
| 6 | [Building grids for FERRE](06-grids.md) | turbogrid workflow, on-the-fly LSF, vmicro and broadening axes, NLTE grids |
| 7 | [Radial velocities](07-radial-velocity.md) | CCF and χ² refinement |
| 8 | [Stellar parameters in local mode](08-local-mode.md) | full-spectrum fit with on-demand syntheses, no grid needed |
| 9 | [Stellar parameters with FERRE (grid mode)](09-ferre.md) | the FERRE driver, resampling, fixed parameters |
| 10 | [Line-by-line abundances](10-line-by-line.md) | modes, masks, fitting, broadening, errors, ¹²C/¹³C, systematics |
| 11 | [Analysing a star](11-analysis.md) | `turbofit.analyse()`: local or grid mode → masks → abundances |
| 12 | [Troubleshooting and known quirks](12-troubleshooting.md) | errors you may see and what to do |
| 13 | [API reference](13-api-reference.md) | every public function and its arguments |

## Conventions used throughout

* **Abundances.** `mh` is [M/H], `am` is [α/M], `cm` is [C/M] and `nm` is [N/M].
    * Individual elements are given as `[X/M]` with `elems=[['Mg', 0.3], ['Al', -0.1]]`.
    * The α elements are O, Ne, Mg, Si, S, Ar, Ca and Ti.
    * Solar abundances are Grevesse et al. (2007), which is also Turbospectrum's default.

* **Wavelengths** are in Å. They follow your line lists, so air or vacuum depending on what you use.
* **Spectra passed to the fitting routines** must be continuum-normalised, with 1σ errors per pixel.
    * Bad pixels can be marked with `err = inf` or `nan`.
* **Radial velocities** are in km/s, with positive meaning receding. Observed spectra stay in the observed frame; turbofit applies the shift.
* **Imports are flat.** Put the turbofit folder on your `PYTHONPATH` (or `sys.path`) and use `import synthesis, linefit, ...`.

## [X/M] or [X/Fe]?

turbofit works in **[X/M]** throughout: in the inputs (`elems`, `am`, `cm`, `nm`) and in the line-by-line results (`xm`). This follows the APOGEE/ASPCAP convention, where [M/H] (FERRE's METALS) is the overall scaling and each element is quoted relative to it.

* `mh` is added to *every* element heavier than He, **Fe included**. So by default [Fe/H] = [M/H] and **[X/M] = [X/Fe]**.
* An override such as `elems=[['Mg', 0.3]]` means A(Mg) = A☉(Mg) + [M/H] + 0.3, which is [Mg/M] = 0.3.
* The two differ as soon as Fe gets its own value. That happens when you set it (`elems=[['Fe', 0.1]]`), or when Fe is fitted line by line and `iterate` feeds the measured [Fe/M] back in.
* The results also give `xh` = [X/H] = [X/M] + [M/H], with [M/H] the star's value (from FERRE).

When Fe has its own value:

**[X/Fe] = [X/M] − [Fe/M] = [X/H] − [Fe/H]**

In code, that is `res.summary['X']['xh'] - res.summary['Fe']['xh']`.
