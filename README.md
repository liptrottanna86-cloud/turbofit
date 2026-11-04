# turbofit

A Python wrapper for Turbospectrum and FERRE: synthetic spectra, FERRE-ready model grids and spectral fitting.

> **Work in progress.** The code is under active development and is not yet validated for science use. Interfaces may change without notice.

turbofit is adapted from David Nidever's packages, in particular
[dnidever/turbospectrum](https://github.com/dnidever/turbospectrum) (GPL-3.0).

## Citing

If you use turbofit, please cite the codes and models it depends on:

- **Turbospectrum:** Alvarez & Plez 1998, A&A, 330, 1109; Plez 2012, ASCL, ascl:1205.004
- **Turbospectrum NLTE (v20):** Gerber et al. 2023, A&A, 669, A43
- **FERRE:** Allende Prieto et al. 2006, ApJ, 636, 804
- **MARCS model atmospheres:** Gustafsson et al. 2008, A&A, 486, 951

If you use the APOGEE line lists, also cite Smith et al. 2021, AJ, 161, 254.

## Licence

GPL-3.0. See `LICENSE`.

## Documentation

The user guide is in [`docs/`](docs/), as one Markdown page per topic, starting with [`docs/README.md`](docs/README.md). A PDF version is in [`docs/turbofit_user_guide.pdf`](docs/turbofit_user_guide.pdf).

Topics: installation, synthesis, NLTE, broadening, continuum normalisation, grids, radial velocities, local mode, FERRE, line-by-line fitting, analysis, troubleshooting and the API reference. Worked examples are in [`examples/`](examples/).
