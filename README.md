# turbofit

A Python toolkit around [Turbospectrum v20](https://github.com/bertrandplez/Turbospectrum_NLTE) for stellar spectroscopy. It can:

* **synthesise spectra**, either one at a time or many abundance patterns for the same star;
* **build grids** for FERRE (Synple-style, Slurm-ready);
* **fit stellar parameters** with FERRE;
* **measure radial velocities**;
* **fit abundances line by line**.

You choose the model atmospheres:

| mode | atmospheres | radiative transfer | use it for |
|---|---|---|---|
| **CNO-LTE** | APOGEE MARCS grid with **[C/M] and [α/M] nodes** | LTE | Teff/log g/[M/H] fitting where the CNO content of the atmosphere matters; C, N, O from molecules; any element without NLTE data |
| **NLTE** | standard MARCS + departure coefficients | NLTE for the chosen elements | abundances of the 17 elements with public departure-coefficient grids (Al, Ba, Ca, Co, Cr, Eu, Fe, H, Mg, Mn, Na, Ni, O, Si, Sr, Ti, Y) |

The CNO-LTE mode never interpolates atmospheres. Spectra are computed at grid nodes and interpolated: by FERRE for grids, and across the 2⁵ bracketing nodes for line fitting. In NLTE mode, Turbospectrum's own interpolator builds the atmosphere and the departure coefficients at the exact parameters.

### What it does that other wrappers don't

* **CNO atmospheres and NLTE in one chain.** Stellar parameters come from FERRE on a CNO-node grid. Abundances are then measured in NLTE wherever departure grids exist and in CNO-LTE otherwise. Every abundance records which mode was used.
* **Automatic line masks for each star and instrument** (`linefit.automask`). The star is synthesised at its parameters (base, [X/M]+Δ, and X removed) and convolved with *your* instrument. Windows are kept where X has a strong response and other species contribute little. This also works for molecular features: C from CO/CN/CH, N from CN, O from OH. Each window records its response and blend fraction.
* **Per-line systematic errors at almost no extra cost.** In CNO-LTE mode the bracketing-node spectra are already computed. `d[X/H]/dTeff`, `d/dlogg`, `d/d[M/H]`, `d/d[α/M]` and `d/d[C/M]` therefore come from re-weighting them, without new syntheses. FERRE's uncertainties are propagated into every line and into the element means, treating the lines as correlated.
* **One LSF treatment everywhere** (`broaden.py`). The same code serves grids, RV and line fits. It accepts a constant R, an R(λ) that varies with wavelength, a callable kernel, or a per-pixel LSF array such as APOGEE's or NIRSpec's. Convolution, Doppler shift and resampling are a single sparse matrix evaluated at each output pixel's own wavelength.

## Files (flat layout: put the folder on `PYTHONPATH` or run from it)

| file | what |
|---|---|
| `synthesis.py` | `synthesize()`: one spectrum. LTE on any atmosphere file, or NLTE via `nlte=`/`nlte_elems=`, optionally restricted to `segments` (windows). `make_abundances()`. |
| `nlte.py` | `NLTEConfig`, bracketing of standard MARCS models, departure-coefficient interpolation, NLTE species file |
| `abund.py` | `synthesize_sets()` / `ladder()`: many abundance patterns for one star |
| `broaden.py` | vmac/vsini, R / R(λ) / LSF convolution, Doppler shift, resampling |
| `turbogrid.py`, `fill_grid_holes.py` | grid building → collection → hole filling → FERRE grid (see `turbogrid_user_guide.pdf`) |
| `ferre.py` | FERRE driver: resample spectra onto the grid's pixels, write inputs, run FERRE, read results |
| `rv.py` | CCF radial velocity (Zucker 2003 errors) + χ² refinement |
| `linefit.py` | line-by-line fitting, automatic masks, systematics, plots |
| `pipeline.py` | RV → FERRE → masks → line-by-line for one star |
| `atmos.py`, `atomic.py`, `utils.py` | atmosphere-file utilities (Holtzman), solar abundances, helpers |
| `examples/` | runnable templates for each of the above |

## Installation

1. **Turbospectrum v20**
   ```
   git clone https://github.com/bertrandplez/Turbospectrum_NLTE
   cd Turbospectrum_NLTE/exec-gf && make          # or exec/ with ifx/ifort
   export TURBOSPEC_DIR=/path/to/Turbospectrum_NLTE
   export PATH=$TURBOSPEC_DIR/exec-gf:$PATH          # babsma_lu, bsyn_lu
   ```
   `csh` must be installed: the wrapper drives Turbospectrum through small csh scripts.
2. **Python**: numpy, scipy, astropy, matplotlib.
3. **Atmospheres**:
   * CNO-LTE: the APOGEE MARCS grid with C/α nodes, read through `turbogrid.make_atmod_lookup`.
   * NLTE: the standard MARCS grid (ASCII `*.mod`, flat directory).
4. **NLTE data** (only for NLTE mode):
   * Compile the interpolator: `cd $TURBOSPEC_DIR/interpolator && gfortran -o interpol_modeles_nlte interpol_modeles_nlte_gfort.f`.
   * From <https://keeper.mpdl.mpg.de/d/6eaecbf95b88448f98a4/>, download the *MARCS* departure grids (`NLTEgrid*MARCS*.bin`, `auxData*MARCS*.dat`) and the model atoms for the elements you need. Lay them out as in `nlte_filenames.cfg`, with the model atoms in `<nlte_dir>/model_atoms/`.
   * NLTE line lists need level labels. The Gaia-ESO NLTE list at the same link has them. VALD lists can be labelled with TSFitPy's utilities.
5. **FERRE** (optional): `git clone https://github.com/callendeprieto/ferre && cd ferre/src && make`.

## Quick start

```python
import synthesis, abund, broaden, turbogrid, nlte

lookup = turbogrid.make_atmod_lookup('/data/marcs_apogee', flat=True)
atm = lookup(4500., 2.0, -0.5, 0.25, -0.25)               # a CNO node
flux, cont, wave = synthesis.synthesize(4500., 2.0, mh=-0.5, am=0.25, cm=-0.25, nm=0.5,
                                        vmicro=1.8, wrange=[15150, 15250], dw=0.01,
                                        atmod=atm, linelists=ll)
w, f = broaden.convolve(wave, flux / cont, R=22500, dlogw_out=6e-6, rv=-45.)

cfg = nlte.NLTEConfig('/data/marcs_std', '/data/nlte_data')
fn, cn, wn = synthesis.synthesize(4500., 2.0, mh=-0.5, vmicro=1.8, wrange=[5500, 5600],
                                  dw=0.01, linelists=ll, nlte=cfg, nlte_elems=['Mg', 'Ca'])
```

Line-by-line fitting is shown in `examples/example_linefit.py` and the full chain in `examples/example_pipeline.py`.

## Things worth knowing

* **babsma is rerun for every abundance set by default.** bsyn takes the electron pressure from babsma's opacity file. Reusing one opacity file across ±0.3 dex in Fe gave normalised-flux errors of 0.09 in a 4500 K giant (0.001–0.003 for Ti/Ni). `rerun_babsma=False` is available for trace elements.
* **Segments.** Line fitting synthesises only the windows, using Turbospectrum's `SEGMENTSFILE` at its default R = 500000 sampling. The `RESOLUTION` keyword is not written, because v20.1 reads it with an `(a)` format, which gives a garbage value.
* **Blends.** `linefit.fit(iterate=1)` refits after feeding the measured abundances back into the synthesis, so blending lines of fitted elements are modelled at the right abundance. Known abundances of other elements can be given as `star['elems']`.
* **vmac** defaults to `'global'`. A first pass fits it per line, then it is fixed to the response-weighted median. A free per-line vmac is degenerate with the abundance for weak lines in narrow windows.
* **FERRE ends with `STOP 1`**, a non-zero exit status, even after a successful run. `ferre.py` therefore checks the output file instead of the exit code. With the current FERRE, `COVPRINT=1` writes the full n×n covariance, not the upper triangle described in the manual. Points exactly on a grid edge are rejected, so starting values are nudged 0.001 inside.
* **FERRE grid axis order.** In grids from `write_ferre_grid`, the last label varies fastest. This was checked by FERRE interpolation at every node of a test grid. Log-λ grids are written with `LOGW = 1`.

## Testing status

Tested here against Turbospectrum v20.1 (gfortran) and FERRE (GitHub head):

* LTE synthesis and segments.
* The NLTE code path in bsyn, run with an all-LTE species file.
* MARCS bracketing, including holes.
* The departure-file NaN clean-up.
* Abundance sets and ladders.
* Broadening and LSF (FWHM and equivalent-width conservation, R(λ), array LSF, RV shift).
* CCF and χ² RV pulls (≈1σ).
* Grid → FERRE grid → FERRE node interpolation and fit.
* Line-by-line injection–recovery. On noise-free data, Fe is recovered to ±0.003 dex.

**Not yet run end to end:** NLTE synthesis with real departure grids. The data could not be downloaded in the build environment. The interpolator inputs and species file follow TSFitPy exactly, but please run one NLTE spectrum against TSFitPy before relying on it.

## Credits and citations

* Turbospectrum: Alvarez & Plez (1998); Plez (2012); NLTE version Gerber et al. (2023, A&A 669, A43). Cite the NLTE data/model-atom papers listed in Turbospectrum's `public-ref-NLTE-TS.bib`.
* NLTE interpolation inputs, species file and departure-file clean-up follow **TSFitPy** (Gerber et al. 2023; Storm & Bergemann 2023), GPL-3.
* Synthesis wrapper based on **David Nidever**'s Turbospectrum Python wrapper; `atmos.py` by **Jon Holtzman** (BSD-3).
* MARCS: Gustafsson et al. (2008). FERRE: Allende Prieto et al. (2006).
* RV errors: Zucker (2003).

Licence: GPL-3.0 (see `LICENSE`).
