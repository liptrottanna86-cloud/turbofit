# turbofit

*Anna Liptrott, Astrophysics Research Institute, Liverpool John Moores University*

> [!WARNING]
> **Work in progress.** turbofit is under active development. Interfaces may change, and parts of it are not yet validated for science use (see [Status](#status)). Check results against an independent code before relying on them.

turbofit is a Python wrapper for the spectral-synthesis code [Turbospectrum v20](https://github.com/bertrandplez/Turbospectrum_NLTE) and the spectral-fitting code [FERRE](https://github.com/callendeprieto/ferre). It does an ASPCAP-style analysis (the APOGEE pipeline's approach) of one spectrum or many, from any spectrograph (e.g. APOGEE, Gaia RVS, FLAMES, JWST NIRSpec). It is adapted from David Nidever's Turbospectrum Python wrapper ([dnidever/turbospectrum](https://github.com/dnidever/turbospectrum)).

What it does:

* **Synthesis**: computes model spectra for one star, for many abundance patterns, or for a series of ¹²C/¹³C ratios.
* **Stellar parameters**: fits Teff, log g, [M/H], [α/M], [C/M], [N/M], microturbulence and macroturbulence/rotation, in one of two modes:
  * **local mode** computes only the model spectra needed around each star, as it goes, so no pre-computed grid is needed;
  * **grid mode** builds a full grid of model spectra once (with `turbogrid`), then fits stars against it with FERRE. This is faster when there are many stars.
* **Radial velocities and continuum normalisation.**
* **Line-by-line abundances and ¹²C/¹³C**, with automatic choice of which pixels to fit, measurement errors, and the extra error from uncertain stellar parameters.

There are no empirical calibrations and no built-in instrument settings. You describe your instrument yourself: a resolving power R, an R that varies with wavelength, or a line-spread function (LSF, the instrument's profile for an unresolved line), for example the APOGEE LSF.

```python
import turbofit, turbogrid
out = turbofit.analyse(wave, flux, err, linelists=ll, inst=dict(R=22500.),
                       elements=['C', 'N', 'O', 'Mg', 'Fe', 'C12C13'], mode='local',
                       nodes=nodes, atmod_lookup=turbogrid.make_atmod_lookup('/path/to/marcs_apogee', flat=True),
                       guess=dict(teff=4600, logg=2.2), n_jobs=16, cache_dir='/path/to/cache')
```

### Two ways of modelling the star

| mode | model atmospheres | how the lines are computed | use it for |
|---|---|---|---|
| **CNO-LTE** | the APOGEE MARCS grid, which has separate grid points in [C/M] and [α/M] | LTE (local thermodynamic equilibrium) | stellar parameters, where the C, N and O content of the atmosphere matters; C, N and O from molecular lines; any element without NLTE data |
| **NLTE** | the standard MARCS grid, plus pre-computed NLTE corrections ("departure coefficients") | non-LTE for the elements you choose, LTE for the rest | abundances of the 17 elements with public NLTE grids: Al, Ba, Ca, Co, Cr, Eu, Fe, H, Mg, Mn, Na, Ni, O, Si, Sr, Ti, Y |

In **CNO-LTE** mode, atmospheres are never interpolated. Spectra are computed at the grid points and the *spectra* are interpolated: by FERRE in grid mode, or between the 2⁵ = 32 surrounding grid points in line fitting. In **NLTE** mode, Turbospectrum's own interpolator builds the atmosphere and the NLTE corrections at the star's exact parameters.

### What it does that other wrappers don't

* **Stellar parameters without a pre-built grid (local mode).** It fits the whole spectrum for Teff, log g, [M/H], [α/M], [C/M], [N/M], microturbulence and broadening. The model spectra are computed only at the atmosphere grid points around each star and cached for re-use. The instrument profile is applied on your own pixels, and the continuum is fitted at the same time, so you don't need to normalise the spectrum first.
* **CNO atmospheres and NLTE in one chain.** Stellar parameters come from the CNO-LTE atmospheres. Abundances are then measured in NLTE where NLTE grids exist, and in CNO-LTE otherwise. Every abundance records which mode was used.
* **Automatic line masks for each star and instrument** (`linefit.automask`). turbofit synthesises the star three ways: as it is, with element X raised, and with X removed. It convolves each with *your* instrument profile. Pixels are kept where X changes the spectrum strongly and other species contribute little. This works for molecules too: C from CO/CN/CH, N from CN, O from OH. Each window records how strongly it responds and how blended it is.
* **A measurement error on every abundance.** It comes from the shape of the χ² curve (including asymmetric errors), with optional refits on noise-added copies. Each element's error is the larger of the propagated error and the line-to-line scatter.
* **The effect of stellar-parameter errors on each line, at almost no extra cost.** In CNO-LTE mode the spectra at the surrounding grid points are already computed, so how each abundance changes with Teff, log g, [M/H], [α/M] and [C/M] is found by re-weighting them, without new syntheses. FERRE's parameter uncertainties are then propagated into every line and into each element's mean.
* **Carbon isotope ratio.** ¹²C/¹³C can be set per star when synthesising, in grids, or as a series of ratios. It can also be fitted line by line like an element (`'C12C13'`).
* **One instrument-profile treatment everywhere** (`broaden.py`). The same code is used for grids, radial velocities and line fits. It accepts a constant R, an R(λ), a function, or a per-pixel LSF array such as APOGEE's or NIRSpec's. Broadening, Doppler shift and resampling onto your pixels are done in one step.

## Documentation

* **[User guide (Markdown)](docs/README.md)**: one page per topic — installation, synthesis, NLTE, broadening, continuum, grids, radial velocities, local mode, FERRE, line-by-line fitting, analysis, troubleshooting and the API reference.
* **[User guide (PDF)](docs/turbofit_user_guide.pdf)**: the same content as one document. Rebuild it with `python docs/build_pdf.py`.
* **[`examples/`](examples/)**: runnable scripts. Each has a `CONFIG` block at the top; edit the paths there before running.

## Installation

turbofit itself is plain Python. It drives three external codes and needs some downloaded data. Which pieces you need depends on what you want to do:

| you want to… | Python + turbofit | Turbospectrum | APOGEE MARCS grid | standard MARCS grid | NLTE data | FERRE |
|---|:---:|:---:|:---:|:---:|:---:|:---:|
| compute LTE spectra | ✓ | ✓ | ✓ | | | |
| compute NLTE spectra | ✓ | ✓ | | ✓ | ✓ | |
| fit stellar parameters, local mode | ✓ | ✓ | ✓ | | | |
| build grids and fit with FERRE | ✓ | ✓ | ✓ | | | ✓ |
| line-by-line abundances | ✓ | ✓ | ✓ | for NLTE | for NLTE | |

All paths below (`/path/to/...`) are placeholders: put things wherever suits your machine. The model grids are large, so use a disk with plenty of space. On a cluster, run long downloads in `screen`/`tmux` or as a batch job.

### 1. Python and turbofit

You need Python ≥ 3.9 with numpy, scipy, astropy and matplotlib. A virtual environment keeps these separate from your system Python:

```bash
python3 -m venv /path/to/envs/turbofit
source /path/to/envs/turbofit/bin/activate
pip install numpy scipy astropy matplotlib
```

turbofit is not installed with `pip`. Its modules are imported directly (`import synthesis`, `import broaden`, …), so clone it and add the folder to your `PYTHONPATH`:

```bash
git clone https://github.com/liptrottanna86-cloud/turbofit /path/to/turbofit
export PYTHONPATH=/path/to/turbofit:$PYTHONPATH
```

To set this up automatically, add the `export` line to the end of `/path/to/envs/turbofit/bin/activate`, so it is set every time you activate the environment.

> [!TIP]
> If you have another Turbospectrum wrapper on your `PYTHONPATH`, Python may load its `synthesis.py` or `atmos.py` instead of turbofit's. Check with
> `python -c "import synthesis; print(synthesis.__file__)"`. It should print a path inside your turbofit folder.

### 2. Turbospectrum v20

Turbospectrum computes the model spectra. It is Fortran, so you need `gfortran` (or the Intel compilers).

```bash
git clone https://github.com/bertrandplez/Turbospectrum_NLTE /path/to/Turbospectrum_NLTE
cd /path/to/Turbospectrum_NLTE/exec-gf && make        # Intel compilers: use exec/ instead
```

Then tell turbofit where it is, e.g. in your `~/.bashrc` or the venv's `activate`:

```bash
export TURBOSPEC_DIR=/path/to/Turbospectrum_NLTE
export PATH=$TURBOSPEC_DIR/exec-gf:$PATH              # puts babsma_lu and bsyn_lu on your PATH
```

turbofit runs Turbospectrum through small `csh` scripts, so `csh` or `tcsh` must be installed. Most clusters have it; on your own Linux machine, run `sudo apt install tcsh` or the equivalent. turbofit also reads Turbospectrum's own data files (molecular data, partition functions) from `$TURBOSPEC_DIR/DATA`.

### 3. Model atmospheres

A model atmosphere describes the temperature and pressure through the star's outer layers. Turbospectrum needs one for every spectrum it computes. turbofit uses two MARCS grids.

#### 3a. APOGEE MARCS grid (for CNO-LTE mode, grids and FERRE)

This is the grid APOGEE's pipeline uses. It has grid points in Teff, log g, [M/H], [α/M] **and [C/M]**, so the atmosphere's carbon content follows the star's. It is on the SDSS DR17 archive:
<https://data.sdss.org/sas/dr17/apogee/spectro/speclib/atmos/marcs/MARCS_v3_2016/>

The archive blocks recursive `wget -r`, so list the directory and download it file by file:

```bash
URL=https://data.sdss.org/sas/dr17/apogee/spectro/speclib/atmos/marcs/MARCS_v3_2016/
OUT=/path/to/marcs_apogee
mkdir -p $OUT/raw && cd $OUT/raw
curl -s $URL | grep -o 'href="[^"?]*"' | cut -d'"' -f2      # first, look at what is there

# download everything in the directory and one level of sub-directories
for d in "" $(curl -s $URL | grep -o 'href="[^"?/]*/"' | cut -d'"' -f2); do
  for f in $(curl -s $URL$d | grep -o 'href="[^"?/]*\.[a-z0-9.]*"' | cut -d'"' -f2); do
    wget -nc -q -P "./$d" "$URL$d$f"
  done
done

# unpack into ONE flat directory of *.mod files (what turbofit expects)
cd $OUT
find raw \( -name '*.tar.gz' -o -name '*.tgz' -o -name '*.tar' \) -exec tar -xf {} -C $OUT \;
find raw -name '*.mod.gz' -exec sh -c 'gunzip -c "$1" > "'$OUT'/$(basename "${1%.gz}")"' _ {} \;
find $OUT -mindepth 2 -name '*.mod' -not -path '*/raw/*' -exec mv -n {} $OUT \;
ls $OUT/*.mod | wc -l                                         # number of models found
```

turbofit reads it with:

```python
lookup = turbogrid.make_atmod_lookup('/path/to/marcs_apogee', flat=True)
atm = lookup(teff, logg, mh, am, cm)        # the model file for one grid point
```

Converted copies are cached in `/path/to/marcs_apogee/_turbo_cache`.

**Holes in the grid.** Some grid points have no model, usually because the MARCS calculation did not converge. You have two options:
* `marcs_fill.py` makes stand-in atmospheres for the missing points, by interpolating between neighbouring models.
* `fill_grid_holes.py` instead fills the missing *spectra* in a finished spectral grid.

See the user guide for both.

#### 3b. Standard MARCS grid (for NLTE mode only)

The NLTE corrections were computed on the standard MARCS models, so NLTE mode uses those. Download `marcs_standard_comp.zip` from the folder of NLTE data shared by the Turbospectrum/TSFitPy developers, and unpack the `*.mod` files into one directory:

```bash
mkdir -p /path/to/marcs_standard && cd /path/to/marcs_standard
wget -O marcs_standard_comp.zip \
  'https://keeper.mpdl.mpg.de/d/6eaecbf95b88448f98a4/files/?p=%2Fatmospheres%2Fmarcs_standard_comp.zip&dl=1'
unzip -j marcs_standard_comp.zip '*.mod*'               # -j: no sub-folders
gunzip -f *.mod.gz 2>/dev/null; ls *.mod | wc -l
```

If the link has changed, browse the folder at <https://keeper.mpdl.mpg.de/d/6eaecbf95b88448f98a4/> (under `atmospheres/`).

### 4. NLTE data (for NLTE mode only)

In NLTE mode, Turbospectrum needs, for each NLTE element:

* a **departure-coefficient grid**: pre-computed NLTE corrections to the level populations, tabulated over Teff, log g, [Fe/H] and abundance. These are the `NLTEgrid*_MARCS_*.bin` and `auxData*_MARCS_*.dat` files.
* a **model atom**: the element's energy levels and transitions, the `atom.*` file.

**Step 1: compile the interpolator** that comes with Turbospectrum. It interpolates the NLTE grid to each star's parameters.

```bash
cd $TURBOSPEC_DIR/interpolator
gfortran -o interpol_modeles_nlte interpol_modeles_nlte_gfort.f    # Intel: interpol_modeles_nlte.f
```

**Step 2: download** the files for the elements you need from <https://keeper.mpdl.mpg.de/d/6eaecbf95b88448f98a4/>. Take the grids computed on **MARCS** models, not 1D/3D STAGGER. The grids are large, so download only what you need. TSFitPy keeps direct links to every file in [`utilities/nlte_grids_links.cfg`](https://github.com/TSFitPy-developers/TSFitPy/tree/main/utilities). For example, for Na:

```bash
NLTE=/path/to/nlte_data
mkdir -p $NLTE/Na $NLTE/model_atoms
wget -O $NLTE/Na/NLTEgrid4TS_Na_MARCS_Jul-14-2023.bin.zip 'https://keeper.mpdl.mpg.de/f/4223b6468d914975896b/?dl=1'
wget -O $NLTE/Na/auxData_Na_MARCS_Jul-14-2023.dat         'https://keeper.mpdl.mpg.de/f/f844a547fe0746c196e8/?dl=1'
wget -O $NLTE/model_atoms/atom.na_qmh                     'https://keeper.mpdl.mpg.de/f/e591b04ed6694306b876/?dl=1'
cd $NLTE/Na && unzip NLTEgrid4TS_Na_MARCS_Jul-14-2023.bin.zip && rm NLTEgrid4TS_Na_MARCS_Jul-14-2023.bin.zip
```

**Step 3: lay the files out** like this. The file names turbofit looks for are listed in `nlte_filenames.cfg`.

```
/path/to/nlte_data/
  Na/NLTEgrid4TS_Na_MARCS_Jul-14-2023.bin
  Na/auxData_Na_MARCS_Jul-14-2023.dat
  Mg/...
  model_atoms/atom.na_qmh
  model_atoms/atom.mg86b
```

**Step 4: check that turbofit can see them.**

```python
import nlte
cfg = nlte.NLTEConfig(marcs_dir='/path/to/marcs_standard', nlte_dir='/path/to/nlte_data')
print(cfg.available_elements())       # e.g. ['Na']
```

**NLTE line lists.** Turbospectrum only treats a line in NLTE if the line list says which energy levels of the model atom it connects (its "level labels"). Lines without labels are computed in LTE, silently. `aplinelists/turbospec.20180901t20.atoms_naNLTE` is the APOGEE atomic list with Na labels added. For other elements or lists, [TSFitPy](https://github.com/TSFitPy-developers/TSFitPy)'s `utilities/convert_lte_to_nlte.py` adds the labels. For the optical, the Gaia-ESO NLTE list in the keeper folder already has them.

### 5. Line lists

A line list gives the wavelength, strength and other data for every atomic and molecular line Turbospectrum should include.

* **APOGEE (H band, 15100–17000 Å):** the APOGEE DR17 list (Smith et al. 2021) is included in `aplinelists/`: atoms, molecules and hydrogen lines. Its wavelengths are in **air**, while APOGEE spectra are in **vacuum**. Convert one to the other before comparing (`turbogrid.write_grid_table(..., air_to_vac=True)` does this for grids).
* **Other wavelengths:** use any Turbospectrum-format list, e.g. from [VALD](https://vald.astro.uu.se) converted with Turbospectrum's tools.

### 6. FERRE (for grid mode only)

FERRE fits observed spectra by interpolating in a pre-computed grid of model spectra.

```bash
git clone https://github.com/callendeprieto/ferre /path/to/ferre
cd /path/to/ferre/src && make          # makes ../bin/ferre.x and ../bin/ascii2bin
```

Give turbofit the path to the executable as `ferre_exe='/path/to/ferre/bin/ferre.x'`.

### 7. APOGEE data and LSFs (for APOGEE spectra only)

* **LSF files.** APOGEE's line-spread function differs between fibres and changes with time, so each observation has its own calibration files. Download the `apLSF-*`/`apWave-*` files (APO, northern telescope) or `asLSF-*`/`asWave-*` files (LCO, southern telescope) you need from <https://data.sdss.org/sas/dr17/apogee/spectro/redux/dr17/cal/lsf/>. Put them all in one directory and set
  ```bash
  export APOGEE_LSF_DIR=/path/to/apogee_lsf
  ```
* **Spectra and catalogue.** The allStar catalogue (`allStar-dr17-synspec_rev1.fits`) and the combined spectra (`apStar-*`/`asStar-*` files; the file name of each star is in the allStar `FILE` column) are under <https://data.sdss.org/sas/dr17/apogee/spectro/>. `examples/apogee_gc_test.py` shows how to download and read them.

### 8. Check that it all works

```python
import synthesis, turbogrid
lookup = turbogrid.make_atmod_lookup('/path/to/marcs_apogee', flat=True)
atm = lookup(4500., 2.0, 0.0, 0.0, 0.0)
ll = ['/path/to/turbofit/aplinelists/turbospec.20180901t20.atoms',
      '/path/to/turbofit/aplinelists/turbospec.20180901t20.molec']
f, c, w = synthesis.synthesize(4500., 2.0, wrange=[15200, 15210], dw=0.01, atmod=atm, linelists=ll)
print(len(w), (f / c).min())       # about 1001 pixels, and a minimum below 1 (absorption lines)
```

If this fails, see [`docs/12-troubleshooting.md`](docs/12-troubleshooting.md).

### Summary: what goes where

| what | how turbofit finds it |
|---|---|
| turbofit code | `export PYTHONPATH=/path/to/turbofit:$PYTHONPATH` |
| Turbospectrum | `export TURBOSPEC_DIR=...`, plus `$TURBOSPEC_DIR/exec-gf` on your `PATH` |
| APOGEE MARCS grid | path given to `turbogrid.make_atmod_lookup(...)` |
| standard MARCS grid | `nlte.NLTEConfig(marcs_dir=...)` |
| NLTE grids and model atoms | `nlte.NLTEConfig(nlte_dir=...)` |
| line lists | list of file paths given as `linelists=` |
| FERRE | `ferre_exe='/path/to/ferre/bin/ferre.x'` |
| APOGEE LSFs | `export APOGEE_LSF_DIR=...` |

## Quick start

```python
import synthesis, broaden, turbogrid, nlte

# LTE spectrum on the APOGEE MARCS grid, at one grid point
lookup = turbogrid.make_atmod_lookup('/path/to/marcs_apogee', flat=True)
atm = lookup(4500., 2.0, -0.5, 0.25, -0.25)               # Teff, logg, [M/H], [a/M], [C/M]
flux, cont, wave = synthesis.synthesize(4500., 2.0, mh=-0.5, am=0.25, cm=-0.25, nm=0.5,
                                        vmicro=1.8, wrange=[15150, 15250], dw=0.01,
                                        atmod=atm, linelists=ll)
# normalise, broaden to R = 22500, resample to APOGEE-like pixels and shift by -45 km/s
w, f = broaden.convolve(wave, flux / cont, R=22500, dlogw_out=6e-6, rv=-45.)

# NLTE spectrum, Mg and Ca in NLTE
cfg = nlte.NLTEConfig('/path/to/marcs_standard', '/path/to/nlte_data')
fn, cn, wn = synthesis.synthesize(4500., 2.0, mh=-0.5, vmicro=1.8, wrange=[5500, 5600],
                                  dw=0.01, linelists=ll, nlte=cfg, nlte_elems=['Mg', 'Ca'])
```

Line-by-line fitting is shown in `examples/example_linefit.py`, and the full chain in `examples/example_pipeline.py`.

## Files

All modules sit in one folder and are imported directly (see Installation).

| file | what it does |
|---|---|
| `synthesis.py` | `synthesize()`: one spectrum, in LTE on any atmosphere, or in NLTE (`nlte=`, `nlte_elems=`); optionally only in chosen windows (`segments=`) |
| `nlte.py` | NLTE set-up: finds the standard MARCS models around a star and interpolates the NLTE grids |
| `abund.py` | `synthesize_sets()` / `ladder()`: many abundance patterns for one star |
| `broaden.py` | macroturbulence/rotation, instrument profile (R, R(λ) or LSF), Doppler shift, resampling |
| `turbogrid.py`, `fill_grid_holes.py` | build a grid of spectra → collect → fill holes → write a FERRE grid |
| `marcs_fill.py` | stand-in atmospheres for missing points of a MARCS grid |
| `continuum.py` | pseudo-continuum normalisation |
| `starfit.py` | **local mode**: full-spectrum stellar parameters with on-demand, cached syntheses |
| `turbofit.py` | `analyse()`: parameters (local or grid mode) → line masks → abundances for one spectrum |
| `ferre.py` | runs FERRE: puts spectra on the grid's pixels, writes the inputs, reads the results |
| `pipeline.py` | the grid-mode chain: RV → FERRE → masks → line-by-line |
| `linefit.py` | line-by-line fitting, automatic masks, errors, plots |
| `rv.py` | radial velocities by cross-correlation, with errors (Zucker 2003) |
| `apogee_lsf.py` | the APOGEE LSF for any fibre or combination of fibres |
| `atmos.py`, `atomic.py`, `utils.py` | atmosphere-file tools, solar abundances, helpers |
| `aplinelists/` | the APOGEE DR17 line list (air wavelengths), incl. a Na NLTE-labelled atomic list |
| `examples/` | runnable templates |
| `docs/` | the user guide |

## Things worth knowing

* **Turbospectrum's opacities are recomputed for every abundance set by default.** Turbospectrum runs in two steps: `babsma` computes the continuous opacity and electron pressure, then `bsyn` computes the lines. Reusing one `babsma` output for a different abundance mixture can be badly wrong for elements that supply electrons: ±0.3 dex in Fe changed the normalised flux by up to 0.09 in a 4500 K giant. `rerun_babsma=False` is available for trace elements, where the error is 0.001–0.003.
* **Fitting only windows.** Line fitting synthesises only the windows it fits, using Turbospectrum's `SEGMENTSFILE` option at its default sampling (R = 500 000). turbofit does not write the `RESOLUTION` keyword, because v20.1 misreads it.
* **Blended lines.** `linefit.fit(iterate=1)` refits after putting the measured abundances back into the synthesis, so lines of fitted elements that blend with each other are modelled at the right abundance. Known abundances of other elements can be given in `star['elems']`.
* **Macroturbulence** defaults to `'global'`: a first pass fits it per line, then it is fixed to a weighted median. A free per-line value is degenerate with the abundance for weak lines.
* **FERRE quirks.** FERRE exits with `STOP 1` (an error status) even after a successful run, so `ferre.py` checks for the output file instead. With the current FERRE, `COVPRINT=1` writes the full covariance matrix, not the upper triangle described in the manual. FERRE rejects starting points exactly on a grid edge, so turbofit moves them 0.001 inside.
* **FERRE grid order.** In grids written by `write_ferre_grid`, the last parameter varies fastest. Grids on a log-λ scale are written with `LOGW = 1`.

## Status

Tested against Turbospectrum v20.1 (gfortran) and FERRE (GitHub version):

* LTE synthesis, windows, abundance sets and ladders, ¹²C/¹³C.
* MARCS model selection (including holes), the NLTE input files and departure-file handling.
* Broadening and LSFs: line widths, equivalent-width conservation, R(λ), LSF arrays, Doppler shifts.
* Radial velocities and their errors.
* Grid building → FERRE grid → FERRE interpolation and fits, including microturbulence and broadening axes.
* Local-mode stellar parameters, line-by-line abundances and their errors, by injecting and recovering known values.

**Not yet validated:**

* NLTE synthesis end to end with real NLTE grids. Compare with TSFitPy before relying on it.
* Teff and log g in local mode with a full atmosphere grid.
* The APOGEE LSF from the real `apLSF`/`apWave` files, against `apogee.spec.lsf`.

## Credits and citations

If you use turbofit, please cite the codes and data it relies on:

* **Turbospectrum:** Alvarez & Plez (1998, A&A 330, 1109); Plez (2012, ASCL ascl:1205.004). **NLTE version (v20):** Gerber et al. (2023, A&A 669, A43). Also cite the papers for each NLTE grid/model atom you use, listed in Turbospectrum's `public-ref-NLTE-TS.bib`.
* **MARCS model atmospheres:** Gustafsson et al. (2008, A&A 486, 951). **APOGEE MARCS grid:** Mészáros et al. (2012, AJ 144, 120).
* **FERRE:** Allende Prieto et al. (2006, ApJ 636, 804).
* **APOGEE line list:** Smith et al. (2021, AJ 161, 254). If you use APOGEE data: Majewski et al. (2017, AJ 154, 94) and Abdurro'uf et al. (2022, ApJS 259, 35).
* **Radial-velocity errors:** Zucker (2003, MNRAS 342, 1291).

Code this builds on:

* The synthesis wrapper is adapted from **David Nidever**'s Turbospectrum Python wrapper, [dnidever/turbospectrum](https://github.com/dnidever/turbospectrum) (GPL-3.0). `atmos.py` is by **Jon Holtzman** (BSD-3).
* The NLTE inputs, species file and departure-file clean-up follow **TSFitPy** (Gerber et al. 2023; Storm & Bergemann 2023), GPL-3.0.
* The APOGEE LSF code is adapted from **Jo Bovy**'s `apogee.spec.lsf` (BSD-3), based on APOGEE pipeline code by J. Holtzman and D. Nidever.

## Licence

GPL-3.0. See [`LICENSE`](LICENSE).