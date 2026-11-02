# 1. Installation

turbofit is plain Python on top of three external pieces:

* **Turbospectrum**: always needed.
* **The NLTE data**: only for NLTE mode.
* **FERRE**: only for stellar-parameter fitting.

## 1.1 Turbospectrum v20

```bash
git clone https://github.com/bertrandplez/Turbospectrum_NLTE
cd Turbospectrum_NLTE/exec-gf && make        # gfortran; or exec/ with ifx/ifort
```

Then set two environment variables, for example in your `~/.bashrc`:

```bash
export TURBOSPEC_DIR=/path/to/Turbospectrum_NLTE
export PATH=$TURBOSPEC_DIR/exec-gf:$PATH      # puts babsma_lu and bsyn_lu on PATH
```

turbofit uses `$TURBOSPEC_DIR/DATA` for Turbospectrum's data files (molecular equilibrium, partition functions, H lines). If `TURBOSPEC_DIR` is not set, it falls back to a `data/` folder or symlink next to the code.

> [!IMPORTANT]
> Turbospectrum is driven through small **csh** scripts, so `csh` (or `tcsh`) must be installed. On clusters it almost always is; on a fresh Linux machine use `apt install tcsh` or similar.

Check the installation by synthesising one short spectrum, as in [chapter 2](02-synthesis.md), before doing anything large.

## 1.2 Python packages

numpy, scipy, astropy and matplotlib. Nothing else is required.

```bash
pip install numpy scipy astropy matplotlib
```

turbofit uses flat imports, so there is nothing to install:

```python
import sys
sys.path.insert(0, '/path/to/turbofit')
import synthesis, broaden, linefit
```

## 1.3 Model atmospheres

| used for | grid | how turbofit reads it |
|---|---|---|
| CNO-LTE mode, grids, FERRE | APOGEE MARCS grid with [C/M] and [α/M] nodes | `turbogrid.make_atmod_lookup(griddir, model='marcs', flat=True)` returns a function `lookup(teff, logg, mh, am, cm)` giving the Turbospectrum-ready file for a node |
| NLTE mode | standard MARCS grid (ASCII web format, `*.mod`, all in one directory) | `nlte.NLTEConfig(marcs_dir=..., ...)` indexes it automatically |

`make_atmod_lookup` converts each needed model once and caches the converted copy in `griddir/_turbo_cache`. It also uses the `*.filled` stand-ins written by `atmos.fill_holes()` for missing nodes.

## 1.4 NLTE data (only for NLTE mode)

**Step 1: compile the NLTE interpolator** that ships with Turbospectrum.

```bash
cd $TURBOSPEC_DIR/interpolator
gfortran -o interpol_modeles_nlte interpol_modeles_nlte_gfort.f
# (Intel compilers: interpol_modeles_nlte.f)
```

**Step 2: download** from <https://keeper.mpdl.mpg.de/d/6eaecbf95b88448f98a4/>, for each element you want:

* the departure-coefficient grid computed on **MARCS** models: `NLTEgrid*MARCS*.bin` and `auxData*MARCS*.dat`;
* the model atom, `atom.*`.

The grids are large, so download only the elements you need.

**Step 3: lay the files out** as listed in `nlte_filenames.cfg` (a copy ships with turbofit):

```
nlte_data/
  Mg/NLTEgrid4TS_Mg_MARCS_....bin
  Mg/auxData_Mg_MARCS_....dat
  Ca/...
  model_atoms/atom.mg86b
  model_atoms/atom.ca105b
  nlte_filenames.cfg          (optional; turbofit's copy is used otherwise)
```

**Step 4: check what turbofit can see.**

```python
import nlte
cfg = nlte.NLTEConfig(marcs_dir='/data/marcs_standard', nlte_dir='/data/nlte_data')
print(cfg.available_elements())     # e.g. ['Ca', 'Fe', 'Mg']
```

NLTE needs line lists whose lines carry **NLTE level labels**. The Gaia-ESO NLTE atomic list on the same keeper page has them (4200–9200 Å). VALD extractions do not, but TSFitPy's `utilities/` can add the labels. Lines without labels are simply computed in LTE.

## 1.5 APOGEE LSF files (only for the APOGEE LSF)

To use the real APOGEE line-spread function (`apogee_lsf.py`, chapter 4), put these six files in one directory:

* `apLSF-a-05440020.fits`, `apLSF-b-05440020.fits`, `apLSF-c-05440020.fits`;
* `apWave-a-02420038.fits`, `apWave-b-02420038.fits`, `apWave-c-02420038.fits`.

Then pass that directory to `ApogeeLSF(...)` or set `export APOGEE_LSF_DIR=/path/to/them`.

## 1.6 FERRE (only for stellar parameters)

```bash
git clone https://github.com/callendeprieto/ferre
cd ferre/src && make          # produces ../bin/ferre.x and ../bin/ascii2bin
```

Pass the path to `ferre.x` as `ferre_exe=` (see [chapter 9](09-ferre.md)).

## 1.7 Quick check

```python
import os, sys
sys.path.insert(0, '/path/to/turbofit')
import synthesis, turbogrid

lookup = turbogrid.make_atmod_lookup('/data/marcs_apogee', flat=True)
atm = lookup(4500., 2.0, 0.0, 0.0, 0.0)
f, c, w = synthesis.synthesize(4500., 2.0, wrange=[15200, 15210], dw=0.01,
                               atmod=atm, linelists=[...])
print(len(w), (f / c).min())          # ~1001 pixels and some absorption
```
