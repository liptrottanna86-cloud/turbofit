# 3. NLTE

## 3.1 When to use it

turbofit has two ways of synthesising an element:

| | **CNO-LTE** | **NLTE** |
|---|---|---|
| atmospheres | APOGEE MARCS grid with [C/M] and [α/M] nodes | standard MARCS grid (α tied to [Fe/H], solar C/N) |
| radiative transfer | LTE | NLTE for the chosen elements, LTE for everything else |
| between grid points | spectra are interpolated (never atmospheres) | Turbospectrum's interpolator builds the atmosphere *and* the departure coefficients at the exact parameters |
| best for | Teff, log g, [M/H] (where CNO atmospheres matter); C, N, O from molecules; elements without NLTE data | atomic abundances with strong NLTE effects |

The two are kept separate on purpose. The public departure-coefficient grids were computed on standard MARCS models, so applying them to the CNO atmospheres would be inconsistent. The usual chain is therefore:

1. stellar parameters from FERRE on a CNO grid (LTE);
2. abundances in NLTE wherever data exist, and CNO-LTE otherwise.

Every abundance turbofit reports records which mode was used.

NLTE data exist for H, O, Na, Mg, Al, Si, Ca, Ti, Cr, Mn, Fe, Co, Ni, Sr, Y, Ba and Eu. Molecules (CO, CN, OH, CH) are always in LTE in Turbospectrum.

## 3.2 Setting it up

```python
import nlte
cfg = nlte.NLTEConfig(marcs_dir='/data/marcs_standard',     # *.mod, ASCII web format
                      nlte_dir='/data/nlte_data')           # bins, aux files, model_atoms/
cfg.available_elements()                                     # what is actually on disk
```

Optional arguments:

* `config_file`: which `nlte_filenames.cfg` to use;
* `model_atom_dir`: default `nlte_dir/model_atoms`;
* `interpol_dir`: default `$TURBOSPEC_DIR/interpolator`;
* `comp`: MARCS composition class, default `'st'` (standard).

## 3.3 One NLTE spectrum

```python
f, c, w = synthesis.synthesize(4500., 2.0, mh=-0.5, vmicro=1.8,
                               wrange=[5500., 5600.], dw=0.01, linelists=ll,
                               nlte=cfg, nlte_elems=['Mg', 'Ca'])
```

`atmod` is not needed, and is ignored if given. `am` and `cm` still change the *synthesis* abundances, but not the atmosphere, which is standard MARCS.

If a requested element has no data, it is computed in LTE with a warning. If none of the requested elements has data, an error is raised.

## 3.4 What happens under the hood

For each call, `nlte.prepare()` does the following:

1. **Picks the model family.** For log g < 3 it uses spherical, 1 M☉, ξ = 2 km/s models. For log g ≥ 3 it uses plane-parallel ξ = 1 km/s models, or ξ = 2 if those are missing. This choice only concerns the atmosphere *structure*: the microturbulence used in the synthesis is still your `vmicro`.
2. **Finds the 8 MARCS models bracketing (Teff, log g, [M/H]).** If a corner model is missing, the face of the cube with the most missing corners is moved out by one grid step until all 8 exist, as TSFitPy does.
3. **Runs `interpol_modeles_nlte` once per NLTE element.** It interpolates the atmosphere and that element's departure coefficients with the same weights, at the element's A(X). It writes `*.interpol` (the model) and `*_<El>_coef.dat` (the departures).
4. **Cleans the departure files.** NaNs in the deepest layer become 1 (LTE); other NaNs take the value of the layer below.
5. **Writes the NLTE species file** that bsyn reads (`'NLTE :' '.true.'`, `'NLTEINFOFILE:'`).

The interpolator inputs, species file and NaN handling follow TSFitPy (Gerber et al. 2023; Storm & Bergemann 2023).

## 3.5 NLTE in grids and line fits

* **Grids:** `turbogrid.write_grid_table(..., nlte=cfg, nlte_elems=[...])` with `atmod_lookup=turbogrid.nlte_atmod_lookup`. See [chapter 5](05-grids.md).
* **Line by line:** give `LineFitter(nlte=cfg, nlte_elems=[...])`. Those elements are then fitted in NLTE when their data are present, and in CNO-LTE otherwise. See [chapter 8](08-line-by-line.md).

> [!CAUTION]
> The NLTE path was built and checked piece by piece (bracketing, NaN clean-up, bsyn reading an NLTE species file), but it has **not yet been run end to end with real departure grids**. Before relying on it, compare one NLTE spectrum with TSFitPy for the same star and line.
