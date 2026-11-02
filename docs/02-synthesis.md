# 2. Synthesising spectra

## 2.1 One spectrum: `synthesis.synthesize()`

```python
import synthesis, turbogrid
lookup = turbogrid.make_atmod_lookup('/data/marcs_apogee', flat=True)
atm = lookup(4500., 2.0, -0.5, 0.25, -0.25)          # teff, logg, [M/H], [a/M], [C/M]

flux, cont, wave = synthesis.synthesize(
    4500., 2.0, mh=-0.5, am=0.25, cm=-0.25, nm=0.5, vmicro=1.8,
    elems=[['Al', -0.2]],                 # optional [X/M] overrides
    wrange=[15150., 15250.], dw=0.01,
    atmod=atm, linelists=['atomic.list', 'molecular.list', 'Hlinedata'])
norm = flux / cont
```

It returns three arrays on Turbospectrum's wavelength grid:

* `flux`: absolute flux;
* `cont`: continuum;
* `wave`: wavelength in Å.

### What each argument controls

| argument | meaning |
|---|---|
| `teff, logg` | used for the geometry choice and the default isotope ratios (the atmosphere itself comes from `atmod`) |
| `mh, am, cm, nm` | abundances: [M/H] scales every element heavier than He; [α/M] adds to O, Ne, Mg, Si, S, Ar, Ca, Ti; [C/M] and [N/M] add to C and N |
| `elems` | `[[element, [X/M]], ...]` sets these elements to exactly that [X/M] (overriding `am`/`cm`/`nm`) |
| `vmicro` | microturbulence in km/s, passed to Turbospectrum as a fixed value (`XIFIX`) |
| `wrange, dw` | wavelength range and step (Å) |
| `atmod, atmos_type` | atmosphere file; type auto-detected (`'marcs'`, `'kurucz'`, …). For Turbospectrum-format interpolated models pass any type other than `'marcs'`, e.g. `'interpol'` |
| `dospherical` | spherical radiative transfer for MARCS models with log g ≤ 3 (default `True`) |
| `solarisotopes` | `None` (default) picks from the model name: `p…` plane-parallel → solar ratios, `s…` spherical → giant ratios (¹²C/¹³C = 15); otherwise from log g |
| `c12c13` | carbon isotope ratio ¹²C/¹³C (number ratio), e.g. `c12c13=8.`; overrides `solarisotopes` for carbon |
| `segments` | list of `(wmin, wmax)`: compute only these windows (see below) |
| `nlte, nlte_elems` | NLTE mode, see [chapter 3](03-nlte.md) |
| `workdir, save` | where Turbospectrum runs (default: a temporary directory, deleted afterwards unless `save=True`) |

The abundance array Turbospectrum receives is built by `synthesis.make_abundances(mh, am, cm, nm, elems)`. Call it yourself to see exactly which A(X) are used.

### Only some windows: `segments`

When you only need a few windows, synthesising the whole range is wasteful:

```python
f, c, w = synthesis.synthesize(..., wrange=[15100, 16900],
                               segments=[(15739.5, 15742.5), (15765.0, 15768.5)])
```

Turbospectrum then computes just those windows, using its `SEGMENTSFILE`. The returned arrays hold the windows one after another.

> [!NOTE]
> Within segments, Turbospectrum samples at R = 500 000, its default, and `dw` is ignored there. turbofit does not set the `RESOLUTION` keyword because Turbospectrum v20.1 reads it incorrectly. 500 000 is finer than any instrument, and everything downstream resamples anyway.

## 2.2 Many abundance patterns for one star: `abund.py`

`abund.synthesize_sets()` uses one atmosphere and returns one spectrum per abundance pattern:

```python
import abund
wave, fluxes, conts = abund.synthesize_sets(
    4500., 2.0, mh=-0.5, am=0.25, cm=-0.25, nm=0.5, vmicro=1.8,
    abund_sets=[[],                                   # the base star
                [['Mg', 0.5]],
                [['Mg', 0.5], ['Al', -0.3]]],
    wrange=[15150., 15250.], dw=0.01, atmod=atm, linelists=ll)
# fluxes.shape == (3, npix)
```

Each entry of `abund_sets` is applied on top of `elems`.

`abund.ladder()` is the common special case: one element, stepped around its starting value.

```python
xm, wave, fluxes, conts = abund.ladder(4500., 2.0, 'Mg', offsets=[-0.4, -0.2, 0, 0.2, 0.4],
                                       mh=-0.5, am=0.25, wrange=[15150., 15250.], dw=0.01,
                                       atmod=atm, linelists=ll)
# xm = [X/M] values actually used: base [Mg/M] (= am here) + offsets
```

`abund.base_xm(elem, am, cm, nm, elems)` gives the starting [X/M] of an element: `am` for α elements, `cm` for C, `nm` for N, the `elems` value if one is given, and 0 otherwise.

### Different ¹²C/¹³C ratios

`abund.isotope_ladder()` gives the same star at several carbon isotope ratios, with everything else fixed:

```python
ratios, wave, fluxes, conts = abund.isotope_ladder(
    4500., 2.0, ratios=[4, 6, 8, 10, 15, 20, 30, 50, 90],
    mh=-0.5, am=0.25, cm=-0.25, nm=0.5, vmicro=1.8,
    wrange=[15150., 15250.], dw=0.01, atmod=atm, linelists=ll)
```

Any abundance set can also carry a ratio. `synthesize_sets(..., c12c13=8.)` uses one ratio for every set, and `c12c13=[...]` gives one ratio per set.

The ratio is written to Turbospectrum's `ISOTOPES` block as the fractions ¹²C = r/(1+r) and ¹³C = 1/(1+r), for both babsma and bsyn. Without `c12c13`, carbon follows `solarisotopes`: Turbospectrum's solar mix (¹²C/¹³C ≈ 89.9) or the "giant" value of 15. Only molecular lines labelled with explicit isotopes (e.g. `'0607.013014'` for ¹³C¹⁴N, `'0608.013016'` for ¹³C¹⁶O) respond to the ratio, which is how standard Turbospectrum molecular line lists are labelled.

> [!WARNING]
> **Continuum opacities are recomputed for every set by default (`rerun_babsma=True`).** bsyn takes the electron pressure from babsma's opacity file, so reusing one file across abundance sets freezes the ionisation balance. This matters most when the varied element is an electron donor (Fe, Mg, Si, Na, Al, Ca). Set `rerun_babsma=False` only for trace elements and small steps; it is somewhat faster.

NLTE works here too. If the varied element is in `nlte_elems`, its departure coefficients are re-interpolated for each abundance step.

