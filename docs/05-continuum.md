# 5. Continuum normalisation

turbofit's fitting routines don't need you to normalise your spectra first, and that is deliberate:

* **Local-mode parameter fits** (chapter 8) multiply the model by a smooth polynomial per chunk of the spectrum, solved at every step together with the stellar parameters. Data and model are therefore treated identically, which avoids the usual bias of a pseudo-continuum placed differently on the two.
* **Line fits** (chapter 10) do the same in each window, with a low-order polynomial.

Normalisation is still needed in two places:

* for **grid mode**, because FERRE grids are continuum-normalised;
* for anything you want to look at, such as RV templates and quick-look plots.

`continuum.py` provides it, without instrument presets.

## 5.1 Automatic chunks

```python
import continuum
continuum.find_chunks(wave)          # [(w0, w1), ...] split at gaps in the data
```

A new chunk starts wherever the wavelength jumps by more than 20 median pixels (`gap=` to change it), or where the flux is not finite. Chip gaps (APOGEE), separate gratings or detectors (NIRSpec MOS), separate setups (FLAMES) and masked regions therefore become separate chunks, and each gets its own continuum.

## 5.2 Polynomial pseudo-continuum

```python
cont, nflux, nerr = continuum.normalise(wave, flux, err, deg=3)
cont, nflux, nerr = continuum.normalise(wave, flux, err, deg=4,
                                        chunks=[(15140., 15810.), (15850., 16440.), (16470., 16960.)])
```

For each chunk, a polynomial is fitted iteratively with **asymmetric clipping**: points more than `lo` = 1σ below or `hi` = 3σ above are rejected. This pushes the fit up onto the continuum and away from absorption lines. Choose the degree for your data:

* low (1–3) for short chunks or crowded spectra, such as cool giants or the H band;
* higher for long, clean chunks.

## 5.3 Running percentile

```python
cont = continuum.running_percentile(wave, flux, width=20., q=90)
```

This takes the q-th percentile in a running window of `width` Å. It is robust in very crowded regions, where a polynomial struggles.

## 5.4 How good is it?

On a simulated R = 22 500, S/N = 150 spectrum with 150 lines per 200 Å, a chip gap and a curved continuum:

* the median deviation from the true continuum was 0.18 %;
* the maximum deviation was 1.3 %.

Part of this is simply that the pseudo-continuum sits slightly below the true one, because of blended line wings. That offset is the reason the fitting routines solve for the continuum together with the model, rather than relying on this normalisation.
