# 12. Troubleshooting and known quirks

## 12.1 Errors you may see

| symptom | likely cause | what to do |
|---|---|---|
| `Turbospectrum failed` with `babsma_lu: Command not found` | Turbospectrum executables not on `PATH` | `export PATH=$TURBOSPEC_DIR/exec-gf:$PATH` |
| `/bin/csh: not found` or the scripts do nothing | no csh | install `tcsh` |
| `FileNotFoundError` for a line list or atmosphere | wrong path | check `linelists`, `atmod`, `griddir`, and `flat=` in `make_atmod_lookup` |
| bsyn stops with `NDP too small` or `too many wavelengths` | model has more layers than compiled for, or the range/step is too large | split the range, or raise `NDP`/`lpoint` in `source/spectrum.inc` and recompile |
| `no atmosphere for corner(s) ...` in `linefit` | a node of the cell around the star is missing | fill holes with `atmos.fill_holes()`, or check the star is inside the nodes |
| `... outside the atmosphere nodes ...` | the star lies beyond the node range | extend the nodes or drop the star |
| `none of the requested nlte_elems have NLTE data available` | bin/aux/model-atom files not found | `cfg.available_elements()`; check the layout against `nlte_filenames.cfg` |
| `interpol_modeles_nlte not found` | the NLTE interpolator is not compiled | see [1.4](01-installation.md) |
| `NLTE interpolation failed for X` | MARCS models missing around the star, or the abundance is outside the departure grid | check the MARCS directory; departure grids span about −2…+1 dex around solar-scaled |
| `FERRE failed, see .../ferre.log` | FERRE never wrote its output file | read the log: grid path, `INTER` larger than the axis sizes allow, or start values outside the grid |
| `LSF narrower than the model sampling` | model step too coarse for the LSF | synthesise on a finer `dw` |
| `wave_out does not fall on the APOGEE detectors` | output grid outside the APOGEE range | use the apStar grid, or wavelengths within ~15100–17000 Å (vacuum) |
| `fitting 12C/13C needs atmosphere nodes` | `C12C13` requested without `nodes`/`atmod_lookup` | isotope fits run in CNO-LTE mode; configure the nodes |
| `no atmosphere for node {...}` in `starfit` | a node of the current cell is missing | fill holes (`atmos.fill_holes`), fix that parameter, or start elsewhere |
| `starfit` flags `teff_at_grid_edge` (etc.) | the best fit wants to go beyond the node range | the parameter is outside your atmosphere grid; treat it as a limit |
| `starfit` flags `not_converged` | the cell kept moving for `max_iter` iterations | give a better `guess`, raise `max_iter`, or fix a poorly constrained parameter |
| `starfit` is slow | 288 syntheses per new cell over a wide range | use `segments=` to synthesise only what you fit, fix parameters you know, share a `cache_dir` between stars |
| grid-mode vsini/vmac biased high | linear resampling of the data | keep `resample='cubic'` (default) in `ferre.fit` |
| `rv_err = nan` | the CCF peak is at the edge of `vmin..vmax`, or the CCF is too flat | widen the range, or check normalisation and template |

## 12.2 Lines flagged in `linefit`

* **`weak`**: a 0.1 dex change moves the window by less than 3σ (for `C12C13`, a 0.3 dex change in the ratio). The line is uninformative at this S/N. Raise `min_response` in `automask` or accept fewer lines.
* **`edge`**: the fit sits at the end of the abundance ladder even after extension. Usually the window is badly blended, has the wrong continuum, or contains a defect. Look at the plot.
* **`poorfit`**: χ²_red > 5. Typical causes are an unmodelled blend, a bad pixel, or a line-list problem. Check `drv` too: large offsets point to misidentified or blended lines.

## 12.3 Known quirks (handled for you)

* **Turbospectrum v20.1's `RESOLUTION` keyword** is read with an `(a)` format and gives a garbage value. turbofit never writes it, so segments use Turbospectrum's default R = 500 000.
* **Reusing babsma opacities** across abundance changes freezes the electron pressure, which matters for electron donors such as Fe. babsma is rerun by default ([2.2](02-synthesis.md)).
* **FERRE ends with `STOP 1`**, a non-zero exit status, even after a successful run. Success is judged by the output file.
* **FERRE's `COVPRINT`** writes the full covariance matrix, not the upper triangle described in the manual.
* **FERRE rejects starting values exactly on a grid edge.** They are nudged 0.001 inside.
* **Departure files sometimes contain NaN.** These are replaced as described in [3.4](03-nlte.md).

## 12.4 Getting help

When reporting a problem, include:

* the call;
* the full traceback;
* for Turbospectrum failures, the run directory (`save=True, workdir=...` in `synthesize()` keeps it, including the `*_babsma.log` and `*_bsyn.log` files).
