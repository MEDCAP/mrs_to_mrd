# `mrd2recon.py` vs `mrd2_recon_to_incorporate.py` — recon logic comparison

Side-by-side check of whether the EPSI recon logic in `mrd2_recon_to_incorporate.py::epsi_recon`
has been adopted into `mrd2recon.py`, and where it was deliberately changed.

| Step | Legacy (`mrd2_recon_to_incorporate.py`) | New (`mrd2recon.py`) | Status |
|---|---|---|---|
| Read acquisitions | `raw_acquisition_list`, sorted by timestamp, split later by `average` flag | Streamed once, routed by `encoding_space_ref` into `buffers` dict | Adopted, improved (structural routing vs timestamp sort + average heuristic) |
| Geometry | hardcoded `necho=64, nro=12` (comment admits it overrides the real read) | read per-acquisition: `nswitches` from header, `kept`/`discard_pre`/`discard_post` from acq head, checked via `EncodingBuffer.check` | **Deliberate departure**, documented in module docstring line 49 |
| Line broadening / k-space fill | `for iecho in range(64): kspace[...] = a.data[0, iecho*totalppswitch+discard_pre : +nro]; *= exp(-tk*lb)` | `apply_line_broadening`: same reshape-and-slice logic, generalized to `nswitches`/`npoints_per_switch`, `tk` computed via `np.arange` | Adopted, same math, vectorized |
| Separating series vs prescan | `average > 1` → phantom; `average == 1` → series (fragile convention) | `IS_NOISE_MEASUREMENT` flag on acquisition selects `is_prescan` | Adopted, made explicit/robust (flag vs magic average count) |
| FFT | `np.fft.fftshift(np.fft.fftn(kspace))` per image, phantom and series separately | `fft_kspace_to_image`, same fftn+fftshift, run per repetition for each buffer | Adopted |
| Find brightest voxel | triple nested loop `for ide/j/k: max(abs(...))` | `np.abs(img).max(axis=-1)` then `np.argmax`/`unravel_index` | Adopted, vectorized, same semantics (first max in C order) |
| Noise floor | `noise = mean(abs(imgset[-1]))`, only for non-phantom; phantom threshold hardcoded 0 | Same: `noise*NOISE_THRESHOLD_MULTIPLIER` for series, `0.0` for prescan, plus a warning added for an all-zero last rep | Adopted, warning added |
| Phase/roll alignment | roll search `range(-15,16)`, `S0`,`Spi2` overlap, `arctan2` phase correction, applied per voxel, `globalspect` accumulated | `phase_align`: identical math (`PHASE_SEARCH_RANGE=15` default), same overlap formula, same accumulation, now returns aligned array + global spectrum instead of mutating module globals | Adopted 1:1 |
| Spectral axis (`xscale`)/BW | `BW = 1/sampletime/totalppswitch`; `xscale = arange(n)/n * BW/centerfreq*1e6` | Same formula in `make_buffer`: `spectral_bw_hz = 1e9/(sample_time_ns*totalppswitch)`, `xscale = arange(nfreq)/nfreq*bw_ppm`, with `nfreq = nswitches*FIDPAD` | Adopted, extended with `FIDPAD` zero-fill (documented rationale, not in legacy) |
| Width guess (FWHM) | manual half-max crossing search on `globalspect` around its argmax | `estimate_width_fwhm` (in `lorentzian_fitter.py`) — same FWHM concept | Adopted via helper module |
| Phantom global fit | single peak at `argmax`, `lornfit` with bounds `[widthguess/2, widthguess*1.5]` only on width, everything else unbounded | `PHANTOM_SPEC` (one peak, offset 0) run through same `fit_global_multipeak`/`LorentzianFitter` as series, pinned same as series voxels | Adopted, **deliberately changed**: phantom voxels now pinned like series (docstring lines 53–55), legacy fitted unbounded |
| Phantom scaling / `phantomscaling` | `A[0]*w[0]*globalspectscaling`, averaged over phantom images, used nowhere downstream except saved to `.mat` | `phantom_global_peak_areas` computed the same way (amplitude×width×scaling) via `fit_and_emit_peaks`, emitted as NdArray, "recorded, never applied" | Adopted, same semantics, output channel changed (mrd NdArray vs `.mat`) |
| Restoring module state after phantom loop | explicit `lornputspect(xscale, hpglobalspect, ...)` restore call needed because `lorn` module globals were clobbered by phantom voxel loop | N/A — `EncodingBuffer`/`LorentzianFitter` are per-encoding objects, no shared mutable state | **Correctly not ported** — this restore existed only because of the legacy's global-state bug (see module docstring lines 349–353) |
| Global multi-peak fit (hypothesis search) | loop over `biggestpeaklist`, try each as "the tallest peak", `centers = (argmax - (offsets-offsets[icg])) % BW`, `minimize(lornfit, x0, bounds)`, pick `argmin(diff)` | `fit_global_multipeak`: loop over `spec.biggest_idx`, `candidate_centers(...)`, `fitter.fit_global(...)`, pick lowest `params.loss` | Adopted 1:1 in structure |
| Global fit bounds | width bound `[widthguess/2, widthguess*1.5]` i.e. (0.5,1.5) scale; **no center bound at all** | `WIDTH_BOUND_SCALE=(0.1,1.9)` (restores original arctan-implied range), `GLOBAL_CENTER_WINDOW_PPM=(0.1,1.0)` (new, absent from legacy) | **Deliberate correction**, extensively justified in comments lines 94–119 — this is the "correction of syntax and field" bug fix, not a straight port |
| Per-voxel fit | `for ide in range(4, hpimgset.shape[0])` — **skips first 4 repetitions** ("debug shortcut") | `for ide in range(nreps)` — all repetitions | **Deliberate correction**, documented line 52 |
| Per-voxel fit bounds | `bnds[ip] = [-fit_df, fit_df]` for **center itself** (not an offset from global fit — looks like a bug, center bound is absolute not relative to `centers[ip]`), width `[widths[ip]-fit_dw, widths[ip]+fit_dw]`, phase `[phases[ip]-fit_dph, phases[ip]+fit_dph]`, amplitude `[0, None]` | `fitter.fit_windowed(spect, df=fit_df, dw=fit_dw, dph=fit_dph)` — windows applied relative to the global fit's centers/widths/phases | Adopted with the apparent legacy bug (`bnds[ip]=[-fit_df,fit_df]` ignoring `centers[ip]`) fixed |
| Voxel amplitude/area output | `metabolites[:, ide,j,k] = abs(A)*scaling`; `metabolites2 = abs(A*w)*scaling` | `amplitudes[:, ide,j,k] = abs(params.amplitudes)`; `areas = abs(params.amplitudes*params.widths)`, scaled by `global_scaling` in `fit_and_emit_peaks` | Adopted 1:1 |
| Voxel starting amplitude | `x0` amplitudes start at 0 implicitly (`x0 = np.zeros(...)`) | starts from the spectrum's own values | **Deliberate departure**, documented line 56–57 |
| `1puls`/`fid` path (`spectra_recon`, kinetics `kABfit`) | present, separate branch keyed on sequence name | **absent** — module docstring says "EPSI only: anything else raises" | **Not ported** — explicit scope narrowing, stated up front |
| Image/plot output | `matplotlib` plots, `.mat` save, mrd `Image` stream items (`generate_epsi_images`, `generate_aux_images`) | mrd `NdArray` stream items only (`emit`), no images, no plotting, no `.mat` | **Not ported** — deliberate output-format change (module docstring lines 44–45) |
| Header record of run params | none | `append_recon_header` records `line_broadening`, fit windows, global center window, width bounds, and the peak list itself | **New**, not in legacy |
| Raw acquisitions passthrough | `writer.write_data(generate_stream(raw_acquisition_list))` | `writer.write_data(mrd.StreamItem.Acquisition(a) for a in acquisitions)` | Adopted 1:1 |
| CLI peak-arg parsing (`-name_mods value`) | manual loop over `sys.argv`, same `_s/_t/_m` modifier convention | `split_peak_args` — same convention, pulled out before argparse for negative-value safety | Adopted, same semantics |

## Net assessment

Every piece of load-bearing EPSI logic from `epsi_recon` is present in `mrd2recon.py`: k-space
assembly with line broadening, FFT, brightest-voxel reference, roll/phase alignment, FWHM width
guess, global multi-hypothesis Lorentzian fit, per-voxel fit, phantom handling. The differences
are all things the new module's docstring already calls out as intentional fixes, not omissions:

- hardcoded `necho=64,nro=12` → read from the acquisition header
- first FID point discarded → kept
- reps 0–3 skipped → all reps fitted
- phantom voxels fitted unbounded → now pinned like the series
- global fit center left totally unbounded (the bug that let a tiny peak slide onto a
  neighbor) → `GLOBAL_CENTER_WINDOW_PPM` added
- width bound narrowed to (0.5,1.5) by the legacy port → restored to (0.1,1.9), matching the
  original arctan-based transform

Not ported, and correctly so given the new module's stated scope: the `1puls`/`fid` spectra path
and kinetic (`kABfit`) modeling, all matplotlib/`.mat` output, and the module-global `lorn` state
(and its explicit "restore after phantom" workaround, which existed only to patch around that
same global state).

## Open question

The legacy's per-voxel center bound `bnds[ip] = [-fit_df, fit_df]` looks like it ignores
`centers[ip]` entirely (an absolute bound near zero rather than a window around the fitted
center) — if that was intentional rather than a bug, the new `fit_windowed(df=fit_df, ...)`
semantics (window relative to global center) would be a behavior change, not just a syntax
correction. Worth verifying against `lorn_to_incorporate.py` if it matters.
