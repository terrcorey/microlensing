# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A playground for analysing real microlensing light curves, to learn the
challenges involved in modeling brown dwarf populations with microlensing --
not an application anyone installs or upgrades. 

## Rules

1. You are never to edit the rules and instructions sections on your own. 
   The rest are safe to touch, but if you want to edit these sections
   file, the answer is always NO. Instead, draft up what needs to be added 
   or changed and send it to the user directly.
2. You are never to make any git commits and pushes. Read-only commands are
   okay. 
3. If you require access to files outside this `microlensing/` folder, you
   MUST ask for permission from the user first.
4. Before you write or edit any code, check if the ponytail skill is
   active. If not, confirm with the user that this is intended before
   proceeding.
5. Prioritize a clean, scalable, legible codebase: prefer one shared
   implementation over copy-pasted logic across scripts, keep one-time/
   dev-script code and output separated from the regular pipeline
   (`scratch/`, see "Output layout" below), and fold single-caller files
   into their one caller instead of leaving them as separate scripts.
6. Your main role in this repository is to act as a guiding role. Give an overview rundown of what needs to be done, give suitable hints and direction to the user but allow the user to write their own code. After they finish, you can simplify using ponytail and tidy up.


## Instructions

1. If the user says "get up to speed": read this file and CHANGELOG.md to
   understand the most recent changes and any stated short/long-term goals.
   If CHANGELOG.md has no goals recorded, use `/grill-me` to establish what
   the user wants to work on this session.
2. If the user says "save state": review the session's changes. If they
   changed the architecture/commands described in this file, you should 
   apply the changes below (but never rules and instructions). Log a summary to
   CHANGELOG.md as a new dated entry (see its own header for the format), then 
   run /graphify --update if major architectural changes have occured to ensure 
   the graph is consistent with current structure. Then suggest possible goals
   for next session and use `/grill-me` to confirm them with the user before 
   adding them to that entry.
3. Whenever the user describes wanting to do or build something that isn't
   already covered in CHANGELOG.md, use `/grill-me` to reach a consensus
   before taking any action.

## Setup and commands

```
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt


python3 download_data.py               # fetch raw photometry into data/ (skips files that already exist)
python3 preprocess_binary_data.py      # required before mcmc_fit_binary.py (see below); absorbs the old joint_fit.py

python3 mcmc_fit.py --stage=raw        # raw_lc/O-05-BLG086.png
python3 mcmc_fit.py --stage=quicklook  # + fit_lc/O-05-BLG086.png (curve_fit point estimate, no MCMC)
python3 mcmc_fit.py                    # default --stage=mcmc: + corner_plots/O-05-BLG086.png
                                        #                        hist_plots/O-05-BLG086.png
                                        #                        data/processed/O-05-BLG086_fit_summary.dat

python3 mcmc_fit_binary.py --stage=raw
python3 mcmc_fit_binary.py --stage=quicklook
python3 mcmc_fit_binary.py             # default --stage=mcmc: runs both fits below, in order:
                                        #   run_ogle_only_diagnostic() -> fit_lc/corner_plots *_ogle_only.png
                                        #   run_joint_fit()            -> fit_lc/O-03-BLG235.png + corner_plots/O-03-BLG235.png (canonical)
```

`mcmc_fit_binary.py` also takes `--dataset` (default `O-03-BLG235`), backed by a
`DATASETS` dict keyed by short name (raw OGLE filename + dataset-specific
`u0`/`tE` guesses for `run_ogle_only_diagnostic()` -- kept per-dataset so a
second entry can't silently reproduce session 1's u0-sign bug).

`plot_lightcurve.py`, `fit_lightcurve.py`, `plot_lightcurve_binary.py`, and
`joint_fit.py` no longer exist as separate files -- their logic was folded
into `mcmc_fit.py`/`mcmc_fit_binary.py` (gated by `--stage`) and
`preprocess_binary_data.py` respectively, to kill single-caller files and
give a fast/no-MCMC path without a second script to keep in sync.

`scratch/` holds every one-time/dev script, not just their output:
`fit_2l1s.py`, `mcmc_fit_2l1s.py`, `compare_2l1s_fits.py`,
`cross_check_mulensmodel.py`, `cross_check_mulensmodel_fit.py`,
`cross_check_mulensmodel_parallax.py`, `derive_binary_quintic.py`,
`fit_2l1s_moa19008.py`. None are covered by the pip line above --
`derive_binary_quintic.py` needs `sympy`, the three `cross_check_mulensmodel*.py`
need `MulensModel`. All are self-documented in their own docstrings and run
manually as `python3 scratch/<name>.py` from the project root (each has its
own `sys.path` line so that works despite living one directory down). If one
of these ever graduates into the pipeline, move both its code and its output
path out of `scratch/` at the same time.

`fit_2l1s_moa19008.py` fits a new event, MOA-2019-BLG-008 (see
dataset_names.txt's `M-19-BLG008`), using only its KMT I-band subset of
`data/MOA-2019-BLG-008L.dat` for now -- **currently produces garbage fits
(chi2/dof ~100-350) because that raw data's own `mag_err` column massively
underestimates the true photometric scatter (confirmed model-independently,
~16-49x); needs an error-bar rescaling step before its output means
anything.** See CHANGELOG.md.


There is no test suite or lint step; each script is run directly and its
correctness is judged by inspecting the plot/printed fit it produces.


## Two datasets, two different pipelines

This repo fits a Paczynski point-source point-lens (PSPL) model
(`lc_models.py`: `trajectory`, `magnification`, `flux`, `magnitude`, plus a
point-source binary-lens (2L1S) extension -- see below) to two
real, published microlensing events, at different levels of complexity:

- **O-05-BLG086** (OGLE-2005-BLG-086): single instrument, single band.
  `mcmc_fit.py` is the whole pipeline, staged via `--stage=raw|quicklook|
  mcmc`; it's also imported as a module (see below) for `fit_pspl_mcmc()`.

- **O-03-BLG235** (OGLE-2003-BLG-235 / MOA-2003-BLG-53): the first
  microlensing planet detection, a real binary-lens/caustic-crossing event,
  observed by two instruments in incompatible units (OGLE: calibrated I
  magnitude; MOA: DIA-style flux *relative to a template*, not a total
  flux -- it's what makes plain "convert flux to magnitude" invalid here,
  since it can be negative and has no shared zero point with OGLE).

Only PSPL is implemented, so fits to O-03-BLG235 are a known-incomplete
sanity check, not a claim that the model is correct -- the caustic
anomaly is visible in the auto-zoomed panel as points sitting above the
smooth fitted curve (most obviously in the MOA data, which has the cadence
to resolve it; OGLE's sparser sampling mostly misses it).

## Annual parallax + robust (Student-t) likelihood for PSPL fits

The plain PSPL model doesn't fit real data well (checked on O-05-BLG086:
chi2/dof=2.14, standardized-residual std=1.46 against the pre-parallax
model -- errors running too tight, moderate excess kurtosis, and the
largest residuals cluster in specific time windows rather than scattering
randomly, consistent with unmodeled physics rather than bad photometry).
Two changes, planned across *every* PSPL-based fit (`mcmc_fit.py`, both
functions in `mcmc_fit_binary.py`, and `preprocess_binary_data.py`'s
`fit_joint_pspl()` calibration fit -- not the `scratch/` 2L1S code, a
separate track) -- **done for O-05-BLG086 (`mcmc_fit.py`), not yet
propagated to O-03-BLG235 (`mcmc_fit_binary.py`/`preprocess_binary_data.py`
-- deliberately deferred)**:

1. **Annual parallax** (Gould 2004 geocentric formalism): two new free
   parameters `piE_N`/`piE_E`, perturbing the trajectory via the target's
   sky position and Earth's orbital motion. `astropy` is a pipeline
   dependency for this (`get_body_barycentric_posvel`, analytic
   position+velocity, no finite-differencing). `t0_par` (the reference
   epoch the geocentric frame is anchored to) is fixed at each dataset's
   preliminary non-parallax best-fit `t0` (posterior median), per
   convention -- not a fitted parameter itself.
2. **Robust likelihood**: the Gaussian log-likelihood is replaced with a
   Student-t likelihood, with both `scale` and `dof` fit as free MCMC
   parameters (motivated by the residual diagnostic above -- a globally
   heavier-than-Gaussian error distribution, not a small population of
   discrete outliers, so a smooth down-weight beats Huber/sigma-clipping).

Both changes **replace the existing model outright** -- no flag/toggle;
the old Gaussian-PSPL behavior stays recoverable via git history only.

**Status (O-05-BLG086)**: both done and validated end to end.
`lc_models.trajectory()`/`flux()`/`magnitude()` take `piE_N`, `piE_E`,
`delta_sN`, `delta_sE` as plain arguments (the latter two precomputed once
by `sun_earth_projection(t, coords, t0_par)` -- never recomputed inside a
per-MCMC-step call, since it's a real astropy ephemeris query).
`plain_flux()`/`plain_magnitude()` (piE_N=piE_E=0) exist for the
`t0_par` bootstrap fit only. Sign/unit convention cross-checked against
MulensModel's own parallax model (`scratch/cross_check_mulensmodel_parallax.py`,
mirroring `cross_check_mulensmodel.py`'s pattern) to 6.8e-7 relative
agreement -- needed two sign fixes to get there (`sun_earth_projection()`'s
overall sign, and `trajectory()`'s `delta_beta` combination sign), found by
testing all 8 combinations of (heliocentric vs. bare-barycentric Earth
position) x (delta_s's own sign) x (delta_beta's combination sign) rather
than guessing. `mcmc_fit.fit_parallax_pspl_mcmc()` runs the real 9-param
(7 light-curve + `scale` + `dof`) MCMC; `mcmc_fit.get_t0_par()` bootstraps
`t0_par` (see "Shared plotting/MCMC helpers" below). Real O-05-BLG086
result: chi2/dof 2.14 -> 1.499, `piE_N`=0.27(+0.04/-0.04), `piE_E`=0.10(1),
`scale`=1.08(5), `dof`=9(+5/-3) -- a real, well-constrained parallax
detection, and `scale`/`dof` closer to Gaussian than the pre-parallax
diagnostic's own MLE (dof~6), consistent with parallax explaining away
real signal rather than the two changes competing for the same residuals.

Target coordinates gathered so far: O-05-BLG086 RA 18h04m45.70s / Dec
-26d59m15.5s (OGLE-III EWS alert page, field BLG234.6 -- a Wikipedia-
sourced pair used briefly during development was off by ~3 degrees in
both RA and Dec and was caught before it reached any fit); O-03-BLG235 RA
18h05m16.35s / Dec -28d53m42.0s (Bond et al. 2004 / NASA Exoplanet
Archive).

## O-03-BLG235's calibrate-once-then-fit pipeline

Because OGLE and MOA can't be compared directly, there's a dedicated
reduction step:

1. **`preprocess_binary_data.py`'s `fit_joint_pspl()`** fits a single shared
   trajectory (t0, u0, tE) to both instruments at once, each with its own
   flux model: `OGLE flux = fs_ogle * A(t) + fb_ogle`,
   `MOA flux = fs_moa * (A(t) - 1)` (no `fb_moa`: MOA's differencing
   already cancels out any constant blend). This is the only place those
   per-instrument flux parameters exist.
2. The rest of `preprocess_binary_data.py` runs that fit once and uses the
   resulting (fs_ogle, fb_ogle, fs_moa) to invert both instruments onto a
   common, physically meaningful, dimensionless scale: magnification A(t)
   (`A = (F_ogle - fb_ogle) / fs_ogle` and `A = 1 + F_moa / fs_moa`). It
   writes `data/processed/O-03-BLG235_{OGLE,MOA}_magnification.dat`.
3. Everything downstream (`mcmc_fit_binary.py --stage=raw` and its
   `run_joint_fit()`) reads only those processed files and fits the plain
   3-parameter PSPL model against the combined dataset -- no
   per-instrument bookkeeping.


**This freezes the flux calibration at step 1's point estimate.** The
(t0, u0, tE) posterior from `run_joint_fit()` is therefore tighter than a
fully joint fit would give, because it doesn't propagate calibration
uncertainty. Fine for exploration; re-run `preprocess_binary_data.py`
whenever the raw data changes, and if a rigorous error budget is ever
needed, that requires one true joint MCMC over all 6 parameters instead.

`mcmc_fit_binary.py`'s other function, `run_ogle_only_diagnostic()`, is a
deliberately naive diagnostic: it fits PSPL to OGLE's data *alone* (no MOA,
no calibration step) to show that a single sparse instrument can look
deceptively well-fit even on a real binary-lens event. Its outputs are
suffixed `_ogle_only` specifically so they don't collide with
`run_joint_fit()`'s canonical outputs. `python3 mcmc_fit_binary.py` runs
both functions in sequence.

`run_ogle_only_diagnostic()` gets its actual 5-parameter PSPL+MCMC fit from
`mcmc_fit.fit_pspl_mcmc()` (shared with the O-05-BLG086 pipeline) rather
than duplicating it. **`u0_guess`/`tE_guess` are required, dataset-specific
arguments to that function, not incidental defaults** -- `curve_fit` can
converge to the mirror-image (`u0 -> -u0`) solution from a bad starting
point, since the model only depends on `u0` squared. Silently reusing one
dataset's guess for another is a real, previously-hit bug (wrong sign,
artificially tight MCMC posterior), not just a style nit.

## O-03-BLG235's PSPL-vs-2L1S model comparison (in progress)

`lc_models.py` also implements a point-source binary-lens (2L1S) model:
`binary_trajectory`, `lens_position`, `binary_images`,
`binary_magnification`. Image positions come from a degree-5 polynomial
whose coefficients were derived symbolically with `sympy`
(`scratch/derive_binary_quintic.py`, a one-time dev script -- not a project
dependency) rather than hand-transcribed from a paper, then hardcoded as
plain numpy. Validated three ways internally (`q->0` recovers PSPL, image
count is always 3 or 5, magnification stays finite/positive off-caustic)
and cross-checked against MulensModel as an independent oracle
(`scratch/cross_check_mulensmodel.py`, also one-time -- needs a separate
`pip install MulensModel`) to ~1e-9 agreement.

`scratch/fit_2l1s.py` fits `(t0, u0, tE, alpha, s, q)` directly against raw
OGLE+MOA flux, not the `data/processed/O-03-BLG235_*_magnification.dat`
files `run_joint_fit()` uses -- `chi2()`/`profile_flux()`/`_profile_fs_moa()`
re-solve each instrument's flux calibration (`fs_ogle`/`fb_ogle`/`fs_moa`)
analytically (closed-form weighted least squares) for every trial
trajectory, since that frozen one-time PSPL-based calibration was shown to
be measurably biased for a real caustic-crossing trajectory (chi2 31,982
vs. 1,953 at Bond et al.'s own published parameters -- see CHANGELOG
session 5/6). This supersedes session 2's planned "confirmatory flux-space
re-fit" step; it's now baked into every trial instead of a one-time check.
`chi2(theta, use_ogle=False)` / `run_fit_moa_only()` /
`mcmc_fit_2l1s.run_mcmc_moa_only()` fit MOA alone (only `fs_moa` profiled,
no cross-instrument calibration at all) as an independent,
calibration-ambiguity-free control on the joint fit. Multi-starts over
`(s, q, alpha)` (seeded from the PSPL joint fit's `t0`/`u0`/`tE`) across
close/resonant/wide topologies, parallelized across seeds with
`ProcessPoolExecutor` using the spawn context -- forking after the parent
process has touched `torch` (which spins up its own thread pool) can
deadlock a second pool; hit for real running `run_fit()` then
`run_fit_moa_only()` back to back. Runnable standalone via
`python3 scratch/fit_2l1s.py` (runs joint then MOA-only in sequence).
`scratch/mcmc_fit_2l1s.py`'s `run_mcmc()` persists its chain
(`samples`/`log_probs`) to `scratch/{name}_2l1s_mcmc_chain[_moa_only].npz`.

**Status**: the calibration bug is fixed and independently validated, but
the search itself is now the open problem. Under the fixed calibration,
multi-start and most MCMC runs converge to non-physical "near-miss cusp"
solutions (no 3->5 image-count transition along the trajectory) that score
*numerically better* than Bond et al.'s own published chi2 -- only
catchable via the existing image-count check, not from chi2 alone. Exactly
one result so far is a confirmed genuine caustic crossing: a joint MCMC run
refined with a tight Nelder-Mead fit of its best sample, chi2=1825.35,
`alpha`~209 deg (vs. Bond's 223.8 deg) -- the current best validated 2L1S
solution, but session 2's BIC-vs-PSPL verdict is still blocked until this
search problem is resolved. `scratch/compare_2l1s_fits.py` overlays every
candidate found so far (both search modes x both instrument scopes, plus
Bond's reference) on one figure for comparison.


## Output layout

- `raw_lc/`, `fit_lc/`, `hist_plots/`, `corner_plots/` hold ONLY the regular
  pipeline's canonical outputs (the ones documented in "Setup and commands").
- `scratch/` holds both the code and the output of every one-time/
  diagnostic/not-yet-graduated script (currently: fit_2l1s.py,
  mcmc_fit_2l1s.py, compare_2l1s_fits.py, cross_check_mulensmodel.py,
  cross_check_mulensmodel_fit.py, cross_check_mulensmodel_parallax.py,
  derive_binary_quintic.py, fit_2l1s_moa19008.py). If a script
  isn't part of the documented pipeline, both it and its plots go in
  scratch/, never in the four dirs above -- when a script graduates into
  the pipeline, move both its code and its output path out at the same
  time.

## Shared plotting/MCMC helpers -- do not re-copy these

- `zoom_utils.plot_fit_panels(ogle, moa, model_fn, fit_label, out_path)`:
  the two-panel (full baseline + auto-zoomed peak) OGLE+MOA data+fit overlay.
  Every magnification-space fit script should call this rather than
  hand-rolling the panel loop again -- it's already been duplicated and
  de-duplicated once.
- `mcmc_fit.save_corner(samples, labels, truths, out_path)`: the
  corner-plot save. Same rule -- call it, don't copy its 4 lines.
- `mcmc_fit.plot_fit_lc(...)`: O-05-BLG086's fit-overlay plot -- full
  baseline + auto-zoomed peak, each paired with a standardized-residual
  scatter panel (+-1 sigma band, symmetric+inverted y-axis matching the
  magnitude panel above it) and a rotated residual histogram sharing that
  panel's y-axis. `_plot_residual_panel()`/`_plot_residual_hist()` are its
  internal per-panel helpers -- call `plot_fit_lc()` itself, don't copy them.
- `mcmc_fit.FitResult`: a `NamedTuple` (`best_fit`, `samples`, `labels`)
  with a `.column(name)` method, returned by `fit_pspl_mcmc()` and
  `fit_parallax_pspl_mcmc()` instead of a bare `(best_fit, samples)` tuple.
  Bundles samples with their own column labels so a caller can't pair a
  samples array with the wrong labels constant (a real bug this session --
  see CHANGELOG). Always pull a column via `.column(name)`, never
  `samples[:, i]` by hand.
- `mcmc_fit.get_t0_par(time, mag, mag_err, cache_path, u0_guess, tE_guess)`:
  the parallax fits' `t0_par` bootstrap (see above) -- reads `cache_path`
  if present, else runs the plain PSPL fit and writes it, always returning
  the posterior median. Self-contained; a caller never touches the cache
  file format. `u0_guess`/`tE_guess` have no default (see `u0_guess`/
  `tE_guess` note below) -- always pass the dataset's own values.
- Before adding a new fit/plot script, check whether its plot shape already
  matches one of these helpers. If it's a genuine structural difference
  (like scratch/fit_2l1s.py's extra caustic-geometry panel), it's fine to stay
  custom -- don't force an abstraction over a real difference.


## Conventions shared across both datasets

- **Time** is `HJD - 2450000` everywhere. Raw files that already use it
  (O-05-BLG086) are loaded as-is; files in full JD/HJD (O-03-BLG235's OGLE
  and MOA tables) get `- 2450000.0` inlined right after loading.
- **NASA Exoplanet Archive table format** (O-03-BLG235's raw files): metadata
  lines start with `\` or `|` before the data rows. Load with
  `np.loadtxt(path, comments=("\\", "|"), unpack=True)` -- no custom parser
  needed.
- **Dataset short names** (used in output filenames) follow
  `dataset_names.txt` -- always add new mappings there rather than
  inventing filename abbreviations inline.
- **`zoom_utils.find_zoom_window`** auto-detects the time window around a
  light curve's peak (by thresholding excess above the median baseline) so
  plots can zoom in without a hardcoded per-dataset window. Fits always run
  on the full dataset; only the plotted view is zoomed.
- Raw data under `data/` (and `data/processed/`) is never edited by hand --
  everything derived is regenerated by re-running the relevant script.
