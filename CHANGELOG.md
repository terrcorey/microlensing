# Changelog

Dated session entries, not semantic versions -- this repo isn't an
application anyone installs. Each entry: what was **Built**, what was
**Learned & open questions** (including negative results -- that's the
actual point of this repo), and **Next session**'s planned goal.

## 2026-09-14 — session 1

### Built
- Two real, published microlensing events set up end to end: **O-05-BLG086**
  (OGLE-2005-BLG-086, single-lens) and **O-03-BLG235**
  (OGLE-2003-BLG-235/MOA-2003-BLG-53, the first microlensing planet
  detection -- binary-lens, two instruments).
- `lc_models.py`: PSPL functional forms (`trajectory`, `magnification`,
  `flux`, `magnitude`).
- Per-dataset pipeline: raw light curve plot -> fit overlay (auto-zoomed on
  the peak via `zoom_utils.find_zoom_window`) -> MCMC + corner plot.
  `mcmc_fit.fit_pspl_mcmc()` is shared between O-05-BLG086's pipeline and
  O-03-BLG235's OGLE-only diagnostic fit.
- For O-03-BLG235: `joint_fit.py` fits a single shared trajectory across
  both instruments' native units, `preprocess_binary_data.py` uses that
  calibration to convert both onto a common physical scale (magnification
  A(t)), written once to `data/processed/`; everything downstream fits the
  plain 3-parameter PSPL model against the combined, calibrated data.
- Ran a ponytail-audit cleanup: deleted `time_utils.py` (inlined a 1-line
  subtraction), replaced a hand-rolled NASA-Exoplanet-Archive table parser
  with `np.loadtxt(..., comments=("\\","|"))`, folded single-caller helper
  functions into closures, de-duplicated the two MCMC-fitting scripts.

### Learned & open questions
- O-05-BLG086's independently-fit PSPL parameters reproduce the published
  (Wyrzykowski et al. 2015) solution to ~2% -- a real, useful validation
  that the pipeline's fitting is correct.
- Even a clean single-lens fit shows real parameter degeneracies in its
  MCMC posterior: u0-tE anti-correlation, f_source-f_blend
  anti-correlation. The best-fit point alone hides this.
- Fitting single-lens PSPL to O-03-BLG235 (a *known* binary-lens event)
  using only OGLE's sparse data gives a deceptively plausible fit and a
  clean-looking unimodal corner plot (chi2/dof=2.0) -- neither the fit
  quality number nor the posterior shape reveals the model is wrong. The
  caustic anomaly only becomes visible by zooming into the peak and
  comparing against MOA's higher-cadence coverage of the same event.
- OGLE and MOA report incompatible units (calibrated magnitude vs.
  template-relative differential flux) with no valid direct conversion
  between them (MOA's column can be negative, and lacks a shared zero
  point). Resolved by fitting a joint calibration once, then converting
  both instruments to a shared *physical* quantity (magnification A(t))
  rather than either instrument's native units.
- Real bug hit while consolidating two MCMC scripts into one shared
  function: silently reused one dataset's hardcoded initial guess for the
  other. Since PSPL's u0 only enters the model squared, `curve_fit`
  converged to the mirror-image (wrong-sign) solution with an artificially
  *tight* posterior -- a wrong answer that still looked converged and
  plausible; chi2 didn't flag it. Only caught by comparing against a
  previous run's known values. Lesson: per-dataset fit configuration
  (initial guesses) is not safe to treat as a shared constant.
- Central open question for the project's actual goal: PSPL's fitted
  parameters alone (t0, u0, tE, fs, fb) cannot give a lens mass or
  distance -- tE conflates the angular Einstein radius and the relative
  proper motion, with no way to separate them from single-band photometry
  alone (needs finite-source effects, parallax, or an independent proper
  motion). This mass-distance-velocity degeneracy is exactly why
  brown-dwarf-mass lenses are hard to identify/characterize from
  microlensing alone.

### Next session
- Propagate the full MCMC posterior (not the best-fit point) through to
  derived physical quantities: peak magnification, blend fraction, source
  baseline magnitude, effective event duration -- plus an explicit,
  prior-based mass/distance estimate (assumed relative proper motion and
  source/lens distance distributions) propagated through the same
  posterior, so the result is a full distribution of plausible lens
  masses rather than one number. The goal is to make the mass-distance
  degeneracy visible directly in the output, not hidden behind a point
  estimate.

## 2026-09-14 — session 2

### Built
- Derived physical quantities from the existing posterior samples, added to
  all three fit entry points (`mcmc_fit.py`, and both
  `run_ogle_only_diagnostic()`/`run_joint_fit()` in `mcmc_fit_binary.py`):
  `A_max` (peak magnification), `t_eff` (effective timescale, `u0*tE`),
  plus `blend_fraction`/`m_source` wherever `f_source`/`f_blend` exist (not
  in the 3-param joint fit, which only has t0/u0/tE).
- `estimate_mass()` (`mcmc_fit.py`): a prior-based lens mass estimate.
  Pairs each posterior `tE` sample with one Monte Carlo draw from assumed
  priors on relative proper motion (log-normal, median 4 mas/yr) and lens
  distance (uniform over 0.1-8 kpc, source fixed at the standard 8 kpc
  bulge distance), via `theta_E = tE*mu_rel`, `M = theta_E^2/(kappa*pi_rel)`.
  Produces a full `M_lens` *distribution* per fit, not a point estimate --
  this directly answers session 1's central open question. Flagged with a
  `ponytail:` comment: the uniform-D_L/fixed-D_S prior is a naive
  line-of-sight prior, not a real Galactic density model; upgrade path is a
  disk/bulge-weighted prior (e.g. Han & Gould 2003) if the mass
  distribution's shape needs to be quantitatively trustworthy rather than
  illustrative.
- Persistent, human-readable output for every fit (previously stdout-only):
  `save_summary()` (`mcmc_fit.py`, shared across all three fit functions)
  writes `data/processed/{name}_fit_summary.dat` -- fixed-width-aligned
  `param p16 p50 p84` rows for every raw + derived quantity, regenerated
  each run, following the existing `data/processed/` convention.
- `plot_histograms()` (`mcmc_fit.py`, also shared): a full histogram grid
  per fit, one panel per raw + derived quantity, saved to new
  `hist_plots/{name}.png`. Auto-switches to log-spaced bins/log x-axis for
  any strictly-positive quantity with >100x dynamic range (needed for
  `M_lens`, whose prior gives it a long tail as drawn `D_L` approaches
  `D_S` and `pi_rel -> 0`) -- a first version with linear bins crushed that
  panel into one unreadable bin near zero.

### Learned & open questions
- Lens mass estimates: O-05-BLG086 ~1.29 (+6.50/-1.08) Msun;
  O-03-BLG235 ogle-only ~0.28 (+1.47/-0.24) Msun;
  O-03-BLG235 joint ~0.25 (+1.26/-0.21) Msun. The order-of-magnitude,
  strongly asymmetric spread is dominated by the assumed mu_rel/D_L
  priors, not by how well tE itself is constrained -- i.e. the
  mass-distance degeneracy is now visible directly in the output, which
  was the actual point of this step.
- Checked every fit's full posterior (raw + derived, via the new histogram
  grids) for hidden multi-modality that a 16/50/84 summary alone could
  mask -- none found in any of the three fits, including
  `run_ogle_only_diagnostic()`'s known-wrong single-lens fit to a real
  binary-lens event. That fit's posterior looks just as clean and
  unimodal as the correct ones; its wrongness is only visible by comparing
  against MOA's independent data, not from its own posterior's shape.
  Reinforces (doesn't just repeat) session 1's finding that neither
  chi2/dof nor posterior shape can flag this kind of model
  misspecification on their own.

### Next session
- Confirmed via `/grill-me`: build a real PSPL-vs-2L1S (binary-lens) model
  comparison for O-03-BLG235, producing a statistically grounded verdict
  rather than a residual eyeball check -- needed as reusable groundwork
  for future brown-dwarf population work, not just this one event.
  - **Model**: point-source 2L1S only for now (no finite-source -- it'll
    underfit the exact caustic peak, accepted as a separate, explicitly
    deferred future step, matching `lc_models.py`'s own docstring).
  - **Implementation**: hand-rolled magnification via lens-equation
    root-finding (quintic polynomial + spurious-image filtering),
    alongside `trajectory`/`magnification`/`flux` in `lc_models.py` --
    chosen deliberately over reusing an existing library (e.g.
    MulensModel, which this project's raw data already comes from) for
    the learning value, despite the real risk that a subtle root-finding
    bug would silently corrupt the comparison. Validated two ways before
    it's trusted: (1) internal limiting-case checks -- `q->0` recovers
    the existing PSPL magnification, image count matches the known
    3-or-5 invariant, magnification stays finite/positive off-caustic;
    (2) a one-time offline cross-check against MulensModel at a handful
    of parameter combinations (including near-caustic) as a trusted
    oracle only, not a pipeline dependency.
  - **Fitting strategy**: small multi-start over qualitatively distinct
    topologies (close/resonant/wide binary regimes) for `(s, q, alpha)`,
    Nelder-Mead refinement at each, best chi2 wins. A full systematic
    grid search is deferred future work. The published Bond et al. 2004
    solution is used only as a sanity check the converged fit should
    land near -- never assumed correct or hardcoded in.
  - **Data**: the 2L1S search runs against the existing PSPL-calibrated
    `data/processed/O-03-BLG235_*_magnification.dat` (6-parameter fit,
    same architecture as today's `run_joint_fit()`). One confirmatory
    flux-space re-fit at the converged solution (no multi-start needed
    there) empirically measures whether `fs`/`fb` actually shift from
    the PSPL calibration -- quantifying, not just flagging, the bias
    from calibrating with the wrong (PSPL) trajectory shape.
  - **Verdict statistic**: BIC (`chi2 + k*ln(N)`, N~1535 combined OGLE+
    MOA points), `deltaBIC > 10` read as decisive (Kass & Raftery). If
    `deltaBIC` comes out only marginally under that threshold, that's
    the documented trigger to escalate to a real Bayes factor via nested
    sampling (e.g. `dynesty`) rather than accept an ambiguous call --
    logged as future work, not built now.
  - **Pipeline validation**: also run the same comparison against
    O-05-BLG086 (known single-lens) as a negative control -- not
    trustworthy as a classifier until confirmed not to falsely flag a
    genuine single-lens event, not just confirmed to catch the one
    binary-lens case already in hand.
  - **Deferred future work, explicitly logged rather than dropped**: a
    full systematic (s,q) grid search; a from-scratch flux-space 2L1S
    calibration (instead of reusing the PSPL-calibrated data); a real
    Bayes-factor comparison via nested sampling; finite-source effects.
  - **Collaboration mode for this build**: user writes the implementation
    directly; Claude's role is to guide/sequence/review (ponytail active
    throughout) rather than write the code, per explicit request.

## 2026-09-14 — session 3

### Built
- 2L1S (point-source binary-lens) magnification implemented in
  `lc_models.py`: `binary_trajectory`, `lens_position`,
  `_quintic_coefficients`, `binary_images`, `binary_magnification`.
  Built across a mid-session collaboration-mode change: user wrote
  `binary_trajectory`/`lens_position`; the quintic/root-finding/
  magnification core was implemented by Claude directly (an explicit,
  scoped override of session 2's "user writes, Claude guides" mode)
  after briefly considering and then reversing a switch to an existing
  package (e.g. MulensModel) via `/grill-me` -- net decision: keep
  hand-rolling the model, Claude takes over this one piece.
- Quintic polynomial coefficients derived symbolically with `sympy`
  rather than hand-transcribed from Witt & Mao (1995), to remove
  transcription-error risk; the derivation (`derive_binary_quintic.py`)
  ships alongside the generated coefficients for auditability. `sympy`
  is a one-time dev tool, not a project dependency.
- Validated `binary_images`/`binary_magnification` three ways before
  trusting them: `q->0` recovers the existing PSPL `magnification()` to
  ~2e-8 relative error; image count matches the known 3-or-5 invariant
  across 300 random `(s,q,zeta)` combinations with zero violations;
  magnification stays finite and >=1 off-caustic across three
  topologies.
- Cross-checked against MulensModel (`cross_check_mulensmodel.py`, a
  one-time oracle check, not a project dependency) across 6 `(s,q)`
  topologies (close/resonant/wide x two mass ratios) x far-field/
  near-caustic points: worst-case relative difference 1.4e-9 after
  fixing a coordinate-convention bug in the *cross-check script itself*
  (see Learned).
- Started `fit_2l1s.py`: a 6-parameter `(t0,u0,tE,alpha,s,q)` fit
  against the existing PSPL-calibrated combined OGLE+MOA magnification
  data, multi-started only over `(s,q,alpha)` (seeded from the PSPL
  joint fit's posterior median for `t0,u0,tE`) across close/resonant/
  wide topologies and four approach angles, parallelized across seeds
  with `ProcessPoolExecutor`, two-stage tolerance (loose for the 24-seed
  exploratory pass, tight for the final refit of the winner). Built
  with Claude back in the session-2 "user writes, Claude reviews" mode,
  via iterative bug-fixing.

### Learned & open questions
- A cross-check against a second, independent implementation isn't
  automatically trustworthy just because it's "the trusted oracle" --
  it needs its own convention audit. The first MulensModel comparison
  pass showed 1-6% disagreement far from caustics and up to ~250x
  disagreement near them, which looked like a real bug in
  `binary_magnification`; it was actually a coordinate-convention bug
  in the comparison script (assumed `alpha=0` meant no rotation, but
  MulensModel's own docs note its `alpha` is shifted 180 deg from the
  convention it otherwise follows -- confirmed empirically, not just
  trusted from the docstring). Same *category* of failure as session
  1's u0-sign bug: a wrong-but-plausible result that superficially
  looks like the thing you were testing for.
- Real bugs hit and fixed one-by-one while building `fit_2l1s.py`'s
  multi-start plumbing: `np.concatenate` called with two positional
  args instead of a list; `binary_magnification` called missing `s,q`;
  a closure-based `chi2`/`_fit_one_seed` that `ProcessPoolExecutor`
  can't pickle (functions handed to worker processes must be defined at
  module scope, not nested); a seed-shape mismatch (`alpha, s, q =
  seed` against a 6-element seed); `min(results, key=...)` misused to
  unpack a single `OptimizeResult` as if it were a `(chi2, theta)`
  tuple (`OptimizeResult` is dict-like, not positionally indexable).
  None were physics bugs -- `binary_magnification` itself was already
  independently validated before this file existed.
- `binary_magnification` is expensive: it root-finds per data point
  (not vectorizable across time samples), so a single tight-tolerance
  Nelder-Mead run took ~3m13s standalone. The planned 24-seed
  multi-start grid needed both parallelization (`ProcessPoolExecutor`,
  one process per seed) and a two-stage tolerance (loose while
  searching for the right basin, tight only on the eventual winner) to
  stay tractable at all.
- The multi-start run's slowness was diagnosed by attaching `py-spy`
  (read-only, no interruption) to the live process rather than guessing:
  the parallel multi-start phase itself had already finished in a few
  minutes (winning seed: chi2~2139 in 647 iterations), and the process
  was stuck for 35+ minutes entirely inside the *final* refit's
  `xatol=fatol=1e-8` Nelder-Mead call -- 3798 iterations / 28,749
  function evaluations spent shaving chi2 from 2139.2092354 to
  2139.2092229, an utterly negligible improvement. Root cause: 1e-8 is
  far tighter than a 6-parameter, expensive-per-call objective needs.
  Killed the run; `fit_2l1s.py`'s final-refit tolerance loosened to
  1e-6. The near-converged parameters pulled from the live process were
  promising -- `s~0.780, q~0.00246, alpha~1.273` -- a small, planet-like
  mass ratio, unlike the single-start run's `q~0.103`.
- A single, non-multi-started Nelder-Mead run from one `(alpha=0,
  s=1.0, q=0.1)` guess converged to `chi2~18000` (double-digit
  chi2/dof against ~1500 combined points) with `s,q` barely moved from
  the seed -- read as "stuck in that seed's own basin" (the expected
  symptom motivating multi-start), not a sign `binary_magnification`
  itself is wrong, since that was independently confirmed via the
  MulensModel cross-check beforehand.

### Next session
- Confirmed via `/grill-me`: re-run the now-faster multi-start fit to
  get its actual final answer, sanity-check the converged `(s,q,alpha)`
  against Bond et al. 2004's published solution, then continue session
  2's already-agreed plan in order: the confirmatory flux-space re-fit
  (quantifying how much `fs`/`fb` shift once the true 2L1S trajectory
  replaces the PSPL one used to calibrate the data), the BIC verdict
  against the existing PSPL joint fit (`deltaBIC > 10` decisive), and
  the O-05-BLG086 negative control.
- Also add a caustic-crossing close-up plot (zoomed on the entry/exit,
  like Bond et al. 2004's Fig. 2) into the regular pipeline, not just a
  one-off diagnostic.

## 2026-09-15 — session 4

### Built
- Re-ran the multi-start fit with progressively wider `q`/`alpha` grids
  (24 -> 36 -> 72 -> 32 seeds, `alpha` down to 22.5 deg spacing) -- every
  run converged to the same `chi2=2139-2161` basin regardless of grid
  density, never the caustic-crossing solution.
- Resolved the `alpha`-vs-Bond et al.'s `phi` convention empirically
  (no MulensModel-style 180 deg offset): plugging their raw published
  values straight into `binary_trajectory`/`binary_magnification`
  reproduces their reported caustic entry/exit days (2835, 2842) to
  within 0.05-0.2 days, and visually matches the real light curve's
  spike timing -- confirming the true solution really does sit near
  their parameters, and that the multi-start's repeated failure is a
  search problem, not a data, convention, or model problem.
- Cross-checked with an independent full fit (not just point evaluations):
  `cross_check_mulensmodel_fit.py` runs MulensModel's own optimizer over
  the same seed grid -- it converged to the *identical* wrong basin
  (matching chi2 and parameters to 5 decimals, `alpha` differing by
  exactly MulensModel's documented 180 deg offset), ruling out an
  implementation bug definitively.
- `lc_models.py`'s `binary_images`/`binary_magnification` rewritten to
  batch the per-point root-finding into one call via a batched
  companion-matrix eigensolve (`torch.linalg.eigvals`, CUDA if
  available else CPU) instead of looping `np.roots` per point --
  requested for faster diagnostic iteration. Validated against the
  pre-rewrite implementation (agrees to ~1e-9 on a realistic batch,
  identical q->0 error, and actually 4/300 fewer image-count-invariant
  violations on the same random draws) and against the independent
  MulensModel oracle (1.73e-9 worst-case, matching the original's
  1.4e-9). ~7x faster per chi2 call on CPU alone (69ms -> 9.7ms for
  1535 points); real GPU benefit unverified -- this sandbox has no
  CUDA device. New `torch` dependency; `requirements.txt` added via
  `pip freeze` (hand-filtered to drop a sandbox-only stray package and
  re-pin `torch` off its CPU-only build string).
- `zoom_utils.find_zoom_window` reworked from a raw-amplitude threshold
  (fraction of the single largest excess point) to a significance cut
  (`excess >= sigma_threshold * error`, default 5sigma) -- the old
  version had no stable threshold and could balloon from a 30-day
  window to the entire multi-year dataset on a small parameter change,
  because MOA's per-point errors vary too much for a fixed amplitude
  fraction to separate real signal from baseline noise. Fixed at all
  four call sites.
- Added `lc_models.caustic_curve(s, q)`: traces the source-plane caustic
  by solving the critical-curve quartic (`np.convolve`-built
  coefficients, no conjugate to clear so no symbolic derivation needed)
  and mapping through the lens equation. `fit_2l1s.py`'s `plot_fit` now
  has a third panel plotting it against the source trajectory, making
  "does this fit actually cross the caustic" a visual check instead of
  an ad hoc image-count-transition query.
- Built `mcmc_fit_2l1s.py`: emcee exploration of the full 2L1S parameter
  space (48 walkers, 1500 steps) as an alternative to grid multi-start,
  walkers initialized from broad priors (log-uniform `s`/`q` across ~4
  decades, full-circle `alpha`) rather than a handful of seed points.

### Learned & open questions
- The caustic-panel check caught a fit that looked deceptively good: a
  multi-start result with `chi2=2046` (a new best at the time) whose
  zoomed light curve tracked the data's general shape and had a sharp
  spike right at the single most extreme MOA outlier -- but the image
  count along its trajectory never left 3 (no 3->5 transition), meaning
  it was a near-miss cusp approach, not a real caustic crossing. Same
  category of trap as session 1/2's chi2-can't-detect-wrongness finding,
  now caught by a geometric check instead of an external reference.
- The wrong basin (`alpha~73 deg`) that every multi-start run keeps
  finding is not just easier to seed into than the correct one
  (`alpha~224 deg`) -- it may have a genuinely larger basin of
  attraction. `mcmc_fit_2l1s.py`'s broad-prior run landed there too as
  its single best sample, plus found a second, *more densely populated*
  mode near `alpha~172 deg, s~1.55` that wasn't even on the multi-start
  grid's radar. The corner plot shows scattered, disconnected point
  clusters rather than one well-mixed posterior -- 1500 steps of plain
  `emcee` likely isn't enough to mix across a landscape this rugged, and
  a small bump near `alpha~229 deg` (close to the confirmed-correct 224)
  is visible in the histogram but unquantified -- the raw chain wasn't
  saved, only the corner plot and single best sample.
- The 2L1S PSPL-vs-model-comparison plan from session 2 (flux-space
  re-fit, BIC verdict, O-05-BLG086 negative control) is still blocked on
  actually finding the true caustic-crossing solution -- not yet reached.
- A tight Nelder-Mead refit seeded exactly at Bond et al.'s converted
  parameters (`alpha=223.8 deg`, the confirmed-correct region) barely
  moved in `alpha` (converged to 225.89 deg) but `q` collapsed to
  `0.00001` and `chi2=2305.82` -- worse than the familiar wrong basin.
  Refining *against our own PSPL-calibrated data* pulls even a
  correctly-seeded caustic solution toward degenerating back to q->0
  (single-lens) rather than sharpening into a real binary fit. Concrete
  evidence for the calibration-bias question the session 2 flux-space
  re-fit step already exists to quantify -- our data's own calibration
  may not actually preserve enough of the caustic amplitude to reward
  keeping q away from zero.

### Next session
- User is restructuring the codebase and will refresh CLAUDE.md
  themselves once that's done -- no goals confirmed via `/grill-me` this
  entry; pick up from here once that's finished. Open threads to return
  to: save the MCMC chain and quantify the `alpha~229 deg` bump, or
  switch to a sampler built for multimodal posteriors (e.g. parallel
  tempering) instead of plain `emcee`; then resume session 2's blocked
  comparison plan once a genuine caustic-crossing fit is in hand.

## 2026-09-15 — session 5

### Built
- Confirmed, with hard numbers rather than the standing suspicion from
  session 4, that O-03-BLG235's frozen PSPL-based flux calibration
  (`preprocess_binary_data.py`) is the actual reason every 2L1S search
  (multi-start, MulensModel's own optimizer, `emcee`) has kept landing in
  the wrong basin. Pulled Bond et al. 2004's real Table 1 parameters
  directly from the paper (arXiv:astro-ph/0404309: `t0=2848.06, u0=0.133,
  tE=61.5, alpha=223.8deg, s=1.120, q=0.0039`) and evaluated chi2 against
  them: 31,982 using the frozen calibration vs. 1,953 once `fs_ogle`/
  `fb_ogle`/`fs_moa` are re-solved (closed-form weighted least squares)
  against Bond's real trajectory instead -- right in the ballpark of
  their own reported 1390.49. The frozen values were off by ~40%+ on
  both flux scales, and `fb_ogle` even had the wrong sign. A fix (profile
  `fs`/`fb` analytically per trial instead of freezing them from a
  one-time PSPL point estimate) was designed and the user selected it
  when asked, but then deferred actually implementing it this session --
  not yet applied.
- Added `--dataset` to `mcmc_fit_binary.py` (default `O-03-BLG235`,
  unchanged from before), backed by a small `DATASETS` dict keyed by
  short name (raw OGLE filename + dataset-specific `u0`/`tE` guesses for
  `run_ogle_only_diagnostic()`, kept per-dataset so a second entry can't
  silently reproduce session 1's u0-sign bug).
- Looked up and onboarded a new event: MOA-2019-BLG-008 (Bachelet et al.
  2022, arXiv:2205.07522) -- a real ~30 Jupiter-mass object at the
  planet/brown-dwarf boundary, directly on this project's central theme.
  Its raw multi-survey photometry (`data/MOA-2019-BLG-008L.dat`) has 12
  source/site codes across 6 bands (MOA, OGLE, 4 KMTNet site+chip
  combos, 5 LCO codes) -- far beyond the existing 2-instrument
  architecture, so scoped down to just its KMT I-band subset for a first
  pass (single instrument, so no cross-survey calibration is needed at
  all). New `scratch/fit_2l1s_moa19008.py`: multi-start 2L1S fit against
  that subset, reusing the same analytically-profiled-fs/fb approach
  designed for the O-03-BLG235 fix above. Added `M-19-BLG008` to
  `dataset_names.txt`.

### Learned & open questions
- The session's central result is that session 4's "it's a search
  problem, not a data/convention/model problem" conclusion was itself
  wrong -- it's a calibration bug, just one that happened to also
  produce a real-looking wrong basin every search kept re-finding. The
  tell that finally forced checking it with real numbers (rather than
  leaving it as a documented hunch): seeding a refit exactly at the
  *known-correct* answer scored worse, not better.
- `fit_2l1s_moa19008.py`'s first run: chi2/dof ~100-350 across every
  multi-start seed, and identically bad even for a trivial 3-parameter
  PSPL pre-fit regardless of `alpha`/`s`/`q` -- suspicious enough
  (uniformly bad regardless of the one thing varying between seeds) to
  investigate rather than accept as "just a hard event." Root cause,
  confirmed model-independently: a flat-line fit (no lensing model at
  all) to the baseline region alone gives chi2/dof=2441, because the raw
  KMT `mag_err` column (median 0.006 mag) understates the real
  point-to-point scatter (~0.1-0.2 mag, robust-vs-plain estimate) by
  roughly 16-49x. This is a known property of raw KMTNet/pySIS pipeline
  photometry -- every published microlensing analysis rescales error
  bars (typically to baseline chi2/dof~1) before fitting. A smaller,
  secondary effect was also found and ruled out as the dominant cause:
  the 6 different KMT site+chip codes have ~0.08 mag baseline offsets
  from each other, dwarfed by the error-underestimation issue.
- Both findings this session share a pattern worth remembering before
  the next debugging session: a catastrophically bad chi2 can come
  entirely from something upstream of the physical model (flux
  calibration; error-bar normalization) rather than indict the model or
  the search algorithm.
- Killed the MOA-2019-BLG-008 multi-start run partway through (stuck on
  one slow seed) once the error-bar finding made its result moot --
  no point letting a fit run against data whose error bars are known
  wrong.

### Next session
- Confirmed via `/grill-me`: chase session 4's still-open
  `alpha~229 deg` MCMC bump next, ahead of either calibration fix above.
  Scope is deliberately capped at two steps -- add persistence
  (`np.savez` of `samples`/`log_probs`) to `mcmc_fit_2l1s.py`'s
  `run_mcmc()`, re-run the same 48-walker/1500-step/seed=42 job, then
  inspect the `alpha` marginal (mode population, competitive log-prob)
  -- and stop there; escalating to a multimodal-aware sampler (parallel
  tempering, nested sampling) is explicitly deferred until that
  inspection actually shows a real, separable, poorly-mixed mode
  justifying it. Collaboration mode for this specific task: Claude
  implements it directly, not the user. The O-03-BLG235 calibration fix
  and MOA-2019-BLG-008's error-bar rescaling remain open but
  deprioritized behind this.

## 2026-09-18 — session 6

### Built
- Applied the O-03-BLG235 flux-calibration fix designed in session 5:
  `fit_2l1s.py`'s `chi2()` now profiles `fs_ogle`/`fb_ogle`/`fs_moa`
  analytically per trial via closed-form weighted least squares
  (`profile_flux()`/`_profile_fs_moa()`) against raw flux
  (`lc_models.mag_to_flux()`, `preprocess_binary_data.load_raw()` factored
  out for reuse), instead of the old frozen one-time PSPL calibration.
  Validated: chi2 at Bond et al. 2004's real published parameters now comes
  out to 1952.90, matching session 5's hand-derived estimate (~1,953)
  almost exactly.
- Added the persistence to `mcmc_fit_2l1s.py`'s `run_mcmc()` confirmed at
  the end of session 5 (`np.savez` of `samples`/`log_probs`), then re-ran
  the joint 2L1S MCMC under the fixed calibration and inspected the
  `alpha`/`t0` marginals quantitatively rather than by eye: confirmed the
  chain is severely poorly-mixed -- a 2068-sample-populated `t0` bin is
  still 106 log-units less probable than the true best sample sitting in a
  258-sample bin; a further 494-sample cluster is 163-169 log-units worse.
  Bar height in these histograms tracks walker traffic, not posterior
  probability.
- Via `/grill-me`: scoped and built a MOA-only mode as a
  calibration-ambiguity-free control on the joint fit
  (`chi2(theta, use_ogle=False)`, `run_fit_moa_only()`,
  `run_mcmc_moa_only()`), mirroring `mcmc_fit_binary.py`'s existing
  `run_joint_fit()`/`run_ogle_only_diagnostic()` split rather than
  duplicating scripts. Ran the full controlled comparison -- joint vs.
  MOA-only x multi-start vs. MCMC (4 searches) -- plus a tight Nelder-Mead
  refit of the joint MCMC's best sample, checking every result against the
  existing image-count (3->5 transition) caustic-crossing test rather than
  trusting chi2 alone.
- Fixed a real fork-after-threading deadlock: `fit_2l1s.py`'s multi-start
  `ProcessPoolExecutor` used the default fork context; running `run_fit()`
  then `run_fit_moa_only()` back to back left `torch`'s thread pool holding
  a lock every worker in the second pool then waited on forever. Confirmed
  via `/proc/<pid>/wchan` showing `futex_do_wait` on all workers at 0% CPU
  for over an hour -- genuinely stuck, not just slow -- before fixing by
  switching to the spawn context, matching `mcmc_fit_2l1s.py`'s existing
  precedent (which already carried a comment warning about exactly this
  class of bug).
- Fixed a real plotting bug: the auto-zoomed panel's y-limits were sized
  from the data points alone, not the fitted model curve, so a narrow model
  spike between data points could get visually clipped. Fixed in both
  `zoom_utils.plot_fit_panels` and `fit_2l1s.py`'s own panel loop.
- Reworked `mcmc_fit.save_corner()` (shared by every fit script) to use
  smoothed, `fill_contours`-based gradient shading instead of raw scatter
  points plus a flat "hide the center" block -- applies project-wide.
- Built `scratch/compare_2l1s_fits.py`: a one-off comparison figure
  overlaying every 2L1S candidate found this session (both search modes x
  both instrument scopes, plus Bond et al.'s published solution) --
  light-curve overlay (event-season + auto-zoomed peak, y-axis capped at a
  multiple of the main hump's peak rather than stretching to fit narrow
  caustic spikes) plus a shared-color, zeta=0-centered caustic-geometry
  inset.
- Began, then parked, an HPC path for a wider multi-start search: diagnosed
  an `sshfs` login-username mismatch (not an SSH-key problem as first
  suspected), started scoping a SLURM array-job approach before deferring
  it to run locally instead.

### Learned & open questions
- The central result this session: unfreezing the flux calibration fixed
  the bug session 5 found, but it also made non-physical "near-miss cusp"
  trajectories *more* attractive to chi2 minimization than before -- 3 of 4
  independent searches (joint multi-start, both MOA-only searches)
  converged to solutions that never cross the caustic yet score
  numerically better than Bond et al.'s own published parameters. Only the
  joint MCMC (then refined) found a real crossing.
- MOA-only, despite being calibration-ambiguity-free by construction
  (`fs_moa` is provably identical whether solved jointly or alone), did not
  independently reproduce the joint fit's correct answer in either search
  mode -- inconclusive on whether that's a remaining calibration issue or
  just search inadequacy (weak seed grid, already-known MCMC mixing
  problems); leans toward the latter given `fs_moa`'s proven identity.
- Widening the multi-start seed grid from 16 to 128 seeds (user's own edit,
  `s` in 4 values x `q` in 2 values x 16 `alpha`s) did not change the
  answer at all -- identical chi2=1909.02 near-miss basin found both times,
  reinforcing that this specific wrong basin has an unusually large basin
  of attraction, not just an unlucky seed gap.
- A user-edited linear `s` prior (0.7-1.8, replacing log-uniform) initially
  still used `log`/`exp` in both `log_prior` and `sample_prior` -- caught
  before running (would have silently restricted `s` to (2.01, 6.05),
  excluding every physically plausible value found this session). Once
  fixed, the resulting MCMC run landed on yet another different
  non-crossing near-miss (`alpha`~221deg, `s`=1.624, chi2=2142.58) rather
  than rediscovering the known-good ~209deg solution -- another data point
  for the search-inadequacy read above.
- The current best validated 2L1S solution for O-03-BLG235: chi2=1825.35,
  t0=2847.57, u0=0.094, tE=82.55, alpha~209.2deg, s=1.085, q=0.0062 (joint
  MCMC best sample, tight-refined) -- better than Bond et al.'s own chi2
  (1952.9) and a genuine caustic crossing, but `alpha` is still ~15deg off
  Bond's 223.8deg, an open discrepancy not yet chased down.

### Next session
- No goals confirmed via `/grill-me` this entry (explicitly deferred to
  next time) -- open threads to pick up from:
  - The near-miss-artifact problem: either build the image-count/
    caustic-crossing check directly into the search (reject/penalize
    non-crossing trials) rather than only catching it after the fact,
    and/or escalate to a multimodal-aware sampler (parallel tempering /
    nested sampling) -- both flagged as options, neither chosen yet.
  - MOA-only still hasn't independently confirmed the joint fit's ~209deg
    solution in either search mode -- worth another attempt once the above
    is addressed.
  - The ~15deg `alpha` gap between our best solution and Bond et al.'s
    published value is unexplained.
  - HPC path parked mid-diagnosis: the `sshfs` username mismatch likely
    needs a `User` line in `~/.ssh/config` or an explicit `user@host` in
    the mount command; SLURM script shape (array-per-seed vs. multi-node
    MPI-pool vs. bare scaffold) not yet decided.
  - MOA-2019-BLG-008's KMT error-bar rescaling remains open and
    deprioritized, untouched since session 5.
