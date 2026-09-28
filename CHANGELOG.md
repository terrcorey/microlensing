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

## 2026-09-21 — session 7

### Built
- Via `/grill-me`: scoped a new project direction after the user found the
  existing PSPL fits weren't fitting the data well. Two changes planned
  across *every* PSPL-based fit (`mcmc_fit.py`, both functions in
  `mcmc_fit_binary.py`, `preprocess_binary_data.py`'s `fit_joint_pspl()`
  calibration fit -- not the `scratch/` 2L1S track): add annual parallax
  (Gould 2004 geocentric `piE_N`/`piE_E`, via `astropy`), then replace the
  Gaussian likelihood with a Student-t likelihood (`scale` and `dof` both
  fit freely). Both replace the existing model outright, no flag/toggle.
  O-05-BLG086 is the proving ground; parallax first, then the likelihood
  change. See CLAUDE.md's new "Annual parallax + robust (Student-t)
  likelihood for PSPL fits" section for the full plan.
- A residual diagnostic on the current (pre-parallax) O-05-BLG086 PSPL fit
  motivated the Student-t choice over Huber/sigma-clipping: chi2/dof=2.14,
  standardized-residual std=1.46 (errors running ~46% too tight), only
  moderate excess kurtosis (1.23), Student-t MLE dof~6 -- a globally
  heavier-than-Gaussian error distribution, not a small distinct outlier
  population. Several of the largest residuals cluster in specific time
  windows rather than scattering randomly, consistent with unmodeled
  physics (parallax) rather than bad photometry.
- Gathered target coordinates for the parallax calculation (O-05-BLG086:
  RA 18h04m45.70s / Dec -26d59m15.5s, OGLE-III EWS alert page, field
  BLG234.6; O-03-BLG235: RA 18h05m16.35s / Dec -28d53m42.0s, Bond et al.
  2004 / NASA Exoplanet Archive) and the Gould (2004) `t0_par` convention
  (fixed at each dataset's preliminary non-parallax best-fit `t0`, not
  jointly fit).
- Built Step 1 of the parallax implementation:
  `lc_models.sun_earth_projection(time, ra_str, dec_str, t0_par)` --
  Earth's sky-projected position relative to the Sun (AU), via astropy's
  `get_body_barycentric_posvel` (analytic position *and* velocity, no
  finite-differencing), projected onto the target's North/East
  tangent-plane basis, with the constant-velocity part at `t0_par`
  subtracted out (already degenerate with (t0, u0, tE); only the
  curvature signal is new information). `astropy` added as a new pipeline
  dependency.
- Got `scratch/scratchpad.ipynb` running against the project's own
  `.venv` in VSCode (`ipykernel` installed, kernel selected via
  Cmd+Shift+G to work around the native file-picker hiding dotfiles), for
  interactive development of the parallax code. `%load_ext autoreload` +
  `%autoreload 2` set up so edits to `lc_models.py` take effect without a
  kernel restart.
- Pruned `requirements.txt` back to pipeline-only dependencies -- a
  `pip freeze`-style regeneration during the session had pulled in the
  full notebook/dev toolchain (`ipykernel`, `jupyter_client`, `ipython`,
  `debugpy`, etc.) alongside `astropy`, inconsistent with this project's
  existing convention of keeping dev-only deps out of the main
  requirements file (mirrors how `MulensModel`/`sympy` stay scoped to
  `scratch/`). `astropy`'s own real runtime deps (`astropy-iers-data`,
  `pyerfa`) were kept; `ipykernel` etc. stay installed locally but
  untracked, same treatment as the scratch-only libraries.

### Learned & open questions
- Real, previously-hit-pattern bugs caught and fixed while building
  `sun_earth_projection()`: a `SkyCoord` can't be unpacked as
  `ra, dec = SkyCoord(...)`; `np.dot()` doesn't correctly project an
  astropy `CartesianRepresentation` time series onto a fixed direction
  (needs explicit x/y/z component combination -- pulled into a shared
  `_project()` helper rather than repeating it six times); a `TimeDelta`
  stripped to `.value` before multiplying against a velocity `Quantity`
  silently drops its unit tag, producing a `UnitConversionError` several
  lines later rather than failing at the actual mistake.
- A more serious near-miss: RA/Dec for O-05-BLG086 initially copied from
  Wikipedia (17h54m19.2s / -30d22m38s) differed from the OGLE-III EWS
  alert-page value by ~3 degrees in both RA and Dec -- caught by
  comparing against the coordinates gathered earlier in the session,
  before it propagated into any fit. Another wrong-but-plausible-looking
  input that wouldn't have thrown an error on its own.
- Local dev machine hit an SSL certificate verification failure running
  `download_data.py` under Homebrew's Python (missing local CA bundle) --
  unrelated to the script's own correctness; worked around with `curl`
  for this session, not a repo-level issue.

### Next session
- Confirmed via `/grill-me`: continue the parallax roadmap in order.
  - Step 2: wire `piE_N`/`piE_E` into `lc_models.trajectory()`/`flux()`/
    `magnitude()`.
  - Step 3 (non-optional per this session's plan): cross-check the new
    parallax trajectory against MulensModel's own parallax-enabled PSPL
    model to lock down the sign/unit convention, mirroring the existing
    `scratch/cross_check_mulensmodel.py` pattern used to validate 2L1S.
  - Step 4-5: wire into `mcmc_fit.py`'s MCMC (new params, priors, `p0`)
    and validate against real O-05-BLG086 data (does chi2/dof actually
    improve? is `piE` well-constrained or prior-dominated noise?).
  - Step 6: propagate the same trajectory change into
    `mcmc_fit_binary.py`/`preprocess_binary_data.py` for O-03-BLG235.
  - Robust (Student-t) weighting is planned after all of the above.

## 2026-09-22 — session 8

### Built
- Continued the parallax roadmap in order. Step 2 (wire `piE_N`/`piE_E`
  into `lc_models.py`'s `trajectory`/`flux`/`magnitude`, plus
  `plain_flux`/`plain_magnitude` for the `t0_par` bootstrap) was written by
  the user directly, per this project's "user writes, Claude reviews" mode
  (Claude's own first attempt was reverted at the user's request early in
  the session, once the convention was reasserted). Claude found and fixed
  two incidental bugs while reviewing: `_project()`'s `.radians` -> `.radian`
  typo, and `sun_earth_projection()` returning astropy `Quantity` objects
  instead of plain floats (needed since `piE_N`/`piE_E` are conventionally
  bare AU^-1 numbers, matching MulensModel's own convention).
- Step 3 (non-optional cross-check): built
  `scratch/cross_check_mulensmodel_parallax.py`, testing all 8 combinations
  of (heliocentric vs. bare-barycentric Earth position) x (`delta_s`'s own
  sign) x (`delta_beta`'s combination sign) against MulensModel's own
  parallax model, rather than guessing. Needed two real, simultaneous sign
  fixes (`sun_earth_projection()`'s overall sign; `trajectory()`'s
  `delta_beta` combination sign) to reach 6.8e-7 relative agreement --
  testing only the `delta_beta` sign alone first gave an ambiguous ~3%
  residual that didn't cleanly discriminate, because both bugs were
  present and partially masking each other.
- Steps 4-5: wired parallax into `mcmc_fit.py`'s real MCMC
  (`fit_parallax_pspl_mcmc()`, 9-param: 7 light-curve + `scale` + `dof`).
  A long iterative debugging pass surfaced and fixed a chain of distinct
  real bugs while wiring this up: missing `plain_flux`/`plain_magnitude`
  calls, an `args.state` typo, missing `delta_sN`/`delta_sE` passthrough,
  `curve_fit`'s `p0` including whole precomputed arrays as if they were
  scalar parameters, a hardcoded `ndim=5` left over from copy-paste, a
  `LABELS`-constant length mismatch between the plain and parallax fits, an
  `np.loadtxt` cache read that could never succeed on a mixed string/float
  file, and a `plot_fit_lc()` grid mismatch (`delta_sN`/`delta_sE` computed
  at the data's own time array, then reused against a different, denser
  plotting grid). Verified end to end on both the cache-hit and
  cache-miss code paths.
- Student-t robust likelihood swap done for O-05-BLG086 (session 7's other
  planned change): `fit_parallax_pspl_mcmc()`'s Gaussian log-likelihood
  replaced outright with a Student-t (`scipy.stats.t.logpdf`), 2 new free
  params `scale`/`dof`, no toggle.
- Real O-05-BLG086 result (both changes together): chi2/dof 2.14 -> 1.499;
  `piE_N`=0.27(+0.04/-0.04), `piE_E`=0.10(1) -- a real, well-constrained
  parallax detection; `scale`=1.08(5), `dof`=9(+5/-3) -- closer to Gaussian
  than the pre-parallax diagnostic's own dof~6 MLE, consistent with
  parallax explaining away real signal rather than the two changes fighting
  over the same residuals.
- Ran the `/improve-codebase-architecture` skill on the session's hot spot
  (`lc_models.py`/`mcmc_fit.py`), producing an HTML report of 5 deepening
  candidates. Grilled and implemented two: (1) `FitResult`, a
  `NamedTuple(best_fit, samples, labels)` with a `.column(name)` method,
  replacing the module-level `PLAIN_LABELS`/`PARALLAX_LABELS`/`LABELS`
  constants and every hand-indexed `samples[:, i]` -- the review's own audit
  found this pattern had already caused a live, reproducible `IndexError`
  in `mcmc_fit_binary.py` (`LABELS` grew from 7 to 9 items but that file's
  own samples array is still 5 columns); (2) `get_t0_par()`, replacing the
  inline `try/except` cache block -- fixes a real, previously-silent
  correctness bug where `t0_par` was computed one of two different ways
  (posterior median on a cache hit, `curve_fit` point estimate on a cache
  miss) depending on whether a file happened to exist. `mcmc_fit_binary.py`
  itself deliberately left untouched (already broken for the separate,
  deferred Step 6 reason).
- Repo-wide plotting convention pass (explicit user request, not part of
  the roadmap): every `errorbar()` call (14 across 7 files) now has
  `capsize` plus `markeredgewidth`/`capthick` matching `elinewidth`; every
  plot's `dpi` doubled (150->300, 300->600, 200->400 across 9 `savefig`
  calls).
- `plot_fit_lc()` gained a standardized-residual panel under each of its
  two light-curve panels (+-1 sigma translucent band, symmetric+inverted
  y-axis matching the panel above), plus a rotated marginal histogram
  sharing each residual panel's y-axis -- iterated through layout bugs
  (`tight_layout()` incompatible with the `sharey`-linked grid, causing
  title/label overlap; the light-curve panel spanning the full figure width
  while its residual panel below only spanned the left column, breaking
  date-axis alignment) before landing on a working `GridSpec` layout.
- `/graphify`: built once scoped to the whole Desktop (no path given),
  updated once after the parallax/Student-t work, then rebuilt from scratch
  scoped to just `microlensing/` per explicit request, since the graph had
  been living in `../graphify-out/` -- one level outside the project,
  which CLAUDE.md's own rule on outside-folder access flagged as needing
  the user's sign-off as a standing instruction. `graphify-out/` now lives
  inside `microlensing/`. The old Desktop-level `graphify-out/` was left in
  place, not deleted.

### Learned & open questions
- The 2L1S cross-check pattern (test the oracle comparison exhaustively
  rather than one hypothesis at a time) generalizes cleanly to parallax: a
  single-hypothesis test gave an ambiguous, inconclusive residual because
  two independent sign bugs were both present and partially canceling in
  some regimes. Only testing the full cross-product of candidates isolated
  both.
- Real design-pattern lesson from the architecture review: pairing a
  samples array with a separate "column labels" constant by convention
  (not enforced by any interface) is exactly the shape of bug that doesn't
  show up until a second consumer (`mcmc_fit_binary.py`) needs the same
  building block under slightly different conditions -- confirmed by this
  session's own live `IndexError`. Bundling data with its own labels in one
  object (`FitResult`) removes the failure mode structurally rather than
  adding a check on top of it.
- A second, quieter version of the same lesson: `t0_par` silently had two
  different estimators depending on which branch of a file-existence check
  ran -- no error, no test, just a slightly different downstream answer
  depending on incidental machine state. Silent-but-wrong is a worse
  failure mode than loud-and-wrong, harder to catch by inspection alone.
- `scale`~1.08/`dof`~9 (vs. the pre-parallax diagnostic's own
  scale~1.46-equivalent/dof~6) is a genuine cross-check that parallax and
  the Student-t correction aren't just two knobs fighting over the same
  residual pattern -- if they were, adding parallax wouldn't have changed
  what Student-t needed to compensate for.
- Known small issues, not yet scheduled: the old Desktop-level
  `graphify-out/` is superseded and undeleted; `requirements.txt` lists
  `sympy`/`networkx`, which CLAUDE.md says shouldn't be pipeline
  dependencies (and `networkx` isn't mentioned in any doc at all);
  `README.md`'s setup section claims "no requirements.txt yet" even though
  one exists with 27 pinned packages.

### Next session
- Confirmed via `/grill-me`-style consensus: Step 6 -- propagate both
  changes (annual parallax + Student-t) to O-03-BLG235
  (`mcmc_fit_binary.py`'s two functions and `preprocess_binary_data.py`'s
  `fit_joint_pspl()`), completing the roadmap CLAUDE.md's "Annual parallax
  + robust likelihood" section describes. The 2L1S search problem and
  MOA-2019-BLG-008's `mag_err` rescaling stay backlog, not scheduled.

## 2026-09-22 — session 9

### Built
- Completed session 8's confirmed goal (Step 6): propagated annual
  parallax + Student-t robust likelihood to O-03-BLG235, finishing
  CLAUDE.md's "Annual parallax + robust likelihood" roadmap for both
  datasets. Claude implemented this session (user requested it directly,
  a deliberate one-off exception to this repo's usual "user writes, Claude
  reviews" mode), with ponytail active per CLAUDE.md's rule 4.
  - `preprocess_binary_data.py`'s `fit_joint_pspl()`: added `piE_N`/`piE_E`
    to the 6-param calibration fit. `t0_par` bootstrapped by reusing the
    same function with `delta_sN`/`delta_sE` fixed at zero rather than
    writing a separate plain-fit function -- those arrays multiply
    `piE_N`/`piE_E` in `trajectory()`, so zeroing them collapses the
    parallax terms regardless of `piE_N`/`piE_E`'s value, an exact
    pre-parallax special case (the same trick `lc_models.plain_flux()`
    already used).
  - `mcmc_fit_binary.py`'s `run_ogle_only_diagnostic()`: replaced its own
    plain-PSPL fit + hand-rolled plot with a direct call to
    `mcmc_fit.fit_parallax_pspl_mcmc()`/`get_t0_par()`/`plot_fit_lc()` --
    turned out to be the exact same single-instrument magnitude-space shape
    O-05-BLG086 already solves, so this was a pure reuse, not new code.
  - `mcmc_fit_binary.py`'s `run_joint_fit()`: genuinely new code (no
    existing analog works in magnification space) -- added `piE_N`/`piE_E`
    and Student-t `scale`/`dof` on top of the existing `(t0, u0, tE)`
    Nelder-Mead-then-emcee fit.
  - Real result: joint fit chi2/dof 1.51; `dof`~5.7-6.3 in both new O-03-BLG235
    fits, consistent with O-05-BLG086's own heavier-than-Gaussian finding.
- User then asked where the residual plots/histograms were for O-03-BLG235
  -- `zoom_utils.plot_fit_panels()` had never gotten session 8's
  residual-panel treatment (only `mcmc_fit.plot_fit_lc()` had). First pass:
  promoted the two per-panel helpers (`_plot_residual_panel`/
  `_plot_residual_hist`) out of `mcmc_fit.py` and into `zoom_utils.py` as
  public, multi-series-capable functions, used by both `plot_fit_lc()` and
  the newly-upgraded `plot_fit_panels()` -- also needed `plot_fit_panels()`'s
  own GridSpec rework to match `plot_fit_lc()`'s 4-row layout. This first
  pass had a real bug, caught by the user (not by inspection): those helpers
  show *standardized* residuals (`(obs-model)/err`), where every point's
  error bar is hardcoded to 1 by construction -- correct for a single
  instrument, but on a panel mixing OGLE and MOA it made OGLE's real ~4x
  better photometric precision (verified: median `A_err` 0.080 vs. 0.324)
  invisible, since both instruments' points got the same-size error bar
  regardless. Fixed by giving `plot_fit_panels()` its own raw
  (non-standardized) residual panel instead -- `obs-model` in physical
  magnification units, each point keeping its own instrument's real error
  bar -- rather than forcing it through the shared standardized-residual
  helpers, which stayed as-is for `plot_fit_lc()`'s single-instrument case
  (see CLAUDE.md's updated "Shared plotting/MCMC helpers" note). The joint
  fit's zoomed residual panel still visibly shows the real caustic-crossing
  bump as an outlying spike, exactly the "known-incomplete PSPL model"
  signature CLAUDE.md already described from the raw overlay alone.
- Hit and fixed a real robustness bug while re-running: `fit_parallax_pspl_mcmc()`'s
  `curve_fit` call was hitting scipy's default `maxfev=1600` cap (200*(N+1)
  for its 7 params) on the OGLE-only diagnostic's deliberately marginal fit
  -- reproducible standalone, not a one-off flake. Fixed with an explicit
  `maxfev=10000`; no effect on O-05-BLG086's already-converging fit.

### Learned & open questions
- The "reuse a function with a degenerate/zeroed input to recover a simpler
  special case" trick (first used for `plain_flux()`) generalizes cleanly
  one level up: the same move worked for an entire *calibration fit*
  (`fit_joint_pspl()` with `delta_s=0`), not just a single model-evaluation
  function. Avoided writing and maintaining a second near-duplicate chi2.
- A repo-shape lesson, corrected mid-session: `plot_fit_lc()` and
  `plot_fit_panels()`'s residual panels *looked* like the same shape
  (scatter + errorbar + histogram beside a fit overlay) but weren't --
  single- vs. two-instrument wasn't the relevant axis; standardized vs. raw
  residuals was, and that distinction only became visible once two
  differently-precise instruments shared one panel. A pattern match on
  visual shape alone (session 8's helpers "looked reusable") isn't the same
  check as "does the abstraction hide something instrument-specific" --
  CLAUDE.md's existing "don't force an abstraction over a real difference"
  note existed for exactly this, but the difference wasn't obvious until a
  human looked at the actual output, not just the code.
- O-03-BLG235's joint fit shows a real discrepancy between its Nelder-Mead
  point estimate and MCMC posterior median for `piE_N`/`piE_E` (point
  ~0.05/-0.11 vs. posterior median ~0.68/-0.39) -- plausible given PSPL is
  a known-incomplete model for this real 2L1S event (unlike O-05-BLG086,
  where point estimate and posterior agreed closely), but not yet
  root-caused. Not treated as a bug this session since the model is already
  known-incomplete here, but flagged rather than silently accepted.
- The `curve_fit` `maxfev` failure was deterministic and reproducible given
  fixed inputs (not a flaky/random failure) -- worth remembering before
  assuming a re-run will "just work" the second time.
- Doc/dependency cleanup, done this session rather than deferred (user
  asked to knock it out immediately once raised): `requirements.txt`
  dropped `sympy`/`mpmath`/`networkx` -- `sympy` is only needed by the
  standalone `scratch/derive_binary_quintic.py` (already documented as
  needing its own throwaway venv, not this project's), `mpmath` was only
  there as `sympy`'s own dependency, and `networkx` isn't imported by any
  script in this repo at all (likely an artifact of `graphify`'s own setup
  in this venv, not a pipeline dependency). `README.md`'s stale "no
  requirements.txt yet" setup line replaced with `pip install -r
  requirements.txt`, matching CLAUDE.md's own documented command. The
  superseded Desktop-level `graphify-out/` flagged last session turned out
  to already be gone -- nothing left to do there.

### Next session
- Confirmed via `/grill-me`: extend annual parallax (`piE_N`/`piE_E`) *and*
  the Student-t robust likelihood into the 2L1S track
  (`scratch/fit_2l1s.py` + `scratch/mcmc_fit_2l1s.py`), explicitly crossing
  the boundary CLAUDE.md previously drew ("not the `scratch/` 2L1S code, a
  separate track") -- both changes now apply everywhere.
  - Both `fit_2l1s.py`'s multi-start Nelder-Mead point estimate *and*
    `mcmc_fit_2l1s.py`'s MCMC get `piE_N`/`piE_E` (not just the MCMC seeded
    from a non-parallax point estimate) -- decided this way specifically
    because the existing multi-start step exists to avoid false minima, and
    seeding an MCMC that includes 2 new free parameters from a point
    estimate that never saw them risked starting in the wrong basin.
  - Applies to both the joint (OGLE+MOA) fit and the MOA-only control fit
    (`run_fit_moa_only()`/`run_mcmc_moa_only()`) -- MOA-only exists
    specifically as an independent check on the joint fit, so it needs to
    stay comparable once the joint fit's model changes.
  - Explicit decision: add parallax (+ Student-t) on top of the *current*
    2L1S search as-is, accepting the existing "near-miss cusp solution"
    risk for now (still only catchable via the image-count check, not chi2
    alone) rather than also trying to fix the search's robustness in the
    same session -- keeps the two problems (new physics vs. search
    reliability) from getting conflated if something looks wrong
    afterward. Search-robustness stays a separate, still-open problem.
  - The "MOA-2019-BLG-008 mag_err rescaling" and "resume the 2L1S search
    problem" backlog items were both considered and explicitly deferred in
    favor of this -- still backlog, not scheduled.

## 2026-09-23 — session 10

### Built
- Completed session 9's confirmed goal: extended annual parallax
  (`piE_N`/`piE_E`) + Student-t robust likelihood into the 2L1S scratch
  track (`scratch/fit_2l1s.py`, `scratch/mcmc_fit_2l1s.py`), crossing the
  boundary CLAUDE.md previously drew around it. `lc_models.binary_trajectory()`
  gained `piE_N`/`piE_E`/`delta_sN`/`delta_sE` params, mirroring
  `trajectory()`'s own `delta_tau`/`delta_beta` math exactly, no
  flag/default. `fit_2l1s.py`'s `chi2()` now wraps a new `residuals()`
  (extracted, returns the raw per-point standardized-residual array);
  `mcmc_fit_2l1s.py`'s real MCMC adds `scale`/`dof` and a Student-t
  log-likelihood mirroring `mcmc_fit.fit_parallax_pspl_mcmc()`'s existing
  pattern (Gaussian chi2 for the point-estimate stage only, Student-t for
  the real posterior).
- User wrote most of this by hand per CLAUDE.md's usual mode; Claude's role
  was reviewing each iteration against a live run and catching bugs before
  they silently corrupted results -- six were caught this way, not by
  inspection alone: `residuals()` returning a bare tuple for the joint
  branch (crashed immediately), a pre-squared scalar for the MOA-only
  branch (chi2 silently squared again), unpacking its own theta in a
  different param order than every other call site used (silently swapped
  `s`/`q` with `piE_N`/`piE_E`), `plot_fit()`'s caustic-panel trajectory
  reusing a stale, wrong-length `delta_sN`/`delta_sE` against a freshly
  resampled plotting grid (shape-mismatch crash, fixed twice -- missed one
  call site the first pass), the Student-t log-likelihood keeping the old
  Gaussian `-0.5*` scaling factor and missing the `-log(scale)` Jacobian
  term after the switch, and a stale 6-item print-label list silently
  mislabeling `piE_N`/`piE_E`'s printed values as `s`/`q` in console output
  (real `s`/`q` never printed at all -- only recoverable from the saved
  plot's caustic-panel title).
- Reshaped `fit_2l1s.py`'s `plot_fit()`: full multi-year baseline panel
  replaced with an "event season" window (HJD 2700-3000, matching
  `compare_2l1s_fits.py`'s existing convention), caustic geometry moved
  from its own subplot into a square inset
  (`mpl_toolkits.axes_grid1.inset_locator.inset_axes`, physical inches --
  `Axes.inset_axes()`'s fraction-based sizing doesn't give a square box
  against a non-square parent) inside the zoomed panel, and a raw-residual
  panel added beneath it (per-instrument real error bars, same convention
  as `zoom_utils.plot_fit_panels()`, hand-rolled here since that helper has
  no MOA-only mode). Gained a `tag` kwarg so a second fit to the same
  `use_ogle` mode doesn't overwrite the first's plot.
- Built a direct chi2-vs-Student-t MCMC comparison for the 2L1S track
  (`mcmc_fit_2l1s.run_mcmc_chi2()`/`log_probability_chi2()`, plain
  `-0.5*chi2`, no `scale`/`dof` -- a deliberate diagnostic, not wired into
  `__main__`) to make concrete a hypothesis raised mid-session: does chi2's
  unbounded quadratic penalty drag the whole fit toward one extreme point.
- Investigated the SLURM cluster path: recommended the CPU partition over
  GPU (the `Pool()` parallelism is per-walker OS processes each doing a
  small, sub-GPU-batch-sized calculation -- many processes sharing one GPU
  context would serialize, not speed up), drafted
  `scratch/submit_mcmc_2l1s.sbatch`, and fixed `run_mcmc()`'s `Pool()` to
  read `SLURM_CPUS_PER_TASK` explicitly rather than the default
  `cpu_count()` (which reads the physical node's core count, not what
  SLURM's cgroup actually allocated -- oversubscribes on a shared node).
- Found and fixed a real, independent bug while building a Bond et al.
  2004 comparison table: `preprocess_binary_data.py`'s hardcoded
  `COORDS = SkyCoord("18h05m16.35s -28d53m42.0s")` had the wrong RA --
  both raw OGLE and MOA files' own `\RA` header lines say `18h01m16.35s`,
  4 arcmin of RA (~1 degree on sky) off, silently corrupting every
  O-03-BLG235 parallax calculation made to date, this session included.
  Ran `/grill-me` to scope the fix rather than guess: replaced the
  hardcode with a new `preprocess_binary_data.load_coords()` that parses
  each raw file's own header and raises if OGLE's and MOA's disagree,
  deliberately scoped to O-03-BLG235 only (O-05-BLG086's `.dat` file has
  no header to parse at all; MOA-2019-BLG-008's raw file isn't downloaded
  yet, so scoping any further was guessing).
- Re-ran the full O-03-BLG235 chain (`preprocess_binary_data.py` ->
  `mcmc_fit_binary.py` -> `fit_2l1s.py` -> `mcmc_fit_2l1s.py`) with the
  corrected coordinates to get valid post-fix numbers.

### Learned & open questions
- The chi2-vs-Student-t comparison's result was concrete, not just
  theoretical: under chi2, the model curve visibly distorted into a
  spurious double-peak trying to reach one extreme real MOA
  caustic-crossing point -- and still missed it by ~8 magnification units
  even so -- while the residual panel's scale grew ~3x across the *rest*
  of the light curve as the cost of that chase. Under Student-t, the model
  didn't bend toward that point at all and fit everything else cleanly.
  Conclusion reached together: a single point neither likelihood can
  explain by parameter tuning alone (even chi2's maximal chase falls short)
  is more consistent with a single-epoch photometric anomaly or an
  unresolved finite-source/cusp effect than with something this trajectory
  model family can capture -- not worth chasing via the image-count check.
- The coordinate bug's fix shifted real numbers, not just cosmetics:
  joint MCMC's `alpha` moved from 74.7 deg (Bond: 223.8 deg -- way off) to
  217.5 deg (within 6 deg) after the fix; `piE_N`/`piE_E` dropped from
  being pinned at the prior's +/-2 boundary to modest `-0.109`/`-0.024`.
  Plausible reading: the boundary-pinned parallax values seen
  pre-fix were the wrong coordinates forcing the fit to explain away
  geometric error as spurious parallax signal, not a real detection --
  though this isn't proven, just newly plausible.
- `alpha`'s posterior bimodality was checked quantitatively (KDE
  peak-finding on the saved chains, not eyeballing the corner plot) and
  turned out to be ~78-108 deg apart across both joint and MOA-only
  chains -- not the ~180 deg mirror-image degeneracy assumed at a glance.
  More consistent with genuinely separate local chi2 minima in the
  still-unresolved search landscape (CLAUDE.md's "near-miss cusp" problem)
  than with one clean, well-understood symmetry.
- The MOA-only control fit still disagrees substantially with the joint
  fit on topology even post-fix (`s=1.098`/`q=0.322` vs.
  `s=1.537`/`q=0.0079`) -- predates and is unrelated to the coordinate
  bug, still open.
- Two bugs this session were caught only by actually running the code and
  reading its output against expectations (a saved plot's title, a
  Bond-comparison table) -- not by code review, and not close calls either
  time. Reinforces this repo's own stated convention (CLAUDE.md, "Setup and
  commands": correctness "judged by inspecting the plot/printed fit it
  produces") over static inspection alone.
- Known naming gap, not fixed this session (see Next session): `run_fit()`
  and `run_mcmc()` both default to `plot_fit(..., tag="")`, so whichever
  ran most recently silently overwrites the other's output at
  `scratch/O-03-BLG235_2l1s.png`. Pre-dates this session; became concrete
  while running the chi2-vs-Student-t comparison side by side.

### Next session
- Confirmed via `/grill-me`, scoped to the 2L1S track only (single-lens
  PSPL stays on Student-t as-is -- already validated, and there's no
  physical reason for a spatially-varying statistic there):
  1. **File naming cleanup**: drop `plot_fit()`'s `tag=""` default (make it
     required) and give every 2L1S call site a method-identifying tag
     (e.g. `_nelder_mead`, `_mcmc_studentt`, `_mcmc_chi2`) so different
     fitting methods stop silently overwriting each other's plots/corner
     plots/chain files.
  2. **Better-statistic exploration**: try Huber loss first -- a
     well-documented, widely-used approach -- before anything custom. A
     "joint statistic combining Student-t and chi2" idea was raised too,
     explicitly deprioritized until Huber (and other established options)
     have been tried.

## 2026-09-24 — session 11

### Built
- **File-naming fix** (closes session 10 goal #1): `fit_2l1s.py`'s
  `plot_fit()` dropped its `tag=""` default -- keyword-only, required now.
  Every 2L1S call site passes its own method-identifying tag:
  `run_fit()` -> `"_nelder_mead"`, `run_mcmc()` -> `"_mcmc_studentt"`
  (`run_mcmc_chi2()`'s `"_chi2"` was already distinct pre-session-11, left
  as is), `run_mcmc_huber()` -> `"_huber"` (see below). No more silent
  overwrite between methods' plots.
- **`TwoL1SParams` refactor** (landed mid-session, in parallel): a
  `NamedTuple` holding the 8 physical params in one canonical order.
  `seeds`, `residuals()`, `_fit_one_seed()`, `run_fit()`, and `plot_fit()`
  in `fit_2l1s.py` all unpack through it now instead of each re-deriving
  the order by hand -- directly closes the bug class session 10 hit
  (`residuals()` silently swapping `s`/`q` with `piE_N`/`piE_E`).
  `mcmc_fit_2l1s.py`'s `LABELS`/`LABELS_CHI2`/`LABELS_HUBER` are now
  derived from `TwoL1SParams._fields` rather than three hand-typed lists
  that could drift apart. `residuals()` also now always returns an
  ndarray (`np.inf`-filled on an unphysical/degenerate trial) instead of
  sometimes a bare scalar -- callers simplified their `np.isscalar(resids)
  or ...` guard down to just the finite check. The three near-identical
  likelihoods' bound-checks and walker-seeding were deduplicated into
  `_physical_log_prior()`/`_sample_physical_prior()`.
- **Better-statistic exploration** (closes session 10 goal #2, Huber tried
  first as planned): implemented a Huber-loss likelihood for the 2L1S
  track (`huber()`, `log_prior_huber()`, `log_probability_huber()`,
  `sample_prior_huber()`, `run_mcmc_huber()`). `DELTA=1.345` fixed, not
  fit, so it needs no normalizing-constant correction of its own; `scale`
  stays free with the same `-log(scale)` Jacobian term `log_probability()`
  already has. Six real debugging passes went into getting this right,
  including a missing `delta` argument to `huber()` and a `sample_prior`
  that didn't include a `scale` column for the sampler's dimensionality.
- Fixed the caustic-geometry inset plot not actually rendering square:
  `inset_axes(width=1.8, height=1.8, ...)` requested a square box, but
  `set_aspect("equal")`'s default `adjustable="box"` was silently
  reshaping that box to match the data's own aspect ratio. Fixed with
  `set_aspect("equal", adjustable="datalim")`.
- Fixed docstring/style inconsistencies across `mcmc_fit_2l1s.py` (missing
  space in tuple-unpack lines, trailing whitespace, inconsistent blank-line
  spacing between defs, a stale copy-pasted docstring on `log_prior_huber`,
  missing docstrings on `log_probability()`/`log_probability_huber()`/
  `sample_prior_chi2()`).
- Cleaned `scratch/`: removed 3 files with zero code references, all last
  touched in the very first pre-parallax/pre-weighting 2L1S commit
  (`O-03-BLG235_2l1s_refined_theta.npy`, `fit_2l1s_run.log`,
  `mcmc_fit_2l1s_run.log`), plus 3 stale plot outputs that no current call
  site produces filenames for any more post-naming-fix
  (`O-03-BLG235_2l1s.png`, `O-03-BLG235_2l1s_moa_only.png`,
  `O-03-BLG235_2l1s_studentt.png`), plus `__pycache__`. Deliberately kept
  `compare_2l1s_fits.py`/`fit_2l1s_moa19008.py` (documented as
  intentionally frozen/incomplete, not dead) and the
  `cross_check_mulensmodel*.py`/`derive_binary_quintic.py` provenance
  scripts.
- Widened the search to reach more extreme topologies: `mcmc_fit_2l1s.py`'s
  `S_RANGE` (0.7-1.8 -> 0.5-2.0) and `LOG_Q_RANGE` (1e-4-1.0 -> 1e-5-1.0)
  priors, and `fit_2l1s.py`'s Nelder-Mead multi-start seed grid. The seed
  grid went through two revisions: first expanded to 6x4x16=384 seeds,
  which on this machine's 6 physical cores (`ProcessPoolExecutor` defaults
  to `os.cpu_count()` workers) was a real ~3x wall-clock jump; a two-stage
  coarse-then-refine search was implemented as a fix, then deliberately
  reverted in favor of simply widening the grid spacing and cutting seed
  count back down -- single-stage, 96 seeds (`s_list=[0.5,1.0,1.5,2.0]`,
  `q_list=[0.001,0.03,1.0]`, 8 `alpha` values), same overall range as the
  384-seed attempt, just sparser.
- Ran all six joint/MOA-only x chi2/Student-t/Huber combinations for real,
  post-refactor -- the first time every one of them has had a genuine
  end-to-end run under current code (previously only unit-level checks
  and one earlier joint-Huber run existed). All six completed without
  exceptions.

### Learned & open questions
- **All three statistics (chi2, Student-t, Huber) fail to fit the obvious
  caustic spike.** None actually captures the sharp caustic-crossing
  feature -- they just differ in how much they tolerate or chase it as an
  apparent outlier. This effectively answers (and closes, for now) the
  three-way "which loss function is best" comparison session 10 planned:
  the answer is none of them, so it's not currently a likelihood-shape
  question.
- **More fundamental finding**: the six runs' fitted parameters disagree
  with each other far more than a likelihood-shape difference alone should
  cause -- joint `alpha` ranges 6.60 rad (Nelder-Mead) -> 2.96 rad
  (Student-t) -> 1.34 rad (chi2 MCMC), with `s`/`q` moving comparably.
  Consistent with the older, still-open "near-miss cusp" search-landscape
  problem (session 2's BIC-vs-PSPL verdict, blocked since before this
  session) being the dominant issue right now, not the loss function.
- Two smaller loose ends, flagged but not investigated: MOA-only
  Nelder-Mead's `piE_N` landed exactly on `PIE_RANGE`'s `2.0` boundary (a
  prior-boundary pin, not a converged interior value); Huber MOA-only's
  MCMC run took 24:33 at only 43% CPU vs. 2-10 min at 400-500% for every
  other run this session -- looks like a resource-contention artifact, not
  genuinely more computation.
- Statistics alternatives discussed but not implemented: whether Poisson
  statistics would suit the heavy tails better than Huber/Student-t
  (concluded no -- Poisson-vs-Gaussian is about noise shape at low photon
  counts, not applicable to this bright-bulge photometry, and not
  implementable anyway without raw counts, which OGLE/MOA's published
  magnitudes/relative fluxes don't provide); a Sivia & Skilling
  good-and-bad-data mixture likelihood; an additive jitter term
  (`sigma_eff^2 = sigma_reported^2 + sigma_jitter^2`) as a physically-closer
  alternative to the current multiplicative `scale` for systematics that
  don't scale with the reported error. Finite-source effects (2L1S+FS)
  were also raised as the physical explanation for why *no* point-source
  loss function can consistently fit a caustic spike -- not evaluated
  against the alternatives above yet.
- Two full `/grill-me` rounds were run this session trying to scope the
  statistic-comparison work methodically; the user found the back-and-forth
  unproductive once the real run batch had already answered the practical
  question. Worth remembering: once actual results are in hand, re-deriving
  next steps from first principles via further questioning has diminishing
  returns -- a direct decision from the person who just saw the results is
  faster and no less valid than another interview round.

### Next session
- Confirmed directly (not via further `/grill-me` -- see note above):
  explore rescaling the reported per-point error bars (e.g. an additive
  jitter term, per "Learned" above, or another recalibration approach) and
  re-fit with plain Gaussian chi2, rather than continuing to tune the loss
  function's shape.
- Huber loss is explicitly parked -- not considered the right next step
  for the caustic-spike problem.
- Not scheduled, but still open and unresolved: the near-miss-cusp
  search-landscape/basin-disagreement question, the `piE_N` boundary pin,
  Huber MOA-only's timing anomaly, and the finite-source-effects idea.

### Planned roadmap (post-session addendum, papers + ordered next steps)
User-supplied reading list and implementation order for the 2L1S track,
recorded ahead of starting the work so intent is on record before code
changes begin. Order as given, error-bar rescaling first (already the
confirmed "Next session" item above):
1. **Error bar rescaling** (in progress next -- see "Next session" above).
2. **Finite source effects (`rho*`)** -- extend the point-source 2L1S model
   to account for finite source size, the deferred step `lc_models.py`'s
   own docstring already flags (see "O-03-BLG235's PSPL-vs-2L1S model
   comparison" above) and the leading candidate explanation raised in
   session 11's "Learned" section for why no point-source loss function
   fits the caustic spike.
3. **Caustic parametrisation** -- explore Cassan (2008)'s curvilinear
   abscissa parametrisation of the caustic curve, as groundwork for the
   genetic-algorithm step below.
4. **(d, q) grid search** -- grid over separation/mass-ratio, chi2-fit
   caustic entry/exit times, entry/exit points, and `rho*` at each grid
   point to determine caustic geometry -- directly targets the still-open
   near-miss-cusp/basin-disagreement problem (session 2 onward) by
   replacing multi-start Nelder-Mead's ad hoc seed grid with a
   geometry-driven search.
5. **Genetic algorithm** (Charbonneau 1995) -- explore entry/exit points on
   the parametrised caustic and entry/exit times, seeded/checked by eye;
   MCMC refinement afterward, following the same point-estimate-then-MCMC
   pattern already used elsewhere in this repo (see "Annual parallax +
   robust likelihood" above).
6. **Anomaly vs. outlier test** (added 2026-09-25, session 13; order
   relative to steps 1-5 not yet decided) -- adopt Dominik et al. (2007)'s
   SIGNALMEN criterion to decide whether a deviation from the ordinary
   (PSPL) model is real signal or a bad point: an anomaly needs >=5 recent
   points deviating to the same side, >=3 of them significantly
   (`DEV_SIG=2` sigma, scatter-adjusted), no contradicting non-deviant
   points in between, and no "zig-zag" (more than one change of direction
   => likely photometry failure). A lone deviant point is only a "check",
   not an anomaly. First use: the single MOA point at HJD~2843 that no
   2L1S fit reaches -- is it signal or an outlier? (vs. the persistent
   departure from HJD~2835, which already looks anomaly-shaped). Test
   against the ordinary model only -- SIGNALMEN's bisquare down-weighting
   rejects exactly the large residuals a caustic fit needs, so it isn't a
   2L1S likelihood. Formulas/defaults were read via a summarised fetch of
   the arXiv HTML; verify against the PDF before implementing.

**References**:
- Charbonneau (1995) -- genetic algorithm for fitting (step 5).
- Cassan (2008) -- curvilinear abscissa parametrisation of the caustic
  curve (step 3).
- Kains (2009) -- a real dataset applying the Cassan (2008)/genetic-
  algorithm methodology end to end; reference implementation to check
  this repo's own approach against once steps 3-5 are underway.
- Dominik et al. (2007), MNRAS 380, 792 (arXiv:0706.2566) -- SIGNALMEN
  anomaly detector (step 6).

### Long-term pipeline direction (post-session addendum, via `/grill-me`)
User described wanting this repo to eventually become a real pipeline, not
just a two-event playground -- a genuine philosophy shift not previously
recorded anywhere, so run through a full `/grill-me` round for consensus per
CLAUDE.md rule 3 before any code changed. No code changed this round; this
is a direction decision only.
- **Vision**: config-driven, not fully automatic. A new event = a config
  (raw file paths, coordinates, instrument format if it matches the existing
  OGLE/MOA-style parsing, model choice, initial guesses) rather than a
  system that auto-detects formats, auto-picks PSPL vs. 2L1S, or
  auto-generates initial guesses -- explicitly ruled out as a much bigger,
  more open-ended project than the config-driven version.
- **Sequencing**: this is a later, separate goal. The physics/logistics
  roadmap logged above (error-bar rescaling -> finite source effects ->
  caustic parametrisation -> (d,q) grid search -> genetic algorithm)
  continues uninterrupted; pipeline conversion doesn't start until that
  work is further along -- generalizing plumbing around a model that
  doesn't yet fit the one hard case in hand (O-03-BLG235's caustic) risks
  building the wrong abstraction before it's clear what "general" needs to
  support.
- **Output bar**: whatever the pipeline eventually produces per event must
  match the full treatment the two current datasets already get -- MCMC
  posterior, corner plots, best-fit light curve overlay, derived physical
  quantities (mass distribution, blend fraction, etc.) -- not a reduced
  version.
- **Claude's role unchanged**: still a guide, not an implementer, per
  CLAUDE.md rule 6 -- but with a standing consideration added going
  forward: favor implementation patterns that stay convertible later
  (parameterized functions over hardcoded per-event constants scattered
  inline, proper module structure) over anything better suited to
  notebook-style exploration, without starting pipeline work itself yet.
- A corresponding note was drafted for CLAUDE.md's "What this is" section
  (not applied without the user's explicit sign-off on wording, per its own
  rule 1).

## 2026-09-24 — session 12

### Built
- Diagnosed why every current 2L1S fit method (Nelder-Mead, Student-t/chi2/
  Huber MCMC) fails to produce a trustworthy fit, by actually looking at the
  fit-overlay plots rather than reasoning from chi2 numbers alone:
  Nelder-Mead's point estimate produces a razor-thin, physically implausible
  A~27 spike chasing a single outlier MOA point at HJD~2843, not a real
  caustic crossing. The three MCMC runs each look locally plausible (small
  residuals, a real double-hump structure) but land on substantially
  different (s,q) geometries (Student-t: s=1.64,q=0.106; chi2:
  s=0.97,q=0.0034; Huber: s=1.64,q=0.085) and none of them reach up to match
  that same extreme point -- confirming this is the search-landscape/
  near-miss-cusp problem (open since session 2), not something error-bar
  rescaling could fix (rescaling only reweights residuals at a given
  trajectory, it doesn't change which basin the search finds).
- User identified, by eye, that the caustic entry point is around HJD~2835
  (where the data persistently departs above the single-lens model) -- a
  concrete anchor the current (t0,u0,tE,alpha,s,q) parametrisation has no
  direct way to use, since entry/exit time is only an indirect, emergent
  property of that parametrisation rather than a fittable quantity. Decided
  directly (not via `/grill-me`, given results were already in hand) to
  implement Cassan (2008)'s curvilinear-abscissa caustic parametrisation
  next -- ahead of session 11's error-bar-rescaling roadmap item, since
  rescaling only matters once there's a trustworthy trajectory to rescale
  around.
- Reorganized `scratch/`'s 2L1S outputs into
  `scratch/2l1s/{nelder_mead,chi2,huber,studentt,stylecheck}/` (previously
  ~25 files flat in `scratch/`), and updated `fit_2l1s.py`'s `plot_fit()`
  and `mcmc_fit_2l1s.py`'s three `run_mcmc*()` functions to write into those
  subfolders going forward (tag -> folder mapping in `plot_fit()`, `KeyError`
  on an unregistered tag -- same "no silent default" reasoning as `tag`
  itself having no default). Already committed by the user (`d369ff1`,
  "added scratch folder layers").
- Started implementing the Cassan parametrisation's prerequisites in
  `lc_models.py` (uncommitted, still in progress):
  - `caustic_curve()` now orders its 4 roots-per-`phi` into continuous
    branches via nearest-neighbor continuation
    (`scipy.optimize.linear_sum_assignment` on the pairwise distance matrix
    between consecutive `phi` steps) -- previously each `phi` step's 4 roots
    were returned in whatever arbitrary order `np.roots()` gave, fine for
    the existing scatter-plot-only caller but meaningless for defining arc
    length along the curve. Also added a closed-curve check (`phi` now
    sampled with `endpoint=True` specifically so `phi=0`/`phi=2*pi` can be
    compared; matched via the same assignment approach, compared with a
    tolerance rather than `==`, since two independent `np.roots()` calls at
    mathematically-identical `phi` values won't be bit-identical).
  - `equidistant_caustic()`: resamples each of the 4 ordered branches to be
    evenly spaced *in arc length* (not `phi`) -- cumulative arc length via
    `np.diff`/`np.cumsum` along the `phi` axis, normalized per-branch
    (deliberately not sharing one normalization constant across branches),
    then interpolated (`np.interp` on real/imaginary parts separately,
    recombined) at a uniform `zeta` grid. This is a numerical answer to the
    concern that a naive chord-based method under/over-resolves
    high-curvature regions like cusps -- raised when the user pointed out
    Cassan's own paper derives an analytic `d(zeta)/d(phi)` via implicit
    differentiation for exactly this reason. Worked out the analytic route
    too (`d(zeta)/d(phi) = dz/dphi + e^(-i*phi) * d(zbar)/dphi`, following
    directly from conjugating the critical curve equation already in this
    file's docstring -- no need for `m1/m2/z1/z2` in that step at all) but
    decided to go numerical instead, deferring the analytic derivative
    rather than implementing it now.

### Learned & open questions
- Confirmed the failure mode across all six existing 2L1S fits is real, not
  a rescaling-fixable likelihood-shape issue -- visual inspection of the
  fit-overlay plots (not just the chi2/parameter numbers already in session
  11's notes) was needed to see this clearly, particularly Nelder-Mead's
  unphysical spike, which isn't visible from its `chi2` value alone.
- Two validation gaps in the new `caustic_curve()`/`equidistant_caustic()`
  code, not yet checked: (1) whether the closed-curve check actually passes
  for all 4 branches at a real `(s,q)` (or whether some need the
  cross-branch permutation-stitching case flagged during development, which
  isn't implemented); (2) whether `equidistant_caustic()`'s output distances
  are actually uniform (the sanity check suggested wasn't run before the
  session ended).
- Deliberately deferred: the analytic `d(zeta)/d(phi)` derivative (formula
  above) in favor of the numerical resampling approach, on the reasoning
  that the numerical route is simpler to get right immediately and should
  already resolve the cusp-clustering concern, given the underlying `phi`
  sampling (n_phi=600) is dense enough. Revisit the analytic route only if
  the numerical resampling turns out not to be accurate enough in practice.
- CLAUDE.md's "O-03-BLG235's PSPL-vs-2L1S model comparison" section
  (documenting `lc_models.py`'s binary-lens functions) was deliberately NOT
  updated to describe `caustic_curve()`'s new behavior or
  `equidistant_caustic()` this session -- both are still unvalidated
  (see gaps above), and documenting unvalidated behavior risks committing to
  a description that changes once actually checked. Fold into CLAUDE.md once
  validated, not before.

### Next session
- Not yet confirmed via `/grill-me` -- session ended on a time constraint,
  picking up exactly where this one stopped rather than re-opening the
  direction question:
  1. Run the two validation checks flagged above (closed-curve check passes
     per-branch; equidistant output is actually evenly spaced) before
     trusting either function further.
  2. Pick which of the (possibly multiple, topology-dependent) disjoint
     caustic branches is the physically relevant one for an existing
     solution (e.g. Student-t's `s=1.641, q=0.106`), cross-checked against
     that solution's own caustic-inset plot.
  3. Build the actual reparametrisation function: `(zeta_entry, zeta_exit,
     t_entry, t_exit)` + fixed `(s,q)` -> `(t0, u0, tE, alpha)`, via the two
     caustic points' implied velocity.
  4. Validate in `scratch/cassan_caustic.py` (created this session, still
     empty): fix `(s,q)`, pick `zeta_entry`/`zeta_exit` by eye,
     `t_entry~=2835` (the user's visual anchor) plus a guessed `t_exit`, run
     one Nelder-Mead refinement, sanity-check the resulting light curve
     before building the full `(d,q)` grid search on top.
- `/graphify --update` not run this session -- the new caustic-geometry code
  is still mid-implementation/unvalidated, so deferred until this track is
  further along rather than rebuilding the graph around unfinished code.

## 2026-09-28 — session 13

### Built
- **Cassan (2008) caustic parametrisation, end to end** (session 11
  roadmap item 3; session 12's "Next session" steps 1-4). Most of the
  code and every output under `scratch/2l1s/cassan/` date from 2026-09-25
  (file timestamps); committed 2026-09-28 in `a261e54` ("cassan mcmc
  implemented") on top of session 12's `a2eb25d` ("mid-update on the
  cassan caustic implementation"). Only the SIGNALMEN roadmap addendum
  below was logged at the time -- this entry backfills the rest.
- `lc_models.py`:
  - `caustic_curve()` now returns a **list of closed caustics** (1
    resonant, 2 wide, 3 close) instead of a raw `(n_phi, 4)` array. The
    critical-curve quartic is solved in Cassan's frame (`m1` at 0, `m2` at
    `-s`) and mapped back to `binary_trajectory()`'s frame via
    `-conj(z) + z1`. Per-`phi` root ordering keeps session 12's
    `linear_sum_assignment` continuation; the 4 resulting phi-pieces are
    then stitched into closed curves by matching each piece's `phi=2pi`
    end to some piece's `phi=0` start (same assignment trick), and it
    **raises `ValueError`** if any join is off by more than `tol` (1e-5)
    instead of printing a check. This closes session 12's validation gap
    #1: a non-closing caustic is now an error, not a silent print, and the
    cross-piece stitching that session flagged as unimplemented is exactly
    what the new code does.
  - `equidistant_caustic()` rewritten to take **one** closed curve (an
    element of `caustic_curve()`'s list) and resample it evenly in arc
    length *including the closing segment* back to the first point --
    the old version worked on all 4 branches at once, each open-ended.
  - New: `cassan_caustic(s, q, idx=0)` (one caustic, rolled to start at
    its rightmost point, forced counter-clockwise via the shoelace sign,
    then equidistantly resampled -- so abscissa `sigma` in [0, 1) is
    Cassan's curvilinear abscissa rescaled, with an origin/direction that
    don't depend on `np.roots`' ordering); `caustic_point(caustic, sigma)`
    (periodic interpolation); `cassan_to_standard()` (`sigma_in`,
    `sigma_out`, `t_in`, `t_out` -> `t0`, `u0`, `tE`, `alpha` for the
    straight, parallax-free trajectory through both caustic points -- every
    trial is a genuine crossing by construction); `standard_to_cassan()`
    (inverse: every `(sigma, t)` where a given trajectory crosses the
    caustic, time-sorted, empty on a miss).
  - `torch.set_num_threads(1)` at import -- stops each spawned emcee
    worker also spinning up a full torch thread pool.
  - Uncommitted (the only `git diff`): dropped the now-unused `SkyCoord`
    import.
- `scratch/cassan_caustic.py` (empty at end of session 12, now 152
  lines): `CassanParams` NamedTuple `(sigma_in, sigma_out, t_in, t_out,
  s, q)`; `__main__` first round-trips Bond et al. 2004's published
  solution standard -> Cassan -> standard and prints both chi2s, then
  runs `fit()` from two starting `(s, q)` -- Bond's (1.120, 0.0039) and
  session 11's Student-t basin (1.641, 0.106) -- each a 20x20 `sigma`
  grid x 5 `t_out` values at `t_in=2835` (the by-eye anchor; `t_out`
  gridded, not assumed), Nelder-Mead on the best 5 grid points with
  `(s, q)` fixed, then a final 6-param Nelder-Mead with `(s, q)` free.
  The better of the two seeds `run_mcmc()`: emcee, 32 walkers x 3000
  steps, spawn-context pool, Gaussian `-0.5*chi2` under flat priors
  (`t_in < t_out` inside `T_WINDOW` 2820-2870, `S_RANGE`/`LOG_Q_RANGE`/
  `TE_RANGE` reused from `mcmc_fit_2l1s.py`, resonant topology only).
  Chain saved to `scratch/2l1s/cassan/O-03-BLG235_2l1s_cassan_mcmc_chain.npz`.
- `scratch/fit_2l1s.py`'s `plot_fit()`: gained a rotated per-instrument
  residual histogram (`ax_hist`, sharing the residual panel's y-axis,
  zoom window only; figure widened to 10.5in), and three new tag ->
  folder entries (`_cassan_from_bond_sq`, `_cassan_from_studentt_sq`,
  `_cassan_mcmc` -> `scratch/2l1s/cassan/`). It and
  `compare_2l1s_fits.py`/`fit_2l1s_moa19008.py` now
  `np.concatenate(caustic_curve(...))` for their scatter-only insets.
- `requirements.txt`: `tqdm==4.70.1` (emcee's `progress=True`).
- Roadmap addendum: item 6, the SIGNALMEN anomaly-vs-outlier test
  (Dominik et al. 2007), was added to session 11's "Planned roadmap" on
  2026-09-25, during this same session -- the "session 13" label there
  refers to this entry.

### Learned & open questions
- **First Cassan-parametrised result is a genuine resonant-caustic
  crossing near Bond et al.'s solution.** MCMC best sample: `s`=1.108,
  `q`=0.0079 (Bond: 1.120, 0.0039 -- `s` close, `q` ~2x higher), which
  maps to `t0`=2847.72, `u0`=0.135, `tE`=60.9 d, `alpha`~214.6deg (Bond:
  2848.06, 0.133, 61.5, 223.8deg). The `alpha` gap to Bond narrows from
  session 6's ~15deg to ~9deg. Best-sample chi2 = -2*max(log_probs) =
  **1765.98**. Posterior (16/50/84 percentiles, 4800 samples after
  discarding 750 steps and thinning by 15):
  `sigma_in` 0.1064 +0.0120/-0.0100, `sigma_out` 0.7370 +0.0073/-0.0068,
  `t_in` 2832.61 +0.61/-0.29, `t_out` 2842.075 +0.012/-0.008,
  `s` 1.1096 +0.0078/-0.0069, `q` 0.00775 +0.00102/-0.00096.
- The Bond-seeded Nelder-Mead fit (`_cassan_from_bond_sq`: s=1.113,
  q=0.0087) is essentially the same solution the MCMC then explored; the
  Student-t-seeded fit (`_cassan_from_studentt_sq`) stayed near its old
  geometry (s=1.640, q=0.0772), grazing one end of an elongated
  caustic, with large structured residuals -- a clearly worse fit.
  Neither run's printed chi2 was saved, and neither was the Bond round-trip
  check's output, so the round trip is implemented but its result isn't on
  record.
- **The caustic-exit spike lands on the single MOA point at HJD~2842**
  (`t_out` = 2842.075) that every earlier 2L1S fit treated as an outlier
  (session 12: Nelder-Mead's implausible A~27 spike chased it, the three
  MCMC runs all missed it). Under this parametrisation it reads as the
  exit caustic crossing, not a bad point -- directly relevant to roadmap
  item 6's stated first use case ("is the HJD~2843 point signal or an
  outlier?"), which this result now argues is signal. SIGNALMEN would
  still be a useful independent check, but the question it was meant to
  settle first has a strong model-based answer already.
- The model's caustic **entry is at ~2832.6**, inside a data gap, slightly
  earlier than the user's by-eye anchor of ~2835 (which was only ever
  `fit()`'s grid starting point, not a constraint). The entry spike is
  therefore unconstrained by any point sitting on it.
- **Point-source spikes are razor-thin** (both crossings in the fit plot
  are near-vertical lines peaking far above any data) -- physically
  implausible for a real source of finite size, and the natural next thing
  to fix.
- **chi2 comparison with session 6's best**: 1765.98 vs. the previous best
  validated genuine crossing, chi2=1825.35 (session 6, joint MCMC best
  sample tight-refined). Comparable on paper: both are the same
  flux-profiled Gaussian chi2 (`fs_ogle`/`fb_ogle`/`fs_moa` solved per
  trial, introduced in session 6), both joint OGLE+MOA against raw flux,
  both parallax-free with 6 free trajectory params (session 6 predates
  parallax; session 10's RA fix only touches parallax terms, which vanish
  at `piE`=0). On that basis this is ~59 lower, and ~187 below Bond's own
  1952.90. Not re-evaluated under current code this session, so the
  1825.35 figure is from the session 6 log, not a fresh recomputation.
- **Session 12's validation gap #2** (whether `equidistant_caustic()`'s
  output is actually evenly spaced): no numerical check found anywhere --
  `scratch/scratchpad.ipynb` only plots resampled caustics as `+` markers
  (including an `ipywidgets` (s, q) explorer), which is at best a visual
  check. **Still unverified.**
- Known limitations of the current Cassan fit:
  - **Parallax-free** (`piE_N`=`piE_E`=0): the two-point -> straight-line
    mapping is only exact without parallax.
  - `cassan_caustic()` picks the caustic by **list index** (`idx=0`):
    unambiguous at fixed topology, but the index can point at a different
    caustic once `(s, q)` crosses a topology boundary -- hence the MCMC
    prior's resonant-only restriction (`len(caustic_curve(...)) == 1`,
    with `ValueError` at a boundary also rejected).
  - The MCMC starts from a **tight ball around the Nelder-Mead best**
    (1e-3 in `sigma`/`s`, 0.05 d in `t`, 1% in `q`), so it explores that
    one solution's posterior -- not a global search, and says nothing
    about other crossing geometries.
  - Likelihood is **Gaussian chi2**, not Student-t -- no `scale`/`dof`.
- `caustic_curve()`'s docstring still describes the old formulation
  (`m1/(z-z1)^2 + m2/(z-z2)^2`, coefficients via `np.convolve`); the code
  now solves the Cassan-frame quartic directly. Stale, not yet fixed.

### Next session
- Roadmap reordered by the user (Cassan parametrisation, old item 3, is
  done):
  1. **Finite-source effects (`rho*`)** -- starting now; motivated
     directly by the razor-thin point-source spikes above.
  2. **Error-bar rescaling** -- pushed to second.
  3. The rest unchanged: **(d, q) grid search**, then the **genetic
     algorithm** (Charbonneau 1995). The SIGNALMEN anomaly test's
     position is still undecided.
- For the user to fold in at "save state" (not fixed here): CLAUDE.md
  doesn't yet document the Cassan functions in `lc_models.py` or
  `scratch/cassan_caustic.py` (deliberately deferred in session 12 until
  validated -- gap #1 is now resolved, gap #2 still isn't); CLAUDE.md lists
  a `scratch/submit_mcmc_2l1s.sbatch` that doesn't exist on disk; and
  `scratch/scratchpad.ipynb` exists but isn't listed.

## 2026-09-28 — session 14

### Built
- Session 13's CHANGELOG entry (written this session, by a subagent, from the
  commits and saved outputs). Roadmap reordered by the user: **finite source
  (`rho*`) first, error-bar rescaling second**, the rest unchanged.
- **Finite-source 2L1S magnification**, `lc_models.binary_magnification_fs(zeta,
  s, q, rho, n=4000, gate=5)`, written by the user with guidance. Went through
  three sampling designs:
  1. Equal-area rings at edge radii `sqrt(k/n)` -- corrected to midpoint
     radii `rho*sqrt((k+0.5)/n)` (samples on the edges bias outward and leave
     the centre unsampled), and `rho` confirmed as the source *radius*, not
     diameter (literature convention; a factor-2 error would break both the
     analytic and MulensModel checks).
  2. Rings with a per-ring angular stagger -- worked, but a sample-point plot
     showed inner rings grossly oversampled (same `n_theta` on every ring)
     and a hole inside the innermost ring.
  3. **Fibonacci/sunflower spiral** (`r_i = rho*sqrt((i+0.5)/n)`, `theta_i =
     i * golden angle`): one radius per point, uniform density, no central
     hole, no spoke alignment. Kept.
- **Gating**: plain point-source everywhere, disk average only for points
  within `gate*rho` of any caustic point (min distance to
  `np.concatenate([equidistant_caustic(c) for c in caustic_curve(s, q)])`,
  caustic spacing ~0.29 rho at `n_phi=1000` -- `caustic_curve`/
  `cassan_caustic` defaults raised from 600 to 1000). Final filter is
  `A[mask] = _disk_average(zeta_arr[mask], ...)`, deliberately not
  `np.where` (which would evaluate the expensive branch for every point).
- **Threaded `rho` through the 2L1S track**: `TwoL1SParams` gains `rho` (9
  fields); `fit_2l1s.residuals()`/`plot_fit()` switch to
  `binary_magnification_fs` and reject `rho <= 0`; seeds use `rho=1e-3`
  (not gridded -- the seed grid is for topology, not source size).
  `mcmc_fit_2l1s.py`: `RHO_RANGE=(1e-5, 0.1)`, `rho` bounded in
  `_physical_log_prior` and drawn log-uniformly in `_sample_physical_prior`,
  and every hardcoded `theta[:8]`/`theta[8]`/`theta[9]` replaced by `N_PHYS =
  len(TwoL1SParams._fields)`. `cassan_caustic.py`: `CassanParams` gains
  `rho`; `fit(s, q, rho, ...)` runs its grid + fixed-`(s, q)` stages at
  `GRID_RHO=1e-6` (~point source, 0.06 s/chi2) and frees `(s, q, rho)` only
  in the final Nelder-Mead (finite source ~0.8 s/chi2); MCMC ball and prior
  gain `rho`; Bond check uses Bond's own `rho=0.00096`.
- `scratch/scratchpad.ipynb`: a cell plotting every disk sample
  `binary_magnification_fs` evaluates, coloured by `log10(A_PS,i / A_FS)`
  (captures `binary_magnification`'s real inputs/outputs by temporarily
  wrapping it, rather than re-deriving the grid).
- Ran both `cassan_caustic.fit()` starts with finite source (MCMC stage not
  run). Outputs overwrite `scratch/2l1s/cassan/*_from_{bond,studentt}_sq.png`
  (point-source versions recoverable from `a261e54`).

### Learned & open questions
- **Bond et al. 2004's Table 1, read from the arXiv PDF** (not previously
  recorded in this repo): best fit q=0.0039(+11/-7), **rho=0.00096(11)**,
  s=1.120(7), alpha=223.8(1.4) deg, u0=0.133(3), t0=2848.06(13),
  tE=61.5(1.8) d, chi2=1390.49 on 1267 dof (their reduction -- not
  comparable to ours in absolute terms). Also an **"early caustic"**
  alternative: q=0.0070, rho=0.00104, s=1.121, alpha=218.9 deg, tE=58.5,
  chi2 +7.4. Bennett et al. 2006 quote t*=0.059(7) d (= rho*tE, consistent).
- **Finite source recovers Bond's best fit.** Cassan fit from Bond's
  `(s, q)`, `rho` free: **chi2=1650.06**, s=1.1223, q=0.00390, rho=0.00097,
  alpha=224.4 deg, t0=2848.03, u0=0.123, tE=65.1 d, caustic entry/exit
  2834.26/2842.07 (entry near the user's by-eye ~2835 anchor). 116 below
  session 13's point-source best (1766), which, in hindsight, had matched
  Bond's *early-caustic* solution (q=0.0079, alpha~214.6, entry ~2832.6)
  rather than the best one. The high MOA points at the exit (HJD~2842) are now
  fitted by a finite-height spike rather than chased or missed. A few MOA
  points just after the entry (~2835) sit ~2-3 below the model.
- At Bond's published geometry, chi2 vs rho: 1952.90 (rho=1e-5, ~point
  source) -> **1776.24 (rho=0.00096)** -> 1889.83 (rho=3e-3) -- a clean,
  independent confirmation of both the finite-source code and Bond's rho.
- Fit from session 11's Student-t `(s, q)`=(1.641, 0.106): stuck,
  chi2=3304.85 (worse than its own point-source stage, 3267.99), tE=214 d
  (outside the MCMC prior's 10-150), rho unmoved from its start. A wrong
  basin, not a competitor (Delta chi2 ~1655).
- **Open**: in the winning fit, q finished at exactly its start (0.00390) and
  rho moved only 0.00096 -> 0.00097. Could be genuine (Bond's rho
  uncertainty, 0.00011, is about Nelder-Mead's initial 5% step) or an
  under-explored simplex -- only an MCMC posterior will tell.
- **Validation of `binary_magnification_fs`** (VBBL check saved to the repo,
  the rest were one-off scratchpad scripts):
  - rho->0 matches `binary_magnification` to ~1e-13.
  - Point lens, source centred (exact `A=sqrt(1+4/rho^2)`): error -0.302/sqrt(n_r)
    with rings -- matches the analytic midpoint-rule error on the 1/sqrt(u)
    singularity to three digits, i.e. the method's limit, not a bug.
  - vs MulensModel VBBL (source stepped along the caustic normal at the
    session-13 fit's exit point), saved as
    `scratch/cross_check_mulensmodel_fs.py`: ~1e-6 away from the caustic;
    straddling (|d|<1 rho) up to 2.6% at n=1000, **0.3-1.6% at n=4000**,
    0.5-0.9% at n=16000 -- ~1/sqrt(n), noisy. Hence the n=4000 default
    (below MOA's few-% errors). The gate edge (d=+5 rho) costs 0.1%.
  - Outside the caustic the finite-source effect vanishes by d~1 rho;
    inside a fold it decays slowly (~9% at 1 rho, ~1% at 2, ~0.3% at 4, ~0.05%
    at 8) -- hence `gate=5` (worst gate miss ~0.15%; `gate=3` measured
    0.57%).
- **Speed**: ungated disk average on 1800 points = 16 s at n=1250 / 97 s at
  n=4000; gated = 1.4 s; point source = 0.13 s; caustic computation inside
  the gate = 0.02 s. One finite-source chi2 on the real data ~0.8 s -> a
  32x3000 Cassan MCMC is hours even on 6 cores.
- Why source-plane brute force converges slowly near caustics: the
  integrand jumps and diverges along the caustic line. The standard
  alternatives work in the image plane or avoid the full integral
  (hexadecapole, Gould 2008 / Pejcha & Heyrovsky 2009; contour integration,
  Gould & Gaucherel 1997 / Dominik 1998 / Bozza 2010 VBBinaryLensing) --
  discussed, not implemented; the brute-force version stays as the
  reference either way.
- Minor: a VS Code Jupyter kernel desync (cell "running", kernel idle in its
  message loop) cost some time; fixed by restarting/reloading, not a code
  issue.
- Still stale/unverified from session 13: `caustic_curve()`'s docstring;
  `equidistant_caustic()` spacing never numerically checked (though the gate
  measured max caustic-point spacing 0.29 rho at `n_phi=1000`, which is
  indirect evidence).

### Next session
Confirmed with the user (quick question round, per session 11's note that
long `/grill-me` rounds add little once results are in hand):
1. **Finite-source Cassan MCMC** (first): `cassan_caustic.run_mcmc()` around
   the chi2=1650.06 fit, for real posteriors on `rho` and `q` -- settles
   whether Nelder-Mead actually explored them (see "Open" above). Hours at
   ~0.8 s/chi2; a cluster candidate. It will overwrite the point-source
   `O-03-BLG235_2l1s_cassan_mcmc_chain.npz` with 7 columns -- the notebook
   cells that unpack 6 values from it need a `, rho` added then.
2. **Error-bar rescaling** (roadmap item 2): additive jitter, Gaussian chi2
   refit, now that there's a trustworthy trajectory to rescale around.
- Not scheduled but discussed: seeding a finite-source fit from Bond's
  early-caustic solution to reproduce their Delta chi2 (+7.4); a
  hexadecapole / contour-integration speed-up; fixing `caustic_curve()`'s
  stale docstring.
