# Graph Report - microlensing  (2026-10-09)

## Corpus Check
- 28 files · ~59,417 words
- Verdict: corpus is large enough that graph structure adds value.
- Unclassified: 30 file(s) not represented in the graph (top: .npz 10, .toml 9, .sbatch 8)

## Summary
- 540 nodes · 1348 edges · 38 communities (34 shown, 4 thin omitted)
- Extraction: 95% EXTRACTED · 5% INFERRED · 0% AMBIGUOUS · INFERRED: 66 edges (avg confidence: 0.87)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `8464027a`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- mcmc_fit_2l1s.py
- fit_2l1s.py
- mcmc_fit.py
- event.py
- search.py
- plot_fit
- CHANGELOG session log
- load_instrument
- CLAUDE.md project guide
- binary_trajectory
- diagnose
- lc_models.py
- numpy
- Event
- ld_test.py
- cassan_caustic.py
- Dataset short-name mapping
- _binary_images_batch
- plot_fit
- run_grid
- README.md
- Planned 2L1S roadmap (rescale -> finite source -> Cassan -> grid -> GA)
- plot_fit
- caustic_curve
- Per-instrument error-bar rescaling k
- log_prior
- distinct_modes
- TwoL1SParams
- derive_binary_quintic.py
- submit_all.sh
- preprocess_binary_data.py
- probe_fs_grid.py
- torch
- chi2
- mdm_systematics.py
- pool
- flat_chain
- download_data.py

## God Nodes (most connected - your core abstractions)
1. `CHANGELOG session log` - 44 edges
2. `CLAUDE.md project guide` - 35 edges
3. `Event` - 29 edges
4. `diagnose()` - 26 edges
5. `plot_fit()` - 21 edges
6. `binary_trajectory()` - 20 edges
7. `magnification()` - 19 edges
8. `TwoL1SParams` - 19 edges
9. `caustic_curve()` - 18 edges
10. `sun_earth_projection()` - 18 edges

## Surprising Connections (you probably didn't know these)
- `u0 sign degeneracy / required u0_guess` --references--> `run_ogle_only_diagnostic()`  [EXTRACTED]
  CLAUDE.md → mcmc_fit_binary.py
- `TwoL1SParams canonical param order` --references--> `TwoL1SParams`  [EXTRACTED]
  CLAUDE.md → scratch/fit_2l1s.py
- `diagnose summary.txt readable output` --references--> `diagnose()`  [EXTRACTED]
  CHANGELOG.md → search.py
- `find_zoom_window significance handling` --references--> `find_zoom_window()`  [EXTRACTED]
  CHANGELOG.md → zoom_utils.py
- `Geocentric JD light-travel correction` --references--> `load_instrument()`  [EXTRACTED]
  CHANGELOG.md → event.py

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **search.py config-driven stages (FSPL -> grid -> refine/MCMC -> BIC -> diagnose)** — search_fit_fspl, search_run_grid, search_refine, search_run_mcmc, search_bic, search_diagnose [EXTRACTED 1.00]
- **Session 22 robustness fixes against VBBL hangs** — changelog_vbbl_hang_pool_fix, changelog_rho_min_guard, changelog_mcmc_checkpointing, changelog_config_hash_cache, changelog_check_search_selfcheck [EXTRACTED 1.00]
- **Likelihood choices (Gaussian chi2, Student-t, Huber)** — claude_student_t_likelihood, claude_huber_likelihood, scratch_fit_2l1s_chi2 [INFERRED 0.85]

## Communities (38 total, 4 thin omitted)

### Community 0 - "mcmc_fit_2l1s.py"
Cohesion: 0.18
Nodes (16): Corner plot of the raw MCMC parameters -- shared by every fit script.…, save_corner(), log_probability(), MCMC (emcee) exploration of the full 2L1S parameter space for…, Draw n walker start columns for the physical params -- shared by all three…, Draw n walker start positions directly from the prior -- broad coverage in…, Draw n walker start positions directly from the prior -- broad coverage in…, Draw n walker start positions directly from the prior -- broad coverage in… (+8 more)

### Community 1 - "fit_2l1s.py"
Cohesion: 0.16
Nodes (15): functools, mag_to_flux(), Calibrated magnitude + error -> flux + error on the ZERO_POINT_MAG scale., mpl_toolkits_axes_grid1_inset_locator, moa_to_magnification(), ogle_to_magnification(), Invert OGLE's flux model to recover A(t): A = (F - fb_ogle) / fs_ogle., Invert MOA's flux model to recover A(t): A = 1 + F / fs_moa. (+7 more)

### Community 2 - "mcmc_fit.py"
Cohesion: 0.05
Nodes (76): argparse, Annual parallax (Gould 2004 geocentric), FitResult NamedTuple (samples + labels), corner, emcee, flux(), magnification(), magnitude() (+68 more)

### Community 3 - "event.py"
Cohesion: 0.12
Nodes (19): concurrent_futures_process, chi2(), error_scale(), flux_residuals(), ld_penalty(), load_event(), ndarray, One event's config + data, and the model-agnostic half of every chi2. Library… (+11 more)

### Community 4 - "search.py"
Cohesion: 0.13
Nodes (25): datetime, hashlib, scipy_ndimage, as_binary(), chi2_binary(), from_mcmc(), in_range(), ld_ok() (+17 more)

### Community 5 - "plot_fit"
Cohesion: 0.12
Nodes (21): Instrument, profile_flux(), NamedTuple, Weighted least-squares (fs, fb) ("mag": flux = fs*A + fb) or (fs,) ("dia": flux…, Invert an instrument's flux onto magnification via its fs/fb profiled at model…, One telescope/band: its [[instruments]] config entry plus its loaded data., to_magnification(), Runs-test z (negative: fewer sign changes than chance), lag-1 autocorrelation,… (+13 more)

### Community 6 - "CHANGELOG session log"
Cohesion: 0.12
Nodes (26): CHANGELOG session log, Autocorrelation (tau) convergence check, Blending degeneracy (model-dependent magnification), scratch/check_search.py self-check, Prior-based lens mass estimate, Fair fit_lc.png in OGLE I magnitude, find_zoom_window significance handling, Gould et al. 2006 (O-05-BLG169 discovery) (+18 more)

### Community 7 - "load_instrument"
Cohesion: 0.50
Nodes (4): Geocentric JD light-travel correction, load_instrument(), One [[instruments]] entry -> Instrument: loadtxt, time shift by time_fmt, K…, SkyCoord

### Community 8 - "CLAUDE.md project guide"
Cohesion: 0.13
Nodes (18): BIC caveats (look-elsewhere, correlated errors, circular K), CLAUDE.md project guide, PSPL vs 2L1S BIC comparison, Blind search, no literature seeds, O-03-BLG235 calibrate-once-then-fit pipeline, Instrument kind mag vs dia flux model, No compute on shared login node (srun/sbatch), MulensModel cross-check oracle (+10 more)

### Community 9 - "binary_trajectory"
Cohesion: 0.17
Nodes (19): itertools, binary_magnification(), binary_trajectory(), Total binary-lens magnification at source position(s) `zeta`. Sums 1/|J| over…, Source-lens separation for binary lenses on the lens plane, given in complex…, plot_panel(), Compare every 2L1S candidate solution found this session (Nelder-Mead multi-…, chi2() (+11 more)

### Community 10 - "diagnose"
Cohesion: 0.09
Nodes (24): Concurrent MCMC modes on shared Pool, Path, acceptance(), bic(), cached(), caustic_crossing(), diagnose(), fmt_params() (+16 more)

### Community 11 - "lc_models.py"
Cohesion: 0.16
Nodes (12): astropy_time, astropy_units, caustic_point(), _project(), Functional forms used to model a microlensing light curve. Point-source point-…, Point at abscissa sigma (wraps with period 1) on a cassan_caustic() output., Project a Cartesian position/velocity vector onto the sky's (N, E) tangent-…, numpy_typing (+4 more)

### Community 12 - "numpy"
Cohesion: 0.22
Nodes (10): mulensmodel, numpy, pathlib, One-time cross-check: fit MulensModel's own binary-lens model to the same…, One-time offline cross-check of lc_models.binary_magnification_fs against…, mm_magnification(), One-time offline cross-check of lc_models.binary_magnification against…, MulensModel's magnification at the same complex zeta values. MulensModel's docs… (+2 more)

### Community 13 - "Event"
Cohesion: 0.18
Nodes (19): Event, Event with the free bands' LD set to `ld` (one per event.ld_prior entry); ()…, with_ld(), fspl_magnification(), Finite-source point-lens A(u) via VBBL's ESPLMag2. `rho` is the source radius…, anomaly_times(), as_fspl(), chi2_fspl() (+11 more)

### Community 14 - "ld_test.py"
Cohesion: 0.17
Nodes (19): polish(), nelder_mead() restarted from its own result with a fresh simplex (`steps(x)`)…, fit(), Free limb darkening test (session 23, one-off): does fitting the linear…, (FSPL, 2L1S overall best, K at FSPL) from a search.py summary.txt (alpha…, Pool task: (model, free LD?, start, event) -> (x, chi2); LD appended to x when…, read_summary(), with_ld() (+11 more)

### Community 15 - "cassan_caustic.py"
Cohesion: 0.21
Nodes (16): Inverse of cassan_to_standard(): every (sigma, t) where the parallax-free…, standard_to_cassan(), CassanParams, chi2_cassan(), fit(), log_probability(), NamedTuple, 2L1S fit of OGLE-2003-BLG-235 in Cassan (2008)'s caustic-crossing… (+8 more)

### Community 16 - "Dataset short-name mapping"
Cohesion: 0.15
Nodes (16): Bennett 2015 re-reduction of O-05-BLG169, Correlated MDM residuals in O-05-BLG169 2L1S, Offset-scaled archive magnitudes (FTN, MDM, Auckland), Grid inner search N_SIGMA/N_POLISH/STD_MAXFEV 16/5/600, Why the coarse grid missed Bond's basin, O-03-BLG235 search failing done-checks, Bond et al. 2004 published O-03-BLG235 solution, Drop-one instrument sensitivity runs (+8 more)

### Community 17 - "_binary_images_batch"
Cohesion: 0.15
Nodes (14): Symbolically derived binary-lens quintic, Torch batched eigvals root solver, binary_images(), _binary_images_batch(), _companion_eigvals(), lens_position(), _quintic_coefficients(), Coefficients (highest degree first) of the degree-5 polynomial in z equivalent… (+6 more)

### Community 18 - "plot_fit"
Cohesion: 0.19
Nodes (13): _data_delta_s(), delta_s(), plot_fit(), model_fn(), plot_panel(), profile_flux(), _profile_fs_moa(), Overlay the 2L1S model on the data: event season (HJD 2700-3000), then a zoomed… (+5 more)

### Community 19 - "run_grid"
Cohesion: 0.13
Nodes (18): Cache keyed on config hash + inner-search settings, MCMC chain checkpointing every 200 steps, Grid checkpoint/resume (grid_partial.npz), cassan_to_standard(), Cassan (2008) caustic-crossing parameters -> binary_trajectory()'s (t0, u0, tE,…, axis(), finite_grid(), grid_cell() (+10 more)

### Community 20 - "README.md"
Cohesion: 0.22
Nodes (8): Mass-distance degeneracy in lens mass estimate, Prior-based lens mass estimate (D_L uniform, mu_rel log-normal), Flat-bound MCMC priors, astropy, corner, emcee, torch, VBBinaryLensing 3.7.0

### Community 21 - "Planned 2L1S roadmap (rescale -> finite source -> Cassan -> grid -> GA)"
Cohesion: 0.20
Nodes (10): Cassan 2008 curvilinear abscissa, Charbonneau 1995 genetic algorithm, Kains 2009 reference application, Planned 2L1S roadmap (rescale -> finite source -> Cassan -> grid -> GA), SIGNALMEN anomaly-vs-outlier test (Dominik et al. 2007), Cassan (2008) caustic-crossing parametrisation, Config-driven pipeline (search.py + event.py + TOML), Huber loss likelihood (DELTA=1.345) (+2 more)

### Community 22 - "plot_fit"
Cohesion: 0.33
Nodes (7): chi2(), mm_magnification(), plot_fit(), Chi2 of MulensModel's binary magnification against the combined OGLE+MOA data., Serial multi-start search, then a tight-tolerance refit of the best seed., Overlay the MulensModel fit on the OGLE+MOA data, full baseline + zoomed on the…, run_fit()

### Community 23 - "caustic_curve"
Cohesion: 0.23
Nodes (12): complex128, binary_magnification_fs(), _disk_average(), _sample_coords(), cassan_caustic(), caustic_curve(), equidistant_caustic(), NDArray (+4 more)

### Community 24 - "Per-instrument error-bar rescaling k"
Cohesion: 0.40
Nodes (5): MOA-2019-BLG-008 error bars underestimated 16-49x, Per-instrument error-bar rescaling k, M-19-BLG008 (MOA-2019-BLG-008, KMT I subset), error_rescaling(), Rescales the error bars per instrument by settings chi2/dof ~ 1. raw_resid, err…

### Community 25 - "log_prior"
Cohesion: 0.20
Nodes (11): log_prior(), log_prior_chi2(), log_prior_huber(), log_probability_chi2(), _log_range_ok(), _physical_log_prior(), Same bounds as log_prior(), minus scale/dof -- this run has no Student-t params., Plain Gaussian likelihood (-0.5*chi2) in place of log_probability()'s Student-t… (+3 more)

### Community 26 - "distinct_modes"
Cohesion: 0.67
Nodes (3): distinct_modes tolerance-based dedup, distinct_modes(), Refined (x, chi2, ...) sorted by chi2, and whether each crosses a caustic ->…

### Community 27 - "TwoL1SParams"
Cohesion: 0.29
Nodes (10): huber(), NamedTuple, ndarray, Per-point standardized residuals of the 2L1S model against raw flux. Always…, The 2L1S track's 9 physical parameters, in one canonical order -- every call…, residuals(), TwoL1SParams, log_probability_huber() (+2 more)

### Community 30 - "preprocess_binary_data.py"
Cohesion: 0.25
Nodes (8): astropy_coordinates, load_coords(), load_raw(), _parse_header_coord(), Reduce the raw OGLE + MOA files for OGLE-2003-BLG-235 / MOA-2003-BLG-53 into…, Extract (RA, Dec) sexagesimal strings from a raw table's own \\RA/\\DEC header…, Target coordinates, parsed from OGLE_PATH/MOA_PATH's own headers rather than…, Raw OGLE (mag) + MOA (differential flux) tables, time-shifted to HJD-2450000.

### Community 31 - "probe_fs_grid.py"
Cohesion: 0.25
Nodes (7): concurrent_futures, Event with each instrument's flux_err multiplied by its k (on top of the…, rescale(), multiprocessing, cell(), Finite-source grid probe (session 23, one-off): what would search.grid_cell()…, time

### Community 33 - "chi2"
Cohesion: 0.33
Nodes (7): chi2(), _fit_one_seed(), Chi2 of the 2L1S model against raw flux, with each instrument's flux…, Nelder-Mead on one seed at loose tolerance; module-level so it can be pickled., Parallel multi-start search, then a tight-tolerance refit of the best seed., run_fit(), run_fit_moa_only()

### Community 34 - "mdm_systematics.py"
Cohesion: 0.33
Nodes (5): astropy_timeseries, matplotlib_pyplot, detrended(), MDM systematics check (session 24, one-off): is O-05-BLG169's MDM night-1 wave…, Standardized residuals left after a weighted linear fit of r against x.

### Community 35 - "pool"
Cohesion: 0.40
Nodes (4): init(), _init(), pool(), Spawned worker pool. ProcessPoolExecutor, not multiprocessing.Pool: when a…

### Community 36 - "flat_chain"
Cohesion: 0.50
Nodes (4): EnsembleSampler, flat_chain(), ndarray, Flat (samples, log_probs) after burn-in/thinning -- shared by every emcee…

## Knowledge Gaps
- **10 isolated node(s):** `submit_all.sh script`, `astropy`, `corner`, `emcee`, `torch` (+5 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 210 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **4 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `CHANGELOG session log` connect `CHANGELOG session log` to `load_instrument`, `CLAUDE.md project guide`, `diagnose`, `Dataset short-name mapping`, `_binary_images_batch`, `run_grid`, `README.md`, `Planned 2L1S roadmap (rescale -> finite source -> Cassan -> grid -> GA)`, `Per-instrument error-bar rescaling k`, `distinct_modes`?**
  _High betweenness centrality (0.099) - this node is a cross-community bridge._
- **Why does `CLAUDE.md project guide` connect `CLAUDE.md project guide` to `mcmc_fit.py`, `CHANGELOG session log`, `Dataset short-name mapping`, `run_grid`, `README.md`, `Planned 2L1S roadmap (rescale -> finite source -> Cassan -> grid -> GA)`, `Per-instrument error-bar rescaling k`?**
  _High betweenness centrality (0.087) - this node is a cross-community bridge._
- **Why does `find_zoom_window()` connect `mcmc_fit.py` to `fit_2l1s.py`, `search.py`, `plot_fit`, `CHANGELOG session log`, `binary_trajectory`, `diagnose`, `plot_fit`?**
  _High betweenness centrality (0.036) - this node is a cross-community bridge._
- **Are the 18 inferred relationships involving `Event` (e.g. with `anomaly_times()` and `bic()`) actually correct?**
  _`Event` has 18 INFERRED edges - model-reasoned connections that need verification._
- **What connects `submit_all.sh script`, `astropy`, `corner` to the rest of the system?**
  _10 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `mcmc_fit.py` be split into smaller, more focused modules?**
  _Cohesion score 0.05063291139240506 - nodes in this community are weakly interconnected._
- **Should `event.py` be split into smaller, more focused modules?**
  _Cohesion score 0.12380952380952381 - nodes in this community are weakly interconnected._