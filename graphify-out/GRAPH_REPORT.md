# Graph Report - microlensing  (2026-10-08)

## Corpus Check
- 11 files · ~51,618 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 493 nodes · 1188 edges · 33 communities (27 shown, 6 thin omitted)
- Extraction: 96% EXTRACTED · 4% INFERRED · 0% AMBIGUOUS · INFERRED: 51 edges (avg confidence: 0.87)
- Token cost: 90,573 input · 0 output

## Community Hubs (Navigation)
- MCMC and Cassan Helpers
- O-03-BLG235 Preprocessing
- Old PSPL Pipeline Scripts
- Event Data and Flux Model
- Plotting and VBBL FSPL
- Cassan Caustic Parametrisation
- Session 22 Findings
- PSPL Flux and FitResult
- Project Guide and BIC
- Binary Lens Trajectory
- 2L1S Chi2 and Diagnose
- lc_models Parallax Core
- MulensModel Cross-checks
- Joint PSPL Fit
- FSPL and Parallax Search
- Search Outputs and Caches
- Grid Miss and Residuals
- Hand-rolled Quintic Solver
- Config Loading
- Pool and Checkpointing
- Mass Estimate and Deps
- Roadmap and References
- MulensModel Fit Script
- O-05-BLG169 Verdict
- Error-bar Rescaling
- Header Coordinates
- Mode Deduplication
- Data Download
- Quintic Derivation
- SLURM Submit Script
- NamedTuple Type
- ndarray Type
- torch Dependency

## God Nodes (most connected - your core abstractions)
1. `CHANGELOG session log` - 44 edges
2. `CLAUDE.md project guide` - 35 edges
3. `Event` - 24 edges
4. `diagnose()` - 22 edges
5. `plot_fit()` - 21 edges
6. `TwoL1SParams` - 19 edges
7. `magnification()` - 19 edges
8. `binary_trajectory()` - 19 edges
9. `plot_fit()` - 18 edges
10. `trajectory()` - 17 edges

## Surprising Connections (you probably didn't know these)
- `TwoL1SParams canonical param order` --references--> `TwoL1SParams`  [EXTRACTED]
  CLAUDE.md → scratch/fit_2l1s.py
- `u0 sign degeneracy / required u0_guess` --references--> `run_ogle_only_diagnostic()`  [EXTRACTED]
  CLAUDE.md → mcmc_fit_binary.py
- `diagnose summary.txt readable output` --references--> `diagnose()`  [EXTRACTED]
  CHANGELOG.md → search.py
- `find_zoom_window significance handling` --references--> `find_zoom_window()`  [EXTRACTED]
  CHANGELOG.md → zoom_utils.py
- `Profiled per-instrument fs/fb (weighted lstsq)` --references--> `profile_flux()`  [EXTRACTED]
  CLAUDE.md → scratch/fit_2l1s.py

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **search.py config-driven stages (FSPL -> grid -> refine/MCMC -> BIC -> diagnose)** — search_fit_fspl, search_run_grid, search_refine, search_run_mcmc, search_bic, search_diagnose [EXTRACTED 1.00]
- **Session 22 robustness fixes against VBBL hangs** — changelog_vbbl_hang_pool_fix, changelog_rho_min_guard, changelog_mcmc_checkpointing, changelog_config_hash_cache, changelog_check_search_selfcheck [EXTRACTED 1.00]
- **Likelihood choices (Gaussian chi2, Student-t, Huber)** — claude_student_t_likelihood, claude_huber_likelihood, scratch_fit_2l1s_chi2 [INFERRED 0.85]

## Communities (33 total, 6 thin omitted)

### Community 0 - "MCMC and Cassan Helpers"
Cohesion: 0.05
Nodes (76): emcee, EnsembleSampler, Inverse of cassan_to_standard(): every (sigma, t) where the parallax-free…, standard_to_cassan(), flat_chain(), ndarray, Corner plot of the raw MCMC parameters -- shared by every fit script.…, Flat (samples, log_probs) after burn-in/thinning -- shared by every emcee… (+68 more)

### Community 1 - "O-03-BLG235 Preprocessing"
Cohesion: 0.14
Nodes (24): mag_to_flux(), Calibrated magnitude + error -> flux + error on the ZERO_POINT_MAG scale., matplotlib_pyplot, numpy, pathlib, load_raw(), moa_to_magnification(), ogle_to_magnification() (+16 more)

### Community 2 - "Old PSPL Pipeline Scripts"
Cohesion: 0.12
Nodes (27): argparse, corner, magnitude(), Returns (delta_s_n, delta_s_e): Earth's sky-projected position relative to the…, Observed magnitude, for overplotting against real photometry., sun_earth_projection(), plot_raw(), MCMC PSPL fits for OGLE-2003-BLG-235 / MOA-2003-BLG-53. Both fits below use the… (+19 more)

### Community 3 - "Event Data and Flux Model"
Cohesion: 0.14
Nodes (25): chi2(), error_scale(), Event, flux_residuals(), Instrument, profile_flux(), NamedTuple, ndarray (+17 more)

### Community 4 - "Plotting and VBBL FSPL"
Cohesion: 0.09
Nodes (23): datetime, hashlib, fspl_magnification(), lens_position(), Finite-source point-lens A(u) via VBBL's ESPLMag2. `rho` is the source radius…, Returns the lens masses and positions from separation s and mass ratio q, mpl_toolkits_axes_grid1_inset_locator, scipy_ndimage (+15 more)

### Community 5 - "Cassan Caustic Parametrisation"
Cohesion: 0.12
Nodes (23): Cassan 2008 curvilinear abscissa, Cassan (2008) caustic-crossing parametrisation, Near-miss cusp solutions problem, complex128, binary_magnification_fs(), _disk_average(), _sample_coords(), cassan_caustic() (+15 more)

### Community 6 - "Session 22 Findings"
Cohesion: 0.14
Nodes (22): CHANGELOG session log, Autocorrelation (tau) convergence check, Blending degeneracy (model-dependent magnification), scratch/check_search.py self-check, Prior-based lens mass estimate, Fair fit_lc.png in OGLE I magnitude, find_zoom_window significance handling, Han & Gould 2003 Galactic prior (+14 more)

### Community 7 - "PSPL Flux and FitResult"
Cohesion: 0.13
Nodes (21): FitResult NamedTuple (samples + labels), flux(), plain_flux(), plain_magnitude(), Observed flux: magnified source plus a constant blend flux., flux() with piE_N=piE_E=0 -- the exact pre-parallax model, used only to…, magnitude() with piE_N=piE_E=0 -- see plain_flux()., fit_parallax_pspl_mcmc() (+13 more)

### Community 8 - "Project Guide and BIC"
Cohesion: 0.12
Nodes (21): BIC caveats (look-elsewhere, correlated errors, circular K), CLAUDE.md project guide, PSPL vs 2L1S BIC comparison, Blind search, no literature seeds, Bond et al. 2004 published O-03-BLG235 solution, O-03-BLG235 calibrate-once-then-fit pipeline, Finite-source binary magnification (Fibonacci disk average), Instrument kind mag vs dia flux model (+13 more)

### Community 9 - "Binary Lens Trajectory"
Cohesion: 0.18
Nodes (18): concurrent_futures, binary_magnification(), binary_trajectory(), Total binary-lens magnification at source position(s) `zeta`. Sums 1/|J| over…, Source-lens separation for binary lenses on the lens plane, given in complex…, plot_panel(), chi2(), _fit_one_seed() (+10 more)

### Community 10 - "2L1S Chi2 and Diagnose"
Cohesion: 0.17
Nodes (17): Binary, binary_A(), binary_steps(), chi2_binary(), diagnose(), log_prob(), plot_parallax_map(), NamedTuple (+9 more)

### Community 11 - "lc_models Parallax Core"
Cohesion: 0.15
Nodes (14): astropy_coordinates, astropy_time, astropy_units, functools, _project(), Functional forms used to model a microlensing light curve. Point-source point-…, Project a Cartesian position/velocity vector onto the sky's (N, E) tangent-…, numpy_typing (+6 more)

### Community 12 - "MulensModel Cross-checks"
Cohesion: 0.14
Nodes (12): itertools, mulensmodel, scipy_optimize, One-time cross-check: fit MulensModel's own binary-lens model to the same…, One-time offline cross-check of lc_models.binary_magnification_fs against…, mm_magnification(), One-time offline cross-check of lc_models.binary_magnification against…, MulensModel's magnification at the same complex zeta values. MulensModel's docs… (+4 more)

### Community 13 - "Joint PSPL Fit"
Cohesion: 0.23
Nodes (17): magnification(), Source-lens separation u(t), in units of the Einstein radius., Paczynski point-lens magnification A(u)., trajectory(), run_joint_fit(), chi2(), log_likelihood(), log_prior() (+9 more)

### Community 14 - "FSPL and Parallax Search"
Cohesion: 0.24
Nodes (16): nelder_mead(), Nelder-Mead from an explicit initial simplex: scipy's default steps 5% of each…, anomaly_times(), chi2_fspl(), fit_fspl(), FSPL, fspl_A(), fspl_cell() (+8 more)

### Community 15 - "Search Outputs and Caches"
Cohesion: 0.14
Nodes (14): Concurrent MCMC modes on shared Pool, Path, acceptance(), cached(), out_dir(), plot_map(), plot_raw(), plot_trace() (+6 more)

### Community 16 - "Grid Miss and Residuals"
Cohesion: 0.18
Nodes (12): Bennett 2015 re-reduction of O-05-BLG169, Correlated MDM residuals in O-05-BLG169 2L1S, Grid inner search N_SIGMA/N_POLISH/STD_MAXFEV 16/5/600, Why the coarse grid missed Bond's basin, O-03-BLG235 search failing done-checks, Drop-one instrument sensitivity runs, Dataset short-name mapping, -fineq variant (log q step 0.1) (+4 more)

### Community 17 - "Hand-rolled Quintic Solver"
Cohesion: 0.18
Nodes (12): Symbolically derived binary-lens quintic, Torch batched eigvals root solver, binary_images(), _binary_images_batch(), _companion_eigvals(), _quintic_coefficients(), Coefficients (highest degree first) of the degree-5 polynomial in z equivalent…, (torch, device), imported on first use: only the hand-rolled solver needs it,… (+4 more)

### Community 18 - "Config Loading"
Cohesion: 0.18
Nodes (10): Geocentric JD light-travel correction, concurrent_futures_process, load_event(), load_instrument(), One [[instruments]] entry -> Instrument: loadtxt, time shift by time_fmt, K…, Parse the TOML at `path` (tomllib, opened "rb") and load every instrument., os, Self-check for search.py's session-22 robustness fixes (~1 min; srun, not the… (+2 more)

### Community 19 - "Pool and Checkpointing"
Cohesion: 0.22
Nodes (10): Cache keyed on config hash + inner-search settings, MCMC chain checkpointing every 200 steps, Grid checkpoint/resume (grid_partial.npz), axis(), _init(), pool(), [lo, hi, step] -> lo..hi inclusive., Spawned worker pool. ProcessPoolExecutor, not multiprocessing.Pool: when a… (+2 more)

### Community 20 - "Mass Estimate and Deps"
Cohesion: 0.22
Nodes (8): Mass-distance degeneracy in lens mass estimate, Prior-based lens mass estimate (D_L uniform, mu_rel log-normal), Flat-bound MCMC priors, astropy, corner, emcee, torch, VBBinaryLensing 3.7.0

### Community 21 - "Roadmap and References"
Cohesion: 0.25
Nodes (8): Charbonneau 1995 genetic algorithm, Kains 2009 reference application, Planned 2L1S roadmap (rescale -> finite source -> Cassan -> grid -> GA), SIGNALMEN anomaly-vs-outlier test (Dominik et al. 2007), Config-driven pipeline (search.py + event.py + TOML), Huber loss likelihood (DELTA=1.345), Student-t robust likelihood (scale, dof free), Event

### Community 22 - "MulensModel Fit Script"
Cohesion: 0.33
Nodes (7): chi2(), mm_magnification(), plot_fit(), Chi2 of MulensModel's binary magnification against the combined OGLE+MOA data., Serial multi-start search, then a tight-tolerance refit of the best seed., Overlay the MulensModel fit on the OGLE+MOA data, full baseline + zoomed on the…, run_fit()

### Community 23 - "O-05-BLG169 Verdict"
Cohesion: 0.33
Nodes (6): Offset-scaled archive magnitudes (FTN, MDM, Auckland), Gould et al. 2006 (O-05-BLG169 discovery), O-05-BLG169 2L1S verdict (Delta BIC ~132, q~3.5e-5), Parallax-only null search (parallax_search), Annual parallax (Gould 2004 geocentric), O-05-BLG169 (OGLE-2005-BLG-169)

### Community 24 - "Error-bar Rescaling"
Cohesion: 0.40
Nodes (5): MOA-2019-BLG-008 error bars underestimated 16-49x, Per-instrument error-bar rescaling k, M-19-BLG008 (MOA-2019-BLG-008, KMT I subset), error_rescaling(), Rescales the error bars per instrument by settings chi2/dof ~ 1. raw_resid, err…

### Community 25 - "Header Coordinates"
Cohesion: 0.50
Nodes (4): load_coords(), _parse_header_coord(), Extract (RA, Dec) sexagesimal strings from a raw table's own \\RA/\\DEC header…, Target coordinates, parsed from OGLE_PATH/MOA_PATH's own headers rather than…

### Community 26 - "Mode Deduplication"
Cohesion: 0.67
Nodes (3): distinct_modes tolerance-based dedup, distinct_modes(), Refined (x, chi2) pairs, sorted by chi2 -> the Binary modes within MODE_DCHI2…

## Knowledge Gaps
- **10 isolated node(s):** `submit_all.sh script`, `astropy`, `corner`, `emcee`, `torch` (+5 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 195 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **6 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `CHANGELOG session log` connect `Session 22 Findings` to `Cassan Caustic Parametrisation`, `Project Guide and BIC`, `Search Outputs and Caches`, `Grid Miss and Residuals`, `Hand-rolled Quintic Solver`, `Config Loading`, `Pool and Checkpointing`, `Mass Estimate and Deps`, `Roadmap and References`, `O-05-BLG169 Verdict`, `Error-bar Rescaling`, `Mode Deduplication`?**
  _High betweenness centrality (0.106) - this node is a cross-community bridge._
- **Why does `CLAUDE.md project guide` connect `Project Guide and BIC` to `Cassan Caustic Parametrisation`, `Session 22 Findings`, `PSPL Flux and FitResult`, `Grid Miss and Residuals`, `Pool and Checkpointing`, `Mass Estimate and Deps`, `Roadmap and References`, `O-05-BLG169 Verdict`, `Error-bar Rescaling`?**
  _High betweenness centrality (0.097) - this node is a cross-community bridge._
- **Why does `find_zoom_window()` connect `O-03-BLG235 Preprocessing` to `MCMC and Cassan Helpers`, `Old PSPL Pipeline Scripts`, `Event Data and Flux Model`, `Plotting and VBBL FSPL`, `Session 22 Findings`, `Binary Lens Trajectory`, `MulensModel Cross-checks`, `Search Outputs and Caches`?**
  _High betweenness centrality (0.041) - this node is a cross-community bridge._
- **Are the 15 inferred relationships involving `Event` (e.g. with `anomaly_times()` and `bic()`) actually correct?**
  _`Event` has 15 INFERRED edges - model-reasoned connections that need verification._
- **What connects `submit_all.sh script`, `astropy`, `corner` to the rest of the system?**
  _10 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `MCMC and Cassan Helpers` be split into smaller, more focused modules?**
  _Cohesion score 0.05194805194805195 - nodes in this community are weakly interconnected._
- **Should `O-03-BLG235 Preprocessing` be split into smaller, more focused modules?**
  _Cohesion score 0.14482758620689656 - nodes in this community are weakly interconnected._