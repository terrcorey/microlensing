# Graph Report - microlensing  (2026-10-09)

## Corpus Check
- 28 files · ~57,967 words
- Verdict: corpus is large enough that graph structure adds value.
- Unclassified: 26 file(s) not represented in the graph (top: .toml 8, .sbatch 8, .npz 7)

## Summary
- 528 nodes · 1314 edges · 28 communities (25 shown, 3 thin omitted)
- Extraction: 95% EXTRACTED · 5% INFERRED · 0% AMBIGUOUS · INFERRED: 63 edges (avg confidence: 0.87)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `8464027a`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- mcmc_fit_2l1s.py
- fit_2l1s.py
- mcmc_fit.py
- Event
- search.py
- plot_fit
- CHANGELOG session log
- load_instrument
- CLAUDE.md project guide
- binary_trajectory
- diagnose
- mdm_systematics.py
- pathlib
- plot_caustic_inset
- ld_test.py
- _project
- Dataset short-name mapping
- lc_models.py
- numpy
- run_grid
- README.md
- Planned 2L1S roadmap (rescale -> finite source -> Cassan -> grid -> GA)
- cross_check_mulensmodel_fit.py
- Per-instrument error-bar rescaling k
- distinct_modes
- derive_binary_quintic.py
- submit_all.sh
- torch

## God Nodes (most connected - your core abstractions)
1. `CHANGELOG session log` - 44 edges
2. `CLAUDE.md project guide` - 35 edges
3. `Event` - 25 edges
4. `diagnose()` - 25 edges
5. `plot_fit()` - 21 edges
6. `binary_trajectory()` - 20 edges
7. `magnification()` - 19 edges
8. `TwoL1SParams` - 19 edges
9. `FSPL` - 19 edges
10. `Binary` - 19 edges

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

## Communities (28 total, 3 thin omitted)

### Community 0 - "mcmc_fit_2l1s.py"
Cohesion: 0.06
Nodes (66): EnsembleSampler, cassan_to_standard(), Cassan (2008) caustic-crossing parameters -> binary_trajectory()'s (t0, u0, tE,…, flat_chain(), ndarray, Corner plot of the raw MCMC parameters -- shared by every fit script.…, Flat (samples, log_probs) after burn-in/thinning -- shared by every emcee…, save_corner() (+58 more)

### Community 1 - "fit_2l1s.py"
Cohesion: 0.08
Nodes (43): mag_to_flux(), magnification(), Calibrated magnitude + error -> flux + error on the ZERO_POINT_MAG scale., Source-lens separation u(t), in units of the Einstein radius., Paczynski point-lens magnification A(u)., trajectory(), chi2(), log_likelihood() (+35 more)

### Community 2 - "mcmc_fit.py"
Cohesion: 0.06
Nodes (63): argparse, Annual parallax (Gould 2004 geocentric), FitResult NamedTuple (samples + labels), corner, emcee, flux(), magnitude(), plain_flux() (+55 more)

### Community 3 - "Event"
Cohesion: 0.19
Nodes (16): chi2(), error_scale(), Event, flux_residuals(), load_event(), ndarray, One event's config + data, and the model-agnostic half of every chi2. Library…, Event with each instrument's parallax offsets precomputed at t0_par -- one… (+8 more)

### Community 4 - "search.py"
Cohesion: 0.08
Nodes (25): concurrent_futures, datetime, hashlib, fspl_magnification(), Finite-source point-lens A(u) via VBBL's ESPLMag2. `rho` is the source radius…, multiprocessing, scipy_ndimage, init() (+17 more)

### Community 5 - "plot_fit"
Cohesion: 0.12
Nodes (21): Instrument, profile_flux(), NamedTuple, Invert an instrument's flux onto magnification via its fs/fb profiled at model…, One telescope/band: its [[instruments]] config entry plus its loaded data., Weighted least-squares (fs, fb) ("mag": flux = fs*A + fb) or (fs,) ("dia": flux…, to_magnification(), binary_magnification_vbbl() (+13 more)

### Community 6 - "CHANGELOG session log"
Cohesion: 0.14
Nodes (23): CHANGELOG session log, Autocorrelation (tau) convergence check, Blending degeneracy (model-dependent magnification), scratch/check_search.py self-check, Prior-based lens mass estimate, Fair fit_lc.png in OGLE I magnitude, find_zoom_window significance handling, Gould et al. 2006 (O-05-BLG169 discovery) (+15 more)

### Community 7 - "load_instrument"
Cohesion: 0.50
Nodes (4): Geocentric JD light-travel correction, load_instrument(), One [[instruments]] entry -> Instrument: loadtxt, time shift by time_fmt, K…, SkyCoord

### Community 8 - "CLAUDE.md project guide"
Cohesion: 0.13
Nodes (19): BIC caveats (look-elsewhere, correlated errors, circular K), CLAUDE.md project guide, PSPL vs 2L1S BIC comparison, Blind search, no literature seeds, O-03-BLG235 calibrate-once-then-fit pipeline, Instrument kind mag vs dia flux model, No compute on shared login node (srun/sbatch), MulensModel cross-check oracle (+11 more)

### Community 9 - "binary_trajectory"
Cohesion: 0.07
Nodes (43): complex128, binary_magnification(), binary_magnification_fs(), _disk_average(), _sample_coords(), binary_trajectory(), cassan_caustic(), caustic_curve() (+35 more)

### Community 10 - "diagnose"
Cohesion: 0.10
Nodes (30): Concurrent MCMC modes on shared Pool, Path, acceptance(), Binary, binary_steps(), cached(), chi2_binary(), diagnose() (+22 more)

### Community 11 - "mdm_systematics.py"
Cohesion: 0.22
Nodes (8): astropy_coordinates, astropy_time, astropy_timeseries, astropy_units, One-time offline cross-check of lc_models' annual-parallax PSPL trajectory…, detrended(), MDM systematics check (session 24, one-off): is O-05-BLG169's MDM night-1 wave…, Standardized residuals left after a weighted linear fit of r against x.

### Community 12 - "pathlib"
Cohesion: 0.18
Nodes (9): Download raw photometry for each event we're working with.…, mulensmodel, pathlib, One-time offline cross-check of lc_models.binary_magnification_fs against…, mm_magnification(), One-time offline cross-check of lc_models.binary_magnification against…, MulensModel's magnification at the same complex zeta values. MulensModel's docs…, report() (+1 more)

### Community 13 - "plot_caustic_inset"
Cohesion: 0.50
Nodes (4): lens_position(), Returns the lens masses and positions from separation s and mass ratio q, plot_caustic_inset(), Square inset (physical inches + adjustable="datalim": a fraction-based or…

### Community 14 - "ld_test.py"
Cohesion: 0.16
Nodes (23): fit(), Free limb darkening test (session 23, one-off): does fitting the linear…, (FSPL, 2L1S overall best, K at FSPL) from a search.py summary.txt (alpha…, Pool task: (model, free LD?, start, event) -> (x, chi2); LD appended to x when…, read_summary(), with_ld(), anomaly_times(), chi2_fspl() (+15 more)

### Community 15 - "_project"
Cohesion: 0.50
Nodes (4): _project(), Project a Cartesian position/velocity vector onto the sky's (N, E) tangent-…, bare_earth_projection(), Same linear-detrending-at-t0_par logic as lc_models.sun_earth_projection(), but…

### Community 16 - "Dataset short-name mapping"
Cohesion: 0.15
Nodes (16): Bennett 2015 re-reduction of O-05-BLG169, Correlated MDM residuals in O-05-BLG169 2L1S, Offset-scaled archive magnitudes (FTN, MDM, Auckland), Grid inner search N_SIGMA/N_POLISH/STD_MAXFEV 16/5/600, Why the coarse grid missed Bond's basin, O-03-BLG235 search failing done-checks, Bond et al. 2004 published O-03-BLG235 solution, Drop-one instrument sensitivity runs (+8 more)

### Community 17 - "lc_models.py"
Cohesion: 0.12
Nodes (20): Symbolically derived binary-lens quintic, Torch batched eigvals root solver, functools, binary_images(), _binary_images_batch(), caustic_point(), _companion_eigvals(), _quintic_coefficients() (+12 more)

### Community 18 - "numpy"
Cohesion: 0.20
Nodes (9): concurrent_futures_process, itertools, numpy, os, scipy_optimize, Self-check for search.py's session-22 robustness fixes (~1 min; srun, not the…, Multi-start Nelder-Mead 2L1S fit for MOA-2019-BLG-008 (Bachelet et al. 2022), a…, sys (+1 more)

### Community 19 - "run_grid"
Cohesion: 0.20
Nodes (11): Cache keyed on config hash + inner-search settings, MCMC chain checkpointing every 200 steps, Grid checkpoint/resume (grid_partial.npz), axis(), finite_grid(), grid_timing(), [lo, hi, step] -> lo..hi inclusive., Does the source size shape the peak (parallax-free FSPL rho >= |u0| / 20)? Then… (+3 more)

### Community 20 - "README.md"
Cohesion: 0.22
Nodes (8): Mass-distance degeneracy in lens mass estimate, Prior-based lens mass estimate (D_L uniform, mu_rel log-normal), Flat-bound MCMC priors, astropy, corner, emcee, torch, VBBinaryLensing 3.7.0

### Community 21 - "Planned 2L1S roadmap (rescale -> finite source -> Cassan -> grid -> GA)"
Cohesion: 0.20
Nodes (10): Cassan 2008 curvilinear abscissa, Charbonneau 1995 genetic algorithm, Kains 2009 reference application, Planned 2L1S roadmap (rescale -> finite source -> Cassan -> grid -> GA), SIGNALMEN anomaly-vs-outlier test (Dominik et al. 2007), Cassan (2008) caustic-crossing parametrisation, Config-driven pipeline (search.py + event.py + TOML), Huber loss likelihood (DELTA=1.345) (+2 more)

### Community 22 - "cross_check_mulensmodel_fit.py"
Cohesion: 0.33
Nodes (8): chi2(), mm_magnification(), plot_fit(), One-time cross-check: fit MulensModel's own binary-lens model to the same…, Chi2 of MulensModel's binary magnification against the combined OGLE+MOA data., Serial multi-start search, then a tight-tolerance refit of the best seed., Overlay the MulensModel fit on the OGLE+MOA data, full baseline + zoomed on the…, run_fit()

### Community 24 - "Per-instrument error-bar rescaling k"
Cohesion: 0.40
Nodes (5): MOA-2019-BLG-008 error bars underestimated 16-49x, Per-instrument error-bar rescaling k, M-19-BLG008 (MOA-2019-BLG-008, KMT I subset), error_rescaling(), Rescales the error bars per instrument by settings chi2/dof ~ 1. raw_resid, err…

### Community 26 - "distinct_modes"
Cohesion: 0.67
Nodes (3): distinct_modes tolerance-based dedup, distinct_modes(), Refined (x, chi2, ...) sorted by chi2, and whether each crosses a caustic ->…

## Knowledge Gaps
- **10 isolated node(s):** `submit_all.sh script`, `astropy`, `corner`, `emcee`, `torch` (+5 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 207 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **3 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `CHANGELOG session log` connect `CHANGELOG session log` to `load_instrument`, `CLAUDE.md project guide`, `diagnose`, `Dataset short-name mapping`, `lc_models.py`, `run_grid`, `README.md`, `Planned 2L1S roadmap (rescale -> finite source -> Cassan -> grid -> GA)`, `Per-instrument error-bar rescaling k`, `distinct_modes`?**
  _High betweenness centrality (0.101) - this node is a cross-community bridge._
- **Why does `CLAUDE.md project guide` connect `CLAUDE.md project guide` to `mcmc_fit.py`, `CHANGELOG session log`, `Dataset short-name mapping`, `run_grid`, `README.md`, `Planned 2L1S roadmap (rescale -> finite source -> Cassan -> grid -> GA)`, `Per-instrument error-bar rescaling k`?**
  _High betweenness centrality (0.089) - this node is a cross-community bridge._
- **Why does `find_zoom_window()` connect `mcmc_fit.py` to `mcmc_fit_2l1s.py`, `fit_2l1s.py`, `search.py`, `plot_fit`, `CHANGELOG session log`, `binary_trajectory`, `ld_test.py`, `numpy`?**
  _High betweenness centrality (0.036) - this node is a cross-community bridge._
- **Are the 16 inferred relationships involving `Event` (e.g. with `anomaly_times()` and `bic()`) actually correct?**
  _`Event` has 16 INFERRED edges - model-reasoned connections that need verification._
- **What connects `submit_all.sh script`, `astropy`, `corner` to the rest of the system?**
  _10 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `mcmc_fit_2l1s.py` be split into smaller, more focused modules?**
  _Cohesion score 0.06321334503950835 - nodes in this community are weakly interconnected._
- **Should `fit_2l1s.py` be split into smaller, more focused modules?**
  _Cohesion score 0.07712765957446809 - nodes in this community are weakly interconnected._