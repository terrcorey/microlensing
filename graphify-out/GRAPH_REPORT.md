# Graph Report - microlensing  (2026-10-12)

## Corpus Check
- 36 files · ~65,187 words
- Verdict: corpus is large enough that graph structure adds value.
- Unclassified: 58 file(s) not represented in the graph (top: .npz 12, .out 12, .toml 9)

## Summary
- 566 nodes · 1446 edges · 32 communities (27 shown, 5 thin omitted)
- Extraction: 95% EXTRACTED · 5% INFERRED · 0% AMBIGUOUS · INFERRED: 68 edges (avg confidence: 0.87)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `9aa909b1`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- mcmc_fit_2l1s.py
- compare_2l1s_fits.py
- mcmc_fit.py
- event.py
- as_binary
- bic
- CHANGELOG session log
- variants2.sh
- CLAUDE.md project guide
- plot_fit
- search.py
- lc_models.py
- cross_check_mulensmodel_fs.py
- chi2_fspl
- ld_test.py
- variants.sh
- Dataset short-name mapping
- _binary_images_batch
- fit_2l1s.py
- chi2_binary
- README.md
- Planned 2L1S roadmap (rescale -> finite source -> Cassan -> grid -> GA)
- plot_fit
- binary_trajectory
- Per-instrument error-bar rescaling k
- derive_binary_quintic.py
- submit_all.sh
- load_coords
- torch
- compare_pspl_2l1s.py
- run_grid
- download_data.py

## God Nodes (most connected - your core abstractions)
1. `CHANGELOG session log` - 44 edges
2. `CLAUDE.md project guide` - 35 edges
3. `Event` - 29 edges
4. `diagnose()` - 26 edges
5. `save_fig()` - 26 edges
6. `plot_fit()` - 22 edges
7. `binary_trajectory()` - 21 edges
8. `magnification()` - 19 edges
9. `TwoL1SParams` - 19 edges
10. `plot_fit()` - 19 edges

## Surprising Connections (you probably didn't know these)
- `u0 sign degeneracy / required u0_guess` --references--> `run_ogle_only_diagnostic()`  [EXTRACTED]
  CLAUDE.md → mcmc_fit_binary.py
- `TwoL1SParams canonical param order` --references--> `TwoL1SParams`  [EXTRACTED]
  CLAUDE.md → scratch/fit_2l1s.py
- `find_zoom_window significance handling` --references--> `find_zoom_window()`  [EXTRACTED]
  CHANGELOG.md → zoom_utils.py
- `MCMC chain checkpointing every 200 steps` --semantically_similar_to--> `Grid checkpoint/resume (grid_partial.npz)`  [INFERRED] [semantically similar]
  CHANGELOG.md → CLAUDE.md
- `Geocentric JD light-travel correction` --references--> `load_instrument()`  [EXTRACTED]
  CHANGELOG.md → event.py

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **search.py config-driven stages (FSPL -> grid -> refine/MCMC -> BIC -> diagnose)** — search_fit_fspl, search_run_grid, search_refine, search_run_mcmc, search_bic, search_diagnose [EXTRACTED 1.00]
- **Session 22 robustness fixes against VBBL hangs** — changelog_vbbl_hang_pool_fix, changelog_rho_min_guard, changelog_mcmc_checkpointing, changelog_config_hash_cache, changelog_check_search_selfcheck [EXTRACTED 1.00]
- **Likelihood choices (Gaussian chi2, Student-t, Huber)** — claude_student_t_likelihood, claude_huber_likelihood, scratch_fit_2l1s_chi2 [INFERRED 0.85]

## Communities (32 total, 5 thin omitted)

### Community 0 - "mcmc_fit_2l1s.py"
Cohesion: 0.06
Nodes (72): EnsembleSampler, Inverse of cassan_to_standard(): every (sigma, t) where the parallax-free…, standard_to_cassan(), flat_chain(), ndarray, Corner plot of the raw MCMC parameters -- shared by every fit script.…, Flat (samples, log_probs) after burn-in/thinning -- shared by every emcee…, save_corner() (+64 more)

### Community 1 - "compare_2l1s_fits.py"
Cohesion: 0.17
Nodes (14): binary_magnification(), binary_magnification_fs(), _disk_average(), _sample_coords(), Total binary-lens magnification at source position(s) `zeta`. Sums 1/|J| over…, Finite-source binary-lens magnification at source position(s) `zeta`. `rho` is…, moa_to_magnification(), Invert MOA's flux model to recover A(t): A = 1 + F / fs_moa. (+6 more)

### Community 2 - "mcmc_fit.py"
Cohesion: 0.06
Nodes (69): argparse, Annual parallax (Gould 2004 geocentric), FitResult NamedTuple (samples + labels), corner, emcee, flux(), magnification(), magnitude() (+61 more)

### Community 3 - "event.py"
Cohesion: 0.06
Nodes (54): Geocentric JD light-travel correction, concurrent_futures_process, chi2(), error_scale(), Event, flux_residuals(), Instrument, ld_penalty() (+46 more)

### Community 4 - "as_binary"
Cohesion: 0.14
Nodes (17): distinct_modes tolerance-based dedup, MCMC chain checkpointing every 200 steps, as_binary(), distinct_modes(), from_mcmc(), ld_steps(), log_prob(), Refined (x, chi2, ...) sorted by chi2, and whether each crosses a caustic ->… (+9 more)

### Community 5 - "bic"
Cohesion: 0.50
Nodes (4): BIC caveats (look-elsewhere, correlated errors, circular K), PSPL vs 2L1S BIC comparison, bic(), chi2 + k ln N, with k counting the profiled flux parameters too.

### Community 6 - "CHANGELOG session log"
Cohesion: 0.15
Nodes (22): CHANGELOG session log, Autocorrelation (tau) convergence check, Blending degeneracy (model-dependent magnification), scratch/check_search.py self-check, Prior-based lens mass estimate, Fair fit_lc.png in OGLE I magnitude, find_zoom_window significance handling, Gould et al. 2006 (O-05-BLG169 discovery) (+14 more)

### Community 7 - "variants2.sh"
Cohesion: 0.67
Nodes (3): PYTHONWARNINGS, run(), variants2.sh script

### Community 8 - "CLAUDE.md project guide"
Cohesion: 0.14
Nodes (18): SIGNALMEN anomaly-vs-outlier test (Dominik et al. 2007), CLAUDE.md project guide, Blind search, no literature seeds, O-03-BLG235 calibrate-once-then-fit pipeline, Huber loss likelihood (DELTA=1.345), No compute on shared login node (srun/sbatch), MulensModel cross-check oracle, Point-source (s,q) grid (BinaryMag0, constant cost) (+10 more)

### Community 9 - "plot_fit"
Cohesion: 0.22
Nodes (11): chi2(), _fit_one_seed(), plot_fit(), plot_panel(), Data + 2L1S fit overlay (full baseline + auto-zoomed peak, magnitude space)…, Closed-form weighted least squares for (fs, fb) at fixed magnification shape --…, Chi2 of the 2L1S model against KMT I-band flux, profiling fs/fb at every trial., Nelder-Mead on one seed at loose tolerance; module-level so it can be pickled. (+3 more)

### Community 10 - "search.py"
Cohesion: 0.08
Nodes (34): Concurrent MCMC modes on shared Pool, diagnose summary.txt readable output, datetime, hashlib, importlib_metadata, fspl_magnification(), Finite-source point-lens A(u) via VBM's ESPLMag2. `rho` is the source radius in…, mpl_toolkits_axes_grid1_inset_locator (+26 more)

### Community 11 - "lc_models.py"
Cohesion: 0.14
Nodes (16): astropy_coordinates, astropy_time, astropy_units, mag_to_flux(), _project(), Functional forms used to model a microlensing light curve. Point-source point-…, Calibrated magnitude + error -> flux + error on the ZERO_POINT_MAG scale., Project a Cartesian position/velocity vector onto the sky's (N, E) tangent-… (+8 more)

### Community 12 - "cross_check_mulensmodel_fs.py"
Cohesion: 0.20
Nodes (8): caustic_point(), Point at abscissa sigma (wraps with period 1) on a cassan_caustic() output., mulensmodel, One-time offline cross-check of lc_models.binary_magnification_fs against…, mm_magnification(), One-time offline cross-check of lc_models.binary_magnification against…, MulensModel's magnification at the same complex zeta values. MulensModel's docs…, report()

### Community 13 - "chi2_fspl"
Cohesion: 0.21
Nodes (18): nelder_mead(), Nelder-Mead from an explicit initial simplex: scipy's default steps 5% of each…, as_fspl(), chi2_fspl(), fit_fspl(), FSPL, fspl_A(), fspl_cell() (+10 more)

### Community 14 - "ld_test.py"
Cohesion: 0.14
Nodes (20): binary_magnification_vbm(), Finite-source 2L1S A(zeta) via VBM's BinaryMag2 -- `zeta` in…, polish(), nelder_mead() restarted from its own result with a fresh simplex (`steps(x)`)…, fit(), Free limb darkening test (session 23, one-off): does fitting the linear…, (FSPL, 2L1S overall best, K at FSPL) from a search.py summary.txt (alpha…, Pool task: (model, free LD?, start, event) -> (x, chi2); LD appended to x when… (+12 more)

### Community 16 - "Dataset short-name mapping"
Cohesion: 0.13
Nodes (18): Bennett 2015 re-reduction of O-05-BLG169, Correlated MDM residuals in O-05-BLG169 2L1S, Offset-scaled archive magnitudes (FTN, MDM, Auckland), Grid inner search N_SIGMA/N_POLISH/STD_MAXFEV 16/5/600, Why the coarse grid missed Bond's basin, O-03-BLG235 search failing done-checks, Bond et al. 2004 published O-03-BLG235 solution, Drop-one instrument sensitivity runs (+10 more)

### Community 17 - "_binary_images_batch"
Cohesion: 0.18
Nodes (12): Symbolically derived binary-lens quintic, Torch batched eigvals root solver, binary_images(), _binary_images_batch(), _companion_eigvals(), _quintic_coefficients(), Coefficients (highest degree first) of the degree-5 polynomial in z equivalent…, (torch, device), imported on first use: only the hand-rolled solver needs it,… (+4 more)

### Community 18 - "fit_2l1s.py"
Cohesion: 0.15
Nodes (14): concurrent_futures, functools, load_raw(), Raw OGLE (mag) + MOA (differential flux) tables, time-shifted to HJD-2450000., _data_delta_s(), delta_s(), plain_pspl(), model_fn() (+6 more)

### Community 19 - "chi2_binary"
Cohesion: 0.16
Nodes (15): cassan_to_standard(), Cassan (2008) caustic-crossing parameters -> binary_trajectory()'s (t0, u0, tE,…, cell(), chi2_binary(), grid_cell(), c2(), cassan_f(), polish() (+7 more)

### Community 20 - "README.md"
Cohesion: 0.22
Nodes (8): Mass-distance degeneracy in lens mass estimate, Prior-based lens mass estimate (D_L uniform, mu_rel log-normal), Flat-bound MCMC priors, astropy, corner, emcee, torch, VBBinaryLensing 3.7.0

### Community 21 - "Planned 2L1S roadmap (rescale -> finite source -> Cassan -> grid -> GA)"
Cohesion: 0.29
Nodes (7): Cassan 2008 curvilinear abscissa, Charbonneau 1995 genetic algorithm, Kains 2009 reference application, Planned 2L1S roadmap (rescale -> finite source -> Cassan -> grid -> GA), Cassan (2008) caustic-crossing parametrisation, Config-driven pipeline (search.py + event.py + TOML), Near-miss cusp solutions problem

### Community 22 - "plot_fit"
Cohesion: 0.33
Nodes (7): chi2(), mm_magnification(), plot_fit(), Chi2 of MulensModel's binary magnification against the combined OGLE+MOA data., Serial multi-start search, then a tight-tolerance refit of the best seed., Overlay the MulensModel fit on the OGLE+MOA data, full baseline + zoomed on the…, run_fit()

### Community 23 - "binary_trajectory"
Cohesion: 0.17
Nodes (16): complex128, binary_trajectory(), cassan_caustic(), caustic_curve(), equidistant_caustic(), lens_position(), NDArray, Resample one closed caustic (an element of caustic_curve()'s list) to points… (+8 more)

### Community 24 - "Per-instrument error-bar rescaling k"
Cohesion: 0.40
Nodes (5): MOA-2019-BLG-008 error bars underestimated 16-49x, Per-instrument error-bar rescaling k, M-19-BLG008 (MOA-2019-BLG-008, KMT I subset), error_rescaling(), Rescales the error bars per instrument by settings chi2/dof ~ 1. raw_resid, err…

### Community 30 - "load_coords"
Cohesion: 0.50
Nodes (4): load_coords(), _parse_header_coord(), Extract (RA, Dec) sexagesimal strings from a raw table's own \\RA/\\DEC header…, Target coordinates, parsed from OGLE_PATH/MOA_PATH's own headers rather than…

### Community 34 - "compare_pspl_2l1s.py"
Cohesion: 0.15
Nodes (20): astropy_timeseries, itertools, matplotlib_pyplot, numpy, pathlib, scipy_optimize, PSPL vs 2L1S model comparison (BIC) for OGLE-2003-BLG-235/MOA-2003-BLG-53. Both…, One-time cross-check: fit MulensModel's own binary-lens model to the same… (+12 more)

### Community 35 - "run_grid"
Cohesion: 0.20
Nodes (10): Cache keyed on config hash + inner-search settings, Grid checkpoint/resume (grid_partial.npz), init(), axis(), grid_timing(), _init(), [lo, hi, step] -> lo..hi inclusive., (ii) delta-chi2(s, q) map, cached at `path`, + the timing of the cells this run… (+2 more)

## Knowledge Gaps
- **11 isolated node(s):** `PYTHONWARNINGS`, `submit_all.sh script`, `astropy`, `corner`, `emcee` (+6 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 219 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **5 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `CHANGELOG session log` connect `CHANGELOG session log` to `run_grid`, `as_binary`, `bic`, `event.py`, `CLAUDE.md project guide`, `search.py`, `Dataset short-name mapping`, `_binary_images_batch`, `README.md`, `Planned 2L1S roadmap (rescale -> finite source -> Cassan -> grid -> GA)`, `Per-instrument error-bar rescaling k`?**
  _High betweenness centrality (0.092) - this node is a cross-community bridge._
- **Why does `CLAUDE.md project guide` connect `CLAUDE.md project guide` to `mcmc_fit.py`, `run_grid`, `bic`, `CHANGELOG session log`, `Dataset short-name mapping`, `README.md`, `Planned 2L1S roadmap (rescale -> finite source -> Cassan -> grid -> GA)`, `Per-instrument error-bar rescaling k`?**
  _High betweenness centrality (0.081) - this node is a cross-community bridge._
- **Why does `save_fig()` connect `compare_pspl_2l1s.py` to `mcmc_fit_2l1s.py`, `compare_2l1s_fits.py`, `mcmc_fit.py`, `event.py`, `plot_fit`, `search.py`, `chi2_fspl`, `ld_test.py`, `fit_2l1s.py`?**
  _High betweenness centrality (0.038) - this node is a cross-community bridge._
- **Are the 18 inferred relationships involving `Event` (e.g. with `anomaly_times()` and `bic()`) actually correct?**
  _`Event` has 18 INFERRED edges - model-reasoned connections that need verification._
- **What connects `PYTHONWARNINGS`, `submit_all.sh script`, `astropy` to the rest of the system?**
  _11 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `mcmc_fit_2l1s.py` be split into smaller, more focused modules?**
  _Cohesion score 0.055905220288781934 - nodes in this community are weakly interconnected._
- **Should `mcmc_fit.py` be split into smaller, more focused modules?**
  _Cohesion score 0.05516431924882629 - nodes in this community are weakly interconnected._