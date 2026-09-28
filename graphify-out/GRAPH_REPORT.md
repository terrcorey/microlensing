# Graph Report - microlensing  (2026-09-25)

## Corpus Check
- 8 files · ~31,296 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 334 nodes · 698 edges · 18 communities (14 shown, 4 thin omitted)
- Extraction: 93% EXTRACTED · 7% INFERRED · 0% AMBIGUOUS · INFERRED: 50 edges (avg confidence: 0.83)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- Lens Model Physics & Parallax
- PSPL Pipeline & Mass Estimation
- 2L1S Search & MOA-2019-BLG-008
- 2L1S Likelihoods (chi2/Student-t/Huber)
- Cassan Caustic Parametrisation
- Caustic Search Roadmap & Calibration Bugs
- Parallax + Student-t PSPL MCMC
- Datasets & Dependencies
- MulensModel Fit Cross-Check
- Target Coordinate Parsing
- Caustic Spike Anomaly Test
- Repo Rules & Collaboration
- Workflow Instructions & Direction
- Raw Table Loading
- Data Download
- requirements.txt (isolated)
- M-19-BLG008 short name (isolated)

## God Nodes (most connected - your core abstractions)
1. `TwoL1SParams` - 20 edges
2. `plot_fit()` - 20 edges
3. `binary_magnification()` - 17 edges
4. `chi2()` - 14 edges
5. `run_joint_fit()` - 13 edges
6. `fit_parallax_pspl_mcmc()` - 13 edges
7. `binary_trajectory()` - 13 edges
8. `caustic_curve()` - 13 edges
9. `residuals()` - 13 edges
10. `run_ogle_only_diagnostic()` - 11 edges

## Surprising Connections (you probably didn't know these)
- `TwoL1SParams canonical param-order refactor` --semantically_similar_to--> `FitResult`  [INFERRED] [semantically similar]
  CLAUDE.md → mcmc_fit.py
- `t0_par reference-epoch convention` --implements--> `get_t0_par()`  [EXTRACTED]
  CLAUDE.md → mcmc_fit.py
- `HJD - 2450000 time convention` --conceptually_related_to--> `find_zoom_window()`  [INFERRED]
  CLAUDE.md → zoom_utils.py
- `Cassan curvilinear-abscissa caustic parametrisation` --implements--> `cassan_caustic()`  [INFERRED]
  CHANGELOG.md → lc_models.py
- `Cassan curvilinear-abscissa caustic parametrisation` --implements--> `cassan_to_standard()`  [INFERRED]
  CHANGELOG.md → lc_models.py

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Recurring Wrong-But-Plausible Convergence Pattern** — changelog_session_1, changelog_session_3, changelog_session_4, changelog_session_5, changelog_session_6, changelog_wrong_but_plausible_pattern [EXTRACTED 1.00]
- **2L1S Search Roadmap (rescale -> rho* -> Cassan -> grid -> GA)** — changelog_error_bar_rescaling, changelog_finite_source_effects, changelog_cassan_parametrisation, changelog_dq_grid_search, changelog_genetic_algorithm, changelog_signalmen_anomaly_test [EXTRACTED 1.00]
- **Likelihood Options for 2L1S (chi2 / Student-t / Huber)** — claude_studentt_likelihood, claude_huber_loss, changelog_chi2_vs_studentt_comparison, changelog_all_statistics_fail_spike [INFERRED 0.85]

## Communities (18 total, 4 thin omitted)

### Community 0 - "Lens Model Physics & Parallax"
Cohesion: 0.06
Nodes (53): astropy_coordinates, astropy_time, astropy_units, MulensModel cross-check oracle, PSPL (Paczynski point-source point-lens) model, functools, itertools, flux() (+45 more)

### Community 1 - "PSPL Pipeline & Mass Estimation"
Cohesion: 0.06
Nodes (51): argparse, BIC PSPL-vs-2L1S verdict (deltaBIC>10), Prior-based lens-mass estimate, Han & Gould 2003 (disk/bulge-weighted Galactic prior), Kass & Raftery (BIC decisive-threshold convention), Mass-distance-velocity degeneracy, Session 1 — Two Real Events Set Up End-to-End, Session 2 — Derived Physical Quantities & Lens-Mass Estimation (+43 more)

### Community 2 - "2L1S Search & MOA-2019-BLG-008"
Cohesion: 0.06
Nodes (48): Bachelet et al. 2022 (arXiv:2205.07522), Batched companion-matrix eigensolve (torch), Error-bar rescaling / additive jitter plan, KMT mag_err underestimate (~16-49x), MOA-2019-BLG-008, Session 3 — 2L1S Magnification Implemented & Validated, Session 4 — 2L1S Search Stuck in Wrong Basin, Session 5 — Flux-Calibration Bug Found (+40 more)

### Community 3 - "2L1S Likelihoods (chi2/Student-t/Huber)"
Cohesion: 0.09
Nodes (43): chi2 vs Student-t 2L1S MCMC comparison, Session 10 — Parallax/Student-t Extended to 2L1S Track, Session 11 — TwoL1SParams Refactor & Huber Loss, Session 9 — Parallax & Student-t Propagated to O-03-BLG235, Huber loss likelihood (2L1S track), Student-t robust likelihood, scratch/submit_mcmc_2l1s.sbatch (SLURM CPU job), TwoL1SParams canonical param-order refactor (+35 more)

### Community 4 - "Cassan Caustic Parametrisation"
Cohesion: 0.10
Nodes (31): Caustic branch ordering via linear_sum_assignment, Arc-length equidistant caustic resampling, scratch/2l1s/{method}/ output layout, Session 12 — Cassan Caustic Parametrisation Started, cassan_caustic(), cassan_to_standard(), caustic_curve(), caustic_point() (+23 more)

### Community 5 - "Caustic Search Roadmap & Calibration Bugs"
Cohesion: 0.11
Nodes (25): Bond et al. 2004 (arXiv:astro-ph/0404309), Cassan 2008 (caustic curvilinear abscissa), Cassan curvilinear-abscissa caustic parametrisation, Charbonneau 1995 (genetic algorithm), (d, q) grid search over caustic geometry, torch fork-after-threading deadlock, Frozen PSPL flux-calibration bias, Genetic-algorithm caustic search (+17 more)

### Community 6 - "Parallax + Student-t PSPL MCMC"
Cohesion: 0.13
Nodes (19): /improve-codebase-architecture review (FitResult, get_t0_par), Session 7 — Annual Parallax + Student-t Roadmap Scoped, Session 8 — Parallax & Student-t Implemented for O-05-BLG086, Annual parallax (Gould 2004 geocentric formalism), Gould 2004 (annual parallax geocentric formalism), t0_par reference-epoch convention, fit_parallax_pspl_mcmc(), log_likelihood() (+11 more)

### Community 7 - "Datasets & Dependencies"
Cohesion: 0.18
Nodes (9): O-03-BLG235 (OGLE-2003-BLG-235 / MOA-2003-BLG-53), O-05-BLG086 (OGLE-2005-BLG-086), astropy, corner, emcee, matplotlib, numpy, scipy (+1 more)

### Community 8 - "MulensModel Fit Cross-Check"
Cohesion: 0.33
Nodes (7): chi2(), mm_magnification(), plot_fit(), Chi2 of MulensModel's binary magnification against the combined OGLE+MOA data., Serial multi-start search, then a tight-tolerance refit of the best seed., Overlay the MulensModel fit on the OGLE+MOA data, full baseline + zoomed on the…, run_fit()

### Community 9 - "Target Coordinate Parsing"
Cohesion: 0.40
Nodes (5): O-03-BLG235 hardcoded RA bug (load_coords fix), load_coords(), _parse_header_coord(), Extract (RA, Dec) sexagesimal strings from a raw table's own \\RA/\\DEC header…, Target coordinates, parsed from OGLE_PATH/MOA_PATH's own headers rather than…

### Community 10 - "Caustic Spike Anomaly Test"
Cohesion: 0.50
Nodes (4): All statistics fail caustic spike finding, Dominik et al. 2007 (SIGNALMEN, arXiv:0706.2566), Finite source effects (rho*), SIGNALMEN anomaly-vs-outlier test

### Community 11 - "Repo Rules & Collaboration"
Cohesion: 0.50
Nodes (4): User-writes/Claude-guides collaboration mode, Output layout convention (scratch/ vs. pipeline dirs), ponytail skill (minimal-solution rule), Repository rules section

### Community 12 - "Workflow Instructions & Direction"
Cohesion: 0.50
Nodes (4): /graphify knowledge-graph skill, /grill-me consensus mechanism, Repository instructions section, Config-driven pipeline long-term direction

### Community 13 - "Raw Table Loading"
Cohesion: 0.67
Nodes (3): NASA Exoplanet Archive table format, load_raw(), Raw OGLE (mag) + MOA (differential flux) tables, time-shifted to HJD-2450000.

## Knowledge Gaps
- **29 isolated node(s):** `M-19-BLG008 (MOA-2019-BLG-008, KMT I-band subset)`, `O-03-BLG235 (OGLE-2003-BLG-235 / MOA-2003-BLG-53)`, `O-05-BLG086 (OGLE-2005-BLG-086)`, `astropy`, `emcee` (+24 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 134 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **4 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `PSPL (Paczynski point-source point-lens) model` connect `Lens Model Physics & Parallax` to `PSPL Pipeline & Mass Estimation`, `Datasets & Dependencies`?**
  _High betweenness centrality (0.067) - this node is a cross-community bridge._
- **Why does `TwoL1SParams` connect `2L1S Likelihoods (chi2/Student-t/Huber)` to `Lens Model Physics & Parallax`, `2L1S Search & MOA-2019-BLG-008`, `Cassan Caustic Parametrisation`, `Caustic Search Roadmap & Calibration Bugs`?**
  _High betweenness centrality (0.062) - this node is a cross-community bridge._
- **Why does `corner` connect `Datasets & Dependencies` to `Lens Model Physics & Parallax`?**
  _High betweenness centrality (0.057) - this node is a cross-community bridge._
- **Are the 3 inferred relationships involving `chi2()` (e.g. with `binary_magnification()` and `_fit_one_seed()`) actually correct?**
  _`chi2()` has 3 INFERRED edges - model-reasoned connections that need verification._
- **Are the 3 inferred relationships involving `run_joint_fit()` (e.g. with `chi2()` and `log_probability()`) actually correct?**
  _`run_joint_fit()` has 3 INFERRED edges - model-reasoned connections that need verification._
- **What connects `M-19-BLG008 (MOA-2019-BLG-008, KMT I-band subset)`, `O-03-BLG235 (OGLE-2003-BLG-235 / MOA-2003-BLG-53)`, `O-05-BLG086 (OGLE-2005-BLG-086)` to the rest of the system?**
  _29 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Lens Model Physics & Parallax` be split into smaller, more focused modules?**
  _Cohesion score 0.060655737704918035 - nodes in this community are weakly interconnected._