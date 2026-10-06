# Graph Report - microlensing  (2026-10-06)

## Corpus Check
- Corpus is ~44,604 words - fits in a single context window. You may not need a graph.

## Summary
- 406 nodes · 972 edges · 10 communities (8 shown, 2 thin omitted)
- Extraction: 95% EXTRACTED · 5% INFERRED · 0% AMBIGUOUS · INFERRED: 46 edges (avg confidence: 0.87)
- Token cost: 80,871 input · 0 output

## Community Hubs (Navigation)
- Event Loading & Flux Calibration
- Binary-Lens Magnification
- 2L1S Scratch Fits
- Data, Parallax & PSPL Model
- Docs, Papers & Concepts
- Cassan & Binary Preprocessing
- Parallax PSPL MCMC
- PSPL Core & Joint Fit
- Quintic Derivation
- SLURM Submission

## God Nodes (most connected - your core abstractions)
1. `plot_fit()` - 21 edges
2. `Event` - 20 edges
3. `magnification()` - 19 edges
4. `binary_trajectory()` - 18 edges
5. `TwoL1SParams` - 18 edges
6. `trajectory()` - 17 edges
7. `caustic_curve()` - 16 edges
8. `flat_chain()` - 16 edges
9. `run_joint_fit()` - 15 edges
10. `binary_magnification_fs()` - 14 edges

## Surprising Connections (you probably didn't know these)
- `plot_panel()` --calls--> `model_fn()`  [INFERRED]
  zoom_utils.py → scratch/fit_2l1s.py
- `_model_magnitude()` --calls--> `magnitude()`  [EXTRACTED]
  mcmc_fit.py → lc_models.py
- `Config-driven pipeline (search.py + event.py + input TOML)` --references--> `emcee`  [INFERRED]
  CLAUDE.md → requirements.txt
- `anomaly_times()` --uses--> `Event`  [INFERRED]
  search.py → event.py
- `bic()` --uses--> `Event`  [INFERRED]
  search.py → event.py

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Config-driven blind search flow (FSPL -> (s,q) grid -> refine -> MCMC -> BIC)** — claude_config_driven_pipeline, changelog_s_q_grid_search, claude_cassan_parametrisation, claude_profiled_flux_calibration, claude_error_bar_rescaling, changelog_autocorrelation_check, claude_bic_comparison [EXTRACTED 1.00]
- **Alternative likelihoods for heavy-tailed residuals** — claude_student_t_likelihood, claude_huber_likelihood, claude_error_bar_rescaling [INFERRED 0.75]
- **Effects that break the tE mass-distance degeneracy** — claude_mass_distance_degeneracy, claude_annual_parallax, claude_finite_source, readme_estimate_mass_priors [INFERRED 0.85]

## Communities (10 total, 2 thin omitted)

### Community 0 - "Event Loading & Flux Calibration"
Cohesion: 0.06
Nodes (76): chi2(), error_scale(), Event, flux_residuals(), Instrument, load_event(), load_instrument(), profile_flux() (+68 more)

### Community 1 - "Binary-Lens Magnification"
Cohesion: 0.05
Nodes (61): astropy_time, astropy_units, complex128, binary_images(), _binary_images_batch(), binary_magnification(), binary_magnification_fs(), _disk_average() (+53 more)

### Community 2 - "2L1S Scratch Fits"
Cohesion: 0.06
Nodes (67): concurrent_futures, functools, itertools, Corner plot of the raw MCMC parameters -- shared by every fit script.…, save_corner(), mpl_toolkits_axes_grid1_inset_locator, multiprocessing, os (+59 more)

### Community 3 - "Data, Parallax & PSPL Model"
Cohesion: 0.06
Nodes (54): argparse, astropy_coordinates, corner, Download raw photometry for each event we're working with.…, emcee, magnitude(), Returns (delta_s_n, delta_s_e): Earth's sky-projected position relative to the…, Observed magnitude, for overplotting against real photometry. (+46 more)

### Community 4 - "Docs, Papers & Concepts"
Cohesion: 0.07
Nodes (49): MCMC autocorrelation convergence check (nsteps/tau > 50), Bennett et al. 2015, Bond et al. 2004, Bozza 2010, Cassan 2008, Charbonneau 1995, Dominik et al. 2007, FTN negative-offset flux scale (free-sign fb) (+41 more)

### Community 5 - "Cassan & Binary Preprocessing"
Cohesion: 0.13
Nodes (24): Inverse of cassan_to_standard(): every (sigma, t) where the parallax-free…, standard_to_cassan(), moa_to_magnification(), ogle_to_magnification(), Invert OGLE's flux model to recover A(t): A = (F - fb_ogle) / fs_ogle., Invert MOA's flux model to recover A(t): A = 1 + F / fs_moa., CassanParams, chi2_cassan() (+16 more)

### Community 6 - "Parallax PSPL MCMC"
Cohesion: 0.11
Nodes (24): EnsembleSampler, flux(), plain_flux(), plain_magnitude(), Observed flux: magnified source plus a constant blend flux., flux() with piE_N=piE_E=0 -- the exact pre-parallax model, used only to…, magnitude() with piE_N=piE_E=0 -- see plain_flux()., fit_parallax_pspl_mcmc() (+16 more)

### Community 7 - "PSPL Core & Joint Fit"
Cohesion: 0.19
Nodes (19): magnification(), Source-lens separation u(t), in units of the Einstein radius., Paczynski point-lens magnification A(u)., trajectory(), run_joint_fit(), chi2(), log_likelihood(), log_prior() (+11 more)

## Knowledge Gaps
- **2 isolated node(s):** `submit_all.sh script`, `corner`
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 155 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **2 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `binary_trajectory()` connect `Binary-Lens Magnification` to `Event Loading & Flux Calibration`, `2L1S Scratch Fits`, `Cassan & Binary Preprocessing`?**
  _High betweenness centrality (0.034) - this node is a cross-community bridge._
- **Why does `flat_chain()` connect `Parallax PSPL MCMC` to `Event Loading & Flux Calibration`, `2L1S Scratch Fits`, `Data, Parallax & PSPL Model`, `Cassan & Binary Preprocessing`, `PSPL Core & Joint Fit`?**
  _High betweenness centrality (0.034) - this node is a cross-community bridge._
- **Why does `caustic_curve()` connect `Binary-Lens Magnification` to `Event Loading & Flux Calibration`, `2L1S Scratch Fits`, `Cassan & Binary Preprocessing`?**
  _High betweenness centrality (0.024) - this node is a cross-community bridge._
- **Are the 11 inferred relationships involving `Event` (e.g. with `anomaly_times()` and `bic()`) actually correct?**
  _`Event` has 11 INFERRED edges - model-reasoned connections that need verification._
- **What connects `submit_all.sh script`, `corner` to the rest of the system?**
  _2 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Event Loading & Flux Calibration` be split into smaller, more focused modules?**
  _Cohesion score 0.05527805527805528 - nodes in this community are weakly interconnected._
- **Should `Binary-Lens Magnification` be split into smaller, more focused modules?**
  _Cohesion score 0.050724637681159424 - nodes in this community are weakly interconnected._