# Graph Report - microlensing  (2026-09-23)

## Corpus Check
- 6 files · ~26,040 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 264 nodes · 521 edges · 14 communities (8 shown, 6 thin omitted)
- Extraction: 96% EXTRACTED · 4% INFERRED · 0% AMBIGUOUS · INFERRED: 23 edges (avg confidence: 0.85)
- Token cost: 0 input · 142,386 output

## Community Hubs (Navigation)
- PSPL Fitting Pipeline
- Lens Model Physics
- 2L1S Trajectory & Parallax
- 2L1S Validation & Comparison
- Session 10 Findings & Bugs
- Early Sessions & Datasets
- MCMC Likelihood Functions
- Data Loading & Calibration
- Grill-Me Workflow
- 2L1S Model (concept)
- Output Layout Convention
- MOA-2019-BLG-008 Dataset
- SLURM Submission Script

## God Nodes (most connected - your core abstractions)
1. `binary_magnification()` - 19 edges
2. `plot_fit()` - 17 edges
3. `run_joint_fit()` - 14 edges
4. `binary_trajectory()` - 14 edges
5. `chi2()` - 13 edges
6. `fit_parallax_pspl_mcmc()` - 11 edges
7. `sun_earth_projection()` - 10 edges
8. `run_ogle_only_diagnostic()` - 9 edges
9. `plot_fit_panels()` - 9 edges
10. `plot_fit()` - 9 edges

## Surprising Connections (you probably didn't know these)
- `O-03-BLG235 calibrate-once-then-fit pipeline` --rationale_for--> `run_joint_fit()`  [EXTRACTED]
  CLAUDE.md → mcmc_fit_binary.py
- `estimate_mass()` --implements--> `Prior-based lens mass estimate (M_lens)`  [EXTRACTED]
  mcmc_fit.py → README.md
- `fit_parallax_pspl_mcmc()` --references--> `emcee`  [INFERRED]
  mcmc_fit.py → CHANGELOG.md
- `binary_images()` --references--> `torch (PyTorch)`  [EXTRACTED]
  lc_models.py → CHANGELOG.md
- `plot_fit() tag="" naming-collision gap` --rationale_for--> `plot_fit()`  [EXTRACTED]
  CLAUDE.md → scratch/fit_2l1s.py

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Shared plotting/MCMC helper functions (do not re-copy)** — zoom_utils_plot_fit_panels, zoom_utils_plot_residual_panel, zoom_utils_plot_residual_hist, mcmc_fit_save_corner, mcmc_fit_plot_fit_lc, mcmc_fit_fitresult, mcmc_fit_get_t0_par [EXTRACTED 1.00]
- **Wrong-but-plausible bugs recurring across sessions** — changelog_u0_sign_bug, changelog_mulensmodel_alpha_convention_bug, changelog_coordinate_bug, changelog_flux_calibration_bug [INFERRED 0.85]
- **Real microlensing events studied for brown-dwarf-mass characterization** — dataset_names_o05blg086, dataset_names_o03blg235, changelog_moa_2019_blg008 [INFERRED 0.75]

## Communities (14 total, 6 thin omitted)

### Community 0 - "PSPL Fitting Pipeline"
Cohesion: 0.06
Nodes (54): argparse, Han & Gould 2003, Kass & Raftery (BIC decisive threshold), Mass-distance-velocity degeneracy, Session 2 (2026-09-14), Folding single-caller scripts into mcmc_fit.py/mcmc_fit_binary.py, plot_fit() tag="" naming-collision gap, corner (+46 more)

### Community 1 - "Lens Model Physics"
Cohesion: 0.07
Nodes (39): astropy_coordinates, astropy_time, astropy_units, Witt & Mao (1995), O-03-BLG235 calibrate-once-then-fit pipeline, binary_images(), _binary_images_batch(), _companion_eigvals() (+31 more)

### Community 2 - "2L1S Trajectory & Parallax"
Cohesion: 0.10
Nodes (39): concurrent_futures, functools, binary_magnification(), binary_trajectory(), Total binary-lens magnification at source position(s) `zeta`. Sums 1/|J| over…, Returns (delta_s_n, delta_s_e): Earth's sky-projected position relative to the…, Source-lens separation for binary lenses on the lens plane, given in complex…, sun_earth_projection() (+31 more)

### Community 3 - "2L1S Validation & Comparison"
Cohesion: 0.09
Nodes (30): MulensModel (oracle library), Download raw photometry for each event we're working with. Sources are…, itertools, caustic_curve(), Caustic curve(s) in the source plane: the image, under the lens equation, of…, matplotlib_pyplot, mulensmodel, numpy (+22 more)

### Community 4 - "Session 10 Findings & Bugs"
Cohesion: 0.11
Nodes (25): alpha posterior bimodality finding, Bachelet et al. 2022 (arXiv:2205.07522), chi2-vs-Student-t MCMC comparison diagnostic, Wrong RA coordinate bug (O-03-BLG235), Frozen PSPL flux-calibration bias bug, torch/ProcessPoolExecutor fork-after-threading deadlock, /graphify skill build/update, KMT mag_err underestimation bug (+17 more)

### Community 5 - "Early Sessions & Datasets"
Cohesion: 0.09
Nodes (21): MulensModel alpha 180deg convention bug, Session 1 (2026-09-14), Session 3 (2026-09-14), sympy, u0 mirror-image sign bug pattern, Wyrzykowski et al. 2015, Bond et al. 2004 (arXiv:astro-ph/0404309), Ponytail-check rule (CLAUDE.md rule 4) (+13 more)

### Community 6 - "MCMC Likelihood Functions"
Cohesion: 0.17
Nodes (14): emcee, fit_parallax_pspl_mcmc(), log_likelihood(), log_prior(), log_probability(), fit_pspl_mcmc(), log_likelihood(), log_prior() (+6 more)

### Community 7 - "Data Loading & Calibration"
Cohesion: 0.18
Nodes (13): joint_fit.py (removed, folded into preprocess_binary_data.py), mag_to_flux(), Calibrated magnitude + error -> flux + error on the ZERO_POINT_MAG scale., load_coords(), load_raw(), ogle_to_magnification(), _parse_header_coord(), Reduce the raw OGLE + MOA files for OGLE-2003-BLG-235 / MOA-2003-BLG-53 into… (+5 more)

## Knowledge Gaps
- **24 isolated node(s):** `M-19-BLG008 (MOA-2019-BLG-008, KMT I-band subset)`, `emcee`, `matplotlib`, `numpy`, `scipy` (+19 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 110 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **6 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `sun_earth_projection()` connect `2L1S Trajectory & Parallax` to `Lens Model Physics`, `Data Loading & Calibration`?**
  _High betweenness centrality (0.075) - this node is a cross-community bridge._
- **Why does `binary_magnification()` connect `2L1S Trajectory & Parallax` to `Lens Model Physics`, `2L1S Validation & Comparison`, `Early Sessions & Datasets`?**
  _High betweenness centrality (0.074) - this node is a cross-community bridge._
- **Why does `Session 2 (2026-09-14)` connect `PSPL Fitting Pipeline` to `Early Sessions & Datasets`?**
  _High betweenness centrality (0.074) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `binary_magnification()` (e.g. with `binary_images()` and `chi2()`) actually correct?**
  _`binary_magnification()` has 2 INFERRED edges - model-reasoned connections that need verification._
- **Are the 3 inferred relationships involving `run_joint_fit()` (e.g. with `chi2()` and `log_probability()`) actually correct?**
  _`run_joint_fit()` has 3 INFERRED edges - model-reasoned connections that need verification._
- **Are the 3 inferred relationships involving `chi2()` (e.g. with `binary_magnification()` and `_fit_one_seed()`) actually correct?**
  _`chi2()` has 3 INFERRED edges - model-reasoned connections that need verification._
- **What connects `M-19-BLG008 (MOA-2019-BLG-008, KMT I-band subset)`, `emcee`, `matplotlib` to the rest of the system?**
  _24 weakly-connected nodes found - possible documentation gaps or missing edges._