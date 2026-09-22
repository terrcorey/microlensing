# Graph Report - microlensing  (2026-09-22)

## Corpus Check
- 8 files · ~23,240 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 245 nodes · 479 edges · 13 communities (9 shown, 4 thin omitted)
- Extraction: 96% EXTRACTED · 4% INFERRED · 0% AMBIGUOUS · INFERRED: 19 edges (avg confidence: 0.87)
- Token cost: 134,651 input · 0 output

## Community Hubs (Navigation)
- O-03-BLG235 Binary Pipeline
- 2L1S Model + Torch Solver
- Magnification Calibration + Multi-start
- CLAUDE.md / CHANGELOG Narrative
- 2L1S Joint-Fit Search
- PSPL Model + Parallax Roadmap
- Parallax + Student-t MCMC
- 2L1S Search Pitfalls + Model Comparison
- MulensModel Cross-Check
- Binary Quintic Derivation (sympy)
- Output Layout Convention
- MOA-2019-BLG-008 Entry

## God Nodes (most connected - your core abstractions)
1. `binary_magnification()` - 15 edges
2. `plot_fit()` - 14 edges
3. `run_joint_fit()` - 14 edges
4. `PSPL (point-source point-lens) model` - 14 edges
5. `binary_trajectory()` - 13 edges
6. `2L1S (point-source binary-lens) model` - 13 edges
7. `fit_parallax_pspl_mcmc()` - 12 edges
8. `chi2()` - 11 edges
9. `find_zoom_window()` - 11 edges
10. `run_ogle_only_diagnostic()` - 10 edges

## Surprising Connections (you probably didn't know these)
- `binary_magnification()` --implements--> `2L1S (point-source binary-lens) model`  [EXTRACTED]
  lc_models.py → CLAUDE.md
- `binary_trajectory()` --implements--> `2L1S (point-source binary-lens) model`  [EXTRACTED]
  lc_models.py → CLAUDE.md
- `binary_images()` --implements--> `2L1S (point-source binary-lens) model`  [EXTRACTED]
  lc_models.py → CLAUDE.md
- `sun_earth_projection()` --implements--> `Annual parallax (Gould 2004 geocentric formalism)`  [EXTRACTED]
  lc_models.py → CLAUDE.md
- `plain_flux()` --implements--> `"Zero out new params to recover old special case" reuse trick`  [EXTRACTED]
  lc_models.py → CHANGELOG.md

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Recurring wrong-but-plausible-result failure pattern** — changelog_wrong_but_plausible_pattern, changelog_u0_sign_bug, claude_near_miss_cusp_problem, changelog_mulensmodel_convention_bug, changelog_upstream_chi2_pattern [EXTRACTED 0.90]
- **Annual parallax implementation roadmap** — claude_annual_parallax, lc_models_sun_earth_projection, mcmc_fit_fit_parallax_pspl_mcmc, preprocess_binary_data_fit_joint_pspl, claude_gould_2004 [EXTRACTED 0.90]
- **Shared plotting/MCMC helper functions** — mcmc_fit_save_corner, mcmc_fit_plot_fit_lc, mcmc_fit_fitresult, mcmc_fit_get_t0_par, zoom_utils_plot_fit_panels, zoom_utils_plot_residual_panel, zoom_utils_plot_residual_hist, claude_shared_helpers_note [EXTRACTED 0.90]

## Communities (13 total, 4 thin omitted)

### Community 0 - "O-03-BLG235 Binary Pipeline"
Cohesion: 0.08
Nodes (40): argparse, FitResult NamedTuple pattern, get_t0_par() cache-bootstrap pattern, Shared plotting/MCMC helpers ("do not re-copy"), corner, emcee, matplotlib_pyplot, plot_raw() (+32 more)

### Community 1 - "2L1S Model + Torch Solver"
Cohesion: 0.08
Nodes (33): astropy_coordinates, astropy_time, astropy_units, Download raw photometry for each event we're working with. Sources are…, binary_images(), _binary_images_batch(), _companion_eigvals(), lens_position() (+25 more)

### Community 2 - "Magnification Calibration + Multi-start"
Cohesion: 0.12
Nodes (28): concurrent_futures, binary_magnification(), binary_trajectory(), caustic_curve(), Total binary-lens magnification at source position(s) `zeta`. Sums 1/|J| over…, Source-lens separation for binary lenses on the lens plane, given in complex…, Caustic curve(s) in the source plane: the image, under the lens equation, of…, moa_to_magnification() (+20 more)

### Community 3 - "CLAUDE.md / CHANGELOG Narrative"
Cohesion: 0.10
Nodes (25): Bachelet et al. 2022 (arXiv:2205.07522), Bond et al. 2004 (arXiv:astro-ph/0404309), "User writes, Claude reviews" collaboration mode, graphify-out/ relocation into microlensing/, MOA-2019-BLG-008, networkx (erroneously listed, unused, dropped), "Bad chi2 can come from upstream issues, not the model" pattern, Wyrzykowski et al. 2015 (+17 more)

### Community 4 - "2L1S Joint-Fit Search"
Cohesion: 0.11
Nodes (26): functools, mag_to_flux(), Calibrated magnitude + error -> flux + error on the ZERO_POINT_MAG scale., multiprocessing, load_raw(), Raw OGLE (mag) + MOA (differential flux) tables, time-shifted to HJD-2450000., chi2(), _fit_one_seed() (+18 more)

### Community 5 - "PSPL Model + Parallax Roadmap"
Cohesion: 0.10
Nodes (25): "Zero out new params to recover old special case" reuse trick, Han & Gould (2003) Galactic density prior, Central open question restated (session 1): tE can't separate mass from distance, Annual parallax (Gould 2004 geocentric formalism), Gould (2004) parallax formalism, PSPL (point-source point-lens) model, flux(), magnification() (+17 more)

### Community 6 - "Parallax + Student-t MCMC"
Cohesion: 0.17
Nodes (14): Student-t robust likelihood, fit_parallax_pspl_mcmc(), log_likelihood(), log_prior(), log_probability(), fit_pspl_mcmc(), log_likelihood(), log_prior() (+6 more)

### Community 7 - "2L1S Search Pitfalls + Model Comparison"
Cohesion: 0.18
Nodes (13): BIC-based 2L1S-vs-PSPL model-comparison verdict plan, Kass & Raftery (BIC decisive threshold), mpmath (transitive sympy dependency, dropped), MulensModel cross-check convention-bug pattern (alpha 180deg offset), u0-sign / mirror-image degeneracy bug pattern, Witt & Mao (1995), "Wrong-but-plausible result" recurring failure pattern, 2L1S (point-source binary-lens) model (+5 more)

### Community 8 - "MulensModel Cross-Check"
Cohesion: 0.25
Nodes (10): itertools, scipy_optimize, chi2(), mm_magnification(), plot_fit(), One-time cross-check: fit MulensModel's own binary-lens model to the same…, Chi2 of MulensModel's binary magnification against the combined OGLE+MOA data., Serial multi-start search, then a tight-tolerance refit of the best seed. (+2 more)

## Ambiguous Edges - Review These
- `/graphify skill` → `networkx (erroneously listed, unused, dropped)`  [AMBIGUOUS]
  CHANGELOG.md · relation: conceptually_related_to

## Knowledge Gaps
- **12 isolated node(s):** `M-19-BLG008 (MOA-2019-BLG-008, KMT I-band subset)`, `ponytail skill`, `Wyrzykowski et al. 2015`, `Bond et al. 2004 (arXiv:astro-ph/0404309)`, `Bachelet et al. 2022 (arXiv:2205.07522)` (+7 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 102 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **4 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **What is the exact relationship between `/graphify skill` and `networkx (erroneously listed, unused, dropped)`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._
- **Why does `2L1S (point-source binary-lens) model` connect `2L1S Search Pitfalls + Model Comparison` to `2L1S Model + Torch Solver`, `Magnification Calibration + Multi-start`, `CLAUDE.md / CHANGELOG Narrative`?**
  _High betweenness centrality (0.120) - this node is a cross-community bridge._
- **Why does `run_joint_fit()` connect `O-03-BLG235 Binary Pipeline` to `CLAUDE.md / CHANGELOG Narrative`, `PSPL Model + Parallax Roadmap`?**
  _High betweenness centrality (0.112) - this node is a cross-community bridge._
- **Why does `PSPL (point-source point-lens) model` connect `PSPL Model + Parallax Roadmap` to `O-03-BLG235 Binary Pipeline`, `CLAUDE.md / CHANGELOG Narrative`, `Parallax + Student-t MCMC`, `2L1S Search Pitfalls + Model Comparison`?**
  _High betweenness centrality (0.098) - this node is a cross-community bridge._
- **Are the 3 inferred relationships involving `run_joint_fit()` (e.g. with `chi2()` and `log_probability()`) actually correct?**
  _`run_joint_fit()` has 3 INFERRED edges - model-reasoned connections that need verification._
- **What connects `M-19-BLG008 (MOA-2019-BLG-008, KMT I-band subset)`, `ponytail skill`, `Wyrzykowski et al. 2015` to the rest of the system?**
  _12 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `O-03-BLG235 Binary Pipeline` be split into smaller, more focused modules?**
  _Cohesion score 0.08309178743961353 - nodes in this community are weakly interconnected._