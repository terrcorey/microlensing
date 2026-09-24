# Graph Report - microlensing  (2026-09-24)

## Corpus Check
- 4 files · ~28,006 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 298 nodes · 562 edges · 12 communities (9 shown, 3 thin omitted)
- Extraction: 96% EXTRACTED · 4% INFERRED · 0% AMBIGUOUS · INFERRED: 23 edges (avg confidence: 0.85)
- Token cost: 131,319 input · 0 output

## Community Hubs (Navigation)
- Lens Model Physics & Parallax
- 2L1S Fitting Track (chi2/Student-t/Huber)
- PSPL Fitting Pipeline
- Session History & Project Conventions
- Binary-Lens Model & MOA-2019-BLG-008
- CLAUDE.md 2L1S Architecture Notes
- PSPL Dataset & Requirements
- MulensModel Cross-Validation
- PSPL MCMC Likelihood Functions
- MOA-2019-BLG-008 short name (isolated)

## God Nodes (most connected - your core abstractions)
1. `TwoL1SParams` - 15 edges
2. `plot_fit()` - 12 edges
3. `run_joint_fit()` - 11 edges
4. `binary_magnification()` - 11 edges
5. `chi2()` - 11 edges
6. `fit_parallax_pspl_mcmc()` - 10 edges
7. `residuals()` - 10 edges
8. `Annual parallax (Gould 2004 geocentric formalism)` - 10 edges
9. `run_ogle_only_diagnostic()` - 9 edges
10. `plot_fit()` - 9 edges

## Surprising Connections (you probably didn't know these)
- `estimate_mass()` --implements--> `Prior-based lens mass estimate (M_lens)`  [EXTRACTED]
  mcmc_fit.py → README.md
- `Wrong-but-plausible convergence pattern` --semantically_similar_to--> `Near-miss-cusp search-landscape problem`  [INFERRED] [semantically similar]
  CHANGELOG.md → CLAUDE.md
- `chi2()` --calls--> `binary_magnification()`  [INFERRED]
  scratch/fit_2l1s.py → lc_models.py
- `run_ogle_only_diagnostic()` --calls--> `fit_parallax_pspl_mcmc()`  [EXTRACTED]
  mcmc_fit_binary.py → mcmc_fit.py
- `run_ogle_only_diagnostic()` --calls--> `fit_pspl_mcmc()`  [EXTRACTED]
  mcmc_fit_binary.py → mcmc_fit.py

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Annual Parallax + Student-t Rollout Across Events** — claude_annual_parallax, claude_studentt_likelihood, claude_o05_blg086, claude_o03_blg235 [INFERRED 0.85]
- **Recurring Wrong-But-Plausible Convergence Pattern** — changelog_session_1, changelog_session_2, changelog_session_3, changelog_session_4, changelog_wrong_but_plausible_pattern [EXTRACTED 1.00]
- **TwoL1SParams Canonical-Order Refactor Group** — claude_fit_2l1s_twol1sparams, claude_fit_2l1s_residuals, claude_fit_2l1s_plot_fit, claude_mcmc_fit_2l1s_run_mcmc [EXTRACTED 1.00]

## Communities (12 total, 3 thin omitted)

### Community 0 - "Lens Model Physics & Parallax"
Cohesion: 0.06
Nodes (51): astropy_coordinates, astropy_time, astropy_units, binary_images(), _binary_images_batch(), _companion_eigvals(), flux(), lens_position() (+43 more)

### Community 1 - "2L1S Fitting Track (chi2/Student-t/Huber)"
Cohesion: 0.08
Nodes (50): functools, mpl_toolkits_axes_grid1_inset_locator, multiprocessing, os, chi2(), _fit_one_seed(), huber(), plot_fit() (+42 more)

### Community 2 - "PSPL Fitting Pipeline"
Cohesion: 0.08
Nodes (40): argparse, corner, Download raw photometry for each event we're working with. Sources are…, emcee, matplotlib_pyplot, plot_raw(), MCMC PSPL fits for OGLE-2003-BLG-235 / MOA-2003-BLG-53. Both fits below use the…, Combined OGLE+MOA calibrated magnification, no fit overlay -- the --stage=raw… (+32 more)

### Community 3 - "Session History & Project Conventions"
Cohesion: 0.07
Nodes (41): Bachelet et al. 2022 (arXiv:2205.07522), Bond et al. 2004 (arXiv:astro-ph/0404309), mcmc_fit.estimate_mass(), Han & Gould 2003 (disk/bulge-weighted Galactic prior), Kass & Raftery (BIC decisive-threshold convention), Mass-distance-velocity degeneracy, MOA-2019-BLG-008, Session 1 — Two Real Events Set Up End-to-End (+33 more)

### Community 4 - "Binary-Lens Model & MOA-2019-BLG-008"
Cohesion: 0.14
Nodes (23): concurrent_futures, binary_magnification(), binary_trajectory(), caustic_curve(), Total binary-lens magnification at source position(s) `zeta`. Sums 1/|J| over…, Source-lens separation for binary lenses on the lens plane, given in complex…, Caustic curve(s) in the source plane: the image, under the lens equation, of…, moa_to_magnification() (+15 more)

### Community 5 - "CLAUDE.md 2L1S Architecture Notes"
Cohesion: 0.12
Nodes (16): scratch/cross_check_mulensmodel*.py scripts, fit_2l1s.TwoL1SParams NamedTuple, lc_models.binary_magnification(), mcmc_fit_2l1s.run_mcmc(), mcmc_fit_binary.run_ogle_only_diagnostic(), mcmc_fit.fit_parallax_pspl_mcmc(), mcmc_fit.FitResult NamedTuple, mcmc_fit.get_t0_par() (+8 more)

### Community 6 - "PSPL Dataset & Requirements"
Cohesion: 0.11
Nodes (13): OGLE-2005-BLG-086 (single-lens event), PSPL (Paczynski point-source point-lens) model, O-03-BLG235 (OGLE-2003-BLG-235 / MOA-2003-BLG-53), O-05-BLG086 (OGLE-2005-BLG-086), Prior-based lens mass estimate (M_lens), Mass-distance-velocity degeneracy (central open question), astropy, corner (+5 more)

### Community 7 - "MulensModel Cross-Validation"
Cohesion: 0.16
Nodes (15): itertools, mulensmodel, chi2(), mm_magnification(), plot_fit(), One-time cross-check: fit MulensModel's own binary-lens model to the same…, Chi2 of MulensModel's binary magnification against the combined OGLE+MOA data., Serial multi-start search, then a tight-tolerance refit of the best seed. (+7 more)

### Community 8 - "PSPL MCMC Likelihood Functions"
Cohesion: 0.18
Nodes (13): fit_parallax_pspl_mcmc(), log_likelihood(), log_prior(), log_probability(), fit_pspl_mcmc(), log_likelihood(), log_prior(), log_probability() (+5 more)

## Knowledge Gaps
- **18 isolated node(s):** `M-19-BLG008 (MOA-2019-BLG-008, KMT I-band subset)`, `astropy`, `O-03-BLG235 (OGLE-2003-BLG-235 / MOA-2003-BLG-53)`, `O-05-BLG086 (OGLE-2005-BLG-086)`, `emcee` (+13 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 118 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **3 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `lc_models.binary_trajectory()` connect `Session History & Project Conventions` to `Binary-Lens Model & MOA-2019-BLG-008`?**
  _High betweenness centrality (0.122) - this node is a cross-community bridge._
- **Why does `MOA-2019-BLG-008` connect `Session History & Project Conventions` to `Binary-Lens Model & MOA-2019-BLG-008`?**
  _High betweenness centrality (0.119) - this node is a cross-community bridge._
- **Why does `Output layout convention (scratch/ vs. pipeline dirs)` connect `CLAUDE.md 2L1S Architecture Notes` to `Binary-Lens Model & MOA-2019-BLG-008`?**
  _High betweenness centrality (0.111) - this node is a cross-community bridge._
- **Are the 3 inferred relationships involving `run_joint_fit()` (e.g. with `chi2()` and `log_probability()`) actually correct?**
  _`run_joint_fit()` has 3 INFERRED edges - model-reasoned connections that need verification._
- **What connects `M-19-BLG008 (MOA-2019-BLG-008, KMT I-band subset)`, `astropy`, `O-03-BLG235 (OGLE-2003-BLG-235 / MOA-2003-BLG-53)` to the rest of the system?**
  _18 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Lens Model Physics & Parallax` be split into smaller, more focused modules?**
  _Cohesion score 0.055218855218855216 - nodes in this community are weakly interconnected._
- **Should `2L1S Fitting Track (chi2/Student-t/Huber)` be split into smaller, more focused modules?**
  _Cohesion score 0.08446455505279035 - nodes in this community are weakly interconnected._