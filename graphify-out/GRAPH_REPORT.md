# Graph Report - microlensing  (2026-09-22)

## Corpus Check
- Corpus is ~19,624 words - fits in a single context window. You may not need a graph.

## Summary
- 215 nodes · 452 edges · 9 communities (8 shown, 1 thin omitted)
- Extraction: 95% EXTRACTED · 5% INFERRED · 1% AMBIGUOUS · INFERRED: 21 edges (avg confidence: 0.85)
- Token cost: 112,288 input · 0 output

## Community Hubs (Navigation)
- 2L1S Model & Calibration Preprocessing
- MCMC Fitting Pipeline & Setup
- Project Documentation & Citations
- PSPL Flux/Magnitude & Parallax Fit
- lc_models Low-Level Geometry & Astropy
- MulensModel Cross-Check & Known Bugs
- PSPL Core & Joint-Fit Calibration
- MOA-2019-BLG-008 Fit & Error-Bar Bug
- Interactive Dev Notebook

## God Nodes (most connected - your core abstractions)
1. `binary_magnification()` - 15 edges
2. `fit_pspl_mcmc()` - 14 edges
3. `plot_fit()` - 14 edges
4. `magnification()` - 13 edges
5. `binary_trajectory()` - 12 edges
6. `run_joint_fit()` - 12 edges
7. `chi2()` - 12 edges
8. `fit_parallax_pspl_mcmc()` - 10 edges
9. `run_ogle_only_diagnostic()` - 10 edges
10. `find_zoom_window()` - 10 edges

## Surprising Connections (you probably didn't know these)
- `binary_magnification()` --references--> `Witt & Mao 1995`  [EXTRACTED]
  lc_models.py → CHANGELOG.md
- `u0 mirror-image sign bug` --rationale_for--> `fit_pspl_mcmc()`  [EXTRACTED]
  CHANGELOG.md → mcmc_fit.py
- `fit_pspl_mcmc()` --shares_data_with--> `emcee`  [EXTRACTED]
  mcmc_fit.py → requirements.txt
- `O-03-BLG235 frozen PSPL flux-calibration bug` --rationale_for--> `profile_flux()`  [EXTRACTED]
  CHANGELOG.md → scratch/fit_2l1s.py
- `Non-physical near-miss cusp solutions` --rationale_for--> `caustic_curve()`  [EXTRACTED]
  CHANGELOG.md → lc_models.py

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Chi2/posterior shape can't detect model misspecification (recurring pattern)** — changelog_u0_sign_bug, changelog_mulensmodel_alpha_convention_bug, changelog_near_miss_cusp_problem, mcmc_fit_binary_run_ogle_only_diagnostic [INFERRED 0.85]
- **O-03-BLG235 flux-calibration pipeline (frozen vs. per-trial profiling)** — preprocess_binary_data_fit_joint_pspl, mcmc_fit_binary_run_joint_fit, scratch_fit_2l1s_profile_flux, changelog_calibration_bug [INFERRED 0.85]
- **2L1S implementation validation trio (symbolic derivation, internal checks, oracle cross-check)** — lc_models_binary_magnification, scratch_cross_check_mulensmodel, scratch_derive_binary_quintic, changelog_witt_mao_1995 [INFERRED 0.85]

## Communities (9 total, 1 thin omitted)

### Community 0 - "2L1S Model & Calibration Preprocessing"
Cohesion: 0.09
Nodes (43): functools, binary_magnification(), binary_trajectory(), caustic_curve(), Total binary-lens magnification at source position(s) `zeta`. Sums 1/|J| over…, Source-lens separation for binary lenses on the lens plane, given in complex…, Caustic curve(s) in the source plane: the image, under the lens equation, of…, multiprocessing (+35 more)

### Community 1 - "MCMC Fitting Pipeline & Setup"
Cohesion: 0.08
Nodes (33): argparse, Han & Gould 2003, corner, Download raw photometry for each event we're working with. Sources are…, emcee, matplotlib_pyplot, plot_raw(), MCMC PSPL fits for OGLE-2003-BLG-235 / MOA-2003-BLG-53. Both fits below use the… (+25 more)

### Community 2 - "Project Documentation & Citations"
Cohesion: 0.09
Nodes (25): Unresolved ~15deg alpha discrepancy vs Bond et al. 2004, Bachelet et al. 2022 (arXiv:2205.07522), BIC-based PSPL-vs-2L1S model comparison verdict, Bond et al. 2004 (arXiv:astro-ph/0404309), Gould 2004 (geocentric parallax formalism), Kass & Raftery (BIC/Bayes-factor interpretation), Non-physical near-miss cusp solutions, Witt & Mao 1995 (+17 more)

### Community 3 - "PSPL Flux/Magnitude & Parallax Fit"
Cohesion: 0.13
Nodes (22): flux(), magnitude(), plain_flux(), plain_magnitude(), Observed flux: magnified source plus a constant blend flux., Observed magnitude, for overplotting against real photometry., flux() with piE_N=piE_E=0 -- the exact pre-parallax model, used only to…, magnitude() with piE_N=piE_E=0 -- see plain_flux(). (+14 more)

### Community 4 - "lc_models Low-Level Geometry & Astropy"
Cohesion: 0.13
Nodes (19): astropy_coordinates, astropy_time, astropy_units, _binary_images_batch(), _companion_eigvals(), lens_position(), _project(), _quintic_coefficients() (+11 more)

### Community 5 - "MulensModel Cross-Check & Known Bugs"
Cohesion: 0.12
Nodes (19): MulensModel package, MulensModel alpha-convention cross-check bug, u0 mirror-image sign bug, itertools, mulensmodel, chi2(), mm_magnification(), plot_fit() (+11 more)

### Community 6 - "PSPL Core & Joint-Fit Calibration"
Cohesion: 0.15
Nodes (19): O-03-BLG235 frozen PSPL flux-calibration bug, mag_to_flux(), magnification(), Calibrated magnitude + error -> flux + error on the ZERO_POINT_MAG scale., Source-lens separation u(t), in units of the Einstein radius., Paczynski point-lens magnification A(u)., trajectory(), run_joint_fit() (+11 more)

### Community 7 - "MOA-2019-BLG-008 Fit & Error-Bar Bug"
Cohesion: 0.21
Nodes (12): MOA-2019-BLG-008 KMT error-bar underestimation, Annual parallax + Student-t likelihood plan, concurrent_futures, chi2(), _fit_one_seed(), Multi-start Nelder-Mead 2L1S fit for MOA-2019-BLG-008 (Bachelet et al. 2022), a…, Closed-form weighted least squares for (fs, fb) at fixed magnification shape --…, Chi2 of the 2L1S model against KMT I-band flux, profiling fs/fb at every trial. (+4 more)

## Ambiguous Edges - Review These
- `derive_binary_quintic.py` → `sympy`  [AMBIGUOUS]
  requirements.txt · relation: conceptually_related_to
- `CLAUDE.md` → `networkx`  [AMBIGUOUS]
  requirements.txt · relation: conceptually_related_to
- `README.md` → `requirements.txt`  [AMBIGUOUS]
  README.md · relation: conceptually_related_to

## Knowledge Gaps
- **4 isolated node(s):** `M-19-BLG008 (MOA-2019-BLG-008, KMT I-band subset)`, `Point-Source Binary-Lens (2L1S) model`, `BIC-based PSPL-vs-2L1S model comparison verdict`, `scratch/scratchpad.ipynb`
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 84 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **1 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **What is the exact relationship between `derive_binary_quintic.py` and `sympy`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._
- **What is the exact relationship between `CLAUDE.md` and `networkx`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._
- **What is the exact relationship between `README.md` and `requirements.txt`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._
- **Why does `sun_earth_projection()` connect `Project Documentation & Citations` to `MCMC Fitting Pipeline & Setup`, `lc_models Low-Level Geometry & Astropy`, `MOA-2019-BLG-008 Fit & Error-Bar Bug`?**
  _High betweenness centrality (0.076) - this node is a cross-community bridge._
- **Why does `binary_magnification()` connect `2L1S Model & Calibration Preprocessing` to `Project Documentation & Citations`, `lc_models Low-Level Geometry & Astropy`, `MulensModel Cross-Check & Known Bugs`, `MOA-2019-BLG-008 Fit & Error-Bar Bug`?**
  _High betweenness centrality (0.073) - this node is a cross-community bridge._
- **Why does `fit_pspl_mcmc()` connect `PSPL Flux/Magnitude & Parallax Fit` to `MCMC Fitting Pipeline & Setup`, `Project Documentation & Citations`, `MulensModel Cross-Check & Known Bugs`, `MOA-2019-BLG-008 Fit & Error-Bar Bug`?**
  _High betweenness centrality (0.061) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `fit_pspl_mcmc()` (e.g. with `plain_magnitude()` and `log_probability()`) actually correct?**
  _`fit_pspl_mcmc()` has 2 INFERRED edges - model-reasoned connections that need verification._