# Graph Report - microlensing  (2026-09-28)

## Corpus Check
- 8 files · ~34,877 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 308 nodes · 640 edges · 10 communities (7 shown, 3 thin omitted)
- Extraction: 95% EXTRACTED · 5% INFERRED · 0% AMBIGUOUS · INFERRED: 29 edges (avg confidence: 0.83)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- Roadmap, Papers & Session Log
- Lens Model Core (PSPL/2L1S)
- PSPL Pipeline & MCMC
- 2L1S Fit & Flux Profiling
- Cross-Checks & Data Download
- Cassan Caustic Parametrisation
- Finite-Source Magnification
- Quintic Derivation (sympy)
- M-19-BLG008 short name (isolated)

## God Nodes (most connected - your core abstractions)
1. `CHANGELOG (session log)` - 27 edges
2. `CLAUDE.md (project guide)` - 26 edges
3. `TwoL1SParams` - 18 edges
4. `plot_fit()` - 18 edges
5. `binary_magnification()` - 14 edges
6. `chi2()` - 14 edges
7. `binary_trajectory()` - 12 edges
8. `caustic_curve()` - 12 edges
9. `binary_magnification_fs()` - 12 edges
10. `residuals()` - 12 edges

## Surprising Connections (you probably didn't know these)
- `estimate_mass()` --implements--> `Prior-based lens mass estimate (M_lens)`  [EXTRACTED]
  mcmc_fit.py → README.md
- `u0 mirror-sign degeneracy bug (per-dataset initial guesses)` --semantically_similar_to--> `MulensModel oracle cross-checks`  [INFERRED] [semantically similar]
  CHANGELOG.md → CLAUDE.md
- `Error-bar rescaling / additive jitter` --semantically_similar_to--> `Student-t robust likelihood (free scale/dof)`  [INFERRED] [semantically similar]
  CHANGELOG.md → CLAUDE.md
- `SIGNALMEN anomaly-vs-outlier test (Dominik et al. 2007)` --semantically_similar_to--> `Huber loss likelihood (DELTA=1.345)`  [INFERRED] [semantically similar]
  CHANGELOG.md → CLAUDE.md
- `chi2()` --calls--> `binary_magnification()`  [INFERRED]
  scratch/fit_2l1s.py → lc_models.py

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **2L1S model validation via limiting cases + MulensModel oracle** — claude_2l1s_model, claude_finite_source, claude_annual_parallax, claude_mulensmodel_oracle [EXTRACTED 1.00]
- **2L1S roadmap: finite source, error bars, Cassan, (d,q) grid, GA, SIGNALMEN** — claude_finite_source, changelog_error_bar_rescaling, claude_cassan_parametrisation, changelog_dq_grid_search, changelog_genetic_algorithm, changelog_signalmen [EXTRACTED 1.00]
- **Likelihood choices explored (chi2 / Student-t / Huber)** — claude_student_t_likelihood, claude_huber_likelihood, changelog_error_bar_rescaling, claude_mcmc_emcee [INFERRED 0.85]

## Communities (10 total, 3 thin omitted)

### Community 0 - "Roadmap, Papers & Session Log"
Cohesion: 0.06
Nodes (59): CHANGELOG (session log), Bachelet et al. 2022, BIC PSPL-vs-2L1S verdict (deltaBIC>10), Bond et al. 2004 (O-03-BLG235 planet solution), Cassan 2008, Charbonneau 1995, Dominik et al. 2007 (MNRAS 380, 792), (d,q) grid search (+51 more)

### Community 1 - "Lens Model Core (PSPL/2L1S)"
Cohesion: 0.05
Nodes (52): astropy_coordinates, astropy_time, astropy_units, binary_images(), _binary_images_batch(), _companion_eigvals(), flux(), lens_position() (+44 more)

### Community 2 - "PSPL Pipeline & MCMC"
Cohesion: 0.06
Nodes (48): argparse, corner, emcee, plot_raw(), MCMC PSPL fits for OGLE-2003-BLG-235 / MOA-2003-BLG-53. Both fits below use the…, Combined OGLE+MOA calibrated magnification, no fit overlay -- the --stage=raw…, run_joint_fit(), log_likelihood() (+40 more)

### Community 3 - "2L1S Fit & Flux Profiling"
Cohesion: 0.09
Nodes (45): os, chi2(), _fit_one_seed(), huber(), plot_fit(), profile_flux(), _profile_fs_moa(), NamedTuple (+37 more)

### Community 4 - "Cross-Checks & Data Download"
Cohesion: 0.10
Nodes (25): concurrent_futures, Download raw photometry for each event we're working with. Sources are…, itertools, matplotlib_pyplot, mulensmodel, numpy, pathlib, Compare every 2L1S candidate solution found this session (Nelder-Mead multi-… (+17 more)

### Community 5 - "Cassan Caustic Parametrisation"
Cohesion: 0.12
Nodes (27): functools, cassan_caustic(), cassan_to_standard(), caustic_curve(), caustic_point(), caustic_curve(s, q)[idx], resampled evenly in arc length, starting at its…, Point at abscissa sigma (wraps with period 1) on a cassan_caustic() output., Cassan (2008) caustic-crossing parameters -> binary_trajectory()'s (t0, u0, tE,… (+19 more)

### Community 6 - "Finite-Source Magnification"
Cohesion: 0.12
Nodes (24): binary_magnification(), binary_magnification_fs(), _disk_average(), _sample_coords(), binary_trajectory(), equidistant_caustic(), Resample one closed caustic (an element of caustic_curve()'s list) to points…, Total binary-lens magnification at source position(s) `zeta`. Sums 1/|J| over… (+16 more)

## Knowledge Gaps
- **6 isolated node(s):** `M-19-BLG008 (MOA-2019-BLG-008, KMT I-band subset)`, `Time = HJD - 2450000 convention`, `scipy`, `numpy`, `matplotlib` (+1 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 117 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **3 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `estimate_mass()` connect `PSPL Pipeline & MCMC` to `Roadmap, Papers & Session Log`?**
  _High betweenness centrality (0.316) - this node is a cross-community bridge._
- **Why does `Prior-based lens mass estimate (M_lens)` connect `Roadmap, Papers & Session Log` to `PSPL Pipeline & MCMC`?**
  _High betweenness centrality (0.305) - this node is a cross-community bridge._
- **What connects `M-19-BLG008 (MOA-2019-BLG-008, KMT I-band subset)`, `Time = HJD - 2450000 convention`, `scipy` to the rest of the system?**
  _6 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Roadmap, Papers & Session Log` be split into smaller, more focused modules?**
  _Cohesion score 0.059887005649717516 - nodes in this community are weakly interconnected._
- **Should `Lens Model Core (PSPL/2L1S)` be split into smaller, more focused modules?**
  _Cohesion score 0.0512987012987013 - nodes in this community are weakly interconnected._
- **Should `PSPL Pipeline & MCMC` be split into smaller, more focused modules?**
  _Cohesion score 0.0632996632996633 - nodes in this community are weakly interconnected._
- **Should `2L1S Fit & Flux Profiling` be split into smaller, more focused modules?**
  _Cohesion score 0.09082125603864734 - nodes in this community are weakly interconnected._