# Graph Report - microlensing  (2026-10-06)

## Corpus Check
- 26 files · ~46,401 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 436 nodes · 1004 edges · 15 communities (10 shown, 5 thin omitted)
- Extraction: 95% EXTRACTED · 5% INFERRED · 0% AMBIGUOUS · INFERRED: 48 edges (avg confidence: 0.89)
- Token cost: 87,484 input · 0 output

## Community Hubs (Navigation)
- Event Loading & Flux Calibration
- Cassan & Binary Preprocessing
- Search Pipeline & Diagnostics
- Search Pipeline & Diagnostics
- Binary-Lens Magnification
- Binary-Lens Magnification
- Event Loading & Flux Calibration
- Parallax PSPL MCMC
- Community 8
- Community 9
- Data, Parallax & PSPL Model
- Quintic Derivation
- SLURM Submission
- Community 13
- Community 14

## God Nodes (most connected - your core abstractions)
1. `plot_fit()` - 21 edges
2. `TwoL1SParams` - 19 edges
3. `magnification()` - 19 edges
4. `binary_trajectory()` - 16 edges
5. `binary_magnification_fs()` - 15 edges
6. `trajectory()` - 15 edges
7. `run_joint_fit()` - 15 edges
8. `caustic_curve()` - 14 edges
9. `residuals()` - 14 edges
10. `flat_chain()` - 14 edges

## Surprising Connections (you probably didn't know these)
- `TwoL1SParams NamedTuple (canonical 2L1S param order)` --references--> `TwoL1SParams`  [INFERRED]
  CLAUDE.md → scratch/fit_2l1s.py
- `mcmc_fit.FitResult NamedTuple (best_fit, samples, labels)` --references--> `FitResult`  [INFERRED]
  CLAUDE.md → mcmc_fit.py
- `Error-bar rescaling (per-instrument K, chi2/dof=1)` --references--> `error_scale()`  [INFERRED]
  CLAUDE.md → event.py
- `Per-trial profiled flux calibration (weighted lstsq fs/fb)` --references--> `flux_residuals()`  [INFERRED]
  CLAUDE.md → event.py
- `Per-event TOML config (input/<short>.toml)` --references--> `load_event()`  [INFERRED]
  CLAUDE.md → event.py

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Config-driven blind search flow** — claude_toml_event_config, claude_blind_search, claude_point_source_grid, claude_cassan_parametrisation, claude_bic_model_comparison, changelog_autocorrelation_check [EXTRACTED 1.00]
- **Likelihood options explored (Student-t, Huber, rescaled Gaussian chi2)** — claude_student_t_likelihood, claude_huber_likelihood, claude_error_bar_rescaling [INFERRED 0.85]
- **2L1S magnification implementations cross-validated** — claude_vbbinarylensing, claude_finite_source_fs, changelog_torch_batched_eigvals, claude_mulensmodel_cross_check, changelog_sympy_quintic [INFERRED 0.85]

## Communities (15 total, 5 thin omitted)

### Community 0 - "Event Loading & Flux Calibration"
Cohesion: 0.05
Nodes (82): argparse, concurrent_futures, corner, functools, mag_to_flux(), magnification(), magnitude(), Calibrated magnitude + error -> flux + error on the ZERO_POINT_MAG scale. (+74 more)

### Community 1 - "Cassan & Binary Preprocessing"
Cohesion: 0.06
Nodes (69): emcee, EnsembleSampler, cassan_to_standard(), Cassan (2008) caustic-crossing parameters -> binary_trajectory()'s (t0, u0, tE,…, flat_chain(), ndarray, Flat (samples, log_probs) after burn-in/thinning -- shared by every emcee…, Corner plot of the raw MCMC parameters -- shared by every fit script.… (+61 more)

### Community 2 - "Search Pipeline & Diagnostics"
Cohesion: 0.06
Nodes (61): Autocorrelation (tau) convergence check, nsteps/tau > 50, Bennett et al. 2015 re-reduction of OGLE-2005-BLG-169, Concurrent MCMC modes on one shared Pool, Prior-based lens mass estimate (estimate_mass, session 2), find_zoom_window significance cut (5 sigma), Offset-scaled archive magnitudes (FTN/MDM/Auckland fb<0), Geocentric JD light-travel-time correction, Gould et al. 2006 (OGLE-2005-BLG-169 data, UID 0300030) (+53 more)

### Community 3 - "Search Pipeline & Diagnostics"
Cohesion: 0.07
Nodes (57): distinct_modes() tolerance-based mode dedup, Event, NamedTuple, ndarray, Path, scipy_ndimage, anomaly_times(), axis() (+49 more)

### Community 4 - "Binary-Lens Magnification"
Cohesion: 0.06
Nodes (34): astropy_coordinates, astropy_time, astropy_units, binary_images(), _binary_images_batch(), binary_magnification_vbbl(), caustic_point(), _companion_eigvals() (+26 more)

### Community 5 - "Binary-Lens Magnification"
Cohesion: 0.08
Nodes (38): complex128, binary_magnification(), binary_magnification_fs(), _disk_average(), _sample_coords(), binary_trajectory(), cassan_caustic(), caustic_curve() (+30 more)

### Community 6 - "Event Loading & Flux Calibration"
Cohesion: 0.14
Nodes (25): chi2(), error_scale(), Event, flux_residuals(), Instrument, load_event(), load_instrument(), profile_flux() (+17 more)

### Community 7 - "Parallax PSPL MCMC"
Cohesion: 0.14
Nodes (20): flux(), plain_flux(), plain_magnitude(), Observed flux: magnified source plus a constant blend flux., flux() with piE_N=piE_E=0 -- the exact pre-parallax model, used only to…, magnitude() with piE_N=piE_E=0 -- see plain_flux()., fit_parallax_pspl_mcmc(), log_likelihood() (+12 more)

### Community 8 - "Community 8"
Cohesion: 0.29
Nodes (9): itertools, chi2(), mm_magnification(), plot_fit(), One-time cross-check: fit MulensModel's own binary-lens model to the same…, Chi2 of MulensModel's binary magnification against the combined OGLE+MOA data., Serial multi-start search, then a tight-tolerance refit of the best seed., Overlay the MulensModel fit on the OGLE+MOA data, full baseline + zoomed on the… (+1 more)

### Community 9 - "Community 9"
Cohesion: 0.50
Nodes (4): load_coords(), _parse_header_coord(), Extract (RA, Dec) sexagesimal strings from a raw table's own \\RA/\\DEC header…, Target coordinates, parsed from OGLE_PATH/MOA_PATH's own headers rather than…

## Knowledge Gaps
- **6 isolated node(s):** `submit_all.sh script`, `astropy`, `corner`, `emcee`, `torch` (+1 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 171 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **5 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `binary_magnification_fs()` connect `Binary-Lens Magnification` to `Event Loading & Flux Calibration`, `Cassan & Binary Preprocessing`, `Search Pipeline & Diagnostics`, `Binary-Lens Magnification`?**
  _High betweenness centrality (0.061) - this node is a cross-community bridge._
- **Why does `TwoL1SParams` connect `Cassan & Binary Preprocessing` to `Event Loading & Flux Calibration`, `Search Pipeline & Diagnostics`?**
  _High betweenness centrality (0.057) - this node is a cross-community bridge._
- **Why does `Config-driven pipeline (input/<short>.toml -> search.py)` connect `Search Pipeline & Diagnostics` to `Search Pipeline & Diagnostics`?**
  _High betweenness centrality (0.056) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `TwoL1SParams` (e.g. with `TwoL1SParams NamedTuple (canonical 2L1S param order)` and `fit()`) actually correct?**
  _`TwoL1SParams` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `submit_all.sh script`, `astropy`, `corner` to the rest of the system?**
  _6 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Event Loading & Flux Calibration` be split into smaller, more focused modules?**
  _Cohesion score 0.05207835642618251 - nodes in this community are weakly interconnected._
- **Should `Cassan & Binary Preprocessing` be split into smaller, more focused modules?**
  _Cohesion score 0.05875251509054326 - nodes in this community are weakly interconnected._