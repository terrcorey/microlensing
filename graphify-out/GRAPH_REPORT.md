# Graph Report - microlensing  (2026-10-05)

## Corpus Check
- 24 files · ~43,051 words
- Verdict: corpus is large enough that graph structure adds value.
- Unclassified: 18 file(s) not represented in the graph (top: .sbatch 8, .npz 7, (none) 1)

## Summary
- 408 nodes · 957 edges · 16 communities (12 shown, 4 thin omitted)
- Extraction: 95% EXTRACTED · 5% INFERRED · 0% AMBIGUOUS · INFERRED: 47 edges (avg confidence: 0.86)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `934b4d49`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- CHANGELOG (session log)
- fit_2l1s.py
- mcmc_fit.py
- mcmc_fit_2l1s.py
- cross_check_mulensmodel.py
- search.py
- binary_trajectory
- derive_binary_quintic.py
- M-19-BLG008 (MOA-2019-BLG-008, KMT I-band subset)
- fit_parallax_pspl_mcmc
- lc_models.py
- caustic_curve
- _binary_images_batch
- cross_check_mulensmodel_fs.py
- submit_all.sh

## God Nodes (most connected - your core abstractions)
1. `CHANGELOG (session log)` - 27 edges
2. `CLAUDE.md (project guide)` - 26 edges
3. `plot_fit()` - 21 edges
4. `Event` - 19 edges
5. `magnification()` - 19 edges
6. `binary_trajectory()` - 18 edges
7. `TwoL1SParams` - 18 edges
8. `trajectory()` - 17 edges
9. `caustic_curve()` - 16 edges
10. `flat_chain()` - 16 edges

## Surprising Connections (you probably didn't know these)
- `estimate_mass()` --implements--> `Prior-based lens mass estimate (M_lens)`  [EXTRACTED]
  mcmc_fit.py → README.md
- `u0 mirror-sign degeneracy bug (per-dataset initial guesses)` --semantically_similar_to--> `MulensModel oracle cross-checks`  [INFERRED] [semantically similar]
  CHANGELOG.md → CLAUDE.md
- `Error-bar rescaling / additive jitter` --semantically_similar_to--> `Student-t robust likelihood (free scale/dof)`  [INFERRED] [semantically similar]
  CHANGELOG.md → CLAUDE.md
- `SIGNALMEN anomaly-vs-outlier test (Dominik et al. 2007)` --semantically_similar_to--> `Huber loss likelihood (DELTA=1.345)`  [INFERRED] [semantically similar]
  CHANGELOG.md → CLAUDE.md
- `_model_magnitude()` --calls--> `magnitude()`  [EXTRACTED]
  mcmc_fit.py → lc_models.py

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **2L1S roadmap: finite source, error bars, Cassan, (d,q) grid, GA, SIGNALMEN** — claude_finite_source, changelog_error_bar_rescaling, claude_cassan_parametrisation, changelog_dq_grid_search, changelog_genetic_algorithm, changelog_signalmen [EXTRACTED 1.00]
- **2L1S model validation via limiting cases + MulensModel oracle** — claude_2l1s_model, claude_finite_source, claude_annual_parallax, claude_mulensmodel_oracle [EXTRACTED 1.00]
- **Likelihood choices explored (chi2 / Student-t / Huber)** — claude_student_t_likelihood, claude_huber_likelihood, changelog_error_bar_rescaling, claude_mcmc_emcee [INFERRED 0.85]

## Communities (16 total, 4 thin omitted)

### Community 0 - "CHANGELOG (session log)"
Cohesion: 0.06
Nodes (59): CHANGELOG (session log), Bachelet et al. 2022, BIC PSPL-vs-2L1S verdict (deltaBIC>10), Bond et al. 2004 (O-03-BLG235 planet solution), Cassan 2008, Charbonneau 1995, Dominik et al. 2007 (MNRAS 380, 792), (d,q) grid search (+51 more)

### Community 1 - "fit_2l1s.py"
Cohesion: 0.06
Nodes (58): concurrent_futures, functools, mag_to_flux(), magnification(), Calibrated magnitude + error -> flux + error on the ZERO_POINT_MAG scale., Source-lens separation u(t), in units of the Einstein radius., Paczynski point-lens magnification A(u)., trajectory() (+50 more)

### Community 2 - "mcmc_fit.py"
Cohesion: 0.07
Nodes (52): argparse, corner, Download raw photometry for each event we're working with. Sources are…, emcee, itertools, magnitude(), Returns (delta_s_n, delta_s_e): Earth's sky-projected position relative to the…, Observed magnitude, for overplotting against real photometry. (+44 more)

### Community 3 - "mcmc_fit_2l1s.py"
Cohesion: 0.07
Nodes (59): EnsembleSampler, flat_chain(), ndarray, Flat (samples, log_probs) after burn-in/thinning -- shared by every emcee…, Corner plot of the raw MCMC parameters -- shared by every fit script.…, save_corner(), multiprocessing, os (+51 more)

### Community 4 - "cross_check_mulensmodel.py"
Cohesion: 0.40
Nodes (5): mulensmodel, mm_magnification(), One-time offline cross-check of lc_models.binary_magnification against…, MulensModel's magnification at the same complex zeta values. MulensModel's docs…, report()

### Community 5 - "search.py"
Cohesion: 0.06
Nodes (68): chi2(), error_scale(), Event, flux_residuals(), Instrument, load_event(), load_instrument(), NamedTuple (+60 more)

### Community 6 - "binary_trajectory"
Cohesion: 0.18
Nodes (16): binary_magnification(), binary_trajectory(), Total binary-lens magnification at source position(s) `zeta`. Sums 1/|J| over…, Source-lens separation for binary lenses on the lens plane, given in complex…, plot_panel(), chi2(), _fit_one_seed(), plot_fit() (+8 more)

### Community 10 - "fit_parallax_pspl_mcmc"
Cohesion: 0.14
Nodes (20): flux(), plain_flux(), plain_magnitude(), Observed flux: magnified source plus a constant blend flux., flux() with piE_N=piE_E=0 -- the exact pre-parallax model, used only to…, magnitude() with piE_N=piE_E=0 -- see plain_flux()., fit_parallax_pspl_mcmc(), log_likelihood() (+12 more)

### Community 11 - "lc_models.py"
Cohesion: 0.13
Nodes (16): astropy_coordinates, astropy_time, astropy_units, _project(), Functional forms used to model a microlensing light curve. Point-source point-…, Inverse of cassan_to_standard(): every (sigma, t) where the parallax-free…, Project a Cartesian position/velocity vector onto the sky's (N, E) tangent-…, standard_to_cassan() (+8 more)

### Community 12 - "caustic_curve"
Cohesion: 0.23
Nodes (12): complex128, binary_magnification_fs(), _disk_average(), _sample_coords(), cassan_caustic(), caustic_curve(), equidistant_caustic(), NDArray (+4 more)

### Community 13 - "_binary_images_batch"
Cohesion: 0.20
Nodes (10): binary_images(), _binary_images_batch(), _companion_eigvals(), lens_position(), _quintic_coefficients(), Coefficients (highest degree first) of the degree-5 polynomial in z equivalent…, Roots of a batch of degree-5 polynomials via a batched companion-matrix…, Shared batched core for binary_images/binary_magnification. Solves the binary… (+2 more)

### Community 14 - "cross_check_mulensmodel_fs.py"
Cohesion: 0.40
Nodes (3): caustic_point(), Point at abscissa sigma (wraps with period 1) on a cassan_caustic() output., One-time offline cross-check of lc_models.binary_magnification_fs against…

## Knowledge Gaps
- **7 isolated node(s):** `submit_all.sh script`, `M-19-BLG008 (MOA-2019-BLG-008, KMT I-band subset)`, `Time = HJD - 2450000 convention`, `matplotlib`, `numpy` (+2 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 159 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **4 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `estimate_mass()` connect `mcmc_fit.py` to `CHANGELOG (session log)`?**
  _High betweenness centrality (0.251) - this node is a cross-community bridge._
- **Why does `Prior-based lens mass estimate (M_lens)` connect `CHANGELOG (session log)` to `mcmc_fit.py`?**
  _High betweenness centrality (0.244) - this node is a cross-community bridge._
- **Are the 10 inferred relationships involving `Event` (e.g. with `anomaly_times()` and `bic()`) actually correct?**
  _`Event` has 10 INFERRED edges - model-reasoned connections that need verification._
- **What connects `submit_all.sh script`, `M-19-BLG008 (MOA-2019-BLG-008, KMT I-band subset)`, `Time = HJD - 2450000 convention` to the rest of the system?**
  _7 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `CHANGELOG (session log)` be split into smaller, more focused modules?**
  _Cohesion score 0.059887005649717516 - nodes in this community are weakly interconnected._
- **Should `fit_2l1s.py` be split into smaller, more focused modules?**
  _Cohesion score 0.059395801331285206 - nodes in this community are weakly interconnected._
- **Should `mcmc_fit.py` be split into smaller, more focused modules?**
  _Cohesion score 0.06957047791893527 - nodes in this community are weakly interconnected._