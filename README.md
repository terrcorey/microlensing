# Microlensing light-curve fitting playground

A playground for analysing real, published microlensing light curves, to
learn the challenges involved in modeling brown dwarf populations with
microlensing. Not an installable application -- each script is run
directly and judged by inspecting the plot/printed fit it produces.

The central open question this project is working towards: a microlensing
timescale (`tE`) alone can't separate a lens's mass from its distance, which
is exactly why brown-dwarf-mass lenses are hard to identify from
microlensing alone. Everything below -- the model, the priors, the derived
quantities -- exists to make that degeneracy visible, not to hide it behind
a single confident number.

## Setup

```
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Running it

Both fit scripts below take `--stage`, so you can go as far as you want:
`raw` (just the raw_lc/ plot), `quicklook` (+ fit_lc/, point-estimate only,
no MCMC), or `mcmc` (default: the full run, + corner_plots/, hist_plots/,
a fit_summary.dat).

```
python3 download_data.py             # fetch raw photometry into data/

python3 mcmc_fit.py --stage=raw        # raw_lc/O-05-BLG086.png
python3 mcmc_fit.py --stage=quicklook  # + fit_lc/O-05-BLG086.png
python3 mcmc_fit.py                    # full run (default)

python3 preprocess_binary_data.py      # required before the O-03-BLG235 commands below
python3 mcmc_fit_binary.py --stage=raw
python3 mcmc_fit_binary.py --stage=quicklook
python3 mcmc_fit_binary.py             # full run (default) -- see "Known limitations" below
```

## The model and its assumptions

- **Point-source point-lens (PSPL) only.** No finite-source, parallax, or
  binary-lens physics (`lc_models.py`) -- applied even to O-03-BLG235, a
  real binary-lens event, deliberately (see "Known limitations" below).
- **Zero-point magnitude (18.0) is arbitrary.** It cancels out of every
  fitted flux ratio; it is not a calibrated zeropoint.
- **Time is `HJD - 2450000` everywhere.** Files already in that convention
  are used as-is; OGLE/MOA's O-03-BLG235 tables (full JD/HJD) get
  `- 2450000` applied right after loading.
- **O-03-BLG235's two instruments get different flux models:**
  `OGLE flux = fs_ogle * A(t) + fb_ogle`, `MOA flux = fs_moa * (A(t) - 1)`
  -- MOA has no blend term because its differencing already cancels any
  constant blend flux.
- **O-03-BLG235's flux calibration is frozen at a point estimate.**
  `preprocess_binary_data.py` runs its joint fit once to get
  `(fs_ogle, fb_ogle, fs_moa)`, then treats that as fixed for every
  downstream fit. The resulting `(t0, u0, tE)` posterior is therefore
  tighter than a fully joint 6-parameter MCMC would give -- it doesn't
  propagate calibration uncertainty. Fine for exploration; would need a
  true joint MCMC for a rigorous error budget.
- **`u0`/`tE` initial guesses are per-dataset, not shared.** The model only
  depends on `u0` squared, so a bad starting point can converge to the
  mirror-image (`u0 -> -u0`) solution. This is a real bug that was hit
  once (see CHANGELOG.md) -- silently reusing one dataset's guess for
  another gave a wrong-sign, artificially tight posterior that still
  looked converged.
- **MCMC priors are flat bounds, not physically motivated:** `t0` within
  the data's time span, `0 < u0 < 5`, `0.1 < tE < 1000` days, `f_source >
  0`. They exist to keep the sampler off unphysical territory, not to
  encode real prior knowledge.
- **MCMC settings** (32 walkers, 4000 steps, 1000 discarded as burn-in,
  thinned by 15) were chosen empirically -- no convergence diagnostic
  (e.g. autocorrelation time, Gelman-Rubin) is run.
- **Reported uncertainties are posterior percentiles, not Gaussian
  errors.** Every `p16/p50/p84` triplet is read directly off the
  (burned-in, thinned) samples, so asymmetric spreads are real signal, not
  a formatting quirk -- see `M_lens` below.
- **The lens mass estimate (`M_lens`) is prior-dominated, not primarily
  data-driven.** `estimate_mass()` in `mcmc_fit.py` pairs each posterior
  `tE` sample with one Monte Carlo draw from assumed priors, since `tE`
  alone can't fix the mass:
  - source distance `D_S` fixed at 8 kpc (standard Galactic bulge distance
    -- both events are OGLE bulge fields);
  - lens distance `D_L` ~ Uniform(0.1, `D_S`) kpc -- a naive line-of-sight
    prior, **not** a real Galactic density model. This is the single
    biggest caveat in the whole pipeline: it's what gives `M_lens` its
    long right tail (`D_L -> D_S` implies `pi_rel -> 0` implies `M -> `
    very large);
  - relative proper motion `mu_rel` ~ log-normal, median 4 mas/yr, ~1.5x
    scatter -- a typical literature value for bulge microlensing, not
    measured for these specific events;
  - physics: `theta_E = tE * mu_rel`, `pi_rel = 1/D_L[kpc] - 1/D_S[kpc]`
    (mas), `M = theta_E^2 / (kappa * pi_rel)`, `kappa = 8.144 mas/Msun`.

  **Treat `M_lens` as illustrating the mass-distance degeneracy, not as a
  real mass measurement of these lenses.**
- **Raw data is never hand-edited.** Everything under `data/processed/` is
  regenerated by re-running the relevant script. Raw files under `data/`
  come from the MulensModel project's redistribution of OGLE/MOA
  photometry (`download_data.py`), not directly from OGLE's own archive.

## How to read the data files

### Raw (`data/`)

| File | Format |
|---|---|
| `OGLE-2005-BLG-086.dat` | plain 3-column text: `time(HJD-2450000)  I_mag  mag_err`. Single instrument, already in project time convention. |
| `OGLE-2003-BLG-235_OGLE.tbl.txt` | NASA Exoplanet Archive format (metadata lines start with `\` or `|`): `JD  relative_I_mag  mag_uncertainty`. Full JD -- needs `-2450000`. |
| `OGLE-2003-BLG-235_MOA.tbl.txt` | same table format: `HJD  relative_flux  flux_uncertainty`. DIA-style flux relative to a fiducial template -- can be negative, no shared zero point with OGLE's magnitudes. Full HJD -- needs `-2450000`. |

### Processed (`data/processed/` -- regenerated by scripts, never hand-edited)

| File | Format |
|---|---|
| `O-03-BLG235_{OGLE,MOA}_magnification.dat` | `time(HJD-2450000)  magnification_A  magnification_A_err`. Both instruments converted onto the same physical scale using the frozen joint-fit calibration. |
| `{name}_fit_summary.dat` (e.g. `O-05-BLG086_fit_summary.dat`, `O-03-BLG235_fit_summary.dat`, `O-03-BLG235_ogle_only_fit_summary.dat`) | One row per parameter: `param  p16  p50  p84` (68% credible interval, see percentiles above). Raw MCMC parameters (`t0`, `u0`, `tE`, and `fs`/`fb` where fit) plus derived quantities: `A_max` (peak magnification), `t_eff` (`u0*tE`), `blend_fraction` (`fb/(fs+fb)`, only where fs/fb were fit), `m_source` (source baseline magnitude, only where fs/fb were fit), `M_lens` (prior-based mass estimate -- read with the caution above). |

### Plots

| Directory | Contents |
|---|---|
| `raw_lc/{name}.png` | light curve as observed, no fit overlay. |
| `fit_lc/{name}.png` | data + PSPL fit, full baseline and auto-zoomed peak panels. |
| `corner_plots/{name}.png` | pairwise posterior correlations + 1D marginals, raw MCMC parameters only. |
| `hist_plots/{name}.png` | 1D histogram per raw *and* derived quantity (including `M_lens`), log-x-axis auto-applied to heavy-tailed positive quantities. Exists to catch bimodality/skew that a `p16/p50/p84` triplet alone could hide. |

`raw_lc/`, `fit_lc/`, `hist_plots/`, `corner_plots/` hold only these regular
pipeline outputs. Everything from the one-time/dev scripts in `scratch/`
(2L1S fit, MulensModel cross-check) writes there instead, not into the
four directories above.

### `dataset_names.txt`

Short-name convention (`O-05-BLG086`, `O-03-BLG235`) used in every output
filename. Always add new mappings there rather than inventing an
abbreviation inline.

## Known limitations of each fit

- **O-05-BLG086**: independently reproduces the published solution
  (Wyrzykowski et al. 2015) to ~2% -- validated.
- **O-03-BLG235**: a real binary-lens event; PSPL is a known-incomplete
  model here (visible as points sitting above the fitted curve in the
  auto-zoomed panel, most clearly in MOA's higher-cadence data). The
  `ogle_only` diagnostic exists specifically to show that a single sparse
  instrument can still produce a clean, plausible-looking, unimodal PSPL
  fit to a real binary-lens event -- neither chi2/dof nor posterior shape
  flags the misspecification on their own.

See `CHANGELOG.md` for the session-by-session history of what was built,
what was learned, and open questions.
