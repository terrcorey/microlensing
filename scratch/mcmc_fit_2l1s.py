"""MCMC (emcee) exploration of the full 2L1S parameter space for
OGLE-2003-BLG-235/MOA-2003-BLG-53, as an alternative to fit_2l1s.py's
multi-start Nelder-Mead. Walkers start broadly across log-uniform s/q and
full-circle alpha instead of a handful of grid seeds -- multi-start kept
finding the same wrong, caustic-missing basin because its grid never
seeded q anywhere near the true ~0.004 decade; log-uniform coverage across
~4 decades fixes that directly rather than guessing a better grid.

Not yet graduated into the regular pipeline (see CLAUDE.md) -- lives in
scratch/ with the project's other one-time/dev scripts, run manually as
`python3 scratch/mcmc_fit_2l1s.py` from the project root.
"""
import sys
from multiprocessing import get_context
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import emcee
import numpy as np
from scipy.stats import t as student_t
from fit_2l1s import SHORT_NAME, chi2, guess, plot_fit, residuals
from mcmc_fit import save_corner

LABELS = ["t0", "u0", "tE", "alpha", "piE_N", "piE_E", "s", "q", "scale", "dof"]

T0_WIDTH = 10.0  # days around the PSPL joint fit's t0
U0_RANGE = (0.01, 1.0)
TE_RANGE = (10.0, 150.0)
S_RANGE = (0.7, 1.8)
LOG_Q_RANGE = (np.log(1e-4), np.log(1.0))
PIE_RANGE = (-2.0, 2.0)
LOG_SCALE_RANGE = (np.log(0.1), np.log(10.0))
LOG_DOF_RANGE = (np.log(0.5), np.log(50.0))


def log_prior(theta):
    """Flat prior; log(s) and log(q) uniform so orders of magnitude get equal weight."""
    t0, u0, tE, alpha, piE_N, piE_E, s, q, scale, dof= theta
    if not (guess["t0"] - T0_WIDTH < t0 < guess["t0"] + T0_WIDTH):
        return -np.inf
    if not (U0_RANGE[0] < u0 < U0_RANGE[1]):
        return -np.inf
    if not (TE_RANGE[0] < tE < TE_RANGE[1]):
        return -np.inf
    if not (0 <= alpha < 2 * np.pi):
        return -np.inf
    if s <= 0 or q <= 0:
        return -np.inf
    if not (S_RANGE[0] < s < S_RANGE[1]):
        return -np.inf
    if not (LOG_Q_RANGE[0] < np.log(q) < LOG_Q_RANGE[1]):
        return -np.inf
    if not (PIE_RANGE[0] < piE_N < PIE_RANGE[1]):
        return -np.inf
    if not (PIE_RANGE[0] < piE_E < PIE_RANGE[1]):
        return -np.inf
    if not (LOG_SCALE_RANGE[0] < np.log(scale) < LOG_SCALE_RANGE[1]):
        return -np.inf
    if not (LOG_DOF_RANGE[0] < np.log(dof) < LOG_DOF_RANGE[1]):
        return -np.inf   
    return 0.0


def log_probability(theta, use_ogle=True):
    t0, u0, tE, alpha, piE_N, piE_E, s, q, scale, dof= theta
    binary_theta = (t0, u0, tE, alpha, piE_N, piE_E, s, q)
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    resids = residuals(binary_theta, use_ogle)
    t = student_t.logpdf(resids/scale, df=dof).sum() - resids.size * np.log(scale)
    return -np.inf if not np.isfinite(t) else lp + t


def sample_prior(rng, n):
    """Draw n walker start positions directly from the prior -- broad coverage in place of grid seeds."""
    t0 = rng.uniform(guess["t0"] - T0_WIDTH, guess["t0"] + T0_WIDTH, n)
    u0 = rng.uniform(*U0_RANGE, n)
    tE = rng.uniform(*TE_RANGE, n)
    alpha = rng.uniform(0, 2 * np.pi, n)
    piE_N = rng.uniform(*PIE_RANGE, n)
    piE_E = rng.uniform(*PIE_RANGE, n)
    s = rng.uniform(*S_RANGE, n)
    q = np.exp(rng.uniform(*LOG_Q_RANGE, n))
    scale = np.exp(rng.uniform(*LOG_SCALE_RANGE, n))
    dof = np.exp(rng.uniform(*LOG_DOF_RANGE, n))
    return np.column_stack([t0, u0, tE, alpha, piE_N, piE_E, s, q, scale, dof])


def run_mcmc(nwalkers=48, nsteps=1500, seed=42, use_ogle=True):
    tag = "mcmc-2l1s" if use_ogle else "mcmc-2l1s-moa-only"
    suffix = "" if use_ogle else "_moa_only"

    rng = np.random.default_rng(seed)
    p0 = sample_prior(rng, nwalkers)

    # spawn (not fork) -- avoids CUDA-in-forked-subprocess issues once binary_magnification runs on GPU
    with get_context("spawn").Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, len(LABELS), log_probability, pool=pool, args=(use_ogle,))
        sampler.run_mcmc(p0, nsteps, progress=True)

    samples = sampler.get_chain(discard=nsteps // 4, thin=15, flat=True)
    log_probs = sampler.get_log_prob(discard=nsteps // 4, thin=15, flat=True)
    best = samples[np.argmax(log_probs)]

    print(f"[{tag}] {samples.shape[0]} posterior samples")
    print(f"[{tag}] best sample: chi2={chi2(best, use_ogle):.2f}")
    for label, value in zip(LABELS, best):
        print(f"[{tag}] {label} = {value:.5f}")

    plot_fit(best, use_ogle=use_ogle)

    save_corner(samples, LABELS, best, f"scratch/{SHORT_NAME}_2l1s_mcmc{suffix}.png")
    np.savez(f"scratch/{SHORT_NAME}_2l1s_mcmc_chain{suffix}.npz", samples=samples, log_probs=log_probs, labels=LABELS)

    return sampler, samples


def run_mcmc_moa_only():
    return run_mcmc(use_ogle=False)


if __name__ == "__main__":
    run_mcmc()
    run_mcmc_moa_only()
