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
import os
import sys
from multiprocessing import get_context
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import emcee
import numpy as np
from scipy.stats import t as student_t
from fit_2l1s import SHORT_NAME, TwoL1SParams, chi2, guess, plot_fit, residuals, huber
from mcmc_fit import save_corner

# Derived from TwoL1SParams._fields rather than hand-typed three times, so the
# physical ordering can't drift between the three likelihood variants below.
LABELS = list(TwoL1SParams._fields) + ["scale", "dof"]
LABELS_CHI2 = list(TwoL1SParams._fields)  # no scale/dof -- Gaussian chi2, not Student-t
LABELS_HUBER = list(TwoL1SParams._fields) + ["scale"]

T0_WIDTH = 10.0  # days around the PSPL joint fit's t0
U0_RANGE = (0.01, 1.0)
TE_RANGE = (10.0, 150.0)
S_RANGE = (0.5, 2.0)
LOG_Q_RANGE = (np.log(1e-5), np.log(1.0))
PIE_RANGE = (-2.0, 2.0)
LOG_SCALE_RANGE = (np.log(0.1), np.log(10.0))
LOG_DOF_RANGE = (np.log(0.5), np.log(50.0))
DELTA = 1.345


def _physical_log_prior(params):
    """Flat-prior bound checks on the 8 physical params, shared by all three 2L1S
    likelihoods below (Student-t, chi2, Huber) -- previously an identical 9-line
    if-chain copy-pasted three times as log_prior/log_prior_chi2/log_prior_huber.
    log(s)/log(q) uniform so orders of magnitude get equal weight."""
    t0, u0, tE, alpha, piE_N, piE_E, s, q = params
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
    return 0.0


def _log_range_ok(x, log_range):
    return log_range[0] < np.log(x) < log_range[1]


def log_prior(theta):
    """Student-t variant: physical bounds plus scale/dof."""
    params = TwoL1SParams(*theta[:8])
    scale, dof = theta[8], theta[9]
    lp = _physical_log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    if not _log_range_ok(scale, LOG_SCALE_RANGE):
        return -np.inf
    if not _log_range_ok(dof, LOG_DOF_RANGE):
        return -np.inf
    return 0.0


def log_probability(theta, use_ogle=True):
    """Student-t log-likelihood (heavier-tailed than Gaussian chi2), with the
    `-log(scale)` Jacobian term for the r -> r/scale change of variables -- dropping it
    would let the sampler inflate `scale` for free (see conversation)."""
    params = TwoL1SParams(*theta[:8])
    scale, dof = theta[8], theta[9]
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    resids = residuals(params, use_ogle)
    if not np.all(np.isfinite(resids)):
        return -np.inf
    t = student_t.logpdf(resids/scale, df=dof).sum() - resids.size * np.log(scale)
    return -np.inf if not np.isfinite(t) else lp + t


def _sample_physical_prior(rng, n):
    """Draw n walker start columns for the 8 physical params -- shared by all three
    likelihoods' sample_prior* below, previously copy-pasted three times."""
    t0 = rng.uniform(guess["t0"] - T0_WIDTH, guess["t0"] + T0_WIDTH, n)
    u0 = rng.uniform(*U0_RANGE, n)
    tE = rng.uniform(*TE_RANGE, n)
    alpha = rng.uniform(0, 2 * np.pi, n)
    piE_N = rng.uniform(*PIE_RANGE, n)
    piE_E = rng.uniform(*PIE_RANGE, n)
    s = rng.uniform(*S_RANGE, n)
    q = np.exp(rng.uniform(*LOG_Q_RANGE, n))
    return [t0, u0, tE, alpha, piE_N, piE_E, s, q]


def sample_prior(rng, n):
    """Draw n walker start positions directly from the prior -- broad coverage in place of grid seeds."""
    scale = np.exp(rng.uniform(*LOG_SCALE_RANGE, n))
    dof = np.exp(rng.uniform(*LOG_DOF_RANGE, n))
    return np.column_stack(_sample_physical_prior(rng, n) + [scale, dof])


def run_mcmc(nwalkers=48, nsteps=1500, seed=42, use_ogle=True):
    tag = "mcmc-2l1s" if use_ogle else "mcmc-2l1s-moa-only"
    suffix = "" if use_ogle else "_moa_only"

    rng = np.random.default_rng(seed)
    p0 = sample_prior(rng, nwalkers)

    # explicit process count -- os.cpu_count()/mp's own default reads the physical node's
    # core count, not what SLURM's cgroup actually allocated to this job, which oversubscribes
    # on a shared cluster node; SLURM_CPUS_PER_TASK is unset when run outside SLURM (falls back
    # to the default). spawn (not fork) avoids CUDA-in-forked-subprocess issues once
    # binary_magnification runs on GPU.
    nprocs = int(os.environ["SLURM_CPUS_PER_TASK"]) if "SLURM_CPUS_PER_TASK" in os.environ else None
    with get_context("spawn").Pool(processes=nprocs) as pool:
        sampler = emcee.EnsembleSampler(nwalkers, len(LABELS), log_probability, pool=pool, args=(use_ogle,))
        sampler.run_mcmc(p0, nsteps, progress=True)

    samples = sampler.get_chain(discard=nsteps // 4, thin=15, flat=True)
    log_probs = sampler.get_log_prob(discard=nsteps // 4, thin=15, flat=True)
    best = samples[np.argmax(log_probs)]
    binary_best = TwoL1SParams(*best[:8])  # drop scale/dof -- chi2()/plot_fit() take only the 8 physical params

    print(f"[{tag}] {samples.shape[0]} posterior samples")
    print(f"[{tag}] best sample: chi2={chi2(binary_best, use_ogle):.2f}")
    for label, value in zip(LABELS, best):
        print(f"[{tag}] {label} = {value:.5f}")

    plot_fit(binary_best, use_ogle=use_ogle, tag="_mcmc_studentt")

    save_corner(samples, LABELS, best, f"scratch/{SHORT_NAME}_2l1s_mcmc{suffix}.png")
    np.savez(f"scratch/{SHORT_NAME}_2l1s_mcmc_chain{suffix}.npz", samples=samples, log_probs=log_probs, labels=LABELS)

    return sampler, samples


def run_mcmc_moa_only():
    return run_mcmc(use_ogle=False)


def log_prior_chi2(theta):
    """Same bounds as log_prior(), minus scale/dof -- this run has no Student-t params."""
    return _physical_log_prior(TwoL1SParams(*theta))


def log_probability_chi2(theta, use_ogle=True):
    """Plain Gaussian likelihood (-0.5*chi2) in place of log_probability()'s Student-t --
    a direct comparison run: does the caustic-crossing spike drag the whole fit toward
    itself under chi2's unbounded quadratic penalty, the way it can't under Student-t's
    logarithmic tail (see conversation)."""
    lp = log_prior_chi2(theta)
    if not np.isfinite(lp):
        return -np.inf
    c2 = chi2(theta, use_ogle)
    return -np.inf if not np.isfinite(c2) else lp - 0.5 * c2


def sample_prior_chi2(rng, n):
    """Draw n walker start positions directly from the prior -- broad coverage in place of grid seeds."""
    return np.column_stack(_sample_physical_prior(rng, n))


def run_mcmc_chi2(nwalkers=48, nsteps=1500, seed=42, use_ogle=True):
    tag = "mcmc-2l1s-chi2" if use_ogle else "mcmc-2l1s-chi2-moa-only"
    suffix = ("" if use_ogle else "_moa_only") + "_chi2"

    rng = np.random.default_rng(seed)
    p0 = sample_prior_chi2(rng, nwalkers)

    nprocs = int(os.environ["SLURM_CPUS_PER_TASK"]) if "SLURM_CPUS_PER_TASK" in os.environ else None
    with get_context("spawn").Pool(processes=nprocs) as pool:
        sampler = emcee.EnsembleSampler(nwalkers, len(LABELS_CHI2), log_probability_chi2, pool=pool, args=(use_ogle,))
        sampler.run_mcmc(p0, nsteps, progress=True)

    samples = sampler.get_chain(discard=nsteps // 4, thin=15, flat=True)
    log_probs = sampler.get_log_prob(discard=nsteps // 4, thin=15, flat=True)
    best = samples[np.argmax(log_probs)]

    print(f"[{tag}] {samples.shape[0]} posterior samples")
    print(f"[{tag}] best sample: chi2={chi2(best, use_ogle):.2f}")
    for label, value in zip(LABELS_CHI2, best):
        print(f"[{tag}] {label} = {value:.5f}")

    plot_fit(best, use_ogle=use_ogle, tag="_chi2")

    save_corner(samples, LABELS_CHI2, best, f"scratch/{SHORT_NAME}_2l1s_mcmc{suffix}.png")
    np.savez(f"scratch/{SHORT_NAME}_2l1s_mcmc_chain{suffix}.npz", samples=samples, log_probs=log_probs, labels=LABELS_CHI2)

    return sampler, samples


def log_prior_huber(theta):
    """Same bounds as log_prior(), minus dof -- Huber has no dof-like shape parameter
    (delta is fixed, not fit)."""
    params = TwoL1SParams(*theta[:8])
    scale = theta[8]
    lp = _physical_log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    if not _log_range_ok(scale, LOG_SCALE_RANGE):
        return -np.inf
    return 0.0


def log_probability_huber(theta, use_ogle=True):
    """Huber loss in place of log_probability()'s Student-t -- quadratic near zero,
    linear past DELTA (see conversation). DELTA is fixed, not fit, so its own
    normalizing constant is a dropped constant; `scale` is still free, so it keeps
    the same `-log(scale)` Jacobian term log_probability() has."""
    params = TwoL1SParams(*theta[:8])
    scale = theta[8]
    lp = log_prior_huber(theta)
    if not np.isfinite(lp):
        return -np.inf
    resids = residuals(params, use_ogle)
    if not np.all(np.isfinite(resids)):
        return -np.inf
    log_prob = - huber(resids/scale, DELTA) - resids.size * np.log(scale)
    return -np.inf if not np.isfinite(log_prob) else lp + log_prob


def sample_prior_huber(rng, n):
    """Draw n walker start positions directly from the prior -- broad coverage in place of grid seeds."""
    scale = np.exp(rng.uniform(*LOG_SCALE_RANGE, n))
    return np.column_stack(_sample_physical_prior(rng, n) + [scale])


def run_mcmc_huber(nwalkers=48, nsteps=1500, seed=42, use_ogle=True):
    tag = "mcmc-2l1s-huber" if use_ogle else "mcmc-2l1s-huber-moa-only"
    suffix = ("" if use_ogle else "_moa_only") + "_huber"

    rng = np.random.default_rng(seed)
    p0 = sample_prior_huber(rng, nwalkers)

    nprocs = int(os.environ["SLURM_CPUS_PER_TASK"]) if "SLURM_CPUS_PER_TASK" in os.environ else None
    with get_context("spawn").Pool(processes=nprocs) as pool:
        sampler = emcee.EnsembleSampler(nwalkers, len(LABELS_HUBER), log_probability_huber, pool=pool, args=(use_ogle,))
        sampler.run_mcmc(p0, nsteps, progress=True)

    samples = sampler.get_chain(discard=nsteps // 4, thin=15, flat=True)
    log_probs = sampler.get_log_prob(discard=nsteps // 4, thin=15, flat=True)
    best = samples[np.argmax(log_probs)]
    binary_best = TwoL1SParams(*best[:8])  # drop scale -- chi2()/residuals()/plot_fit() take only the 8 physical params

    print(f"[{tag}] {samples.shape[0]} posterior samples")
    print(f"[{tag}] best sample: huber_loss={huber(residuals(binary_best, use_ogle), DELTA):.2f}")
    for label, value in zip(LABELS_HUBER, best):
        print(f"[{tag}] {label} = {value:.5f}")

    plot_fit(binary_best, use_ogle=use_ogle, tag="_huber")

    save_corner(samples, LABELS_HUBER, best, f"scratch/{SHORT_NAME}_2l1s_mcmc{suffix}.png")
    np.savez(f"scratch/{SHORT_NAME}_2l1s_mcmc_chain{suffix}.npz", samples=samples, log_probs=log_probs, labels=LABELS_HUBER)

    return sampler, samples


if __name__ == "__main__":
    run_mcmc()
    run_mcmc_moa_only()
