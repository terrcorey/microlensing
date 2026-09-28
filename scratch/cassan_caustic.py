"""2L1S fit of OGLE-2003-BLG-235 in Cassan (2008)'s caustic-crossing parametrisation.

Fits (sigma_in, sigma_out, t_in, t_out, s, q) instead of (t0, u0, tE, alpha, s, q):
sigma_in/sigma_out are the entry/exit points' curvilinear abscissa along the caustic
(lc_models.cassan_caustic()), t_in/t_out the times the source crosses them. Every
trial trajectory passes through the caustic by construction, so the search can't
fall into the near-miss-cusp basins that multi-start Nelder-Mead kept finding (see
CHANGELOG). Parallax-free (piE_N = piE_E = 0) for now -- the mapping to a straight
trajectory is only exact without it.

Run from the project root: `python3 scratch/cassan_caustic.py`.
"""
import sys
from pathlib import Path
from typing import NamedTuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import emcee
import numpy as np
from itertools import product
from multiprocessing import get_context
from scipy.optimize import minimize

from lc_models import caustic_curve, cassan_caustic, cassan_to_standard, standard_to_cassan
from fit_2l1s import SHORT_NAME, TwoL1SParams, chi2, plot_fit
from mcmc_fit import save_corner
from mcmc_fit_2l1s import S_RANGE, LOG_Q_RANGE, TE_RANGE

T_WINDOW = (2820.0, 2870.0)  # t_in/t_out prior range: the anomaly with room either side
OUT_DIR = Path("scratch/2l1s/cassan")


class CassanParams(NamedTuple):
    sigma_in: float
    sigma_out: float
    t_in: float
    t_out: float
    s: float
    q: float


def to_standard(theta, caustic=None):
    """CassanParams -> the parallax-free TwoL1SParams residuals()/plot_fit() take."""
    sigma_in, sigma_out, t_in, t_out, s, q = CassanParams(*theta)
    caustic = cassan_caustic(s, q) if caustic is None else caustic
    t0, u0, tE, alpha = cassan_to_standard(caustic, sigma_in, sigma_out, t_in, t_out)
    return TwoL1SParams(t0, u0, tE, alpha, 0.0, 0.0, s, q)


def chi2_cassan(theta, caustic=None):
    """chi2 in Cassan parameters; pass `caustic` to hold (s, q) fixed and skip recomputing it."""
    _, _, t_in, t_out, s, q = CassanParams(*theta)
    if t_out <= t_in or s <= 0 or q <= 0:
        return np.inf
    return chi2(to_standard(theta, caustic))


def fit(s, q, t_in, t_out_list, n_sigma=20):
    """Grid over (sigma_in, sigma_out, t_out) at fixed (s, q) and t_in, Nelder-Mead the best
    few grid points in the 4 crossing parameters, then free (s, q) for a final 6-param refit."""
    caustic = cassan_caustic(s, q)
    sigmas = np.arange(n_sigma) / n_sigma
    grid = [(si, so, t_in, to) for si, so, to in product(sigmas, sigmas, t_out_list) if si != so]
    grid_chi2 = [chi2_cassan((*g, s, q), caustic) for g in grid]
    starts = [grid[i] for i in np.argsort(grid_chi2)[:5]]

    fixed = [minimize(lambda x: chi2_cassan((*x, s, q), caustic), x0=x0, method="Nelder-Mead",
                      options={"xatol": 1e-4, "fatol": 1e-3, "maxiter": 4000}) for x0 in starts]
    best_fixed = min(fixed, key=lambda r: r.fun)
    print(f"  fixed s={s}, q={q}: grid best chi2={min(grid_chi2):.2f} -> refined {best_fixed.fun:.2f}")

    result = minimize(chi2_cassan, x0=(*best_fixed.x, s, q), method="Nelder-Mead",
                      options={"xatol": 1e-5, "fatol": 1e-4, "maxiter": 20000})
    best = CassanParams(*result.x)
    print(f"  free (s, q): chi2={result.fun:.2f}")
    for label, value in zip(CassanParams._fields + TwoL1SParams._fields[:4], best + to_standard(best)[:4]):
        print(f"    {label} = {value:.5f}")
    return best, result.fun


def log_probability(theta):
    """Gaussian likelihood (-0.5*chi2) under flat priors: sigma_in/sigma_out uniform on the
    (periodic) caustic, t_in < t_out uniform in T_WINDOW, s uniform and log(q) uniform in
    mcmc_fit_2l1s's ranges, derived tE inside its TE_RANGE, and resonant topology only (one
    caustic -- Bond et al.'s is, and it keeps cassan_caustic()'s idx=0 unambiguous)."""
    _, _, t_in, t_out, s, q = CassanParams(*theta)
    if not (T_WINDOW[0] < t_in < t_out < T_WINDOW[1]):
        return -np.inf
    if not (S_RANGE[0] < s < S_RANGE[1]) or q <= 0 or not (LOG_Q_RANGE[0] < np.log(q) < LOG_Q_RANGE[1]):
        return -np.inf
    try:
        if len(caustic_curve(s, q, n_phi=200)) != 1:
            return -np.inf
    except ValueError:  # pieces didn't join up -- right at a topology boundary
        return -np.inf
    standard = to_standard(theta)
    if not (TE_RANGE[0] < standard.tE < TE_RANGE[1]):
        return -np.inf
    c2 = chi2(standard)
    return -0.5 * c2 if np.isfinite(c2) else -np.inf


def run_mcmc(best, nwalkers=32, nsteps=3000, seed=42):
    """emcee from a tight ball around the Nelder-Mead best fit -- explores the posterior of
    the solution already found, not a global search."""
    rng = np.random.default_rng(seed)
    ball = np.array([1e-3, 1e-3, 0.05, 0.05, 1e-3, 1e-2 * best.q])
    p0 = np.array(best) + ball * rng.standard_normal((nwalkers, len(best)))

    with get_context("spawn").Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, len(best), log_probability, pool=pool)
        sampler.run_mcmc(p0, nsteps, progress=True)

    samples = sampler.get_chain(discard=nsteps // 4, thin=15, flat=True)
    samples[:, :2] %= 1  # sigma is periodic
    log_probs = sampler.get_log_prob(discard=nsteps // 4, thin=15, flat=True)
    top = CassanParams(*samples[np.argmax(log_probs)])
    print(f"[mcmc-cassan] {samples.shape[0]} samples, acceptance={np.mean(sampler.acceptance_fraction):.2f}, "
          f"best sample chi2={chi2_cassan(top):.2f}")
    lo, med, hi = np.percentile(samples, [16, 50, 84], axis=0)
    for label, l, m, h in zip(CassanParams._fields, lo, med, hi):
        print(f"[mcmc-cassan] {label} = {m:.5f} +{h - m:.5f} -{m - l:.5f}")

    plot_fit(to_standard(top), tag="_cassan_mcmc")
    save_corner(samples, list(CassanParams._fields), top, str(OUT_DIR / f"{SHORT_NAME}_2l1s_cassan_mcmc_corner.png"))
    np.savez(OUT_DIR / f"{SHORT_NAME}_2l1s_cassan_mcmc_chain.npz", samples=samples, log_probs=log_probs,
             labels=CassanParams._fields)
    return sampler, samples


if __name__ == "__main__":
    # validation: Bond et al. 2004's published solution, standard -> Cassan -> standard
    bond = TwoL1SParams(2848.06, 0.133, 61.5, np.radians(223.8), 0.0, 0.0, 1.120, 0.0039)
    caustic = cassan_caustic(bond.s, bond.q)
    sigma, t = standard_to_cassan(caustic, *bond[:4])
    bond_cassan = CassanParams(sigma[0], sigma[-1], t[0], t[-1], bond.s, bond.q)
    print(f"Bond et al. as Cassan params: {np.round(bond_cassan, 4)}")
    print(f"  chi2 standard={chi2(bond):.2f}, via Cassan round trip={chi2_cassan(bond_cassan):.2f}")

    # fits: only the by-eye entry anchor t_in~2835 is given; t_out is gridded, not assumed.
    # Starting (s, q): Bond's, then session 11's Student-t MCMC geometry (a wrong basin before).
    t_out_list = [2838, 2842, 2846, 2850, 2855]
    fits = []
    for s, q, tag in [(1.120, 0.0039, "_cassan_from_bond_sq"), (1.641, 0.106, "_cassan_from_studentt_sq")]:
        print(f"fit starting from s={s}, q={q}:")
        best, best_chi2 = fit(s, q, t_in=2835.0, t_out_list=t_out_list)
        plot_fit(to_standard(best), tag=tag)
        fits.append((best_chi2, best))

    run_mcmc(min(fits)[1])
