"""Config-driven FSPL vs 2L1S search -- the new pipeline's CLI.

    python3 search.py --config input/O-03-BLG235.toml

Stages (session 18 plan, see CHANGELOG):
  (i)   automatic FSPL fit from the data alone (no parallax, then parallax);
        per-instrument K derived here too
  (ii)  (s, q) grid from event.grid: Cassan (every caustic) + standard
        (n_alpha alphas) in each cell, point source, parallax off; each cell
        is checkpointed to grid_partial.npz, so a killed run resumes
  (iii) distinct local minima of the delta-chi2(s, q) map
  (iv)  refine each minimum with everything free (both u0 signs), MCMC on
        those within delta-chi2 <~ 10 of the best
  (v)   delta-chi2 map plot + BIC verdict (K from 2L1S and from FSPL)
  (vi)  diagnostics on the saved outputs (also --stage=diagnose, nothing refit):
        MCMC convergence + trace plots, Nelder-Mead vs MCMC chi2, fs/fb, fit plot

Outputs go to results/<short_name>/. The grid is cached there as grid.npz
(delete it to recompute). Heavy -- run via slurm/search.sbatch, never on the
login node.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from itertools import combinations, product
from multiprocessing import get_context
from pathlib import Path
from typing import NamedTuple

import emcee
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import median_filter, minimum_filter
from scipy.optimize import minimize

from event import (Event, chi2, error_scale, flux_residuals, load_event, profile_flux, rescale, to_magnification,
                   with_t0_par)
from lc_models import (binary_magnification_vbbl, binary_trajectory, cassan_caustic, cassan_to_standard,
                       caustic_curve, fspl_magnification, sun_earth_projection, trajectory)
from mcmc_fit import flat_chain, save_corner
from zoom_utils import find_zoom_window

N_SIGMA = 8          # Cassan inner grid: entry/exit abscissae per caustic
N_REFINE = 10        # lowest local minima of the (s, q) map that get refined
MODE_DCHI2 = 10.0    # refined minima within this of the best get an MCMC


class FSPL(NamedTuple):
    t0: float
    u0: float
    tE: float
    rho: float
    piE_N: float
    piE_E: float


class Binary(NamedTuple):
    t0: float
    u0: float
    tE: float
    rho: float
    piE_N: float
    piE_E: float
    s: float
    q: float
    alpha: float


def fspl_A(event: Event, p: FSPL) -> list[np.ndarray]:
    return [fspl_magnification(trajectory(i.time, p.t0, p.u0, p.tE, p.piE_N, p.piE_E, i.dsN, i.dsE), p.rho, i.ld)
            for i in event.instruments]


def binary_A(event: Event, p: Binary) -> list[np.ndarray]:
    return [binary_magnification_vbbl(binary_trajectory(i.time, p.t0, p.u0, p.tE, p.alpha, p.piE_N, p.piE_E,
                                                        i.dsN, i.dsE), p.s, p.q, p.rho, i.ld)
            for i in event.instruments]


def chi2_fspl(theta, event: Event) -> float:
    p = FSPL(*theta)
    if p.tE <= 0 or p.rho <= 0:
        return np.inf
    return chi2(event, fspl_A(event, p))


def chi2_binary(theta, event: Event) -> float:
    p = Binary(*theta)
    if p.tE <= 0 or p.rho < 0 or p.s <= 0 or not 0 < p.q <= 1:  # rho = 0: point source
        return np.inf
    return chi2(event, binary_A(event, p))


def fspl_steps(p: FSPL):
    return (0.05 * p.tE, 0.1 * abs(p.u0) + 1e-3, 0.1 * p.tE, 0.5 * p.rho, 0.05, 0.05)


def binary_steps(p: Binary):
    return (0.01 * p.tE, 0.05 * abs(p.u0) + 1e-3, 0.05 * p.tE, 0.2 * p.rho, 0.05, 0.05, 0.01 * p.s, 0.1 * p.q, 0.02)


def nelder_mead(f, x0, step, maxfev=20000):
    """Nelder-Mead from an explicit initial simplex: scipy's default steps 5% of each
    value, ~140 d on t0 ~ 2848 (see CHANGELOG session 17)."""
    x0 = np.asarray(x0, float)
    r = minimize(f, x0, method="Nelder-Mead", options={
        "initial_simplex": np.vstack([x0, x0 + np.diag(step)]), "xatol": 1e-6, "fatol": 1e-3, "maxfev": maxfev})
    return r.x, float(r.fun)


# ---- (i) FSPL ------------------------------------------------------------------------------

def fit_fspl(event: Event):
    """Blind FSPL fit: t0 candidates at each instrument's (median-smoothed) brightest point,
    a coarse u0 x tE grid, Nelder-Mead from the best 3 without parallax (-> `plain`, whose t0
    becomes t0_par), then with parallax from both u0 signs. Returns (plain, fspl, event with
    parallax offsets)."""
    t0s = {float(i.time[np.argmax(median_filter(i.flux, size=5))]) for i in event.instruments}
    grid = [FSPL(t0, u0, tE, 1e-4, 0, 0) for t0, u0, tE in product(t0s, [0.01, 0.03, 0.1, 0.3, 1.0], [3, 10, 30, 100])]
    starts = sorted(grid, key=lambda p: chi2_fspl(p, event))[:3]
    fits = [nelder_mead(lambda x: chi2_fspl((*x, 0, 0), event), p[:4], fspl_steps(p)[:4]) for p in starts]
    x, c2 = min(fits, key=lambda r: r[1])
    plain = FSPL._make([*x, 0.0, 0.0])
    print(f"[fspl] no parallax: chi2={c2:.2f} {plain}")

    event = with_t0_par(event, plain.t0)
    fits = [nelder_mead(lambda x: chi2_fspl(x, event), p, fspl_steps(p))
            for p in (plain, plain._replace(u0=-plain.u0))]
    x, c2 = min(fits, key=lambda r: r[1])
    fspl = FSPL(*x)
    print(f"[fspl] parallax (t0_par={plain.t0:.4f}): chi2={c2:.2f} {fspl}")
    return plain, fspl, event


# ---- (ii) grid -----------------------------------------------------------------------------

def anomaly_times(event: Event, plain: FSPL, n=6):
    """Times of the n largest |FSPL residuals| within 2 tE of t0, >= 1 d apart: data-driven
    candidates for Cassan's t_in/t_out (where caustic crossings could sit)."""
    t = np.concatenate([i.time for i in event.instruments])
    r = np.abs(flux_residuals(event, fspl_A(event, plain)))
    r[np.abs(t - plain.t0) > 2 * plain.tE] = 0
    picks: list[float] = []
    for k in np.argsort(r)[::-1]:
        if all(abs(t[k] - p) >= 1 for p in picks):
            picks.append(float(t[k]))
        if len(picks) == n:
            break
    return picks


# Pool worker state, set once per worker by _init (spawned workers re-import this module)
EVENT: Event
CTX: tuple = ()


def _init(event, *ctx):
    global EVENT, CTX
    EVENT, CTX = event, ctx


def grid_cell(k_sq):
    """Best parallax-free point-source 2L1S at fixed (s, q), over both families:
    Cassan (inner grid over entry/exit abscissa x candidate times, every caustic, Nelder-Mead
    the best 2) and standard (n_alpha alphas from FSPL's t0/u0/tE, short Nelder-Mead each).
    Point source keeps every chi2 the same cost: finite source near a cusp took up to minutes
    per call (sessions 20-21). The returned Binary carries FSPL's rho to start refinement from.
    Takes (cell index, (s, q)), returns (cell index, chi2, Binary) for imap_unordered."""
    k, (s, q) = k_sq
    plain, times, n_alpha = CTX
    rho = plain.rho

    def c2(t0, u0, tE, alpha):
        return chi2_binary((t0, u0, tE, 0.0, 0, 0, s, q, alpha), EVENT)

    results = []
    try:
        caustics = [cassan_caustic(s, q, idx) for idx in range(len(caustic_curve(s, q)))]
    except ValueError:  # pieces didn't join: right at a topology boundary -- standard family only
        caustics = []
    # half-step offset: sigma = 0 and 0.5 are the on-axis cusps, and pairing them sends the
    # trajectory straight down the binary axis through both (session 21)
    sigmas = (np.arange(N_SIGMA) + 0.5) / N_SIGMA
    for caustic in caustics:
        def cassan(x, caustic=caustic):
            if x[3] <= x[2]:
                return np.inf
            return c2(*cassan_to_standard(caustic, *x))
        grid = [(si, so, ti, to) for si, so in product(sigmas, sigmas) if si != so for ti, to in combinations(times, 2)]
        for x0 in [x for x in sorted(grid, key=cassan)[:2] if np.isfinite(cassan(x))]:
            x, _ = nelder_mead(cassan, x0, (0.02, 0.02, 0.01 * plain.tE, 0.01 * plain.tE), maxfev=600)
            t0, u0, tE, alpha = cassan_to_standard(caustic, *x)
            results.append((c2(t0, u0, tE, alpha), Binary(t0, u0, tE, rho, 0, 0, s, q, alpha)))

    for alpha in np.linspace(0, 2 * np.pi, n_alpha, endpoint=False):
        x0 = (plain.t0, abs(plain.u0), plain.tE, alpha)
        x, c = nelder_mead(lambda x: c2(*x), x0, (0.01 * plain.tE, 0.05 * abs(plain.u0) + 1e-3, 0.05 * plain.tE, 0.05),
                           maxfev=300)
        results.append((c, Binary(x[0], x[1], x[2], rho, 0, 0, s, q, x[3])))
    return (k, *min(results, key=lambda r: r[0]))


def axis(spec):
    """[lo, hi, step] -> lo..hi inclusive."""
    lo, hi, step = spec
    return np.linspace(lo, hi, round((hi - lo) / step) + 1)


def run_grid(event: Event, plain: FSPL, path: Path):
    """(ii) delta-chi2(s, q) map, cached at `path`."""
    if path.exists():
        print(f"[grid] loading cached {path}")
        return dict(np.load(path))
    log_s, log_q = axis(event.grid["log_s"]), axis(event.grid["log_q"])
    times = sorted(anomaly_times(event, plain) + list(plain.t0 + plain.tE * np.linspace(-1, 1, 5)))
    print(f"[grid] {log_s.size} x {log_q.size} cells, Cassan times {np.round(times, 2)}")
    cells = [(10 ** a, 10 ** b) for a, b in product(log_s, log_q)]
    chi2_map, theta = np.empty(len(cells)), np.empty((len(cells), len(Binary._fields)))
    done = np.zeros(len(cells), bool)
    partial = path.with_name("grid_partial.npz")  # every finished cell, so a killed run resumes
    if partial.exists():
        d = np.load(partial)
        if np.array_equal(d["log_s"], log_s) and np.array_equal(d["log_q"], log_q):
            chi2_map, theta, done = d["chi2"], d["theta"], d["done"]
            print(f"[grid] resuming: {done.sum()}/{len(cells)} cells from {partial}")
    todo = [(k, c) for k, c in enumerate(cells) if not done[k]]
    with get_context("spawn").Pool(initializer=_init, initargs=(event, plain, times, event.grid["n_alpha"])) as pool:
        for k, c, p in pool.imap_unordered(grid_cell, todo):
            chi2_map[k], theta[k], done[k] = c, p, True
            with open(partial.with_suffix(".tmp"), "wb") as f:  # write-then-rename: a kill can't corrupt it
                np.savez(f, log_s=log_s, log_q=log_q, chi2=chi2_map, theta=theta, done=done)
            partial.with_suffix(".tmp").replace(partial)
            if done.sum() % 50 == 0:
                print(f"[grid] {done.sum()}/{len(cells)}")
    chi2_map, theta = chi2_map.reshape(log_s.size, log_q.size), theta.reshape(log_s.size, log_q.size, -1)
    np.savez(path, log_s=log_s, log_q=log_q, chi2=chi2_map, theta=theta)
    partial.unlink()
    return dict(log_s=log_s, log_q=log_q, chi2=chi2_map, theta=theta)


# ---- (iii)-(iv) minima, refinement, MCMC ------------------------------------------------------

def local_minima(chi2_map, n=N_REFINE):
    """(iii) Indices of the n lowest cells that are minima of their 3x3 neighbourhood."""
    is_min = (chi2_map == minimum_filter(chi2_map, size=3, mode="nearest")) & np.isfinite(chi2_map)
    return np.argwhere(is_min)[np.argsort(chi2_map[is_min])][:n]


def distinct_modes(refined, tol=(0.02, 0.1)):
    """Refined (x, chi2) pairs, sorted by chi2 -> the Binary modes within MODE_DCHI2 of the best,
    one per (log s, log q, sign u0) cluster: a fit within `tol` in (log s, log q) of a kept mode
    with the same u0 sign is a duplicate. A tolerance, not rounding: rounding split s = 0.751
    and 0.748 into two "modes" across a bin edge (session 21)."""
    modes = []
    for x, c2 in refined:
        p = Binary(*x)
        near = lambda m: (np.sign(m.u0) == np.sign(p.u0) and abs(np.log10(m.s / p.s)) < tol[0]
                          and abs(np.log10(m.q / p.q)) < tol[1])
        if c2 - refined[0][1] <= MODE_DCHI2 and not any(near(m) for m in modes):
            modes.append(p)
    return modes


def refine_starts(theta):
    """Everything free, parallax on (from zero), from both (u0, alpha) and its mirror (-u0, -alpha):
    identical without parallax, distinct with it (ecliptic degeneracy)."""
    p = Binary(*theta)._replace(piE_N=0.0, piE_E=0.0)
    return [p, p._replace(u0=-p.u0, alpha=-p.alpha)]


def refine(x0: Binary):
    """One Nelder-Mead per Pool task, so both starts of a minimum run on separate cores."""
    return nelder_mead(lambda x: chi2_binary(x, EVENT), x0, binary_steps(x0), maxfev=5000)


def log_prob(theta):
    """Flat priors (chi2_binary's bounds, |piE| < 5), Gaussian likelihood."""
    if np.hypot(theta[4], theta[5]) > 5:
        return -np.inf
    c2 = chi2_binary(theta, EVENT)
    return -0.5 * c2 if np.isfinite(c2) else -np.inf


def run_mcmc(pool, best: Binary, tag: str, nwalkers=32, nsteps=12000, seed=42):
    """emcee around a refined minimum. nsteps: the Cassan chain had tau ~ 230 steps (session 19),
    so 3000 was only ~13 tau; 12000 aims at the ~50 tau emcee recommends. `pool` is shared
    by every mode at once (see __main__); save_mcmc() does the rest, in the main thread."""
    rng = np.random.default_rng(seed)
    p0 = np.array(best) + 0.01 * np.array(binary_steps(best)) * rng.standard_normal((nwalkers, len(best)))
    sampler = emcee.EnsembleSampler(nwalkers, len(best), log_prob, pool=pool)
    sampler.run_mcmc(p0, nsteps, progress=True, progress_kwargs={"desc": tag})
    return sampler


def save_mcmc(sampler, out: Path, tag: str):
    nsteps = sampler.iteration
    tau = sampler.get_autocorr_time(quiet=True)  # warns (not raises) when nsteps < 50 tau
    print(f"[mcmc {tag}] acceptance={np.mean(sampler.acceptance_fraction):.2f}, "
          f"nsteps/max(tau)={nsteps / np.max(tau):.1f} (want > 50), tau={np.round(tau, 1)}")
    samples, log_probs = flat_chain(sampler, discard=nsteps // 4)
    lo, med, hi = np.percentile(samples, [16, 50, 84], axis=0)
    for label, l, m, h in zip(Binary._fields, lo, med, hi):
        print(f"[mcmc {tag}] {label} = {m:.6g} +{h - m:.3g} -{m - l:.3g}")
    top = samples[np.argmax(log_probs)]
    save_corner(samples, list(Binary._fields), top, str(out / f"mcmc_{tag}_corner.png"))
    chain = sampler.get_chain()
    assert chain is not None  # emcee's getters are inferred Optional, see flat_chain()
    np.savez(out / f"mcmc_{tag}_chain.npz", samples=samples, log_probs=log_probs, chain=chain,
             labels=Binary._fields)


# ---- plots -----------------------------------------------------------------------------------

def plot_raw(event: Event, plain: FSPL, path: Path):
    """Data only, every instrument inverted onto magnification A via its own fs/fb profiled
    at the parallax-free FSPL (no common zero point exists otherwise: "dia" fluxes, FTN's
    offset mags). Full baseline + peak zoom from the longest-baseline instrument; log A only
    for high-magnification events (on a linear-A event, noisy A ~ 0 baseline points stretch it)."""
    calibrated = [to_magnification(i, A_i) for i, A_i in zip(event.instruments, fspl_A(event, plain))]
    k = max(range(len(event.instruments)), key=lambda k: np.ptp(event.instruments[k].time))
    zoom = find_zoom_window(event.instruments[k].time, *calibrated[k])
    yscale = "log" if max(np.nanmax(A) for A, _ in calibrated) > 50 else "linear"
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, xlim in zip(axes, (None, zoom)):
        for i, (A, A_err) in zip(event.instruments, calibrated):
            ax.errorbar(i.time, A, A_err, fmt=".", ms=3, elinewidth=0.5, label=f"{i.name} ({i.band})")
        ax.set(xlabel="HJD - 2450000", ylabel="magnification A (calibrated at FSPL)", xlim=xlim, yscale=yscale)
    axes[0].legend(fontsize=8)
    axes[0].set_title(f"{event.short_name}: full baseline")
    axes[1].set_title("peak")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"saved {path}")


# ---- (vi) diagnostics: read the saved outputs, nothing refit ----------------------------------

def on_grid(event: Event, inst, t, t0_par):
    """One-instrument Event at times t (inst's kind/ld), parallax offsets at t0_par: model curves."""
    dsN, dsE = sun_earth_projection(t, event.coords, t0_par)
    return event._replace(instruments=[inst._replace(time=t, dsN=dsN, dsE=dsE)])


def plot_fit(event: Event, fspl: FSPL, best: Binary, t0_par: float, path: Path):
    """Data in A via each instrument's fs/fb at the best 2L1S, 2L1S and FSPL curves, then one raw
    residual row per model (data in that model's own calibration minus it). Left: peak; right:
    the anomaly: centred on the point where 2L1S lowers chi2 most, half-width 1.5x the median
    distance of the points gaining > 10% as much (median: a few wing points also gain and would
    stretch a min-max span) -- not where the curves differ, which is everywhere.
    FSPL's curve is mapped onto the 2L1S scale through the longest-baseline instrument: with
    blending, that mapping differs per instrument, so other instruments can sit off it."""
    A_f, A_b = fspl_A(event, fspl), binary_A(event, best)
    k = max(range(len(event.instruments)), key=lambda k: np.ptp(event.instruments[k].time))
    ref = event.instruments[k]
    cf, cb = profile_flux(ref, A_f[k])[0], profile_flux(ref, A_b[k])[0]

    def fspl_on_2l1s(A):  # FSPL's predicted flux for `ref`, in the 2L1S calibration's A units
        return (cf[0] * A + cf[1] - cb[1]) / cb[0] if ref.kind == "mag" else 1 + cf[0] * (A - 1) / cb[0]

    t = np.concatenate([i.time for i in event.instruments])
    gain = flux_residuals(event, A_f) ** 2 - flux_residuals(event, A_b) ** 2
    tc = t[np.argmax(gain)]
    half = max(1.5 * np.median(np.abs(t[gain > 0.1 * gain.max()] - tc)), 0.05)
    windows = (find_zoom_window(ref.time, *to_magnification(ref, A_b[k])), (tc - half, tc + half))
    yscale = "log" if max(A.max() for A in A_b) > 50 else "linear"
    fig, axes = plt.subplots(3, 2, figsize=(15, 11), height_ratios=[4, 1.5, 1.5], sharex="col")
    for col, xlim in enumerate(windows):
        for n, (i, Af, Ab) in enumerate(zip(event.instruments, A_f, A_b)):
            m = (i.time >= xlim[0]) & (i.time <= xlim[1])
            for row, Am in ((0, Ab), (1, Af), (2, Ab)):
                A, err = to_magnification(i, Am)
                y = A if row == 0 else A - Am
                axes[row, col].errorbar(i.time[m], y[m], err[m], fmt=".", ms=3, elinewidth=0.5, color=f"C{n}",
                                        label=f"{i.name} ({i.band})")
        tg = np.linspace(*xlim, 2000)
        g = on_grid(event, ref, tg, t0_par)
        axes[0, col].plot(tg, binary_A(g, best)[0], "k-", lw=1, label="2L1S")
        axes[0, col].plot(tg, fspl_on_2l1s(fspl_A(g, fspl)[0]), "k--", lw=1, label=f"FSPL ({ref.name} calibration)")
        axes[0, col].set(xlim=xlim, yscale=yscale, ylabel="magnification A (2L1S calibration)")
        for row, name in ((1, "FSPL"), (2, "2L1S")):
            axes[row, col].axhline(0, color="k", lw=0.5)
            axes[row, col].set_ylabel(f"data - {name} (A)")
        axes[2, col].set_xlabel("HJD - 2450000")
        axes[2, col].ticklabel_format(axis="x", useOffset=False)
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].set_title(f"{event.short_name}: peak")
    axes[0, 1].set_title("anomaly (where 2L1S gains most chi2)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"saved {path}")


def plot_trace(chain, path: Path):
    """Every walker's path per parameter: burn-in, stuck walkers, drift."""
    fig, axes = plt.subplots(chain.shape[2], 1, figsize=(10, 1.4 * chain.shape[2]), sharex=True)
    for n, (ax, label) in enumerate(zip(axes, Binary._fields)):
        ax.plot(chain[:, :, n], lw=0.3, alpha=0.4)
        ax.set_ylabel(label)
    axes[-1].set_xlabel("step")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"saved {path}")


def diagnose(event: Event, fspl: FSPL, t0_par: float, out: Path):
    """(vi) Session 20's done-checks: per MCMC mode nsteps/tau > 50, acceptance 0.2-0.5, its best
    sample within ~1 of the Nelder-Mead chi2 it started from, trace plot; per instrument fs/fb
    at FSPL and the best 2L1S (fb < 0 flagged: expected for FTN's offset scale, suspicious
    elsewhere); plot_fit() at the best refined 2L1S."""
    ok = lambda c: "ok" if c else "FAIL"
    refined = np.load(out / "refined.npz")
    for f in sorted(out.glob("mcmc_mode*_chain.npz")):
        tag = f.stem.removeprefix("mcmc_").removesuffix("_chain")
        d = np.load(f)
        chain = d["chain"]  # (nsteps, nwalkers, ndim)
        tau = emcee.autocorr.integrated_time(chain, quiet=True)
        acc = np.mean(np.any(np.diff(chain, axis=0) != 0, axis=2))  # a walker that moved accepted
        start = chain[0].mean(axis=0)  # walkers start in a 1%-of-step ball around their refined minimum
        nm = refined["chi2"][np.argmin(np.linalg.norm((refined["theta"] - start) / binary_steps(Binary(*start)), axis=1))]
        best = -2 * d["log_probs"].max()
        print(f"[diag {tag}] nsteps/tau={len(chain) / tau.max():.1f} {ok(len(chain) / tau.max() > 50)}, "
              f"acceptance={acc:.2f} {ok(0.2 <= acc <= 0.5)}, chi2 Nelder-Mead={nm:.2f} vs MCMC best={best:.2f} "
              f"{ok(abs(nm - best) <= 1)}")
        plot_trace(chain, out / f"mcmc_{tag}_trace.png")
    best = Binary(*refined["theta"][0])
    for name, A in (("FSPL", fspl_A(event, fspl)), ("2L1S", binary_A(event, best))):
        for i, A_i in zip(event.instruments, A):
            c = profile_flux(i, A_i)[0]
            fb = f", fb={c[1]:.4g}{' (< 0)' if c[1] < 0 else ''}" if i.kind == "mag" else ""
            print(f"[diag {name}] {i.name}: fs={c[0]:.4g}{fb}")
    plot_fit(event, fspl, best, t0_par, out / "fit_lc.png")


# ---- (v) map + BIC ----------------------------------------------------------------------------

def plot_map(grid, minima, path: Path):
    d = np.ma.masked_invalid(grid["chi2"] - np.nanmin(grid["chi2"][np.isfinite(grid["chi2"])]))
    fig, ax = plt.subplots(figsize=(7, 6))
    mesh = ax.pcolormesh(grid["log_s"], grid["log_q"], np.log10(1 + d).T, shading="nearest", cmap="viridis_r")
    fig.colorbar(mesh, label=r"$\log_{10}(1 + \Delta\chi^2)$")
    ax.plot(grid["log_s"][minima[:, 0]], grid["log_q"][minima[:, 1]], "r+", ms=10, label="refined minima")
    ax.set(xlabel=r"$\log_{10} s$", ylabel=r"$\log_{10} q$", title=r"2L1S $\Delta\chi^2(s, q)$, parallax off")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"saved {path}")


def bic(c2, n_model, event: Event):
    """chi2 + k ln N, with k counting the profiled flux parameters too."""
    n = sum(i.time.size for i in event.instruments)
    k = n_model + sum(2 if i.kind == "mag" else 1 for i in event.instruments)
    return c2 + k * np.log(n)


def out_dir(event: Event) -> Path:
    path = Path("results") / event.short_name
    path.mkdir(parents=True, exist_ok=True)
    return path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Config-driven FSPL vs 2L1S search.")
    parser.add_argument("--config", required=True, help="path to the event's TOML")
    parser.add_argument("--stage", choices=["raw", "search", "diagnose"], default="search",
                        help="raw: FSPL fit + calibrated light curve only (~minutes; srun, not the login node). "
                             "diagnose: FSPL fit + (vi) on a finished search's saved outputs")
    args = parser.parse_args()
    base = load_event(args.config)
    out = out_dir(base)
    plain, fspl, base = fit_fspl(base)  # base now carries the parallax offsets at t0_par
    plot_raw(base, plain, out / "raw_lc.png")
    if args.stage == "raw":
        raise SystemExit
    k_fspl = error_scale(base, fspl_A(base, fspl))
    event = rescale(base, k_fspl)  # the search runs on FSPL-rescaled errors
    print(f"[fspl] error rescaling k = {np.round(k_fspl, 4)}")
    if args.stage == "diagnose":
        diagnose(event, fspl, plain.t0, out)
        raise SystemExit

    grid = run_grid(event, plain, out / "grid.npz")
    minima = local_minima(grid["chi2"])
    plot_map(grid, minima, out / "delta_chi2_map.png")

    with get_context("spawn").Pool(initializer=_init, initargs=(event,)) as pool:
        refined = sorted(pool.map(refine, [x0 for i, j in minima for x0 in refine_starts(grid["theta"][i, j])]),
                         key=lambda r: r[1])
        for x, c2 in refined:
            print(f"[refine] chi2={c2:.2f} {Binary(*x)}")
        np.savez(out / "refined.npz", theta=[x for x, _ in refined], chi2=[c for _, c in refined])
        # every mode at once, one thread each, on the one pool: emcee maps only nwalkers/2 walkers at
        # a time and each step waits for its slowest chi2, so a lone mode leaves most cores idle
        modes = distinct_modes(refined)
        with ThreadPoolExecutor(len(modes)) as threads:
            samplers = list(threads.map(lambda kp: run_mcmc(pool, kp[1], f"mode{kp[0]}"), enumerate(modes)))
    for k, sampler in enumerate(samplers):  # plotting stays in the main thread: pyplot isn't thread-safe
        save_mcmc(sampler, out, f"mode{k}")

    best = Binary(*refined[0][0])
    k_2l1s = error_scale(base, binary_A(base, best))
    print(f"[2l1s] error rescaling k = {np.round(k_2l1s, 4)}")
    # ponytail: both models are re-scored at their fixed best fits under each K, not refit --
    # a K change only reweights instruments, shifting the best fit slightly
    for name, k in (("K at 2L1S", k_2l1s), ("K at FSPL", k_fspl)):
        ev = rescale(base, k)
        bic_f, bic_b = bic(chi2_fspl(fspl, ev), 6, ev), bic(chi2_binary(best, ev), 9, ev)
        print(f"[bic {name}] FSPL {bic_f:.2f}, 2L1S {bic_b:.2f}, Delta BIC (FSPL - 2L1S) = {bic_f - bic_b:.2f}")
    diagnose(event, fspl, plain.t0, out)
