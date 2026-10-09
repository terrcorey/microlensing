"""Config-driven FSPL vs 2L1S search -- the new pipeline's CLI.

    python3 search.py --config input/O-03-BLG235.toml

Stages (session 18 plan, see CHANGELOG):
  (i)   automatic FSPL fit from the data alone (no parallax, then a (piE_N, piE_E) grid x both
        u0 signs + refinement: the parallax-only model 2L1S is judged against, as thorough as the
        planet search); per-instrument K derived here too
  (ii)  (s, q) grid from event.grid: Cassan (every caustic) + standard
        (n_alpha alphas) in each cell, point source, parallax off; each cell
        is checkpointed to grid_partial.npz, so a killed run resumes
  (iii) distinct local minima of the delta-chi2(s, q) map
  (iv)  refine each minimum with everything free (both u0 signs, both start families: Cassan and
        standard), MCMC on those within delta-chi2 <~ 10 of the best, crossing and near-miss
        solutions kept apart (session 23)
  (v)   delta-chi2 map plot + BIC verdict (K from 2L1S and from FSPL)
  (vi)  diagnostics on the saved outputs (also --stage=diagnose, nothing refit):
        MCMC convergence + trace plots, Nelder-Mead vs MCMC chi2, fs/fb, fit plot

Outputs go to results/<short_name>/. The grid is cached there as grid.npz
(delete it to recompute). Heavy -- run via slurm/search.sbatch, never on the
login node.
"""
import argparse
import hashlib
import subprocess
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from datetime import datetime
from itertools import combinations, product
from multiprocessing import get_context
from pathlib import Path
from typing import NamedTuple

import emcee
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.ndimage import median_filter, minimum_filter

from event import (Event, chi2, error_scale, flux_residuals, load_event, profile_flux, rescale, to_magnification,
                   with_t0_par)
from lc_models import (ZERO_POINT_MAG, binary_magnification_vbbl, binary_trajectory, cassan_caustic, cassan_to_standard,
                       caustic_curve, fspl_magnification, lens_position, sun_earth_projection, trajectory)
from mcmc_fit import nelder_mead, polish, save_corner
from zoom_utils import find_zoom_window

FAMILIES = ("Cassan", "standard")  # grid_cell's two start families, refined separately (session 23)
# Cassan inner grid per caustic: N_SIGMA entry/exit abscissae, best N_POLISH coarse starts get a Nelder-Mead;
# STD_MAXFEV per standard-family alpha start. 16/5/600 cost 2.8x per cell and changed no real grid cell on
# O-03-BLG235 (Bond's (s, q) sits between cells: the log q spacing, not the inner search, loses it; session 22)
# 8/2/300 refined O-03-BLG235's cell next to Bond (log q -2.4) to 1447; 16/5/600 found 1283 there and refined
# into Bond's basin (1176) -- denser inner search only pays together with a finer log q step (session 22)
N_SIGMA = 16
N_POLISH = 5
STD_MAXFEV = 600
N_REFINE = 10        # lowest local minima of the (s, q) map that get refined (30 found nothing new on O-03-BLG235, session 22)
MODE_DCHI2 = 10.0    # refined minima within this of the best get an MCMC
# Finite-source floor: rho -> 0 collapsed FSPL to 1e-15..1e-24 on O-05-BLG169-noOGLE (session 22) and
# every refinement then ran VBBL BinaryMag2 at that rho -- the heavy tail / segfault regime. Physical
# rho = theta_* / theta_E is >~ 1e-4 for bulge sources; rho = 0 (exact point source) stays allowed in 2L1S.
RHO_MIN = 1e-5
PIE_MAX = 5.0        # |piE| bound for every chi2 (was MCMC-only: a refined |piE| > 5 froze its walkers)


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
    # isfinite first: VBBL loops forever on a NaN/inf input, and NaN slips past every `<=` test
    if not np.all(np.isfinite(theta)) or p.tE <= 0 or p.rho < RHO_MIN or np.hypot(p.piE_N, p.piE_E) > PIE_MAX:
        return np.inf
    return chi2(event, fspl_A(event, p))


def in_range(x, spec):
    """x > 0 with log10(x) inside a config [lo, hi, step] axis."""
    return x > 0 and spec[0] - 1e-9 <= np.log10(x) <= spec[1] + 1e-9  # slack: 10**lo round-trips inexactly


def chi2_binary(theta, event: Event) -> float:
    """(s, q) confined to the config's grid box: q -> 0 is FSPL again, and walkers drifting
    down it (q ~ 1e-6) burned MCMC time at acceptance ~0.1 (session 21)."""
    p = Binary(*theta)
    if (not np.all(np.isfinite(theta)) or p.tE <= 0 or not (p.rho == 0 or p.rho >= RHO_MIN)  # rho = 0: point source
            or np.hypot(p.piE_N, p.piE_E) > PIE_MAX or not in_range(p.s, event.grid["log_s"])
            or not in_range(p.q, event.grid["log_q"]) or p.q > 1):
        return np.inf
    return chi2(event, binary_A(event, p))


def fspl_steps(p: FSPL):
    return (0.05 * p.tE, 0.1 * abs(p.u0) + 1e-3, 0.1 * p.tE, 0.5 * p.rho, 0.05, 0.05)


def binary_steps(p: Binary):
    return (0.01 * p.tE, 0.05 * abs(p.u0) + 1e-3, 0.05 * p.tE, 0.2 * p.rho + RHO_MIN, 0.05, 0.05, 0.01 * p.s, 0.1 * p.q, 0.02)


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
    assert np.isfinite(c2), "every FSPL start scored inf: check each instrument's data/kind (fs <= 0?)"
    plain = FSPL._make([*x, 0.0, 0.0])
    print(f"[fspl] no parallax: chi2={c2:.2f} {plain}")

    event = with_t0_par(event, plain.t0)
    fits = [nelder_mead(lambda x: chi2_fspl(x, event), p, fspl_steps(p))
            for p in (plain, plain._replace(u0=-plain.u0))]
    x, c2 = min(fits, key=lambda r: r[1])
    fspl = FSPL(*x)
    print(f"[fspl] parallax (t0_par={plain.t0:.4f}): chi2={c2:.2f} {fspl}")
    return plain, fspl, event


PIE_STEP = 0.25  # parallax-only grid spacing in piE_N, piE_E (1-sigma on piE_E ~0.3 on O-05-BLG169)


def fspl_cell(k_pie):
    """Pool task: FSPL at fixed (piE_N, piE_E), (t0, u0, tE, rho) fitted from `plain`, once per u0 sign.
    Returns (cell index, [(chi2 per instrument, x) for u0 > 0, u0 < 0]) -- per instrument so the map
    can be re-expressed at any error rescaling K without refitting."""
    k, (pN, pE) = k_pie
    plain, = CTX
    sizes = np.cumsum([i.time.size for i in EVENT.instruments])[:-1]
    out = []
    for sign in (1, -1):
        x, c = nelder_mead(lambda x: chi2_fspl((*x, pN, pE), EVENT), (plain.t0, sign * abs(plain.u0), plain.tE, plain.rho),
                           fspl_steps(plain)[:4], maxfev=1000)
        r = flux_residuals(EVENT, fspl_A(EVENT, FSPL(*x, pN, pE)))
        out.append(([float(np.sum(ri ** 2)) for ri in np.split(r, sizes)] if np.isfinite(c) else
                    [np.inf] * len(EVENT.instruments), x))
    return k, out


def refine_fspl(x0):
    """Pool task: FSPL with everything free (parallax too) from x0."""
    return nelder_mead(lambda x: chi2_fspl(x, EVENT), x0, fspl_steps(FSPL(*x0)))


def parallax_search(event: Event, plain: FSPL, quick: FSPL, ex):
    """(i) Parallax-only model, searched as thoroughly as the planet: FSPL on a (piE_N, piE_E) grid
    inside |piE| <= PIE_MAX for each u0 sign (ex runs fspl_cell), the 3 lowest local minima refined
    with everything free; `quick` (fit_fspl's two-start fit) competes too, so this can't do worse.
    `event` carries the parallax offsets, errors not yet rescaled. Returns (fspl, pgrid for the map)."""
    pie = np.arange(-PIE_MAX, PIE_MAX + PIE_STEP / 2, PIE_STEP)
    cells = [(k, (a, b)) for k, (a, b) in enumerate(product(pie, pie)) if np.hypot(a, b) <= PIE_MAX]
    n_inst = len(event.instruments)
    c2i, theta = np.full((2, pie.size ** 2, n_inst), np.inf), np.full((2, pie.size ** 2, 4), np.nan)
    for k, out in ex.map(fspl_cell, cells):
        for sgn, (c, x) in enumerate(out):
            c2i[sgn, k], theta[sgn, k] = c, x
    c2i, theta = c2i.reshape(2, pie.size, pie.size, n_inst), theta.reshape(2, pie.size, pie.size, 4)
    total = c2i.sum(axis=-1)
    minima = sorted(((total[g][i, j], g, i, j) for g in range(2) for i, j in local_minima(total[g], n=pie.size ** 2)))
    starts = [FSPL(*theta[g, i, j], pie[i], pie[j]) for _, g, i, j in minima[:3]]
    fits = list(ex.map(refine_fspl, starts)) + [(np.array(quick), chi2_fspl(quick, event))]
    x, c2 = min(fits, key=lambda r: r[1])
    fspl = FSPL(*x)
    on_bound = np.hypot(fspl.piE_N, fspl.piE_E) > PIE_MAX - PIE_STEP / 2
    print(f"[fspl] parallax grid ({len(cells)} cells x 2 u0 signs, {len(minima)} local minima): chi2={c2:.2f} {fspl}"
          + (f"\n[fspl] WARNING: best |piE| = {np.hypot(fspl.piE_N, fspl.piE_E):.2f} sits on the |piE| <= {PIE_MAX} bound"
             if on_bound else ""))
    return fspl, dict(pie=pie, c2i=c2i, on_bound=on_bound)


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
    Takes (cell index, (s, q)), returns (cell index, [chi2], [Binary]), one per FAMILIES entry: both
    get refined at every local minimum (session 23: crossing vs near miss)."""
    k, (s, q) = k_sq
    plain, times, n_alpha = CTX
    rho = plain.rho

    def c2(t0, u0, tE, alpha):
        return chi2_binary((t0, u0, tE, 0.0, 0, 0, s, q, alpha), EVENT)

    cassan_r, std_r = [], []  # (chi2, Binary)
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
        for x0 in [x for x in sorted(grid, key=cassan)[:N_POLISH] if np.isfinite(cassan(x))]:
            x, _ = nelder_mead(cassan, x0, (0.02, 0.02, 0.01 * plain.tE, 0.01 * plain.tE), maxfev=600)
            t0, u0, tE, alpha = cassan_to_standard(caustic, *x)
            cassan_r.append((c2(t0, u0, tE, alpha), Binary(t0, u0, tE, rho, 0, 0, s, q, alpha)))

    for alpha in np.linspace(0, 2 * np.pi, n_alpha, endpoint=False):
        x0 = (plain.t0, abs(plain.u0), plain.tE, alpha)
        if not np.isfinite(c2(*x0)):  # an all-inf simplex never meets fatol: burns maxfev for nothing
            continue
        x, c = nelder_mead(lambda x: c2(*x), x0, (0.01 * plain.tE, 0.05 * abs(plain.u0) + 1e-3, 0.05 * plain.tE, 0.05),
                           maxfev=STD_MAXFEV)
        std_r.append((c, Binary(x[0], x[1], x[2], rho, 0, 0, s, q, x[3])))
    none = (np.inf, Binary(plain.t0, plain.u0, plain.tE, rho, 0, 0, s, q, 0.0))
    best = [min(r, key=lambda r: r[0], default=none) for r in (cassan_r, std_r)]
    return k, [c for c, _ in best], [p for _, p in best]


def axis(spec):
    """[lo, hi, step] -> lo..hi inclusive."""
    lo, hi, step = spec
    return np.linspace(lo, hi, round((hi - lo) / step) + 1)


def cached(path: Path, key: str):
    """np.load(path) if it exists and was written under this config's `key`, else None."""
    if path.exists():
        d = np.load(path)
        if "key" in d and str(d["key"]) == key:
            return dict(d)
        print(f"[cache] {path} is from another config: recomputing")
    return None


def pool(*initargs):
    """Spawned worker pool. ProcessPoolExecutor, not multiprocessing.Pool: when a worker dies
    (VBBL segfault), Pool respawns it, loses its task and waits forever (sessions 21-22, two
    jobs hung ~40 h); the executor raises BrokenProcessPool instead."""
    return ProcessPoolExecutor(mp_context=get_context("spawn"), initializer=_init, initargs=initargs)


def run_grid(event: Event, plain: FSPL, path: Path, key: str, k_fspl):
    """(ii) delta-chi2(s, q) map, cached at `path`."""
    # code-side knobs the TOML key can't see, and the K the cells were scored at (a new FSPL moves it)
    settings = np.array([N_SIGMA, N_POLISH, STD_MAXFEV, *k_fspl])
    fresh = lambda d: (d is not None and "plain" in d and np.allclose(d["plain"], plain)  # cells start from plain
                       and "settings" in d and d["settings"].shape == settings.shape
                       and np.allclose(d["settings"], settings) and "fam_theta" in d)  # pre-session-23: one family
    if fresh(d := cached(path, key)):
        print(f"[grid] loading cached {path}")
        return d
    log_s, log_q = axis(event.grid["log_s"]), axis(event.grid["log_q"])
    times = sorted(anomaly_times(event, plain) + list(plain.t0 + plain.tE * np.linspace(-1, 1, 5)))
    print(f"[grid] {log_s.size} x {log_q.size} cells, Cassan times {np.round(times, 2)}")
    cells = [(10 ** a, 10 ** b) for a, b in product(log_s, log_q)]
    chi2_map, theta = np.empty((len(cells), 2)), np.empty((len(cells), 2, len(Binary._fields)))
    done = np.zeros(len(cells), bool)
    partial = path.with_name("grid_partial.npz")  # every finished cell, so a killed run resumes
    if fresh(d := cached(partial, key)):
        chi2_map, theta, done = d["fam_chi2"], d["fam_theta"], d["done"]
        print(f"[grid] resuming: {done.sum()}/{len(cells)} cells from {partial}")
    todo = [(k, c) for k, c in enumerate(cells) if not done[k]]
    path.with_name("refined.npz").unlink(missing_ok=True)  # refined from the old grid's minima: stale
    with pool(event, plain, times, event.grid["n_alpha"]) as ex:
        futures = [ex.submit(grid_cell, kc) for kc in todo]
        for fut in as_completed(futures):
            try:
                k, c, p = fut.result()
            except BaseException:  # else __exit__ waits for every queued cell, then discards them
                ex.shutdown(cancel_futures=True)
                raise
            chi2_map[k], theta[k], done[k] = c, p, True
            with open(partial.with_suffix(".tmp"), "wb") as f:  # write-then-rename: a kill can't corrupt it
                np.savez(f, log_s=log_s, log_q=log_q, fam_chi2=chi2_map, fam_theta=theta, done=done, key=key,
                         plain=plain, settings=settings)
            partial.with_suffix(".tmp").replace(partial)
            if done.sum() % 50 == 0:
                print(f"[grid] {done.sum()}/{len(cells)}")
    d = dict(log_s=log_s, log_q=log_q, chi2=chi2_map.min(1).reshape(log_s.size, log_q.size),
             fam_chi2=chi2_map.reshape(log_s.size, log_q.size, 2), fam_theta=theta.reshape(log_s.size, log_q.size, 2, -1))
    np.savez(path, **d, key=key, plain=plain, settings=settings)
    partial.unlink(missing_ok=True)
    return d


# ---- (iii)-(iv) minima, refinement, MCMC ------------------------------------------------------

def local_minima(chi2_map, n=N_REFINE):
    """(iii) Indices of the n lowest cells that are minima of their 3x3 neighbourhood."""
    is_min = (chi2_map == minimum_filter(chi2_map, size=3, mode="nearest")) & np.isfinite(chi2_map)
    return np.argwhere(is_min)[np.argsort(chi2_map[is_min])][:n]


def distinct_modes(refined, crosses, tol=(0.02, 0.1)):
    """Refined (x, chi2, ...) sorted by chi2, and whether each crosses a caustic -> the Binary modes
    within MODE_DCHI2 of the best, one per (log s, log q, sign u0, crossing) cluster: a fit within
    `tol` in (log s, log q) of a kept mode with the same u0 sign and class is a duplicate. A tolerance,
    not rounding: rounding split s = 0.751 and 0.748 into two "modes" across a bin edge (session 21).
    The class keeps a near miss within MODE_DCHI2 of a crossing at the same (s, q) (session 23)."""
    modes, kept = [], []
    for (x, c2, *_), cross in zip(refined, crosses):
        p = Binary(*x)
        near = lambda m, mc: (mc == cross and np.sign(m.u0) == np.sign(p.u0) and abs(np.log10(m.s / p.s)) < tol[0]
                              and abs(np.log10(m.q / p.q)) < tol[1])
        if c2 - refined[0][1] <= MODE_DCHI2 and not any(near(m, mc) for m, mc in zip(modes, kept)):
            modes.append(p)
            kept.append(cross)
    return modes


CLOSE = 0.05  # thetaE: caustic_crossing() measures near misses out to this, farther ones print "> CLOSE/rho"


def track(event: Event, t0_par: float):
    """caustic_crossing()'s polyline: every 0.1 d from the first to the last observation, with the
    parallax offsets (one astropy query per event, not per solution; a straight 0.1 d chord errs by
    ~4e-7 piE thetaE), plus the sorted observation times."""
    t_obs = np.sort(np.concatenate([i.time for i in event.instruments]))
    t = np.arange(t_obs[0], t_obs[-1] + 0.1, 0.1)
    return (t, *sun_earth_projection(t, event.coords, t0_par), t_obs)


def caustic_crossing(trk, p: Binary):
    """Crossing vs near miss (session 23): does the source centre's parallax-curved trajectory cross a
    caustic between the first and last observation? -> (crosses, closest approach / rho, unobserved:
    no data point within max(2 rho tE, 0.05 d) of any crossing), or None if the caustic can't be
    built (pieces don't join: right at a topology boundary)."""
    try:
        caustics = caustic_curve(p.s, p.q)
    except ValueError:
        return None
    t, dsN, dsE, t_obs = trk
    z = binary_trajectory(t, p.t0, p.u0, p.tE, p.alpha, p.piE_N, p.piE_E, dsN, dsE)
    cross = lambda u, v: (np.conj(u) * v).imag
    t_x, d = [], np.inf
    for c in caustics:
        c = np.append(c, c[0])
        a, b = z[:-1], z[1:]  # only the segments whose bounding box comes within CLOSE of this caustic's
        k = np.flatnonzero((np.minimum(a.real, b.real) <= c.real.max() + CLOSE)
                           & (np.maximum(a.real, b.real) >= c.real.min() - CLOSE)
                           & (np.minimum(a.imag, b.imag) <= c.imag.max() + CLOSE)
                           & (np.maximum(a.imag, b.imag) >= c.imag.min() - CLOSE))
        if not k.size:
            continue
        a, r, e, f = a[k, None], (b - a)[k, None], c[None, :-1], np.diff(c)[None]  # trajectory x caustic segments
        with np.errstate(divide="ignore", invalid="ignore"):  # parallel segments: never a hit
            lam, mu = cross(e - a, f) / cross(r, f), cross(e - a, r) / cross(r, f)
        hit = (lam >= 0) & (lam <= 1) & (mu >= 0) & (mu <= 1)
        t_x += list((t[k, None] + lam * 0.1)[hit])
        w = np.clip((np.conj(r) * (c[None] - a)).real / np.abs(r) ** 2, 0, 1)  # caustic vertex -> nearest segment point
        d = min(d, np.abs(c[None] - a - w * r).min())
    seen = any(np.abs(t_obs - tx).min() < max(2 * p.rho * p.tE, 0.05) for tx in t_x)
    with np.errstate(divide="ignore"):
        return bool(t_x), 0.0 if t_x else d / p.rho, bool(t_x) and not seen  # d: nearest caustic *vertex*


def refine_starts(theta):
    """Everything free, parallax on (from zero), from both (u0, alpha) and its mirror (-u0, -alpha):
    identical without parallax, distinct with it (ecliptic degeneracy)."""
    p = Binary(*theta)._replace(piE_N=0.0, piE_E=0.0)
    return [p, p._replace(u0=-p.u0, alpha=-p.alpha)]


def refine(x0: Binary):
    """One restarted Nelder-Mead per Pool task, so both starts of a minimum run on separate cores."""
    return polish(lambda x: chi2_binary(x, EVENT), x0, lambda x: binary_steps(Binary(*x)))


def to_mcmc(p: Binary):
    """MCMC coordinates (scalars or arrays): t_eff = u0 tE and t_star = rho tE, what a high-magnification
    peak measures (corr(u0, tE) = -0.94 on O-05-BLG169, session 23), and log s, log q (q spanned two
    decades). Flat in log s, log q = log-uniform priors on s, q; log_prob undoes t_eff/t_star's Jacobian."""
    return np.stack([p.t0, p.u0 * p.tE, p.tE, p.rho * p.tE, p.piE_N, p.piE_E, np.log10(p.s), np.log10(p.q),
                     p.alpha], -1)


def from_mcmc(x) -> Binary:
    t0, t_eff, tE, t_star, piE_N, piE_E, log_s, log_q, alpha = np.moveaxis(x, -1, 0)
    return Binary(t0, t_eff / tE, tE, t_star / tE, piE_N, piE_E, 10**log_s, 10**log_q, alpha)


def log_prob(x):
    """In to_mcmc() coordinates: Gaussian likelihood; priors flat in (t0, u0, tE, rho, piE, alpha) --
    -2 log tE is the Jacobian of (u0, rho) -> (t_eff, t_star), else the prior goes as tE^2 -- and in
    log s, log q (chi2_binary's bounds)."""
    p = from_mcmc(x)
    c2 = chi2_binary(p, EVENT)
    return -0.5 * c2 - 2 * np.log(p.tE) if np.isfinite(c2) else -np.inf


def run_mcmc(ex, best: Binary, path: Path, tag: str, nwalkers=32, nsteps=12000, every=200, seed=42):
    """emcee around a refined minimum, checkpointed to `path` (chain + log_prob, unflattened) every
    `every` steps and resumed from it: a 48 h wall-clock kill used to lose every mode (session 22).
    nsteps: the Cassan chain had tau ~ 230 steps (session 19), so 3000 was only ~13 tau; 12000 aims
    at the ~50 tau emcee recommends. `ex` is shared by every mode at once (see __main__)."""
    chain, lp = np.empty((0, nwalkers, len(best))), np.empty((0, nwalkers))
    if path.exists() and "start" in (d := np.load(path)) and np.allclose(d["start"], best):  # same mode, not index
        chain, lp = d["chain"], d["log_prob"]
        print(f"[mcmc {tag}] resuming at step {len(chain)}")
    rng = np.random.default_rng(seed)
    p0 = chain[-1] if len(chain) else (
        np.array(best) + 0.01 * np.array(binary_steps(best)) * rng.standard_normal((nwalkers, len(best))))
    # DE moves (emcee docs' mix): stretch alone left tau ~500-1200 on *every* parameter of O-03-BLG235 with no
    # burn-in drift or stuck walkers -- a curved posterior it can only cross along walker-to-walker lines (session 23)
    sampler = emcee.EnsembleSampler(nwalkers, len(best), log_prob, pool=ex,
                                    moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])

    def save():  # saved physical, log_prob = -chi2/2 (diagnose reads chi2 off it): the MCMC coordinates stay in here
        c, l = sampler.get_chain(), sampler.get_log_prob()
        assert c is not None and l is not None  # emcee's getters are inferred Optional
        c = np.stack(from_mcmc(c), -1)
        l = l + 2 * np.log(c[..., Binary._fields.index("tE")])
        with open(path.with_suffix(".tmp"), "wb") as f:  # write-then-rename, as run_grid()
            np.savez(f, chain=np.concatenate([chain, c]), log_prob=np.concatenate([lp, l]), labels=Binary._fields,
                     start=np.array(best))
        path.with_suffix(".tmp").replace(path)

    for _ in sampler.sample(to_mcmc(Binary(*p0.T)), iterations=nsteps - len(chain), skip_initial_state_check=len(chain) > 0,
                            progress=True, progress_kwargs={"desc": tag}):
        if sampler.iteration % every == 0:
            save()
    if sampler.iteration:
        save()


def acceptance(chain):
    """Fraction of steps a walker moved (a moved walker accepted) -- works on a saved chain."""
    return np.mean(np.any(np.diff(chain, axis=0) != 0, axis=2))


def save_mcmc(path: Path, out: Path, tag: str):
    """Summary + corner from a finished chain file (main thread: pyplot isn't thread-safe)."""
    d = np.load(path)
    chain, lp = d["chain"], d["log_prob"]
    tau = emcee.autocorr.integrated_time(chain, quiet=True)  # warns (not raises) when nsteps < 50 tau
    print(f"[mcmc {tag}] acceptance={acceptance(chain):.2f}, "
          f"nsteps/max(tau)={len(chain) / np.max(tau):.1f} (want > 50), tau={np.round(tau, 1)}")
    burn, thin = len(chain) // 4, 15
    samples, log_probs = chain[burn::thin].reshape(-1, chain.shape[2]), lp[burn::thin].ravel()
    lo, med, hi = np.percentile(samples, [16, 50, 84], axis=0)
    for label, l, m, h in zip(Binary._fields, lo, med, hi):
        print(f"[mcmc {tag}] {label} = {m:.6g} +{h - m:.3g} -{m - l:.3g}")
    save_corner(samples, list(Binary._fields), samples[np.argmax(log_probs)], str(out / f"mcmc_{tag}_corner.png"))


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


def fmt_params(p, err=None, sig=3):
    """'name = value +hi -lo' per field, shared by fit_lc.png (sig figures) and summary.txt (sig=None:
    full precision). alpha in degrees [0, 360); t0 to 1e-4 d (3 sig figs of 3491.9 say nothing)."""
    out = []
    for name, v in zip(p._fields, p):
        e = err[name] if err else None
        if name == "alpha":
            v, e = np.degrees(v) % 360, e and np.degrees(e)
        f = ((lambda x: repr(float(x))) if sig is None else (lambda x: f"{x:.4f}") if name == "t0"
             else (lambda x: f"{x:#.{sig}g}".rstrip(".")))  # "#": 48.0 not 48; rstrip: 120 not "120."
        out.append(f"{name} = {f(v)}" + (f" +{f(e[0])} -{f(e[1])}" if e is not None else "")
                   + (" deg" if name == "alpha" else ""))
    return out


def plot_caustic_inset(ax, event: Event, best: Binary, t0_par: float, tc: float):
    """Square inset (physical inches + adjustable="datalim": a fraction-based or adjustable="box" inset
    isn't square, CLAUDE.md session 11): the caustics, the fitted trajectory (parallax-curved) with a
    direction arrow, the source disk at tc, both lenses (only if inside the zoom); zoomed on the caustic
    nearest the source at tc, wide enough to hold the source there."""
    t = np.linspace(best.t0 - 3 * best.tE, best.t0 + 3 * best.tE, 20000)
    z = binary_trajectory(t, best.t0, best.u0, best.tE, best.alpha, best.piE_N, best.piE_E,
                          *sun_earth_projection(t, event.coords, t0_par))
    zc = np.interp(tc, t, z.real) + 1j * np.interp(tc, t, z.imag)
    try:
        caustics = caustic_curve(best.s, best.q)
    except ValueError:  # pieces didn't join (topology boundary): trajectory only
        caustics = []
    near = min(caustics, key=lambda c: np.abs(c - zc).min()) if caustics else np.array([zc])
    centre = near.mean()
    half = max(0.75 * max(np.ptp(near.real), np.ptp(near.imag)), 1.2 * abs(zc - centre), 5 * best.rho)
    ins = inset_axes(ax, width=2.2, height=2.2, loc="upper left", borderpad=3.5)  # clear of ax's y labels
    for c in caustics:
        ins.plot(c.real, c.imag, "r-", lw=0.8)
    ins.plot(z.real, z.imag, "k-", lw=0.6)
    ins.add_patch(plt.Circle((zc.real, zc.imag), best.rho, color="C1", alpha=0.7))
    za = np.interp(tc + 0.3 * half * best.tE, t, z.real) + 1j * np.interp(tc + 0.3 * half * best.tE, t, z.imag)
    ins.annotate("", xy=(za.real, za.imag), xytext=(zc.real, zc.imag), arrowprops={"arrowstyle": "->", "lw": 1})
    _, _, z1, z2 = lens_position(best.s, best.q)
    ins.plot([z1, z2], [0, 0], "b+", ms=8)
    ins.set(xlim=(centre.real - half, centre.real + half), ylim=(centre.imag - half, centre.imag + half))
    ins.set_aspect("equal", adjustable="datalim")
    ins.tick_params(labelsize=6)
    ins.set_title(r"source plane ($\theta_E$)", fontsize=7)


def plot_fit(event: Event, fspl: FSPL, best: Binary, t0_par: float, path: Path, text):
    """Every instrument in the reference instrument's magnitude system (the longest-baseline "mag"
    instrument -- OGLE's calibrated I for both current events), the convention of published
    microlensing light curves. Other instruments have their own zero points and bands, so they can
    only be put there through a model: F_ref = fs_ref * A_i + fb_ref with A_i = (F_i - fb_i) / fs_i
    (session 22: magnification itself is model-dependent through blending, FSPL and 2L1S explain the
    same fluxes with A differing ~4x on O-05-BLG169). Top row: data aligned with the 2L1S fit, both
    models as predicted reference-instrument magnitudes (exact for the reference instrument's own
    points). Residual rows: data aligned with that row's own model minus it, in mag, shared axis --
    the fair per-model comparison. Left: peak; right: the anomaly, centred on the point where 2L1S
    lowers chi2 most, half-width 1.5x the median distance of the points gaining > 10% as much (median:
    a few wing points also gain and would stretch a min-max span). The anomaly panel carries the
    caustic inset; `text` (three lists of lines, from diagnose()) fills a strip below."""
    A_f, A_b = fspl_A(event, fspl), binary_A(event, best)
    mags = [n for n, i in enumerate(event.instruments) if i.kind == "mag"]
    k = max(mags or range(len(event.instruments)), key=lambda k: np.ptp(event.instruments[k].time))
    ref = event.instruments[k]
    t = np.concatenate([i.time for i in event.instruments])
    gain = flux_residuals(event, A_f) ** 2 - flux_residuals(event, A_b) ** 2
    tc = t[np.argmax(gain)]
    half = max(1.5 * np.median(np.abs(t[gain > 0.1 * gain.max()] - tc)), 0.05)
    windows = (find_zoom_window(ref.time, ref.flux, ref.flux_err), (tc - half, tc + half))  # model-free window

    def to_mag(F):  # ref-instrument flux -> its magnitude (mag_to_flux's inverse); F <= 0 -> NaN
        return ZERO_POINT_MAG - 2.5 * np.log10(np.where(F > 0, F, np.nan))

    def aligned(i, A_i, A_ref):  # instrument i's data and the model, as ref-instrument flux, under one model
        cr = profile_flux(ref, A_ref)[0]
        fs_r, fb_r = cr[0], (cr[1] if ref.kind == "mag" else 0.0)
        A_data, A_err = to_magnification(i, A_i)
        return fs_r * A_data + fb_r, fs_r * A_err, fs_r * A_i + fb_r, (fs_r, fb_r)

    fig = plt.figure(figsize=(15, 13))
    gs = fig.add_gridspec(4, 2, height_ratios=[4, 1.5, 1.5, 1.7])
    axes = np.empty((3, 2), dtype=object)
    for col in range(2):
        axes[0, col] = fig.add_subplot(gs[0, col])
        axes[1, col] = fig.add_subplot(gs[1, col], sharex=axes[0, col])
        axes[2, col] = fig.add_subplot(gs[2, col], sharex=axes[0, col], sharey=axes[1, col])  # same residual scale
    for col, xlim in enumerate(windows):
        for n, (i, Af, Ab) in enumerate(zip(event.instruments, A_f, A_b)):
            m = (i.time >= xlim[0]) & (i.time <= xlim[1])
            F, F_err, _, _ = aligned(i, Ab, A_b[k])
            axes[0, col].errorbar(i.time[m], to_mag(F)[m], (1.0857 * F_err / F)[m], fmt=".", ms=3, elinewidth=0.5,
                                  color=f"C{n}", label=f"{i.name} ({i.band})" + ("" if n == k else ", aligned by 2L1S"))
            for row, Am, Aref in ((1, Af, A_f[k]), (2, Ab, A_b[k])):
                F, F_err, F_mod, _ = aligned(i, Am, Aref)
                axes[row, col].errorbar(i.time[m], (to_mag(F) - to_mag(F_mod))[m], (1.0857 * F_err / F)[m], fmt=".",
                                        ms=3, elinewidth=0.5, color=f"C{n}")
        tg = np.linspace(*xlim, 2000)
        g = on_grid(event, ref, tg, t0_par)
        for name, Ag, Aref, style in (("2L1S", binary_A(g, best)[0], A_b[k], "k-"),
                                      ("FSPL + parallax", fspl_A(g, fspl)[0], A_f[k], "r--")):
            _, _, _, (fs_r, fb_r) = aligned(ref, Aref, Aref)
            axes[0, col].plot(tg, to_mag(fs_r * Ag + fb_r), style, lw=1.2, label=f"{name} best fit")
        axes[0, col].set(xlim=xlim, ylabel=f"{ref.name} {ref.band} magnitude")
        axes[0, col].invert_yaxis()
        for row, (name, color) in ((1, ("FSPL + parallax", "r")), (2, ("2L1S", "k"))):
            axes[row, col].axhline(0, color=color, lw=0.8)
            axes[row, col].set_ylabel(f"data - {name} (mag)")
        axes[2, col].set_xlabel("HJD - 2450000")
        axes[2, col].ticklabel_format(axis="x", useOffset=False)
    axes[1, 0].invert_yaxis()  # brighter up, as in the light curve (shared by both residual rows)
    axes[1, 1].invert_yaxis()
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].set_title(f"{event.short_name}: peak")
    axes[0, 1].set_title("anomaly (where 2L1S gains most chi2)")
    plot_caustic_inset(axes[0, 1], event, best, t0_par, tc)
    strip = fig.add_subplot(gs[3, :])
    strip.axis("off")
    for x, lines in zip((0.0, 0.36, 0.62), text):
        strip.text(x, 1.0, "\n".join(lines), va="top", family="monospace", fontsize=9, transform=strip.transAxes)
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


def diagnose(base: Event, k_fspl, fspl: FSPL, pgrid, t0_par: float, out: Path, key: str, config: str):
    """(vi) Session 20's done-checks: per MCMC mode nsteps/tau > 50, acceptance 0.2-0.5, its best
    sample within ~1 of the Nelder-Mead chi2 it started from, trace plot; the overall best 2L1S
    (refined[0] or any chain's best sample: the MCMC regularly beats Nelder-Mead, session 21), its +/-
    from the 16-84% of the chain of the mode nearest it (modes aren't pooled: close/wide would merge);
    BIC under K at FSPL (conservative) and K at 2L1S; per instrument fs/fb (fb < 0 flagged: expected
    for FTN's offset scale, suspicious elsewhere). Writes fit_lc.png and summary.txt (session 22)."""
    event = rescale(base, k_fspl)  # the search's footing
    ok = lambda c: "ok" if c else "FAIL"
    refined = cached(out / "refined.npz", key)
    if refined is None or "family" not in refined:
        raise SystemExit(f"{out / 'refined.npz'} missing, from another config/FSPL or pre-session-23: rerun the search")
    trk = track(base, t0_par)
    modes = []  # (tag, best chi2, best sample, chain, done-check line, converged, caustic_crossing(), crossing fraction)
    for f in sorted(out.glob("mcmc_mode*_chain.npz"), key=lambda f: int(f.stem.split("_")[1].removeprefix("mode"))):
        tag = f.stem.removeprefix("mcmc_").removesuffix("_chain")
        d = np.load(f)
        chain, lp = d["chain"], d["log_prob"]  # (nsteps, nwalkers, ndim), (nsteps, nwalkers)
        ntau = len(chain) / emcee.autocorr.integrated_time(chain, quiet=True).max()
        acc = acceptance(chain)
        start = chain[0].mean(axis=0)  # walkers start in a 1%-of-step ball around their refined minimum
        nm = refined["chi2"][np.argmin(np.linalg.norm((refined["theta"] - start) / binary_steps(Binary(*start)), axis=1))]
        i = np.unravel_index(np.argmax(lp), lp.shape)  # the full chain, not the burned-in thinned subset
        best = -2 * lp[i]
        conv = ntau > 50 and 0.2 <= acc <= 0.5 and abs(nm - best) <= 1
        check = (f"nsteps/tau={ntau:.1f} {ok(ntau > 50)}, acceptance={acc:.2f} {ok(0.2 <= acc <= 0.5)}, "
                 f"chi2 Nelder-Mead={nm:.2f} vs MCMC best={best:.2f} {ok(abs(nm - best) <= 1)}")
        print(f"[diag {tag}] {check}")
        # ponytail: 200 random post-burn-in samples (binomial error <~ 3.5%): one caustic_curve() each
        smp = chain[len(chain) // 4:].reshape(-1, chain.shape[2])
        smp = smp[np.random.default_rng(0).choice(len(smp), min(200, len(smp)), replace=False)]
        frac = np.mean([(cc := caustic_crossing(trk, Binary(*x))) is not None and cc[0] for x in smp])
        modes.append((tag, best, Binary(*chain[i]), chain, check, conv, caustic_crossing(trk, Binary(*chain[i])), frac))
        plot_trace(chain, out / f"mcmc_{tag}_trace.png")

    cands = [(refined["chi2"][0], Binary(*refined["theta"][0]))] + [(m[1], m[2]) for m in modes]
    best = min(cands, key=lambda c: c[0])[1]
    c2_f, c2_b = chi2_fspl(fspl, event), chi2_binary(best, event)
    print(f"[diag] overall best 2L1S: chi2={c2_b:.2f} {best}")
    if not np.isfinite(c2_b):
        print("[diag] WARNING: best is outside the current bounds (pre-session-22 outputs?): rerun the search")
    err, err_note = None, "no MCMC chains"
    if modes:
        tag, _, _, chain, _, conv, *_ = min(modes, key=lambda m: np.linalg.norm((np.array(m[2]) - best) / binary_steps(best)))
        lo, med, hi = np.percentile(chain[len(chain) // 4:].reshape(-1, chain.shape[2]), [16, 50, 84], axis=0)
        err = {name: (h - m, m - l) for name, l, m, h in zip(Binary._fields, lo, med, hi)}
        err_note = f"+/- from {tag} (16-84%)" + ("" if conv else ", MCMC NOT CONVERGED")

    k_2l1s = error_scale(base, binary_A(base, best))
    print(f"[2l1s] error rescaling k = {np.round(k_2l1s, 4)}")
    n = sum(i.time.size for i in event.instruments)
    bic_lines = []
    # ponytail: both models are re-scored at their fixed best fits under each K, not refit --
    # a K change only reweights instruments, shifting the best fit slightly
    for name, k in (("K at FSPL, conservative", k_fspl), ("K at 2L1S, optimistic", k_2l1s)):
        ev = rescale(base, k)
        bic_f, bic_b = bic(chi2_fspl(fspl, ev), 6, ev), bic(chi2_binary(best, ev), 9, ev)
        bic_lines.append(f"BIC ({name}): FSPL {bic_f:.2f}, 2L1S {bic_b:.2f}, Delta (FSPL - 2L1S) = {bic_f - bic_b:.2f}")
        print(f"[bic] {bic_lines[-1]}")
    flux_lines = []
    for name, A in (("FSPL", fspl_A(event, fspl)), ("2L1S", binary_A(event, best))):
        for i, A_i in zip(event.instruments, A):
            c = profile_flux(i, A_i)[0]
            fb = f", fb={float(c[1])!r}{' (< 0)' if c[1] < 0 else ''}" if i.kind == "mag" else ""
            flux_lines.append(f"{name} {i.name}: fs={float(c[0])!r}{fb}")
            print(f"[diag] {flux_lines[-1]}")

    # crossing vs near miss (session 23): (label, chi2, Binary, caustic_crossing()) candidates
    ref = [(f"{FAMILIES[f]} start", c2, Binary(*x), caustic_crossing(trk, Binary(*x)), tuple(cell))
           for x, c2, cell, f in zip(refined["theta"], refined["chi2"], refined["cell"], refined["family"])]

    def versus(cands):
        side = lambda c, what: (f"no {what} found" if c is None else
                                f"{what}: chi2 {c[1]:.1f} ({c[0]}, s {c[2].s:.3g}, q {c[2].q:.3g}"
                                + ((", unobserved" if c[3][2] else "") if c[3][0] else
                                   f", closest {c[3][1]:.1f} rho" if np.isfinite(c[3][1]) else
                                   f", closest > {CLOSE / c[2].rho:.0f} rho") + ")")
        x, n = [min((c for c in cands if c[3] is not None and c[3][0] == k), key=lambda c: c[1], default=None)
                for k in (True, False)]
        line = f"{side(x, 'crossing')} | {side(n, 'near miss')}"
        if x is not None and n is not None:
            dc = n[1] - x[1]
            line += f" | Delta = {dc:.1f} -> " + ("close" if abs(dc) <= MODE_DCHI2 else
                                                  "crossing preferred" if dc > 0 else "near miss preferred")
        return line

    overall = versus([(f"refined, {c[0]}", *c[1:4]) for c in ref] + [(m[0], m[1], m[2], m[6]) for m in modes])
    overall += " | MCMC crossing fraction: " + (", ".join(f"{m[0]} {m[7]:.2f}" for m in modes) or "no chains")
    cells = sorted({c[4] for c in ref}, key=lambda cell: min(c[1] for c in ref if c[4] == cell))
    per_min = [f"log s {cell[0]:+.2f}, log q {cell[1]:+.2f}: " + versus([c for c in ref if c[4] == cell]) for cell in cells]
    print(f"[diag] {overall}", *(f"[diag] {r}" for r in per_min), sep="\n")

    ranked = sorted(modes, key=lambda m: m[1])
    kind = lambda m: "?" if m[6] is None else "cross" if m[6][0] else "miss"
    row = lambda m: f"{m[0]:7s} {m[1]:8.1f} {m[2].s:7.3g} {m[2].q:9.3g} {m[2].u0:+10.3g} {kind(m):>5s} {m[7]:5.2f} {ok(m[5])}"
    header = f"{'mode':7s} {'chi2':>8s} {'s':>7s} {'q':>9s} {'u0':>10s} {'class':>5s} {'xfrac':>5s} conv"
    text = (["2L1S overall best", f"({err_note})"] + fmt_params(best, err),
            ["FSPL (parallax)"] + fmt_params(fspl) + ["", f"chi2 FSPL  = {c2_f:.1f}", f"chi2 2L1S  = {c2_b:.1f}",
                                                       f"Delta chi2 = {c2_f - c2_b:.1f}", f"N = {n}"],
            [header] + [row(m) for m in ranked[:6]] + ([f"+{len(ranked) - 6} more, see summary.txt"] if len(ranked) > 6 else []))
    plot_fit(event, fspl, best, t0_par, out / "fit_lc.png", text)
    verdict = plot_parallax_map(pgrid, k_fspl, fspl, c2_f, best, c2_b, out / "fspl_parallax_map.png")

    commit = subprocess.run(["git", "describe", "--always", "--dirty"], capture_output=True, text=True).stdout.strip()
    k_of = lambda ks: ", ".join(f"{i.name}={float(k)!r}" for i, k in zip(event.instruments, ks))
    lines = [f"# {event.short_name} -- search.py summary", f"config: {config}",
             f"written: {datetime.now():%Y-%m-%d %H:%M}, git: {commit or 'unknown'}", "",
             "## FSPL (parallax)", *fmt_params(fspl, sig=None),
             f"parallax grid: step {PIE_STEP}, |piE| <= {PIE_MAX}, both u0 signs; {verdict}"
             + ("; WARNING: best on the |piE| bound" if pgrid["on_bound"] else ""), "",
             f"## 2L1S overall best ({err_note})", *fmt_params(best, err, sig=None), "",
             "## scores (errors rescaled by K at FSPL)", f"chi2 FSPL = {c2_f!r}", f"chi2 2L1S = {c2_b!r}",
             f"N = {n}, k FSPL = {6 + sum(2 if i.kind == 'mag' else 1 for i in event.instruments)}, "
             f"k 2L1S = {9 + sum(2 if i.kind == 'mag' else 1 for i in event.instruments)}", *bic_lines, "",
             "## error rescaling", f"K at FSPL: {k_of(k_fspl)}", f"K at 2L1S: {k_of(k_2l1s)}", "",
             "## flux (fs, fb) at the errors above", *flux_lines, "",
             "## crossing vs near miss (source centre crosses a caustic in the data span; Delta = near miss - crossing)",
             overall, "per refined local minimum (grid cell):", *per_min, "",
             "## modes (best sample per MCMC chain, by chi2)", header, *[row(m) for m in ranked], "",
             "## done-checks", *[f"{m[0]}: {m[4]}" for m in modes]]
    (out / "summary.txt").write_text("\n".join(lines) + "\n")
    print(f"saved {out / 'summary.txt'}")


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


def plot_parallax_map(pgrid, k_fspl, fspl: FSPL, c2_f: float, best: Binary, c2_b: float, path: Path):
    """Parallax-only Delta chi2(piE_N, piE_E), one panel per u0 sign, at the K-at-FSPL errors (cells'
    per-instrument chi2 / k^2: fit parameters stay at their base-error optimum, a small approximation).
    Contours at the 1/2/3-sigma levels for 2 parameters; local minima within Delta chi2 < 25 marked;
    the verdict counts separate minima within Delta chi2 < 11.8: more than one -> MCMC recommended.
    Returns the verdict line (summary.txt reuses it)."""
    pie = pgrid["pie"]
    total = (pgrid["c2i"] / np.asarray(k_fspl) ** 2).sum(axis=-1)  # (2, n, n)
    d = total - min(np.min(total), c2_f)
    mins = [(g, i, j) for g in range(2) for i, j in local_minima(total[g], n=pie.size ** 2) if d[g, i, j] < 25]
    n_sep = sum(d[g, i, j] < 11.8 for g, i, j in mins)
    verdict = ("single minimum within Delta chi2 < 11.8: MCMC optional" if n_sep <= 1 else
               f"{n_sep} separate minima within Delta chi2 < 11.8: MCMC recommended")
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.8), sharey=True)
    for g, (ax, label) in enumerate(zip(axes, ("u0 > 0", "u0 < 0"))):
        dg = np.ma.masked_invalid(d[g])
        mesh = ax.pcolormesh(pie, pie, np.log10(1 + dg).T, shading="nearest", cmap="viridis_r",
                             vmin=0, vmax=np.log10(1 + np.nanmax(np.where(np.isfinite(d), d, np.nan))))
        cs = ax.contour(pie, pie, dg.filled(np.nan).T, levels=[2.3, 6.2, 11.8], colors=["w", "w", "w"],
                        linewidths=[1.5, 1, 0.6])
        ax.clabel(cs, fmt={2.3: "1σ", 6.2: "2σ", 11.8: "3σ"}, fontsize=8)
        ax.add_patch(plt.Circle((0, 0), PIE_MAX, fill=False, ls=":", color="grey"))
        for gg, i, j in mins:
            if gg == g:
                ax.plot(pie[i], pie[j], "rx", ms=8, mew=1.5)
        if (fspl.u0 > 0) == (g == 0):
            ax.plot(fspl.piE_N, fspl.piE_E, "r*", ms=16, mec="k", label="parallax-only best")
        ax.plot(best.piE_N, best.piE_E, "D", color="orange", mec="k", ms=8,
                label=f"2L1S best (u0 {'>' if best.u0 > 0 else '<'} 0)")
        ax.set(title=f"FSPL + parallax, {label}", xlabel=r"$\pi_{E,N}$", aspect="equal")
        ax.legend(loc="lower left", fontsize=8)
    axes[0].set_ylabel(r"$\pi_{E,E}$")
    fig.colorbar(mesh, ax=axes, label=r"$\log_{10}(1 + \Delta\chi^2)$ vs parallax-only best", shrink=0.8)
    fig.text(0.02, 0.02, f"chi2 parallax-only = {c2_f:.1f}    chi2 2L1S = {c2_b:.1f}    Delta = {c2_f - c2_b:.1f}"
             f"    (errors rescaled by K at FSPL)\nlocal minima: x (Delta chi2 < 25)    {verdict}"
             + ("    WARNING: best on the |piE| bound" if pgrid["on_bound"] else ""), family="monospace", fontsize=9)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {path}")
    return verdict


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
    if args.stage == "raw":  # quick look: fit_fspl's two-start parallax fit is enough here
        raise SystemExit
    with pool(base, plain) as ex:
        fspl, pgrid = parallax_search(base, plain, fspl, ex)
    k_fspl = error_scale(base, fspl_A(base, fspl))
    assert np.all(np.isfinite(k_fspl)), f"error rescaling k = {k_fspl}"
    event = rescale(base, k_fspl)  # the search runs on FSPL-rescaled errors
    print(f"[fspl] error rescaling k = {np.round(k_fspl, 4)}")
    # one key for every cache below: a config edit (instruments, K, n_alpha, grid...) invalidates them all;
    # run_grid() also checks `plain` (code changes move it without touching the TOML: RHO_MIN, session 22).
    # ponytail: hashes the TOML, not the data files it points to -- edit data in place => delete the caches
    key = hashlib.sha256(Path(args.config).read_bytes()).hexdigest()
    if args.stage == "diagnose":
        diagnose(base, k_fspl, fspl, pgrid, plain.t0, out, key, args.config)
        raise SystemExit

    grid = run_grid(event, plain, out / "grid.npz", key, k_fspl)
    minima = local_minima(grid["chi2"])
    assert len(minima), "no finite cell in the (s, q) grid"
    plot_map(grid, minima, out / "delta_chi2_map.png")

    with pool(event) as ex:
        r = cached(out / "refined.npz", key)  # delete refined.npz to re-refine
        if r is not None and ("plain" not in r or not np.allclose(r["plain"], plain)  # grid was rebuilt under it
                              or "family" not in r):  # pre-session-23: one family per minimum
            r = None
        if r is not None:
            refined = list(zip(r["theta"], r["chi2"], r["cell"], r["family"]))
            print(f"[refine] loaded cached {out / 'refined.npz'}")
        else:  # both families' best start at every minimum, each with its mirror (session 23)
            starts = [((grid["log_s"][i], grid["log_q"][j]), f, x0) for i, j in minima for f in range(len(FAMILIES))
                      if np.isfinite(grid["fam_chi2"][i, j, f]) for x0 in refine_starts(grid["fam_theta"][i, j, f])]
            refined = sorted(((x, c2, cell, f) for (cell, f, _), (x, c2) in zip(starts, ex.map(refine, [s[2] for s in starts]))),
                             key=lambda r: r[1])
            for f in out.glob("mcmc_mode*"):  # new modes: old chains (resume points, plots) are stale
                f.unlink()
            np.savez(out / "refined.npz", theta=[r[0] for r in refined], chi2=[r[1] for r in refined],
                     cell=[r[2] for r in refined], family=[r[3] for r in refined], key=key, plain=plain)
        trk = track(base, plain.t0)
        crosses = [(cc := caustic_crossing(trk, Binary(*x))) is not None and cc[0] for x, *_ in refined]
        for (x, c2, _, f), cross in zip(refined, crosses):
            print(f"[refine] chi2={c2:.2f} {FAMILIES[f]} start, {'crossing' if cross else 'near miss'} {Binary(*x)}")
        modes = distinct_modes(refined, crosses)
        assert modes, "every refined chi2 is inf"
        for f in out.glob("mcmc_mode*"):  # orphans from a run with more modes (globbed by diagnose)
            if int(f.name.split("_")[1].removeprefix("mode")) >= len(modes):
                f.unlink()
        # every mode at once, one thread each, on the one pool: emcee maps only nwalkers/2 walkers at
        # a time and each step waits for its slowest chi2, so a lone mode leaves most cores idle
        with ThreadPoolExecutor(len(modes)) as threads:
            list(threads.map(lambda kp: run_mcmc(ex, kp[1], out / f"mcmc_mode{kp[0]}_chain.npz", f"mode{kp[0]}"),
                             enumerate(modes)))
    for k in range(len(modes)):
        save_mcmc(out / f"mcmc_mode{k}_chain.npz", out, f"mode{k}")

    diagnose(base, k_fspl, fspl, pgrid, plain.t0, out, key, args.config)
