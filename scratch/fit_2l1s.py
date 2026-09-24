"""Multi-start Nelder-Mead 2L1S fit for OGLE-2003-BLG-235/MOA-2003-BLG-53.

Not yet graduated into the regular pipeline (see CLAUDE.md) -- lives in
scratch/ with the project's other one-time/dev scripts, run manually as
`python3 scratch/fit_2l1s.py` from the project root.
"""
import sys
from pathlib import Path
from typing import NamedTuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes as make_square_inset_axes
from scipy.optimize import minimize
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from itertools import product
from multiprocessing import get_context

from lc_models import binary_magnification, binary_trajectory, caustic_curve, mag_to_flux, sun_earth_projection
from preprocess_binary_data import load_raw, moa_to_magnification, ogle_to_magnification, fit_joint_pspl, COORDS
from zoom_utils import find_zoom_window

class TwoL1SParams(NamedTuple):
    """The 2L1S track's 8 physical parameters, in one canonical order -- every
    call site (seeds, residuals, plot_fit, MCMC theta) unpacks through this
    instead of re-deriving the order by hand. A previously real bug (see
    CHANGELOG session 10): residuals() was once unpacked in a different order
    than seeds/plot_fit/LABELS, silently swapping s/q with piE_N/piE_E."""
    t0: float
    u0: float
    tE: float
    alpha: float
    piE_N: float
    piE_E: float
    s: float
    q: float


SHORT_NAME = "O-03-BLG235"

# Raw flux, not the old data/processed/*_magnification.dat -- those froze fs/fb from a
# one-time PSPL point estimate that's since been shown biased for a real caustic-crossing
# trajectory (see CHANGELOG session 5). chi2() below re-solves fs/fb per trial instead.
ogle_time, ogle_mag, ogle_mag_err, moa_time, moa_flux, moa_flux_err = load_raw()
ogle_flux, ogle_flux_err = mag_to_flux(ogle_mag, ogle_mag_err)

zeros_ogle, zeros_moa = np.zeros_like(ogle_time), np.zeros_like(moa_time)
plain_fit, _ = fit_joint_pspl(ogle_time, ogle_mag, ogle_mag_err, moa_time, moa_flux, moa_flux_err,
                                zeros_ogle, zeros_ogle, zeros_moa, zeros_moa)
t0_par = plain_fit[0]
print(f"t0_par (parallax reference epoch, from plain pre-parallax fit) = {t0_par:.5f}")

delta_sN_ogle, delta_sE_ogle = sun_earth_projection(ogle_time, COORDS, t0_par)
delta_sN_moa, delta_sE_moa = sun_earth_projection(moa_time, COORDS, t0_par)

# t0/u0/tE seeded from the existing PSPL joint fit's posterior median.
param, val = np.loadtxt(f"data/processed/{SHORT_NAME}_fit_summary.dat", unpack=True, usecols=(0, 2), dtype=str)
guess = dict(zip(param, val.astype(float)))

# Multi-start only over (s, q, alpha): close/resonant/wide topologies x approach angle.
# Same overall range as before (s: 0.4-1.9, q: 0.001-1.0 log-spaced, alpha: full circle),
# just fewer/wider-spaced points -- 4x3x8=96 seeds instead of 384, on a 6-core machine
# (ProcessPoolExecutor defaults to os.cpu_count() workers) that's 16/core instead of 64/core.
s_list = [0.5, 1.0, 1.5, 2.0]
q_list = [0.001, 0.03, 1.0]
alpha_list = [i * np.pi / 4 for i in range(8)]
seeds = [TwoL1SParams(t0=guess["t0"], u0=guess["u0"], tE=guess["tE"], alpha=alpha,
                       piE_N=guess["piE_N"], piE_E=guess["piE_E"], s=s, q=q)
         for s, q, alpha in product(s_list, q_list, alpha_list)]


def _profile_fs_moa(A_moa):
    """Closed-form weighted-least-squares fs_moa alone (MOA has no additive blend term)."""
    x_moa = A_moa - 1.0
    w_moa = 1.0 / moa_flux_err ** 2
    return np.sum(w_moa * x_moa * moa_flux) / np.sum(w_moa * x_moa ** 2)


def profile_flux(A_ogle, A_moa):
    """Closed-form weighted-least-squares fs_ogle/fb_ogle/fs_moa for one trial's model
    magnification -- replaces the old frozen-PSPL-calibration approach (see CHANGELOG).
    ponytail: OGLE side is solved in (linearized) flux space, not the exact magnitude-space
    likelihood fit_joint_pspl uses -- fine while mag errors are small; revisit if this ever
    needs to match a mag-space chi2 exactly.
    """
    w = 1.0 / ogle_flux_err ** 2
    X = np.column_stack([A_ogle, np.ones_like(A_ogle)])
    fs_ogle, fb_ogle = np.linalg.solve((X * w[:, None]).T @ X, (X * w[:, None]).T @ ogle_flux)
    return fs_ogle, fb_ogle, _profile_fs_moa(A_moa)

def residuals(theta, use_ogle=True):
    """Per-point standardized residuals of the 2L1S model against raw flux.
    Always returns an array shaped like the data -- an unphysical/degenerate
    trial fills it with np.inf rather than returning a bare scalar, so callers
    check np.isfinite() alone instead of also guarding for np.isscalar()."""
    t0, u0, tE, alpha, piE_N, piE_E, s, q = TwoL1SParams(*theta)
    invalid = np.full(len(moa_time) + (len(ogle_time) if use_ogle else 0), np.inf)

    if tE <= 0 or s <= 0 or q <= 0:
        return invalid  # unphysical
    if np.abs(piE_N) > 2 or np.abs(piE_E) > 2:
        return invalid
    A_moa = binary_magnification(binary_trajectory(moa_time, t0, u0, tE, alpha, piE_N, piE_E, delta_sN_moa, delta_sE_moa), s, q)
    if not use_ogle:
        fs_moa = _profile_fs_moa(A_moa)
        if fs_moa <= 0:
            return invalid
        return (moa_flux - fs_moa * (A_moa - 1.0)) / moa_flux_err

    A_ogle = binary_magnification(binary_trajectory(ogle_time, t0, u0, tE, alpha, piE_N, piE_E, delta_sN_ogle, delta_sE_ogle), s, q)
    try:
        fs_ogle, fb_ogle, fs_moa = profile_flux(A_ogle, A_moa)
    except np.linalg.LinAlgError:
        return invalid  # degenerate trial (e.g. flat A_ogle) -- singular normal equations
    if fs_ogle <= 0 or fs_moa <= 0:
        return invalid
    resid_ogle = (ogle_flux - (fs_ogle * A_ogle + fb_ogle)) / ogle_flux_err
    resid_moa = (moa_flux - fs_moa * (A_moa - 1.0)) / moa_flux_err
    return np.concatenate([resid_ogle, resid_moa])

def chi2(theta, use_ogle=True):
    """Chi2 of the 2L1S model against raw flux, with each instrument's flux calibration
    (fs/fb) profiled analytically for this specific trial trajectory. `use_ogle=False`
    fits MOA alone (only fs_moa profiled) -- a calibration-ambiguity-free control on the
    joint fit, since there's no cross-instrument scale to get wrong (see CHANGELOG)."""
    return np.sum(residuals(theta, use_ogle)**2)

def huber(residual, delta):
    loss = np.where(np.abs(residual) <= delta, 0.5 * residual ** 2, delta * (np.abs(residual) - 0.5 * delta))
    return np.sum(loss)

def _fit_one_seed(seed, use_ogle=True):
    """Nelder-Mead on one seed at loose tolerance; module-level so it can be pickled."""
    t0, u0, tE, alpha, piE_N, piE_E, s, q = TwoL1SParams(*seed)
    tag = "2l1s" if use_ogle else "2l1s-moa-only"
    result = minimize(chi2, x0=seed, args=(use_ogle,), method="Nelder-Mead",
                       options={"xatol": 1e-2, "fatol": 1e-2, "maxiter": 2000})
    print(f"[{tag}] seed s={s:.2f} q={q:.3f} alpha={alpha:.2f} -> chi2={result.fun:.2f}")
    return result


def run_fit(use_ogle=True):
    """Parallel multi-start search, then a tight-tolerance refit of the best seed."""
    tag = "2l1s" if use_ogle else "2l1s-moa-only"
    print(f"[{tag}] starting multi-start search: {len(seeds)} seeds")
    # spawn (not fork) -- forking after the parent process has touched torch (which spins up
    # its own internal thread pool) can deadlock the child on a lock none of its threads hold;
    # hit exactly this running run_fit() then run_fit_moa_only() back to back (see CHANGELOG).
    with ProcessPoolExecutor(mp_context=get_context("spawn")) as exec:
        results = list(exec.map(partial(_fit_one_seed, use_ogle=use_ogle), seeds))
    best = min(results, key=lambda r: r.fun)
    print(f"[{tag}] multi-start done, best chi2={best.fun:.2f}, refining...")

    result = minimize(chi2, x0=best.x, args=(use_ogle,), method="Nelder-Mead",
                       options={"xatol": 1e-6, "fatol": 1e-6, "maxiter": 20000})
    best_params = TwoL1SParams(*result.x)
    print(f"[{tag}] final chi2={result.fun:.2f}")
    for label, value in zip(TwoL1SParams._fields, best_params):
        print(f"[{tag}] {label} = {value:.5f}")

    plot_fit(best_params, use_ogle=use_ogle, tag="_nelder_mead")
    return result


def run_fit_moa_only():
    return run_fit(use_ogle=False)


def plot_fit(theta, use_ogle=True, *, tag):
    """Overlay the 2L1S model on the data: event season (HJD 2700-3000), then a
    zoomed peak panel carrying a caustic-geometry inset and a raw-residual panel
    beneath it.

    Magnification-space view only, for plotting -- inverts the raw flux using this
    theta's own profiled fs/fb (chi2's calibration), not the frozen PSPL one.

    tag: required, method-identifying suffix appended to the output filename
    (e.g. "_nelder_mead", "_mcmc_studentt", "_mcmc_chi2") -- no default, so a new
    call site can't silently collide with an existing method's output the way
    run_fit()/run_mcmc() used to (see CHANGELOG).
    """
    t0, u0, tE, alpha, piE_N, piE_E, s, q = TwoL1SParams(*theta)

    def model_fn(t):
        delta_sN, delta_sE = sun_earth_projection(t, COORDS, t0_par)
        return binary_magnification(binary_trajectory(t, t0, u0, tE, alpha, piE_N, piE_E, delta_sN, delta_sE), s, q)

    A_moa_best = model_fn(moa_time)
    if use_ogle:
        A_ogle_best = model_fn(ogle_time)
        fs_ogle, fb_ogle, fs_moa = profile_flux(A_ogle_best, A_moa_best)
        ogle_A, ogle_A_err = ogle_to_magnification(ogle_mag, ogle_mag_err, fs_ogle, fb_ogle)
        ogle_resid = ogle_A - A_ogle_best
    else:
        fs_moa = _profile_fs_moa(A_moa_best)
    moa_A, moa_A_err = moa_to_magnification(moa_flux, moa_flux_err, fs_moa)
    moa_resid = moa_A - A_moa_best
    time = np.concatenate([ogle_time, moa_time]) if use_ogle else moa_time
    A_obs = np.concatenate([ogle_A, moa_A]) if use_ogle else moa_A
    zoom_start, zoom_end = find_zoom_window(moa_time, moa_A, moa_A_err, padding_fraction=0.3)

    def plot_panel(ax, xlim=None):
        t_grid = np.linspace(*(xlim if xlim else (time.min(), time.max())), 3000)
        A_model = model_fn(t_grid)

        if use_ogle:
            ax.errorbar(ogle_time, ogle_A, yerr=ogle_A_err, fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5,
                        color="black", label="OGLE")
        ax.errorbar(moa_time, moa_A, yerr=moa_A_err, fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5,
                    color="tab:orange", label="MOA")
        ax.plot(t_grid, A_model, color="crimson", lw=1.5, label="2L1S fit")
        ax.set_ylabel("Magnification A(t)")

        if xlim is not None:
            ax.set_xlim(*xlim)
            in_zoom = (time >= xlim[0]) & (time <= xlim[1])
            # include the model curve, not just the data -- a sharp model peak/dip
            # between data points would otherwise get clipped by the y-limits.
            y_vals = np.concatenate([A_obs[in_zoom], A_model])
            if y_vals.size:
                pad = 0.1 * (y_vals.max() - y_vals.min())
                ax.set_ylim(y_vals.min() - pad, y_vals.max() + pad)

    fig = plt.figure(figsize=(9, 12))
    gs = fig.add_gridspec(3, 1, height_ratios=[4, 4, 1.5])
    ax_season = fig.add_subplot(gs[0])
    ax_zoom = fig.add_subplot(gs[1])
    ax_resid = fig.add_subplot(gs[2], sharex=ax_zoom)

    # event season, not the multi-year full baseline -- same HJD 2700-3000 window
    # scratch/compare_2l1s_fits.py already uses for O-03-BLG235.
    plot_panel(ax_season, xlim=(2700, 3000))
    ax_season.set_title("2L1S fit, event season (HJD 2700-3000)")
    ax_season.legend(loc="upper right")

    plot_panel(ax_zoom, xlim=(zoom_start, zoom_end))
    ax_zoom.set_title("zoomed on peak (auto-detected)")
    ax_zoom.tick_params(labelbottom=False)

    # caustic geometry as an inset in the zoomed panel's corner, rather than its own subplot
    caustic = caustic_curve(s, q)
    t_traj = np.linspace(zoom_start, zoom_end, 3000)
    delta_sN_traj, delta_sE_traj = sun_earth_projection(t_traj, COORDS, t0_par)
    traj = binary_trajectory(t_traj, t0, u0, tE, alpha, piE_N, piE_E, delta_sN_traj, delta_sE_traj)
    ax_caustic = make_square_inset_axes(ax_zoom, width=1.8, height=1.8, loc="upper right", borderpad=2.0)
    ax_caustic.scatter(caustic.real, caustic.imag, s=0.5, color="crimson")
    ax_caustic.plot(traj.real, traj.imag, color="black", lw=1)
    # adjustable="datalim" (not the default "box") -- keeps the inset's box the square
    # shape we just requested from inset_axes(), and instead pads the *data* limits on
    # whichever axis is narrower so x/y data units stay equal-scaled.
    ax_caustic.set_aspect("equal", adjustable="datalim")
    ax_caustic.set_title(f"s={s:.3f}, q={q:.4f}", fontsize=8)
    ax_caustic.tick_params(labelsize=6)

    # raw (non-standardized) residuals, real per-instrument error bars -- matches
    # zoom_utils.plot_fit_panels()'s convention (see CLAUDE.md) rather than the
    # standardized plot_residual_panel()/plot_residual_hist() helpers, which hide
    # OGLE's real ~4x-better precision when both instruments share one panel.
    in_zoom_moa = (moa_time >= zoom_start) & (moa_time <= zoom_end)
    ax_resid.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    resid_parts = [moa_resid[in_zoom_moa]]
    if use_ogle:
        in_zoom_ogle = (ogle_time >= zoom_start) & (ogle_time <= zoom_end)
        ax_resid.errorbar(ogle_time[in_zoom_ogle], ogle_resid[in_zoom_ogle], yerr=ogle_A_err[in_zoom_ogle],
                           fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5, color="black")
        resid_parts.append(ogle_resid[in_zoom_ogle])
    ax_resid.errorbar(moa_time[in_zoom_moa], moa_resid[in_zoom_moa], yerr=moa_A_err[in_zoom_moa],
                       fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5, color="tab:orange")
    all_resid_z = np.concatenate(resid_parts)
    if all_resid_z.size:
        ylim = np.abs(all_resid_z).max() * 1.1
        ax_resid.set_ylim(-ylim, ylim)
    ax_resid.set_ylabel("residual (A(t))")
    ax_resid.set_xlabel("HJD - 2450000")

    fig.tight_layout()

    Path("scratch").mkdir(exist_ok=True)
    suffix = ("" if use_ogle else "_moa_only") + tag
    out_path = f"scratch/{SHORT_NAME}_2l1s{suffix}.png"
    fig.savefig(out_path, dpi=600)
    print(f"saved {out_path}")


if __name__ == "__main__":
    run_fit()
    run_fit_moa_only()


