"""Multi-start Nelder-Mead 2L1S fit for OGLE-2003-BLG-235/MOA-2003-BLG-53.

Not yet graduated into the regular pipeline (see CLAUDE.md) -- lives in
scratch/ with the project's other one-time/dev scripts, run manually as
`python3 scratch/fit_2l1s.py` from the project root.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from itertools import product
from multiprocessing import get_context

from lc_models import binary_magnification, binary_trajectory, caustic_curve, mag_to_flux
from preprocess_binary_data import load_raw, moa_to_magnification, ogle_to_magnification
from zoom_utils import find_zoom_window

SHORT_NAME = "O-03-BLG235"

# Raw flux, not the old data/processed/*_magnification.dat -- those froze fs/fb from a
# one-time PSPL point estimate that's since been shown biased for a real caustic-crossing
# trajectory (see CHANGELOG session 5). chi2() below re-solves fs/fb per trial instead.
ogle_time, ogle_mag, ogle_mag_err, moa_time, moa_flux, moa_flux_err = load_raw()
ogle_flux, ogle_flux_err = mag_to_flux(ogle_mag, ogle_mag_err)

# t0/u0/tE seeded from the existing PSPL joint fit's posterior median.
param, val = np.loadtxt(f"data/processed/{SHORT_NAME}_fit_summary.dat", unpack=True, usecols=(0, 2), dtype=str)
guess = dict(zip(param, val.astype(float)))

# Multi-start only over (s, q, alpha): close/resonant/wide topologies x approach angle.
s_list = [0.7, 1, 1.3, 1.6]
q_list = [0.01, 0.1]
alpha_list = [i * np.pi / 8 for i in range(16)]
seeds = [[guess["t0"], guess["u0"], guess["tE"], alpha, s, q] for s, q, alpha in product(s_list, q_list, alpha_list)]


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


def chi2(theta, use_ogle=True):
    """Chi2 of the 2L1S model against raw flux, with each instrument's flux calibration
    (fs/fb) profiled analytically for this specific trial trajectory. `use_ogle=False`
    fits MOA alone (only fs_moa profiled) -- a calibration-ambiguity-free control on the
    joint fit, since there's no cross-instrument scale to get wrong (see CHANGELOG)."""
    t0, u0, tE, alpha, s, q = theta
    if tE <= 0 or s <= 0 or q <= 0:
        return np.inf  # unphysical
    A_moa = binary_magnification(binary_trajectory(moa_time, t0, u0, tE, alpha), s, q)
    if not use_ogle:
        fs_moa = _profile_fs_moa(A_moa)
        if fs_moa <= 0:
            return np.inf
        return np.sum(((moa_flux - fs_moa * (A_moa - 1.0)) / moa_flux_err) ** 2)

    A_ogle = binary_magnification(binary_trajectory(ogle_time, t0, u0, tE, alpha), s, q)
    try:
        fs_ogle, fb_ogle, fs_moa = profile_flux(A_ogle, A_moa)
    except np.linalg.LinAlgError:
        return np.inf  # degenerate trial (e.g. flat A_ogle) -- singular normal equations
    if fs_ogle <= 0 or fs_moa <= 0:
        return np.inf
    resid_ogle = (ogle_flux - (fs_ogle * A_ogle + fb_ogle)) / ogle_flux_err
    resid_moa = (moa_flux - fs_moa * (A_moa - 1.0)) / moa_flux_err
    return np.sum(resid_ogle ** 2) + np.sum(resid_moa ** 2)


def _fit_one_seed(seed, use_ogle=True):
    """Nelder-Mead on one seed at loose tolerance; module-level so it can be pickled."""
    t0, u0, tE, alpha, s, q = seed
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
    print(f"[{tag}] final chi2={result.fun:.2f}")
    for label, value in zip(["t0", "u0", "tE", "alpha", "s", "q"], result.x):
        print(f"[{tag}] {label} = {value:.5f}")

    plot_fit(result.x, use_ogle=use_ogle)
    return result


def run_fit_moa_only():
    return run_fit(use_ogle=False)


def plot_fit(theta, use_ogle=True):
    """Overlay the 2L1S model on the data, full baseline + zoomed on the peak.

    Magnification-space view only, for plotting -- inverts the raw flux using this
    theta's own profiled fs/fb (chi2's calibration), not the frozen PSPL one.
    """
    t0, u0, tE, alpha, s, q = theta
    A_moa_best = binary_magnification(binary_trajectory(moa_time, t0, u0, tE, alpha), s, q)
    if use_ogle:
        A_ogle_best = binary_magnification(binary_trajectory(ogle_time, t0, u0, tE, alpha), s, q)
        fs_ogle, fb_ogle, fs_moa = profile_flux(A_ogle_best, A_moa_best)
        ogle_A, ogle_A_err = ogle_to_magnification(ogle_mag, ogle_mag_err, fs_ogle, fb_ogle)
    else:
        fs_moa = _profile_fs_moa(A_moa_best)
    moa_A, moa_A_err = moa_to_magnification(moa_flux, moa_flux_err, fs_moa)
    time = np.concatenate([ogle_time, moa_time]) if use_ogle else moa_time
    A_obs = np.concatenate([ogle_A, moa_A]) if use_ogle else moa_A

    zoom_start, zoom_end = find_zoom_window(moa_time, moa_A, moa_A_err, padding_fraction=0.3)

    def plot_panel(ax, xlim=None):
        t_grid = np.linspace(*(xlim if xlim else (time.min(), time.max())), 3000)
        A_model = binary_magnification(binary_trajectory(t_grid, t0, u0, tE, alpha), s, q)

        if use_ogle:
            ax.errorbar(ogle_time, ogle_A, yerr=ogle_A_err, fmt="+", ms=3, elinewidth=0.5,
                        color="black", label="OGLE")
        ax.errorbar(moa_time, moa_A, yerr=moa_A_err, fmt="+", ms=3, elinewidth=0.5,
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

    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(2, 2)
    ax_full = fig.add_subplot(gs[0, :])
    ax_zoom = fig.add_subplot(gs[1, 0])
    ax_caustic = fig.add_subplot(gs[1, 1])

    plot_panel(ax_full)
    ax_full.set_title("2L1S fit, full baseline")
    ax_full.legend(loc="upper right")
    plot_panel(ax_zoom, xlim=(zoom_start, zoom_end))
    ax_zoom.set_title("zoomed on peak (auto-detected)")
    ax_zoom.set_xlabel("HJD - 2450000")

    caustic = caustic_curve(s, q)
    t_traj = np.linspace(zoom_start, zoom_end, 3000)
    traj = binary_trajectory(t_traj, t0, u0, tE, alpha)
    ax_caustic.scatter(caustic.real, caustic.imag, s=0.5, color="crimson", label="caustic")
    ax_caustic.plot(traj.real, traj.imag, color="black", lw=1, label="source trajectory")
    ax_caustic.set_aspect("equal")
    ax_caustic.set_xlabel("Re(zeta)")
    ax_caustic.set_ylabel("Im(zeta)")
    ax_caustic.set_title(f"caustic geometry (s={s:.3f}, q={q:.4f})")
    ax_caustic.legend(loc="upper right", fontsize=8)

    fig.tight_layout()

    Path("scratch").mkdir(exist_ok=True)
    suffix = "" if use_ogle else "_moa_only"
    out_path = f"scratch/{SHORT_NAME}_2l1s{suffix}.png"
    fig.savefig(out_path, dpi=300)
    print(f"saved {out_path}")


if __name__ == "__main__":
    run_fit()
    run_fit_moa_only()


