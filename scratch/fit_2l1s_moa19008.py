"""Multi-start Nelder-Mead 2L1S fit for MOA-2019-BLG-008 (Bachelet et al. 2022),
a real object at the planet/brown-dwarf boundary (~30 M_Jup) -- directly on
this project's central brown-dwarf-mass-lens theme.

Uses only its KMT I-band data for now (source contains "KMT", band == "I"),
filtered straight out of the raw multi-survey file rather than a notebook-derived
copy. Single instrument -- unlike O-03-BLG235's OGLE+MOA case, there's no
cross-survey unit mismatch to jointly calibrate, so fs/fb (the flux scale and
blend) are solved in closed form (linear least squares) at every trial theta
instead of frozen from a separate preprocessing step. This also sidesteps the
frozen-calibration bias found in fit_2l1s.py's O-03-BLG235 pipeline (see
CHANGELOG.md) -- there's nothing to freeze here.

t0/u0/tE guesses below come from this dataset's own light curve (peak time,
rough half-depth width), NOT reused from O-03-BLG235 -- see CHANGELOG.md's
u0-sign bug for why that's not safe to share across datasets.

Not yet graduated into the regular pipeline (see CLAUDE.md) -- lives in
scratch/ with the project's other one-time/dev scripts, run manually as
`python3 scratch/fit_2l1s_moa19008.py` from the project root.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from concurrent.futures import ProcessPoolExecutor
from itertools import product

from lc_models import ZERO_POINT_MAG, binary_magnification, binary_trajectory, caustic_curve
from zoom_utils import find_zoom_window

SHORT_NAME = "M-19-BLG008"  # see dataset_names.txt

# raw file has no header: source  band  HJD(full)  mag  mag_err, all surveys mixed together
source, band, hjd, mag, mag_err = np.loadtxt("data/MOA-2019-BLG-008L.dat", dtype=str, unpack=True)
kmt_i = np.char.startswith(source, "KMT") & (band == "I")

time = hjd[kmt_i].astype(float) - 2450000.0
mag = mag[kmt_i].astype(float)
mag_err = mag_err[kmt_i].astype(float)
order = np.argsort(time)
time, mag, mag_err = time[order], mag[order], mag_err[order]

flux_obs = 10 ** (-0.4 * (mag - ZERO_POINT_MAG))
flux_err = flux_obs * mag_err * (np.log(10) / 2.5)

# seeded from this dataset's own peak (argmin mag) and half-depth width, not another dataset's fit
t0_guess = float(time[np.argmin(mag)])
u0_guess = 0.2
tE_guess = 50.0

s_list = [1]
q_list = [0.01]
alpha_list = [i * np.pi / 8 for i in range(16)]
seeds = [[t0_guess, u0_guess, tE_guess, alpha, s, q] for s, q, alpha in product(s_list, q_list, alpha_list)]


def solve_flux_calibration(A_model):
    """Closed-form weighted least squares for (fs, fb) at fixed magnification shape --
    single instrument, so no joint cross-survey calibration is needed, just this."""
    w = 1 / flux_err**2
    X = np.column_stack([A_model, np.ones_like(A_model)])
    return np.linalg.solve(X.T @ (X * w[:, None]), X.T @ (w * flux_obs))


def chi2(theta):
    """Chi2 of the 2L1S model against KMT I-band flux, profiling fs/fb at every trial."""
    t0, u0, tE, alpha, s, q = theta
    if tE <= 0 or s <= 0 or q <= 0:
        return np.inf
    A_model = binary_magnification(binary_trajectory(time, t0, u0, tE, alpha), s, q)
    fs, fb = solve_flux_calibration(A_model)
    return np.sum(((flux_obs - (fs * A_model + fb)) / flux_err) ** 2)


def _fit_one_seed(seed):
    """Nelder-Mead on one seed at loose tolerance; module-level so it can be pickled."""
    t0, u0, tE, alpha, s, q = seed
    result = minimize(chi2, x0=seed, method="Nelder-Mead",
                       options={"xatol": 1e-2, "fatol": 1e-2, "maxiter": 2000})
    print(f"[2l1s-moa19008] seed s={s:.2f} q={q:.3f} alpha={alpha:.2f} -> chi2={result.fun:.2f}")
    return result


def run_fit():
    """Parallel multi-start search, then a tight-tolerance refit of the best seed."""
    print(f"[2l1s-moa19008] starting multi-start search: {len(seeds)} seeds")
    with ProcessPoolExecutor() as exec:
        results = list(exec.map(_fit_one_seed, seeds))
    best = min(results, key=lambda r: r.fun)
    print(f"[2l1s-moa19008] multi-start done, best chi2={best.fun:.2f}, refining...")

    result = minimize(chi2, x0=best.x, method="Nelder-Mead",
                       options={"xatol": 1e-6, "fatol": 1e-6, "maxiter": 20000})
    print(f"[2l1s-moa19008] final chi2={result.fun:.2f}  (dof={len(time) - 6})")
    for label, value in zip(["t0", "u0", "tE", "alpha", "s", "q"], result.x):
        print(f"[2l1s-moa19008] {label} = {value:.5f}")

    plot_fit(result.x)
    return result


def plot_fit(theta):
    """Data + 2L1S fit overlay (full baseline + auto-zoomed peak, magnitude space)
    plus the caustic-geometry panel -- same shape as fit_2l1s.py's plot_fit."""
    t0, u0, tE, alpha, s, q = theta
    A_model = binary_magnification(binary_trajectory(time, t0, u0, tE, alpha), s, q)
    fs, fb = solve_flux_calibration(A_model)

    zoom_start, zoom_end = find_zoom_window(time, -mag, mag_err, padding_fraction=0.3)

    def plot_panel(ax, xlim=None):
        window = slice(None) if xlim is None else (time >= xlim[0]) & (time <= xlim[1])
        t_grid = np.linspace(*(xlim if xlim else (time.min(), time.max())), 3000)
        A_grid = binary_magnification(binary_trajectory(t_grid, t0, u0, tE, alpha), s, q)
        mag_grid = ZERO_POINT_MAG - 2.5 * np.log10(fs * A_grid + fb)

        ax.errorbar(time, mag, yerr=mag_err, fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5, color="black", label="KMT I")
        ax.plot(t_grid, mag_grid, color="crimson", lw=1.5, label="2L1S fit")
        ax.set_ylabel("I magnitude")

        y_vals = mag[window]
        pad = 0.1 * (y_vals.max() - y_vals.min())
        ax.set_ylim(y_vals.max() + pad, y_vals.min() - pad)  # inverted: brighter (lower mag) at top
        if xlim is not None:
            ax.set_xlim(*xlim)

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
    out_path = f"scratch/{SHORT_NAME}_2l1s.png"
    fig.savefig(out_path, dpi=600)
    print(f"saved {out_path}")


if __name__ == "__main__":
    run_fit()
