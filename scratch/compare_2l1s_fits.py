"""Compare every 2L1S candidate solution found this session (Nelder-Mead multi-start,
tight refit, and MCMC best samples, joint + MOA-only) on one figure: light-curve
overlay (full baseline + zoomed peak) plus one caustic-geometry panel per candidate.
Bond et al. 2004's published solution is included as a reference.

Data points are shown once, converted to magnification via the tight refit's profiled
calibration (our best validated solution) -- model curves themselves are purely
geometric (A(t) doesn't depend on fs/fb) so every candidate overlays on that same
fixed backdrop regardless of which chi2 (joint or MOA-only) found it.

One-time comparison script -- lives in scratch/, run manually as
`python3 scratch/compare_2l1s_fits.py` from the project root.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 14})

from fit_2l1s import ogle_time, ogle_mag, ogle_mag_err, moa_time, moa_flux, moa_flux_err, profile_flux
from lc_models import binary_magnification, binary_trajectory, caustic_curve
from preprocess_binary_data import moa_to_magnification, ogle_to_magnification
from zoom_utils import find_zoom_window

# label, theta=(t0,u0,tE,alpha,s,q), chi2, crosses_caustic
CANDIDATES = [
    ("joint multi-start",    (2848.97616, 0.12824, 66.13835, 6.66082, 0.92931, 0.00257), 1909.02, False),
    ("joint MCMC best",      (2848.18893, 0.11076, 76.77174, 3.62710, 1.08204, 0.00632), 1932.17, True),
    ("joint tight refit",    (2847.56905, 0.09376, 82.54863, 3.65096, 1.08476, 0.00620), 1825.35, True),
    ("moa-only multi-start", (2850.48656, 0.20045, 42.22517, 6.17412, 0.90842, 0.01909), 1452.74, False),
    ("moa-only MCMC best",   (2847.20542, 0.19141, 85.50410, 1.65220, 1.96359, 0.08877), 1653.42, False),
    ("Bond et al. 2004",     (2848.06, 0.133, 61.5, np.radians(223.8), 1.120, 0.0039), 1952.90, True),
]
COLORS = plt.get_cmap("tab10").colors  # one distinct color per candidate, shared by every panel

# Fixed data calibration for display -- the tight refit, our best validated solution.
refit_theta = CANDIDATES[2][1]
t0, u0, tE, alpha, s, q = refit_theta
A_ogle_ref = binary_magnification(binary_trajectory(ogle_time, t0, u0, tE, alpha), s, q)
A_moa_ref = binary_magnification(binary_trajectory(moa_time, t0, u0, tE, alpha), s, q)
fs_ogle, fb_ogle, fs_moa = profile_flux(A_ogle_ref, A_moa_ref)
ogle_A, ogle_A_err = ogle_to_magnification(ogle_mag, ogle_mag_err, fs_ogle, fb_ogle)
moa_A, moa_A_err = moa_to_magnification(moa_flux, moa_flux_err, fs_moa)

zoom_start, zoom_end = find_zoom_window(moa_time, moa_A, moa_A_err, padding_fraction=0.3)


def plot_panel(ax, xlim=None, cap_mult=2.0):
    t_grid = np.linspace(*(xlim if xlim else (ogle_time.min(), moa_time.max())), 4000)
    ax.errorbar(ogle_time, ogle_A, yerr=ogle_A_err, fmt="+", ms=6, mew=1.2, elinewidth=1.2, capsize=3, capthick=1.2,
                color="black", alpha=0.5, zorder=1, label="OGLE")
    ax.errorbar(moa_time, moa_A, yerr=moa_A_err, fmt="+", ms=6, mew=1.2, elinewidth=1.2, capsize=3, capthick=1.2,
                color="tab:orange", alpha=0.5, zorder=1, label="MOA")

    model_y_vals = []
    for (label, theta, chi2, crosses), color in zip(CANDIDATES, COLORS):
        t0, u0, tE, alpha, s, q = theta
        A_model = binary_magnification(binary_trajectory(t_grid, t0, u0, tE, alpha), s, q)
        ls = "-" if crosses else "--"
        ax.plot(t_grid, A_model, color=color, lw=1.5, ls=ls, zorder=2,
                 label=f"{label} (chi2={chi2:.0f}{'' if crosses else ', no crossing'})")
        model_y_vals.append(A_model)  # t_grid already spans xlim, so no masking needed

    ax.set_ylabel("Magnification A(t)")
    if xlim is not None:
        ax.set_xlim(*xlim)
    in_ogle = (ogle_time >= xlim[0]) & (ogle_time <= xlim[1]) if xlim else slice(None)
    in_moa = (moa_time >= xlim[0]) & (moa_time <= xlim[1]) if xlim else slice(None)
    data_in_view = np.concatenate([ogle_A[in_ogle], moa_A[in_moa]])
    if data_in_view.size:
        # cap at cap_mult x the main (broad-hump) peak, not the caustic-crossing spike --
        # the single highest MOA point IS the real caustic spike (not noise), but it's a
        # 1-point-wide feature; the 99th percentile robustly lands on the broad hump
        # instead, which is what should set the scale.
        main_peak = np.percentile(data_in_view, 99)
        y_vals = np.concatenate([data_in_view] + model_y_vals)
        y_min, y_max = y_vals.min(), min(y_vals.max(), cap_mult * main_peak)
        pad = 0.1 * (y_max - y_min)
        ax.set_ylim(y_min - pad, y_max + pad)


fig = plt.figure(figsize=(13, 18))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.3])

ax_full = fig.add_subplot(gs[0])
plot_panel(ax_full, xlim=(2700, 3000), cap_mult=1.3)
ax_full.set_title("2L1S candidates, event season (HJD 2700-3000)")
ax_full.legend(loc="upper right", fontsize=11)

ax_zoom = fig.add_subplot(gs[1])
plot_panel(ax_zoom, xlim=(zoom_start, zoom_end))
ax_zoom.set_title("zoomed on peak (auto-detected)")
ax_zoom.set_xlabel("HJD - 2450000")

# Caustic geometry as an inset in the zoom panel's top-right corner -- that corner is
# clear of any curve (they've all decayed back down by the right edge of the window).
ax_caustic = ax_zoom.inset_axes([0.548, 0.475, 0.432, 0.495])
t_traj = np.linspace(zoom_start, zoom_end, 3000)
for (label, theta, chi2, crosses), color in zip(CANDIDATES, COLORS):
    t0, u0, tE, alpha, s, q = theta
    caustic = caustic_curve(s, q)
    if label == "moa-only MCMC best":
        caustic = caustic[np.abs(caustic.real) < 0.6]  # wide topology's 2nd island sits far from zeta=0
    traj = binary_trajectory(t_traj, t0, u0, tE, alpha)
    ax_caustic.scatter(caustic.real, caustic.imag, s=1, color=color)
    ls = "-" if crosses else "--"
    ax_caustic.plot(traj.real, traj.imag, color=color, lw=1, ls=ls)
ax_caustic.set_aspect("equal")
ax_caustic.set_xlim(-0.4, 0.4)
ax_caustic.set_ylim(-0.3, 0.3)
ax_caustic.set_title("caustic geometry", fontsize=10)
ax_caustic.tick_params(labelsize=7)

fig.tight_layout()

Path("scratch").mkdir(exist_ok=True)
out_path = "scratch/O-03-BLG235_2l1s_comparison.png"
fig.savefig(out_path, dpi=400)
print(f"saved {out_path}")
