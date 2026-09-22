"""Automatically find the time window around a light curve's peak, so plots
can zoom in without hardcoding a window per dataset."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def find_zoom_window(time, brightness, error, sigma_threshold=5.0, padding_fraction=0.2):
    """Window around points whose excess over baseline is >= sigma_threshold * error (significance, not raw amplitude -- see CHANGELOG)."""
    baseline = np.median(brightness)
    excess = brightness - baseline
    in_event = excess >= sigma_threshold * error
    event_times = time[in_event]
    start, end = event_times.min(), event_times.max()
    pad = padding_fraction * (end - start)
    return start - pad, end + pad


def plot_residual_panel(ax, series, invert=True):
    """Standardized-residual scatter (yerr=1, since residuals are already
    divided by their own error) with a +-1 sigma translucent band and a
    symmetric y-axis around 0 -- shared by every fit-overlay plot.

    series: list of (x, residuals[, color]) tuples, one per instrument/dataset
    sharing this panel. invert=True matches a magnitude panel above it
    (brighter/lower-mag residuals plot upward); pass False for a
    magnification-space panel, where higher values already plot upward.
    """
    all_residuals = np.concatenate([np.asarray(s[1]) for s in series])
    ylim = max(np.abs(all_residuals).max() * 1.1, 1.5) if all_residuals.size else 1.5
    ax.axhspan(-1, 1, color="gray", alpha=0.15, linewidth=0)
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    for x, residuals, *color in series:
        ax.errorbar(x, residuals, yerr=1, fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5,
                     color=color[0] if color else None)
    ax.set_ylim(-ylim, ylim)
    if invert:
        ax.invert_yaxis()
    ax.set_ylabel("residual (σ)")


def plot_residual_hist(ax, series):
    """Residual histogram rotated 90deg (count on x, residual (sigma) on y) so
    it sits directly beside its matching plot_residual_panel() scatter,
    sharing that panel's y-axis (a sharey ax needs no ylim/invert of its own).

    series: list of residuals arrays, matching plot_residual_panel()'s series.
    """
    all_residuals = np.concatenate([np.asarray(r) for r in series])
    ax.axhspan(-1, 1, color="gray", alpha=0.15, linewidth=0)
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax.hist(all_residuals, bins=30, orientation="horizontal", color="steelblue", edgecolor="white")
    ax.tick_params(labelleft=False)
    ax.set_xlabel("count")


def plot_fit_panels(ogle, moa, model_fn, fit_label, out_path):
    """Two-panel (full baseline + auto-zoomed peak) OGLE+MOA magnification plot with
    a model overlay, each paired with a RAW (non-standardized) residual scatter
    panel -- obs-model in physical magnification units, each point keeping its
    own instrument's real error bar -- and a rotated residual histogram sharing
    that panel's y-axis. model_fn(t_grid) -> A_model.

    Deliberately NOT built on plot_residual_panel()/plot_residual_hist() (those
    show standardized (obs-model)/err residuals, where every point's error bar
    is 1 by construction): with two instruments of genuinely different
    precision on one panel, a uniform yerr=1 hides that difference. Raw
    residuals with each instrument's real error bar show it directly -- see
    CHANGELOG. A real structural difference from plot_fit_lc()'s
    single-instrument case, not a candidate for forcing through one shared
    helper (see CLAUDE.md's "don't force an abstraction over a real
    difference" note).
    """
    ogle_time, ogle_A, ogle_A_err = ogle
    moa_time, moa_A, moa_A_err = moa
    time = np.concatenate([ogle_time, moa_time])
    A_obs = np.concatenate([ogle_A, moa_A])
    zoom_start, zoom_end = find_zoom_window(moa_time, moa_A, moa_A_err, padding_fraction=0.3)

    ogle_resid = ogle_A - model_fn(ogle_time)
    moa_resid = moa_A - model_fn(moa_time)

    def plot_panel(ax, ax_resid, ax_hist, xlim=None):
        t_grid = np.linspace(*(xlim if xlim else (time.min(), time.max())), 3000)
        A_model = model_fn(t_grid)
        ax.errorbar(ogle_time, ogle_A, yerr=ogle_A_err, fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5,
                    color="black", label="OGLE")
        ax.errorbar(moa_time, moa_A, yerr=moa_A_err, fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5,
                    color="tab:orange", label="MOA")
        ax.plot(t_grid, A_model, color="crimson", lw=1.5, label=fit_label)
        ax.set_ylabel("Magnification A(t)")
        ax.sharex(ax_resid)

        in_zoom_ogle = (ogle_time >= xlim[0]) & (ogle_time <= xlim[1]) if xlim else np.ones_like(ogle_time, dtype=bool)
        in_zoom_moa = (moa_time >= xlim[0]) & (moa_time <= xlim[1]) if xlim else np.ones_like(moa_time, dtype=bool)
        if xlim is not None:
            ax.set_xlim(*xlim)
            in_zoom = (time >= xlim[0]) & (time <= xlim[1])
            # include the model curve, not just the data -- a sharp model peak/dip
            # between data points would otherwise get clipped by the y-limits.
            y_vals = np.concatenate([A_obs[in_zoom], A_model])
            if y_vals.size:
                pad = 0.1 * (y_vals.max() - y_vals.min())
                ax.set_ylim(y_vals.min() - pad, y_vals.max() + pad)
            ax_resid.set_xlim(*xlim)

        ogle_resid_z, moa_resid_z = ogle_resid[in_zoom_ogle], moa_resid[in_zoom_moa]
        ax_resid.axhline(0, color="gray", linestyle="--", linewidth=0.8)
        ax_resid.errorbar(ogle_time[in_zoom_ogle], ogle_resid_z, yerr=ogle_A_err[in_zoom_ogle],
                           fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5, color="black")
        ax_resid.errorbar(moa_time[in_zoom_moa], moa_resid_z, yerr=moa_A_err[in_zoom_moa],
                           fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5, color="tab:orange")
        all_resid_z = np.concatenate([ogle_resid_z, moa_resid_z])
        if all_resid_z.size:
            ylim = np.abs(all_resid_z).max() * 1.1
            ax_resid.set_ylim(-ylim, ylim)
        ax_resid.set_ylabel("residual (A(t))")
        ax_resid.set_xlabel("HJD - 2450000")

        ax_hist.axhline(0, color="gray", linestyle="--", linewidth=0.8)
        ax_hist.hist(all_resid_z, bins=30, orientation="horizontal", color="steelblue", edgecolor="white")
        ax_hist.tick_params(labelleft=False)
        ax_hist.set_xlabel("count")

    fig = plt.figure(figsize=(9, 13))
    gs = fig.add_gridspec(4, 2, height_ratios=[3, 2, 3, 2], width_ratios=[4, 1], hspace=0.6, wspace=0.05)
    ax_full = fig.add_subplot(gs[0, 0])
    ax_full_resid = fig.add_subplot(gs[1, 0])
    ax_full_hist = fig.add_subplot(gs[1, 1], sharey=ax_full_resid)
    ax_zoom = fig.add_subplot(gs[2, 0])
    ax_zoom_resid = fig.add_subplot(gs[3, 0])
    ax_zoom_hist = fig.add_subplot(gs[3, 1], sharey=ax_zoom_resid)

    plot_panel(ax_full, ax_full_resid, ax_full_hist)
    ax_full.set_title(f"{fit_label}, full baseline")
    ax_full.legend(loc="upper right")

    plot_panel(ax_zoom, ax_zoom_resid, ax_zoom_hist, xlim=(zoom_start, zoom_end))
    ax_zoom.set_title("zoomed on peak (auto-detected)")

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=600)
    plt.close(fig)
    print(f"saved {out_path}")
