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


def plot_fit_panels(ogle, moa, model_fn, fit_label, out_path):
    """Two-panel (full baseline + auto-zoomed peak) OGLE+MOA magnification plot with
    a model overlay -- shared by every O-03-BLG235 magnification-space fit script.
    model_fn(t_grid) -> A_model."""
    ogle_time, ogle_A, ogle_A_err = ogle
    moa_time, moa_A, moa_A_err = moa
    time = np.concatenate([ogle_time, moa_time])
    A_obs = np.concatenate([ogle_A, moa_A])
    zoom_start, zoom_end = find_zoom_window(moa_time, moa_A, moa_A_err, padding_fraction=0.3)

    def plot_panel(ax, xlim=None):
        t_grid = np.linspace(*(xlim if xlim else (time.min(), time.max())), 3000)
        A_model = model_fn(t_grid)
        ax.errorbar(ogle_time, ogle_A, yerr=ogle_A_err, fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5,
                    color="black", label="OGLE")
        ax.errorbar(moa_time, moa_A, yerr=moa_A_err, fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5,
                    color="tab:orange", label="MOA")
        ax.plot(t_grid, A_model, color="crimson", lw=1.5, label=fit_label)
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

    fig, (ax_full, ax_zoom) = plt.subplots(2, 1, figsize=(8, 8))
    plot_panel(ax_full)
    ax_full.set_title(f"{fit_label}, full baseline")
    ax_full.legend(loc="upper right")
    plot_panel(ax_zoom, xlim=(zoom_start, zoom_end))
    ax_zoom.set_title("zoomed on peak (auto-detected)")
    ax_zoom.set_xlabel("HJD - 2450000")
    fig.tight_layout()

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=600)
    plt.close(fig)
    print(f"saved {out_path}")
