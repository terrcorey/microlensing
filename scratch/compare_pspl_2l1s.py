"""PSPL vs 2L1S model comparison (BIC) for OGLE-2003-BLG-235/MOA-2003-BLG-53.

Both models are scored on the same raw OGLE+MOA flux, the same rescaled errors
(fit_2l1s's K_OGLE/K_MOA) and the same profiled flux calibration
(fit_2l1s.flux_residuals), all parallax-free to match the Cassan 2L1S fit.
plot_comparison() overlays both best fits on one light curve.
Run as `python3 scratch/compare_pspl_2l1s.py` from the project root.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

plt.rcParams.update({"font.size": 14})

from fit_2l1s import (flux_residuals, plain_fit, plain_chi2_fn, profile_flux,
                      ogle_time, ogle_mag, ogle_mag_err, moa_time, moa_flux, moa_flux_err)
from lc_models import binary_magnification_fs, binary_trajectory, magnification, trajectory
from preprocess_binary_data import moa_to_magnification, ogle_to_magnification
from zoom_utils import find_zoom_window
from cassan_caustic import OUT_DIR, SHORT_NAME, CassanParams, chi2_cassan, to_standard

def run_pspl():
    def chi2(theta):
        t0, u0, tE = theta
        A_ogle = magnification(trajectory(ogle_time, t0, u0, tE, 0, 0, 0, 0))
        A_moa = magnification(trajectory(moa_time, t0, u0, tE, 0, 0, 0, 0))
        return np.sum(flux_residuals(A_ogle, A_moa) ** 2)
    result = minimize(chi2, x0=plain_fit[:3], method="Nelder-Mead",
                    options={"xatol": 1e-6, "fatol": 1e-6, "maxiter": 20000})
    return result.x, result.fun

def run_2l1s():
    chain = np.load(OUT_DIR / f"{SHORT_NAME}_2l1s_cassan_mcmc_chain.npz")
    best_sample = chain["samples"][np.argmax(chain["log_probs"])]
    step = np.std(chain["samples"], axis = 0)
    simplex = np.vstack([best_sample, best_sample + np.diag(step)])
    result = minimize(chi2_cassan, x0=best_sample, method="Nelder-Mead",
                      options={"xatol": 1e-3, "fatol": 1e-3, "maxiter": 5000, 
                               "initial_simplex": simplex})
    return result.x, result.fun

def plot_comparison(pspl_theta, cassan_theta):
    """Both best fits over the data: event season, zoomed peak, then one raw-residual row
    per model. Data are converted to magnification with the 2L1S fit's profiled fs/fb.
    PSPL's own fs/fb differ, so its predicted *flux* is re-expressed on that same scale --
    one curve per instrument, since OGLE's blend term makes the two mappings differ."""
    t0, u0, tE = pspl_theta
    std = to_standard(cassan_theta)
    models = {
        "PSPL": lambda t: magnification(trajectory(t, t0, u0, tE, 0, 0, 0, 0)),
        "2L1S": lambda t: binary_magnification_fs(binary_trajectory(t, *std[:4], 0, 0, 0, 0), std.s, std.q, std.rho),
    }
    cal = {name: profile_flux(A(ogle_time), A(moa_time)) for name, A in models.items()}
    fs_ogle, fb_ogle, fs_moa = cal["2L1S"]
    ogle_A, ogle_A_err = ogle_to_magnification(ogle_mag, ogle_mag_err, fs_ogle, fb_ogle)
    moa_A, moa_A_err = moa_to_magnification(moa_flux, moa_flux_err, fs_moa)

    def predicted(name, t, instrument):
        """Model `name`'s predicted flux at t, in the 2L1S calibration's magnification units."""
        fs_o, fb_o, fs_m = cal[name]
        A = models[name](t)
        return (fs_o * A + fb_o - fb_ogle) / fs_ogle if instrument == "ogle" else 1 + fs_m / fs_moa * (A - 1)

    data = {"ogle": (ogle_time, ogle_A, ogle_A_err, "black"), "moa": (moa_time, moa_A, moa_A_err, "tab:orange")}
    curves = [("2L1S", "moa", "crimson", "-"), ("PSPL", "ogle", "tab:blue", "--"), ("PSPL", "moa", "tab:cyan", "--")]
    style = dict(fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5)
    zoom = find_zoom_window(moa_time, moa_A, moa_A_err, padding_fraction=0.3)

    fig, axes = plt.subplots(4, 1, figsize=(10.5, 14), height_ratios=[4, 4, 1.5, 1.5])
    for ax, xlim in zip(axes[:2], [(2700, 3000), zoom]):
        t_grid = np.linspace(*xlim, 3000)
        y_vals = []
        for inst, (t, A, A_err, color) in data.items():
            ax.errorbar(t, A, yerr=A_err, color=color, label=inst.upper(), **style)
            y_vals.append(A[(t >= xlim[0]) & (t <= xlim[1])])
        for name, inst, color, ls in curves:
            curve = predicted(name, t_grid, inst)  # 2L1S maps to itself, so one curve covers both instruments
            label = name if name == "2L1S" else f"{name} ({inst.upper()} calibration)"
            ax.plot(t_grid, curve, color=color, ls=ls, lw=1.5, label=label)
            y_vals.append(curve)
        y_vals = np.concatenate(y_vals)
        pad = 0.1 * np.ptp(y_vals)
        ax.set(xlim=xlim, ylim=(y_vals.min() - pad, y_vals.max() + pad), ylabel="Magnification A(t)")
    axes[0].set_title("PSPL vs 2L1S, event season (HJD 2700-3000)")
    axes[0].legend(loc="upper right", fontsize=10)
    axes[1].set_title("zoomed on peak (auto-detected)")

    # raw residuals (data - model) in the same A units, real per-instrument error bars --
    # plot_fit()'s convention; shared y-limits so the two models' misfits compare directly
    resid_lim = 0
    for ax, name in zip(axes[2:], ["PSPL", "2L1S"]):
        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
        for inst, (t, A, A_err, color) in data.items():
            in_zoom = (t >= zoom[0]) & (t <= zoom[1])
            resid = A[in_zoom] - predicted(name, t[in_zoom], inst)
            ax.errorbar(t[in_zoom], resid, yerr=A_err[in_zoom], color=color, **style)
            resid_lim = max(resid_lim, np.abs(resid).max() * 1.1)
        ax.set(xlim=zoom, ylabel=f"{name} resid")
    for ax in axes[2:]:
        ax.set_ylim(-resid_lim, resid_lim)
    axes[3].set_xlabel("HJD - 2450000")

    fig.tight_layout()
    out_dir = OUT_DIR.parent / "compare"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{SHORT_NAME}_pspl_vs_2l1s.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"saved {out_path}")

if __name__ == "__main__":
    best_fit, chi2 = run_pspl()
    t0, u0, tE = best_fit
    best_fit_2l1s, cassan_chi2 = run_2l1s()
    print(f"PSPL chi2 (flux space) = {chi2:.2f} | t0 = {t0:.5f}, u0 = {u0:.5f}, tE = {tE:.5f}")
    print(f"fit_joint_pspl chi2 (OGLE in mag space) at plain_fit = {plain_chi2_fn(plain_fit):.2f}")
    print(f"2L1S chi2 (flux space) = {cassan_chi2:.2f}")
    for label, value in zip(CassanParams._fields, best_fit_2l1s):
        print(f"    {label} = {value:.5f}")

    n = len(ogle_time) + len(moa_time)
    k_pspl, k_2l1s = 6, 10
    bic_pspl = chi2 + k_pspl * np.log(n)
    bic_2l1s = cassan_chi2 + k_2l1s * np.log(n)
    print(f"N = {n} | BIC_PSPL = {bic_pspl:.2f}, BIC_2L1S = {bic_2l1s:.2f}, "
          f"Delta BIC (PSPL - 2L1S, > 0 favours 2L1S) = {bic_pspl - bic_2l1s:.2f}")

    plot_comparison(best_fit, best_fit_2l1s)
