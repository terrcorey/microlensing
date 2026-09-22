"""MCMC PSPL fits for OGLE-2003-BLG-235 / MOA-2003-BLG-53.

Both fits below use the single-lens PSPL model, even though this event is
a real binary-lens/caustic-crossing detection (the first microlensing
planet found) -- see the caustic bump in the joint fit's zoomed panel.

- run_ogle_only_diagnostic(): fits OGLE's data ALONE, no calibration. A
  deliberate sanity check showing a single sparse instrument can look
  deceptively well-fit even on a real binary-lens event.
- run_joint_fit(): fits the combined OGLE+MOA magnification data written
  by preprocess_binary_data.py -- run that first. This is the canonical
  fit for this event.

Staged via --stage: raw (raw_lc/ only, both fits skipped), quicklook (both
fits' point-estimate + fit_lc/, no MCMC), mcmc (default: full run, both
fits' corner_plots/, hist_plots/, fit_summary.dat too).
"""

from pathlib import Path

import emcee
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from lc_models import ZERO_POINT_MAG, magnification, magnitude, trajectory
from mcmc_fit import LABELS as PSPL_LABELS
from mcmc_fit import estimate_mass, fit_pspl_mcmc, plot_histograms, save_corner, save_summary
from zoom_utils import plot_fit_panels

# short name (see dataset_names.txt) -> (raw OGLE filename, u0_guess, tE_guess)
# for run_ogle_only_diagnostic(). u0/tE guesses are dataset-specific -- reusing
# another dataset's blindly reproduces the u0-sign bug from CHANGELOG.md
# (session 1): curve_fit can converge to the mirror-image (u0 -> -u0) solution
# from a bad starting point, since the model only depends on u0 squared.
DATASETS = {
    "O-03-BLG235": ("OGLE-2003-BLG-235_OGLE.tbl.txt", 0.2, 30.0),
}

SHORT_NAME = "O-03-BLG235"  # overwritten from --dataset in __main__


def plot_raw():
    """Combined OGLE+MOA calibrated magnification, no fit overlay -- the --stage=raw output."""
    ogle_time, ogle_A, ogle_A_err = np.loadtxt(f"data/processed/{SHORT_NAME}_OGLE_magnification.dat", unpack=True)
    moa_time, moa_A, moa_A_err = np.loadtxt(f"data/processed/{SHORT_NAME}_MOA_magnification.dat", unpack=True)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(ogle_time, ogle_A, yerr=ogle_A_err, fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5,
                color="black", label="OGLE")
    ax.errorbar(moa_time, moa_A, yerr=moa_A_err, fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5,
                color="tab:orange", label="MOA")
    ax.set_xlabel("HJD - 2450000")
    ax.set_ylabel("Magnification A(t)")
    ax.set_title("OGLE-2003-BLG-235 / MOA-2003-BLG-53 (raw, joint-calibrated)")
    ax.legend()
    fig.tight_layout()

    Path("raw_lc").mkdir(exist_ok=True)
    out_path = f"raw_lc/{SHORT_NAME}.png"
    fig.savefig(out_path, dpi=600)
    plt.close(fig)
    print(f"saved {out_path}")


def run_ogle_only_diagnostic(stage="mcmc"):
    ogle_file, u0_guess, tE_guess = DATASETS[SHORT_NAME]
    time, mag, mag_err = np.loadtxt(
        f"data/{ogle_file}", comments=("\\", "|"), unpack=True
    )
    time = time - 2450000.0

    best_fit, samples = fit_pspl_mcmc(time, mag, mag_err, u0_guess=u0_guess, tE_guess=tE_guess,
                                       run_mcmc=(stage == "mcmc"))

    chi2 = np.sum(((mag - magnitude(time, *best_fit)) / mag_err) ** 2)
    print(f"[ogle_only] PSPL curve_fit chi2/dof = {chi2 / (len(time) - len(best_fit)):.1f}  (bad fit expected)")

    t_model = np.linspace(time.min(), time.max(), 3000)
    fig_fit, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(time, mag, yerr=mag_err, fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5, color="black", label="OGLE data")
    ax.plot(t_model, magnitude(t_model, *best_fit), color="crimson", lw=1.5, label="PSPL fit")
    ax.invert_yaxis()
    ax.set_xlabel("HJD - 2450000")
    ax.set_ylabel("I magnitude")
    ax.set_title("OGLE-2003-BLG-235 (PSPL fit onto a real binary-lens event)")
    ax.legend()
    fig_fit.tight_layout()

    Path("fit_lc").mkdir(exist_ok=True)
    fig_fit.savefig(f"fit_lc/{SHORT_NAME}_ogle_only.png", dpi=600)
    plt.close(fig_fit)
    print(f"saved fit_lc/{SHORT_NAME}_ogle_only.png")

    if stage != "mcmc":
        return

    print(f"[ogle_only] {samples.shape[0]} posterior samples after burn-in/thinning")
    _, u0_s, tE_s, fs_s, fb_s = samples.T
    derived = {
        "A_max": magnification(u0_s),
        "t_eff": u0_s * tE_s,
        "blend_fraction": fb_s / (fs_s + fb_s),
        "m_source": ZERO_POINT_MAG - 2.5 * np.log10(fs_s),
        "M_lens": estimate_mass(tE_s),
    }
    results = {label: samples[:, i] for i, label in enumerate(PSPL_LABELS)} | derived
    save_summary(results, f"data/processed/{SHORT_NAME}_ogle_only_fit_summary.dat", prefix="[ogle_only] ")
    plot_histograms(results, f"hist_plots/{SHORT_NAME}_ogle_only.png")

    save_corner(samples, PSPL_LABELS, best_fit, f"corner_plots/{SHORT_NAME}_ogle_only.png")


def run_joint_fit(stage="mcmc"):
    labels = ["t0", "u0", "tE"]
    ogle_time, ogle_A, ogle_A_err = np.loadtxt(f"data/processed/{SHORT_NAME}_OGLE_magnification.dat", unpack=True)
    moa_time, moa_A, moa_A_err = np.loadtxt(f"data/processed/{SHORT_NAME}_MOA_magnification.dat", unpack=True)

    time = np.concatenate([ogle_time, moa_time])
    A_obs = np.concatenate([ogle_A, moa_A])
    A_err = np.concatenate([ogle_A_err, moa_A_err])

    def chi2(theta):
        t0, u0, tE = theta
        if tE <= 0:
            return np.inf
        A_model = magnification(trajectory(time, t0, u0, tE))
        return np.sum(((A_obs - A_model) / A_err) ** 2)

    t0_guess = time[np.argmax(A_obs)]
    p0_guess = [t0_guess, 0.2, 30.0]
    result = minimize(chi2, x0=p0_guess, method="Nelder-Mead",
                       options={"xatol": 1e-8, "fatol": 1e-8, "maxiter": 20000})
    best_fit = result.x
    best_fit[1], best_fit[2] = abs(best_fit[1]), abs(best_fit[2])  # u0, tE sign is arbitrary
    t0, u0, tE = best_fit

    n_dof = len(time) - 3
    print(f"[joint] combined PSPL fit: chi2/dof = {chi2(best_fit) / n_dof:.2f}")
    for label, value in zip(labels, best_fit):
        print(f"[joint] {label} = {value:.5f}")

    plot_fit_panels(
        (ogle_time, ogle_A, ogle_A_err), (moa_time, moa_A, moa_A_err),
        model_fn=lambda t: magnification(trajectory(t, t0, u0, tE)),
        fit_label="PSPL fit", out_path=f"fit_lc/{SHORT_NAME}.png",
    )

    if stage != "mcmc":
        return

    def log_prior(theta):
        t0, u0, tE = theta
        if not (time.min() < t0 < time.max()):
            return -np.inf
        if not (0 < u0 < 5):
            return -np.inf
        if not (0.1 < tE < 1000):
            return -np.inf
        return 0.0

    def log_probability(theta):
        lp = log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp - 0.5 * chi2(theta)

    ndim, nwalkers = 3, 32
    rng = np.random.default_rng(42)
    spread = np.array([0.1, 0.02, 1.0])
    p0 = best_fit + spread * rng.standard_normal((nwalkers, ndim))

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability)
    sampler.run_mcmc(p0, 4000, progress=False)

    samples = sampler.get_chain(discard=1000, thin=15, flat=True)
    print(f"[joint] {samples.shape[0]} posterior samples after burn-in/thinning")
    _, u0_s, tE_s = samples.T
    derived = {"A_max": magnification(u0_s), "t_eff": u0_s * tE_s, "M_lens": estimate_mass(tE_s)}
    results = {label: samples[:, i] for i, label in enumerate(labels)} | derived
    save_summary(results, f"data/processed/{SHORT_NAME}_fit_summary.dat", prefix="[joint] ")
    plot_histograms(results, f"hist_plots/{SHORT_NAME}.png")

    save_corner(samples, labels, best_fit, f"corner_plots/{SHORT_NAME}.png")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", default="O-03-BLG235", choices=sorted(DATASETS),
                         help="short name from dataset_names.txt. Must already have processed "
                              "magnification files from preprocess_binary_data.py.")
    parser.add_argument("--stage", choices=["raw", "quicklook", "mcmc"], default="mcmc",
                         help="raw: raw_lc/ only. quicklook: + fit_lc/ for both fits (no MCMC). "
                              "mcmc (default): full run, + corner_plots/, hist_plots/, fit_summary.dat.")
    args = parser.parse_args()

    SHORT_NAME = args.dataset

    if args.stage == "raw":
        plot_raw()
    else:
        run_ogle_only_diagnostic(stage=args.stage)
        run_joint_fit(stage=args.stage)
