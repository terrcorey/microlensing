from pathlib import Path

import corner
import emcee
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from lc_models import ZERO_POINT_MAG, flux, magnification, magnitude
from zoom_utils import find_zoom_window

LABELS = ["t0", "u0", "tE", "f_source", "f_blend"]


def fit_pspl_mcmc(time, mag, mag_err, u0_guess=0.5, tE_guess=50.0, run_mcmc=True):
    """Fit the 5-param PSPL model: curve_fit for a starting point, then
    emcee. Returns (best_fit, samples) -- samples is the flat, burned-in,
    thinned chain, or None if run_mcmc=False (quicklook: point estimate only).

    u0_guess/tE_guess matter: curve_fit can converge to the mirror-image
    (u0 -> -u0) solution from a bad starting point, since the model only
    depends on u0 squared. Tune these per dataset if the fit looks off.
    """
    t0_guess = time[np.argmin(mag)]
    f_source_guess = 10 ** (-0.4 * (np.median(mag) - ZERO_POINT_MAG))
    p0_guess = [t0_guess, u0_guess, tE_guess, f_source_guess, 0.0]
    best_fit, _ = curve_fit(magnitude, time, mag, sigma=mag_err, p0=p0_guess)
    if not run_mcmc:
        return best_fit, None

    def log_prior(theta):
        t0, u0, tE, f_source, f_blend = theta
        if not (time.min() < t0 < time.max()):
            return -np.inf
        if not (0 < u0 < 5):
            return -np.inf
        if not (0.1 < tE < 1000):
            return -np.inf
        if not (f_source > 0):
            return -np.inf
        return 0.0

    def log_likelihood(theta):
        model_flux = flux(time, *theta)
        if np.any(model_flux <= 0):
            return -np.inf
        model_mag = ZERO_POINT_MAG - 2.5 * np.log10(model_flux)
        return -0.5 * np.sum(((mag - model_mag) / mag_err) ** 2)

    def log_probability(theta):
        lp = log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + log_likelihood(theta)

    ndim, nwalkers = 5, 32
    rng = np.random.default_rng(42)
    spread = np.array([0.1, 0.02, 1.0, 0.01 * best_fit[3], 0.01 * max(abs(best_fit[3]), 1e-3)])
    p0 = best_fit + spread * rng.standard_normal((nwalkers, ndim))

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability)
    sampler.run_mcmc(p0, 4000, progress=False)
    return best_fit, sampler.get_chain(discard=1000, thin=15, flat=True)


KAPPA = 8.144  # mas / Msun; theta_E[mas]^2 = KAPPA * M[Msun] * pi_rel[mas]


def estimate_mass(tE_days, D_S_kpc=8.0, mu_rel_median=4.0, seed=42):
    """Monte Carlo lens mass from tE alone via assumed priors on relative
    proper motion (mu_rel) and lens distance (D_L) -- tE alone can't
    separate mass from distance, so this pairs each posterior tE sample
    with one prior draw to show the resulting mass *distribution*, not a
    point estimate.

    D_S fixed at the standard Galactic bulge distance (both datasets are
    OGLE bulge fields). mu_rel: log-normal around a typical bulge-lensing
    value. D_L: uniform along the line of sight to the source.
    """
    # ponytail: D_L uniform + fixed D_S is a naive line-of-sight prior, not
    # a Galactic density model; upgrade to a disk/bulge-weighted prior
    # (e.g. Han & Gould 2003) if the mass distribution's shape needs to be
    # quantitatively trustworthy rather than illustrative.
    rng = np.random.default_rng(seed)
    n = len(tE_days)
    mu_rel = rng.lognormal(mean=np.log(mu_rel_median), sigma=0.4, size=n)  # mas/yr
    D_L_kpc = rng.uniform(0.1, D_S_kpc, size=n)
    pi_rel = 1.0 / D_L_kpc - 1.0 / D_S_kpc  # mas
    theta_E = (np.asarray(tE_days) / 365.25) * mu_rel  # mas
    return theta_E**2 / (KAPPA * pi_rel)  # Msun


def plot_histograms(results, out_path, bins=50):
    """1D histogram per entry in results -- a percentile triplet alone can
    hide bimodality (e.g. the u0 mirror-image degeneracy); this shows the
    actual shape."""
    ncols = min(len(results), 3)
    nrows = -(-len(results) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows), squeeze=False)
    axes = axes.flatten()
    for ax, (label, values) in zip(axes, results.items()):
        values = np.asarray(values)
        if values.min() > 0 and values.max() / values.min() > 100:
            # heavy-tailed (e.g. M_lens, whose prior can send pi_rel -> 0):
            # linear bins crush the whole shape into one bin near zero.
            ax.hist(values, bins=np.logspace(np.log10(values.min()), np.log10(values.max()), bins), color="steelblue")
            ax.set_xscale("log")
        else:
            ax.hist(values, bins=bins, color="steelblue")
        ax.set_title(label)
    for ax in axes[len(results):]:
        ax.axis("off")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"saved {out_path}")


def plot_raw(time, mag, mag_err, title, out_path):
    """Raw light curve, no fit overlay -- the --stage=raw output."""
    fig, ax = plt.subplots()
    ax.errorbar(time, mag, yerr=mag_err, fmt="+", ms=3, elinewidth=0.5)
    ax.invert_yaxis()
    ax.set_xlabel("HJD - 2450000")
    ax.set_ylabel("I magnitude")
    ax.set_title(title)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"saved {out_path}")


def plot_fit_lc(time, mag, mag_err, best_fit, dataset_label, out_path):
    """Two-panel (full baseline + auto-zoomed peak) magnitude-space fit overlay
    -- the --stage=quicklook/mcmc output."""
    t_model = np.linspace(time.min(), time.max(), 2000)
    mag_model = magnitude(t_model, *best_fit)

    zoom_start, zoom_end = find_zoom_window(time, -mag, mag_err)
    in_zoom = (time >= zoom_start) & (time <= zoom_end)
    t_model_zoom = np.linspace(zoom_start, zoom_end, 2000)
    mag_model_zoom = magnitude(t_model_zoom, *best_fit)

    fig, (ax_full, ax_zoom) = plt.subplots(2, 1, figsize=(8, 8))

    ax_full.errorbar(time, mag, yerr=mag_err, fmt="+", ms=3, elinewidth=0.5, label="data")
    ax_full.plot(t_model, mag_model, color="crimson", label="PSPL fit")
    ax_full.invert_yaxis()
    ax_full.set_ylabel("I magnitude")
    ax_full.set_title(f"{dataset_label} (PSPL fit, full baseline)")
    ax_full.legend()

    ax_zoom.errorbar(time, mag, yerr=mag_err, fmt="+", ms=3, elinewidth=0.5, label="data")
    ax_zoom.plot(t_model_zoom, mag_model_zoom, color="crimson", label="PSPL fit")
    ax_zoom.set_xlim(zoom_start, zoom_end)
    ax_zoom.set_ylim(mag[in_zoom].max() + 0.05, mag[in_zoom].min() - 0.05)
    ax_zoom.set_xlabel("HJD - 2450000")
    ax_zoom.set_ylabel("I magnitude")
    ax_zoom.set_title("zoomed on peak (auto-detected)")
    fig.tight_layout()

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"saved {out_path}")


def save_corner(samples, labels, truths, out_path):
    """Corner plot of the raw MCMC parameters -- shared by every fit script.

    plot_datapoints=False + fill_contours drop the raw-point/flat-center-block
    look in favor of a smoothed, graduated-alpha density gradient.
    """
    fig = corner.corner(samples, labels=labels, truths=truths, bins=30, smooth=1.0,
                         plot_datapoints=False, fill_contours=True, color="mediumblue")
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"saved {out_path}")


def save_summary(results, out_path, prefix=""):
    """Print each param's 16/50/84 percentiles and write them to out_path."""
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    label_width = max(len(label) for label in results) + 2
    with open(out_path, "w") as f:
        f.write(f"{'# param':<{label_width}}{'p16':>14}{'p50':>14}{'p84':>14}\n")
        for label, values in results.items():
            lo, mid, hi = np.percentile(values, [16, 50, 84])
            print(f"{prefix}{label} = {mid:.5f} (+{hi - mid:.5f} / -{mid - lo:.5f})")
            f.write(f"{label:<{label_width}}{lo:>14.6f}{mid:>14.6f}{hi:>14.6f}\n")
    print(f"{prefix}saved {out_path}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["raw", "quicklook", "mcmc"], default="mcmc",
                         help="raw: raw_lc/ only. quicklook: + fit_lc/ (curve_fit, no MCMC). "
                              "mcmc (default): full run, + corner_plots/, hist_plots/, fit_summary.dat.")
    args = parser.parse_args()

    SHORT_NAME = "O-05-BLG086"  # see dataset_names.txt
    time, mag, mag_err = np.loadtxt("data/OGLE-2005-BLG-086.dat", unpack=True)

    if args.stage == "raw":
        plot_raw(time, mag, mag_err, "OGLE-2005-BLG-086 (raw)", f"raw_lc/{SHORT_NAME}.png")
    else:
        best_fit, samples = fit_pspl_mcmc(time, mag, mag_err, run_mcmc=(args.stage == "mcmc"))

        chi2 = np.sum(((mag - magnitude(time, *best_fit)) / mag_err) ** 2)
        print(f"chi2/dof = {chi2 / (len(time) - len(best_fit)):.3f}")
        for label, value in zip(LABELS, best_fit):
            print(f"{label} = {value:.5f}")
        plot_fit_lc(time, mag, mag_err, best_fit, "OGLE-2005-BLG-086", f"fit_lc/{SHORT_NAME}.png")

        if args.stage == "mcmc":
            print(f"{samples.shape[0]} posterior samples after burn-in/thinning")
            _, u0_s, tE_s, fs_s, fb_s = samples.T
            derived = {
                "A_max": magnification(u0_s),
                "t_eff": u0_s * tE_s,
                "blend_fraction": fb_s / (fs_s + fb_s),
                "m_source": ZERO_POINT_MAG - 2.5 * np.log10(fs_s),
                "M_lens": estimate_mass(tE_s),
            }
            results = {label: samples[:, i] for i, label in enumerate(LABELS)} | derived
            save_summary(results, f"data/processed/{SHORT_NAME}_fit_summary.dat")
            plot_histograms(results, f"hist_plots/{SHORT_NAME}.png")

            save_corner(samples, LABELS, best_fit, f"corner_plots/{SHORT_NAME}.png")
