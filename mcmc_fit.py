from pathlib import Path
from typing import NamedTuple

import corner
import emcee
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import t as student_t
from astropy.coordinates import SkyCoord

from lc_models import ZERO_POINT_MAG, flux, plain_flux, magnification, magnitude, plain_magnitude, sun_earth_projection
from zoom_utils import find_zoom_window


class FitResult(NamedTuple):
    """One fit's (best_fit, samples, labels), kept together so a caller
    can never pass a mismatched pair -- .column(name) is the only way to
    pull a column out, and it always reads its own samples/labels rather
    than ones handed in separately. samples is None in quicklook mode
    (run_mcmc=False)."""
    best_fit: np.ndarray
    samples: np.ndarray | None
    labels: list

    def column(self, name):
        return self.samples[:, self.labels.index(name)]


def fit_pspl_mcmc(time, mag, mag_err, u0_guess=0.5, tE_guess=50.0, run_mcmc=True):
    """Fit the 5-param PSPL model: curve_fit for a starting point, then
    emcee. Returns a FitResult -- .samples is the flat, burned-in, thinned
    chain, or None if run_mcmc=False (quicklook: point estimate only).

    u0_guess/tE_guess matter: curve_fit can converge to the mirror-image
    (u0 -> -u0) solution from a bad starting point, since the model only
    depends on u0 squared. Tune these per dataset if the fit looks off.
    """
    labels = ["t0", "u0", "tE", "f_source", "f_blend"]
    t0_guess = time[np.argmin(mag)]
    f_source_guess = 10 ** (-0.4 * (np.median(mag) - ZERO_POINT_MAG))
    p0_guess = [t0_guess, u0_guess, tE_guess, f_source_guess, 0.0]
    best_fit, _ = curve_fit(plain_magnitude, time, mag, sigma=mag_err, p0=p0_guess)
    if not run_mcmc:
        return FitResult(best_fit, None, labels)

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
        model_flux = plain_flux(time, *theta)
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
    return FitResult(best_fit, sampler.get_chain(discard=1000, thin=15, flat=True), labels)


def fit_parallax_pspl_mcmc(time, mag, mag_err, delta_sN, delta_sE, u0_guess=0.5, tE_guess=50.0, piE_N_guess=0.0, piE_E_guess=0.0, run_mcmc=True):
    """Fit the 7-param parallax PSPL model against a Student-t robust
    likelihood (2 more free params: scale, dof -- replaces the Gaussian
    likelihood outright, motivated by session 7's residual diagnostic: a
    globally heavier-than-Gaussian error distribution, not a small outlier
    population). curve_fit finds a starting point for the 7 light-curve
    params only (least-squares has no scale/dof analog), then emcee samples
    all 9. Returns a FitResult -- .best_fit is the 7-param curve_fit point
    estimate, .samples is the flat, burned-in, thinned 9-column chain, or
    None if run_mcmc=False (quicklook: point estimate only).

    u0_guess/tE_guess matter: curve_fit can converge to the mirror-image
    (u0 -> -u0) solution from a bad starting point, since the model only
    depends on u0 squared. Tune these per dataset if the fit looks off.
    """
    labels = ["t0", "u0", "tE", "f_source", "f_blend", "piE_N", "piE_E", "scale", "dof"]
    t0_guess = time[np.argmin(mag)]
    f_source_guess = 10 ** (-0.4 * (np.median(mag) - ZERO_POINT_MAG))
    def _model_magnitude(t, t0, u0, tE, f_source, f_blend, piE_N, piE_E):
        return magnitude(t, t0, u0, tE, f_source, f_blend, piE_N, piE_E, delta_sN, delta_sE)
    best_fit, _ = curve_fit(_model_magnitude, time, mag, sigma=mag_err,
                            p0=[t0_guess, u0_guess, tE_guess, f_source_guess, 0.0, piE_N_guess, piE_E_guess])
    if not run_mcmc:
        return FitResult(best_fit, None, labels)

    def log_prior(theta):
        t0, u0, tE, f_source, f_blend, piE_N, piE_E, scale, dof = theta
        if not (time.min() < t0 < time.max()):
            return -np.inf
        if not (0 < u0 < 5):
            return -np.inf
        if not (0.1 < tE < 1000):
            return -np.inf
        if not (f_source > 0):
            return -np.inf
        if not (np.abs(piE_N) < 2):
            return -np.inf
        if not (np.abs(piE_E) < 2):
            return -np.inf
        if not (0.1 < scale < 10.0):
            return -np.inf
        if not (0.5 < dof < 50.0):
            return -np.inf
        return 0.0

    def log_likelihood(theta):
        t0, u0, tE, f_source, f_blend, piE_N, piE_E, scale, dof = theta
        model_flux = flux(time, t0, u0, tE, f_source, f_blend, piE_N, piE_E, delta_sN, delta_sE)
        if np.any(model_flux <= 0):
            return -np.inf
        model_mag = ZERO_POINT_MAG - 2.5 * np.log10(model_flux)
        standardized_resid = (mag - model_mag) / mag_err
        return np.sum(student_t.logpdf(standardized_resid / scale, df=dof) - np.log(scale))

    def log_probability(theta):
        lp = log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + log_likelihood(theta)

    ndim, nwalkers = 9, 32
    rng = np.random.default_rng(42)
    spread = np.array([0.1, 0.02, 1.0, 0.01 * best_fit[3], 0.01 * max(abs(best_fit[3]), 1e-3), 0.5, 0.5, 0.1, 1.0])
    p0_center = np.concatenate([best_fit, [1.0, 5.0]])  # scale=1, dof=5 initial guesses
    p0 = p0_center + spread * rng.standard_normal((nwalkers, ndim))

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability)
    sampler.run_mcmc(p0, 4000, progress=False)
    return FitResult(best_fit, sampler.get_chain(discard=1000, thin=15, flat=True), labels)


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


def get_t0_par(time, mag, mag_err, cache_path, u0_guess, tE_guess):
    """Preliminary non-parallax t0, frozen as the reference epoch
    sun_earth_projection() anchors its geocentric frame to (Gould 2004
    convention -- not itself a fitted parameter of the real parallax fit).

    Always the MCMC posterior median, whether read from an existing
    cache_path or freshly fit -- consistent with every other parameter
    this codebase reports (save_summary()'s own p16/p50/p84 convention),
    and avoids the drift this used to have between two different
    estimators depending on whether cache_path already existed.

    u0_guess/tE_guess have no default: silently reusing one dataset's
    guess for another previously produced a real, hard-to-notice bug (the
    u0 mirror-image solution -- see mcmc_fit_binary.py's DATASETS dict).
    """
    try:
        data = np.genfromtxt(cache_path, dtype=None, names=True, encoding=None)
        return float(data["p50"][data["param"] == "t0"][0])
    except FileNotFoundError:
        pass

    fit_result = fit_pspl_mcmc(time, mag, mag_err, u0_guess=u0_guess, tE_guess=tE_guess, run_mcmc=True)
    u0_s, tE_s, fs_s, fb_s = (fit_result.column(name) for name in ("u0", "tE", "f_source", "f_blend"))
    derived = {
        "A_max": magnification(u0_s),
        "t_eff": u0_s * tE_s,
        "blend_fraction": fb_s / (fs_s + fb_s),
        "m_source": ZERO_POINT_MAG - 2.5 * np.log10(fs_s),
        "M_lens": estimate_mass(tE_s),
    }
    results = {label: fit_result.column(label) for label in fit_result.labels} | derived
    save_summary(results, cache_path)
    return float(np.percentile(fit_result.column("t0"), 50))


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
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"saved {out_path}")


def plot_raw(time, mag, mag_err, title, out_path):
    """Raw light curve, no fit overlay -- the --stage=raw output."""
    fig, ax = plt.subplots()
    ax.errorbar(time, mag, yerr=mag_err, fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5)
    ax.invert_yaxis()
    ax.set_xlabel("HJD - 2450000")
    ax.set_ylabel("I magnitude")
    ax.set_title(title)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"saved {out_path}")


def _plot_residual_panel(ax, x, residuals):
    """Standardized-residual scatter (yerr=1, since residuals are already
    divided by mag_err) with a +-1 sigma translucent band and a symmetric
    y-axis around 0 -- shared by both panels in plot_fit_lc(). Y-axis
    inverted to match the magnitude panel above it (brighter/lower-mag
    residuals plot upward here too, instead of the opposite direction)."""
    ylim = max(np.abs(residuals).max() * 1.1, 1.5)
    ax.axhspan(-1, 1, color="gray", alpha=0.15, linewidth=0)
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax.errorbar(x, residuals, yerr=1, fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5)
    ax.set_ylim(-ylim, ylim)
    ax.invert_yaxis()
    ax.set_ylabel("residual (σ)")


def _plot_residual_hist(ax, residuals):
    """Residual histogram rotated 90deg (count on x, residual (sigma) on y) so
    it sits directly beside its matching _plot_residual_panel scatter,
    sharing that panel's y-axis (sharey in plot_fit_lc() propagates its
    limits/inversion, so this needs no ylim/invert of its own)."""
    ax.axhspan(-1, 1, color="gray", alpha=0.15, linewidth=0)
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax.hist(residuals, bins=30, orientation="horizontal", color="steelblue", edgecolor="white")
    ax.tick_params(labelleft=False)
    ax.set_xlabel("count")


def plot_fit_lc(time, mag, mag_err, best_fit, coords, t0_par, dataset_label, out_path):
    """Two-panel (full baseline + auto-zoomed peak) magnitude-space fit overlay,
    each with a standardized-residual ((data-model)/mag_err) panel underneath
    and a rotated residual histogram beside that panel, sharing its y-axis --
    the --stage=quicklook/mcmc output."""

    t_model = np.linspace(time.min(), time.max(), 2000)
    delta_sN_model, delta_sE_model = sun_earth_projection(t_model, coords, t0_par)
    mag_model = magnitude(t_model, *best_fit, delta_sN_model, delta_sE_model)

    zoom_start, zoom_end = find_zoom_window(time, -mag, mag_err)
    in_zoom = (time >= zoom_start) & (time <= zoom_end)
    t_model_zoom = np.linspace(zoom_start, zoom_end, 2000)
    delta_sN_zoom, delta_sE_zoom = sun_earth_projection(t_model_zoom, coords, t0_par)
    mag_model_zoom = magnitude(t_model_zoom, *best_fit, delta_sN_zoom, delta_sE_zoom)

    delta_sN_data, delta_sE_data = sun_earth_projection(time, coords, t0_par)
    residuals = (mag - magnitude(time, *best_fit, delta_sN_data, delta_sE_data)) / mag_err

    fig = plt.figure(figsize=(9, 13))
    gs = fig.add_gridspec(4, 2, height_ratios=[3, 2, 3, 2], width_ratios=[4, 1], hspace=0.6, wspace=0.05)
    ax_full = fig.add_subplot(gs[0, 0])
    ax_full_resid = fig.add_subplot(gs[1, 0])
    ax_full_hist = fig.add_subplot(gs[1, 1], sharey=ax_full_resid)
    ax_zoom = fig.add_subplot(gs[2, 0])
    ax_zoom_resid = fig.add_subplot(gs[3, 0])
    ax_zoom_hist = fig.add_subplot(gs[3, 1], sharey=ax_zoom_resid)

    ax_full.errorbar(time, mag, yerr=mag_err, fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5, label="data")
    ax_full.plot(t_model, mag_model, color="crimson", label="PSPL fit")
    ax_full.invert_yaxis()
    ax_full.set_ylabel("I magnitude")
    ax_full.set_title(f"{dataset_label} (PSPL fit, full baseline)")
    ax_full.legend()
    ax_full.sharex(ax_full_resid)

    _plot_residual_panel(ax_full_resid, time, residuals)
    ax_full_resid.set_xlabel("HJD - 2450000")
    _plot_residual_hist(ax_full_hist, residuals)

    ax_zoom.errorbar(time, mag, yerr=mag_err, fmt="+", ms=3, elinewidth=0.5, capsize=2, markeredgewidth=0.5, capthick=0.5, label="data")
    ax_zoom.plot(t_model_zoom, mag_model_zoom, color="crimson", label="PSPL fit")
    ax_zoom.set_xlim(zoom_start, zoom_end)
    ax_zoom.set_ylim(mag[in_zoom].max() + 0.05, mag[in_zoom].min() - 0.05)
    ax_zoom.set_ylabel("I magnitude")
    ax_zoom.set_title("zoomed on peak (auto-detected)")
    ax_zoom.sharex(ax_zoom_resid)

    _plot_residual_panel(ax_zoom_resid, time[in_zoom], residuals[in_zoom])
    ax_zoom_resid.set_xlim(zoom_start, zoom_end)
    ax_zoom_resid.set_xlabel("HJD - 2450000")
    _plot_residual_hist(ax_zoom_hist, residuals[in_zoom])

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
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
    fig.savefig(out_path, dpi=600)
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
    coords = SkyCoord("18h04m45.70s -26d59m15.5s")
    time, mag, mag_err = np.loadtxt("data/OGLE-2005-BLG-086.dat", unpack=True)

    if args.stage == "raw":
        plot_raw(time, mag, mag_err, "OGLE-2005-BLG-086 (raw)", f"raw_lc/{SHORT_NAME}.png")
    else:
        t0_par = get_t0_par(time, mag, mag_err, f"data/processed/{SHORT_NAME}_plain_fit_summary.dat",
                             u0_guess=0.5, tE_guess=50.0)
        delta_sN, delta_sE = sun_earth_projection(time, coords, t0_par)

        fit_result = fit_parallax_pspl_mcmc(time, mag, mag_err, delta_sN, delta_sE, run_mcmc=(args.stage == "mcmc"))

        chi2 = np.sum(((mag - magnitude(time, *fit_result.best_fit, delta_sN, delta_sE)) / mag_err) ** 2)
        print(f"chi2/dof = {chi2 / (len(time) - len(fit_result.best_fit)):.3f}")
        for label, value in zip(fit_result.labels, fit_result.best_fit):
            print(f"{label} = {value:.5f}")
        plot_fit_lc(time, mag, mag_err, fit_result.best_fit, coords, t0_par, "OGLE-2005-BLG-086", f"fit_lc/{SHORT_NAME}.png")

        if args.stage == "mcmc":
            print(f"{fit_result.samples.shape[0]} posterior samples after burn-in/thinning")
            u0_s, tE_s, fs_s, fb_s = (fit_result.column(name) for name in ("u0", "tE", "f_source", "f_blend"))
            derived = {
                "A_max": magnification(u0_s),
                "t_eff": u0_s * tE_s,
                "blend_fraction": fb_s / (fs_s + fb_s),
                "m_source": ZERO_POINT_MAG - 2.5 * np.log10(fs_s),
                "M_lens": estimate_mass(tE_s),
            }
            results = {label: fit_result.column(label) for label in fit_result.labels} | derived
            save_summary(results, f"data/processed/{SHORT_NAME}_fit_summary.dat")
            plot_histograms(results, f"hist_plots/{SHORT_NAME}.png")

            # fit_result.best_fit is curve_fit's 7-param point estimate --
            # scale/dof have no curve_fit analog (least-squares has no
            # robust-likelihood concept), so pad with None to skip their
            # truth line in the plot.
            truths = list(fit_result.best_fit) + [None, None]
            save_corner(fit_result.samples, fit_result.labels, truths, f"corner_plots/{SHORT_NAME}.png")
