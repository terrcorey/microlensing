"""Reduce the raw OGLE + MOA files for OGLE-2003-BLG-235 / MOA-2003-BLG-53
into one calibrated product: magnification A(t) on a common scale.

Runs a joint PSPL fit ONCE to get the flux calibration -- shared trajectory
(t0, u0, tE), per-instrument flux scaling:
  OGLE (calibrated I mag, total flux) : model_flux = fs_ogle * A(t) + fb_ogle
  MOA  (DIA-style differential flux)  : model_flux = fs_moa * (A(t) - 1)
(MOA has no additive term because its differencing already cancels out any
constant blend) -- then converts both instruments' data to magnification and
writes them out. Downstream scripts (raw plot, fit overlay, MCMC) just load
these processed files and work with the plain 3-parameter PSPL model -- no
per-instrument flux bookkeeping needed once this has run.

Note: this freezes the flux calibration at this fit's point estimate, so
downstream (t0, u0, tE) uncertainties won't include calibration uncertainty.
Fine for exploration; re-run this if the raw data changes.
"""

from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from lc_models import ZERO_POINT_MAG, mag_to_flux, magnification, trajectory

SHORT_NAME = "O-03-BLG235"  # see dataset_names.txt
LABELS = ["t0", "u0", "tE", "fs_ogle", "fb_ogle", "fs_moa"]


def load_raw():
    """Raw OGLE (mag) + MOA (differential flux) tables, time-shifted to HJD-2450000."""
    ogle_time, ogle_mag, ogle_err = np.loadtxt(
        "data/OGLE-2003-BLG-235_OGLE.tbl.txt", comments=("\\", "|"), unpack=True
    )
    moa_time, moa_flux, moa_err = np.loadtxt(
        "data/OGLE-2003-BLG-235_MOA.tbl.txt", comments=("\\", "|"), unpack=True
    )
    ogle_time -= 2450000.0
    moa_time -= 2450000.0
    return ogle_time, ogle_mag, ogle_err, moa_time, moa_flux, moa_err


def fit_joint_pspl(ogle_time, ogle_mag, ogle_err, moa_time, moa_flux, moa_err):
    """Quick point-estimate joint fit. Returns best_fit array in LABELS order."""

    def total_chi2(theta):
        t0, u0, tE, fs_ogle, fb_ogle, fs_moa = theta
        if tE <= 0 or fs_moa <= 0:
            return np.inf

        ogle_model_flux = fs_ogle * magnification(trajectory(ogle_time, t0, u0, tE)) + fb_ogle
        if np.any(ogle_model_flux <= 0):
            return np.inf
        ogle_model_mag = ZERO_POINT_MAG - 2.5 * np.log10(ogle_model_flux)
        ogle_chi2 = np.sum(((ogle_mag - ogle_model_mag) / ogle_err) ** 2)

        moa_model_relflux = fs_moa * (magnification(trajectory(moa_time, t0, u0, tE)) - 1.0)
        moa_chi2 = np.sum(((moa_flux - moa_model_relflux) / moa_err) ** 2)

        return ogle_chi2 + moa_chi2

    t0_guess = ogle_time[np.argmin(ogle_mag)]
    fs_ogle_guess = 10 ** (-0.4 * (np.median(ogle_mag) - ZERO_POINT_MAG))
    fs_moa_guess = moa_flux.max() / 4.0
    p0_guess = [t0_guess, 0.2, 30.0, fs_ogle_guess, 0.0, fs_moa_guess]

    result = minimize(total_chi2, x0=p0_guess, method="Nelder-Mead",
                       options={"xatol": 1e-6, "fatol": 1e-6, "maxiter": 20000})
    best_fit = result.x
    best_fit[1], best_fit[2] = abs(best_fit[1]), abs(best_fit[2])  # u0, tE sign is arbitrary
    return best_fit, total_chi2


def ogle_to_magnification(mag, mag_err, fs_ogle, fb_ogle):
    """Invert OGLE's flux model to recover A(t): A = (F - fb_ogle) / fs_ogle."""
    flux, flux_err = mag_to_flux(mag, mag_err)
    return (flux - fb_ogle) / fs_ogle, flux_err / fs_ogle


def moa_to_magnification(flux, flux_err, fs_moa):
    """Invert MOA's flux model to recover A(t): A = 1 + F / fs_moa."""
    return 1.0 + flux / fs_moa, flux_err / fs_moa


if __name__ == "__main__":
    ogle_time, ogle_mag, ogle_err, moa_time, moa_flux, moa_err = load_raw()

    best_fit, _ = fit_joint_pspl(ogle_time, ogle_mag, ogle_err, moa_time, moa_flux, moa_err)
    print("calibration fit (used only to derive fs_ogle/fb_ogle/fs_moa below):")
    for label, value in zip(LABELS, best_fit):
        print(f"  {label} = {value:.5f}")
    _, _, _, fs_ogle, fb_ogle, fs_moa = best_fit

    ogle_A, ogle_A_err = ogle_to_magnification(ogle_mag, ogle_err, fs_ogle, fb_ogle)
    moa_A, moa_A_err = moa_to_magnification(moa_flux, moa_err, fs_moa)

    Path("data/processed").mkdir(parents=True, exist_ok=True)

    ogle_out = f"data/processed/{SHORT_NAME}_OGLE_magnification.dat"
    np.savetxt(ogle_out, np.column_stack([ogle_time, ogle_A, ogle_A_err]),
               header="HJD-2450000  magnification  magnification_err")
    print(f"saved {ogle_out}")

    moa_out = f"data/processed/{SHORT_NAME}_MOA_magnification.dat"
    np.savetxt(moa_out, np.column_stack([moa_time, moa_A, moa_A_err]),
               header="HJD-2450000  magnification  magnification_err")
    print(f"saved {moa_out}")
