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
from astropy.coordinates import SkyCoord
from scipy.optimize import minimize

from lc_models import ZERO_POINT_MAG, mag_to_flux, magnification, sun_earth_projection, trajectory

SHORT_NAME = "O-03-BLG235"  # see dataset_names.txt
OGLE_PATH = "data/OGLE-2003-BLG-235_OGLE.tbl.txt"
MOA_PATH = "data/OGLE-2003-BLG-235_MOA.tbl.txt"
LABELS = ["t0", "u0", "tE", "fs_ogle", "fb_ogle", "fs_moa", "piE_N", "piE_E"]


def _parse_header_coord(path):
    """Extract (RA, Dec) sexagesimal strings from a raw table's own \\RA/\\DEC
    header lines (NASA Exoplanet Archive format -- the same \\-prefixed
    metadata load_raw() already skips via comments=("\\", "|")). Stops at the
    first non-header line rather than scanning the whole file."""
    ra = dec = None
    with open(path) as f:
        for line in f:
            if not line.startswith("\\"):
                break
            if line.startswith("\\RA "):
                ra = line.split("=", 1)[1].strip().strip('"')
            elif line.startswith("\\DEC "):
                dec = line.split("=", 1)[1].strip().strip('"')
    return ra, dec


def load_coords():
    """Target coordinates, parsed from OGLE_PATH/MOA_PATH's own headers rather
    than hardcoded -- catches a transcription error against either file, or a
    disagreement between the two (this is how a real bug was caught: a
    previously hardcoded RA of 18h05m16.35s was 4 arcmin of RA off from what
    both raw files actually say, 18h01m16.35s)."""
    ogle_coord = _parse_header_coord(OGLE_PATH)
    moa_coord = _parse_header_coord(MOA_PATH)
    if ogle_coord != moa_coord:
        raise ValueError(f"OGLE and MOA header coordinates disagree: OGLE={ogle_coord} vs MOA={moa_coord}")
    ra, dec = ogle_coord
    return SkyCoord(f"{ra} {dec}")


COORDS = load_coords()


def load_raw():
    """Raw OGLE (mag) + MOA (differential flux) tables, time-shifted to HJD-2450000."""
    ogle_time, ogle_mag, ogle_err = np.loadtxt(
        OGLE_PATH, comments=("\\", "|"), unpack=True
    )
    moa_time, moa_flux, moa_err = np.loadtxt(
        MOA_PATH, comments=("\\", "|"), unpack=True
    )
    ogle_time -= 2450000.0
    moa_time -= 2450000.0
    return ogle_time, ogle_mag, ogle_err, moa_time, moa_flux, moa_err


def fit_joint_pspl(ogle_time, ogle_mag, ogle_err, moa_time, moa_flux, moa_err,
                    delta_sN_ogle, delta_sE_ogle, delta_sN_moa, delta_sE_moa):
    """Quick point-estimate joint fit, now including annual parallax
    (piE_N, piE_E) -- see CLAUDE.md's "Annual parallax + robust likelihood"
    roadmap, step 6. delta_s*=0 collapses trajectory()'s parallax terms to
    zero regardless of piE, so this same function also serves as the plain
    (pre-parallax) fit used to bootstrap t0_par -- see __main__. Returns
    best_fit array in LABELS order."""

    def total_chi2(theta):
        t0, u0, tE, fs_ogle, fb_ogle, fs_moa, piE_N, piE_E = theta
        if tE <= 0 or fs_moa <= 0:
            return np.inf

        ogle_A = magnification(trajectory(ogle_time, t0, u0, tE, piE_N, piE_E, delta_sN_ogle, delta_sE_ogle))
        ogle_model_flux = fs_ogle * ogle_A + fb_ogle
        if np.any(ogle_model_flux <= 0):
            return np.inf
        ogle_model_mag = ZERO_POINT_MAG - 2.5 * np.log10(ogle_model_flux)
        ogle_chi2 = np.sum(((ogle_mag - ogle_model_mag) / ogle_err) ** 2)

        moa_A = magnification(trajectory(moa_time, t0, u0, tE, piE_N, piE_E, delta_sN_moa, delta_sE_moa))
        moa_model_relflux = fs_moa * (moa_A - 1.0)
        moa_chi2 = np.sum(((moa_flux - moa_model_relflux) / moa_err) ** 2)

        return ogle_chi2 + moa_chi2

    t0_guess = ogle_time[np.argmin(ogle_mag)]
    fs_ogle_guess = 10 ** (-0.4 * (np.median(ogle_mag) - ZERO_POINT_MAG))
    fs_moa_guess = moa_flux.max() / 4.0
    p0_guess = [t0_guess, 0.2, 30.0, fs_ogle_guess, 0.0, fs_moa_guess, 0.0, 0.0]

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

    zeros_ogle, zeros_moa = np.zeros_like(ogle_time), np.zeros_like(moa_time)
    plain_fit, _ = fit_joint_pspl(ogle_time, ogle_mag, ogle_err, moa_time, moa_flux, moa_err,
                                   zeros_ogle, zeros_ogle, zeros_moa, zeros_moa)
    t0_par = plain_fit[0]
    print(f"t0_par (parallax reference epoch, from plain pre-parallax fit) = {t0_par:.5f}")

    delta_sN_ogle, delta_sE_ogle = sun_earth_projection(ogle_time, COORDS, t0_par)
    delta_sN_moa, delta_sE_moa = sun_earth_projection(moa_time, COORDS, t0_par)

    best_fit, _ = fit_joint_pspl(ogle_time, ogle_mag, ogle_err, moa_time, moa_flux, moa_err,
                                  delta_sN_ogle, delta_sE_ogle, delta_sN_moa, delta_sE_moa)
    print("calibration fit (used only to derive fs_ogle/fb_ogle/fs_moa below):")
    for label, value in zip(LABELS, best_fit):
        print(f"  {label} = {value:.5f}")
    _, _, _, fs_ogle, fb_ogle, fs_moa, piE_N, piE_E = best_fit

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
