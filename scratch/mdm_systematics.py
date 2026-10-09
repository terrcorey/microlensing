"""MDM systematics check (session 24, one-off): is O-05-BLG169's MDM night-1 wave photometric rather than
astrophysical? MDM's file has no seeing/sky/airmass columns, so: fractional residuals at the session-22 fine-q 2L1S
best vs airmass (computed for Kitt Peak, HJD used as JD -- the <= 8 min difference is irrelevant for airmass) and
vs the reported per-point error (a seeing/sky proxy), wave stats before/after removing a linear trend in each, and
a Lomb-Scargle periodogram of the residuals. Writes scratch/ld_test/O-05-BLG169_mdm_systematics.png. Never on the
login node (runs the blind FSPL fit for the parallax offsets, ~1 min):

    srun --partition=small-short --cpus-per-task=2 --mem=4G --time=00:15:00 .venv/bin/python scratch/mdm_systematics.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import AltAz, EarthLocation
from astropy.time import Time
from astropy.timeseries import LombScargle

from event import load_event, profile_flux, rescale
from ld_test import CONFIG, OUT, read_summary, wave
from search import binary_A, fit_fspl

MDM_SITE = EarthLocation(lat=31.9499 * u.deg, lon=-111.6167 * u.deg, height=1938 * u.m)  # MDM, Kitt Peak


def detrended(r, x, err):
    """Standardized residuals left after a weighted linear fit of r against x."""
    c = np.polyfit(x, r, 1, w=1 / err)
    return (r - np.polyval(c, x)) / err


if __name__ == "__main__":
    _, best, k = read_summary()
    _, _, base = fit_fspl(load_event(CONFIG))
    event = rescale(base, k)
    mdm = [i.name for i in event.instruments].index("MDM")
    inst = event.instruments[mdm]
    A = binary_A(event, best)[mdm]
    fs, fb = profile_flux(inst, A)[0]
    model = fs * A + fb
    night = inst.time < 3492.4
    t, r, err = inst.time[night], (inst.flux / model - 1)[night], (inst.flux_err / model)[night]
    X = AltAz(obstime=Time(t + 2450000, format="jd"), location=MDM_SITE)
    airmass = np.asarray(event.coords.transform_to(X).secz)

    dt = np.diff(t) * 1440
    print(f"MDM night 1: {t.size} points, {t[0]:.4f}-{t[-1]:.4f}, cadence median {np.median(dt):.2f} min, "
          f"max gap {dt.max():.1f} min; airmass {airmass.min():.2f}-{airmass.max():.2f}")
    print(f"fractional residual rms {np.std(r) * 100:.2f}%, median reported error {np.median(err) * 100:.2f}%")
    print(f"corr(residual, airmass) {np.corrcoef(r, airmass)[0, 1]:+.2f}, corr(residual, error) {np.corrcoef(r, err)[0, 1]:+.2f}")
    print("wave (runs z, lag-1, chi2/pt):")
    print("  as is                {:+.2f}, {:+.2f}, {:.2f}".format(*wave(r / err)))
    print("  minus linear airmass {:+.2f}, {:+.2f}, {:.2f}".format(*wave(detrended(r, airmass, err))))
    print("  minus linear error   {:+.2f}, {:+.2f}, {:.2f}".format(*wave(detrended(r, err, err))))
    ls = LombScargle(t, r, err)
    freq, power = ls.autopower(minimum_frequency=1440 / 180, maximum_frequency=1440 / 10, samples_per_peak=20)
    peak = power.argmax()
    print(f"Lomb-Scargle peak: {1440 / freq[peak]:.1f} min, power {power[peak]:.2f}, "
          f"false-alarm prob {ls.false_alarm_probability(power[peak]):.2g} (white-noise assumption)")

    fig, ax = plt.subplots(2, 2, figsize=(11, 7))
    ax[0, 0].errorbar(t, r * 100, err * 100, fmt=".", ms=3, lw=0.5, color="k")
    tw = ax[0, 0].twinx()
    tw.plot(t, airmass, "C1-")
    tw.set_ylabel("airmass", color="C1")
    ax[0, 0].set(xlabel="HJD - 2450000", ylabel="residual (%)", title="residual and airmass")
    ax[0, 1].plot(airmass, r * 100, ".", ms=3, color="k")
    ax[0, 1].set(xlabel="airmass", ylabel="residual (%)", title=f"r = {np.corrcoef(r, airmass)[0, 1]:+.2f}")
    ax[1, 0].plot(err * 100, r * 100, ".", ms=3, color="k")
    ax[1, 0].set(xlabel="reported error (%)", ylabel="residual (%)", title=f"r = {np.corrcoef(r, err)[0, 1]:+.2f}")
    ax[1, 1].plot(1440 / freq, power, "k-", lw=0.8)
    ax[1, 1].set(xlabel="period (min)", ylabel="Lomb-Scargle power", title=f"peak {1440 / freq[peak]:.1f} min")
    for a in ax.flat:
        a.axhline(0, color="k", lw=0.4) if a is not ax[1, 1] else None
    fig.suptitle("O-05-BLG169 MDM night 1: residuals at the session-22 2L1S best")
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "O-05-BLG169_mdm_systematics.png", dpi=200)
    print(f"saved {OUT / 'O-05-BLG169_mdm_systematics.png'}")
