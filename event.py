"""One event's config + data, and the model-agnostic half of every chi2.

Library only -- search.py is the CLI. Nothing is loaded at import time:
search.py's __main__ calls load_event() once and passes the Event down.

Every instrument ends up as (time, flux, flux_err), whatever its raw format;
after loading, `kind` only decides the flux-calibration design matrix:
    mag -> F = fs * A + fb    (columns [A, 1])
    dia -> F = fs * (A - 1)   (column  [A - 1]; difference imaging has no blend)
"""
from typing import NamedTuple

import numpy as np
import tomllib
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time


from lc_models import mag_to_flux, sun_earth_projection


class Instrument(NamedTuple):
    """One telescope/band: its [[instruments]] config entry plus its loaded data."""
    name: str
    kind: str        # "mag" | "dia"
    band: str
    ld: float        # linear limb-darkening coefficient, VBBL's a1 (u convention, see lc_models)
    time: np.ndarray       # HJD - 2450000
    flux: np.ndarray
    flux_err: np.ndarray   # already multiplied by the config's K
    # Earth's projected offset (lc_models.sun_earth_projection) -- zero until with_t0_par();
    # only ever multiplied by piE, so zero is exact for parallax-free models
    dsN: np.ndarray | float = 0.0
    dsE: np.ndarray | float = 0.0


class Event(NamedTuple):
    short_name: str        # output filenames, see dataset_names.txt
    coords: SkyCoord
    grid: dict             # the config's [grid] table, as parsed
    instruments: list[Instrument]
    # (band, mean, sigma) per band whose LD is free (an instrument's optional `ld_sigma`), sorted by band:
    # a fit's LD coefficients follow this order; () = every LD fixed at the config's (session 24)
    ld_prior: tuple = ()
    tE_max: float = np.inf  # 2L1S tE ceiling, search.TE_FACTOR x the FSPL tE (session 25)


def load_instrument(cfg: dict, coords: SkyCoord) -> Instrument:
    """One [[instruments]] entry -> Instrument: loadtxt, time shift by time_fmt,
    K applied to the raw error (before mag_to_flux for kind="mag")."""
    if cfg["kind"] not in ("mag", "dia"):
        raise ValueError(f"unknown kind {cfg['kind']!r}")
    time, val, err = np.loadtxt(cfg["path"], comments=("\\", "|"), unpack=True)
    good = np.isfinite(time) & np.isfinite(val) & np.isfinite(err) & (err > 0)  # err = 0 / NaN -> NaN chi2
    if not good.all():
        print(f"[load] {cfg['name']}: dropped {np.sum(~good)} non-finite / zero-error rows")
    order = np.argsort(time[good])  # fit_fspl's median-filtered t0 guess assumes time order
    time, val, err = time[good][order], val[good][order], err[good][order]
    if cfg["time_fmt"] == "HJD":
        t = time - 2450000
    elif cfg["time_fmt"] == "HJD-2450000":
        t = time
    elif cfg["time_fmt"] == "Geocentric JD":
        jd = Time(time, format="jd", scale="utc",
                 location=EarthLocation.from_geocentric(0, 0, 0, unit="m"))
        hjd = np.asarray((jd + jd.light_travel_time(coords, kind="heliocentric")).jd)
        t = hjd - 2450000
    else:
        raise ValueError(f"unknown time_fmt {cfg['time_fmt']!r}")
    if not 0 < t.min() <= t.max() < 20000:  # HJD - 2450000 spans ~1995-2050
        raise ValueError(f"{cfg['name']}: times {t.min():.1f}..{t.max():.1f} after the shift -- check time_fmt")
    flux, flux_err = mag_to_flux(val, cfg["K"] * err) if cfg["kind"] == "mag" else (val, cfg["K"] * err)
    return Instrument(cfg["name"], cfg["kind"], cfg["band"], cfg["ld"], t, flux, flux_err)



def load_event(path: str) -> Event:
    """Parse the TOML at `path` (tomllib, opened "rb") and load every instrument."""
    with open(path, "rb") as f:
        cfg = tomllib.load(f)
    coords = SkyCoord(f'{cfg["ra"]} {cfg["dec"]}')
    free = sorted({inst["band"] for inst in cfg["instruments"] if "ld_sigma" in inst})
    prior = {}
    for inst in cfg["instruments"]:  # one coefficient per band: every instrument in it must agree on the prior
        if inst["band"] in free and prior.setdefault(inst["band"], (inst["ld"], inst.get("ld_sigma"))) != (
                inst["ld"], inst.get("ld_sigma")):
            raise ValueError(f"band {inst['band']!r}: free LD needs the same ld and ld_sigma on every instrument in it")
    return Event(cfg["short_name"], coords, cfg["grid"],
                 [load_instrument(inst, coords) for inst in cfg["instruments"]], tuple((b, *prior[b]) for b in free))


def with_ld(event: Event, ld) -> Event:
    """Event with the free bands' LD set to `ld` (one per event.ld_prior entry); () leaves the config's."""
    if not len(ld):
        return event
    band = {b: float(a) for (b, _, _), a in zip(event.ld_prior, ld, strict=True)}
    return event._replace(instruments=[i._replace(ld=band.get(i.band, i.ld)) for i in event.instruments])


def ld_penalty(event: Event, ld) -> float:
    """Gaussian LD prior as a chi2 term: sum ((a - mean) / sigma)^2 over the free bands; 0 for ld = ()."""
    return float(sum(((a - m) / sd) ** 2 for (_, m, sd), a in zip(event.ld_prior, ld)))


def flux_residuals(event: Event, A: list[np.ndarray]) -> np.ndarray:
    """Any model's magnification per instrument (same order as event.instruments)
    -> profile each instrument's fs/fb by weighted least squares -> concatenated
    standardized residuals. Generalises scratch/fit_2l1s.py's flux_residuals()."""
    residuals = []
    for inst, A_i in zip(event.instruments, A, strict=True):
        if not np.all(np.isfinite(A_i)):
            return np.full(sum(i.time.size for i in event.instruments), np.inf)
        coeffs, Xw, Fw = profile_flux(inst, A_i)
        if coeffs[0] <= 0:  # negative source flux: unphysical trial
            return np.full(sum(i.time.size for i in event.instruments), np.inf)
        residuals.append(Fw - Xw @ coeffs)
    return np.concatenate(residuals)


def profile_flux(inst: Instrument, A_i: np.ndarray):
    """Weighted least-squares (fs, fb) ("mag": flux = fs*A + fb) or (fs,) ("dia": flux =
    fs*(A - 1)) at this magnification, plus the weighted design matrix and data."""
    X = np.column_stack([A_i, np.ones_like(A_i)]) if inst.kind == "mag" else (A_i - 1.0)[:, None]
    Xw, Fw = X / inst.flux_err[:, None], inst.flux / inst.flux_err  # rows / sigma -> weighted LSQ
    return np.linalg.lstsq(Xw, Fw)[0], Xw, Fw


def to_magnification(inst: Instrument, A_i: np.ndarray):
    """Invert an instrument's flux onto magnification via its fs/fb profiled at model A_i:
    (A, A_err), the common scale every instrument can share in one plot."""
    fs, *fb = profile_flux(inst, A_i)[0]
    A = (inst.flux - fb[0]) / fs if inst.kind == "mag" else 1 + inst.flux / fs
    return A, inst.flux_err / fs


def chi2(event: Event, A: list[np.ndarray]) -> float:
    return float(np.sum(flux_residuals(event, A) ** 2))


def with_t0_par(event: Event, t0_par: float) -> Event:
    """Event with each instrument's parallax offsets precomputed at t0_par -- one astropy
    ephemeris query per instrument, never per chi2 call."""
    return event._replace(instruments=[
        inst._replace(**dict(zip(("dsN", "dsE"), sun_earth_projection(inst.time, event.coords, t0_par))))
        for inst in event.instruments])


def error_scale(event: Event, A: list[np.ndarray]) -> list[float]:
    """Per-instrument k making chi2_i / (N_i - n_flux_i) = 1 at this model (no error floor,
    see CLAUDE.md's error-bar rescaling note). ponytail: the shared model parameters aren't
    split across instruments' dof -- negligible while N_i >> n_params."""
    sizes = [inst.time.size for inst in event.instruments]
    chunks = np.split(flux_residuals(event, A), np.cumsum(sizes)[:-1])
    return [float(np.sqrt(np.sum(r ** 2) / (n - (2 if inst.kind == "mag" else 1))))
            for r, n, inst in zip(chunks, sizes, event.instruments)]


def rescale(event: Event, k: list[float]) -> Event:
    """Event with each instrument's flux_err multiplied by its k (on top of the config's K)."""
    return event._replace(instruments=[inst._replace(flux_err=inst.flux_err * ki)
                                       for inst, ki in zip(event.instruments, k, strict=True)])
