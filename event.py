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


def load_instrument(cfg: dict, coords: SkyCoord) -> Instrument:
    """One [[instruments]] entry -> Instrument: loadtxt, time shift by time_fmt,
    K applied to the raw error (before mag_to_flux for kind="mag")."""
    if cfg["kind"] not in ("mag", "dia"):
        raise ValueError(f"unknown kind {cfg['kind']!r}")
    time, val, err = np.loadtxt(cfg["path"], comments=("\\", "|"), unpack=True)
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
    flux, flux_err = mag_to_flux(val, cfg["K"] * err) if cfg["kind"] == "mag" else (val, cfg["K"] * err)
    return Instrument(cfg["name"], cfg["kind"], cfg["band"], cfg["ld"], t, flux, flux_err)



def load_event(path: str) -> Event:
    """Parse the TOML at `path` (tomllib, opened "rb") and load every instrument."""
    with open(path, "rb") as f:
        cfg = tomllib.load(f)
    coords = SkyCoord(f'{cfg["ra"]} {cfg["dec"]}')
    return Event(cfg["short_name"], coords, cfg["grid"],
                 [load_instrument(inst, coords) for inst in cfg["instruments"]])


def flux_residuals(event: Event, A: list[np.ndarray]) -> np.ndarray:
    """Any model's magnification per instrument (same order as event.instruments)
    -> profile each instrument's fs/fb by weighted least squares -> concatenated
    standardized residuals. Generalises scratch/fit_2l1s.py's flux_residuals()."""
    residuals = []
    for inst, A_i in zip(event.instruments, A, strict=True):
        if not np.all(np.isfinite(A_i)):
            return np.full(sum(i.time.size for i in event.instruments), np.inf)
        X = np.column_stack([A_i, np.ones_like(A_i)]) if inst.kind == "mag" else (A_i - 1.0)[:, None]
        Xw, Fw = X / inst.flux_err[:, None], inst.flux / inst.flux_err  # rows / sigma -> weighted LSQ
        coeffs = np.linalg.lstsq(Xw, Fw)[0]  # (fs, fb) or (fs,)
        if coeffs[0] <= 0:  # negative source flux: unphysical trial
            return np.full(sum(i.time.size for i in event.instruments), np.inf)
        residuals.append(Fw - Xw @ coeffs)
    return np.concatenate(residuals)


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
