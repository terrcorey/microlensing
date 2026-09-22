"""One-time offline cross-check of lc_models' annual-parallax PSPL trajectory
against MulensModel (an independent, published implementation), as a trusted
oracle -- not a pipeline dependency. Run manually in the project venv with
`pip install MulensModel` (scratch-only, not in requirements.txt, same
treatment as the existing cross_check_mulensmodel.py for 2L1S).

Two independent unknowns are tangled together here, so this tests all four
combinations rather than guessing which one matters:

1. delta_tau/delta_beta sign convention (flagged unverified in
   lc_models.trajectory()'s docstring). MulensModel's own
   Trajectory._project_delta() computes:
       delta_tau  =  delta_N * pi_E_N + delta_E * pi_E_E
       delta_beta = -delta_N * pi_E_E + delta_E * pi_E_N
   lc_models.trajectory() currently computes the same delta_tau but the
   mirror-sign delta_beta = piE_E * delta_sN - piE_N * delta_sE.

2. delta_s(t) definition: lc_models.sun_earth_projection() projects Earth's
   position *relative to the Sun* (heliocentric), while MulensModel's
   _get_delta_annual() projects Earth's bare barycentric position linearly
   detrended at t_0_par, with no Sun subtraction -- physically similar
   (both isolate the curvature of Earth's orbit) but not identical, since
   the Sun's own barycentric wobble (Jupiter etc.) is not exactly zero.

A first pass comparing only the two sign candidates (both using our
heliocentric delta_s) gave a few-percent residual for BOTH signs -- far
short of the 2L1S cross-check's ~1e-9 agreement, and too close between the
two candidates to call the sign convention from that alone. This version
adds the bare-barycentric delta_s variant (matching MulensModel's own
definition exactly) crossed with both signs, to isolate which of the two
unknowns is actually responsible for the residual.

Lives in scratch/ like the project's other one-time cross-checks; run
manually as `python3 scratch/cross_check_mulensmodel_parallax.py` from the
project root.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from astropy.coordinates import SkyCoord, get_body_barycentric_posvel
from astropy.time import Time
import astropy.units as u
import MulensModel as mm

import lc_models as lm

# O-05-BLG086's own coordinates (this session's parallax proving ground),
# OGLE-III EWS alert page, field BLG234.6 -- see CLAUDE.md.
RA_STR = "18h04m45.70s"
DEC_STR = "-26d59m15.5s"
COORDS = SkyCoord(f"{RA_STR} {DEC_STR}")

# Arbitrary but physically reasonable test parameters, project convention
# (HJD - 2450000). t0_par fixed at t0 here, matching the Gould (2004)
# convention this project follows (t0_par = the dataset's own preliminary
# non-parallax best-fit t0, not a fitted value itself).
T0 = 3600.0
U0 = 0.1
TE = 60.0
PIE_N = 0.15
PIE_E = -0.08
T0_PAR = T0

JD_OFFSET = 2450000.0  # project convention subtracts this; MulensModel wants full JD

times = np.linspace(T0 - 200, T0 + 200, 41)


def bare_earth_projection(time, coords, t0_par):
    """Same linear-detrending-at-t0_par logic as lc_models.sun_earth_projection(),
    but on Earth's bare barycentric position (no Sun subtraction) -- matching
    MulensModel's own Trajectory._get_delta_annual() definition exactly, to
    isolate whether the heliocentric-vs-barycentric choice (rather than the
    tau/beta sign) is the source of the residual seen with the heliocentric
    version.
    """
    time_hjd = Time(time + JD_OFFSET, format="jd")
    t0 = Time(t0_par + JD_OFFSET, format="jd")

    earth_pos, _ = get_body_barycentric_posvel("earth", time_hjd)
    earth_t0pos, earth_t0vel = get_body_barycentric_posvel("earth", t0)

    s_N, s_E = lm._project(earth_pos, coords)
    s0_N, s0_E = lm._project(earth_t0pos, coords)
    v0_N, v0_E = lm._project(earth_t0vel, coords)

    s_N, s_E = s_N.to_value(u.au), s_E.to_value(u.au)
    s0_N, s0_E = s0_N.to_value(u.au), s0_E.to_value(u.au)
    v0_N, v0_E = v0_N.to_value(u.au / u.day), v0_E.to_value(u.au / u.day)

    dt_days = (time_hjd - t0).to_value(u.day)
    delta_sN = s_N - s0_N - dt_days * v0_N
    delta_sE = s_E - s0_E - dt_days * v0_E
    return delta_sN, delta_sE


delta_sN_helio, delta_sE_helio = lm.sun_earth_projection(times, COORDS, T0_PAR)
delta_sN_bary, delta_sE_bary = bare_earth_projection(times, COORDS, T0_PAR)


def our_magnification(delta_sN, delta_sE, delta_s_sign, beta_sign):
    """lc_models' trajectory formula, with the delta_s definition, delta_s's
    own overall sign (position(t)-position(t0) vs MulensModel's
    position(t0)-position(t)), and delta_beta's combination sign all left as
    free switches so every candidate combination can be tested without
    editing lc_models.py itself."""
    delta_sN = delta_s_sign * delta_sN
    delta_sE = delta_s_sign * delta_sE
    delta_tau = PIE_N * delta_sN + PIE_E * delta_sE
    delta_beta = beta_sign * (PIE_E * delta_sN - PIE_N * delta_sE)
    tau = (times - T0) / TE + delta_tau
    beta = U0 + delta_beta
    u_val = np.sqrt(tau**2 + beta**2)
    return lm.magnification(u_val)


def mulensmodel_magnification():
    params = mm.ModelParameters({
        "t_0": T0 + JD_OFFSET,
        "u_0": U0,
        "t_E": TE,
        "pi_E_N": PIE_N,
        "pi_E_E": PIE_E,
        "t_0_par": T0_PAR + JD_OFFSET,
    })
    model = mm.Model(params, coords=f"{RA_STR} {DEC_STR}")
    return model.get_magnification(times + JD_OFFSET)


theirs = mulensmodel_magnification()

definitions = {"helio": (delta_sN_helio, delta_sE_helio), "bary": (delta_sN_bary, delta_sE_bary)}
candidates = {}
for def_label, (dN, dE) in definitions.items():
    for delta_s_sign, sign_label in ((+1, "s(t)-s(t0)"), (-1, "s(t0)-s(t)")):
        for beta_sign, beta_label in ((+1, "as-written"), (-1, "flipped")):
            label = f"{def_label}, {sign_label}, beta {beta_label}"
            candidates[label] = our_magnification(dN, dE, delta_s_sign, beta_sign)

worst = {}
for label, curve in candidates.items():
    rel_diff = np.abs(curve - theirs) / theirs
    worst[label] = rel_diff.max()

for label, max_rel in sorted(worst.items(), key=lambda kv: kv[1]):
    print(f"max rel_diff, {label:35s}: {max_rel:.2e}")

best_label = min(worst, key=worst.get)
print(f"\n=> best match: {best_label} (max rel_diff {worst[best_label]:.2e})")
if worst[best_label] > 1e-3:
    print("   Still far short of the 2L1S cross-check's ~1e-9 -- do not treat this as locked "
          "down yet; something beyond sign/delta_s-definition may still be off.")
else:
    print("   Close enough to treat this combination as the correct convention.")
