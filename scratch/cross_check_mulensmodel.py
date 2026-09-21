"""One-time offline cross-check of lc_models.binary_magnification against
MulensModel (an independent, published binary-lens implementation), as a
trusted oracle -- not a pipeline dependency. Run manually in a scratch venv
with `pip install MulensModel` (not a project dependency).

Both codes place the origin at the binary's center of mass with the higher
mass (q<1) on the negative real axis -- confirmed by reading MulensModel's
source. But its Trajectory docstring warns its *alpha* is shifted 180 deg
from the standard convention it otherwise follows (Skowron et al. 2011,
Appendix A); confirmed empirically that at alpha=0, MulensModel's
(traj.x, traj.y) = (-t, -u0). So t=-zeta.real, u_0=-zeta.imag reproduces a
target zeta exactly, letting us compare magnifications at matched, arbitrary
zeta directly without going through any fitting/trajectory code of our own.

Result (last run): across 6 (s, q) topologies (close/resonant/wide x two
mass ratios) x far-field and near-caustic points, worst-case relative
difference was 1.4e-9 -- everything else at floating-point precision.

Lives in scratch/ with the project's other one-time/dev scripts, run
manually as `python3 scratch/cross_check_mulensmodel.py` from the project
root.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import MulensModel as mm

import lc_models as lm


def mm_magnification(zeta_values, s, q):
    """MulensModel's magnification at the same complex zeta values.

    MulensModel's docs warn its alpha is shifted 180 deg from the standard
    convention; confirmed empirically that at alpha=0, traj.x=-t and
    traj.y=-u0, so t=-zeta.real, u_0=-zeta.imag reproduces zeta exactly.
    """
    out = []
    for zeta in zeta_values:
        params = mm.ModelParameters(
            {"t_0": 0.0, "u_0": -zeta.imag, "t_E": 1.0, "alpha": 0.0, "s": s, "q": q}
        )
        model = mm.Model(params)
        out.append(model.get_magnification(np.array([-zeta.real]))[0])
    return np.array(out)


def report(label, s, q, zeta_values):
    mine = np.array([lm.binary_magnification(z, s, q) for z in zeta_values])
    theirs = mm_magnification(zeta_values, s, q)
    rel_diff = np.abs(mine - theirs) / theirs
    print(f"\n--- {label} (s={s}, q={q}) ---")
    for z, a, b, d in zip(zeta_values, mine, theirs, rel_diff):
        flag = "  <-- MISMATCH" if d > 1e-3 else ""
        print(f"zeta={z:+.4f}  mine={a:12.6f}  MulensModel={b:12.6f}  rel_diff={d:.2e}{flag}")
    print(f"max rel_diff: {rel_diff.max():.2e}")
    return rel_diff.max()


rng = np.random.default_rng(42)
worst = 0.0

topologies = [
    ("close", 0.6, 0.1),
    ("close, near-equal mass", 0.6, 0.5),
    ("resonant", 1.0, 0.1),
    ("resonant (O-03-BLG235-like)", 1.0, 0.05),
    ("wide", 1.8, 0.1),
    ("wide, near-equal mass", 1.8, 0.5),
]

for label, s, q in topologies:
    # far-from-caustic random points
    far_points = np.array(
        [complex(rng.uniform(-2.5, 2.5), rng.uniform(-2.5, 2.5)) for _ in range(5)]
    )
    worst = max(worst, report(f"{label}, random far points", s, q, far_points))

    # near-caustic points: perturb actual caustic points slightly
    caustics = mm.CausticsBinary(q=q, s=s)
    cx, cy = caustics.get_caustics(n_points=2000)
    idx = rng.choice(len(cx), size=5, replace=False)
    near_points = np.array(
        [complex(cx[i], cy[i]) + complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) * 1e-3
         for i in idx]
    )
    worst = max(worst, report(f"{label}, near-caustic points", s, q, near_points))

print(f"\n=== overall worst-case relative difference across all tests: {worst:.2e} ===")
