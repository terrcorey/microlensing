"""One-time cross-check: fit MulensModel's own binary-lens model to the same
OGLE+MOA magnification data fit_2l1s.py uses, as an independent check on
whether our hand-rolled optimizer/grid -- not the underlying physics -- is
what's failing to find the caustic-crossing solution.

Not a pipeline dependency; run manually in a scratch venv (see
cross_check_mulensmodel.py for the `pip install MulensModel` setup).
MulensModel's own magnification is ~1000x faster than our root-finding one
(no per-point Python-level root solve), so a plain serial multi-start is
fine here -- no need for fit_2l1s.py's parallelization/two-stage tolerance.
Lives in scratch/ with the project's other one-time/dev scripts, run
manually as `python3 scratch/cross_check_mulensmodel_fit.py` from the
project root.
"""
import sys
from itertools import product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import MulensModel as mm
import numpy as np
from scipy.optimize import minimize

from zoom_utils import plot_fit_panels

SHORT_NAME = "O-03-BLG235"
ogle_time, ogle_A, ogle_A_err = np.loadtxt(f"data/processed/{SHORT_NAME}_OGLE_magnification.dat", unpack=True)
moa_time, moa_A, moa_A_err = np.loadtxt(f"data/processed/{SHORT_NAME}_MOA_magnification.dat", unpack=True)

time = np.concatenate([ogle_time, moa_time])
A_obs = np.concatenate([ogle_A, moa_A])
A_err = np.concatenate([ogle_A_err, moa_A_err])

# t0/u0/tE seeded from the existing PSPL joint fit's posterior median.
param, val = np.loadtxt(f"data/processed/{SHORT_NAME}_fit_summary.dat", unpack=True, usecols=(0, 2), dtype=str)
guess = dict(zip(param, val.astype(float)))

# Same grid as fit_2l1s.py's multi-start, for an apples-to-apples comparison.
s_list = [0.7, 1.0, 1.5]
q_list = [0.01, 0.05, 0.1]
alpha_list = [0.0, 90.0, 180.0, 270.0]  # degrees -- MulensModel's native alpha unit
seeds = [[guess["t0"], guess["u0"], guess["tE"], alpha, s, q] for s, q, alpha in product(s_list, q_list, alpha_list)]


def mm_magnification(t, t0, u0, tE, alpha, s, q):
    params = mm.ModelParameters({"t_0": t0, "u_0": u0, "t_E": tE, "alpha": alpha, "s": s, "q": q})
    return mm.Model(params).get_magnification(t)


def chi2(theta):
    """Chi2 of MulensModel's binary magnification against the combined OGLE+MOA data."""
    t0, u0, tE, alpha, s, q = theta
    if tE <= 0 or s <= 0 or q <= 0:
        return np.inf  # unphysical
    try:
        A_model = mm_magnification(time, t0, u0, tE, alpha, s, q)
    except ValueError:
        return np.inf  # MulensModel rejects some degenerate (s, q) combos
    return np.sum(((A_obs - A_model) / A_err) ** 2)


def run_fit():
    """Serial multi-start search, then a tight-tolerance refit of the best seed."""
    print(f"[mm] starting multi-start search: {len(seeds)} seeds")
    results = []
    for seed in seeds:
        result = minimize(chi2, x0=seed, method="Nelder-Mead",
                           options={"xatol": 1e-4, "fatol": 1e-4, "maxiter": 10000})
        t0, u0, tE, alpha, s, q = seed
        print(f"[mm] seed s={s:.2f} q={q:.3f} alpha={alpha:.0f} -> chi2={result.fun:.2f}")
        results.append(result)
    best = min(results, key=lambda r: r.fun)
    print(f"[mm] multi-start done, best chi2={best.fun:.2f}, refining...")

    result = minimize(chi2, x0=best.x, method="Nelder-Mead",
                       options={"xatol": 1e-6, "fatol": 1e-6, "maxiter": 20000})
    print(f"[mm] final chi2={result.fun:.2f}")
    for label, value in zip(["t0", "u0", "tE", "alpha", "s", "q"], result.x):
        print(f"[mm] {label} = {value:.5f}")

    plot_fit(result.x)
    return result


def plot_fit(theta):
    """Overlay the MulensModel fit on the OGLE+MOA data, full baseline + zoomed on the peak."""
    t0, u0, tE, alpha, s, q = theta
    plot_fit_panels(
        (ogle_time, ogle_A, ogle_A_err), (moa_time, moa_A, moa_A_err),
        model_fn=lambda t: mm_magnification(t, t0, u0, tE, alpha, s, q),
        fit_label="MulensModel fit", out_path=f"scratch/{SHORT_NAME}_mulensmodel.png",
    )


if __name__ == "__main__":
    run_fit()
