"""One-time offline cross-check of lc_models.binary_magnification_fs against
MulensModel's finite-source binary-lens magnification (VBBL, contour
integration), as a trusted oracle -- not a pipeline dependency. Needs a separate
`pip install MulensModel`. Same zeta mapping as cross_check_mulensmodel.py
(t=-zeta.real, u_0=-zeta.imag at alpha=0).

Steps the source along the caustic's normal through the best Cassan chain
sample's exit point, d = -5..+5 rho (d > 0 inside the caustic), and prints
binary_magnification_fs's relative error vs VBBL for increasing n.

Result (session 14, s=1.1075, q=0.0079, rho=1e-3, default gate=5): ~1e-6 for
d <= -1 (outside) and 2 <= d < 5 (inside); where the disk straddles the fold
(|d| < 1) up to 2.6% at n=1000, 0.3-1.6% at n=4000, 0.5-0.9% at n=16000 --
~1/sqrt(n) and noisy, because the point-source magnification jumps and
diverges across the caustic line. d=+5 shows the gate edge (point source
used there): -0.1%. Basis for the n=4000 default (below MOA's few-% errors).

Run manually as `python3 scratch/cross_check_mulensmodel_fs.py` from the
project root.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import MulensModel as mm

from lc_models import cassan_caustic, caustic_point, binary_magnification, binary_magnification_fs

RHO = 1e-3
N_LIST = [1000, 4000, 16000]

chain = np.load("scratch/2l1s/cassan/O-03-BLG235_2l1s_cassan_mcmc_chain.npz")
best = chain["samples"][np.argmax(chain["log_probs"])]
sigma_out, s, q = best[1], best[4], best[5]
caustic = cassan_caustic(s, q)
z0 = caustic_point(caustic, sigma_out)
tangent = caustic_point(caustic, sigma_out + 1e-4) - caustic_point(caustic, sigma_out - 1e-4)
normal = 1j * tangent / abs(tangent)


def mm_fs(zeta):
    params = mm.ModelParameters({"t_0": 0.0, "u_0": -zeta.imag, "t_E": 1.0, "alpha": 0.0, "s": s, "q": q, "rho": RHO})
    model = mm.Model(params)
    model.set_magnification_methods([-zeta.real - 1, "VBBL", -zeta.real + 1])
    return model.get_magnification(np.array([-zeta.real]))[0]


print(f"s={s:.4f} q={q:.5f} rho={RHO}")
print(f"{'d/rho':>6} {'A_PS':>10} {'VBBL':>10} " + " ".join(f"{'n=' + str(n):>10}" for n in N_LIST))
for d in [-5, -2, -1, -0.5, 0, 0.5, 1, 2, 5]:
    z = z0 + d * RHO * normal
    ref = mm_fs(z)
    errs = [binary_magnification_fs(z, s, q, RHO, n) / ref - 1 for n in N_LIST]
    print(f"{d:>6} {binary_magnification(z, s, q):>10.3f} {ref:>10.4f} " + " ".join(f"{e:>+10.1e}" for e in errs))
