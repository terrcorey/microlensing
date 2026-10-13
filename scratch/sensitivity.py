"""Planet sensitivity map of a real event (session 27): Gaudi & Sackett (2000) / OGLE-III (arXiv:1602.02519)
style, but Delta chi2 itself instead of a detected / not cut. For every (log s, log q) of the grid and N_ALPHA
source angles: a noise-free 2L1S light curve on the event's timestamps and errors (inject.base_event's
footing: the real run's best single-lens part, fs/fb, LD, errors x K at 2L1S), refit by FSPL + parallax
(restarted Nelder-Mead from the truth's single-lens part, LD fixed at the truth's). The 2L1S is exact (chi2 0),
so Delta chi2 = that FSPL's chi2. Saved per alpha (scratch/sensitivity/<short>.npz), plotted as the median
over alpha (<short>.png) with the injections of scratch/injections/<short>/ on top.
  sbatch slurm/sensitivity.sbatch ../input/<short>.toml      (from slurm/; ~1 h on 24 CPUs for O-03-BLG235)
  .venv/bin/python scratch/sensitivity.py scratch/sensitivity/<short>.npz   (re-plot only: login node is fine)"""
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from multiprocessing import get_context
from pathlib import Path

sys.path.insert(0, "/sharedscratch/gyc1/microlensing")
sys.path.insert(0, "/sharedscratch/gyc1/microlensing/scratch")
import matplotlib.pyplot as plt
import numpy as np
from inject import base_event, flux
from mcmc_fit import polish
from search import FSPL, NB, NF, axis, binary_A, chi2_fspl, fspl_A, fspl_steps
from zoom_utils import save_fig

LOG_S, LOG_Q, N_ALPHA = [-0.7, 0.7, 0.05], [-6.0, -1.0, 0.1], 24  # fineq's grid box, q <= 0.1
SKIP = 10.0  # planet vs planet removed below this: a refit can only lower it, so it's kept as is (and is cheap)
OUT = Path("scratch/sensitivity")
CTX: tuple = ()


def _init(*ctx):
    global CTX
    CTX = ctx


def cell(ls_lq):
    """Pool task: one (log s, log q), every alpha -> [(Delta chi2 vs the planet removed, vs the FSPL refit)]."""
    real, best, fsfb, sigma = CTX
    null = np.array([*best[:NF], *best[NB:]])
    F_null = [flux(i, c, a) for i, c, a in zip(real.instruments, fsfb, fspl_A(real, null))]
    out = []
    for alpha in np.linspace(0, 2 * np.pi, N_ALPHA, endpoint=False):
        th = best.copy()
        th[6], th[7], th[8] = 10 ** ls_lq[0], 10 ** ls_lq[1], alpha
        F = [flux(i, c, a) for i, c, a in zip(real.instruments, fsfb, binary_A(real, th))]
        up = float(sum(np.sum(((f - fn) / sg) ** 2) for f, fn, sg in zip(F, F_null, sigma)))
        if up < SKIP:
            out.append((up, up))
            continue
        syn = real._replace(instruments=[i._replace(flux=f, flux_err=sg) for i, f, sg in zip(real.instruments, F, sigma)])
        # ponytail: one start (the truth's single lens); a second FSPL basin elsewhere would only lower Delta chi2
        _, c = polish(lambda x: chi2_fspl([*x, *best[NB:]], syn), null[:NF], lambda x: fspl_steps(FSPL(*x)))
        out.append((up, min(c, up)))
    return out


def plot(path: Path):
    d = np.load(path)
    med = np.median(d["dchi2"], axis=2)
    fig, ax = plt.subplots(figsize=(9, 6.5))
    pc = ax.pcolormesh(d["log_s"], d["log_q"], np.log10(np.maximum(med, 0.1)).T, shading="nearest", cmap="viridis")
    cs = ax.contour(d["log_s"], d["log_q"], med.T, levels=[10, 25, 50, 100, 1000, 10000], colors="w", linewidths=0.8)
    ax.clabel(cs, fmt="%g", fontsize=8)
    fig.colorbar(pc, label=r"log$_{10}$ median over $\alpha$ of $\Delta\chi^2$ (FSPL + parallax refit $-$ 2L1S)")
    b = d["best"]
    ax.plot(np.log10(b[6]), np.log10(b[7]), "r*", ms=14, mec="k", label="real event's best 2L1S")
    inj = [np.load(f) for f in sorted(Path(f"scratch/injections/{path.stem}").glob("inj*/truth.npz"))]
    inj = [t["theta"] for t in inj if str(t["kind"]) == "planet"]
    if inj:
        ax.plot(np.log10([t[6] for t in inj]), np.log10([t[7] for t in inj]), "wx", ms=7, label="injections")
    ax.set(xlabel=r"log$_{10}$ s", ylabel=r"log$_{10}$ q",
           title=f"{path.stem}: planet sensitivity (noise-free, {d['dchi2'].shape[2]} alphas per cell)")
    ax.legend(loc="lower left", fontsize=8)
    save_fig(fig, path.with_suffix(".png"), dpi=150)
    print(f"saved {path.with_suffix('.png')}")


if __name__ == "__main__":
    if sys.argv[1].endswith(".npz"):
        plot(Path(sys.argv[1]))
        raise SystemExit
    real, best, t0_par, fsfb, sigma = base_event(sys.argv[1])
    ls, lq = axis(LOG_S), axis(LOG_Q)
    cells = list(product(ls, lq))
    res = []
    with ProcessPoolExecutor(mp_context=get_context("spawn"), initializer=_init, initargs=(real, best, fsfb, sigma)) as ex:
        for n, r in enumerate(ex.map(cell, cells, chunksize=2)):
            res.append(r)
            if n % 100 == 0:
                print(f"[sens] {n}/{len(cells)} cells", flush=True)
    res = np.array(res).reshape(len(ls), len(lq), N_ALPHA, 2)
    OUT.mkdir(exist_ok=True)
    path = OUT / f"{real.short_name}.npz"
    np.savez(path, log_s=ls, log_q=lq, alpha=np.linspace(0, 2 * np.pi, N_ALPHA, endpoint=False),
             dchi2_removed=res[..., 0], dchi2=res[..., 1], best=best, t0_par=t0_par)
    plot(path)
