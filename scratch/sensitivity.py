"""Planet sensitivity map of a real event (session 27): Gaudi & Sackett (2000) / OGLE-III (arXiv:1602.02519)
style, but Delta chi2 itself instead of a detected / not cut. For every (log s, log q) of the grid and N_ALPHA
source angles: a noise-free 2L1S light curve on the event's timestamps and errors (inject.base_event's
footing: the real run's best single-lens part, fs/fb, LD, errors x K at 2L1S), refit by FSPL + parallax
(inject.dchi2(), shared with the injections' cut). Saved per alpha (scratch/sensitivity/<short>.npz), plotted as the median
over alpha (<short>.png) with the injections of scratch/injections/<short>/ on top: by outcome, one panel per grid,
once `inject.py tally` has written its tally_<grid>.npz (session 28).
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
from inject import ALPHA_TOL, OUTCOMES, Q_TOL, S_TOL, SKIP, U0_TOL, base_event, dchi2
from search import axis
from zoom_utils import save_fig

LOG_S, LOG_Q, N_ALPHA = [-0.7, 0.7, 0.05], [-6.0, -1.0, 0.1], 24  # fineq's grid box, q <= 0.1
OUT = Path("scratch/sensitivity")
CTX: tuple = ()


def _init(*ctx):
    global CTX
    CTX = ctx


def cell(ls_lq):
    """Pool task: one (log s, log q), every alpha -> [(Delta chi2 vs the planet removed, vs the FSPL refit)]."""
    real, best, fsfb, sigma = CTX
    out = []
    for alpha in np.linspace(0, 2 * np.pi, N_ALPHA, endpoint=False):
        th = best.copy()
        th[6], th[7], th[8] = 10 ** ls_lq[0], 10 ** ls_lq[1], alpha
        out.append(dchi2(real, th, best, fsfb, sigma))
    return out


# outcome -> marker (inject.outcome()); wrong fits also get a line from the injected to the found (s, q)
STYLE = {"found": dict(marker="o", color=OUTCOMES["found"][0], mec="k"),
         "wrong fit": dict(marker="^", color=OUTCOMES["wrong fit"][0], mec="k"),
         "not claimed": dict(marker="o", color="none", mec="0.85", mew=1.5),  # light: inject's grey is lost on the map
         "missed": dict(marker="X", color=OUTCOMES["missed"][0], mec="k")}


def plot(path: Path):
    """Median-over-alpha map, one panel per scratch/injections/<short>/tally_<grid>.npz (inject.py tally) with each
    injection marked by its outcome; without a tally, one panel with the injected (s, q) as plain crosses."""
    d = np.load(path)
    med = np.median(d["dchi2"], axis=2)
    root = Path(f"scratch/injections/{path.stem}")
    tallies = sorted(root.glob("tally_*.npz"))
    fig, axes = plt.subplots(1, max(len(tallies), 1), figsize=(9 * max(len(tallies), 1), 7), squeeze=False,
                             layout="constrained")
    extent = (d["log_s"][0], d["log_s"][-1]), (d["log_q"][0], d["log_q"][-1])
    for n, ax in enumerate(axes[0]):
        # below SKIP a cell is the un-refit planet-vs-removed value (an upper bound, inject.dchi2()): one floor colour,
        # not a gradient -- session 27's map showed a fake edge where the refit starts
        pc = ax.pcolormesh(d["log_s"], d["log_q"], np.log10(np.maximum(med, SKIP)).T, shading="nearest", cmap="viridis")
        cs = ax.contour(d["log_s"], d["log_q"], med.T, levels=[10, 25, 50, 100, 1000, 10000], colors="w", linewidths=0.8)
        ax.clabel(cs, fmt="%g", fontsize=8)
        b = d["best"]
        ax.plot(np.log10(b[6]), np.log10(b[7]), "r*", ms=14, mec="k", label="real event's best 2L1S")
        title = ""  # the suptitle says what the map is
        if tallies:
            t = np.load(tallies[n])
            pl, o = t["planet"], t["outcome"]
            lt, lf = np.log10([t["s_true"][pl], t["q_true"][pl]]), np.log10([t["s_found"][pl], t["q_found"][pl]])
            for c, st in STYLE.items():
                m = o[pl] == c
                ax.plot(*lt[:, m], ls="none", ms=9, label=f"{c} ({m.sum()})", **st)
                if c == "wrong fit":  # off-map found (s, q) run to the edge, clipped
                    ax.plot([lt[0, m], np.clip(lf[0, m], *extent[0])], [lt[1, m], np.clip(lf[1, m], *extent[1])],
                            color=OUTCOMES["wrong fit"][0], lw=1)
            nulls = o[~pl]
            title = (f"injections searched on {tallies[n].stem[6:]}\n"
                     f"nulls: {np.sum(nulls == 'false positive')} false positives of {nulls.size}")
        else:
            inj = [np.load(f) for f in sorted(root.glob("inj*/truth.npz"))]
            inj = [t["theta"] for t in inj if str(t["kind"]) == "planet"]
            if inj:
                ax.plot(np.log10([t[6] for t in inj]), np.log10([t[7] for t in inj]), "wx", ms=7, label="injections")
        ax.set(xlabel="s (log scale)", ylabel=r"log$_{10}$ q", xlim=extent[0], ylim=extent[1], title=title)
        ticks = [x for x in (0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5) if extent[0][0] <= np.log10(x) <= extent[0][1]]
        ax.set_xticks(np.log10(ticks), labels=[f"{x:g}" for x in ticks])  # axis stays log10 s: markers unchanged
        ax.legend(loc="lower left", fontsize=8)
    fig.colorbar(pc, ax=axes[0], label=r"log$_{10}$ median over $\alpha$ of $\Delta\chi^2$ (FSPL + parallax refit $-$ 2L1S)")
    fig.suptitle(f"{path.stem}: planet sensitivity (noise-free, {d['dchi2'].shape[2]} alphas per cell)"
                 + (f"\nfound = Delta BIC > 10 and the injected planet: s (or 1/s) within {S_TOL} dex, q {Q_TOL} dex, "
                    f"|u0| {U0_TOL} dex, alpha {ALPHA_TOL:g} deg (vs -alpha if u0 flipped)" if tallies else ""))
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
