"""Free limb darkening test (session 23, one-off): does fitting the linear coefficient per band
(I: OGLE + MDM, clear: Auckland, R: FTN; flat on [0, 1]) remove O-05-BLG169's correlated MDM
residuals? FSPL and 2L1S polished from the session-22 fine-q bests (summary.txt) with LD fixed
(the config's) and free, at that run's K at FSPL; prints chi2, the coefficients and, for MDM's
first night (the anomaly), runs-test z, lag-1 autocorrelation and chi2/pt. Writes
scratch/ld_test/O-05-BLG169_mdm_residuals.png. Never on the login node:

    sbatch --partition=small-short --cpus-per-task=4 --mem=4G --time=06:00:00 --job-name=ld-test \\
        --output=slurm/output/slurm-ld-test-%j.out --export=ALL,PYTHONUNBUFFERED=1,OMP_NUM_THREADS=1 \\
        --wrap ".venv/bin/python scratch/ld_test.py"
"""
import sys
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib.pyplot as plt
import numpy as np

from event import flux_residuals, load_event, rescale
from mcmc_fit import polish
from search import FSPL, Binary, binary_A, binary_steps, chi2_binary, chi2_fspl, fit_fspl, fspl_A, fspl_steps

CONFIG, SUMMARY = "input/O-05-BLG169-fineq.toml", Path("results/O-05-BLG169-fineq/session22/summary.txt")
OUT = Path("scratch/ld_test")
BANDS = ("I", "clear", "R")


def read_summary():
    """(FSPL, 2L1S overall best, K at FSPL) from a search.py summary.txt (alpha printed in degrees)."""
    blocks, k = {}, None
    for block in SUMMARY.read_text().split("\n## "):
        title, *lines = block.splitlines()
        vals = {name: float(l.split(" = ")[1].split()[0]) for l in lines
                if (name := l.split(" = ")[0]) in Binary._fields}  # FSPL's fields are a subset
        blocks[title.split(" (")[0]] = vals
        k = next((l for l in lines if l.startswith("K at FSPL:")), k)
    b = blocks["2L1S overall best"]
    b["alpha"] = np.radians(b["alpha"])
    kk = [float(kv.split("=")[1]) for kv in k.removeprefix("K at FSPL: ").split(", ")]
    return FSPL(**{f: blocks["FSPL"][f] for f in FSPL._fields}), Binary(**{f: b[f] for f in Binary._fields}), kk


def with_ld(event, ld):
    return event._replace(instruments=[i._replace(ld=ld[BANDS.index(i.band)]) for i in event.instruments])


def fit(args):
    """Pool task: (model, free LD?, start, event) -> (x, chi2); LD appended to x when free."""
    model, free, p0, event = args
    chi2, steps = (chi2_fspl, fspl_steps) if model is FSPL else (chi2_binary, binary_steps)
    n = len(model._fields)
    if not free:
        return polish(lambda x: chi2(x, event), p0, lambda x: steps(model(*x)))
    ld0 = [next(i.ld for i in event.instruments if i.band == b) for b in BANDS]
    f = lambda x: chi2(x[:n], with_ld(event, x[n:])) if np.all((0 <= x[n:]) & (x[n:] <= 1)) else np.inf
    return polish(f, [*p0, *ld0], lambda x: [*steps(model(*x[:n])), 0.05, 0.05, 0.05])


def wave(r):
    """Runs-test z (negative: fewer sign changes than chance), lag-1 autocorrelation, mean r^2."""
    pos, n = np.sum(r > 0), r.size
    runs = 1 + np.count_nonzero(np.diff(r > 0))
    mu = 2 * pos * (n - pos) / n + 1
    return (runs - mu) / np.sqrt((mu - 1) * (mu - 2) / (n - 1)), np.corrcoef(r[:-1], r[1:])[0, 1], np.mean(r ** 2)


if __name__ == "__main__":  # spawned pool workers re-import this file
    fspl0, best0, k = read_summary()
    _, _, base = fit_fspl(load_event(CONFIG))  # parallax offsets at the run's t0_par (same blind fit)
    event = rescale(base, k)
    print(f"start chi2: FSPL {chi2_fspl(fspl0, event):.2f}, 2L1S {chi2_binary(best0, event):.2f} "
          "(session 22 summary: 575.00, 424.14 -- must match)")
    jobs = [(m, free, p0, event) for m, p0 in ((FSPL, fspl0), (Binary, best0)) for free in (False, True)]
    with ProcessPoolExecutor(len(jobs), mp_context=get_context("spawn")) as ex:
        results = list(ex.map(fit, jobs))

    mdm = [i.name for i in event.instruments].index("MDM")
    t = event.instruments[mdm].time
    night = t < 3492.4
    OUT.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 4))
    for (model, free, _, _), (x, c2) in zip(jobs, results):
        n = len(model._fields)
        ev = with_ld(event, x[n:]) if free else event
        A = (fspl_A if model is FSPL else binary_A)(ev, model(*x[:n]))
        r = np.split(flux_residuals(ev, A), np.cumsum([i.time.size for i in ev.instruments])[:-1])[mdm][night]
        z, ac, c2pt = wave(r)
        name = f"{'FSPL' if model is FSPL else '2L1S'} LD {'free' if free else 'fixed'}"
        ld = dict(zip(BANDS, np.round(x[n:], 3))) if free else "config"
        print(f"{name:15s} chi2 {c2:8.2f}  LD {ld}  MDM night 1: runs z {z:+.2f}, lag-1 {ac:+.2f}, chi2/pt {c2pt:.2f}")
        print(f"{'':15s} {model(*x[:n])}")
        ax.plot(t[night], r, ".-" if free else "o", ms=3, lw=0.6, label=f"{name} (chi2 {c2:.1f})")
    ax.axhline(0, color="k", lw=0.5)
    ax.set(xlabel="HJD - 2450000", ylabel="MDM standardized residual", title="O-05-BLG169 MDM, anomaly night")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "O-05-BLG169_mdm_residuals.png", dpi=200)
    print(f"saved {OUT / 'O-05-BLG169_mdm_residuals.png'}")
