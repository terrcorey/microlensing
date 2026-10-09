"""Free limb darkening test (session 23, one-off): does fitting the linear coefficient per band
(I: OGLE + MDM, clear: Auckland, R: FTN; flat on [0, 1]) remove O-05-BLG169's correlated MDM
residuals? FSPL and 2L1S polished from the session-22 fine-q bests (summary.txt) with LD fixed
(the config's) and free, at that run's K at FSPL; prints chi2, the coefficients and, for MDM's
first night (the anomaly), runs-test z, lag-1 autocorrelation and chi2/pt. Writes
scratch/ld_test/O-05-BLG169_mdm_residuals.png. Never on the login node:

    sbatch --partition=small-short --cpus-per-task=4 --mem=4G --time=06:00:00 --job-name=ld-test \\
        --output=slurm/output/slurm-ld-test-%j.out --export=ALL,PYTHONUNBUFFERED=1,OMP_NUM_THREADS=1 \\
        --wrap ".venv/bin/python scratch/ld_test.py"

--free0: the 2L1S + per-band LD job alone, LD started from 0 (does session 23's I-band 0.125 depend on the start?).

--scan instead: the 2L1S geometry fixed at the session-22 best, one coefficient shared by every instrument,
scanned over [0, 1] (seconds, plus the blind FSPL fit; scratch/ld_test/O-05-BLG169_ld_scan.png):

    srun --partition=small-short --cpus-per-task=2 --mem=4G --time=00:15:00 .venv/bin/python scratch/ld_test.py --scan
"""
import sys
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib.pyplot as plt
import numpy as np

from event import flux_residuals, load_event, profile_flux, rescale
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


def scan(event, best):
    """Shared LD coefficient a over [0, 1], geometry fixed (fs/fb still profiled). With the geometry fixed each
    instrument's chi2 depends only on its own fs/fb, so the same scan gives MDM's own best a too."""
    names = [i.name for i in event.instruments]
    mdm = names.index("MDM")
    night = event.instruments[mdm].time < 3492.4
    a = np.linspace(0, 1, 21)
    chi2, mdm_r = [], []
    for ai in a:
        ev = event._replace(instruments=[i._replace(ld=ai) for i in event.instruments])
        r = np.split(flux_residuals(ev, binary_A(ev, best)), np.cumsum([i.time.size for i in ev.instruments])[:-1])
        chi2.append([np.sum(x ** 2) for x in r])
        mdm_r.append(r[mdm][night])
    chi2 = np.array(chi2)
    print(f"{'a':>5s} {'total':>9s} " + " ".join(f"{n:>9s}" for n in names) + "   MDM night 1: runs z, lag-1, chi2/pt")
    for ai, c, r in zip(a, chi2, mdm_r):
        print(f"{ai:5.2f} {c.sum():9.2f} " + " ".join(f"{x:9.2f}" for x in c) + "   {:+.2f}, {:+.2f}, {:.2f}".format(*wave(r)))
    best_all, best_mdm = chi2.sum(1).argmin(), chi2[:, mdm].argmin()
    print(f"best shared a = {a[best_all]:.2f}, best for MDM alone = {a[best_mdm]:.2f} (grid step 0.05)")
    # MDM's anomaly night: data in its own flux, each a's model at its own profiled fs/fb (on all MDM points)
    inst = event.instruments[mdm]
    t, F, err = inst.time[night], inst.flux[night], inst.flux_err[night]
    td = np.linspace(t.min(), t.max(), 2000)
    dense = inst._replace(time=td, dsN=np.interp(td, inst.time, inst.dsN), dsE=np.interp(td, inst.time, inst.dsE))

    def model(ai):  # MDM model flux on the data's times and on td
        ev = event._replace(instruments=[i._replace(ld=ai) for i in (inst, dense)])
        A_data, A_dense = binary_A(ev, best)
        fs, fb = profile_flux(ev.instruments[0], A_data)[0]
        return fs * A_data + fb, fs * A_dense + fb

    ref_data, ref_dense = model(a[best_all])
    sig = np.median(err)
    fig, (top, bot) = plt.subplots(2, 1, figsize=(9, 7), sharex=True, height_ratios=[2, 1])
    top.errorbar(t, F, err, fmt=".", color="k", ms=3, lw=0.5, label="MDM")
    bot.plot(t, (F - ref_data[night]) / err, ".", color="k", ms=3)
    for ai in (0, 0.25, a[best_all], 0.75, 1):
        k = np.abs(a - ai).argmin()
        m_dense = model(ai)[1]
        label = f"a = {ai:.2f}{' (best)' if k == best_all else ''}: total chi2 {chi2[k].sum():.1f}, MDM {chi2[k, mdm]:.1f}"
        top.plot(td, m_dense, lw=1, label=label)
        bot.plot(td, (m_dense - ref_dense) / sig, lw=1)
    bot.axhline(0, color="k", lw=0.5)
    top.set(ylabel="MDM flux (instrument scale)",
            title="O-05-BLG169 MDM, anomaly night: 2L1S geometry fixed, shared LD a")
    bot.set(xlabel="HJD - 2450000", ylabel=f"minus best model / sigma\n(curves: median sigma)")
    top.legend(fontsize=8)
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "O-05-BLG169_ld_scan.png", dpi=200)
    print(f"saved {OUT / 'O-05-BLG169_ld_scan.png'}")


if __name__ == "__main__":  # spawned pool workers re-import this file
    fspl0, best0, k = read_summary()
    _, _, base = fit_fspl(load_event(CONFIG))  # parallax offsets at the run's t0_par (same blind fit)
    event = rescale(base, k)
    print(f"start chi2: FSPL {chi2_fspl(fspl0, event):.2f}, 2L1S {chi2_binary(best0, event):.2f} "
          "(session 22 summary: 575.00, 424.14 -- must match)")
    if "--scan" in sys.argv:
        scan(event, best0)
        sys.exit()
    jobs = [(m, free, p0, event) for m, p0 in ((FSPL, fspl0), (Binary, best0)) for free in (False, True)]
    tag = "_free0" if "--free0" in sys.argv else ""
    if tag:  # 2L1S only, per-band LD free from a uniform disk instead of the config's values
        jobs = [(Binary, True, best0, with_ld(event, [0.0] * len(BANDS)))]
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
    fig.savefig(OUT / f"O-05-BLG169_mdm_residuals{tag}.png", dpi=200)
    print(f"saved {OUT / f'O-05-BLG169_mdm_residuals{tag}.png'}")
