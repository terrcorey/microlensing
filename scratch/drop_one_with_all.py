"""One-off: a drop-one run's best fits (fixed), scored and plotted against *all* instruments. The dropped
instrument's fs/fb are profiled against the fixed model -- it never influenced the geometry. Errors use the
full run's K at FSPL, so the total chi2 compares directly with the full run's summary. A band only the
dropped instrument has (e.g. Auckland's clear) gets its LD prior mean. (Was noftn_with_ftn.py, session 23.)
usage: drop_one_with_all.py <short, e.g. O-05-BLG169-noAuckland> [--dir results/<short>] [--xlim lo,hi] [--t0-par T] [--polish]
(--polish: then search.refine() on all data from there -- everything + every band's LD free, chi2 + LD prior --,
scored/plotted as <short>_with_all_polished.png, theta in ..._polished.npz)
(--xlim: the anomaly panel's time range, else plot_fit's automatic window; --t0-par: for pre-session-23
refined.npz files without `plain`, from the run's log line "[fspl] no parallax: ... t0=...")
-> scratch/drop_one/<short>_with_all.png"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/sharedscratch/gyc1/microlensing")
from event import flux_residuals, load_event, rescale, with_t0_par
import search
from search import FSPL, NB, Binary, binary_A, chi2_binary, fspl_A, labels, param_dir, plot_fit

sys.path.insert(0, "/sharedscratch/gyc1/microlensing/scratch")
from ld_test import wave

ap = argparse.ArgumentParser()
ap.add_argument("short")
ap.add_argument("--dir")
ap.add_argument("--xlim", type=lambda x: tuple(map(float, x.split(","))))
ap.add_argument("--t0-par", type=float)
ap.add_argument("--polish", action="store_true")
args = ap.parse_args()
short = args.short
full = short.split("-no")[0]
out = Path(args.dir or f"results/{short}")
t0_par = args.t0_par or float(np.load(out / "refined.npz")["plain"][0])  # the drop-one run's parallax-free t0
drop = load_event(f"input/{short}.toml")
event = with_t0_par(load_event(f"input/{full}.toml"), t0_par)


def summary(path, header):
    """`name = value` lines of one summary.txt section (+/- and units dropped), and its K at FSPL."""
    lines = path.read_text().splitlines()
    i = lines.index(header) + 1
    sec = {l.split(" = ")[0]: float(l.split(" = ")[1].split()[0]) for l in lines[i:lines.index("", i)] if " = " in l}
    k = dict(kv.split("=") for kv in next(l for l in lines if l.startswith("K at FSPL:"))[11:].split(", "))
    return sec, {n: float(v) for n, v in k.items()}


def on_all(model, d):
    """Fit vector for `event`: model fields, then one LD per event.ld_prior band (prior mean if not in d)."""
    return [*(d[f] for f in model._fields), *(d.get(f"ld_{b}", m) for b, m, _ in event.ld_prior)]


fspl, _ = summary(param_dir(out) / "summary.txt", "## FSPL (parallax)")
_, k_full = summary(param_dir(Path(f"results/{full}")) / "summary.txt", "## FSPL (parallax)")
event = rescale(event, [k_full[i.name] for i in event.instruments])

# overall best = refined[0] or any chain's best sample (as diagnose() does)
r = np.load(out / "refined.npz")
cands = [(r["chi2"][0], dict(zip(labels(drop, Binary), r["theta"][0])))]
for f in param_dir(out).glob("mcmc_mode*_chain.npz"):
    d = np.load(f)
    i = np.unravel_index(np.argmax(d["log_prob"]), d["log_prob"].shape)
    cands.append((-2 * d["log_prob"][i], dict(zip(d["labels"], d["chain"][i]))))
best = on_all(Binary, min(cands, key=lambda c: c[0])[1])
fspl = on_all(FSPL, fspl)
tag = ""
if args.polish:  # ponytail: one start, no mirror and no tE ceiling (production's refine() has both): a local polish
    search._init(event)
    print(f"[polish] start chi2 {chi2_binary(best, event):.2f}")
    best, c = search.refine(best)
    tag = "_polished"
    np.savez(f"scratch/drop_one/{short}_with_all{tag}.npz", theta=best, labels=labels(event, Binary), chi2_plus_prior=c)
    print(f"[polish] chi2 + LD prior {c:.2f}")

lines, edges = [], np.cumsum([0] + [i.time.size for i in event.instruments])
for name, A in (("FSPL", fspl_A(event, fspl)), ("2L1S", binary_A(event, best))):
    res = flux_residuals(event, A)
    lines.append(f"{name}: " + ", ".join(f"{i.name} {np.sum(res[a:b] ** 2):.1f}/{b - a}"
                                         for i, a, b in zip(event.instruments, edges, edges[1:])))
    print(lines[-1])
c2 = chi2_binary(best, event)
print(f"2L1S chi2 all = {c2:.2f} (full run's best: see results/{full}/summary.txt)  {Binary(*best[:NB])}")
print(f"LD {[f'{b}={a:.3f}' for (b, _, _), a in zip(event.ld_prior, best[NB:])]}")

# runs test on the 2L1S residuals, per instrument (time-sorted at load) + MDM's first night (the anomaly)
res = np.split(flux_residuals(event, binary_A(event, best)), edges[1:-1])
sets = [(i.name, x) for i, x in zip(event.instruments, res)]
sets += [("MDM night 1", x[i.time < 3492.4]) for i, x in zip(event.instruments, res) if i.name == "MDM"]
print(f"{'2L1S residuals':15s} {'N':>4s} {'runs z':>7s} {'lag-1':>6s} {'chi2/pt':>7s}")
for name, x in sets:
    print(f"{name:15s} {x.size:4d} {wave(x)[0]:+7.2f} {wave(x)[1]:+6.2f} {wave(x)[2]:7.2f}")

text = ([f"{short} best fits, all data" + (", 2L1S polished on all data" if args.polish else ""),
         "(FSPL geometry fixed; dropped fs/fb profiled)", f"K from {full}"],
        [f"chi2 / N per instrument (2L1S total {c2:.2f})"] + lines, [])
plot_fit(event, fspl, best, t0_par, Path(f"scratch/drop_one/{short}_with_all{tag}.png"), text, args.xlim)
