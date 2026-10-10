"""One-off: O-05-BLG169-noFTN's best fits (fixed), plotted against all four instruments incl. FTN.
FTN's fs/fb are profiled against the fixed noFTN model -- FTN never influenced the geometry."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/sharedscratch/gyc1/microlensing")
from event import flux_residuals, load_event, rescale, with_t0_par
from search import FSPL, Binary, as_binary, binary_A, chi2_binary, fit_fspl, fspl_A, plot_fit

sys.path.insert(0, "/sharedscratch/gyc1/microlensing/scratch")
from ld_test import wave

out = Path("results/O-05-BLG169-noFTN/pre24")  # pre-session-24 (fixed LD) backup
plain, _, _ = fit_fspl(load_event("input/O-05-BLG169-noFTN.toml"))  # t0_par = noFTN's parallax-free t0
event = with_t0_par(load_event("input/O-05-BLG169.toml"), plain.t0)
# K at FSPL from noFTN's summary.txt; FTN unscaled
event = rescale(event, [1.133951068605987, 1.2499401086722035, 0.987000534018946, 1.0])
assert [i.name for i in event.instruments] == ["OGLE", "MDM", "Auckland", "FTN"]

# overall best = refined[0] or any chain's best sample (as diagnose() does); summary.txt FSPL
r = np.load(out / "refined.npz")
cands = [(r["chi2"][0], r["theta"][0])]
for f in out.glob("mcmc_mode*_chain.npz"):
    d = np.load(f)
    i = np.unravel_index(np.argmax(d["log_prob"]), d["log_prob"].shape)
    cands.append((-2 * d["log_prob"][i], d["chain"][i]))
best = Binary(*as_binary(min(cands, key=lambda c: c[0])[1]))
fspl = FSPL(3491.881602616428, 0.00114597018663635, 45.941966222521856, 0.00032427499848147143,
            1.972553276386617, -0.012306872851097439)

lines = []
for name, A in (("FSPL", fspl_A(event, fspl)), ("2L1S", binary_A(event, best))):
    r = flux_residuals(event, A)
    edges = np.cumsum([0] + [i.time.size for i in event.instruments])
    per = {i.name: (np.sum(r[a:b] ** 2), b - a) for i, a, b in zip(event.instruments, edges, edges[1:])}
    lines.append(f"{name}: " + ", ".join(f"{k} {c:.1f}/{n}" for k, (c, n) in per.items()))
    print(lines[-1])
print(f"2L1S chi2 all = {chi2_binary(best, event):.2f}  {best}")

# runs test on the 2L1S residuals, per instrument (time-sorted at load) + MDM's first night (the anomaly)
r = np.split(flux_residuals(event, binary_A(event, best)), np.cumsum([i.time.size for i in event.instruments])[:-1])
mdm = event.instruments[1]
sets = [(i.name, x) for i, x in zip(event.instruments, r)] + [("MDM night 1", r[1][mdm.time < 3492.4])]
print(f"{'2L1S residuals':15s} {'N':>4s} {'runs z':>7s} {'lag-1':>6s} {'chi2/pt':>7s}")
for name, x in sets:
    print(f"{name:15s} {x.size:4d} {wave(x)[0]:+7.2f} {wave(x)[1]:+6.2f} {wave(x)[2]:7.2f}")

text = (["noFTN best fits, FTN overplotted", "(FTN fs/fb profiled, geometry fixed)"],
        ["chi2 / N per instrument"] + lines, [])
plot_fit(event, fspl, best, plain.t0, Path(sys.argv[1]), text)
