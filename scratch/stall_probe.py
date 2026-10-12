"""Run one grid cell in-process, logging every chi2_binary call's theta *before* it runs:
the last 'start' line without a matching 'done' is the stalling call.
usage: stall_probe.py <config> <log_s> <log_q> <k1,k2,...> <log>
`plain` comes from results/<short>/refined.npz: the production run's own (a refit lands elsewhere -- rho is unconstrained, session 26)."""
import sys, time
sys.path.insert(0, "/sharedscratch/gyc1/microlensing")
import numpy as np
import search
from event import load_event, rescale, with_t0_par

config, ls, lq, ks, logf = sys.argv[1:]
plain = search.FSPL(*np.load(f"results/{config.split('/')[-1][:-5]}/refined.npz")["plain"])
base = with_t0_par(load_event(config), plain.t0)
event = rescale(base._replace(tE_max=search.TE_FACTOR * plain.tE), [float(k) for k in ks.split(",")])
search._init(event, plain, sorted(search.anomaly_times(event, plain) + list(plain.t0 + plain.tE * np.linspace(-1, 1, 5))),
             event.grid["n_alpha"])
log, orig = open(logf, "w", buffering=1), search.chi2_binary

def logged(theta, ev):
    log.write(f"start {list(map(float, theta))}\n")
    t = time.perf_counter()
    c = orig(theta, ev)
    log.write(f"done {time.perf_counter() - t:.2f} s chi2 {c}\n")
    return c

search.chi2_binary = logged
print(search.grid_cell((0, (10 ** float(ls), 10 ** float(lq))))[1:])
