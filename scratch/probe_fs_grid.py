"""Finite-source grid probe (session 23, one-off): what would search.grid_cell() cost, and score, on
O-05-BLG169-fineq with a finite source instead of a point source? rho = t_star / tE per trial, t_star
from the parallax-free FSPL (`plain`) -- the peak measures t_star, not rho. Runs grid_cell() unchanged
on 16 random cells + the cell nearest the resonant crossing the point-source grid missed (s 0.995,
q 4.4e-6) + the control (the full run's best); prints per cell wall time, VBBL call-time tail and
finite- vs point-source score (the saved grid.npz). Never on the login node:

    sbatch --partition=small-short --cpus-per-task=8 --mem=8G --time=03:00:00 --job-name=probe-fs \\
        --output=slurm/output/slurm-probe-fs-%j.out --export=ALL,PYTHONUNBUFFERED=1,OMP_NUM_THREADS=1 \\
        --wrap ".venv/bin/python scratch/probe_fs_grid.py"
"""
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np

import ld_test  # noqa: E402  (scratch/ is sys.path[0] when run as a script)
import search
from event import load_event, rescale

GRID = "results/O-05-BLG169-fineq/grid.npz"
CALLS = []


def init(event, plain, times, n_alpha):
    search._init(event, plain, times, n_alpha)
    point = search.chi2_binary
    t_star = plain.rho * plain.tE

    def timed(theta, ev):  # grid_cell's rho = 0 -> t_star / tE
        theta = tuple(theta)
        if theta[3] == 0.0 and theta[2] > 0:
            theta = (*theta[:3], max(t_star / theta[2], search.RHO_MIN), *theta[4:])
        t = time.perf_counter()
        c = point(theta, ev)
        CALLS.append(time.perf_counter() - t)
        return c
    search.chi2_binary = timed  # grid_cell looks it up at call time


def cell(k_sq):
    CALLS.clear()
    t = time.perf_counter()
    k, c2, _ = search.grid_cell(k_sq)
    return k_sq, min(c2), time.perf_counter() - t, np.array(CALLS)


if __name__ == "__main__":
    ld_test.SUMMARY = Path("results/O-05-BLG169-fineq/session22/summary.txt")
    _, _, k = ld_test.read_summary()
    plain, _, base = search.fit_fspl(load_event("input/O-05-BLG169-fineq.toml"))
    event = rescale(base, k)
    print(f"plain rho {plain.rho:.3g}, tE {plain.tE:.3g} -> t_star {plain.rho * plain.tE:.4f} d")
    times = sorted(search.anomaly_times(event, plain) + list(plain.t0 + plain.tE * np.linspace(-1, 1, 5)))
    g = np.load(GRID)
    ls, lq = g["log_s"], g["log_q"]
    near = lambda s, q: (int(np.argmin(abs(ls - np.log10(s)))), int(np.argmin(abs(lq - np.log10(q)))))
    rng = np.random.default_rng(0)
    picks = [near(0.9948, 4.37e-6), near(1.228, 3.35e-5)] + [tuple(rng.integers(len(a)) for a in (ls, lq)) for _ in range(16)]
    tasks = [((i, j), (10 ** ls[i], 10 ** lq[j])) for i, j in picks]
    with ProcessPoolExecutor(8, mp_context=get_context("spawn"), initializer=init,
                             initargs=(event, plain, times, event.grid["n_alpha"])) as ex:
        for fut in as_completed([ex.submit(cell, t) for t in tasks]):
            ((i, j), _), c2, wall, calls = fut.result()
            print(f"log s {ls[i]:+.2f} log q {lq[j]:+.2f}: finite {c2:9.1f} vs point {g['chi2'][i, j]:9.1f} | "
                  f"{wall:6.0f} s, {calls.size} calls, median {1e3 * np.median(calls):.0f} ms, "
                  f"99% {1e3 * np.percentile(calls, 99):.0f} ms, max {calls.max():.1f} s, "
                  f"top 1% = {calls[calls >= np.percentile(calls, 99)].sum() / calls.sum():.0%} of time")
