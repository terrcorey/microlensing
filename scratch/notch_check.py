"""Session 30 goal 3: is the O-05-BLG169 sensitivity map's dip at log s = 0 smooth or exactly-s=1? (one-off)
  srun --partition=debug --cpus-per-task=18 --mem=18G --time=00:45:00 .venv/bin/python scratch/notch_check.py"""
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from multiprocessing import get_context

sys.path[:0] = ["/sharedscratch/gyc1/microlensing", "/sharedscratch/gyc1/microlensing/scratch"]
import numpy as np
import sensitivity
from inject import base_event

if __name__ == "__main__":
    real, best, _, fsfb, sigma = base_event("input/O-05-BLG169-fineq.toml")
    cells = list(product([-0.05, -0.02, -0.01, -0.001, 0.0, 0.001, 0.01, 0.02, 0.05], [-4.5, -4.0]))
    with ProcessPoolExecutor(mp_context=get_context("spawn"), initializer=sensitivity._init,
                             initargs=(real, best, fsfb, sigma)) as ex:
        res = list(ex.map(sensitivity.cell, cells))
    print(f"{'log s':>7s} {'log q':>6s}  Delta chi2 vs refit: median [16, 84]   vs planet removed: median")
    for (ls, lq), r in sorted(zip(cells, res), key=lambda c: (c[0][1], c[0][0])):
        refit, removed = np.array(r)[:, 1], np.array(r)[:, 0]
        print(f"{ls:7.3f} {lq:6.1f}  {np.median(refit):8.1f} [{np.percentile(refit, 16):7.1f}, "
              f"{np.percentile(refit, 84):7.1f}]   {np.median(removed):10.1f}")
