"""Session-26 gate for the VBBinaryLensing 3.7 -> VBMicrolensing 5.6.1 swap (run once, while VBBL is still installed):
(1) the session-25 hang geometry (BinaryMag2, source ~0.7 rho from the primary) vs s, both backends, each call in a
    subprocess with a timeout; (2) point-source frame: BinaryMag0 vs lc_models.binary_magnification (pure numpy);
(3) chi2 of the 3 best refined solutions of every finished search under VBM and VBBL (+ refined.npz's VBBL-era
    value, which includes the LD prior), with time per chi2. Pass: |VBM - VBBL| well under 1.
usage: srun --partition=debug --cpus-per-task=1 --mem=4G --time=00:30:00 .venv/bin/python scratch/vbm_check.py"""
import subprocess, sys, time
sys.path.insert(0, "/sharedscratch/gyc1/microlensing")
import numpy as np
import VBBinaryLensing
import lc_models, search
from event import ld_penalty, load_event, rescale, with_t0_par

HANG = ("import time, {m} as V; v = V.{m}(); v.RelTol = 1e-3; v.a1 = 0.53; t = time.perf_counter(); "
        "a = v.BinaryMag2({s}, 5.623413e-4, -0.00557321, 3.94096e-05, 9.682674e-05); "
        "print(f'{{time.perf_counter() - t:.3f}} s  A = {{a:.6g}}')")
for s in (5, 7, 8, 9, 10):
    for m in ("VBBinaryLensing", "VBMicrolensing"):
        try:
            out = subprocess.run([sys.executable, "-c", HANG.format(m=m, s=s)], capture_output=True, text=True,
                                 timeout=60).stdout.strip()
        except subprocess.TimeoutExpired:
            out = "HANG > 60 s"
        print(f"[hang] s = {s:>2}  {m:<16} {out}")

z = np.random.default_rng(1).normal(0, 0.5, (300, 2)) @ np.array([1, 1j])
for s, q in ((1.12, 4e-3), (0.8, 4.4e-5), (1.0, 3e-6), (0.3, 0.03)):
    a, b = lc_models.binary_magnification_vbm(z, s, q, 0, 0), lc_models.binary_magnification(z, s, q)
    print(f"[frame] s {s} q {q}: max |VBM / numpy - 1| = {np.max(np.abs(a / b - 1)):.2e}")

vbm, vbbl = lc_models._VBM, VBBinaryLensing.VBBinaryLensing()
vbbl.LoadESPLTable(str(__import__("pathlib").Path(VBBinaryLensing.__file__).parent / "data" / "ESPL.tbl"))
vbbl.RelTol = 1e-3
for short in ("O-03-BLG235", "O-03-BLG235-fineq", "O-05-BLG169", "O-05-BLG169-noAuckland", "O-05-BLG169-noMDM"):
    r = np.load(f"results/{short}/refined.npz")
    k_line = next(l for l in open(f"results/{short}/summary.txt") if l.startswith("K at FSPL:"))
    k = [float(kv.split("=")[1]) for kv in k_line.split(":", 1)[1].split(",")]
    plain = search.FSPL(*r["plain"])
    ev = with_t0_par(load_event(f"input/{short}.toml"), plain.t0)
    ev = rescale(ev._replace(tE_max=search.TE_FACTOR * plain.tE), k)
    for theta, stored in zip(r["theta"][:3], r["chi2"][:3]):
        p = search.as_binary(theta)
        res = {}
        for name, backend in (("VBM", vbm), ("VBBL", vbbl)):
            lc_models._VBM = backend
            t = time.perf_counter()
            res[name] = (search.chi2_binary(theta, ev), time.perf_counter() - t)
        lc_models._VBM = vbm
        print(f"[chi2] {short:<23} s {p.s:.4f} q {p.q:.3g} rho {p.rho:.3g}: VBM {res['VBM'][0]:.3f} "
              f"({res['VBM'][1]:.2f} s)  VBBL {res['VBBL'][0]:.3f} ({res['VBBL'][1]:.2f} s)  "
              f"diff {res['VBM'][0] - res['VBBL'][0]:+.3f}  | stored {stored:.3f} "
              f"(VBBL + LD prior {res['VBBL'][0] + ld_penalty(ev, theta[search.NB:]):.3f})")
