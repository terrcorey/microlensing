"""One-off: time VBBL BinaryMag2 point by point along one 2L1S theta (as logged by stall_probe.py),
logging each source position *before* the call, so a hang shows which zeta it is.
usage: vbbl_point_probe.py <config> <log> <theta as python list>"""
import sys, time
sys.path.insert(0, "/sharedscratch/gyc1/microlensing")
import numpy as np
from event import load_event
from lc_models import _VBBL, binary_trajectory
from search import as_binary, fit_fspl

config, logf, theta = sys.argv[1], sys.argv[2], eval(sys.argv[3])
p = as_binary(theta)
_, _, ev = fit_fspl(load_event(config))  # parallax offsets (zero here: piE = 0)
log = open(logf, "w", buffering=1)
for inst in ev.instruments:
    _VBBL.a1 = inst.ld
    zeta = binary_trajectory(inst.time, p.t0, p.u0, p.tE, p.alpha, p.piE_N, p.piE_E, inst.dsN, inst.dsE)
    for t, z in zip(inst.time, zeta):
        log.write(f"start {inst.name} t={t:.5f} zeta=({z.real:.6g}, {z.imag:.6g})\n")
        t1 = time.perf_counter()
        a = _VBBL.BinaryMag2(p.s, p.q, z.real, z.imag, p.rho)
        log.write(f"done {time.perf_counter() - t1:.4f} s A={a:.6g}\n")
