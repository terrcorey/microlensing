import sys, time, VBBinaryLensing as V
s, q, x, y, rho, rel, tol = map(float, sys.argv[1:])
v = V.VBBinaryLensing(); v.RelTol = rel; v.Tol = tol; v.a1 = 0.53
t = time.perf_counter(); a = v.BinaryMag2(s, q, x, y, rho); print(f"{time.perf_counter()-t:.3f}s A={a:.6g}")
