"""Injection-recovery for search.py (session 27): synthetic events on a real event's timestamps and error
bars, run through the pipeline untouched (--stage=refine: no MCMC), scored against the truth.
  make (compute node: one light curve per draw):
    srun --partition=small-short --mem=4G --time=00:15:00 .venv/bin/python scratch/inject.py make \
        input/O-03-BLG235-fineq.toml --grids input/O-03-BLG235.toml input/O-03-BLG235-fineq.toml
    -> scratch/injections/<short>/inj<NN>/: one data file per instrument, one TOML per --grids config
       (its [grid], everything else the base's), truth.npz. Submit every TOML from slurm/:
    for f in ../scratch/injections/<short>/inj*/*.toml; do sbatch search.sbatch $f refine; done
  tally (light: two chi2 per run): .venv/bin/python scratch/inject.py tally scratch/injections/<short>
Base event: the real run's best refined 2L1S (results/<short>/refined.npz) -- t0, u0, tE, rho, piE, LD and
fs/fb (profiled on the real data) kept, s, q, alpha drawn uniform in log s, log q, alpha, kept if the
noise-free planet signal (vs the same with the planet removed, lens at the COM) has Delta chi2 >= MIN_DCHI2.
Nulls: the planet removed. Noise: Gaussian, the real errors x the run's K at 2L1S (K at FSPL carries the
real planet's misfit). The truth's chi2 is taken at its own t0_par, under the search's K and with its LD
prior term: the footing of the refined chi2 it's compared with."""
import argparse
import sys
import tomllib
from pathlib import Path

sys.path.insert(0, "/sharedscratch/gyc1/microlensing")
import numpy as np
from event import ld_penalty, load_event, profile_flux, rescale, with_t0_par
from lc_models import ZERO_POINT_MAG
from search import NB, NF, bic, binary_A, chi2_binary, chi2_fspl, fspl_A

MIN_DCHI2 = 50.0
LOG_S, LOG_Q = (-0.5, 0.5), (-5.0, -1.5)  # inside every config's grid box (log s +-0.7, log q -6..0)
FOUND, CLOSE, DETECT = 1.0, 10.0, 10.0  # refined - truth chi2; Delta BIC (FSPL - 2L1S)


def toml_value(v):
    if isinstance(v, str):
        return f'"{v}"'
    return "[" + ", ".join(map(toml_value, v)) + "]" if isinstance(v, list) else repr(v)


def write_toml(path: Path, cfg: dict):  # the shape of a search config: flat keys, [grid], [[instruments]]
    lines = [f"{k} = {toml_value(v)}" for k, v in cfg.items() if k not in ("grid", "instruments")]
    lines += ["", "[grid]", *(f"{k} = {toml_value(v)}" for k, v in cfg["grid"].items())]
    for i in cfg["instruments"]:
        lines += ["", "[[instruments]]", *(f"{k} = {toml_value(v)}" for k, v in i.items())]
    path.write_text("\n".join(lines) + "\n")


def flux(inst, c, A):
    return c[0] * A + c[1] if inst.kind == "mag" else c[0] * (A - 1)


def base_event(path):
    """The real run's footing (also scratch/sensitivity.py's): (event at its t0_par, best refined 2L1S, t0_par,
    per-instrument fs/fb profiled on the real data, noise sigma = the errors x the run's K at 2L1S)."""
    real = load_event(path)
    res = Path("results") / real.short_name
    r = np.load(res / "refined.npz")
    best, t0_par = r["theta"][np.argmin(r["chi2"])], float(r["plain"][0])
    line = next(l for l in (res / "summary.txt").read_text().splitlines() if l.startswith("K at 2L1S:"))
    k2 = dict((n, float(k)) for n, k in (x.strip().split("=") for x in line.split(":", 1)[1].split(",")))
    real = with_t0_par(real, t0_par)
    fsfb = [profile_flux(i, A)[0] for i, A in zip(real.instruments, binary_A(real, best))]
    return real, best, t0_par, fsfb, [i.flux_err * k2[i.name] for i in real.instruments]


def make(args):
    cfg = tomllib.loads(Path(args.base).read_text())
    real, best, t0_par, fsfb, sigma = base_event(args.base)
    null = np.array([*best[:NF], *best[NB:]])
    A_null = fspl_A(real, null)
    root = Path("scratch/injections") / real.short_name
    assert not root.exists(), f"{root} exists: move it away first (its outputs would mix with the new draws)"
    rng = np.random.default_rng(args.seed)
    planets, draws = [], 0
    while len(planets) < args.n:
        th = best.copy()
        th[6], th[7], th[8] = 10 ** rng.uniform(*LOG_S), 10 ** rng.uniform(*LOG_Q), rng.uniform(0, 2 * np.pi)
        A, draws = binary_A(real, th), draws + 1
        d = sum(np.sum(((flux(i, c, a) - flux(i, c, an)) / sg) ** 2)
                for i, c, a, an, sg in zip(real.instruments, fsfb, A, A_null, sigma))
        if d >= MIN_DCHI2:
            planets.append((th, A, d))
    print(f"[make] {args.n} planets with Delta chi2 >= {MIN_DCHI2} from {draws} draws")
    grids = {Path(g).stem: tomllib.loads(Path(g).read_text())["grid"] for g in args.grids}
    for n, (th, A, d) in enumerate(planets + [(null, A_null, 0.0)] * args.nulls):
        kind = "planet" if n < args.n else "null"
        dr = root / f"inj{n:02d}"
        dr.mkdir(parents=True)
        ins = []
        for i, c, a, sg, ic in zip(real.instruments, fsfb, A, sigma, cfg["instruments"]):
            F = flux(i, c, a) + sg * rng.standard_normal(a.size)
            if i.kind == "mag":  # loader: mag_to_flux gives back F and sigma exactly
                assert np.all(F > 0), f"{i.name}: a noisy flux <= 0 has no magnitude"
                cols = (i.time, ZERO_POINT_MAG - 2.5 * np.log10(F), 2.5 / np.log(10) * sg / F)
            else:
                cols = (i.time, F, sg)
            np.savetxt(dr / f"{i.name}.dat", np.column_stack(cols), fmt="%.10g")
            ins.append({**ic, "path": str(dr / f"{i.name}.dat"), "K": 1.0, "time_fmt": "HJD-2450000"})
        for g, grid in grids.items():
            top = {k: v for k, v in cfg.items() if k not in ("grid", "instruments")}
            write_toml(dr / f"{g}.toml", {**top, "short_name": f"../{dr}/out_{g}", "grid": grid, "instruments": ins})
        back = load_event(str(next(dr.glob("*.toml"))))  # the TOML writer + the mag round trip
        assert all(np.allclose(b.time, i.time) and np.allclose(b.flux_err, sg, rtol=1e-6)
                   for b, i, sg in zip(back.instruments, real.instruments, sigma, strict=True)), "data don't load back"
        np.savez(dr / "truth.npz", theta=th, kind=kind, t0_par=t0_par, dchi2=d, seed=args.seed)
        print(f"[make] {dr.name} {kind}" + (f": s {th[6]:.3f} q {th[7]:.2e} alpha {np.degrees(th[8]):.0f} deg, "
                                            f"Delta chi2 {d:.0f}" if n < args.n else ""))


def tally(args):
    print(f"{'inj':5s} {'grid':22s} {'kind':6s} {'dchi2':>7s} {'s':>6s} {'q':>9s} {'s found':>7s} {'q found':>9s} "
          f"{'found-true':>10s} {'search':>7s} {'dBIC':>7s}")
    score = {}
    for dr in sorted(Path(args.root).glob("inj*")):
        t = np.load(dr / "truth.npz")
        th, planet = t["theta"], str(t["kind"]) == "planet"
        for toml in sorted(dr.glob("*.toml")):
            out = dr / f"out_{toml.stem}"
            if not (out / "fspl.npz").exists():
                print(f"{dr.name:5s} {toml.stem:22s} pending")
                continue
            f, r = np.load(out / "fspl.npz"), np.load(out / "refined.npz")
            ev = rescale(load_event(str(toml)), f["k"])
            truth = (chi2_binary(th, with_t0_par(ev, float(t["t0_par"]))) + ld_penalty(ev, th[NB:]) if planet else
                     chi2_fspl(th, with_t0_par(ev, float(t["t0_par"]))) + ld_penalty(ev, th[NF:]))
            k = int(np.argmin(r["chi2"]))
            x, gap = r["theta"][k], r["chi2"][k] - truth
            evp = with_t0_par(ev, float(r["plain"][0]))  # the search's own t0_par
            dbic = bic(float(f["chi2"]), len(f["theta"]), ev) - bic(chi2_binary(x, evp), len(x), ev)
            verdict = "found" if gap <= FOUND else "close" if gap <= CLOSE else "missed"
            s = score.setdefault(toml.stem, {"found": 0, "close": 0, "missed": 0, "detected": 0, "planets": 0,
                                             "false+": 0, "nulls": 0})
            s["planets" if planet else "nulls"] += 1
            s["detected" if planet else "false+"] += dbic > DETECT
            if planet:
                s[verdict] += 1
            print(f"{dr.name:5s} {toml.stem:22s} {t['kind']!s:6s} {float(t['dchi2']):7.0f} "
                  + (f"{th[6]:6.3f} {th[7]:9.2e}" if planet else f"{'-':>6s} {'-':>9s}")
                  + f" {x[6]:7.3f} {x[7]:9.2e} {gap:10.1f} {verdict if planet else '-':>7s} {dbic:7.1f}")
    for g, s in score.items():
        print(f"[{g}] planets: {s['found']} found, {s['close']} close, {s['missed']} missed of {s['planets']}, "
              f"{s['detected']} detected (dBIC > {DETECT:g}); nulls: {s['false+']} false positives of {s['nulls']}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("make")
    m.add_argument("base", help="the real event's TOML (its results/<short>/ must hold refined.npz + summary.txt)")
    m.add_argument("--grids", nargs="+", required=True, help="configs whose [grid] each injection is searched on")
    m.add_argument("--n", type=int, default=16)
    m.add_argument("--nulls", type=int, default=4)
    m.add_argument("--seed", type=int, default=0)
    sub.add_parser("tally").add_argument("root")
    a = p.parse_args()
    make(a) if a.cmd == "make" else tally(a)
