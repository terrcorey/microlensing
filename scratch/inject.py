"""Injection-recovery for search.py (session 27): synthetic events on a real event's timestamps and error
bars, run through the pipeline untouched (--stage=refine: no MCMC), scored against the truth.
  make (compute node: one light curve per draw):
    srun --partition=debug --mem=4G --time=00:15:00 .venv/bin/python scratch/inject.py make \
        input/O-03-BLG235-fineq.toml --grids input/O-03-BLG235.toml input/O-03-BLG235-fineq.toml
    -> scratch/injections/<short>/inj<NN>/: one data file per instrument, one TOML per --grids config
       (its [grid], everything else the base's), truth.npz. Submit every TOML from slurm/:
    for f in ../scratch/injections/<short>/inj*/*.toml; do sbatch search.sbatch $f refine; done
  tally (light: two chi2 per run): .venv/bin/python scratch/inject.py tally scratch/injections/<short>
    --plot (srun: two plot_fit() per run, session 28): out_<grid>/fit_lc.png = the search's best 2L1S vs its FSPL
    (what diagnose would draw), fit_lc_truth.png = the injected model, planet-removed convention (its signal);
    the text strip carries both parameter sets, chi2 found / truth / FSPL and Delta BIC
    tally also writes <short>/tally_<grid>.npz (one row per run: truth / found s, q, Delta BIC of the found 2L1S and
    of the truth, outcome()) for scratch/sensitivity.py's map
  montage (srun, after tally --plot): .venv/bin/python scratch/inject.py montage scratch/injections/<short>
    -> <short>/montage_<grid>.png per grid: every run's two plots side by side at 1/3 size, 4 injections a row,
    framed in its outcome's colour (OUTCOMES) once tally has written tally_<grid>.npz
Base event: the real run's best refined 2L1S (results/<short>/refined.npz) -- t0, u0, tE, rho, piE, LD and
fs/fb (profiled on the real data) kept, s, q, alpha drawn uniform in log s, log q, alpha, kept if the
noise-free planet signal has Delta chi2 >= MIN_DCHI2 vs an FSPL + parallax refit (dchi2(); session 28 -- was vs
the planet removed, lens at the COM: the refit absorbs ~94% of that, all of it for the wide session-27 draws,
which no search could then detect). truth.npz: dchi2 (vs the refit; session 27's: vs the planet removed) and
dchi2_removed.
Nulls: the planet removed. Noise: Gaussian, the real errors x max(1, the run's K at 2L1S) (K at FSPL carries the
real planet's misfit; base_event() for the floor). The truth's chi2 is taken at its own t0_par, under the search's K and with its LD
prior term: the footing of the refined chi2 it's compared with."""
import argparse
import sys
import tomllib
from pathlib import Path

sys.path.insert(0, "/sharedscratch/gyc1/microlensing")
import numpy as np
from event import ld_penalty, load_event, profile_flux, rescale, with_t0_par
from lc_models import ZERO_POINT_MAG
from mcmc_fit import polish
from search import (FSPL, NB, NF, Binary, bic, binary_A, chi2_binary, chi2_fspl, fmt_params, fspl_A, fspl_steps,
                    labels, plot_fit, REFINERS, refined_name)

MIN_DCHI2 = 50.0
# default draw box (make --log-s / --log-q), inside every config's grid box (log s +-0.7, log q -6..0)
LOG_S, LOG_Q = (-0.5, 0.5), (-5.0, -1.5)
DETECT = 10.0  # Delta BIC (FSPL - 2L1S)
SKIP = 10.0  # dchi2(): planet vs planet removed below this gets no refit (a refit can only lower it)
# outcome(): the found planet is the injected one if s (or its close/wide twin 1/s), q, |u0| (dex) and alpha
# (degrees; vs -alpha when u0 flipped sign -- the (-u0, -alpha) mirror) all lie within these. Session 28: good fits
# sit <= 0.11 / 0.17 / 0.08 dex, <= 4 deg (6 deg off the mirror); Q_TOL was 0.5 (set to pass inj12 on the 0.25
# grid: q 0.39 dex off, a mirror at 0.0 deg) -- tightened, it and inj03 (|u0| 0.35 dex off) are wrong fits
S_TOL, Q_TOL, U0_TOL, ALPHA_TOL = 0.15, 0.2, 0.2, 20.0
# outcome -> (colour, what it means), shared by montage() frames and scratch/sensitivity.py's markers
OUTCOMES = {"found": ("limegreen", "planet claimed, right (s, q)"),
            "wrong fit": ("orange", "planet claimed, wrong (s, q)"),
            "false positive": ("magenta", "planet claimed on a null"),
            "missed": ("red", "no claim, injected planet detectable"),
            "not claimed": ("0.6", "no claim, injected planet undetectable"),
            "null ok": ("royalblue", "no claim on a null")}


def mismatch(true, found):
    """true / found = (u0, s, q, alpha) -> the failed right-planet checks (S_TOL ...), e.g. ["q -0.39 dex"];
    [] = the injected planet. s: the nearer of s and its close/wide twin 1/s; alpha: vs -alpha when u0 flipped."""
    (u_t, s_t, q_t, a_t), (u_f, s_f, q_f, a_f) = true, found
    ls_t, ls_f = np.log10(s_t), np.log10(s_f)
    da = np.degrees(a_f - (a_t if np.sign(u_f) == np.sign(u_t) else -a_t))
    checks = (("s", min(ls_f - ls_t, ls_f + ls_t, key=abs), S_TOL, "dex"), ("q", np.log10(q_f / q_t), Q_TOL, "dex"),
              ("|u0|", np.log10(abs(u_f / u_t)), U0_TOL, "dex"), ("alpha", (da + 180) % 360 - 180, ALPHA_TOL, "deg"))
    return [f"{n} {v:+.{2 if u == 'dex' else 0}f} {u}" for n, v, tol, u in checks if abs(v) > tol]


def outcome(planet, true, found, dbic, dbic_true):
    """true / found = (u0, s, q, alpha). Detected (Delta BIC > DETECT) and the injected planet (no mismatch()) ->
    found, else wrong fit; not detected -> not claimed if the truth itself isn't detectable (its own Delta BIC vs the
    FSPL <= DETECT), else missed."""
    if not planet:
        return "false positive" if dbic > DETECT else "null ok"
    if dbic > DETECT:
        return "wrong fit" if mismatch(true, found) else "found"
    return "missed" if dbic_true > DETECT else "not claimed"


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
    per-instrument fs/fb profiled on the real data, noise sigma = the errors x max(1, the run's K at 2L1S)).
    The floor (session 28): O-05-BLG169's MDM has K 0.64 at 2L1S (the 2L1S overfits it, or its errors are too big)
    -- below 1 would make injections easier to detect than the real errors claim; O-03-BLG235's MOA 0.9965 -> 1."""
    real = load_event(path)
    res = Path("results") / real.short_name
    r = np.load(res / "refined.npz")
    best, t0_par = r["theta"][np.argmin(r["chi2"])], float(r["plain"][0])
    line = next(l for l in (res / "summary.txt").read_text().splitlines() if l.startswith("K at 2L1S:"))
    k2 = dict((n, float(k)) for n, k in (x.strip().split("=") for x in line.split(":", 1)[1].split(",")))
    real = with_t0_par(real, t0_par)
    fsfb = [profile_flux(i, A)[0] for i, A in zip(real.instruments, binary_A(real, best))]
    return real, best, t0_par, fsfb, [i.flux_err * max(1.0, k2[i.name]) for i in real.instruments]


def dchi2(real, th, best, fsfb, sigma, floor=SKIP):
    """Noise-free Delta chi2 of the 2L1S `th` on base_event()'s footing: (vs the planet removed, vs an FSPL +
    parallax refit -- restarted Nelder-Mead from the truth's single-lens part, LD fixed at the truth's). The 2L1S
    is exact (chi2 0), so the latter is that FSPL's chi2. Below `floor` vs the planet removed: no refit, both that."""
    null = np.array([*best[:NF], *best[NB:]])
    F = [flux(i, c, a) for i, c, a in zip(real.instruments, fsfb, binary_A(real, th))]
    F_null = [flux(i, c, a) for i, c, a in zip(real.instruments, fsfb, fspl_A(real, null))]
    up = float(sum(np.sum(((f - fn) / sg) ** 2) for f, fn, sg in zip(F, F_null, sigma)))
    if up < floor:
        return up, up
    syn = real._replace(instruments=[i._replace(flux=f, flux_err=sg) for i, f, sg in zip(real.instruments, F, sigma)])
    # ponytail: one start (the truth's single lens); a second FSPL basin elsewhere would only lower Delta chi2
    _, c = polish(lambda x: chi2_fspl([*x, *best[NB:]], syn), null[:NF], lambda x: fspl_steps(FSPL(*x)))
    return up, min(c, up)


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
        th[6], th[7], th[8] = 10 ** rng.uniform(*args.log_s), 10 ** rng.uniform(*args.log_q), rng.uniform(0, 2 * np.pi)
        (up, d), draws = dchi2(real, th, best, fsfb, sigma, floor=MIN_DCHI2), draws + 1
        if d >= MIN_DCHI2:
            planets.append((th, binary_A(real, th), d, up))
    print(f"[make] {args.n} planets with Delta chi2 >= {MIN_DCHI2} (vs an FSPL refit) from {draws} draws")
    grids = {Path(g).stem: tomllib.loads(Path(g).read_text())["grid"] for g in args.grids}
    for n, (th, A, d, up) in enumerate(planets + [(null, A_null, 0.0, 0.0)] * args.nulls):
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
        np.savez(dr / "truth.npz", theta=th, kind=kind, t0_par=t0_par, dchi2=d, dchi2_removed=up, seed=args.seed,
                 log_s=args.log_s, log_q=args.log_q)  # seed + draw box: the batch is reproducible
        print(f"[make] {dr.name} {kind}" + (f": s {th[6]:.3f} q {th[7]:.2e} alpha {np.degrees(th[8]):.0f} deg, "
                                            f"Delta chi2 {d:.0f} ({up:.0f} vs removed)" if n < args.n else ""))


def tally(args):
    print(f"{'inj':5s} {'grid':22s} {'kind':6s} {'dchi2':>7s} {'s':>6s} {'q':>9s} {'s found':>7s} {'q found':>9s} "
          f"{'found-true':>10s} {'dBIC':>7s} {'truth dBIC':>10s} outcome")
    rows, sfx = {}, "" if args.refine == "nm" else f"_{args.refine}"  # non-nm: own tally / fit_lc files
    for dr in sorted(Path(args.root).glob("inj*")):
        t = np.load(dr / "truth.npz")
        th, planet = t["theta"], str(t["kind"]) == "planet"
        for toml in sorted(dr.glob("*.toml")):
            out = dr / f"out_{toml.stem}"
            if not (out / refined_name(args.refine)).exists():
                print(f"{dr.name:5s} {toml.stem:22s} pending")
                continue
            f, r = np.load(out / "fspl.npz"), np.load(out / refined_name(args.refine))
            ev = rescale(load_event(str(toml)), f["k"])
            evt = with_t0_par(ev, float(t["t0_par"]))
            c2t = chi2_binary(th, evt) if planet else chi2_fspl(th, evt)  # data only: BIC's footing
            truth = c2t + ld_penalty(ev, th[NB:] if planet else th[NF:])  # + LD prior: the refined chi2's footing
            k = int(np.argmin(r["chi2"]))
            x, gap = r["theta"][k], r["chi2"][k] - truth
            evp = with_t0_par(ev, float(r["plain"][0]))  # the search's own t0_par
            dbic = bic(float(f["chi2"]), len(f["theta"]), ev) - bic(chi2_binary(x, evp), len(x), ev)
            dbic_t = bic(float(f["chi2"]), len(f["theta"]), ev) - bic(c2t, len(th), ev) if planet else 0.0
            ut, st, qt, at = (th[1], *th[6:9]) if planet else (np.nan,) * 4
            o = outcome(planet, (ut, st, qt, at), (x[1], *x[6:9]), dbic, dbic_t)
            rows.setdefault(toml.stem, []).append((dr.name, planet, st, qt, at, ut, x[6], x[7], x[8], x[1], dbic,
                                                   dbic_t, o))
            print(f"{dr.name:5s} {toml.stem:22s} {t['kind']!s:6s} {float(t['dchi2']):7.0f} "
                  + (f"{th[6]:6.3f} {th[7]:9.2e}" if planet else f"{'-':>6s} {'-':>9s}")
                  + f" {x[6]:7.3f} {x[7]:9.2e} {gap:10.1f} {dbic:7.1f} "
                  + (f"{dbic_t:10.1f}" if planet else f"{'-':>10s}") + f" {o}")
            if args.plot:
                text = ([f"search best 2L1S ({o})"] + fmt_params(labels(ev, Binary), x),
                        [f"injected {'planet' if planet else 'null (no planet)'}"]
                        + fmt_params(labels(ev, Binary if planet else FSPL), th),
                        [f"{dr.name} on {toml.stem}", f"chi2 found {r['chi2'][k]:.1f}, truth {truth:.1f}, "
                         f"FSPL {float(f['chi2']):.1f}", f"found - truth = {gap:.1f}, Delta BIC = {dbic:.1f}",
                         f"noise-free Delta chi2 {float(t['dchi2']):.0f}, K at FSPL {np.round(f['k'], 2)}",
                         f"FSPL chi2 at the truth's t0_par {chi2_fspl(f['theta'], evt):.1f} (fit_lc_truth.png)"])
                plot_fit(evp, f["theta"], x, float(r["plain"][0]), out / f"fit_lc{sfx}.png", text)
                if planet:
                    plot_fit(evt, f["theta"], th, float(t["t0_par"]), out / "fit_lc_truth.png", text, planet_removed=True)
    for g in rows:
        outs = [r[-1] for r in rows[g]]
        print(f"[{g}] outcomes: " + ", ".join(f"{outs.count(c)} {c}" for c in OUTCOMES if c in outs))
        name, planet, *num, o = zip(*rows[g])
        np.savez(Path(args.root) / f"tally_{g}{sfx}.npz", name=name, planet=planet, s_true=num[0], q_true=num[1],
                 alpha_true=num[2], u0_true=num[3], s_found=num[4], q_found=num[5], alpha_found=num[6],
                 u0_found=num[7], dbic=num[8], dbic_true=num[9], outcome=o)


def montage(args):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from PIL import Image
    from zoom_utils import save_fig
    root, per_row = Path(args.root), 4
    drs = sorted(root.glob("inj*"))
    rows = -(-len(drs) // per_row)
    for toml in sorted(drs[0].glob("*.toml")):
        tf = root / f"tally_{toml.stem}.npz"
        par = [f"{p}_{w}" for w in ("true", "found") for p in ("u0", "s", "q", "alpha")]
        t = np.load(tf) if tf.exists() else {k: [] for k in ("name", "outcome", *par)}
        verdict = dict(zip(t["name"], t["outcome"]))
        usqa = {nm: (v[:4], v[4:]) for nm, *v in zip(t["name"], *(t[k] for k in par))}  # (injected, found) u0 s q alpha
        fig, axes = plt.subplots(rows, 2 * per_row, figsize=(10 * per_row, 5.6 * rows))
        for ax in axes.flat:
            ax.axis("off")
        for n, dr in enumerate(drs):
            null = str(np.load(dr / "truth.npz")["kind"]) == "null"
            o = verdict.get(dr.name)
            for j, name in enumerate(("fit_lc.png", "fit_lc_truth.png")):
                p, ax = dr / f"out_{toml.stem}" / name, axes[n // per_row, 2 * (n % per_row) + j]
                if p.exists():
                    ax.imshow(np.asarray(Image.open(p).reduce(3)))  # 3000 px -> 1000: full size is ~5 GB of RGBA
                what = "search best" if j == 0 else "injected: null (no planet)" if null else "injected"
                title = f"{dr.name} {what}" + ("" if p.exists() or (j and null) else " -- not run yet")
                if o is not None:  # a frame in the outcome's colour round both tiles of the injection
                    ax.axis("on")
                    ax.set(xticks=[], yticks=[])
                    for sp in ax.spines.values():
                        sp.set(color=OUTCOMES[o][0], linewidth=8)
                    title += f" [{o}]"
                    if not (j and null):  # search tile: the found (s, q, alpha); injected tile: the injected one
                        uv, sv, qv, av = usqa[dr.name][1 - j]
                        title += f"\ns {sv:.3f}, q {qv:.2e}, " + r"$\alpha$" + f" {np.degrees(av) % 360:.0f}" + r"$^\circ$"
                    if j == 0 and o == "wrong fit":  # which right-planet checks it failed
                        title += "\nwrong: " + ", ".join(mismatch(*usqa[dr.name]))
                ax.set_title(title, fontsize=14, color="k" if o is None else OUTCOMES[o][0], fontweight="bold")
        step = tomllib.loads(toml.read_text())["grid"]["log_q"][2]
        fig.suptitle(f"{root.name} injections searched on {toml.stem} (log q step {step})", fontsize=22)
        if verdict:
            fig.legend(handles=[Patch(color=c, label=f"{k}: {m} ({list(verdict.values()).count(k)})")
                                for k, (c, m) in OUTCOMES.items()], loc="upper center", ncol=6, fontsize=16,
                       bbox_to_anchor=(0.5, 0.975), frameon=False)
        fig.tight_layout(rect=(0, 0, 1, 0.95))  # room for the suptitle + legend
        save_fig(fig, root / f"montage_{toml.stem}.png", dpi=200)
        plt.close(fig)
        print(f"saved {root / f'montage_{toml.stem}.png'}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("make")
    m.add_argument("base", help="the real event's TOML (its results/<short>/ must hold refined.npz + summary.txt)")
    m.add_argument("--grids", nargs="+", required=True, help="configs whose [grid] each injection is searched on")
    m.add_argument("--n", type=int, default=16)
    m.add_argument("--nulls", type=int, default=4)
    m.add_argument("--seed", type=int, default=0)
    m.add_argument("--log-s", type=float, nargs=2, default=LOG_S, help="uniform draw range of log s")
    m.add_argument("--log-q", type=float, nargs=2, default=LOG_Q, help="uniform draw range of log q")
    t = sub.add_parser("tally")
    t.add_argument("root")
    t.add_argument("--plot", action="store_true", help="+ fit_lc.png / fit_lc_truth.png per run (srun)")
    t.add_argument("--refine", choices=REFINERS, default="nm", help="score search.py --refine=<method>'s results")
    sub.add_parser("montage").add_argument("root")
    a = p.parse_args()
    {"make": make, "tally": tally, "montage": montage}[a.cmd](a)
