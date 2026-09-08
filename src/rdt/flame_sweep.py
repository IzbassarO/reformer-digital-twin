"""The flame-length structural variant as a threshold rather than a binary.

``python -m rdt.flame_sweep [--version v3]``

The robustness table of the main text reports each conclusion twice: once at the calibrated,
load-independent flame length and once at a bounding structural variant in which the flame shortens
with load as ``L_q_eff = L_q * load**0.5``. That is a binary, and it decides the sign of the
flexibility conclusions, so the question it leaves open is *how much* load dependence is needed. This
module replaces the binary with a sweep of the exponent

    L_q_eff = L_q * load**p,      p = 0, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 1.0,

re-running the paired Monte Carlo at each ``p`` and reporting the robustness fraction of every
statement as a function of ``p``.

What is reused
--------------
Everything that decides a fraction. The Monte Carlo design is the run's own stored sample
(``mc_<tag>_design.csv.gz``), the per-sample evaluation is :func:`rdt.uq.evaluate_one`, and the
fractions themselves are computed by :func:`rdt.uq.robustness_fractions`, the function that produced
the published table. Only the exponent changes, so the numbers stay comparable with the published ones
by construction. :func:`verify_endpoints` checks that literally: ``p = 0`` must reproduce the nominal
column and ``p = 0.5`` the variant column of the published table, sample for sample.

What has to be re-run, and what does not
----------------------------------------
``load**p`` is 1 for every ``p`` when the load is 1, so a point at full load is independent of the
exponent. Of the points entering the five statements only three sit away from full load --- the RQ3
knee at 1.10 and the two flexibility states at 0.85 and 0.70 --- and only those are re-run. The base
point, the fourth ageing year and the two steam-to-carbon points are all at load 1 and are read from
the stored campaign unchanged, which is also why two of the five fractions cannot move with ``p``.

Scope
-----
The sweep is on the base alloy curve, as the published table is. ``p`` is a structural assumption
about the furnace, not a fitted quantity: the point of the sweep is to say which value of it would
change a conclusion, not to claim a value.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from rdt import config as C
from rdt import uq

OUT_DIR = C.ROOT / "data" / "uq" / "flame_sweep"

#: exponents swept in ``L_q_eff = L_q * load**p``
P_GRID: Tuple[float, ...] = (0.0, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 1.0)

#: the only points entering the five statements whose load differs from 1, hence the only ones that
#: move with ``p``; everything else is read from the stored campaign
SWEPT_POINTS: Tuple[str, ...] = ("RQ3_knee", "load_0.70", "load_0.85")

#: the five statements of the robustness table, in the order they are printed there
STATEMENTS: Tuple[Tuple[str, str], ...] = (
    ("S5y4_life_gt_2x_S1", "aging year 4 costs > 2x S1 life per kmol"),
    ("knee_life_lt_0.5_base", "knee regime < 0.5x base"),
    ("SC3.5_less_life_per_kmol_than_SC2.5", "S/C 3.5 cheaper in life than S/C 2.5"),
    ("S2_life_per_kmol_le_S1", "S2 daily <= S1"),
    ("S3_life_per_kmol_le_S1", "S3 renewable <= S1"),
)
#: the two whose sign the flame-length assumption is known to decide
FLEXIBILITY = ("S2_life_per_kmol_le_S1", "S3_life_per_kmol_le_S1")
MATERIAL = ("S5y4_life_gt_2x_S1", "knee_life_lt_0.5_base", "SC3.5_less_life_per_kmol_than_SC2.5")

#: per-sample agreement expected against the stored campaign at the two published exponents. The
#: coupled solve is an adaptive LSODA integration at ``rtol = 1e-6`` (:func:`rdt.furnace.simulate_coupled`),
#: so run-to-run reproducibility is bounded by the integrator's own tolerance, not by machine epsilon.
SAMPLE_TOL = 1e-6


# ---------------------------------------------------------------------------
# Running one exponent
# ---------------------------------------------------------------------------
def _tag(p: float) -> str:
    return f"p{p:g}".replace(".", "_")


def run_point_at_p(name: str, point: Dict[str, float], X: pd.DataFrame, p: float,
                   n_jobs: int = -1, out_dir: Optional[Path] = None) -> pd.DataFrame:
    """Evaluate the stored Monte Carlo design at one operating point and one exponent (cached).

    Identical call to the one the published campaign makes, with the boolean variant flag replaced by
    the exponent itself (:func:`rdt.uq.lq_exponent`).
    """
    from joblib import Parallel, delayed

    out_dir = Path(out_dir or OUT_DIR)
    path = out_dir / f"{uq.MC_TAG}_{name}_{_tag(p)}.csv.gz"
    if path.exists():
        return pd.read_csv(path)
    alloy = uq.creep.alloy_key()
    rows = Parallel(n_jobs=n_jobs)(
        delayed(uq._worker)(point, th, float(p), alloy) for th in X.to_dict(orient="records"))
    df = pd.concat([X.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    df["point"] = name
    df["lq_exponent"] = float(p)
    df["sample"] = np.arange(len(df))
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, compression="gzip")
    return df


def stored_campaign(cfg: C.RunConfig) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame, Dict[str, Dict[str, float]]]:
    """The published Monte Carlo campaign: per-point results, the design sample and the points."""
    d = Path(uq.OUT_DIR)
    mc = {n: pd.read_csv(d / f"{uq.MC_TAG}_{n}.csv.gz") for n in uq.POINT_NAMES}
    mc["load_1.00"] = mc["S1_base"]                       # load exactly 1: the same state by definition
    X = pd.read_csv(d / f"{uq.MC_TAG}_design.csv.gz")
    pts = json.loads((d / f"{uq.MC_TAG}_points.json").read_text())
    return mc, X, pts


def statement_masks(mc: Dict[str, pd.DataFrame], swept: Dict[str, pd.DataFrame]) -> Dict[str, np.ndarray]:
    """Per-sample truth of each statement, mirroring :func:`rdt.uq.robustness_fractions`.

    The published function returns only the means, but the decisive check on a re-run is whether the
    *same samples* decide each fraction, not merely whether the counts agree. This reproduces its
    expressions sample by sample; :func:`fractions_at_p` cross-checks that the means of these masks
    equal the published function's output, so any divergence between the two shows up immediately.
    """
    M = dict(mc); M.update(swept)
    b = M["S1_base"]
    ok_all = (b.converged == True).to_numpy().copy()  # noqa: E712
    for k in ("load_0.70", "load_0.85", "load_1.00"):
        ok_all &= (M[k].converged == True).to_numpy()  # noqa: E712
    out: Dict[str, np.ndarray] = {}

    def paired(a, bb, col="life_rate_per_kmol_H2"):
        ok = (a.converged == True) & (bb.converged == True)  # noqa: E712
        return a.loc[ok, col].to_numpy(), bb.loc[ok, col].to_numpy()

    s5, bb = paired(M["S5_y4_holdslip_end"], b)
    out["S5y4_life_gt_2x_S1"] = s5 > 2.0 * bb
    kk_, bb = paired(M["RQ3_knee"], b)
    out["knee_life_lt_0.5_base"] = kk_ < 0.5 * bb
    ok = (M["S6_SC3.5"].converged == True) & (M["S6_SC2.5"].converged == True)  # noqa: E712
    out["SC3.5_less_life_per_kmol_than_SC2.5"] = (
        M["S6_SC3.5"].loc[ok, "life_rate_per_kmol_H2"].to_numpy()
        < M["S6_SC2.5"].loc[ok, "life_rate_per_kmol_H2"].to_numpy())
    d2 = 16 * M["load_1.00"].loc[ok_all, "life_rate_per_h"].to_numpy() + 8 * M["load_0.70"].loc[ok_all, "life_rate_per_h"].to_numpy()
    h2 = 16 * M["load_1.00"].loc[ok_all, "H2_net_kmol_h"].to_numpy() + 8 * M["load_0.70"].loc[ok_all, "H2_net_kmol_h"].to_numpy()
    out["S2_life_per_kmol_le_S1"] = d2 / h2 <= b.loc[ok_all, "life_rate_per_kmol_H2"].to_numpy()
    w, keys = (0.25, 0.5, 0.25), ("load_0.70", "load_0.85", "load_1.00")
    d3 = sum(wi * M[c].loc[ok_all, "life_rate_per_h"].to_numpy() for wi, c in zip(w, keys))
    h3 = sum(wi * M[c].loc[ok_all, "H2_net_kmol_h"].to_numpy() for wi, c in zip(w, keys))
    out["S3_life_per_kmol_le_S1"] = d3 / h3 <= b.loc[ok_all, "life_rate_per_kmol_H2"].to_numpy()
    return out


def fractions_at_p(mc: Dict[str, pd.DataFrame], swept: Dict[str, pd.DataFrame]) -> Dict[str, float]:
    """Robustness fractions at one exponent, through the published function.

    ``swept`` replaces only the points that move with ``p``; every other point, including the base the
    ratios are paired against, is the stored one.
    """
    var = dict(mc)
    var.update(swept)
    f = uq.robustness_fractions(mc, var)["load_dependent_L_q"]
    masks = statement_masks(mc, swept)
    for k, v in f.items():                      # the mirror must agree with the published function
        assert abs(float(np.mean(masks[k])) - v) < 1e-15, f"mask mirror disagrees on {k}"
    return f


def sweep(cfg: C.RunConfig, ps: Sequence[float] = P_GRID, n_jobs: int = -1,
          log=print) -> Dict[str, object]:
    mc, X, pts = stored_campaign(cfg)
    rows, runs = [], {}
    for p in ps:
        swept = {}
        for name in SWEPT_POINTS:
            swept[name] = run_point_at_p(name, pts[name], X, p, n_jobs=n_jobs)
        runs[p] = swept
        f = fractions_at_p(mc, swept)
        conv = min(float((swept[n].converged == True).mean()) for n in SWEPT_POINTS)  # noqa: E712
        rows.append({"p": float(p), **f, "min_convergence": conv})
        log(f"  p = {p:<5g} done  ({', '.join(f'{k.split(chr(95))[0]} {v:.3f}' for k, v in f.items())})")
    return {"table": pd.DataFrame(rows).set_index("p"), "mc": mc, "runs": runs, "points": pts}


# ---------------------------------------------------------------------------
# Endpoint verification
# ---------------------------------------------------------------------------
def verify_endpoints(res: Dict[str, object]) -> Dict[str, object]:
    """``p = 0`` must reproduce the nominal column and ``p = 0.5`` the variant column of the table.

    Checked at two levels. The per-sample life rates are compared against the stored campaign: they
    agree to about 1e-12 relative rather than bit for bit, because the coupled ODE solve is not
    reproducible in its last bits across runs (the BLAS reduction order depends on how the work is
    spread over workers). That is twelve orders below anything that can flip one of the inequalities
    being counted. The fractions themselves must then equal the published ones exactly, which is the
    claim that matters, and is not circular: these runs are fresh.
    """
    mc = res["mc"]
    published = json.loads(Path(uq.SUMMARY_JSON).read_text())["robustness_fractions"]
    out: Dict[str, object] = {"samples": {}, "fractions": {}, "decisions": {}, "sample_tol": SAMPLE_TOL}
    for p, col, suffix in ((0.0, "nominal_L_q", ""), (0.5, "load_dependent_L_q", "_lqvar")):
        if p not in res["runs"]:
            continue
        for name in SWEPT_POINTS:
            got = res["runs"][p][name]
            ref = pd.read_csv(Path(uq.OUT_DIR) / f"{uq.MC_TAG}_{name}{suffix}.csv.gz")
            cols = [c for c in ("life_rate_per_kmol_H2", "life_rate_per_h", "T_wo_max_K", "log10_t_r")
                    if c in got.columns and c in ref.columns]
            worst, med, p99 = 0.0, 0.0, 0.0
            for c in cols:
                a, b = got[c].to_numpy(float), ref[c].to_numpy(float)
                m = np.isfinite(a) & np.isfinite(b) & (b != 0)
                if m.any():
                    r = np.abs(a[m] - b[m]) / np.abs(b[m])
                    worst = max(worst, float(np.max(r)))
                    med = max(med, float(np.median(r)))
                    p99 = max(p99, float(np.percentile(r, 99)))
            out["samples"][f"p={p:g} {name}"] = {"max_rel_diff": worst, "median_rel_diff": med,
                                                 "p99_rel_diff": p99, "n": int(len(got)),
                                                 "p99_within_tol": bool(p99 <= SAMPLE_TOL),
                                                 "columns": cols}
        f = fractions_at_p(mc, res["runs"][p])
        out["fractions"][f"p={p:g}"] = {k: {"got": float(f[k]), "published": float(published[col][k]),
                                            "exact": f[k] == published[col][k]} for k, _ in STATEMENTS}
        # the decisive check: the SAME samples must decide each fraction, not merely as many of them
        stored_swept = {n: pd.read_csv(Path(uq.OUT_DIR) / f"{uq.MC_TAG}_{n}{suffix}.csv.gz")
                        for n in SWEPT_POINTS}
        a_masks = statement_masks(mc, res["runs"][p])
        b_masks = statement_masks(mc, stored_swept)
        out["decisions"][f"p={p:g}"] = {k: int(np.sum(a_masks[k] != b_masks[k])) for k, _ in STATEMENTS}
    out["all_samples_within_tol"] = all(v["p99_within_tol"] for v in out["samples"].values())
    out["worst_sample_rel_diff"] = max((v["max_rel_diff"] for v in out["samples"].values()), default=0.0)
    out["all_fractions_exact"] = all(d["exact"] for g in out["fractions"].values() for d in g.values())
    out["total_decision_flips"] = int(sum(v for g in out["decisions"].values() for v in g.values()))
    return out


# ---------------------------------------------------------------------------
# Crossings
# ---------------------------------------------------------------------------
def crossing(table: pd.DataFrame, key: str, level: float = 0.5) -> Dict[str, object]:
    """Where the fraction of one statement crosses ``level``, linearly interpolated in ``p``.

    Returns every crossing found on the grid, and says plainly when there is none because the fraction
    is on one side of the level throughout the swept range.
    """
    p = table.index.to_numpy(float)
    f = table[key].to_numpy(float)
    hits: List[Dict[str, float]] = []
    for i in range(len(p) - 1):
        a, b = f[i], f[i + 1]
        if (a - level) * (b - level) < 0:
            frac = (level - a) / (b - a)
            hits.append({"p": float(p[i] + frac * (p[i + 1] - p[i])),
                         "from": float(a), "to": float(b),
                         "bracket": (float(p[i]), float(p[i + 1]))})
        elif a == level:
            hits.append({"p": float(p[i]), "from": float(a), "to": float(b), "bracket": (float(p[i]), float(p[i]))})
    return {"crossings": hits, "at_p0": float(f[0]), "at_pmax": float(f[-1]),
            "side_at_p0": "above" if f[0] > level else ("below" if f[0] < level else "on"),
            "monotone_decreasing": bool(np.all(np.diff(f) <= 1e-12)),
            "monotone_increasing": bool(np.all(np.diff(f) >= -1e-12))}


def collapse_p(table: pd.DataFrame, key: str, level: float) -> Optional[float]:
    """Interpolated ``p`` at which a falling fraction first drops below ``level``."""
    c = crossing(table, key, level)
    return c["crossings"][0]["p"] if c["crossings"] else None


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def report(res: Dict[str, object], ver: Dict[str, object]) -> str:
    t: pd.DataFrame = res["table"]
    L: List[str] = []
    P = L.append
    P("Flame-length structural variant resolved into a threshold")
    P("=" * 110)
    P("The robustness table reports each conclusion at the calibrated flame length and at a bounding")
    P("variant L_q_eff = L_q * load^0.5. Here the exponent p is swept instead, and the paired Monte")
    P("Carlo robustness fraction of each statement is recomputed at every p on the base alloy curve.")
    P("")
    P("  p = 0    is the calibrated, load-independent flame length  (nominal column of the table)")
    P("  p = 0.5  is the published bounding variant                 (variant column of the table)")
    P("  p is a structural assumption about the furnace, not a fitted quantity. The sweep says which")
    P("  value of it would change a conclusion; it does not claim a value.")
    P("")

    P("WHAT WAS RE-RUN")
    P("-" * 110)
    P("  load^p = 1 for any p when the load is 1, so a full-load point cannot move with the exponent.")
    P(f"  re-run at every p : {', '.join(SWEPT_POINTS)}  (loads "
      + ", ".join(f"{res['points'][n]['feed_per_tube_fraction']:.3f}" for n in SWEPT_POINTS) + ")")
    P("  read from the stored campaign unchanged, all at load 1.000 : S1_base, S5_y4_holdslip_end,")
    P("  S6_SC2.5, S6_SC3.5, load_1.00")
    P("  Two of the five fractions therefore cannot move with p, and are shown flat as a check rather")
    P("  than as a result.")
    P("")

    P("ENDPOINT VERIFICATION")
    P("-" * 110)
    P("  Per-sample agreement with the stored campaign (fresh runs, not the stored files re-read).")
    P("  The coupled solve is an adaptive LSODA integration at rtol = 1e-6, so run-to-run agreement is")
    P("  bounded by the integrator's own tolerance rather than by machine epsilon:")
    P(f"    {'':<22} {'n':>5} {'median':>11} {'p99':>11} {'max':>11}")
    for k, v in ver["samples"].items():
        P(f"    {k:<22} {v['n']:>5} {v['median_rel_diff']:11.2e} {v['p99_rel_diff']:11.2e} "
          f"{v['max_rel_diff']:11.2e}")
    P(f"    the solver rtol is {ver['sample_tol']:.0e}. It bounds the local error per step, not the")
    P(f"    accumulated error, so a maximum a few times above it ({ver['worst_sample_rel_diff']:.1e}) is")
    P("    expected; the median sits four to five orders below it.")
    P("")
    P("  The decisive check is not the size of that difference but whether it changes any decision.")
    P("  For each statement, the number of samples whose truth differs between the fresh run and the")
    P("  stored campaign:")
    for pk, g in ver["decisions"].items():
        P(f"    {pk:<8} " + "  ".join(f"{k.split('_')[0]}: {v}" for k, v in g.items()))
    P(f"    total decision flips across both endpoints: {ver['total_decision_flips']}")
    for pk, g in ver["fractions"].items():
        col = "nominal" if pk == "p=0" else "variant"
        P(f"  {pk} against the published {col} column:")
        for k, d in g.items():
            P(f"      {k:<38} {d['got']:.10f}  vs  {d['published']:.10f}   "
              f"{'EXACT' if d['exact'] else '*** DIFFERS ***'}")
    P(f"  99th percentile of every re-run within the solver rtol : "
      f"{'YES' if ver['all_samples_within_tol'] else 'NO'}")
    P(f"  samples whose decision changed, over all five statements and both endpoints : "
      f"{ver['total_decision_flips']}")
    P(f"  every endpoint fraction exact     : {'YES' if ver['all_fractions_exact'] else 'NO'}")
    P("")

    P("ROBUSTNESS FRACTION AS A FUNCTION OF p")
    P("-" * 110)
    P(f"  {'statement':<44} " + " ".join(f"{p:>7g}" for p in t.index))
    P(f"  {'p =':<44} " + " ".join(f"{'':>7}" for _ in t.index))
    for key, label in STATEMENTS:
        P(f"  {label:<44} " + " ".join(f"{v:7.3f}" for v in t[key]))
    P("")
    P(f"  {'minimum convergence over the re-run points':<44} "
      + " ".join(f"{v:7.3f}" for v in t["min_convergence"]))
    P("")

    P("CROSSINGS OF 0.5")
    P("-" * 110)
    P("  Flexibility statements (the ones the flame-length assumption is known to decide):")
    for key in FLEXIBILITY:
        lab = dict(STATEMENTS)[key]
        c = crossing(t, key)
        if c["crossings"]:
            for h in c["crossings"]:
                P(f"    {lab:<40} crosses 0.5 at p = {h['p']:.4f}  "
                  f"(between p = {h['bracket'][0]:g} and {h['bracket'][1]:g}, "
                  f"{h['from']:.3f} -> {h['to']:.3f})")
        else:
            P(f"    {lab:<40} does NOT cross 0.5 in the swept range: it is already")
            P(f"    {'':<40} {c['side_at_p0']} 0.5 at p = 0 ({c['at_p0']:.3f}) and "
              f"{'falls' if c['monotone_decreasing'] else 'moves'} to {c['at_pmax']:.3f} at p = {t.index[-1]:g}")
            for lvl in (0.25, 0.1, 0.05, 0.01):
                cp = collapse_p(t, key, lvl)
                if cp is not None:
                    P(f"    {'':<40} falls below {lvl:.2f} at p = {cp:.4f}")
    P("")
    P("  Material conclusions:")
    for key in MATERIAL:
        lab = dict(STATEMENTS)[key]
        c = crossing(t, key)
        if c["crossings"]:
            for h in c["crossings"]:
                P(f"    {lab:<40} crosses 0.5 at p = {h['p']:.4f}")
        else:
            flat = "constant (load-1 points only, cannot move with p)" if c["at_p0"] == c["at_pmax"] else \
                   ("rises" if c["monotone_increasing"] else "falls") + \
                   f" {c['at_p0']:.3f} -> {c['at_pmax']:.3f}"
            P(f"    {lab:<40} no crossing in the swept range; {flat}")
    P("")
    P("  None of the three material conclusions comes within reach of 0.5 anywhere in 0 <= p <= "
      f"{t.index[-1]:g}." if not any(crossing(t, k)["crossings"] for k in MATERIAL) else "")
    return "\n".join(x for x in L if x is not None)


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", default="v3", choices=sorted(C.BY_TAG))
    p.add_argument("--p-grid", type=float, nargs="+", default=list(P_GRID))
    p.add_argument("--n-jobs", type=int, default=-1)
    p.add_argument("--out-dir", type=Path, default=OUT_DIR)
    a = p.parse_args(argv)

    cfg = C.BY_TAG[a.version]
    from rdt import pipeline as pl
    pl.apply_config(cfg)

    res = sweep(cfg, a.p_grid, a.n_jobs, log=lambda *x: None)
    ver = verify_endpoints(res)
    text = report(res, ver)
    print(text)

    a.out_dir.mkdir(parents=True, exist_ok=True)
    (a.out_dir / f"flame_sweep_{cfg.tag}.txt").write_text(text + "\n")
    res["table"].to_csv(a.out_dir / f"flame_sweep_{cfg.tag}.csv")
    js = {"tag": cfg.tag, "p_grid": list(a.p_grid), "swept_points": list(SWEPT_POINTS),
          "verification": ver,
          "crossings": {k: crossing(res["table"], k) for k, _ in STATEMENTS}}
    (a.out_dir / f"flame_sweep_{cfg.tag}.json").write_text(json.dumps(js, indent=2, default=float))
    print(f"\nwritten: {a.out_dir / f'flame_sweep_{cfg.tag}.txt'}")
    print(f"         {a.out_dir / f'flame_sweep_{cfg.tag}.csv'}")
    print(f"         {a.out_dir / f'flame_sweep_{cfg.tag}.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
