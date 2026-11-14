"""Revision R1 of the paper, calculation package 1.

``python -m rdt.revision_r1 [--calc A B C D] [--n-d 512] [--n-jobs -1]``

Everything here reads the published v3 results and writes only to ``data/revision_r1/``; no ``*_v3``
file is touched. The v3 run configuration is applied first, so the creep curve, the surrogates, the
Monte Carlo design (the stored output of ``uq.sample(spec, N)`` at ``uq.py:367``) and the seeds are the
published ones.

A  surrogate_cv_full.csv        5-fold CV (mean +/- sd over folds) and Saltelli test metrics of the
                                selected model of every target (Table S4)
B  flex_magnitude_vs_p.csv      magnitude, not fraction, of the paired S2/S1 and S3/S1 life cost per
                                kmol H2 at L_q_eff = L_q * load**p
C  pareto_across_curves.csv     the 198 physics-feasible Pareto members re-scored on the other three
                                creep curves of RQ4 (+ pareto_regimes_across_curves.csv)
D  scenario_intervals.csv       Monte Carlo 5-95 % intervals of the Table 2 scenarios, each scenario
                                represented by <= 6 duration-weighted operating points; figure
                                fig11_scenarios_summary.pdf with error bars
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from rdt import config as C

OUT = C.ROOT / "data" / "revision_r1"
CFG = C.V3


def setup() -> None:
    from rdt import pipeline as pl
    pl.apply_config(CFG)


def _dump(name: str, obj) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(obj, indent=2, default=lambda v: v.item() if hasattr(v, "item") else str(v)))


def mc_design(n: Optional[int] = None) -> pd.DataFrame:
    """The stored v3 Monte Carlo design (written by ``uq.run_campaign`` from ``uq.sample(spec, N)``)."""
    from rdt import uq
    X = pd.read_csv(Path(uq.OUT_DIR) / f"{uq.MC_TAG}_design.csv.gz")
    return X if n is None else X.iloc[:n].reset_index(drop=True)


def check_design() -> float:
    """Max abs difference between the stored design and a fresh ``uq.sample(build_spec(), N)``."""
    from rdt import uq
    Xs = mc_design(); Xr = uq.sample(uq.build_spec(save=False), CFG.n_mc)
    return float(np.nanmax(np.abs(Xs[uq.PARAMS].to_numpy(float) - Xr[uq.PARAMS].to_numpy(float))))


def _ok(df: pd.DataFrame) -> np.ndarray:
    return (df["converged"] == True).to_numpy()  # noqa: E712


def _q(v: np.ndarray) -> Dict[str, float]:
    return {"median": float(np.median(v)), "p05": float(np.percentile(v, 5)), "p95": float(np.percentile(v, 95))}


# ---------------------------------------------------------------------------
# A. Surrogate cross-validation
# ---------------------------------------------------------------------------
def _metrics(y, yhat) -> Dict[str, float]:
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    return {"r2": float(r2_score(y, yhat)), "rmse": float(np.sqrt(mean_squared_error(y, yhat))),
            "mae": float(mean_absolute_error(y, yhat))}


def _fold(name: str, X: np.ndarray, y: np.ndarray, tri: np.ndarray, tei: np.ndarray, seed: int) -> Dict[str, float]:
    from rdt import surrogate as sg
    m = sg._fit(name, X[tri], y[tri], seed)          # GP: random subset of <= GP_SUBSET points of the training fold
    return _metrics(y[tei], m.predict(X[tei]))


def calc_A(n_jobs: int = 5, log=print) -> pd.DataFrame:
    import joblib
    from joblib import Parallel, delayed
    from sklearn.model_selection import KFold
    from rdt import surrogate as sg

    met = json.loads(Path(CFG.surrogate_metrics).read_text())
    tr, te = sg.load_data()
    tr = sg.add_physics_features(tr); te = sg.add_physics_features(te)
    folds = list(KFold(5, shuffle=True, random_state=sg.SEED).split(tr))
    rows = []
    # the 8 targets of Table S4, plus the hoop-stress GP the scenario/optimisation surrogate also uses
    # (scenarios.TwinSurrogate._sigma_model, gp on the 7 base inputs), which Table S4 does not list
    targets = {**{t: met["targets"][t]["best"] for t in sg.TARGETS}, "sigma_hot_MPa": "gp__base"}
    for tgt, best in targets.items():
        name, fs = best.split("__")
        feats = list(sg.INPUTS) + (sg.PHYSICS_FEATURES if fs == "physics" else [])
        X, y = tr[feats].to_numpy(float), tr[tgt].to_numpy(float)
        t0 = time.perf_counter()
        cv = Parallel(n_jobs=n_jobs)(delayed(_fold)(name, X, y, a, b, sg.SEED) for a, b in folds)
        bundle = joblib.load(Path(sg.MODELS_DIR) / f"{tgt}__{best}.joblib")
        yt = te[tgt].to_numpy(float)
        test = _metrics(yt, bundle["pipeline"].predict(te[bundle["features"]].to_numpy(float)))
        r = {"target": tgt, "unit": sg.TARGETS.get(tgt, "MPa"), "in_table_S4": tgt in sg.TARGETS, "model": best, "n_lhs": len(tr), "n_test": len(te),
             "gp_train_subset_per_fold": sg.GP_SUBSET if name == "gp" else None}
        for k in ("r2", "rmse", "mae"):
            v = np.array([c[k] for c in cv])
            r[f"cv_{k}_mean"] = float(v.mean()); r[f"cv_{k}_sd"] = float(v.std(ddof=1))
        for k in ("r2", "rmse", "mae"):
            r[f"test_{k}"] = test[k]
        pub = met["targets"].get(tgt, {}).get("models", {}).get(best, {}).get("test", {})
        r["published_test_rmse"] = pub.get("rmse", np.nan); r["published_test_r2"] = pub.get("r2", np.nan)
        rows.append(r)
        log(f"A {tgt:18s} {best:12s} cv R2 {r['cv_r2_mean']:.6f}+/-{r['cv_r2_sd']:.1e}  test RMSE {test['rmse']:.4g} "
            f"(published {r['published_test_rmse']:.4g})  {time.perf_counter() - t0:.0f}s")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "surrogate_cv_full.csv", index=False)
    _dump("calc_A.json", {"rows": rows, "folds": "KFold(5, shuffle=True, random_state=0) over all converged LHS points"})
    return df


# ---------------------------------------------------------------------------
# B. Magnitude of the load-following effect vs flame-length exponent p
# ---------------------------------------------------------------------------
P_B = (0.0, 0.1, 0.2, 0.25, 0.5, 1.0)
S3_WEIGHTS = (0.25, 0.5, 0.25)


def flex_ratios(b: pd.DataFrame, l07: pd.DataFrame, l085: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Paired S2/S1 and S3/S1 per kmol H2, exactly as ``uq.robustness_fractions`` builds S2 and S3.

    S2 = 16 h at load 1.00 + 8 h at 0.70 per day; S3 = 0.25 / 0.50 / 0.25 weights on loads 0.70 / 0.85 / 1.00;
    per kmol = sum(rate) / sum(H2); the load-1.00 state is the S1 base itself.
    """
    ok = _ok(b) & _ok(l07) & _ok(l085)
    r = lambda d, c: d.loc[ok, c].to_numpy(float)  # noqa: E731
    base = r(b, "life_rate_per_kmol_H2")
    s2 = (16 * r(b, "life_rate_per_h") + 8 * r(l07, "life_rate_per_h")) / (16 * r(b, "H2_net_kmol_h") + 8 * r(l07, "H2_net_kmol_h"))
    w = S3_WEIGHTS
    s3 = (w[0] * r(l07, "life_rate_per_h") + w[1] * r(l085, "life_rate_per_h") + w[2] * r(b, "life_rate_per_h")) / \
         (w[0] * r(l07, "H2_net_kmol_h") + w[1] * r(l085, "H2_net_kmol_h") + w[2] * r(b, "H2_net_kmol_h"))
    return s2 / base, s3 / base, ok


def calc_B(log=print) -> pd.DataFrame:
    from rdt import flame_sweep as fsw
    from rdt import uq
    mc, X, pts = fsw.stored_campaign(CFG)
    nom = uq.nominal(uq.build_spec(save=False))
    published_fsw = pd.read_csv(fsw.OUT_DIR / f"flame_sweep_{CFG.tag}.csv", index_col=0)
    b = mc["S1_base"]
    nb = uq.evaluate_one(pts["S1_base"], nom)
    rows = []
    for p in P_B:
        f = fsw.OUT_DIR / f"{uq.MC_TAG}_{{}}_{fsw._tag(p)}.csv.gz"
        l07 = pd.read_csv(str(f).format("load_0.70")); l085 = pd.read_csv(str(f).format("load_0.85"))
        assert (l07[uq.PARAMS].to_numpy() == X[uq.PARAMS].to_numpy()).all(), "flame-sweep run is not on the v3 design"
        s2, s3, ok = flex_ratios(b, l07, l085)
        n07 = uq.evaluate_one(pts["load_0.70"], nom, float(p)); n085 = uq.evaluate_one(pts["load_0.85"], nom, float(p))
        one = lambda d: pd.DataFrame([d])  # noqa: E731
        ns2, ns3, _ = flex_ratios(one(nb), one(n07), one(n085))
        for lab, v, nv, key in (("S2 daily vs S1", s2, ns2[0], "S2_life_per_kmol_le_S1"), ("S3 renewable vs S1", s3, ns3[0], "S3_life_per_kmol_le_S1")):
            frac = float(np.mean(v <= 1.0))
            ref = float(published_fsw.loc[p, key]) if p in published_fsw.index else np.nan
            rows.append({"statement": lab, "p": p, "nominal": float(nv), **_q(v), "n_ok": int(ok.sum()),
                         "fraction_le_1": frac, "fraction_le_1_flame_sweep_v3": ref,
                         "abs_change_pct_median": 100 * (float(np.median(v)) - 1)})
        log(f"B p={p:<5g} S2 nom {ns2[0]:.4f} med {np.median(s2):.4f} [{np.percentile(s2, 5):.4f},{np.percentile(s2, 95):.4f}]  "
            f"S3 nom {ns3[0]:.4f} med {np.median(s3):.4f}")
    df = pd.DataFrame(rows)
    assert np.allclose(df.fraction_le_1, df.fraction_le_1_flame_sweep_v3, atol=0, rtol=0), "fractions disagree with flame_sweep_v3"
    df.to_csv(OUT / "flex_magnitude_vs_p.csv", index=False)
    loads = {n: pts[n]["feed_per_tube_fraction"] for n in ("load_0.70", "load_0.85", "S1_base")}
    firings = {n: pts[n]["specific_firing_factor"] for n in ("load_0.70", "load_0.85", "S1_base")}
    _dump("calc_B.json", {"rows": rows, "N": int(len(X)), "loads": loads, "firings": firings,
                          "table2_hourly": {"S2": 1.010, "S3": 1.013}})
    return df


# ---------------------------------------------------------------------------
# C. Pareto set re-scored on the other creep curves
# ---------------------------------------------------------------------------
def nondominated(F: np.ndarray) -> np.ndarray:
    """Boolean mask of the non-dominated rows of F (all objectives minimised)."""
    n = len(F); nd = np.ones(n, bool)
    for i in range(n):
        le = np.all(F <= F[i], axis=1); lt = np.any(F < F[i], axis=1)
        nd[i] = not np.any(le & lt)
    return nd


TABLE3_REGIMES = ("min_life_iso_H2", "knee", "max_H2", "min_fuel_overall", "min_life_overall")


def calc_C(log=print) -> pd.DataFrame:
    from scipy.stats import spearmanr
    from rdt import alloy_comparison as ac
    from rdt import creep
    from rdt import operate as op
    from rdt import optimize as opt

    par = pd.read_csv(CFG.pareto_csv)
    feas = par[par.feasible_phys == True].copy()  # noqa: E712
    feas.index.name = "pareto_row"
    meta = json.loads(Path(CFG.pareto_meta).read_text()); bref = meta["base_reference"]
    rb = op.run_case(bref["inputs"], op.base_case())
    assert abs(rb["T_wo_max_K"] - bref["T_wo_max"]) < 1e-6 and abs(rb["H2_net_kmol_h"] - bref["H2"]) < 1e-9
    T, S, H = (feas[c].to_numpy(float) for c in ("T_wo_max_K_phys", "sigma_hot_MPa_phys", "H2_net_kmol_h_phys"))
    keys = list(ac.ALLOYS); base_key = CFG.creep_alloy
    assert keys[0] == base_key
    rel = {}
    for k in keys:
        cur = creep.load_alloy(k).minimum
        rate = 1.0 / np.asarray(cur.time_to_rupture(T, S), float) / H
        rate_b = 1.0 / float(cur.time_to_rupture(rb["T_wo_max_K"], rb["sigma_hot_MPa"])) / rb["H2_net_kmol_h"]
        rel[k] = rate / rate_b
        feas[f"rel_life_{k}"] = rel[k]
    chk = float(np.max(np.abs(rel[base_key] / feas.life_rate_per_kmol_H2_rel_base_phys.to_numpy() - 1)))
    assert chk < 1e-9, f"base curve does not reproduce the published relative life ({chk})"
    nd = {}
    for k in keys:
        F = np.column_stack([-H, feas.fuel_MJ_per_kmol_H2_phys.to_numpy(float), rel[k]])
        nd[k] = nondominated(F); feas[f"nondominated_{k}"] = nd[k]
    summary = {}
    for k in keys:
        rho = float(spearmanr(rel[base_key], rel[k]).statistic)
        summary[k] = {"label": ac.SHORT[k], "C": creep.load_alloy(k).C, "spearman_vs_base": rho,
                      "n_nondominated": int(nd[k].sum()),
                      "lost_nondominance": int(np.sum(nd[base_key] & ~nd[k])),
                      "gained_nondominance": int(np.sum(~nd[base_key] & nd[k]))}
    # Table 3 regimes: the published member, re-scored, and the member re-selected on each curve
    reg = pd.read_csv(CFG.regimes_csv, index_col=0)
    ins = list(op.INPUT_NAMES)
    def member(name):
        m = np.all(np.isclose(feas[ins].to_numpy(float), reg.loc[name, ins].to_numpy(float), rtol=0, atol=1e-12), axis=1)
        assert m.sum() == 1, name
        return feas.index[m][0]
    iso = feas[np.abs(feas.H2_net_kmol_h_phys / bref["H2"] - 1.0) <= 0.01]
    rrows = []
    for name in TABLE3_REGIMES:
        i = member(name)
        row = {"regime": name, "pareto_row": int(i), "published_rel_life": float(reg.loc[name, "life_rate_per_kmol_H2_rel_base_phys"])}
        for k in keys:
            row[f"fixed_{k}"] = float(feas.loc[i, f"rel_life_{k}"])
            if name == "min_life_iso_H2":
                j = iso[f"rel_life_{k}"].idxmin()
            elif name == "min_life_overall":
                j = feas[f"rel_life_{k}"].idxmin()
            elif name == "knee":
                j = opt.knee_point(feas, cols=("H2_net_kmol_h_phys", "fuel_MJ_per_kmol_H2_phys", f"rel_life_{k}"))
            else:                                   # max H2, min fuel: objectives do not depend on the curve
                j = i
            row[f"reselected_{k}"] = float(feas.loc[j, f"rel_life_{k}"]); row[f"reselected_row_{k}"] = int(j)
        rrows.append(row)
        log("C " + name + "  " + "  ".join(f"{ac.SHORT[k][:12]} {row[f'fixed_{k}']:.3f}/{row[f'reselected_{k}']:.3f}" for k in keys))
    rdf = pd.DataFrame(rrows).set_index("regime")
    # knee re-selection uses the same normalised-distance rule as optimize.knee_point
    keep = ins + ["H2_net_kmol_h_phys", "fuel_MJ_per_kmol_H2_phys", "T_wo_max_K_phys", "sigma_hot_MPa_phys"] + \
        [f"rel_life_{k}" for k in keys] + [f"nondominated_{k}" for k in keys]
    feas[keep].to_csv(OUT / "pareto_across_curves.csv")
    rdf.to_csv(OUT / "pareto_regimes_across_curves.csv")
    for k in keys:
        log(f"C {ac.SHORT[k]:26s} rho {summary[k]['spearman_vs_base']:.4f}  nd {summary[k]['n_nondominated']}  lost {summary[k]['lost_nondominance']}")
    _dump("calc_C.json", {"n_members": int(len(feas)), "base_reference": {"T_wo_max_K": rb["T_wo_max_K"], "sigma_hot_MPa": rb["sigma_hot_MPa"],
                                                                          "H2": rb["H2_net_kmol_h"], "firing": bref["inputs"]["specific_firing_factor"]},
                          "base_curve_reproduction_max_rel_err": chk, "curves": summary, "regimes": rdf.reset_index().to_dict(orient="records")})
    return feas


# ---------------------------------------------------------------------------
# D. Monte Carlo intervals of the RQ2 scenarios
# ---------------------------------------------------------------------------
TABLE2 = ("S1_steady", "S2_daily", "S3_renewable", "S4a_mild_overfire", "S4b_severe_overfire", "S4c_hot_band", "S6_SC2.5", "S6_SC3.5",
          "S5_hold_CH4_slip_ageing_y1", "S5_hold_CH4_slip_ageing_y4", "S5_hold_T_out_ageing_y1", "S5_hold_T_out_ageing_y4")
MAX_BINS = 6
BIN_TOL_SURR = 1e-3     # binning error allowed on the surrogate (relative), tighter than the 1 % physics criterion
PHYS_TOL = 1e-2


def hourly_with_add(name: str) -> pd.DataFrame:
    """Hourly v3 history with the additive hot-band perturbation of T_wo,max reconstructed (S4c only)."""
    from rdt import scenarios as sc
    from rdt import operate as op
    h = pd.read_csv(CFG.scenario_hourly(name))
    add = np.zeros(len(h))
    if name == "S4c_hot_band":
        add = np.asarray(sc.make_scenarios()["S4c_hot_band"].T_wo_add_K, float)
        s1 = pd.read_csv(CFG.scenario_hourly("S1_steady")).T_wo_max_K.to_numpy()
        assert np.allclose(h.T_wo_max_K.to_numpy() - add, s1, atol=1e-9), "hot band reconstruction failed"
    h["T_wo_add_K"] = add
    return h


def make_bins(h: pd.DataFrame, k: int) -> pd.DataFrame:
    """At most k duration-weighted representative points: 1-D k-means on log10 of the hourly damage rate,
    each bin represented by the hour-weighted mean of its inputs (and of the hot-band offset)."""
    from sklearn.cluster import KMeans
    from rdt import operate as op
    ld = np.log10(h.dD.to_numpy(float))
    nu = len(np.unique(np.round(ld, 9)))          # float noise (~1e-17 in dD) is not a distinct state
    k = min(k, nu)
    lab = np.zeros(len(h), int) if k == 1 else KMeans(k, n_init=10, random_state=0).fit(ld.reshape(-1, 1)).labels_
    cols = list(op.INPUT_NAMES) + ["T_wo_add_K"]
    g = h.assign(_b=lab).groupby("_b")
    out = g[cols].mean(); out["hours"] = g.size()
    return out.reset_index(drop=True)


def bins_value_surrogate(bins: pd.DataFrame, tw) -> Tuple[float, float]:
    p = tw.predict(bins[list(tw_inputs())], physics_fallback=True)
    t_r = tw.curve.time_to_rupture(p.T_wo_max_K.to_numpy(float) + bins.T_wo_add_K.to_numpy(float), p.sigma_hot_MPa.to_numpy(float))
    n = bins.hours.to_numpy(float)
    return float(np.sum(n / np.asarray(t_r, float))), float(np.sum(n * p.H2_net_kmol_h.to_numpy(float)))


def tw_inputs():
    from rdt import operate as op
    return op.INPUT_NAMES


def _theta_with_add(theta: Dict[str, float], add: float) -> Dict[str, float]:
    t = dict(theta); t["dT_wo_meas_K"] = float(t["dT_wo_meas_K"]) + float(add); return t


def bins_value_physics(bins: pd.DataFrame, theta: Dict[str, float]) -> Tuple[float, float, bool]:
    from rdt import uq
    D = H = 0.0
    for _, b in bins.iterrows():
        r = uq.evaluate_one({k: float(b[k]) for k in tw_inputs()}, _theta_with_add(theta, b.T_wo_add_K))
        if not r["converged"]:
            return np.nan, np.nan, False
        D += b.hours * r["life_rate_per_h"]; H += b.hours * r["H2_net_kmol_h"]
    return D, H, True


def _key(b: pd.Series) -> str:
    s = json.dumps({k: round(float(b[k]), 10) for k in list(tw_inputs()) + ["T_wo_add_K"]}, sort_keys=True)
    return hashlib.sha1(s.encode()).hexdigest()[:12]


def _mc_worker(point, theta, add):
    from rdt import uq
    r = uq.evaluate_one(point, _theta_with_add(theta, add))
    return {k: r.get(k, np.nan) for k in ("converged", "life_rate_per_h", "H2_net_kmol_h", "T_wo_max_K", "life_rate_per_kmol_H2")}


def run_bin_mc(bins_all: Dict[str, pd.DataFrame], X: pd.DataFrame, n_jobs: int = -1, log=print) -> Dict[str, pd.DataFrame]:
    from joblib import Parallel, delayed
    d = OUT / "scenario_mc"; d.mkdir(parents=True, exist_ok=True)
    uniq = {}
    for bins in bins_all.values():
        for _, b in bins.iterrows():
            uniq.setdefault(_key(b), b)
    todo = [k for k in uniq if not (d / f"{k}.csv.gz").exists()]
    log(f"D unique representative points {len(uniq)}, to run {len(todo)} x {len(X)} samples")
    th = X.to_dict(orient="records")
    jobs = [(k, i) for k in todo for i in range(len(th))]
    if jobs:
        t0 = time.perf_counter()
        res = Parallel(n_jobs=n_jobs, batch_size=8)(
            delayed(_mc_worker)({c: float(uniq[k][c]) for c in tw_inputs()}, th[i], float(uniq[k].T_wo_add_K)) for k, i in jobs)
        log(f"D MC done in {time.perf_counter() - t0:.0f}s")
        for j, k in enumerate(todo):
            df = pd.DataFrame(res[j * len(th):(j + 1) * len(th)]); df["sample"] = np.arange(len(th))
            for c in list(tw_inputs()) + ["T_wo_add_K"]:
                df[c] = float(uniq[k][c])
            df.to_csv(d / f"{k}.csv.gz", index=False, compression="gzip")
    return {k: pd.read_csv(d / f"{k}.csv.gz") for k in uniq}


def scenario_ratio(bins: pd.DataFrame, s1bins: pd.DataFrame, runs: Dict[str, pd.DataFrame]) -> Tuple[np.ndarray, np.ndarray]:
    def acc(bb):
        D = H = 0.0; ok = True
        for _, b in bb.iterrows():
            r = runs[_key(b)]; ok = ok & _ok(r)
            D = D + b.hours * r.life_rate_per_h.to_numpy(float); H = H + b.hours * r.H2_net_kmol_h.to_numpy(float)
        return D, H, ok
    D, H, ok = acc(bins); D1, H1, ok1 = acc(s1bins); ok = ok & ok1
    return ((D / H) / (D1 / H1))[ok], ok


LABELS = {"S1_steady": "S1 steady", "S2_daily": "S2 daily", "S3_renewable": "S3 renewable", "S4a_mild_overfire": "S4a mild", "S4b_severe_overfire": "S4b severe",
          "S4c_hot_band": "S4c hot band", "S6_SC2.5": "S6 S/C 2.5", "S6_SC3.5": "S6 S/C 3.5", "S5_hold_CH4_slip_ageing_y1": "S5 slip y1",
          "S5_hold_CH4_slip_ageing_y4": "S5 slip y4", "S5_hold_T_out_ageing_y1": "S5 $T_{out}$ y1", "S5_hold_T_out_ageing_y4": "S5 $T_{out}$ y4"}


def calc_D(n: int = 512, n_jobs: int = -1, log=print) -> pd.DataFrame:
    from rdt import scenarios as sc
    from rdt import uq
    tw = sc.TwinSurrogate()
    nom = uq.nominal(uq.build_spec(save=False))
    summ = pd.read_csv(CFG.scenarios_summary).set_index("scenario")
    H = {s: hourly_with_add(s) for s in TABLE2}
    # --- choose the representation ------------------------------------------------------------
    chosen, info = {}, {}
    s1 = make_bins(H["S1_steady"], 1)
    Ds1, Hs1 = bins_value_surrogate(s1, tw); Dp1, Hp1, _ = bins_value_physics(s1, nom)
    for s in TABLE2:
        tab = float(summ.loc[s, "life_consumption_per_kmol_H2_rel_S1"]); hist = []
        for k in range(1, MAX_BINS + 1):
            bins = make_bins(H[s], k)
            Ds, Hs = bins_value_surrogate(bins, tw); rs = (Ds / Hs) / (Ds1 / Hs1)
            hist.append({"k": k, "n_bins": len(bins), "rel_surrogate": rs, "err_surrogate": rs / tab - 1})
            if abs(rs / tab - 1) <= BIN_TOL_SURR or len(bins) < k or k == MAX_BINS:
                Dp, Hp, okp = bins_value_physics(bins, nom); rp = (Dp / Hp) / (Dp1 / Hp1)
                hist[-1].update({"rel_physics": rp, "err_physics": rp / tab - 1})
                if abs(rp / tab - 1) <= PHYS_TOL or len(bins) < k or k == MAX_BINS:
                    break
        chosen[s] = bins
        info[s] = {"table2": tab, "history": hist, **hist[-1]}
        log(f"D {s:28s} bins {len(bins)}  table {tab:.4f}  surr {hist[-1]['rel_surrogate']:.4f} ({100 * hist[-1]['err_surrogate']:+.3f} %)"
            f"  phys {hist[-1]['rel_physics']:.4f} ({100 * hist[-1]['err_physics']:+.3f} %)")
    pd.concat({s: b for s, b in chosen.items()}, names=["scenario", "bin"]).to_csv(OUT / "scenario_bins.csv")
    # --- Monte Carlo ---------------------------------------------------------------------------
    X = mc_design(n)
    runs = run_bin_mc(chosen, X, n_jobs, log)
    # sanity: the S1 point must reproduce the stored v3 campaign sample for sample
    st = pd.read_csv(Path(uq.OUT_DIR) / f"{uq.MC_TAG}_S1_base.csv.gz").iloc[:n]
    r1 = runs[_key(chosen["S1_steady"].iloc[0])]; okb = _ok(st) & _ok(r1)
    s1_dev = float(np.max(np.abs(r1.life_rate_per_kmol_H2.to_numpy()[okb] / st.life_rate_per_kmol_H2.to_numpy()[okb] - 1)))
    rows = []
    for s in TABLE2:
        v, ok = scenario_ratio(chosen[s], chosen["S1_steady"], runs)
        i = info[s]
        rows.append({"scenario": s, "label": LABELS[s], "table2_rel": i["table2"], "n_bins": i["n_bins"],
                     "hours_per_bin": ";".join(str(int(x)) for x in chosen[s].hours),
                     "rep_rel_surrogate_nominal": i["rel_surrogate"], "rep_err_surrogate_pct": 100 * i["err_surrogate"],
                     "rep_rel_physics_nominal": i["rel_physics"], "rep_err_physics_pct": 100 * i["err_physics"],
                     "mc_n_ok": int(ok.sum()), **{f"mc_{k}": x for k, x in _q(v).items()},
                     "mc_frac_gt_1": float(np.mean(v > 1.0))})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "scenario_intervals.csv", index=False)
    _dump("calc_D.json", {"N": n, "rows": rows, "representation": info, "S1_vs_stored_mc_v3_max_rel_dev": s1_dev,
                          "s1_firing": float(chosen["S1_steady"].specific_firing_factor.iloc[0])})
    log(f"D S1 point vs stored mc_v3_S1_base: max rel dev {s1_dev:.2e}")
    fig_D(df, summ)
    return df


def fig_D(df: pd.DataFrame, summ: pd.DataFrame) -> None:
    """fig11_scenarios_summary in the style of rdt.figures.fig11_scenario_summary, with MC 5-95 % bars."""
    from rdt import plotting as P
    order = list(df.scenario)
    D1, H1 = summ.loc["S1_steady", "annual_damage_D"], summ.loc["S1_steady", "annual_H2_kmol_per_tube"]
    rel = summ.loc[order, "life_consumption_per_kmol_H2_rel_S1"].to_numpy()
    rel_avg = ((summ.loc[order, "D_avg_condition"] / summ.loc[order, "annual_H2_kmol_per_tube"]) / (D1 / H1)).to_numpy()
    lo, hi, med = df.mc_p05.to_numpy(), df.mc_p95.to_numpy(), df.mc_median.to_numpy()
    with P.style():
        fig, ax = P.figure(P.DOUBLE, 6 * P.CM); x = np.arange(len(order))
        ax.bar(x, rel, color=P.PALETTE[0], label="hourly accumulation"); ax.plot(x, rel_avg, "D", color="k", ms=3, label="annual-mean inputs")
        m = np.array([s != "S1_steady" for s in order])
        ax.errorbar(x[m] + 0.22, med[m], yerr=[med[m] - lo[m], hi[m] - med[m]], fmt="o", ms=2.5, color=P.PALETTE[1], ecolor=P.PALETTE[1],
                    elinewidth=0.8, capsize=2, label="Monte Carlo median, 5–95 %")
        ax.axhline(1.0, color="grey", lw=0.5)
        ax.set_xticks(x); ax.set_xticklabels([LABELS[o] for o in order], rotation=40, ha="right"); ax.set_yscale("log")
        P.tidy(ax, None, "Life consumption per kmol H$_2$, rel. S1 (–)", True)
        fig.tight_layout(); P.save(fig, OUT / "fig11_scenarios_summary")


# ---------------------------------------------------------------------------
# Report tables
# ---------------------------------------------------------------------------
def _md(df: pd.DataFrame, fmt: Dict[str, str]) -> str:
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(format(r[c], fmt.get(c, "")) if not isinstance(r[c], str) else r[c] for c in cols) + " |")
    return "\n".join(out)


def report_tables() -> Dict[str, str]:
    from rdt import alloy_comparison as ac
    t: Dict[str, str] = {}
    a = pd.read_csv(OUT / "surrogate_cv_full.csv")
    pm = lambda m, s, f: f"{m:{f}} ± {s:.1e}"  # noqa: E731
    ta = pd.DataFrame({"target": [x + ("" if s4 else " *") for x, s4 in zip(a.target, a.in_table_S4)], "model": a.model,
                       "CV R²": [pm(m, s, ".6f") for m, s in zip(a.cv_r2_mean, a.cv_r2_sd)],
                       "CV RMSE": [f"{m:.4g} ± {s:.2g}" for m, s in zip(a.cv_rmse_mean, a.cv_rmse_sd)],
                       "CV MAE": [f"{m:.4g} ± {s:.2g}" for m, s in zip(a.cv_mae_mean, a.cv_mae_sd)],
                       "test R²": [f"{v:.6f}" for v in a.test_r2], "test RMSE": [f"{v:.4g}" for v in a.test_rmse],
                       "test MAE": [f"{v:.4g}" for v in a.test_mae], "unit": a.unit})
    t["A"] = _md(ta, {})
    b = pd.read_csv(OUT / "flex_magnitude_vs_p.csv")
    tb = b[["statement", "p", "nominal", "median", "p05", "p95", "fraction_le_1"]].rename(columns={"p05": "5 %", "p95": "95 %", "fraction_le_1": "доля ≤ 1 (Table 5)"})
    t["B"] = _md(tb, {"p": "g", "nominal": ".4f", "median": ".4f", "5 %": ".4f", "95 %": ".4f", "доля ≤ 1 (Table 5)": ".3f"})
    cj = json.loads((OUT / "calc_C.json").read_text())
    tc1 = pd.DataFrame([{"кривая": v["label"], "C": v["C"], "Spearman ρ vs G 4852": v["spearman_vs_base"],
                         "недоминируемых из 198": v["n_nondominated"], "потеряли статус": v["lost_nondominance"],
                         "приобрели статус": v["gained_nondominance"]} for v in cj["curves"].values()])
    t["C1"] = _md(tc1, {"C": ".1f", "Spearman ρ vs G 4852": ".4f", "недоминируемых из 198": "d", "потеряли статус": "d", "приобрели статус": "d"})
    reg = pd.DataFrame(cj["regimes"]).set_index("regime")
    rows = []
    for name in TABLE3_REGIMES:
        r = {"режим (Table 3)": name}
        for k in ac.ALLOYS:
            same = reg.loc[name, f"reselected_row_{k}"] == reg.loc[name, "pareto_row"]
            r[ac.SHORT[k]] = f"{reg.loc[name, f'fixed_{k}']:.3f}" + ("" if same else f" (пересел. {reg.loc[name, f'reselected_{k}']:.3f})")
        rows.append(r)
    t["C2"] = _md(pd.DataFrame(rows), {})
    d = pd.read_csv(OUT / "scenario_intervals.csv")
    td = d[["label", "table2_rel", "n_bins", "hours_per_bin", "rep_rel_physics_nominal", "rep_err_physics_pct", "mc_median", "mc_p05", "mc_p95", "mc_frac_gt_1", "mc_n_ok"]]
    td = td.rename(columns={"label": "сценарий", "table2_rel": "Table 2", "n_bins": "бинов", "hours_per_bin": "часы по бинам",
                            "rep_rel_physics_nominal": "номинал (физика, бины)", "rep_err_physics_pct": "откл. от Table 2, %",
                            "mc_median": "MC медиана", "mc_p05": "5 %", "mc_p95": "95 %", "mc_frac_gt_1": "доля > 1", "mc_n_ok": "N ok"})
    t["D"] = _md(td, {"Table 2": ".3f", "бинов": "d", "номинал (физика, бины)": ".3f", "откл. от Table 2, %": "+.2f", "MC медиана": ".3f",
                      "5 %": ".3f", "95 %": ".3f", "доля > 1": ".3f", "N ok": "d"})
    return t


# ===========================================================================
# Package 2
# ===========================================================================
# ---------------------------------------------------------------------------
# E. Catalyst-ageing histories (S5)
# ---------------------------------------------------------------------------
HOURS = 8760
E_HISTORIES = {"E0": "0.30 -> 0.10 linear over 4 y (paper base)", "E1": "0.30 -> 0.15 linear over 4 y",
               "E2": "0.25 -> 0.10 linear over 4 y",
               "E3": "0.30 -> 0.18 geometric over the first 6 months, then linear to 0.10 at the end of year 4",
               "E4": "0.30 -> 0.10 linear over 3 y", "E5": "0.30 -> 0.10 linear over 5 y"}
CONTROLS = ("hold_CH4_slip", "hold_T_out")
A_CAL = 0.20            # calibrated activity (Latham f_prx)


def activity_profile(name: str) -> np.ndarray:
    """Hourly activity over the whole campaign; linear profiles use the same t grid as scenarios.ageing_campaign."""
    lin = {"E0": (0.30, 0.10, 4), "E1": (0.30, 0.15, 4), "E2": (0.25, 0.10, 4), "E4": (0.30, 0.10, 3), "E5": (0.30, 0.10, 5)}
    if name in lin:
        a0, a1, y = lin[name]
        t = np.arange(y * HOURS) / (y * HOURS)
        return a0 + (a1 - a0) * t
    if name == "E3":
        ty = np.arange(4 * HOURS) / HOURS
        return np.where(ty < 0.5, 0.30 * (0.18 / 0.30) ** (ty / 0.5), 0.18 + (0.10 - 0.18) * (ty - 0.5) / 3.5)
    raise KeyError(name)


def year_histories(name: str, control: str):
    from rdt import scenarios as sc
    act = activity_profile(name); n_y = len(act) // HOURS; out = []
    for y in range(n_y):
        d = sc.base_inputs(); a = act[y * HOURS:(y + 1) * HOURS]
        out.append(sc.History(name=f"{name}_{control}_y{y + 1}", load=d["load"], steam_to_carbon=d["steam_to_carbon"], inlet_T=d["inlet_T"],
                              inlet_P=d["inlet_P"], excess_air=d["excess_air"], activity=a, control=control))
    return out


def _phys_point(x: Dict[str, float], eta_mult: float) -> Dict[str, float]:
    """Physics model at the nominal (calibrated) parameters, effectiveness factors of both zones x eta_mult."""
    from rdt import uq
    th = dict(_NOM); th["log_eta_mult"] = float(np.log10(eta_mult))
    return uq.evaluate_one(x, th)


_NOM: Dict[str, float] = {}


def solve_firing_phys(x: Dict[str, float], eta_mult: float, mode: str, target: float, lo: float = 0.85, hi: float = 1.15) -> Dict[str, float]:
    """Physics-model firing that holds T_out or the dry CH4 slip (same bracket as TwinSurrogate.solve_firing)."""
    from scipy.optimize import brentq
    key = "T_out_K" if mode == "hold_T_out" else "CH4_slip_dry_pct"
    cache: Dict[float, Dict[str, float]] = {}
    def run(f):
        if f not in cache:
            cache[f] = _phys_point({**x, "specific_firing_factor": float(f)}, eta_mult)
        return cache[f]
    g = lambda f: run(f)[key] - target  # noqa: E731
    glo, ghi = g(lo), g(hi); sat = None
    if glo * ghi > 0:                                   # target not reachable inside the bracket: clamp
        f = lo if abs(glo) < abs(ghi) else hi; sat = "lo" if f == lo else "hi"
    else:
        f = brentq(g, lo, hi, xtol=1e-6)
    r = dict(run(f)); r["firing"] = float(f); r["firing_saturated"] = sat
    return r


def _eta_bin_task(x, a, eta_on, mode, target):
    _ensure_setup()
    m = (a / A_CAL) ** -0.5 if eta_on else 1.0
    r = solve_firing_phys({**x, "catalyst_activity": float(a)}, m, mode, target)
    r.update({"activity": float(a), "eta_mult": float(m)}); return r


_SETUP_DONE = False


def _ensure_setup() -> None:
    """Workers (loky) start fresh interpreters: re-apply the v3 configuration and the nominal theta."""
    global _SETUP_DONE, _NOM
    if not _SETUP_DONE:
        from rdt import uq
        setup(); _NOM.clear(); _NOM.update(uq.nominal(uq.build_spec(save=False))); _SETUP_DONE = True


def calc_E(n_bins: int = 6, log=print) -> pd.DataFrame:
    from joblib import Parallel, delayed
    from rdt import scenarios as sc
    _ensure_setup()
    tw = sc.TwinSurrogate()
    summ = pd.read_csv(CFG.scenarios_summary).set_index("scenario")
    D1, H1 = float(summ.loc["S1_steady", "annual_damage_D"]), float(summ.loc["S1_steady", "annual_H2_kmol_per_tube"])
    slip0 = sc.steady_state_closed_form(tw)["CH4_slip_dry_pct"]
    rows = []
    # --- (i) surrogate hourly histories, exactly the v3 scenario path --------------------------
    for name in E_HISTORIES:
        for control in CONTROLS:
            for h in year_histories(name, control):
                r = sc.evaluate(h, tw, slip_target=slip0); hh = r["hourly"]
                f = hh.specific_firing_factor.to_numpy()
                y = int(h.name.rsplit("_y", 1)[1])
                rows.append({"history": name, "variant": "surrogate_hourly", "control": control, "year": y,
                             "activity_start": float(h.activity[0]), "activity_end": float(h.activity[-1]),
                             "rel_life_per_kmol_H2_vs_S1": (r["annual_damage"] / r["annual_H2_kmol"]) / (D1 / H1),
                             "T_wo_max_mean_K": r["T_wo_max_mean_K"], "T_wo_max_max_K": r["T_wo_max_max_K"],
                             "firing_start": float(f[0]), "firing_end": float(f[-1]),
                             "firing_at_bound_h": int(np.sum((f < 0.85 + 1e-4) | (f > 1.15 - 1e-4))),
                             "annual_D": r["annual_damage"], "annual_H2_kmol": r["annual_H2_kmol"],
                             "CH4_slip_mean_pct": r["CH4_slip_mean_pct"], "T_out_mean_K": r["T_out_mean_K"],
                             "verify_max_abs_dT_wo_K": r["verify_max_abs_dT_wo_K"]})
                log(f"E {name} {control:14s} y{y}  rel {rows[-1]['rel_life_per_kmol_H2_vs_S1']:.4f}  T_wo {r['T_wo_max_mean_K']:.1f}/{r['T_wo_max_max_K']:.1f}")
    # --- (ii) physics on a binned E0 history, with and without the eta(a) correction -----------
    b = sc.base_inputs(1); x0 = {"steam_to_carbon": float(b["steam_to_carbon"][0]), "feed_per_tube_fraction": 1.0, "inlet_T": float(b["inlet_T"][0]),
                                 "inlet_P": float(b["inlet_P"][0]), "catalyst_activity": A_CAL, "excess_air": float(b["excess_air"][0])}
    s1 = solve_firing_phys(x0, 1.0, "hold_T_out", sc.T_OUT_TARGET)
    slip_phys = float(s1["CH4_slip_dry_pct"])
    target = {"hold_T_out": sc.T_OUT_TARGET, "hold_CH4_slip": slip_phys}
    act = activity_profile("E0"); n_y = len(act) // HOURS
    edges = np.linspace(0, HOURS, n_bins + 1).astype(int)
    tasks = []
    for y in range(n_y):
        for j in range(n_bins):
            a = float(act[y * HOURS + edges[j]:y * HOURS + edges[j + 1]].mean())
            for control in CONTROLS:
                for eta_on in (False, True):
                    tasks.append((y + 1, j, int(edges[j + 1] - edges[j]), control, eta_on, a))
    t0 = time.perf_counter()
    res = Parallel(n_jobs=-1)(delayed(_eta_bin_task)(x0, a, eta_on, c, target[c]) for (_, _, _, c, eta_on, a) in tasks)
    log(f"E physics bins: {len(tasks)} firing solves in {time.perf_counter() - t0:.0f}s")
    bt = pd.DataFrame([{"year": t[0], "bin": t[1], "hours": t[2], "control": t[3], "eta_correction": t[4], **r} for t, r in zip(tasks, res)])
    bt.to_csv(OUT / "ageing_eta_bins.csv", index=False)
    d1p = s1["life_rate_per_kmol_H2"]
    for (control, eta_on, y), g in bt.groupby(["control", "eta_correction", "year"]):
        n = g.hours.to_numpy(float)
        rel = (np.sum(n * g.life_rate_per_h) / np.sum(n * g.H2_net_kmol_h)) / d1p
        rows.append({"history": "E0-eta" if eta_on else "E0", "variant": "physics_binned", "control": control, "year": int(y),
                     "activity_start": float(act[(y - 1) * HOURS]), "activity_end": float(act[y * HOURS - 1]),
                     "rel_life_per_kmol_H2_vs_S1": float(rel), "T_wo_max_mean_K": float(np.average(g.T_wo_max_K, weights=n)),
                     "T_wo_max_max_K": float(g.T_wo_max_K.max()), "firing_start": float(g.firing.iloc[0]), "firing_end": float(g.firing.iloc[-1]),
                     "firing_at_bound_h": int(n[g.firing_saturated.notna().to_numpy()].sum()),
                     "annual_H2_kmol": float(np.sum(n * g.H2_net_kmol_h)), "CH4_slip_mean_pct": float(np.average(g.CH4_slip_dry_pct, weights=n)),
                     "T_out_mean_K": float(np.average(g.T_out_K, weights=n)), "eta_mult_range": f"{g.eta_mult.min():.3f}-{g.eta_mult.max():.3f}"})
    df = pd.DataFrame(rows)
    # check E0 against Table 2 / summary_v3
    chk = {}
    for control in CONTROLS:
        for y in range(1, 5):
            ref = float(summ.loc[f"S5_{control}_ageing_y{y}", "life_consumption_per_kmol_H2_rel_S1"])
            s = df[(df.history == "E0") & (df.variant == "surrogate_hourly") & (df.control == control) & (df.year == y)].rel_life_per_kmol_H2_vs_S1.iloc[0]
            p = df[(df.history == "E0") & (df.variant == "physics_binned") & (df.control == control) & (df.year == y)].rel_life_per_kmol_H2_vs_S1.iloc[0]
            chk[f"{control}_y{y}"] = {"summary_v3": ref, "surrogate_hourly": float(s), "physics_binned": float(p),
                                      "surrogate_rel_err": float(s / ref - 1), "physics_rel_err": float(p / ref - 1)}
    df.to_csv(OUT / "ageing_histories.csv", index=False)
    _dump("calc_E.json", {"histories": E_HISTORIES, "slip_target_surrogate": slip0, "slip_target_physics": slip_phys,
                          "S1_physics": {k: s1[k] for k in ("firing", "T_wo_max_K", "CH4_slip_dry_pct", "T_out_K")}, "n_bins_per_year": n_bins,
                          "check_E0_vs_table2": chk})
    return df


# ---------------------------------------------------------------------------
# F. NSGA-II with narrowed decision bounds
# ---------------------------------------------------------------------------
F_VARIANTS = {"F0": {}, "F1": {"steam_to_carbon": (None, 3.5)}, "F2": {"steam_to_carbon": (None, 3.0)},
              "F3": {"feed_per_tube_fraction": (None, 1.0)}, "F4": {"excess_air": (10.0, None)},
              "F5": {"steam_to_carbon": (None, 3.5), "feed_per_tube_fraction": (None, 1.0), "excess_air": (10.0, None)},
              "F6": {"steam_to_carbon": (None, 3.0), "feed_per_tube_fraction": (None, 1.0), "excess_air": (10.0, None)}}


def calc_F(log=print) -> pd.DataFrame:
    from rdt import optimize as opt
    from rdt import scenarios as sc
    meta = json.loads(Path(CFG.pareto_meta).read_text()); st = meta["settings"]
    tw = sc.TwinSurrogate(); ref = opt.base_reference(tw)
    rows, info = [], {}
    for v, lim in F_VARIANTS.items():
        prob = opt.ReformerProblem(tw, ref)
        xl, xu = prob.xl.copy(), prob.xu.copy()
        for k, (lo, hi) in lim.items():
            i = opt.DECISION.index(k)
            if lo is not None: xl[i] = lo
            if hi is not None: xu[i] = hi
        prob.xl, prob.xu = xl, xu
        res, dt = opt.run_nsga2(prob, st["pop"], st["gens"], st["seed"])
        pareto = prob.evaluate_table(res.X)
        ver = opt.physics_verify(pareto, tw, ref)
        ver.to_csv(OUT / f"pareto_{v}.csv", index=False)
        feas = ver[ver.feasible_phys].copy()
        iso = opt.iso_production(feas, ref); dom = opt.base_domination(feas, ref); reg = opt.regimes_table(feas, ref, iso)
        info[v] = {"bounds": {k: [float(a), float(b)] for k, a, b in zip(opt.DECISION, xl, xu)}, "n_pareto": int(len(pareto)),
                   "n_feasible_phys": int(len(feas)), "n_dominating_base": dom["n_dominating"], "n_iso_band": int(iso.get("n_in_band", 0)),
                   "wall_time_s": dt, "max_disc_life": float(np.nanmax(np.abs(ver.disc_life_rate_per_kmol_H2_rel_base)))}
        if v == "F0":
            p3 = pd.read_csv(CFG.pareto_csv)
            a = np.sort(ver[list(opt.DECISION)].to_numpy(), axis=0); b = np.sort(p3[list(opt.DECISION)].to_numpy(), axis=0)
            info[v]["max_abs_diff_vs_pareto_v3"] = float(np.max(np.abs(a - b))) if a.shape == b.shape else float("nan")
        for rg in ("knee", "min_life_iso_H2"):
            r = {"variant": v, "regime": rg, "constraints": "; ".join(f"{k} {'<=' if hi is not None else '>='} {hi if hi is not None else lo:g}"
                                                                    for k, (lo, hi) in lim.items()) or "v3 bounds",
                 "n_pareto": len(pareto), "n_feasible_phys": len(feas), "n_dominating_base": dom["n_dominating"], "n_iso_band": info[v]["n_iso_band"]}
            if rg in reg.index:
                q = reg.loc[rg]
                r.update({k: float(q[k]) for k in opt.DECISION})
                r.update({"H2_kmol_h": float(q.H2_net_kmol_h_phys), "fuel_MJ_per_kmol_H2": float(q.fuel_MJ_per_kmol_H2_phys),
                          "rel_life_cost": float(q.life_rate_per_kmol_H2_rel_base_phys), "T_wo_max_K": float(q.T_wo_max_K_phys),
                          "CH4_slip_pct": float(q.CH4_slip_dry_pct_phys), "H2_rel_base_pct": float(q.H2_rel_base_pct), "fuel_rel_base_pct": float(q.fuel_rel_base_pct)})
            rows.append(r)
        log(f"F {v}: pareto {len(pareto)} feasible {len(feas)} dominating {dom['n_dominating']} iso {info[v]['n_iso_band']}  "
            + "  ".join(f"{rr['regime']} life {rr.get('rel_life_cost', np.nan):.3f} H2 {rr.get('H2_rel_base_pct', np.nan):+.1f}% fuel {rr.get('fuel_rel_base_pct', np.nan):+.1f}%" for rr in rows[-2:])
            + f"  ({dt:.0f}s)")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "pareto_constrained.csv", index=False)
    _dump("calc_F.json", {"settings": st, "base_reference": {"H2": ref.H2, "fuel": ref.fuel, "slip": ref.slip, "T_wo_max": ref.T_wo_max}, "variants": info})
    return df


# ---------------------------------------------------------------------------
# G. Effective-stress multiplier
# ---------------------------------------------------------------------------
K_STRESS = (1.0, 1.5, 2.0, 3.0)


def calc_G(log=print) -> pd.DataFrame:
    from rdt import creep
    from rdt import operate as op
    cur = creep.load_alloy(CFG.creep_alloy).minimum
    meta = json.loads(Path(CFG.pareto_meta).read_text()); bref = meta["base_reference"]
    rb = op.run_case(bref["inputs"], op.base_case())
    par = pd.read_csv(CFG.pareto_csv); feas = par[par.feasible_phys == True]  # noqa: E712
    reg = pd.read_csv(CFG.regimes_csv, index_col=0); ins = list(op.INPUT_NAMES)
    states = {}
    for name in TABLE3_REGIMES:
        m = np.all(np.isclose(feas[ins].to_numpy(float), reg.loc[name, ins].to_numpy(float), rtol=0, atol=1e-12), axis=1)
        q = feas[m].iloc[0]; states[name] = (float(q.T_wo_max_K_phys), float(q.sigma_hot_MPa_phys), float(q.H2_net_kmol_h_phys), float(q.life_rate_per_kmol_H2_rel_base_phys))
    scen = {s: pd.read_csv(CFG.scenario_hourly(s), usecols=["T_wo_max_K", "sigma_hot_MPa", "H2_net_kmol_h"])
            for s in ("S1_steady", "S5_hold_CH4_slip_ageing_y4", "S6_SC2.5", "S6_SC3.5")}
    summ = pd.read_csv(CFG.scenarios_summary).set_index("scenario")
    lmp_lo, lmp_hi = cur.lmp_range
    s_lo, s_hi = float(10 ** cur.log10_stress(lmp_hi)), float(10 ** cur.log10_stress(lmp_lo))
    rows = []
    for k in K_STRESS:
        sig_all = np.concatenate([[rb["sigma_hot_MPa"]], [s[1] for s in states.values()]] + [d.sigma_hot_MPa.to_numpy() for d in scen.values()])
        L_all = np.vectorize(cur.lmp_of_stress)(np.unique(np.round(k * sig_all, 6)))
        inside = bool(np.all((L_all >= lmp_lo) & (L_all <= lmp_hi)))
        rate = lambda T, s, H: 1.0 / np.asarray(cur.time_to_rupture(T, k * np.asarray(s)), float) / H  # noqa: E731
        base = float(rate(rb["T_wo_max_K"], rb["sigma_hot_MPa"], rb["H2_net_kmol_h"]))
        def slope(T, s):
            L = cur.lmp_of_stress(k * s); return -L / (cur.scale * T ** 2)
        r = {"k": k, "sigma_eff_base_MPa": k * rb["sigma_hot_MPa"], "LMP_base": cur.lmp_of_stress(k * rb["sigma_hot_MPa"]),
             "t_r_base_y": float(cur.time_to_rupture(rb["T_wo_max_K"], k * rb["sigma_hot_MPa"])) / 8760,
             "slope_dlog10tr_dT_base_per_K": slope(rb["T_wo_max_K"], rb["sigma_hot_MPa"]),
             "life_halving_dT_base_K": float(np.log10(2) / abs(slope(rb["T_wo_max_K"], rb["sigma_hot_MPa"]))),
             "sigma_eff_range_MPa": f"{k * sig_all.min():.2f}-{k * sig_all.max():.2f}", "curve_stress_range_MPa": f"{s_lo:.2f}-{s_hi:.2f}",
             "LMP_min": float(L_all.min()), "LMP_max": float(L_all.max()), "inside_digitised_range": inside}
        for name, (T, s, H, pub) in states.items():
            r[f"{name}_vs_base"] = float(rate(T, s, H)) / base
            if k == 1.0:
                assert abs(r[f"{name}_vs_base"] / pub - 1) < 1e-9
        h1 = scen["S1_steady"]; s1 = np.sum(rate(h1.T_wo_max_K.to_numpy(), h1.sigma_hot_MPa.to_numpy(), 1.0)) / h1.H2_net_kmol_h.sum()
        T1, sg1 = float(h1.T_wo_max_K.iloc[0]), float(h1.sigma_hot_MPa.iloc[0])
        r["slope_dlog10tr_dT_S1_per_K"] = slope(T1, sg1)
        for sname, d in scen.items():
            if sname == "S1_steady":
                continue
            v = np.sum(rate(d.T_wo_max_K.to_numpy(), d.sigma_hot_MPa.to_numpy(), 1.0)) / d.H2_net_kmol_h.sum() / s1
            r[f"{sname}_vs_S1"] = float(v)
            if k == 1.0:
                assert abs(v / summ.loc[sname, "life_consumption_per_kmol_H2_rel_S1"] - 1) < 1e-9
        rows.append(r)
        log(f"G k={k:g}: inside {inside}  slope {r['slope_dlog10tr_dT_base_per_K']:.4f}/K  knee {r['knee_vs_base']:.3f}  S5y4 {r['S5_hold_CH4_slip_ageing_y4_vs_S1']:.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "stress_multiplier.csv", index=False)
    _dump("calc_G.json", {"rows": rows, "base_reference": {"T_wo_max_K": rb["T_wo_max_K"], "sigma_hot_MPa": rb["sigma_hot_MPa"]},
                          "curve": {"alloy": CFG.creep_alloy, "kind": "minimum", "C": cur.C, "lmp_range": list(cur.lmp_range)}})
    return df


# ---------------------------------------------------------------------------
# H. Plant C1 as a held-out case
# ---------------------------------------------------------------------------
FEED_ATOMS = {"H2": (0, 2), "CO": (1, 0), "C1": (1, 4), "C2": (2, 6), "C3": (3, 8), "C4": (4, 10), "C5": (5, 12), "C6plus": (6, 14),
              "N2": (0, 0), "CO2": (1, 0), "H2O": (0, 2)}
OUT_ATOMS = {"H2": (0, 2), "CO": (1, 0), "CH4": (1, 4), "N2": (0, 0), "CO2": (1, 0), "H2O": (0, 2)}
H_QUANT = ("out_T_K", "TWT_upper_K", "TWT_lower_K", "flue_out_T_K", "out_P_Pa", "out_x_H2_wet_molpct", "out_x_CH4_wet_molpct",
           "out_x_CO_wet_molpct", "out_x_CO2_wet_molpct")
FREE3 = ("F_gt", "L_q", "f_htg")


def atom_balance(row: pd.Series, feed_kmol_h: Optional[float] = None) -> Dict[str, float]:
    n_in = float(row.feed_total_kmol_h if feed_kmol_h is None else feed_kmol_h); n_out = float(row.out_flow_kmol_h)
    cf = sum(row[f"feed_x_{k}_molpct"] / 100 * a[0] for k, a in FEED_ATOMS.items()); hf = sum(row[f"feed_x_{k}_molpct"] / 100 * a[1] for k, a in FEED_ATOMS.items())
    co = sum(row[f"out_x_{k}_wet_molpct"] / 100 * a[0] for k, a in OUT_ATOMS.items()); ho = sum(row[f"out_x_{k}_wet_molpct"] / 100 * a[1] for k, a in OUT_ATOMS.items())
    return {"feed_kmol_h": n_in, "out_kmol_h": n_out, "C_per_mol_feed": cf, "H_per_mol_feed": hf,
            "C_closure_out_over_in": n_out * co / (n_in * cf), "H_closure_out_over_in": n_out * ho / (n_in * hf),
            "feed_from_C_balance_kmol_h": n_out * co / cf, "feed_from_H_balance_kmol_h": n_out * ho / hf,
            "feed_sum_molpct": float(sum(row[f"feed_x_{k}_molpct"] for k in FEED_ATOMS)), "out_sum_molpct": float(sum(row[f"out_x_{k}_wet_molpct"] for k in OUT_ATOMS))}


def _adopted() -> "object":
    import yaml
    from rdt import latham_cases as lc
    doc = yaml.safe_load(lc.FIT_YAML.read_text())
    return lc.LathamParams(**doc["calibration_alpha_top_fixed"]["best_fit"]["full_parameters"]), doc["calibration_alpha_top_fixed"]


def _fit_task(rows, x0, base, names):
    from rdt import latham_calibration as cal
    f = cal.fit(rows, x0, base, names)
    return {k: f[k] for k in ("x", "params", "chi2", "n_obs", "jac", "fun", "status", "message", "n_evals", "x0", "at_bound", "names")}


def _loco_task(rows, k, x0, base, names):
    from rdt import latham_calibration as cal
    held = rows[k]; train = [r for j, r in enumerate(rows) if j != k]
    f = cal.fit(train, x0, base, names)
    err = cal.errors_table(f["x"], [held], base, names)[str(held.case)]
    return str(held.case), {"params": f["params"], "chi2_train": f["chi2"], "held_out_errors": err, "at_bound": f["at_bound"]}


def calc_H(log=print) -> pd.DataFrame:
    from dataclasses import replace
    from joblib import Parallel, delayed
    from rdt import latham_calibration as cal
    from rdt import latham_cases as lc
    cases = lc.load_cases().set_index("case", drop=False)
    # (1) atom balances, (2) C1 feed from the carbon balance
    bal = {c: atom_balance(cases.loc[c]) for c in cases.index}
    n_c1 = bal["Plant_C1"]["feed_from_C_balance_kmol_h"]
    bal["Plant_C1_restored"] = atom_balance(cases.loc["Plant_C1"], n_c1)
    for c, v in bal.items():
        log(f"H balance {c:18s} feed {v['feed_kmol_h']:8.1f}  C out/in {v['C_closure_out_over_in']:.4f}  H out/in {v['H_closure_out_over_in']:.4f}  "
            f"feed(C) {v['feed_from_C_balance_kmol_h']:.1f}")
    c1 = cases.loc["Plant_C1"].copy(); c1["feed_total_kmol_h"] = n_c1; c1["feed_per_tube_kmol_h"] = n_c1 / float(c1.n_tubes)
    # (3) prediction with the adopted parameters, nothing refitted
    P, cal_doc = _adopted()
    rows_out, ea = [], {}
    def predict(label, row, kind):
        tube, bed, feed, flue, geom = lc.build_case(row, P); ea[label] = float(flue.excess_air_pct)
        res = lc.run_case(row, P); c = lc.compare(row, res)
        r = {"case": label, "kind": kind, "feed_kmol_h": float(row.feed_total_kmol_h), "excess_air_pct": ea[label],
             "excess_air_pct_csv": float(row.excess_air_pct_derived)}
        for q in H_QUANT:
            sc_ = 1e-3 if q == "out_P_Pa" else 1.0
            r[f"{q}_meas"] = c[q]["measured"] * sc_; r[f"{q}_pred"] = c[q]["predicted"] * sc_; r[f"{q}_err"] = c[q]["difference"] * sc_
        r["TWT_peak_frac"] = c["TWT_peak"]["z_frac"]; r["energy_closure"] = c["energy"]["closure_rel_error"]
        return r
    for cname in ("Plant_A", "Plant_B", "Plant_C2"):
        rows_out.append(predict(cname, cases.loc[cname], "calibration case (in-sample)"))
    rows_out.append(predict("Plant_C1_restored", c1, "held out, feed from C balance"))
    rows_out.append(predict("Plant_C1_printed", cases.loc["Plant_C1"], "held out, printed feed (for contrast)"))
    for r in rows_out:
        log(f"H {r['case']:18s} dT_out {r['out_T_K_err']:+6.1f}  dTw,u {r['TWT_upper_K_err']:+6.1f}  dTw,l {r['TWT_lower_K_err']:+6.1f}  "
            f"dT_fg {r['flue_out_T_K_err']:+6.1f}  dP {r['out_P_Pa_err']:+6.1f} kPa  dH2 {r['out_x_H2_wet_molpct_err']:+5.2f}  EA {r['excess_air_pct']:.2f}")
    # (4) sensitivity: refit on A, B, C1 restored, C2 (alpha_top fixed, same multistart rule as the adopted fit)
    base = replace(cal.BASE, alpha_top=0.182)
    rows4 = [cases.loc["Plant_A"], cases.loc["Plant_B"], c1, cases.loc["Plant_C2"]]
    rows3 = [cases.loc["Plant_A"], cases.loc["Plant_B"], cases.loc["Plant_C2"]]
    chi3 = float(np.sum(cal.residuals([P.F_gt, P.L_q, P.f_htg], rows3, base, names=FREE3) ** 2))
    starts = [np.array([0.35, 0.48, 1.68])] + list(cal.lhs_starts(5, 7, FREE3)) + [np.array([P.F_gt, P.L_q, P.f_htg])]
    t0 = time.perf_counter()
    fits = Parallel(n_jobs=len(starts))(delayed(_fit_task)(rows4, s, base, FREE3) for s in starts)
    fits.sort(key=lambda f: f["chi2"]); best = fits[0]
    cov = cal.covariance(best)
    log(f"H 4-case fit: {best['params']}  chi2 {best['chi2']:.2f}/{cov['dof']}  ({time.perf_counter() - t0:.0f}s)")
    loco = dict(Parallel(n_jobs=4)(delayed(_loco_task)(rows4, k, best["x"], base, FREE3) for k in range(4)))
    err4 = cal.errors_table(best["x"], rows4, base, FREE3)
    for c, d in loco.items():
        log(f"H LOCO {c}: dTw,u {d['held_out_errors']['TWT_upper_K']:+.1f}  dTw,l {d['held_out_errors']['TWT_lower_K']:+.1f}  params {d['params']}")
    df = pd.DataFrame(rows_out)
    df.to_csv(OUT / "c1_heldout.csv", index=False)
    pd.DataFrame(bal).T.to_csv(OUT / "c1_atom_balances.csv")
    fit4 = {"params": best["params"], "chi2": best["chi2"], "n_obs": best["n_obs"], "dof": cov["dof"], "reduced_chi2": cov["reduced_chi2"],
            "ci95": cov["ci95"], "se": cov["se"], "corr": cov["corr_table"], "at_bound": best["at_bound"],
            "multistart": [{"x0": list(map(float, f["x0"])), "params": f["params"], "chi2": f["chi2"]} for f in fits],
            "in_sample_errors": err4, "loco": loco}
    pub = {"params": cal_doc["best_fit"]["params"], "chi2": cal_doc["best_fit"]["chi2"], "ci95": cal_doc["best_fit"]["covariance"]["ci95"],
           "reduced_chi2": cal_doc["best_fit"]["covariance"]["reduced_chi2"], "dof": cal_doc["best_fit"]["covariance"]["dof"],
           "chi2_recomputed": chi3,
           "loco": {c: {q: v["held_out_errors"][q] for q in ("TWT_upper_K", "TWT_lower_K")} for c, v in cal_doc["leave_one_out_cv"].items()}}
    _dump("calc_H.json", {"balances": bal, "C1_feed_restored_kmol_h": n_c1, "excess_air_pct": ea, "adopted_fit_3case": pub, "fit_4case": fit4})
    return df


# ===========================================================================
# Package 3
# ===========================================================================
TABLE2_LABELS = {"S1_steady": "S1 steady", "S2_daily": "S2 daily", "S3_renewable": "S3 renewable", "S4a_mild_overfire": "S4a mild",
                 "S4b_severe_overfire": "S4b severe", "S4c_hot_band": "S4c hot band", "S6_SC2.5": "S6 S/C 2.5", "S6_SC3.5": "S6 S/C 3.5",
                 "S5_hold_CH4_slip_ageing_y1": "S5 slip y1", "S5_hold_CH4_slip_ageing_y4": "S5 slip y4",
                 "S5_hold_T_out_ageing_y1": "S5 T_out y1", "S5_hold_T_out_ageing_y4": "S5 T_out y4"}


# ---------------------------------------------------------------------------
# I. Wall-temperature error propagated to the ratios (re-scoring only)
# ---------------------------------------------------------------------------
LOCO_U = ((0.99, -12.8), (0.71, +18.2))          # (plant rate, LOCO upper-wall error, K): Plant A, Plant C2
W_HOT = (0.682 - 0.43) / (0.682 - 0.293)


def e_upper(load) -> np.ndarray:
    (x1, y1), (x2, y2) = LOCO_U
    return y1 + (y2 - y1) * (np.asarray(load, float) - x1) / (x2 - x1)


def calc_I(log=print) -> pd.DataFrame:
    from rdt import creep
    from rdt import operate as op
    from rdt import optimize as opt
    cur = creep.load_alloy(CFG.creep_alloy).minimum
    meta = json.loads(Path(CFG.pareto_meta).read_text()); bref = meta["base_reference"]
    rb = op.run_case(bref["inputs"], op.base_case())
    par = pd.read_csv(CFG.pareto_csv); feas = par[par.feasible_phys == True].copy()  # noqa: E712
    reg = pd.read_csv(CFG.regimes_csv, index_col=0); ins = list(op.INPUT_NAMES)
    reg_row = {}
    for name in TABLE3_REGIMES:
        m = np.all(np.isclose(feas[ins].to_numpy(float), reg.loc[name, ins].to_numpy(float), rtol=0, atol=1e-12), axis=1)
        reg_row[name] = feas.index[m][0]
    knee0 = reg_row["knee"]
    hours = {s: pd.read_csv(CFG.scenario_hourly(s), usecols=["feed_per_tube_fraction", "T_wo_max_K", "sigma_hot_MPa", "H2_net_kmol_h"]) for s in TABLE2}
    summ = pd.read_csv(CFG.scenarios_summary).set_index("scenario")
    variants = {"a: dT = -18 K": ("uniform", -18.0), "a: dT = -12 K": ("uniform", -12.0), "a: dT = +12 K": ("uniform", 12.0), "a: dT = +18 K": ("uniform", 18.0),
                "b: w*e_u(load), sign as is": ("load", W_HOT), "b: -w*e_u(load), reversed": ("load", -W_HOT), "b: 1*e_u(load), w = 1": ("load", 1.0)}
    def dT(kind, par_, load):
        return np.full(np.shape(load), par_, float) if kind == "uniform" else par_ * e_upper(load)
    rate = lambda T, s, H: 1.0 / np.asarray(cur.time_to_rupture(T, s), float) / H  # noqa: E731
    rows, nd_info = [], {}
    F0 = np.column_stack([-feas.H2_net_kmol_h_phys, feas.fuel_MJ_per_kmol_H2_phys, feas.life_rate_per_kmol_H2_rel_base_phys])
    nd0 = nondominated(F0)
    for vname, (kind, pv) in [("none (v3)", ("uniform", 0.0))] + list(variants.items()):
        # Table 2: hourly histories, each hour shifted by its own load
        D = {}
        for s, h in hours.items():
            T = h.T_wo_max_K.to_numpy() + dT(kind, pv, h.feed_per_tube_fraction.to_numpy())
            D[s] = np.sum(rate(T, h.sigma_hot_MPa.to_numpy(), 1.0)) / h.H2_net_kmol_h.sum()
        for s in TABLE2:
            v = D[s] / D["S1_steady"]
            if vname == "none (v3)":
                assert abs(v / summ.loc[s, "life_consumption_per_kmol_H2_rel_S1"] - 1) < 1e-9
            rows.append({"variant": vname, "table": "Table 2", "item": TABLE2_LABELS[s], "rel_life_cost": float(v),
                         "dT_item_K": float(np.mean(dT(kind, pv, hours[s].feed_per_tube_fraction.to_numpy()))),
                         "dT_ref_K": float(dT(kind, pv, np.array([1.0]))[0])})
        # Table 3 regimes and the 198 members, base_reference at load 1.0
        dT_b = float(dT(kind, pv, np.array([float(bref["inputs"]["feed_per_tube_fraction"])]))[0])
        rb_rate = float(rate(rb["T_wo_max_K"] + dT_b, rb["sigma_hot_MPa"], rb["H2_net_kmol_h"]))
        ld = feas.feed_per_tube_fraction.to_numpy()
        life = rate(feas.T_wo_max_K_phys.to_numpy() + dT(kind, pv, ld), feas.sigma_hot_MPa_phys.to_numpy(), feas.H2_net_kmol_h_phys.to_numpy()) / rb_rate
        lf = pd.Series(life, index=feas.index)
        for name in TABLE3_REGIMES:
            i = reg_row[name]
            if vname == "none (v3)":
                assert abs(lf[i] / feas.loc[i, "life_rate_per_kmol_H2_rel_base_phys"] - 1) < 1e-9
            rows.append({"variant": vname, "table": "Table 3", "item": name, "rel_life_cost": float(lf[i]),
                         "dT_item_K": float(dT(kind, pv, np.array([feas.loc[i, "feed_per_tube_fraction"]]))[0]), "dT_ref_K": dT_b})
        if kind == "load" or vname == "none (v3)":
            F = np.column_stack([-feas.H2_net_kmol_h_phys, feas.fuel_MJ_per_kmol_H2_phys, life])
            nd = nondominated(F)
            tmp = feas.assign(_life=life)
            knee = opt.knee_point(tmp, cols=("H2_net_kmol_h_phys", "fuel_MJ_per_kmol_H2_phys", "_life"))
            nd_info[vname] = {"n_nondominated": int(nd.sum()), "lost": int(np.sum(nd0 & ~nd)), "gained": int(np.sum(~nd0 & nd)),
                              "knee_row": int(knee), "knee_same_as_v3": bool(knee == knee0),
                              "knee": {k: float(feas.loc[knee, k]) for k in list(opt.DECISION)},
                              "knee_H2": float(feas.loc[knee, "H2_net_kmol_h_phys"]), "knee_fuel": float(feas.loc[knee, "fuel_MJ_per_kmol_H2_phys"]),
                              "knee_rel_life": float(lf[knee]), "v3_knee_rel_life_this_variant": float(lf[knee0]),
                              "knee_decision_distance_norm": float(np.linalg.norm((feas.loc[knee, list(opt.DECISION)].to_numpy(float) - feas.loc[knee0, list(opt.DECISION)].to_numpy(float)) / (opt.XU - opt.XL)))}
        log(f"I {vname:28s} " + "  ".join(f"{r['item']} {r['rel_life_cost']:.3f}" for r in rows[-17:] if r["item"] in ("S2 daily", "S5 slip y4", "knee", "min_life_iso_H2"))
            + (f"  nd {nd_info[vname]['n_nondominated']} lost {nd_info[vname]['lost']} knee row {nd_info[vname]['knee_row']}" if vname in nd_info else ""))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "walltemp_error_propagation.csv", index=False)
    _dump("calc_I.json", {"e_u_line": {"points": LOCO_U, "slope_K_per_load": (LOCO_U[1][1] - LOCO_U[0][1]) / (LOCO_U[1][0] - LOCO_U[0][0]),
                                       "e_u_at_load": {f"{x:.2f}": float(e_upper(x)) for x in (0.6, 0.7, 0.85, 0.97, 1.0, 1.1)}},
                          "w": W_HOT, "nondominance_and_knee": nd_info, "v3_knee_row": int(knee0)})
    return df


# ---------------------------------------------------------------------------
# J/K. Flame length and bed-side coefficient: refit, then physics ratios
# ---------------------------------------------------------------------------
_BASE_CACHE: Dict[str, object] = {}


def _run_P(x: Dict[str, float], P) -> Dict[str, float]:
    from dataclasses import replace
    from rdt import creep
    from rdt import operate as op
    if "base" not in _BASE_CACHE:
        _BASE_CACHE["base"] = op.base_case()
    return op.run_case(x, replace(_BASE_CACHE["base"], params=P), creep.curve_for(CFG.creep_alloy))


def _solve_P(x: Dict[str, float], P, key: str, target: float, lo: float = 0.75, hi: float = 1.25) -> Dict[str, float]:
    from scipy.optimize import brentq
    c: Dict[float, Dict[str, float]] = {}
    def run(f):
        if f not in c:
            c[f] = _run_P({**x, "specific_firing_factor": float(f)}, P)
        return c[f]
    g = lambda f: run(f)[key] - target  # noqa: E731
    a, b = g(lo), g(hi)
    if a * b > 0:
        f = lo if abs(a) < abs(b) else hi; sat = True
    else:
        f = brentq(g, lo, hi, xtol=1e-6); sat = False
    r = dict(run(f)); r["firing"] = float(f); r["firing_saturated"] = sat
    return r


def _sweep_task(fixed: str, value: float, free: Tuple[str, ...], log_prefix: str = "") -> Dict[str, object]:
    """Refit `free` on A, B, C2 with `fixed` = value (alpha_top = 0.182), then the physics ratios."""
    from dataclasses import replace
    from rdt import latham_calibration as cal
    from rdt import latham_cases as lc
    from rdt import operate as op
    from rdt import scenarios as sc
    setup()
    P0, _ = _adopted()
    base_cal = replace(cal.BASE, alpha_top=0.182, **{fixed: value})
    rows3 = cal.case_rows(("Plant_A", "Plant_B", "Plant_C2"))
    starts = [[getattr(P0, k) for k in free], [{"F_gt": 0.35, "L_q": 0.48, "f_htg": 1.68}[k] for k in free]]
    fits = [cal.fit(rows3, s, base_cal, free) for s in starts]; f = min(fits, key=lambda q: q["chi2"])
    cov = cal.covariance(f)
    err3 = cal.errors_table(f["x"], rows3, base_cal, free)
    P = replace(P0, **{fixed: value}, **dict(zip(free, map(float, f["x"]))))
    cases = lc.load_cases().set_index("case", drop=False)
    c1 = cases.loc["Plant_C1"].copy(); n1 = atom_balance(c1)["feed_from_C_balance_kmol_h"]
    c1["feed_total_kmol_h"] = n1; c1["feed_per_tube_kmol_h"] = n1 / float(c1.n_tubes)
    cc1 = lc.compare(c1, lc.run_case(c1, P))
    out = {"fixed": fixed, "value": value, "free": list(free), "params": dict(zip(free, map(float, f["x"]))), "chi2": f["chi2"], "dof": cov["dof"],
           "reduced_chi2": cov["reduced_chi2"], "ci95": cov["ci95"], "at_bound": f["at_bound"],
           "err_wall": {c: {"TWT_upper_K": e["TWT_upper_K"], "TWT_lower_K": e["TWT_lower_K"], "out_T_K": e["out_T_K"]} for c, e in err3.items()},
           "C1_restored": {q: cc1[q]["difference"] for q in ("out_T_K", "TWT_upper_K", "TWT_lower_K", "flue_out_T_K")}}
    # Plant A base point (firing 1.0, base_reference inputs)
    meta = json.loads(Path(CFG.pareto_meta).read_text()); bx = meta["base_reference"]["inputs"]
    rb = _run_P(bx, P)
    out["base"] = {"T_wo_max_K": rb["T_wo_max_K"], "z_frac_T_wo_max": rb["z_frac_T_wo_max"], "T_out_K": rb["T_out_K"], "converged": rb["converged"]}
    rate_b = rb["life_consumption_rate_per_h"] / rb["H2_net_kmol_h"]
    reg = pd.read_csv(CFG.regimes_csv, index_col=0)
    out["regimes"] = {}
    for name in TABLE3_REGIMES:
        r = _run_P({k: float(reg.loc[name, k]) for k in op.INPUT_NAMES}, P)
        out["regimes"][name] = (r["life_consumption_rate_per_h"] / r["H2_net_kmol_h"]) / rate_b if r["converged"] else float("nan")
    # S1 (hold T_out, physics) and its slip; S5 hold slip year 4 in 6 bins; S6 at hold T_out
    b = sc.base_inputs(1)
    x0 = {"steam_to_carbon": float(b["steam_to_carbon"][0]), "feed_per_tube_fraction": 1.0, "inlet_T": float(b["inlet_T"][0]),
          "inlet_P": float(b["inlet_P"][0]), "catalyst_activity": A_CAL, "excess_air": float(b["excess_air"][0])}
    s1 = _solve_P(x0, P, "T_out_K", sc.T_OUT_TARGET)
    r1_ = s1["life_consumption_rate_per_h"] / s1["H2_net_kmol_h"]
    act = activity_profile("E0")[3 * HOURS:4 * HOURS]; edges = np.linspace(0, HOURS, 7).astype(int)
    D = H = 0.0; sat = 0
    for j in range(6):
        a = float(act[edges[j]:edges[j + 1]].mean()); n = edges[j + 1] - edges[j]
        r = _solve_P({**x0, "catalyst_activity": a}, P, "CH4_slip_dry_pct", s1["CH4_slip_dry_pct"])
        D += n * r["life_consumption_rate_per_h"]; H += n * r["H2_net_kmol_h"]; sat += int(r["firing_saturated"])
    out["S5_slip_y4"] = (D / H) / r1_
    for scv in (2.5, 3.5):
        r = _solve_P({**x0, "steam_to_carbon": scv}, P, "T_out_K", sc.T_OUT_TARGET)
        out[f"S6_SC{scv}"] = (r["life_consumption_rate_per_h"] / r["H2_net_kmol_h"]) / r1_; sat += int(r["firing_saturated"])
    out["S1"] = {"firing": s1["firing"], "T_wo_max_K": s1["T_wo_max_K"], "slip": s1["CH4_slip_dry_pct"]}
    out["n_firing_saturated"] = sat + int(s1["firing_saturated"])
    return out


def _sweep_table(res: List[Dict[str, object]]) -> pd.DataFrame:
    rows = []
    for o in res:
        r = {"fixed": o["fixed"], "value": o["value"], **{f"fit_{k}": v for k, v in o["params"].items()},
             **{f"ci95_{k}": f"{v[0]:.4f}-{v[1]:.4f}" for k, v in o["ci95"].items()},
             "chi2": o["chi2"], "dof": o["dof"], "reduced_chi2": o["reduced_chi2"], "any_at_bound": any(o["at_bound"].values())}
        for c, e in o["err_wall"].items():
            r[f"{c}_dTw_u"] = e["TWT_upper_K"]; r[f"{c}_dTw_l"] = e["TWT_lower_K"]; r[f"{c}_dT_out"] = e["out_T_K"]
        for q, v in o["C1_restored"].items():
            r[f"C1_{q}_err"] = v
        r.update({"base_T_wo_max_K": o["base"]["T_wo_max_K"], "base_z_frac_T_wo_max": o["base"]["z_frac_T_wo_max"], "base_T_out_K": o["base"]["T_out_K"]})
        r.update({f"rel_{k}": v for k, v in o["regimes"].items()})
        r.update({"rel_S5_slip_y4": o["S5_slip_y4"], "rel_S6_SC2.5": o["S6_SC2.5"], "rel_S6_SC3.5": o["S6_SC3.5"],
                  "S1_firing_phys": o["S1"]["firing"], "n_firing_saturated": o["n_firing_saturated"]})
        rows.append(r)
    return pd.DataFrame(rows)


L_Q_GRID = (0.35, 0.40, 0.456, 0.50, 0.55, 0.60)
F_HTG_GRID = (0.5, 1.0, 1.83, 3.0)


def calc_J(log=print) -> pd.DataFrame:
    from joblib import Parallel, delayed
    t0 = time.perf_counter()
    res = Parallel(n_jobs=len(L_Q_GRID))(delayed(_sweep_task)("L_q", v, ("F_gt", "f_htg")) for v in L_Q_GRID)
    df = _sweep_table(res)
    df.to_csv(OUT / "flame_length_sweep.csv", index=False)
    _dump("calc_J.json", {"runs": res})
    for o in res:
        log(f"J L_q {o['value']:.3f}: {o['params']} chi2/dof {o['reduced_chi2']:.3f}  T_wo {o['base']['T_wo_max_K']:.1f} z/L {o['base']['z_frac_T_wo_max']:.3f}  "
            f"knee {o['regimes']['knee']:.3f} S5y4 {o['S5_slip_y4']:.3f}  C1 dTw,u {o['C1_restored']['TWT_upper_K']:+.1f}")
    log(f"J done in {time.perf_counter() - t0:.0f}s")
    fig_J(df)
    return df


def fig_J(df: pd.DataFrame) -> None:
    from rdt import plotting as P
    x = df.value.to_numpy()
    with P.style():
        fig, ax = P.figure(P.DOUBLE, 5.5 * P.CM, 1, 3)
        ax[0].plot(x, df.reduced_chi2, "o-", color=P.PALETTE[0], ms=3); ax[0].axvline(0.456, color="grey", lw=0.5, ls="--")
        P.tidy(ax[0], "$L_q$ (–)", "$\\chi^2$/dof (A, B, C2)")
        ax[1].plot(x, df.base_T_wo_max_K, "o-", color=P.PALETTE[1], ms=3, label="$T_{wo,max}$"); ax[1].axvline(0.456, color="grey", lw=0.5, ls="--")
        P.tidy(ax[1], "$L_q$ (–)", "$T_{wo,max}$ at Plant A base (K)")
        a2 = ax[1].twinx(); a2.plot(x, df.base_z_frac_T_wo_max, "s--", color=P.PALETTE[2], ms=3, label="$z/L$"); a2.set_ylabel("$z/L$ of $T_{wo,max}$ (–)")
        h1, l1 = ax[1].get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels(); ax[1].legend(h1 + h2, l1 + l2, loc="upper center")
        series = (("rel_knee", "knee / base", 3), ("rel_min_life_iso_H2", "min life iso-H$_2$ / base", 2), ("rel_max_H2", "max H$_2$ / base", 4),
                  ("rel_S5_slip_y4", "S5 slip y4 / S1", 1), ("rel_S6_SC2.5", "S/C 2.5 / S1", 5), ("rel_S6_SC3.5", "S/C 3.5 / S1", 0))
        for c, lab, k in series:
            ax[2].plot(x, df[c], "o-", ms=3, color=P.PALETTE[k], label=lab)
        ax[2].axvline(0.456, color="grey", lw=0.5, ls="--"); ax[2].axhline(1.0, color="grey", lw=0.5); ax[2].set_yscale("log"); ax[2].set_ylim(3e-3, 8)
        P.tidy(ax[2], "$L_q$ (–)", "Life consumption per kmol H$_2$ (rel.)", True, fontsize=5, ncol=2, loc="lower center")
        for a, l in zip(ax, "abc"):
            P.panel_label(a, l)
        fig.tight_layout(); P.save(fig, OUT / "flame_length_sweep")


def htc_chain_terms() -> pd.DataFrame:
    """Term-by-term bed-side heat-transfer chain at mid-tube, Xu-Froment 1989 case and Plant A."""
    import inspect
    from rdt import reactor1d as r1
    from rdt import latham_cases as lc
    from rdt import validate_xf1989 as xf
    src = Path(r1.__file__); lines = src.read_text().splitlines()
    def ln(pat, f=lines, path="src/rdt/reactor1d.py", start=0):
        i = next(k for k, s in enumerate(f) if k >= start and pat in s); return f"{path}:{i + 1}"
    xl = Path(xf.__file__).read_text().splitlines(); ll = Path(lc.__file__).read_text().splitlines()
    # --- states -------------------------------------------------------------------------------
    tidy = xf.load_tidy(); cv = xf.curves(tidy); wall = xf.wall_from_digitised(tidy)
    prof = xf.alpha_i_profile_from_digitised(cv)
    res_xf = xf.run_matched(prof["profile"], wall)
    tube_xf = r1.tube_from_xu_froment()
    bed_xf = r1.CatalystBed(rho_bed=1.0, voidage=xf.STAGE_A_BED["voidage"], d_p=xf.STAGE_A_BED["d_p_m"], eta=(0.1,) * 3, activity=1.0)
    P_ad, _ = _adopted()
    rowA = lc.load_cases().set_index("case").loc["Plant_A"]
    resA = lc.run_case(rowA, P_ad)
    tube_A, bed_A, *_ = lc.build_case(rowA, P_ad)
    def state(res, tube, z_mid):
        tr = res.tube if hasattr(res, "tube") else res
        k = int(np.argmin(np.abs(tr.z - z_mid)))
        F = np.array([tr.F[s][k] for s in r1.SPECIES]); X = {s: F[i] / F.sum() for i, s in enumerate(r1.SPECIES)}
        props = r1._gas_props(tr.T[k], tr.P[k], X); _, _, mw = r1._gas()
        m_dot = sum(F[i] * mw[s] for i, s in enumerate(r1.SPECIES)) / 3600.0
        return tr, k, props, m_dot / tube.A_cs
    out = []
    ring = r1.ring_equivalent_diameter(0.0173, 0.0084, 0.010)
    for label, res, tube, bed, z_mid, f_cal in (("Xu-Froment 1989", res_xf, tube_xf, bed_xf, 11.12 / 2, None),
                                               ("Plant A", resA, tube_A, bed_A, P_ad.heated_length_m / 2, P_ad.f_htg)):
        tr, k, pr, G = state(res, tube, z_mid)
        lam, mu, cp = pr["lam"], pr["mu"], pr["cp_mass"]
        lam0 = r1.kunii_smith_lambda_er0(lam, bed.lambda_s, bed.voidage, float(tr.T[k]), bed.emissivity, bed.d_p)
        ch = r1.xu_froment_alpha_i(lam, mu, cp, G, bed.d_p, tube.d_i, lam0)
        lg1 = r1.leva_grummer_alpha(lam, mu, G, bed.d_p, tube.d_i, 1.0)
        Re, Pr = ch["Re"], ch["Pr"]
        if label == "Plant A":
            backcalc = (float("nan"), "not applicable (no inner-wall measurement); the calibrated value is the next row", "")
            calv = (float(tr.alpha_i[k]), f"alpha_i of the calibrated run = f_htg * LG = {f_cal:.4f} * {lg1:.1f}; check {f_cal * lg1:.1f}", ln("alpha_i = leva_grummer_alpha("))
            dp_def = f"thesis pellet equivalent diameter 5.40 mm (LathamParams.d_p_m)"; dp_ref = f"src/rdt/latham_cases.py:{next(i for i, s in enumerate(ll) if 'd_p_m: float' in s) + 1}"
        else:
            j = int(np.argmin(np.abs(prof["z"] - z_mid)))
            backcalc = (float(prof["alpha_i"][j]), f"q_i/(T_wi - T_gas) from digitised Fig. 3, lambda_tube 29.6; at z = {prof['z'][j]:.2f} m (median over 0.5-11 m: {prof['median']:.0f})",
                        f"src/rdt/validate_xf1989.py:{next(i for i, s in enumerate(xl) if 'def alpha_i_profile_from_digitised' in s) + 1}")
            calv = (float("nan"), "not applicable (Xu-Froment case is not calibrated with f_htg)", "")
            dp_def = (f"9.44 mm from the stage-A pressure-drop fit (STAGE_A_BED, used in the runs); the 6V/S ring definition "
                      f"(ring_equivalent_diameter, 17.3/8.4/10 mm) gives {1e3 * ring:.2f} mm")
            dp_ref = f"src/rdt/validate_xf1989.py:{next(i for i, s in enumerate(xl) if 'STAGE_A_BED =' in s) + 1}; " + ln("def ring_equivalent_diameter")
        R = [("z (mid heated length)", float(tr.z[k]), "m", "node nearest L_heated/2", ""),
             ("T_gas", float(tr.T[k]), "K", "1-D model state", ""), ("P", float(tr.P[k]), "bar", "1-D model state", ""),
             ("d_t inner", tube.d_i, "m", "tube internal diameter", ln("def A_cs")),
             ("d_p equivalent", bed.d_p, "m", dp_def, dp_ref),
             ("voidage", bed.voidage, "-", "bed voidage", ""),
             ("G (superficial mass velocity)", G, "kg m-2 s-1", "G = m_dot / (pi d_t^2 / 4)", ln("G_s = m_dot / tube.A_cs")),
             ("mu_g", mu, "Pa s", "Cantera gri30 mixture-averaged viscosity", ln('"mu": gas.viscosity')),
             ("lambda_g", lam, "W m-1 K-1", "Cantera gri30 mixture-averaged conductivity", ln('"lam": gas.thermal_conductivity')),
             ("cp_g", cp, "J kg-1 K-1", "Cantera cp_mass", ln('"cp_mass": gas.cp_mass')),
             ("Re_p", Re, "-", "Re_p = d_p G / mu", ln("Re = d_p * G_s / mu", start=next(k for k, q in enumerate(lines) if "def xu_froment_alpha_i" in q))),
             ("Pr", Pr, "-", "Pr = cp mu / lambda_g", ln("Pr = cp_mass * mu / lam_g")),
             ("lambda_er0 (static, Kunii-Smith)", lam0, "W m-1 K-1",
              f"Kunii-Smith 1960 with lambda_s = {bed.lambda_s}, emissivity {bed.emissivity}, beta 1, gamma 2/3; radiation h_rv, h_rs ~ 4 sigma T^3", ln("def kunii_smith_lambda_er0")),
             ("lambda_er dynamic part", 0.14 * lam * Re * Pr, "W m-1 K-1", "0.14 lambda_g Re_p Pr (Yagi-Kunii form)", ln("lam_er = lam_er0 + 0.14 * lam_g")),
             ("lambda_er total", ch["lam_er"], "W m-1 K-1", "lambda_er0 + 0.14 lambda_g Re_p Pr", ln("lam_er = lam_er0 + 0.14 * lam_g")),
             ("alpha_w0 (static, De Wasch-Froment)", ch["alpha_w0"], "W m-2 K-1", "8.694 lambda_er0 / d_t^(4/3)", ln("alpha_w0 = 8.694 * lam_er0")),
             ("alpha_w dynamic part", 0.444 * Re * Pr * lam / bed.d_p, "W m-2 K-1", "0.444 Re_p Pr lambda_g / d_p", ln("alpha_w = alpha_w0 + 0.444 * Re")),
             ("alpha_w total", ch["alpha_w"], "W m-2 K-1", "alpha_w0 + 0.444 Re_p Pr lambda_g / d_p", ln("alpha_w = alpha_w0 + 0.444 * Re")),
             ("alpha_i, Xu-Froment chain", ch["alpha_i"], "W m-2 K-1", "1/alpha_i = 1/alpha_w + d_t/(8 lambda_er)", ln("alpha_i = 8.0 * lam_er * alpha_w")),
             ("alpha_i, Leva-Grummer (f_htg = 1)", lg1, "W m-2 K-1", "0.813 (lambda_g/d_t) exp(-6 d_p/d_t) Re_p^0.9", ln("return f_htg * 0.813")),
             ("alpha_i, back-calculated", backcalc[0], "W m-2 K-1", backcalc[1], backcalc[2]),
             ("alpha_i, calibrated", calv[0], "W m-2 K-1", calv[1], calv[2]),
             ("ratio chain / Leva-Grummer", ch["alpha_i"] / lg1, "-", "", "")]
        out += [{"case": label, "term": a, "value": v, "unit": u, "formula_or_definition": fo, "code_ref": rf} for a, v, u, fo, rf in R]
    return pd.DataFrame(out)


def calc_K(log=print) -> pd.DataFrame:
    import yaml
    from joblib import Parallel, delayed
    from rdt import latham_cases as lc
    t0 = time.perf_counter()
    res = Parallel(n_jobs=len(F_HTG_GRID))(delayed(_sweep_task)("f_htg", v, ("F_gt", "L_q")) for v in F_HTG_GRID)
    df = _sweep_table(res)
    df.to_csv(OUT / "htc_sweep.csv", index=False)
    for o in res:
        log(f"K f_htg {o['value']:.2f}: {o['params']} chi2/dof {o['reduced_chi2']:.3f} at_bound {o['at_bound']}  T_wo {o['base']['T_wo_max_K']:.1f} z/L {o['base']['z_frac_T_wo_max']:.3f}  "
            f"knee {o['regimes']['knee']:.3f} S5y4 {o['S5_slip_y4']:.3f}")
    doc = yaml.safe_load(lc.FIT_YAML.read_text())["calibration_alpha_top_fixed"]["best_fit"]["covariance"]
    corr = {"f_htg~F_gt": doc["corr_table"]["f_htg"]["F_gt"], "f_htg~L_q": doc["corr_table"]["f_htg"]["L_q"], "F_gt~L_q": doc["corr_table"]["F_gt"]["L_q"]}
    chain = htc_chain_terms(); chain.to_csv(OUT / "htc_chain_terms.csv", index=False)
    _dump("calc_K.json", {"runs": res, "correlations_adopted_fit": corr, "high_correlation_pairs": doc.get("high_correlation_pairs")})
    log(f"K correlations {corr}  ({time.perf_counter() - t0:.0f}s)")
    return df


# ---------------------------------------------------------------------------
# L. Load-following magnitude with T_out held by physics at every p
# ---------------------------------------------------------------------------
def _flex_point_task(point, p):
    from scipy.optimize import brentq
    from rdt import uq
    _ensure_setup()
    g = lambda f: uq.evaluate_one({**point, "specific_firing_factor": float(f)}, _NOM, float(p))["T_out_K"] - 1105.55  # noqa: E731
    f = brentq(g, 0.75, 1.25, xtol=1e-7)
    r = uq.evaluate_one({**point, "specific_firing_factor": f}, _NOM, float(p))
    return {"firing": float(f), "T_wo_max_K": r["T_wo_max_K"], "T_out_K": r["T_out_K"], "nominal": r}


def _flex_mc_worker(point, theta, p):
    from rdt import uq
    _ensure_setup()
    r = uq.evaluate_one(point, theta, float(p))
    return {k: r.get(k, np.nan) for k in ("converged", "life_rate_per_h", "H2_net_kmol_h", "life_rate_per_kmol_H2", "T_wo_max_K", "T_out_K")}


def calc_L(n: int = 512, n_jobs: int = -1, log=print) -> pd.DataFrame:
    from joblib import Parallel, delayed
    from rdt import uq
    _ensure_setup()
    pts = json.loads((Path(uq.OUT_DIR) / f"{uq.MC_TAG}_points.json").read_text())
    X = mc_design(n); th = X.to_dict(orient="records")
    b = pd.read_csv(Path(uq.OUT_DIR) / f"{uq.MC_TAG}_S1_base.csv.gz").iloc[:n].reset_index(drop=True)
    nb = uq.evaluate_one(pts["S1_base"], _NOM)
    sol = Parallel(n_jobs=n_jobs)(delayed(_flex_point_task)({k: v for k, v in pts[nm].items()}, p) for p in P_B for nm in ("load_0.70", "load_0.85"))
    solved = {(p, nm): s for (p, nm), s in zip([(p, nm) for p in P_B for nm in ("load_0.70", "load_0.85")], sol)}
    d = OUT / "flex_holdTout_mc"; d.mkdir(exist_ok=True)
    jobs = [(p, nm) for p in P_B for nm in ("load_0.70", "load_0.85") if not (d / f"{nm}_{_tag_p(p)}.csv.gz").exists()]
    if jobs:
        t0 = time.perf_counter()
        res = Parallel(n_jobs=n_jobs, batch_size=8)(delayed(_flex_mc_worker)({**pts[nm], "specific_firing_factor": solved[(p, nm)]["firing"]}, th[i], p)
                                                    for p, nm in jobs for i in range(n))
        for j, (p, nm) in enumerate(jobs):
            pd.DataFrame(res[j * n:(j + 1) * n]).to_csv(d / f"{nm}_{_tag_p(p)}.csv.gz", index=False, compression="gzip")
        log(f"L MC {len(jobs)} x {n} in {time.perf_counter() - t0:.0f}s")
    rows = []
    one = lambda dct: pd.DataFrame([dct])  # noqa: E731
    for p in P_B:
        l07 = pd.read_csv(d / f"load_0.70_{_tag_p(p)}.csv.gz"); l085 = pd.read_csv(d / f"load_0.85_{_tag_p(p)}.csv.gz")
        s2, s3, ok = flex_ratios(b, l07, l085)
        ns2, ns3, _ = flex_ratios(one(nb), one(solved[(p, "load_0.70")]["nominal"]), one(solved[(p, "load_0.85")]["nominal"]))
        for lab, v, nv in (("S2 daily vs S1", s2, ns2[0]), ("S3 renewable vs S1", s3, ns3[0])):
            rows.append({"statement": lab, "p": p, "firing_0.70": solved[(p, "load_0.70")]["firing"], "firing_0.85": solved[(p, "load_0.85")]["firing"],
                         "T_wo_max_0.70_K": solved[(p, "load_0.70")]["T_wo_max_K"], "T_wo_max_0.85_K": solved[(p, "load_0.85")]["T_wo_max_K"],
                         "T_out_check_0.70": solved[(p, "load_0.70")]["T_out_K"], "nominal": float(nv), **_q(v), "fraction_le_1": float(np.mean(v <= 1.0)), "n_ok": int(ok.sum())})
        log(f"L p={p:<5g} f07 {rows[-1]['firing_0.70']:.4f} f085 {rows[-1]['firing_0.85']:.4f}  S2 nom {ns2[0]:.4f} med {np.median(s2):.4f} frac {np.mean(s2 <= 1):.3f}  "
            f"S3 nom {ns3[0]:.4f} med {np.median(s3):.4f} frac {np.mean(s3 <= 1):.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "flex_magnitude_vs_p_holdTout.csv", index=False)
    _dump("calc_L.json", {"N": n, "rows": rows, "S1_base_firing": pts["S1_base"]["specific_firing_factor"],
                          "published_firing_p0": {"load_0.70": pts["load_0.70"]["specific_firing_factor"], "load_0.85": pts["load_0.85"]["specific_firing_factor"]}})
    return df


def _tag_p(p: float) -> str:
    return f"p{p:g}".replace(".", "_")


# ---------------------------------------------------------------------------
# M. Alternative creep-uncertainty models (re-scoring the stored v3 Monte Carlo)
# ---------------------------------------------------------------------------
M_POINTS = ("RQ3_knee", "RQ3_min_life_iso_H2", "S5_y4_holdslip_end", "S6_SC2.5", "S6_SC3.5", "load_0.70", "load_0.85")
T_REF_SHEET_K = 1373.15      # 1100 degC: steam-reforming max operating temperature printed on the sheet (see report)
T_REF_NEAR_K = 1147.0
T_DF = 4


class _CurveFamily:
    """t_r(T, sigma; C) for one creep parametrisation, cached by rounded C."""

    def __init__(self, kind: str, T_ref: Optional[float] = None):
        from rdt import creep
        self.kind, self.T_ref = kind, T_ref
        self.c0 = creep.load_alloy(CFG.creep_alloy).minimum
        self.C0 = self.c0.C
        pts = pd.read_csv(C.ROOT / "data" / "creep_derived" / "centralloy_g_4852_rupture_curve.csv")
        g = pts[pts.curve == "minimum"].sort_values("lmp")
        self.lmp_pts, self.ls_pts = g.lmp.to_numpy(float), np.log10(g.stress_mpa.to_numpy(float))
        self.cache: Dict[float, object] = {}
        self.max_shift_dev = 0.0

    def curve(self, C_: float):
        from rdt import creep
        key = round(float(C_), 3)
        if key not in self.cache:
            if self.kind == "v3":
                c = creep.LarsonMillerCurve(C=key, scale=self.c0.scale, degree=self.c0.degree, points=[], alloy=self.c0.alloy, source=self.c0.source, curve_kind="minimum")
                c.coeffs = self.c0.coeffs.copy(); c.lmp_range = self.c0.lmp_range
            else:   # consistent: digitised (LMP, sigma) -> (sigma, t_r) at T_ref with C0 -> LMP with C -> refit P3
                s = self.c0.scale
                log_tr = self.lmp_pts / (s * self.T_ref) - self.C0
                lmp_new = s * self.T_ref * (key + log_tr)
                c = creep.LarsonMillerCurve(C=key, scale=s, degree=self.c0.degree, points=[], alloy=self.c0.alloy, source=self.c0.source, curve_kind="minimum")
                c.coeffs = np.polyfit(lmp_new, self.ls_pts, self.c0.degree); c.lmp_range = (float(lmp_new.min()), float(lmp_new.max()))
                # the refit must equal the base polynomial shifted by s T_ref (C - C0)
                d = s * self.T_ref * (key - self.C0)
                xx = np.linspace(*self.c0.lmp_range, 50)
                self.max_shift_dev = max(self.max_shift_dev, float(np.max(np.abs(np.polyval(c.coeffs, xx + d) - np.polyval(np.polyfit(self.lmp_pts, self.ls_pts, self.c0.degree), xx)))))
            self.cache[key] = c
        return self.cache[key]

    def log10_tr(self, T: np.ndarray, sig: np.ndarray, Cs: np.ndarray) -> np.ndarray:
        out = np.full(len(T), np.nan)
        for i, (t, s_, c_) in enumerate(zip(T, sig, Cs)):
            if np.isfinite(t) and np.isfinite(s_):
                out[i] = np.log10(float(self.curve(c_).time_to_rupture(t, s_)))
        return out


def _scatter(x_norm: np.ndarray, sigma_n: float, kind: str) -> np.ndarray:
    from scipy import stats
    if kind == "normal":
        return x_norm
    U = stats.norm.cdf(x_norm / sigma_n)
    s_t = (1.645 * sigma_n) / stats.t.ppf(0.95, T_DF)             # 5 % quantile = lower band (-1.645 sigma_n)
    return stats.t.ppf(U, T_DF) * s_t


def calc_M(log=print) -> pd.DataFrame:
    from rdt import uq
    d = Path(uq.OUT_DIR)
    mc = {n: pd.read_csv(d / f"{uq.MC_TAG}_{n}.csv.gz") for n in uq.POINT_NAMES}
    grp = {g: pd.read_csv(d / f"{uq.MC_TAG}_S1_base_group_{g}.csv.gz") for g in uq.GROUPS}
    spec = uq.build_spec(save=False); sig_n = float(spec["groups"]["creep"]["log10_tr_scatter"]["sigma"])
    summ = json.loads(Path(uq.SUMMARY_JSON).read_text())
    variants = {"M0 (v3)": ("v3", None, "normal"), f"M1 consistent C, T_ref {T_REF_SHEET_K:.0f} K": ("consistent", T_REF_SHEET_K, "normal"),
                f"M1 consistent C, T_ref {T_REF_NEAR_K:.0f} K": ("consistent", T_REF_NEAR_K, "normal"), "M2 Student-t (4 dof) scatter": ("v3", None, "t")}
    rows, info = [], {}
    for vname, (kind, Tref, sk) in variants.items():
        fam = _CurveFamily(kind, Tref)
        def rescore(df):
            df = df.copy()
            lt = fam.log10_tr(df.T_wo_max_K.to_numpy(float), df.sigma_hot_MPa.to_numpy(float), df.C_LM.to_numpy(float)) + _scatter(df.log10_tr_scatter.to_numpy(float), sig_n, sk)
            df["log10_t_r"] = lt; df["life_rate_per_h"] = 10.0 ** -lt; df["life_rate_per_kmol_H2"] = df.life_rate_per_h / df.H2_net_kmol_h
            return df
        M = {n: rescore(v) for n, v in mc.items()}; M["load_1.00"] = M["S1_base"]
        G = {g: rescore(v) for g, v in grp.items()}
        if vname == "M0 (v3)":
            ok = _ok(mc["S1_base"])
            dev = float(np.nanmax(np.abs(M["S1_base"].log10_t_r.to_numpy()[ok] - mc["S1_base"].log10_t_r.to_numpy()[ok])))
            info["M0_reproduction_max_abs_dev_log10tr"] = dev
        vs = uq.variance_shares(M["S1_base"], G, cols=("log10_t_r",))["log10_t_r"]
        iv = uq.intervals(M["S1_base"], M["S1_base"])
        rf = uq.robustness_fractions(M, M)["nominal_L_q"]
        r = {"variant": vname, "creep_share_log10tr": vs["shares"]["creep"], "kinetics_share": vs["shares"]["kinetics"], "heat_transfer_share": vs["shares"]["heat_transfer"],
             "measurement_share": vs["shares"]["measurement"], "S1_log10tr_median": iv["log10_t_r"]["median"], "S1_log10tr_p05": iv["log10_t_r"]["p05"], "S1_log10tr_p95": iv["log10_t_r"]["p95"]}
        for n in M_POINTS:
            q = uq.intervals(M[n], M["S1_base"])["life_rate_per_kmol_H2_rel_base"]
            r.update({f"{n}_median": q["median"], f"{n}_p05": q["p05"], f"{n}_p95": q["p95"]})
        r.update({f"frac_{k}": v for k, v in rf.items()})
        if kind == "consistent":
            r["refit_vs_shift_max_dev_log10sigma"] = fam.max_shift_dev
        rows.append(r)
        log(f"M {vname:34s} creep share {r['creep_share_log10tr']:.3f}  S1 log tr [{r['S1_log10tr_p05']:.2f},{r['S1_log10tr_p95']:.2f}]  knee {r['RQ3_knee_median']:.3f} "
            f"[{r['RQ3_knee_p05']:.3f},{r['RQ3_knee_p95']:.3f}]  S2<=S1 {rf['S2_life_per_kmol_le_S1']:.3f}")
    df = pd.DataFrame(rows)
    pub = summ["robustness_fractions"]["nominal_L_q"]
    info["M0_vs_published_fractions"] = {k: [float(df.loc[0, f"frac_{k}"]), float(v)] for k, v in pub.items()}
    info["M0_vs_published_creep_share"] = [float(df.loc[0, "creep_share_log10tr"]), float(summ["variance_shares_S1_base"]["log10_t_r"]["shares"]["creep"])]
    info["scatter"] = {"normal_sigma": sig_n, "lower_band_decades": 1.645 * sig_n, "t_df": T_DF, "t_scale": 1.645 * sig_n / float(__import__("scipy").stats.t.ppf(0.95, T_DF))}
    info["T_ref_K"] = {"sheet": T_REF_SHEET_K, "near_hot_spot": T_REF_NEAR_K}
    df.to_csv(OUT / "creep_uq_alternatives.csv", index=False)
    _dump("calc_M.json", {"rows": rows, **info})
    log(f"M M0 reproduction: {info['M0_reproduction_max_abs_dev_log10tr']:.2e}; fractions vs published {info['M0_vs_published_fractions']}")
    return df


# ---------------------------------------------------------------------------
# Report tables, package 3
# ---------------------------------------------------------------------------
def report_tables_pkg3() -> Dict[str, str]:
    t: Dict[str, str] = {}
    w = pd.read_csv(OUT / "walltemp_error_propagation.csv")
    for tab, key in (("Table 2", "I2"), ("Table 3", "I3")):
        g = w[w.table == tab].pivot(index="item", columns="variant", values="rel_life_cost")
        order = [x for x in (TABLE2_LABELS.values() if tab == "Table 2" else TABLE3_REGIMES) if x in g.index]
        cols = ["none (v3)", "a: dT = -18 K", "a: dT = -12 K", "a: dT = +12 K", "a: dT = +18 K", "b: w*e_u(load), sign as is", "b: -w*e_u(load), reversed", "b: 1*e_u(load), w = 1"]
        g = g.loc[order, cols].reset_index().rename(columns={"item": "строка"})
        t[key] = _md(g, {c: ".3f" for c in cols})
    ij = json.loads((OUT / "calc_I.json").read_text())
    t["I_nd"] = _md(pd.DataFrame([{"вариант": k, "недоминируемых из 198": v["n_nondominated"], "потеряли статус": v["lost"], "приобрели": v["gained"],
                                   "knee — тот же член, что в v3": "да" if v["knee_same_as_v3"] else "нет", "knee: load / S/C / EA / firing":
                                   f"{v['knee']['feed_per_tube_fraction']:.3f} / {v['knee']['steam_to_carbon']:.2f} / {v['knee']['excess_air']:.1f} / {v['knee']['specific_firing_factor']:.3f}",
                                   "rel life knee": v["knee_rel_life"]} for k, v in ij["nondominance_and_knee"].items()]), {"rel life knee": ".3f"})
    t["I_e"] = _md(pd.DataFrame([{"load": k, "e_u (K)": v, "dT_hot = w e_u (K)": W_HOT * v} for k, v in ij["e_u_line"]["e_u_at_load"].items()]), {"e_u (K)": "+.1f", "dT_hot = w e_u (K)": "+.1f"})
    for name, key, free in (("flame_length_sweep.csv", "J", ("F_gt", "f_htg")), ("htc_sweep.csv", "K", ("F_gt", "L_q"))):
        d = pd.read_csv(OUT / name)
        fit = pd.DataFrame({"фикс.": d.value, **{k: [f"{a:.4f} [{b}]" for a, b in zip(d[f"fit_{k}"], d[f"ci95_{k}"])] for k in free},
                            "χ²/dof": [f"{a:.1f}/{b} = {c:.2f}" for a, b, c in zip(d.chi2, d.dof, d.reduced_chi2)], "на границе": d.any_at_bound.map({True: "да", False: "нет"}),
                            "A ΔT_w,u/l": [f"{a:+.1f}/{b:+.1f}" for a, b in zip(d.Plant_A_dTw_u, d.Plant_A_dTw_l)],
                            "B ΔT_w,u/l": [f"{a:+.1f}/{b:+.1f}" for a, b in zip(d.Plant_B_dTw_u, d.Plant_B_dTw_l)],
                            "C2 ΔT_w,u/l": [f"{a:+.1f}/{b:+.1f}" for a, b in zip(d.Plant_C2_dTw_u, d.Plant_C2_dTw_l)],
                            "C1 восст. ΔT_out/ΔT_w,u/ΔT_w,l": [f"{a:+.1f}/{b:+.1f}/{c:+.1f}" for a, b, c in zip(d.C1_out_T_K_err, d.C1_TWT_upper_K_err, d.C1_TWT_lower_K_err)]})
        t[key + "1"] = _md(fit, {"фикс.": "g"})
        rat = pd.DataFrame({"фикс.": d.value, "T_wo,max база (K)": d.base_T_wo_max_K, "z/L": d.base_z_frac_T_wo_max, "T_out база (K)": d.base_T_out_K,
                            **{n: d[f"rel_{n}"] for n in TABLE3_REGIMES}, "S5 slip y4/S1": d.rel_S5_slip_y4, "S/C 2.5/S1": d["rel_S6_SC2.5"], "S/C 3.5/S1": d["rel_S6_SC3.5"],
                            "firing S1": d.S1_firing_phys, "firing на границе": d.n_firing_saturated})
        t[key + "2"] = _md(rat, {"фикс.": "g", "T_wo,max база (K)": ".1f", "z/L": ".3f", "T_out база (K)": ".1f", **{n: ".3f" for n in TABLE3_REGIMES},
                                 "S5 slip y4/S1": ".3f", "S/C 2.5/S1": ".3f", "S/C 3.5/S1": ".3f", "firing S1": ".4f", "firing на границе": ".0f"})
    kj = json.loads((OUT / "calc_K.json").read_text())
    t["K_corr"] = _md(pd.DataFrame([{"пара": k, "корреляция": v} for k, v in kj["correlations_adopted_fit"].items()]), {"корреляция": "+.3f"})
    ch = pd.read_csv(OUT / "htc_chain_terms.csv")
    piv = ch.pivot(index="term", columns="case", values="value").reindex(ch.term.drop_duplicates())
    meta = ch.drop_duplicates("term").set_index("term")
    tab = pd.DataFrame({"член": piv.index, "Xu-Froment 1989": [("—" if pd.isna(v) else f"{v:.5g}") for v in piv["Xu-Froment 1989"]],
                        "Plant A": [("—" if pd.isna(v) else f"{v:.5g}") for v in piv["Plant A"]], "ед.": meta.loc[piv.index, "unit"].fillna("").to_numpy(),
                        "формула / определение": [str(ch[(ch.term == k)].formula_or_definition.iloc[-1] if k.startswith(("alpha_i, back", "alpha_i, cal", "d_p")) and False else meta.loc[k, "formula_or_definition"]) for k in piv.index],
                        "код": meta.loc[piv.index, "code_ref"].fillna("").to_numpy()})
    t["K_chain"] = _md(tab, {})
    t["K_chain_notes"] = "\n".join(f"- {r.case}, {r.term}: {r.formula_or_definition}" + (f" ({r.code_ref})" if isinstance(r.code_ref, str) and r.code_ref else "")
                                   for r in ch.itertuples() if r.term in ("d_p equivalent", "alpha_i, back-calculated", "alpha_i, calibrated"))
    l = pd.read_csv(OUT / "flex_magnitude_vs_p_holdTout.csv")
    t["L"] = _md(l[["statement", "p", "firing_0.70", "firing_0.85", "T_wo_max_0.70_K", "T_wo_max_0.85_K", "nominal", "median", "p05", "p95", "fraction_le_1"]]
                 .rename(columns={"p05": "5 %", "p95": "95 %", "fraction_le_1": "доля ≤ 1", "T_wo_max_0.70_K": "T_wo,max 0.70 (K)", "T_wo_max_0.85_K": "T_wo,max 0.85 (K)"}),
                 {"p": "g", "firing_0.70": ".4f", "firing_0.85": ".4f", "T_wo,max 0.70 (K)": ".1f", "T_wo,max 0.85 (K)": ".1f", "nominal": ".4f", "median": ".4f", "5 %": ".4f", "95 %": ".4f", "доля ≤ 1": ".3f"})
    b = pd.read_csv(OUT / "flex_magnitude_vs_p.csv")
    cmp_ = l[["statement", "p", "median", "fraction_le_1"]].merge(b[["statement", "p", "median", "fraction_le_1"]], on=["statement", "p"], suffixes=(" (T_out удерживается)", " (Table 6, firing p=0)"))
    t["L_cmp"] = _md(cmp_, {"p": "g", "median (T_out удерживается)": ".4f", "fraction_le_1 (T_out удерживается)": ".3f", "median (Table 6, firing p=0)": ".4f", "fraction_le_1 (Table 6, firing p=0)": ".3f"})
    m = pd.read_csv(OUT / "creep_uq_alternatives.csv")
    m1 = pd.DataFrame({"вариант": m.variant, "доля creep (log10 t_r)": m.creep_share_log10tr, "kinetics": m.kinetics_share, "heat transfer": m.heat_transfer_share,
                       "measurement": m.measurement_share, "S1 log10 t_r: медиана [5, 95 %]": [f"{a:.2f} [{b:.2f}, {c:.2f}]" for a, b, c in zip(m.S1_log10tr_median, m.S1_log10tr_p05, m.S1_log10tr_p95)]})
    t["M1"] = _md(m1, {c: ".3f" for c in ("доля creep (log10 t_r)", "kinetics", "heat transfer", "measurement")})
    labs = {"RQ3_knee": "knee", "RQ3_min_life_iso_H2": "min life iso-H2", "S5_y4_holdslip_end": "S5 y4", "S6_SC2.5": "S/C 2.5", "S6_SC3.5": "S/C 3.5", "load_0.70": "load 0.70", "load_0.85": "load 0.85"}
    m2 = pd.DataFrame({"вариант": m.variant, **{lab: [f"{a:.3f} [{b:.3f}, {c:.3f}]" for a, b, c in zip(m[f"{n}_median"], m[f"{n}_p05"], m[f"{n}_p95"])] for n, lab in labs.items()}})
    t["M2"] = _md(m2, {})
    fr = [c for c in m.columns if c.startswith("frac_")]
    t["M3"] = _md(m[["variant"] + fr].rename(columns={c: c[5:] for c in fr}), {c[5:]: ".3f" for c in fr})
    return t


# ---------------------------------------------------------------------------
# Report tables, package 2
# ---------------------------------------------------------------------------
def report_tables_pkg2() -> Dict[str, str]:
    t: Dict[str, str] = {}
    e = pd.read_csv(OUT / "ageing_histories.csv")
    j = lambda s, f: " / ".join(format(v, f) for v in s)  # noqa: E731
    rows = []
    for (var, hist, ctl), g in e.groupby(["variant", "history", "control"], sort=False):
        g = g.sort_values("year")
        rows.append({"история": hist, "метод": "surrogate, почасово" if var == "surrogate_hourly" else "физика, 6 бинов/год", "управление": ctl,
                     "rel life cost по годам": j(g.rel_life_per_kmol_H2_vs_S1, ".3f"), "T̄_wo,max по годам (K)": j(g.T_wo_max_mean_K, ".1f"),
                     "max T_wo,max по годам (K)": j(g.T_wo_max_max_K, ".1f"),
                     "кампания (ΣD/ΣH2 vs S1)": float((g.rel_life_per_kmol_H2_vs_S1 * g.annual_H2_kmol).sum() / g.annual_H2_kmol.sum())})
    t["E"] = _md(pd.DataFrame(rows), {"кампания (ΣD/ΣH2 vs S1)": ".3f"})
    ej = json.loads((OUT / "calc_E.json").read_text())
    t["E_chk"] = _md(pd.DataFrame([{"год": k, "Table 2": v["summary_v3"], "E0 surrogate": v["surrogate_hourly"], "откл. surrogate": f"{v['surrogate_rel_err']:+.1e}",
                                    "E0 физика (бины)": v["physics_binned"], "откл. физика, %": 100 * v["physics_rel_err"]} for k, v in ej["check_E0_vs_table2"].items()]),
                     {"Table 2": ".4f", "E0 surrogate": ".4f", "E0 физика (бины)": ".4f", "откл. физика, %": "+.2f"})
    f = pd.read_csv(OUT / "pareto_constrained.csv")
    cols = {"variant": "вар.", "constraints": "ограничения", "regime": "режим", "feed_per_tube_fraction": "load", "steam_to_carbon": "S/C", "inlet_T": "T_in (K)",
            "excess_air": "EA (%)", "specific_firing_factor": "firing", "H2_kmol_h": "H2 (kmol/h)", "fuel_MJ_per_kmol_H2": "fuel (MJ/kmol)",
            "rel_life_cost": "rel life", "T_wo_max_K": "T_wo,max (K)", "CH4_slip_pct": "slip (%)", "H2_rel_base_pct": "ΔH2 %", "fuel_rel_base_pct": "Δfuel %",
            "n_dominating_base": "домин. базу", "n_iso_band": "в iso-полосе", "n_feasible_phys": "допуст."}
    ft = f[list(cols)].rename(columns=cols)
    t["F"] = _md(ft, {"load": ".3f", "S/C": ".3f", "T_in (K)": ".1f", "EA (%)": ".1f", "firing": ".3f", "H2 (kmol/h)": ".2f", "fuel (MJ/kmol)": ".1f",
                      "rel life": ".3f", "T_wo,max (K)": ".1f", "slip (%)": ".2f", "ΔH2 %": "+.1f", "Δfuel %": "+.1f", "домин. базу": "d", "в iso-полосе": "d", "допуст.": "d"})
    g = pd.read_csv(OUT / "stress_multiplier.csv")
    gc = {"k": "k", "sigma_eff_base_MPa": "σ_eff база (MPa)", "LMP_base": "LMP база", "t_r_base_y": "t_r база (лет)",
          "slope_dlog10tr_dT_base_per_K": "d log10 t_r/dT база (1/K)", "life_halving_dT_base_K": "ΔT удвоения скорости (K)",
          "sigma_eff_range_MPa": "диапазон kσ (MPa)", "inside_digitised_range": "внутри кривой"}
    t["G1"] = _md(g[list(gc)].rename(columns=gc).assign(**{"внутри кривой": lambda d: d["внутри кривой"].map({True: "да", False: "НЕТ"})}),
                  {"k": "g", "σ_eff база (MPa)": ".2f", "LMP база": ".3f", "t_r база (лет)": ".4g", "d log10 t_r/dT база (1/K)": ".4f", "ΔT удвоения скорости (K)": ".1f"})
    g2 = {"k": "k", **{f"{n}_vs_base": n for n in TABLE3_REGIMES}, "S5_hold_CH4_slip_ageing_y4_vs_S1": "S5 slip y4 / S1",
          "S6_SC2.5_vs_S1": "S/C 2.5 / S1", "S6_SC3.5_vs_S1": "S/C 3.5 / S1"}
    t["G2"] = _md(g[list(g2)].rename(columns=g2), {c: ".3f" for c in g2.values() if c != "k"} | {"k": "g"})
    h = pd.read_csv(OUT / "c1_atom_balances.csv", index_col=0)
    hb = h[["feed_kmol_h", "out_kmol_h", "C_closure_out_over_in", "H_closure_out_over_in", "feed_from_C_balance_kmol_h", "feed_from_H_balance_kmol_h"]].reset_index()
    hb.columns = ["случай", "подача (kmol/h)", "выход (kmol/h)", "C вых/вх", "H вых/вх", "подача из баланса C", "подача из баланса H"]
    t["H1"] = _md(hb, {"подача (kmol/h)": ".1f", "выход (kmol/h)": ".1f", "C вых/вх": ".4f", "H вых/вх": ".4f", "подача из баланса C": ".1f", "подача из баланса H": ".1f"})
    p = pd.read_csv(OUT / "c1_heldout.csv")
    hc = {"case": "случай", "kind": "роль", "feed_kmol_h": "подача", "excess_air_pct": "EA (%)", "out_T_K_err": "ΔT_out", "TWT_upper_K_err": "ΔT_w,u", "TWT_lower_K_err": "ΔT_w,l",
          "flue_out_T_K_err": "ΔT_fg", "out_P_Pa_err": "Δp (kPa)", "out_x_H2_wet_molpct_err": "ΔH2", "out_x_CH4_wet_molpct_err": "ΔCH4",
          "out_x_CO_wet_molpct_err": "ΔCO", "out_x_CO2_wet_molpct_err": "ΔCO2"}
    t["H2"] = _md(p[list(hc)].rename(columns=hc), {"подача": ".1f", "EA (%)": ".2f", **{c: "+.2f" if c.startswith("ΔC") or c == "ΔH2" else "+.1f" for c in list(hc.values())[4:]}})
    hj = json.loads((OUT / "calc_H.json").read_text()); a3, a4 = hj["adopted_fit_3case"], hj["fit_4case"]
    pr = []
    for k in FREE3:
        pr.append({"параметр": k, "A, B, C2 (принят)": f"{a3['params'][k]:.4f} [{a3['ci95'][k][0]:.3f}, {a3['ci95'][k][1]:.3f}]",
                   "A, B, C1 восст., C2": f"{a4['params'][k]:.4f} [{a4['ci95'][k][0]:.3f}, {a4['ci95'][k][1]:.3f}]"})
    pr.append({"параметр": "χ² / dof", "A, B, C2 (принят)": f"{a3['chi2']:.2f} / {a3['dof']} = {a3['reduced_chi2']:.3f}",
               "A, B, C1 восст., C2": f"{a4['chi2']:.2f} / {a4['dof']} = {a4['reduced_chi2']:.3f}"})
    t["H3"] = _md(pd.DataFrame(pr), {})
    lr = []
    for c in ("Plant_A", "Plant_B", "Plant_C1", "Plant_C2"):
        r4 = a4["loco"][c]["held_out_errors"]; r3 = a3["loco"].get(c)
        lr.append({"отложен": c, "LOCO 3 случая ΔT_w,u": f"{r3['TWT_upper_K']:+.1f}" if r3 else "—", "LOCO 3 случая ΔT_w,l": f"{r3['TWT_lower_K']:+.1f}" if r3 else "—",
                   "LOCO 4 случая ΔT_w,u": f"{r4['TWT_upper_K']:+.1f}", "LOCO 4 случая ΔT_w,l": f"{r4['TWT_lower_K']:+.1f}",
                   "параметры фолда (F_gt, L_q, f_htg)": ", ".join(f"{a4['loco'][c]['params'][k]:.3f}" for k in FREE3)})
    t["H4"] = _md(pd.DataFrame(lr), {})
    return t


# ---------------------------------------------------------------------------
def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--calc", nargs="+", default=["A", "B", "C", "D"], choices=list("ABCDEFGHIJKLM"))
    ap.add_argument("--n-d", type=int, default=512)
    ap.add_argument("--n-jobs", type=int, default=-1)
    a = ap.parse_args(argv)
    setup(); OUT.mkdir(parents=True, exist_ok=True)
    log = lambda m: print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)  # noqa: E731
    dev = check_design(); log(f"stored mc_v3 design vs fresh uq.sample: max abs diff {dev:.3g}"); _dump("design_check.json", {"max_abs_diff": dev})
    for c in a.calc:
        {"A": lambda: calc_A(log=log), "B": lambda: calc_B(log=log), "C": lambda: calc_C(log=log),
         "D": lambda: calc_D(a.n_d, a.n_jobs, log=log), "E": lambda: calc_E(log=log), "F": lambda: calc_F(log=log),
         "G": lambda: calc_G(log=log), "H": lambda: calc_H(log=log), "I": lambda: calc_I(log=log), "J": lambda: calc_J(log=log),
         "K": lambda: calc_K(log=log), "L": lambda: calc_L(a.n_d, a.n_jobs, log=log), "M": lambda: calc_M(log=log)}[c]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
