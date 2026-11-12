"""Revision R1 of HE-D-26-14265, calculation package 1.

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
# Report tables (the prose of REPORT_pkg1.md is written around these)
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
    ap.add_argument("--calc", nargs="+", default=["A", "B", "C", "D"], choices=list("ABCDEFGH"))
    ap.add_argument("--n-d", type=int, default=512)
    ap.add_argument("--n-jobs", type=int, default=-1)
    a = ap.parse_args(argv)
    setup(); OUT.mkdir(parents=True, exist_ok=True)
    log = lambda m: print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)  # noqa: E731
    dev = check_design(); log(f"stored mc_v3 design vs fresh uq.sample: max abs diff {dev:.3g}"); _dump("design_check.json", {"max_abs_diff": dev})
    for c in a.calc:
        {"A": lambda: calc_A(log=log), "B": lambda: calc_B(log=log), "C": lambda: calc_C(log=log),
         "D": lambda: calc_D(a.n_d, a.n_jobs, log=log), "E": lambda: calc_E(log=log), "F": lambda: calc_F(log=log),
         "G": lambda: calc_G(log=log), "H": lambda: calc_H(log=log)}[c]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
