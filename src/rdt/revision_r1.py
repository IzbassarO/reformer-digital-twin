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


# ---------------------------------------------------------------------------
def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--calc", nargs="+", default=["A", "B", "C", "D"], choices=list("ABCD"))
    ap.add_argument("--n-d", type=int, default=512)
    ap.add_argument("--n-jobs", type=int, default=-1)
    a = ap.parse_args(argv)
    setup(); OUT.mkdir(parents=True, exist_ok=True)
    log = lambda m: print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)  # noqa: E731
    dev = check_design(); log(f"stored mc_v3 design vs fresh uq.sample: max abs diff {dev:.3g}"); _dump("design_check.json", {"max_abs_diff": dev})
    for c in a.calc:
        {"A": lambda: calc_A(log=log), "B": lambda: calc_B(log=log), "C": lambda: calc_C(log=log),
         "D": lambda: calc_D(a.n_d, a.n_jobs, log=log)}[c]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
