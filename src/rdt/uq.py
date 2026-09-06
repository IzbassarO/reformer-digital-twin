"""RQ4: uncertainty propagation through the physics model and robustness of the RQ2/RQ3 conclusions.

Uncertainty groups (see ``data/uq/uncertainty_spec.yaml``): kinetics (Xu-Froment 95 % intervals, effectiveness
multiplier, activity), heat transfer / furnace (calibration covariance of F_gt, L_q, f_htg; structural variant
``L_q_eff = L_q * load**0.5``), creep (heat-to-heat scatter, Larson-Miller constant, wall thickness, tube
conductivity) and measurement (tube-wall temperature offset). A scrambled Sobol design (N = 4096) is propagated
through the coupled physics model at several fixed operating points (control is not re-solved per sample); the
base point is evaluated with the same samples so that life ratios are paired.
"""

from __future__ import annotations

import json
import math
import time
from dataclasses import replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import yaml
from scipy import stats
from scipy.stats import qmc

from rdt import creep
from rdt import kinetics as kin
from rdt import latham_cases as lc
from rdt import operate as op
from rdt import scenarios as sc

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "data" / "uq"
SPEC_YAML = OUT_DIR / "uncertainty_spec.yaml"
N_MC = 4096
SEED = 0

GROUPS = {
    "kinetics": ["E1", "E2", "E3", "logk1", "logk2", "logk3", "logK_CO", "logK_H2", "logK_CH4", "logK_H2O", "log_eta_mult", "d_activity"],
    "heat_transfer": ["F_gt", "L_q", "f_htg"],
    "creep": ["log10_tr_scatter", "C_LM", "wall_thickness_m", "lambda_tube_factor"],
    "measurement": ["dT_wo_meas_K"],
}
PARAMS = [p for g in GROUPS.values() for p in g]


# ---------------------------------------------------------------------------
# Specification
# ---------------------------------------------------------------------------
def build_spec(save: bool = True) -> Dict[str, object]:
    kp = kin.default_params(); ci = kp.ci95
    doc = yaml.safe_load(lc.FIT_YAML.read_text()); cal = doc["calibration_alpha_top_fixed"]["best_fit"]
    cov_meta = cal["covariance"]; names = list(cal["params"])
    se = np.array([cov_meta["se"][k] for k in names]); corr = np.array([[cov_meta["corr_table"][a][b] for b in names] for a in names])
    cov = corr * np.outer(se, se)
    spec = {"sampling": {"design": "scipy.stats.qmc.Sobol, scrambled, seed 0", "N": N_MC, "note": "N = 4096 (power of 2) instead of 4000 for Sobol balance"},
            "control": "operating points are evaluated at their nominal (solved) firing; control loops are not re-solved per sample",
            "groups": {}}
    kin_spec = {}
    for i in (1, 2, 3):
        v, (ul, ll) = kp.E[i], ci[f"E{i}"]
        kin_spec[f"E{i}"] = {"dist": "normal", "mean": v, "sigma": (ul - ll) / 2 / 1.96, "unit": "kJ/mol", "source": "Xu & Froment 1989 Table 5 95 % CI"}
        v, (ul, ll) = kp.k_ref[i], ci[f"k{i}_648"]
        kin_spec[f"logk{i}"] = {"dist": "normal", "mean": math.log10(v), "sigma": (math.log10(ul) - math.log10(ll)) / 2 / 1.96, "unit": "log10 of k_i at 648 K", "source": "Table 5 95 % CI, normal on log k"}
    for j, key in (("CO", "K_CO_648"), ("H2", "K_H2_648"), ("CH4", "K_CH4_823"), ("H2O", "K_H2O_823")):
        v, (ul, ll) = kp.K_ref[j], ci[key]
        if j == "H2O":
            ll = 2 * v - ul   # symmetric bound (printed LL 0.0317 suspected misprint)
        kin_spec[f"logK_{j}"] = {"dist": "normal", "mean": math.log10(v), "sigma": (math.log10(ul) - math.log10(max(ll, 1e-6))) / 2 / 1.96, "unit": f"log10 of K_{j} at T_ref",
                                 "source": "Table 5 95 % CI, normal on log K" + (" (symmetric K_H2O bound)" if j == "H2O" else "")}
    kin_spec["log_eta_mult"] = {"dist": "uniform", "low": math.log10(0.5), "high": math.log10(2.0), "unit": "log10 multiplier on eta (0.1 / 0.05 top)", "source": "log-uniform [0.5, 2], assumption"}
    kin_spec["d_activity"] = {"dist": "uniform", "low": -0.03, "high": 0.03, "unit": "additive on activity 0.20", "source": "assumption, +/-0.03 around the calibrated value"}
    spec["groups"]["kinetics"] = kin_spec
    spec["groups"]["heat_transfer"] = {"F_gt,L_q,f_htg": {"dist": "multivariate_normal", "names": names, "mean": [cal["params"][k] for k in names], "cov": cov.tolist(),
                                                           "source": "calibration_alpha_top_fixed Jacobian covariance (se x corr)"},
                                       "structural_variant": {"L_q_eff": "L_q * load**0.5", "note": "flame length shrinks with load; ASSUMPTION used only for the robustness check"}}
    spec["groups"]["creep"] = {"log10_tr_scatter": {"dist": "normal", "mean": 0.0, "sigma": 0.3, "unit": "log10 multiplier on t_r", "source": "heat-to-heat scatter, factor ~2, placeholder until NIMS heats"},
                               "C_LM": {"dist": "uniform", "low": 21.96, "high": 23.96, "source": "Larson-Miller constant 22.96 +/- 1.0; master-curve points refitted with each C"},
                               "wall_thickness_m": {"dist": "uniform", "low": 0.012, "high": 0.018, "source": "wall thickness not given by Latham; 12-18 mm"},
                               "lambda_tube_factor": {"dist": "uniform", "low": 0.85, "high": 1.15, "source": "tube conductivity +/-15 % around 29.6 W/(m K)"}}
    spec["groups"]["measurement"] = {"dT_wo_meas_K": {"dist": "uniform", "low": -10.0, "high": 10.0, "unit": "K, additive offset on T_wo,max before creep evaluation", "source": "pyrometer / background-correction uncertainty"}}
    spec["_cov_names"] = names; spec["_cov"] = cov.tolist(); spec["_mean_ht"] = [cal["params"][k] for k in names]
    if save:
        OUT_DIR.mkdir(exist_ok=True)
        yaml.safe_dump({k: v for k, v in spec.items() if not k.startswith("_")}, open(SPEC_YAML, "w"), sort_keys=False, width=120)
    return spec


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------
def sample(spec: Dict[str, object], n: int = N_MC, seed: int = SEED, active_groups: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Sobol design mapped to the marginal distributions; inactive groups are held at nominal."""
    d = len(PARAMS)
    U = qmc.Sobol(d=d, scramble=True, seed=seed).random(n)
    X = pd.DataFrame(index=range(n))
    ks = spec["groups"]["kinetics"]
    col = 0
    for p in GROUPS["kinetics"]:
        s = ks[p]
        X[p] = stats.norm.ppf(U[:, col], s["mean"], s["sigma"]) if s["dist"] == "normal" else s["low"] + (s["high"] - s["low"]) * U[:, col]
        col += 1
    # heat transfer: MVN via Cholesky of the Sobol normals
    mean = np.array(spec["_mean_ht"]); cov = np.array(spec["_cov"]); Lc = np.linalg.cholesky(cov)
    Z = stats.norm.ppf(U[:, col:col + 3]); col += 3
    H = mean + Z @ Lc.T
    for i, p in enumerate(spec["_cov_names"]):
        X[p] = H[:, i]
    cs = spec["groups"]["creep"]
    for p in GROUPS["creep"]:
        s = cs[p]; X[p] = stats.norm.ppf(U[:, col], s["mean"], s["sigma"]) if s["dist"] == "normal" else s["low"] + (s["high"] - s["low"]) * U[:, col]; col += 1
    ms = spec["groups"]["measurement"]["dT_wo_meas_K"]; X["dT_wo_meas_K"] = ms["low"] + (ms["high"] - ms["low"]) * U[:, col]
    nom = nominal(spec)
    if active_groups is not None:
        for g, ps in GROUPS.items():
            if g not in active_groups:
                for p in ps:
                    X[p] = nom[p]
    return X


def nominal(spec: Dict[str, object]) -> Dict[str, float]:
    kp = kin.default_params()
    nom = {f"E{i}": kp.E[i] for i in (1, 2, 3)}
    nom.update({f"logk{i}": math.log10(kp.k_ref[i]) for i in (1, 2, 3)})
    nom.update({f"logK_{j}": math.log10(kp.K_ref[j]) for j in ("CO", "H2", "CH4", "H2O")})
    nom.update({"log_eta_mult": 0.0, "d_activity": 0.0})
    nom.update(dict(zip(spec["_cov_names"], spec["_mean_ht"])))
    nom.update({"log10_tr_scatter": 0.0, "C_LM": 22.96, "wall_thickness_m": 0.015, "lambda_tube_factor": 1.0, "dT_wo_meas_K": 0.0})
    return nom


# ---------------------------------------------------------------------------
# Evaluation of one (operating point, parameter sample)
# ---------------------------------------------------------------------------
def kinetic_params_from(theta: Dict[str, float]) -> kin.XuFromentParams:
    kp = kin.default_params()
    return replace(kp, E={i: float(theta[f"E{i}"]) for i in (1, 2, 3)}, k_ref={i: 10 ** float(theta[f"logk{i}"]) for i in (1, 2, 3)},
                   K_ref={j: 10 ** float(theta[f"logK_{j}"]) for j in ("CO", "H2", "CH4", "H2O")})


def params_from(theta: Dict[str, float], base_params: lc.LathamParams, load: float, lq_variant: bool = False) -> lc.LathamParams:
    m = 10 ** float(theta["log_eta_mult"])
    L_q = float(theta["L_q"]) * (load ** 0.5 if lq_variant else 1.0)
    return replace(base_params, F_gt=float(theta["F_gt"]), L_q=L_q, f_htg=float(theta["f_htg"]), eta=base_params.eta * m, eta_top=base_params.eta_top * m,
                   wall_thickness_m=float(theta["wall_thickness_m"]), lambda_tube=29.6 * float(theta["lambda_tube_factor"]), kinetic_params=kinetic_params_from(theta))


_CURVE_CACHE: Dict[float, creep.LarsonMillerCurve] = {}


def curve_with_C(C: float) -> creep.LarsonMillerCurve:
    key = round(float(C), 3)
    if key not in _CURVE_CACHE:
        c0 = creep.LarsonMillerCurve.from_yaml()
        _CURVE_CACHE[key] = creep.LarsonMillerCurve(C=key, scale=c0.scale, degree=c0.degree, points=c0.points, alloy=c0.alloy, curve_kind=c0.curve_kind).fit()
    return _CURVE_CACHE[key]


def evaluate_one(point: Dict[str, float], theta: Dict[str, float], lq_variant: bool = False) -> Dict[str, float]:
    base = op.base_case()
    x = dict(point); x["catalyst_activity"] = float(x["catalyst_activity"]) + float(theta["d_activity"])
    p = params_from(theta, base.params, float(x["feed_per_tube_fraction"]), lq_variant)
    b2 = replace(base, params=p)
    r = op.run_case(x, b2, creep.LarsonMillerCurve.from_yaml())
    if not r["converged"]:
        return {"converged": False}
    T_wo = r["T_wo_max_K"] + float(theta["dT_wo_meas_K"])
    curve = curve_with_C(theta["C_LM"])
    t_r = float(curve.time_to_rupture(T_wo, r["sigma_hot_MPa"])) * 10 ** float(theta["log10_tr_scatter"])
    return {"converged": True, "T_wo_max_K": T_wo, "T_wo_max_model_K": r["T_wo_max_K"], "T_out_K": r["T_out_K"], "CH4_slip_dry_pct": r["CH4_slip_dry_pct"],
            "H2_net_kmol_h": r["H2_net_kmol_h"], "sigma_hot_MPa": r["sigma_hot_MPa"], "log10_t_r": math.log10(t_r), "life_rate_per_h": 1.0 / t_r,
            "life_rate_per_kmol_H2": 1.0 / t_r / r["H2_net_kmol_h"]}


def _worker(point, theta, lq_variant):
    return evaluate_one(point, theta, lq_variant)


def run_point(name: str, point: Dict[str, float], X: pd.DataFrame, lq_variant: bool = False, n_jobs: int = -1,
              out_dir: Path = OUT_DIR, tag: str = "mc_v1") -> pd.DataFrame:
    """Evaluate all samples for one operating point (resumable: cached to CSV)."""
    from joblib import Parallel, delayed
    path = out_dir / f"{tag}_{name}{'_lqvar' if lq_variant else ''}.csv.gz"
    if path.exists():
        return pd.read_csv(path)
    t0 = time.perf_counter()
    rows = Parallel(n_jobs=n_jobs)(delayed(_worker)(point, th, lq_variant) for th in X.to_dict(orient="records"))
    df = pd.concat([X.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    df["point"] = name; df["lq_variant"] = lq_variant; df["sample"] = np.arange(len(df)); df.attrs["wall_time_s"] = time.perf_counter() - t0
    df.to_csv(path, index=False, compression="gzip")
    return df


# ---------------------------------------------------------------------------
# Operating points
# ---------------------------------------------------------------------------
def operating_points() -> Dict[str, Dict[str, float]]:
    base = op.base_case(); b = base.inputs(); b["excess_air"] = sc.EXCESS_AIR_SCEN
    summ = pd.read_csv(sc.OUT_DIR / "summary_v1.csv").set_index("scenario")
    meta = json.loads((sc.OUT_DIR / "summary_v1_meta.json").read_text())
    reg = pd.read_csv(lc.ROOT / "data" / "optimization" / "regimes_v1.csv", index_col=0)
    pts = {}
    pts["S1_base"] = {**b, "specific_firing_factor": float(summ.loc["S1_steady", "firing_mean"])}
    for r in ("knee", "min_life_iso_H2"):
        pts[f"RQ3_{r}"] = {k: float(reg.loc[r, k]) for k in op.INPUT_NAMES}
    y4 = meta["S5_campaigns"]["hold_CH4_slip_ageing"][-1]
    pts["S5_y4_holdslip_end"] = {**b, "catalyst_activity": float(y4["activity_end"]), "specific_firing_factor": float(y4["firing_end"])}
    pts["S6_SC2.5"] = {**b, "steam_to_carbon": 2.5, "specific_firing_factor": float(summ.loc["S6_SC2.5", "firing_mean"])}
    pts["S6_SC3.5"] = {**b, "steam_to_carbon": 3.5, "specific_firing_factor": float(summ.loc["S6_SC3.5", "firing_mean"])}
    return pts


def flexibility_points() -> Dict[str, Dict[str, float]]:
    """Representative load states for S2 (0.70 night) and S3 (quantiles 0.70/0.85/1.00) with hold_T_out firing."""
    tw = sc.TwinSurrogate(); base = op.base_case(); b = base.inputs(); b["excess_air"] = sc.EXCESS_AIR_SCEN
    pts = {}
    for load in (0.70, 0.85, 1.00):
        X = pd.DataFrame([{**b, "feed_per_tube_fraction": load}])
        X["specific_firing_factor"] = tw.solve_firing(X, "hold_T_out", sc.T_OUT_TARGET)
        pts[f"load_{load:.2f}"] = {k: float(X[k].iloc[0]) for k in op.INPUT_NAMES}
    return pts


# ---------------------------------------------------------------------------
# Analyses
# ---------------------------------------------------------------------------
def intervals(df: pd.DataFrame, base_df: pd.DataFrame, cols=("T_wo_max_K", "T_out_K", "CH4_slip_dry_pct", "log10_t_r")) -> Dict[str, object]:
    ok = (df.converged == True) & (base_df.converged == True)  # noqa: E712
    out = {"n": int(ok.sum())}
    for c in cols:
        v = df.loc[ok, c].to_numpy(float)
        out[c] = {"median": float(np.median(v)), "p05": float(np.percentile(v, 5)), "p95": float(np.percentile(v, 95)), "width_90": float(np.percentile(v, 95) - np.percentile(v, 5)), "std": float(v.std())}
    ratio = df.loc[ok, "life_rate_per_h"].to_numpy() / base_df.loc[ok, "life_rate_per_h"].to_numpy()
    ratio_kmol = df.loc[ok, "life_rate_per_kmol_H2"].to_numpy() / base_df.loc[ok, "life_rate_per_kmol_H2"].to_numpy()
    out["life_rate_rel_base"] = {"median": float(np.median(ratio)), "p05": float(np.percentile(ratio, 5)), "p95": float(np.percentile(ratio, 95))}
    out["life_rate_per_kmol_H2_rel_base"] = {"median": float(np.median(ratio_kmol)), "p05": float(np.percentile(ratio_kmol, 5)), "p95": float(np.percentile(ratio_kmol, 95))}
    return out


def variance_shares(joint: pd.DataFrame, per_group: Dict[str, pd.DataFrame], cols=("log10_t_r", "T_wo_max_K")) -> Dict[str, object]:
    out = {}
    for c in cols:
        vj = float(joint.loc[joint.converged == True, c].var())  # noqa: E712
        shares = {g: float(d.loc[d.converged == True, c].var()) / vj for g, d in per_group.items()}  # noqa: E712
        out[c] = {"joint_variance": vj, "shares": shares, "sum_of_shares": float(sum(shares.values())), "interaction_residual": float(1.0 - sum(shares.values()))}
    return out


def sobol_on_surrogate(joint: pd.DataFrame, cols=("log10_t_r", "T_wo_max_K"), n_fit: int = 1000, N: int = 512, seed: int = 0) -> Dict[str, object]:
    """Group total/first-order Sobol indices from a GP surrogate of the joint MC sample (SALib Saltelli)."""
    from SALib.analyze import sobol
    from SALib.sample import saltelli
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    d = joint[joint.converged == True]  # noqa: E712
    X = d[PARAMS].to_numpy(float); idx = np.random.default_rng(seed).choice(len(X), min(n_fit, len(X)), replace=False)
    problem = {"num_vars": len(PARAMS), "names": PARAMS, "bounds": [[float(d[p].min()), float(d[p].max())] for p in PARAMS]}
    Xs = saltelli.sample(problem, N, calc_second_order=False)
    out = {}
    for c in cols:
        y = d[c].to_numpy(float)
        gp = Pipeline([("s", StandardScaler()), ("gp", GaussianProcessRegressor(ConstantKernel() * Matern(np.ones(len(PARAMS)), (1e-2, 1e3), nu=2.5) + WhiteKernel(1e-4, (1e-8, 1e-1)), normalize_y=True, random_state=seed))])
        gp.fit(X[idx], y[idx])
        r2 = float(gp.score(X, y))
        Si = sobol.analyze(problem, gp.predict(Xs), calc_second_order=False, print_to_console=False)
        S1 = dict(zip(PARAMS, map(float, Si["S1"]))); ST = dict(zip(PARAMS, map(float, Si["ST"])))
        out[c] = {"surrogate_R2_on_joint": r2, "S1": S1, "ST": ST,
                  "group_S1": {g: float(sum(S1[p] for p in ps)) for g, ps in GROUPS.items()}, "group_ST": {g: float(sum(ST[p] for p in ps)) for g, ps in GROUPS.items()}}
    return out


def robustness_fractions(mc: Dict[str, pd.DataFrame], mc_var: Dict[str, pd.DataFrame], weights_s3=(0.25, 0.5, 0.25)) -> Dict[str, object]:
    """Paired fractions of samples supporting the RQ2/RQ3 conclusions, with nominal and load-dependent L_q."""
    def paired(a, b, col="life_rate_per_kmol_H2"):
        ok = (a.converged == True) & (b.converged == True)  # noqa: E712
        return a.loc[ok, col].to_numpy(), b.loc[ok, col].to_numpy()
    out = {}
    for tag, M in (("nominal_L_q", mc), ("load_dependent_L_q", mc_var)):
        res = {}
        b = M["S1_base"]
        if "S5_y4_holdslip_end" in M:
            s5, bb = paired(M["S5_y4_holdslip_end"], b); res["S5y4_life_gt_2x_S1"] = float(np.mean(s5 > 2.0 * bb))
        if "RQ3_knee" in M:
            k, bb = paired(M["RQ3_knee"], b); res["knee_life_lt_0.5_base"] = float(np.mean(k < 0.5 * bb))
        if "S6_SC3.5" in M and "S6_SC2.5" in M:
            ok = (M["S6_SC3.5"].converged == True) & (M["S6_SC2.5"].converged == True)  # noqa: E712
            res["SC3.5_less_life_per_kmol_than_SC2.5"] = float(np.mean(M["S6_SC3.5"].loc[ok, "life_rate_per_kmol_H2"].to_numpy() < M["S6_SC2.5"].loc[ok, "life_rate_per_kmol_H2"].to_numpy()))
        if all(k in M for k in ("load_0.70", "load_1.00")):
            ok = np.ones(len(b), bool)
            for k in ("load_0.70", "load_0.85", "load_1.00"):
                if k in M: ok &= (M[k].converged == True).to_numpy()  # noqa: E712
            ok &= (b.converged == True).to_numpy()  # noqa: E712
            # S2: 16 h at 1.0 + 8 h at 0.70 per day (ramps ignored); per-kmol = sum(rate)/sum(H2)
            d2 = 16 * M["load_1.00"].loc[ok, "life_rate_per_h"].to_numpy() + 8 * M["load_0.70"].loc[ok, "life_rate_per_h"].to_numpy()
            h2 = 16 * M["load_1.00"].loc[ok, "H2_net_kmol_h"].to_numpy() + 8 * M["load_0.70"].loc[ok, "H2_net_kmol_h"].to_numpy()
            res["S2_life_per_kmol_le_S1"] = float(np.mean(d2 / h2 <= b.loc[ok, "life_rate_per_kmol_H2"].to_numpy()))
            if "load_0.85" in M:
                w = weights_s3; keys = ("load_0.70", "load_0.85", "load_1.00")
                d3 = sum(wi * M[k].loc[ok, "life_rate_per_h"].to_numpy() for wi, k in zip(w, keys)); h3 = sum(wi * M[k].loc[ok, "H2_net_kmol_h"].to_numpy() for wi, k in zip(w, keys))
                res["S3_life_per_kmol_le_S1"] = float(np.mean(d3 / h3 <= b.loc[ok, "life_rate_per_kmol_H2"].to_numpy()))
        out[tag] = res
    return out
