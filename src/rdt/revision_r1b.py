"""Revision R1 of HE-D-26-14265, final checks (CLAUDE_CODE_TASK.md).

``python -m rdt.revision_r1b [--calc chain curves checks]``

Reads the published v3 results and the ``data/revision_r1`` outputs; writes only to
``data/revision_r1b/``. No model code is changed: the bed-side chain of ``rdt.reactor1d`` was
compared with p. 100 of Xu & Froment (1989, Part II) and is identical to the print.

chain   chain_plantA_profile.csv   alpha_i(chain) / alpha_i(calibrated) along the Plant A tube
        chain_terms_midheight.csv  Table S17 recomputed (same function as revision_r1.htc_chain_terms)
        xf_chain_run.json          stand-alone Xu-Froment run with the chain; Leva-Grummer residual ranges
curves  four_curves.csv            the four rupture curves, one row per curve (Table S22)
checks  checks.json                timing, adopted row of Tables S18/S19, physics Table 2, over-firing, S20
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd

from rdt import config as C
from rdt import revision_r1 as R1

OUT = C.ROOT / "data" / "revision_r1b"
CFG = C.V3
Z_LO, Z_HI = 0.3, 11.12      # heated length of the Xu-Froment tube compared with the digitised Fig. 3


def _dump(name: str, obj) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(obj, indent=2, default=lambda v: v.item() if hasattr(v, "item") else str(v)))


# ---------------------------------------------------------------------------
# Task 1: the Xu-Froment chain (identical to the print; 1.4 b, c, d)
# ---------------------------------------------------------------------------
def quoted_alternative_form(lam_g, Re, Pr, d_p, d_ti, lam_er0) -> float:
    """The form quoted in CLAUDE_CODE_TASK.md from papers citing De Wasch & Froment (1972).

    NOT what Xu & Froment print (p. 100 of Part II) and not used by the package; evaluated here only
    to reproduce the 1640 / 1170 W m-2 K-1 of the task description from the Table S17 inputs.
    """
    lam_er = lam_er0 + 0.111 * lam_g * Re * Pr ** (1 / 3) / (1 + 46 * (d_p / d_ti) ** 2)
    alpha_w = 8.694 * lam_er0 / d_ti ** (4 / 3) + 0.512 * lam_g * d_ti * Re * Pr ** (1 / 3) / d_p
    return 8 * lam_er * alpha_w / (8 * lam_er + alpha_w * d_ti)


def chain_plantA_profile() -> pd.DataFrame:
    """alpha_i of the printed chain against the calibrated alpha_i at every node of the Plant A run."""
    from rdt import latham_cases as lc
    from rdt import reactor1d as r1
    P, _ = R1._adopted()
    row = lc.load_cases().set_index("case").loc["Plant_A"]
    res = lc.run_case(row, P); tr = res.tube if hasattr(res, "tube") else res
    tube, bed, *_ = lc.build_case(row, P)
    _, _, mw = r1._gas(); rows = []
    for k in range(len(tr.z)):
        if not (tr.z[k] <= P.heated_length_m + 1e-9 and tr.alpha_i[k] > 0):
            continue
        F = np.array([tr.F[s][k] for s in r1.SPECIES]); X = {s: F[i] / F.sum() for i, s in enumerate(r1.SPECIES)}
        pr = r1._gas_props(tr.T[k], tr.P[k], X)
        G = sum(F[i] * mw[s] for i, s in enumerate(r1.SPECIES)) / 3600.0 / tube.A_cs
        lam0 = r1.kunii_smith_lambda_er0(pr["lam"], bed.lambda_s, bed.voidage, float(tr.T[k]), bed.emissivity, bed.d_p)
        ch = r1.xu_froment_alpha_i(pr["lam"], pr["mu"], pr["cp_mass"], G, bed.d_p, tube.d_i, lam0)
        rows.append({"z_m": float(tr.z[k]), "z_frac": float(tr.z[k] / P.heated_length_m), "T_gas_K": float(tr.T[k]), "Re_p": ch["Re"], "Pr": ch["Pr"],
                     "alpha_i_chain": ch["alpha_i"], "alpha_i_calibrated": float(tr.alpha_i[k]), "ratio": ch["alpha_i"] / float(tr.alpha_i[k])})
    return pd.DataFrame(rows)


def _residuals(res, cv, curve: str) -> np.ndarray:
    from rdt import validate_xf1989 as xf
    z, v = cv[curve]; m = (z >= Z_LO) & (z <= Z_HI)
    return xf.model_curve(res, curve, z[m]) - v[m]


def xf_runs() -> Dict[str, object]:
    """Stand-alone Xu-Froment case, settings of the Leva-Grummer row of Table S5 (f_htg 1, eta 0.1,
    voidage 0.526, d_p 9.44 mm, Latham inlet rule), once with each correlation."""
    from rdt import validate_xf1989 as xf
    tidy = xf.load_tidy(); cv = xf.curves(tidy); wall = xf.wall_from_digitised(tidy)
    out = {}
    for ht in ("xu_froment", "leva_grummer"):
        res = xf.run_case(ht, "latham", wall)
        m = xf.metrics_zmin(res, cv)
        eT, ex = _residuals(res, cv, "T_gas"), _residuals(res, cv, "x_CH4")
        out[ht] = {"rmse_x_CH4": m["x_CH4"]["rmse"], "rmse_dx_increment": xf.increment_metric(res, cv)["rmse"], "rmse_T_gas_K": m["T_gas"]["rmse"],
                   "rmse_T_wi_K": m["T_wall_inner"]["rmse"], "rmse_p_bar": m["p_t"]["rmse"],
                   "T_gas_model_minus_digitised_K": {"min": float(eT.min()), "max": float(eT.max()), "n": int(len(eT))},
                   "x_CH4_model_minus_digitised": {"min": float(ex.min()), "max": float(ex.max()), "n": int(len(ex))},
                   "z_range_m": [Z_LO, Z_HI],
                   "alpha_i_W_m2K": {"z_0": float(res.alpha_i[0]), "z_5.56": float(np.interp(5.56, res.z, res.alpha_i)), "z_11": float(np.interp(11.0, res.z, res.alpha_i))}}
    return out


def calc_chain(log=print) -> Dict[str, object]:
    R1.setup()
    prof = chain_plantA_profile(); prof.to_csv(OUT / "chain_plantA_profile.csv", index=False)
    terms = R1.htc_chain_terms(); terms.to_csv(OUT / "chain_terms_midheight.csv", index=False)
    old = pd.read_csv(R1.OUT / "htc_chain_terms.csv")
    dev = float(np.nanmax(np.abs(terms.value.to_numpy(float) / old.value.to_numpy(float) - 1)))
    t = terms.set_index(["case", "term"]).value
    alt = {}
    for case in ("Xu-Froment 1989", "Plant A"):
        g = lambda k: float(t[(case, k)])  # noqa: E731
        alt[case] = quoted_alternative_form(g("lambda_g"), g("Re_p"), g("Pr"), g("d_p equivalent"), g("d_t inner"), g("lambda_er0 (static, Kunii-Smith)"))
    i, j = prof.ratio.idxmin(), prof.ratio.idxmax()
    res = {"verdict": "identical", "chain_terms_max_rel_dev_vs_revision_r1": dev,
           "plantA_ratio_chain_over_calibrated": {"min": float(prof.ratio[i]), "z_frac_min": float(prof.z_frac[i]), "max": float(prof.ratio[j]),
                                                  "z_frac_max": float(prof.z_frac[j]), "n_nodes": int(len(prof))},
           "midheight": {c: {"chain": float(t[(c, "alpha_i, Xu-Froment chain")]), "leva_grummer_f1": float(t[(c, "alpha_i, Leva-Grummer (f_htg = 1)")]),
                             "back_calculated": float(t[(c, "alpha_i, back-calculated")]), "calibrated": float(t[(c, "alpha_i, calibrated")])} for c in alt},
           "quoted_alternative_form_not_in_XF_W_m2K": alt, "xf_runs": xf_runs()}
    _dump("xf_chain_run.json", res)
    log(f"chain: Plant A ratio {res['plantA_ratio_chain_over_calibrated']}; XF chain run {res['xf_runs']['xu_froment']}")
    return res


# ---------------------------------------------------------------------------
# Task 2: four rupture curves (re-scoring only)
# ---------------------------------------------------------------------------
CURVE_LABEL = {"centralloy_g_4852": "Centralloy G 4852, lower band (base)", "centralloy_g_4852_micro": "Centralloy G 4852 Micro",
               "centralloy_et_45_micro": "Centralloy ET 45 Micro", "legacy_yeh_manaurite_xm": "Manaurite XM, minimum curve"}
SCEN = {"S2": "S2_daily", "S3": "S3_renewable", "S6_SC2.5": "S6_SC2.5", "S6_SC3.5": "S6_SC3.5", "S5_slip_y4": "S5_hold_CH4_slip_ageing_y4"}
MC_ITEMS = {"knee": "RQ3_knee", "S5_y4_end": "S5_y4_holdslip_end", "SC3.5": "S6_SC3.5"}


def scatter_sigma(key: str) -> float:
    """sigma of the log10 t_r scatter of one alloy, by the rule of ``uq.build_spec``."""
    from rdt import creep, uq
    a = creep.load_alloy(key); s = a.scatter.summary(uq.SCATTER_REF_T_K)["median_decades"]
    return float(s if a.scatter.assumed else s / 1.645)


def rescore_mc(df: pd.DataFrame, key: str, base_key: str) -> pd.DataFrame:
    """One stored Monte Carlo file re-scored on another curve, sample by sample at equal quantiles:
    C = C_alloy + (C_LM - C_base) (both uniform +/- 1), scatter scaled by sigma_alloy / sigma_base."""
    from rdt import creep, uq
    dC = creep.load_alloy(key).C - creep.load_alloy(base_key).C
    k = scatter_sigma(key) / scatter_sigma(base_key)
    lt = np.full(len(df), np.nan)
    for i, (T, s, c) in enumerate(zip(df.T_wo_max_K.to_numpy(float), df.sigma_hot_MPa.to_numpy(float), df.C_LM.to_numpy(float))):
        if np.isfinite(T) and np.isfinite(s):
            lt[i] = np.log10(float(uq.curve_with_C(c + dC, key).time_to_rupture(T, s)))
    df = df.copy()
    df["log10_t_r"] = lt + k * df.log10_tr_scatter.to_numpy(float)
    df["life_rate_per_h"] = 10.0 ** -df.log10_t_r; df["life_rate_per_kmol_H2"] = df.life_rate_per_h / df.H2_net_kmol_h
    return df


def calc_curves(log=print) -> pd.DataFrame:
    from rdt import alloy_comparison as ac
    from rdt import creep, uq
    R1.setup()
    keys = list(ac.ALLOYS); base = CFG.creep_alloy
    cj = json.loads((R1.OUT / "calc_C.json").read_text()); hot = cj["base_reference"]
    reg = pd.read_csv(R1.OUT / "pareto_regimes_across_curves.csv", index_col=0)
    states = ac.load_hourly_states(CFG, ("S1_steady", *SCEN.values()))
    summ = pd.read_csv(CFG.scenarios_summary).set_index("scenario").life_consumption_per_kmol_H2_rel_S1
    d = Path(uq.OUT_DIR)
    mc = {n: pd.read_csv(d / f"{uq.MC_TAG}_{n}.csv.gz") for n in uq.POINT_NAMES}
    grp = {g: pd.read_csv(d / f"{uq.MC_TAG}_S1_base_group_{g}.csv.gz") for g in uq.GROUPS}
    pub = json.loads(Path(uq.SUMMARY_JSON).read_text())
    rows = []
    for k in keys:
        cur = creep.load_alloy(k).minimum
        T, s = hot["T_wo_max_K"], hot["sigma_hot_MPa"]
        lt = lambda x: float(np.log10(cur.time_to_rupture(x, s)))  # noqa: E731
        slope = (lt(T + 0.5) - lt(T - 0.5)) / 1.0
        r = {"curve": CURVE_LABEL[k], "key": k, "C": creep.load_alloy(k).C, "scatter_sigma_decades": scatter_sigma(k), "T_hot_K": T, "sigma_hot_MPa": s,
             "log10_tr_hot": lt(T), "slope_decades_per_K": slope, "halving_dT_K": float(np.log10(2.0) / abs(slope))}
        # nominal ratios: regimes (base = calibrated point), scenarios (base = S1)
        for name in R1.TABLE3_REGIMES:
            r[f"nom_{name}"] = float(reg.loc[name, f"fixed_{k}"])
            r[f"same_member_{name}"] = bool(int(reg.loc[name, f"reselected_row_{k}"]) == int(reg.loc[name, "pareto_row"]))
        sl = ac.scenario_life_per_kmol(states, cur)
        for lab, sc_ in SCEN.items():
            r[f"nom_{lab}"] = sl[sc_]
            if k == base:
                assert abs(sl[sc_] / summ[sc_] - 1) < 1e-9, sc_
        # Monte Carlo, re-scored
        M = {n: rescore_mc(v, k, base) for n, v in mc.items()}; M["load_1.00"] = M["S1_base"]
        if k == base:
            ok = R1._ok(mc["S1_base"])
            dev = float(np.nanmax(np.abs(M["S1_base"].log10_t_r.to_numpy()[ok] - mc["S1_base"].log10_t_r.to_numpy()[ok])))
            assert dev < 1e-3, dev
            r["base_reproduction_max_abs_dev_log10tr"] = dev
        for lab, n in MC_ITEMS.items():
            q = uq.intervals(M[n], M["S1_base"])["life_rate_per_kmol_H2_rel_base"]
            r.update({f"mc_{lab}_median": q["median"], f"mc_{lab}_p05": q["p05"], f"mc_{lab}_p95": q["p95"]})
        q = uq.intervals(M["S1_base"], M["S1_base"])["log10_t_r"]
        r.update({"mc_S1_log10tr_median": q["median"], "mc_S1_log10tr_p05": q["p05"], "mc_S1_log10tr_p95": q["p95"]})
        fr = uq.robustness_fractions(M, M)["nominal_L_q"]
        r.update({f"frac_{a}": b for a, b in fr.items()})
        if k == base:
            for a, b in pub["robustness_fractions"]["nominal_L_q"].items():
                assert abs(fr[a] - b) < 1e-12, (a, fr[a], b)
        # variance shares: T_wo is not re-scored (the curve does not enter the thermal problem), log10 t_r is
        G = {g: rescore_mc(v, k, base) for g, v in grp.items()}
        vs = uq.variance_shares(M["S1_base"], G)
        for g in uq.GROUPS:
            r[f"share_Two_{g}"] = vs["T_wo_max_K"]["shares"][g]; r[f"share_log10tr_{g}"] = vs["log10_t_r"]["shares"][g]
        c = cj["curves"][k]
        r.update({"spearman_vs_base": c["spearman_vs_base"], "n_nondominated": c["n_nondominated"], "left_nondominated_set": c["lost_nondominance"],
                  "regimes_same_members": all(r[f"same_member_{n}"] for n in R1.TABLE3_REGIMES)})
        rows.append(r)
        log(f"curves {CURVE_LABEL[k][:28]:28s} C {r['C']:.2f} knee {r['nom_knee']:.3f} mc {r['mc_knee_median']:.3f}  frac "
            + " ".join(f"{fr[a]:.3f}" for a in fr))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "four_curves.csv", index=False)
    return df


def latex_table(df: pd.DataFrame) -> str:
    """Booktabs Table S22 in two panels (nominal ratios; Monte Carlo and front), one row per curve."""
    f3 = lambda r, a: f"{r[f'mc_{a}_median']:.2f} [{r[f'mc_{a}_p05']:.2f}, {r[f'mc_{a}_p95']:.2f}]"  # noqa: E731
    short = {"centralloy_g_4852": "G 4852 (base)", "centralloy_g_4852_micro": "G 4852 Micro", "centralloy_et_45_micro": "ET 45 Micro",
             "legacy_yeh_manaurite_xm": "Manaurite XM"}
    A = ["\\begin{tabular}{lcccc ccccc ccccc}", "\\toprule",
         " & & & & & \\multicolumn{5}{c}{Regimes (base = calibrated point)} & \\multicolumn{5}{c}{Scenarios (base = S1)} \\\\",
         "\\cmidrule(lr){6-10}\\cmidrule(lr){11-15}",
         "Curve & $C$ & $\\log_{10} t_r$ (h) & $\\mathrm{d}\\log_{10} t_r/\\mathrm{d}T$ (K$^{-1}$) & $\\Delta T_{1/2}$ (K) & min life, iso-H$_2$ & knee & max H$_2$ & "
         "min fuel & min life & S2 & S3 & S/C 2.5 & S/C 3.5 & aging y4 \\\\", "\\midrule"]
    B = ["\\begin{tabular}{l ccc ccccc ccc}", "\\toprule",
         " & \\multicolumn{3}{c}{Paired Monte Carlo, median [5--95\\,\\%]} & \\multicolumn{5}{c}{Robustness fractions, $p = 0$} & \\multicolumn{3}{c}{Pareto front (198 members)} \\\\",
         "\\cmidrule(lr){2-4}\\cmidrule(lr){5-9}\\cmidrule(lr){10-12}",
         "Curve & knee & end of aging year 4 & S/C 3.5 & y4 $> 2\\times$ S1 & knee $< 0.5\\times$ & S/C 3.5 $<$ 2.5 & S2 $\\le$ S1 & S3 $\\le$ S1 & "
         "Spearman & left front & same regimes \\\\", "\\midrule"]
    for _, r in df.iterrows():
        n = short[r["key"]]
        A.append(" & ".join([n, f"{r['C']:.2f}", f"{r['log10_tr_hot']:.2f}", f"${r['slope_decades_per_K']:.4f}$", f"{r['halving_dT_K']:.1f}",
                             *(f"{r[f'nom_{k}']:.3f}" for k in ("min_life_iso_H2", "knee", "max_H2", "min_fuel_overall", "min_life_overall")),
                             *(f"{r[f'nom_{k}']:.3f}" for k in ("S2", "S3", "S6_SC2.5", "S6_SC3.5", "S5_slip_y4"))]) + " \\\\")
        B.append(" & ".join([n, f3(r, "knee"), f3(r, "S5_y4_end"), f3(r, "SC3.5"),
                             *(f"{r[f'frac_{k}']:.3f}" for k in ("S5y4_life_gt_2x_S1", "knee_life_lt_0.5_base", "SC3.5_less_life_per_kmol_than_SC2.5",
                                                                 "S2_life_per_kmol_le_S1", "S3_life_per_kmol_le_S1")),
                             f"{r['spearman_vs_base']:.4f}", f"{int(r['left_nondominated_set'])}", "yes" if r["regimes_same_members"] else "no"]) + " \\\\")
    end = ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(A + end + ["", "\\medskip", ""] + B + end)


# ---------------------------------------------------------------------------
# Task 3: small consistency checks
# ---------------------------------------------------------------------------
def check_timing(n_phys: int = 20) -> Dict[str, object]:
    import joblib
    from rdt import operate as op
    from rdt import surrogate as sg
    M = json.loads(Path(CFG.surrogate_metrics).read_text())
    meta = json.loads(Path(CFG.pareto_meta).read_text()); x = meta["base_reference"]["inputs"]; base = op.base_case()
    op.run_case(x, base)                                   # warm-up (Cantera objects, caches)
    t = []
    for _ in range(n_phys):
        t0 = time.perf_counter(); op.run_case(x, base); t.append(time.perf_counter() - t0)
    best = M["targets"]["T_wo_max_K"]["best"]
    b = joblib.load(Path(sg.MODELS_DIR) / f"T_wo_max_K__{best}.joblib")
    tr, _ = sg.load_data(); tr = sg.add_physics_features(tr); X = tr[b["features"]].to_numpy(float)
    Xb = X[np.random.default_rng(0).integers(0, len(X), 100_000)]
    t0 = time.perf_counter(); b["pipeline"].predict(Xb); batch = (time.perf_counter() - t0) / 1e5
    phys = float(np.median(t))
    return {"stored": {**M["speed"], "file": str(Path(CFG.surrogate_metrics).relative_to(C.ROOT)), "key": "speed (written by pipeline.py:72)", "surrogate": best,
                       "all_predict_us": {k: {m: v["predict_time_per_1e5_s"] * 10 for m, v in d["models"].items() if m == d["best"]} for k, d in M["targets"].items()}},
            "measured_now": {"physics_s_per_run_median": phys, "physics_s_min": float(min(t)), "physics_s_max": float(max(t)), "n": n_phys,
                             "surrogate_s_per_prediction_batch_1e5": batch, "speedup": phys / batch}}


def check_adopted_row() -> Dict[str, object]:
    """Tables S18/S19 at the exact adopted fit: f_htg fixed at its unrounded value, F_gt and L_q refitted."""
    P, _ = R1._adopted()
    o = R1._sweep_task("f_htg", float(P.f_htg), ("F_gt", "L_q"))
    row = R1._sweep_table([o]).iloc[0].to_dict()
    old = {"S18 (f_htg = 1.83)": pd.read_csv(R1.OUT / "htc_sweep.csv").query("value == 1.83").iloc[0].to_dict(),
           "S19 (L_q = 0.456)": pd.read_csv(R1.OUT / "flame_length_sweep.csv").query("value == 0.456").iloc[0].to_dict()}
    return {"adopted": {"f_htg": float(P.f_htg), "L_q": float(P.L_q), "F_gt": float(P.F_gt)}, "exact_row": row, "printed_rows": old}


def check_table2_physics() -> Dict[str, object]:
    """S1 and S6 with the physics model: firing solved to hold T_out, as revision_r1._sweep_task does."""
    from rdt import scenarios as sc
    P, _ = R1._adopted(); b = sc.base_inputs(1)
    x0 = {"steam_to_carbon": float(b["steam_to_carbon"][0]), "feed_per_tube_fraction": 1.0, "inlet_T": float(b["inlet_T"][0]),
          "inlet_P": float(b["inlet_P"][0]), "catalyst_activity": R1.A_CAL, "excess_air": float(b["excess_air"][0])}
    summ = pd.read_csv(CFG.scenarios_summary).set_index("scenario")
    out = {}
    for name, scv in (("S1_steady", x0["steam_to_carbon"]), ("S6_SC2.5", 2.5), ("S6_SC3.5", 3.5)):
        r = R1._solve_P({**x0, "steam_to_carbon": scv}, P, "T_out_K", sc.T_OUT_TARGET)
        out[name] = {"steam_to_carbon": scv, "firing": r["firing"], "H2_net_kmol_h": r["H2_net_kmol_h"], "annual_H2_kmol_per_tube": 8760 * r["H2_net_kmol_h"],
                     "annual_damage_D": 8760 * r["life_consumption_rate_per_h"], "T_wo_max_K": r["T_wo_max_K"],
                     "rate_per_kmol": r["life_consumption_rate_per_h"] / r["H2_net_kmol_h"],
                     "surrogate": {"annual_H2_kmol_per_tube": float(summ.loc[name, "annual_H2_kmol_per_tube"]), "annual_damage_D": float(summ.loc[name, "annual_damage_D"]),
                                   "rel_life_cost": float(summ.loc[name, "life_consumption_per_kmol_H2_rel_S1"]), "firing": float(summ.loc[name, "firing_mean"])}}
    for v in out.values():
        v["rel_life_cost"] = v["rate_per_kmol"] / out["S1_steady"]["rate_per_kmol"]
    return out


def check_overfiring() -> Dict[str, object]:
    from rdt import creep
    cur = creep.load_alloy(CFG.creep_alloy).minimum
    summ = pd.read_csv(CFG.scenarios_summary).set_index("scenario")
    rate = lambda h: 1.0 / np.asarray(cur.time_to_rupture(h.T_wo_max_K.to_numpy(float), h.sigma_hot_MPa.to_numpy(float)), float)  # noqa: E731
    s1 = rate(pd.read_csv(CFG.scenario_hourly("S1_steady")))
    out = {}
    for s in ("S4a_mild_overfire", "S4b_severe_overfire", "S4c_hot_band"):
        d = rate(pd.read_csv(CFG.scenario_hourly(s))); ev = d > 1.0001 * s1
        out[s] = {"event_hours_summary": float(summ.loc[s, "event_hours"]), "hours_above_S1_rate": int(ev.sum()),
                  "damage_rate_factor_summary": float(summ.loc[s, "dD_event_to_base_ratio"]),
                  "damage_rate_factor_mean_over_elevated_hours": float(d[ev].mean() / s1.mean()), "damage_rate_factor_peak": float(d.max() / s1.mean()),
                  "share_of_annual_damage_in_elevated_hours_pct": float(100 * d[ev].sum() / d.sum()),
                  "excess_over_S1_share_of_annual_damage_pct": float(100 * (d.sum() - s1.sum()) / d.sum()),
                  "rel_life_cost": float(summ.loc[s, "life_consumption_per_kmol_H2_rel_S1"])}
    return out


def check_s20() -> Dict[str, object]:
    ci = json.loads((R1.OUT / "calc_I.json").read_text())
    df = pd.read_csv(R1.OUT / "walltemp_error_propagation.csv")
    (x1, y1), (x2, y2) = R1.LOCO_U
    return {"line_points": R1.LOCO_U, "slope_K_per_unit_load": (y2 - y1) / (x2 - x1), "w": R1.W_HOT, "w_formula": "(0.682 - 0.43) / (0.682 - 0.293)",
            "e_u_at_load": {f"{v:.2f}": float(R1.e_upper(v)) for v in (0.60, 0.62, 0.70, 0.71, 0.85, 0.97, 0.99, 1.00, 1.10)},
            "stored_calc_I": {"e_u_line": ci["e_u_line"], "w": ci["w"]},
            "table3_shift_range_K": {v: [float(g.dT_item_K.min()), float(g.dT_item_K.max())] for v, g in df[(df.table == "Table 3") & df.variant.str.startswith("b")].groupby("variant")}}


def calc_checks(log=print) -> Dict[str, object]:
    R1.setup()
    out = {}
    for name, f in (("timing", check_timing), ("overfiring", check_overfiring), ("table_S20", check_s20), ("table2_physics", check_table2_physics),
                    ("adopted_row_S18_S19", check_adopted_row)):
        t0 = time.perf_counter(); out[name] = f(); log(f"checks {name} ({time.perf_counter() - t0:.0f}s)")
        _dump("checks.json", out)
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--calc", nargs="+", default=["chain", "curves", "checks"], choices=["chain", "curves", "checks"])
    a = p.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    log = lambda s: print(time.strftime("[%H:%M:%S] ") + s, flush=True)  # noqa: E731
    if "chain" in a.calc:
        calc_chain(log)
    if "curves" in a.calc:
        df = calc_curves(log); (OUT / "four_curves_table.tex").write_text(latex_table(df) + "\n")
    if "checks" in a.calc:
        calc_checks(log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
