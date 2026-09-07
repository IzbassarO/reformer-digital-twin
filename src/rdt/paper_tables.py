"""Generate the booktabs tables and number macros of the *archived* sectioned draft from the v2 results.

The manuscript under ``paper/`` is now the consolidated Overleaf pair, which carries its numbers
inline rather than through ``\\input{numbers}``. This generator therefore writes into
``paper/archive/draft_step21/`` -- it still reproduces the draft it was written for, and it can no
longer drop stale ``tables/`` and ``numbers.tex`` files next to the current source of truth. Nothing
under ``data/`` is modified.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import yaml

from rdt import config as C

DRAFT = C.ROOT / "paper" / "archive" / "draft_step21"
TAB = DRAFT / "tables"
NUMBERS = DRAFT / "numbers.tex"
cfg = C.V2
LV = C.DATA / "literature_validation"


def _tex(s: str) -> str:
    return str(s).replace("_", "\\_").replace("%", "\\%").replace("&", "\\&")


def _table(body_rows: List[List[str]], header: List[str], caption: str, label: str, colspec: str = None, note: str = None, small=True) -> str:
    colspec = colspec or "l" + "r" * (len(header) - 1)
    wide = len(header) > 5 or label in ("tab:uq", "tab:var")  # shrink wide tables to the text width
    lines = ["\\begin{table}[htbp]", "\\centering", f"\\caption{{{caption}}}", f"\\label{{{label}}}"]
    if small: lines.append("\\small")
    if wide: lines.append("\\resizebox{\\textwidth}{!}{%")
    lines += [f"\\begin{{tabular}}{{{colspec}}}", "\\toprule", " & ".join(header) + " \\\\", "\\midrule"]
    lines += [" & ".join(r) + " \\\\" for r in body_rows]
    lines += ["\\bottomrule", "\\end{tabular}" + ("}" if wide else "")]
    if note: lines.append(f"\\par\\smallskip\\footnotesize {note}")
    lines.append("\\end{table}")
    return "\n".join(lines) + "\n"


def f(x, nd=2):
    try:
        return f"{float(x):.{nd}f}"
    except Exception:  # noqa: BLE001
        return str(x)


def g(x, nd=3):
    try:
        return f"{float(x):.{nd}g}"
    except Exception:  # noqa: BLE001
        return str(x)


def table_xf_verification() -> str:
    d = yaml.safe_load((LV / "xf1989_fit.yaml").read_text())
    v = d["verification_matched_alpha"]; runs = v["runs"]; lg = d["xu_froment_correlation"]["combinations"]["leva_grummer__inlet_latham"]
    rows = []
    def row(lab, m, inc):
        return [lab, f(m["x_CH4"]["rmse"], 3), f(inc["rmse"], 3), f(m["T_gas"]["rmse"], 1), f(m["T_wall_inner"]["rmse"], 1), f(m["p_t"]["rmse"], 2)]
    rows.append(row("Leva--Grummer, $f_{htg}=1$, $\\eta=0.1$", lg["metrics"], {"rmse": float("nan")}))
    for k, lab in (("alpha_const_median__eta_0.1", f"matched $\\alpha_i$ = {v['alpha_i_back_calculated']['median_W_m2K']:.0f} W m$^{{-2}}$ K$^{{-1}}$, $\\eta=0.1$"),
                   ("alpha_const_median__eta_latham_inlet", "matched $\\alpha_i$, $\\eta$ = 0.05 / 0.1 (adopted)"), ("alpha_profile_pchip__eta_latham_inlet", "matched $\\alpha_i(z)$ profile, $\\eta$ = 0.05 / 0.1")):
        r = runs[k]; rows.append(row(lab, r["metrics"], r["x_CH4_increment"]))
    cap = ("Verification of the tube model against the digitised Fig.~3 of Xu and Froment \\cite{xu1989b}: RMSE of methane conversion $x_{\\mathrm{CH_4}}$, of the conversion increment "
           "$x(z)-x(0.5\\,\\mathrm{m})$, of the process-gas and inner-wall temperatures ($z\\ge0.3$~m) and of pressure. Source: \\texttt{data/literature\\_validation/xf1989\\_fit.yaml} (verification\\_matched\\_alpha).")
    return _table(rows, ["Configuration", "RMSE $x_{\\mathrm{CH_4}}$", "RMSE $\\Delta x$", "RMSE $T_g$ (K)", "RMSE $T_{wi}$ (K)", "RMSE $p$ (bar)"], cap, "tab:xf", note="Back-calculated $\\alpha_i$: median 348, interquartile range 313--363 W m$^{-2}$ K$^{-1}$ over 0.5--11 m.")


def table_latham_cases() -> str:
    df = pd.read_csv(LV / "latham2008_plant_cases.csv"); L = yaml.safe_load((LV / "latham_fit.yaml").read_text())
    err = L["calibration_alpha_top_fixed"]["best_fit"]["errors_pred_minus_meas"]; cv = L["calibration_alpha_top_fixed"]["leave_one_out_cv"]
    rows = []
    for _, r in df.iterrows():
        c = r.case; e = err.get(c); h = cv.get(c, {}).get("held_out_errors")
        rows.append([c.replace("_", " "), f"{int(r.plant_rate_pct)}", f(r.T_in_K, 1), f(r.P_in_Pa / 1e5, 2), f(r.feed_per_tube_kmol_h, 1), f(r.steam_to_carbon_derived, 2), f(r.out_T_K, 1), f(r.TWT_upper_K, 1), f(r.TWT_lower_K, 1),
                     (f"{e['out_T_K']:+.1f}" if e else "--"), (f"{e['TWT_upper_K']:+.1f}" if e else "--"), (f"{e['TWT_lower_K']:+.1f}" if e else "--"),
                     (f"{h['TWT_upper_K']:+.1f}" if h else "--"), (f"{h['TWT_lower_K']:+.1f}" if h else "--")])
    cap = ("Latham plant cases (thesis Appendix H, absolute temperatures) and calibration results: measured inlet and outlet conditions, model errors (predicted minus measured) of the "
           "all-cases fit with $\\alpha_{top}$ fixed, and leave-one-case-out held-out errors. Plant C1 is excluded from the fit (inconsistent feed flow). "
           "Source: \\texttt{latham2008\\_plant\\_cases.csv}, \\texttt{latham\\_fit.yaml} (calibration\\_alpha\\_top\\_fixed).")
    hdr = ["Case", "Rate (\\%)", "$T_{in}$ (K)", "$p_{in}$ (bar)", "Feed (kmol h$^{-1}$)", "S/C", "$T_{out}$ (K)", "$T_{w,u}$ (K)", "$T_{w,l}$ (K)", "$\\Delta T_{out}$", "$\\Delta T_{w,u}$", "$\\Delta T_{w,l}$", "CV $\\Delta T_{w,u}$", "CV $\\Delta T_{w,l}$"]
    return _table(rows, hdr, cap, "tab:latham", colspec="lrrrrrrrrrrrrr", note="Per-tube feed; S/C on total hydrocarbon carbon; upper/lower wall temperatures at 3.66 and 8.53 m of 12.5 m; thesis uncertainties: 2~K outlet, 3~K wall, 8~K flue gas.")


def table_calibration() -> str:
    L = yaml.safe_load((LV / "latham_fit.yaml").read_text())
    rows = []
    for key, lab in (("calibration", "four parameters"), ("calibration_alpha_top_fixed", "$\\alpha_{top}$ fixed (adopted)")):
        b = L[key]["best_fit"]; cov = b["covariance"]
        for k in b["params"]:
            rows.append([lab, f"${{{k}}}$".replace("F_gt", "F_{gt}").replace("L_q", "L_q").replace("alpha_top", "\\alpha_{top}").replace("f_htg", "f_{htg}"), g(b["params"][k], 3), f"[{cov['ci95'][k][0]:.3f}, {cov['ci95'][k][1]:.3f}]", g(cov["se"][k], 2)])
        rows.append([lab, "$\\chi^2$ / dof", g(b["chi2"], 3), f"{cov['dof']}", f"reduced {cov['reduced_chi2']:.2f}"])
    cap = ("Calibrated furnace and heat-transfer parameters on Plants A, B and C2 with 95\\,\\% $t$-intervals from the Jacobian covariance. The four-parameter fit has $\\alpha_{top}$ at its bound with "
           "$r(L_q,\\alpha_{top})=0.985$; the adopted fit fixes $\\alpha_{top}=0.182$. Source: \\texttt{latham\\_fit.yaml}.")
    return _table(rows, ["Fit", "Parameter", "Value", "95\\,\\% interval", "s.e."], cap, "tab:cal", colspec="llrlr")


def table_sobol() -> str:
    s = json.loads(cfg.sobol_json.read_text()); names = list(s["problem"]["names"])
    short = {"steam_to_carbon": "steam-to-carbon", "feed_per_tube_fraction": "load", "specific_firing_factor": "firing", "inlet_T": "inlet $T$", "inlet_P": "inlet $p$", "catalyst_activity": "activity", "excess_air": "excess air"}
    rows = []
    for k in names:
        r = [short[k]]
        for y in ("T_wo_max_K", "log10_t_r_hot", "CH4_slip_dry_pct"):
            o = s["outputs"][y]; r += [f"{o['S1'][k]:.2f} $\\pm$ {o['S1_conf'][k]:.2f}", f"{o['ST'][k]:.2f} $\\pm$ {o['ST_conf'][k]:.2f}"]
        rows.append(r)
    cap = (f"Sobol first-order ($S_1$) and total-effect ($S_T$) indices with 95\\,\\% confidence half-widths from a Saltelli design ($N={s['N']}$, {s['n_runs']} physics runs) over the v2 operating space. "
           f"Source: \\texttt{{{_tex(str(cfg.sobol_json.relative_to(C.ROOT)))}}} ({cfg.tag}).")
    hdr = ["Input", "$S_1$ ($T_{wo,max}$)", "$S_T$ ($T_{wo,max}$)", "$S_1$ ($\\log_{10}t_r$)", "$S_T$ ($\\log_{10}t_r$)", "$S_1$ (CH$_4$ slip)", "$S_T$ (CH$_4$ slip)"]
    return _table(rows, hdr, cap, "tab:sobol")


def table_scenarios() -> str:
    d = pd.read_csv(cfg.scenarios_summary).set_index("scenario"); meta = json.loads(cfg.scenarios_meta.read_text())
    order = ["S1_steady", "S2_daily", "S3_renewable", "S4a_mild_overfire", "S4b_severe_overfire", "S4c_hot_band", "S6_SC2.5", "S6_SC3.5", "S5_hold_CH4_slip_ageing_y1", "S5_hold_CH4_slip_ageing_y4", "S5_hold_T_out_ageing_y1", "S5_hold_T_out_ageing_y4"]
    labels = {"S1_steady": "S1 steady base", "S2_daily": "S2 daily load-following", "S3_renewable": "S3 renewable-following", "S4a_mild_overfire": "S4a mild over-firing (48 h)", "S4b_severe_overfire": "S4b severe over-firing (6 h)", "S4c_hot_band": "S4c hot band +40 K (24 h)",
              "S6_SC2.5": "S6 steam-to-carbon 2.5", "S6_SC3.5": "S6 steam-to-carbon 3.5", "S5_hold_CH4_slip_ageing_y1": "S5 ageing, hold slip, year 1", "S5_hold_CH4_slip_ageing_y4": "S5 ageing, hold slip, year 4", "S5_hold_T_out_ageing_y1": "S5 ageing, hold $T_{out}$, year 1", "S5_hold_T_out_ageing_y4": "S5 ageing, hold $T_{out}$, year 4"}
    rows = []
    for k in order:
        r = d.loc[k]; rows.append([labels[k], f(r.annual_H2_kmol_per_tube / 1e3, 1), g(r.annual_damage_D, 3), f(r.life_consumption_per_kmol_H2_rel_S1, 3), f(r.avg_condition_error, 3), f(r.T_wo_max_mean_K, 0), f(r.T_wo_max_max_K, 0), f(r.firing_mean, 3)])
    cap = (f"RQ2 scenario summary for one year (8760 h) of the Plant A twin at {meta['excess_air_pct']:.1f}\\,\\% excess air: hydrogen per tube, Robinson damage $D$, life consumption per kmol H$_2$ relative to S1, "
           "ratio of the annual-mean-condition damage estimate to the hourly-accumulated damage, mean and maximum hot-spot temperature and mean firing factor. Life uses the Yeh Manaurite XM placeholder curve. "
           f"Source: \\texttt{{{_tex(str(cfg.scenarios_summary.relative_to(C.ROOT)))}}}.")
    hdr = ["Scenario", "H$_2$ (10$^3$ kmol)", "$D$ (--)", "Rel.\\ life cost", "Avg.-cond.\\ ratio", "$\\bar T_{wo,max}$ (K)", "max $T_{wo,max}$ (K)", "Firing"]
    return _table(rows, hdr, cap, "tab:scen")


def table_regimes() -> str:
    reg = pd.read_csv(cfg.regimes_csv, index_col=0); m = json.loads(cfg.pareto_meta.read_text())
    rows = []
    for k in ("base", "min_life_iso_H2", "knee", "max_H2", "min_fuel_overall", "min_life_overall"):
        if k not in reg.index: continue
        r = reg.loc[k]
        rows.append([k.replace("_", " "), f(r.feed_per_tube_fraction, 2), f(r.steam_to_carbon, 2), f(r.inlet_T, 0), f(r.excess_air, 1), f(r.specific_firing_factor, 3), f(r.H2_net_kmol_h_phys, 2), f(r.fuel_MJ_per_kmol_H2_phys, 1), f(r.life_rate_per_kmol_H2_rel_base_phys, 3), f(r.T_wo_max_K_phys, 0), f(r.CH4_slip_dry_pct_phys, 2)])
    cap = (f"RQ3 representative regimes (physics-verified) from the NSGA-II Pareto set ({m['settings']['pop']} $\\times$ {m['settings']['gens']} generations; {m['n_feasible_phys']} of {m['n_pareto']} members feasible under the physics model). "
           "Fuel = firebox fuel per kmol H$_2$; life = consumption rate per kmol H$_2$ relative to the base (Yeh placeholder curve). "
           f"Source: \\texttt{{{_tex(str(cfg.regimes_csv.relative_to(C.ROOT)))}}}.")
    hdr = ["Regime", "Load", "S/C", "$T_{in}$ (K)", "Excess air (\\%)", "Firing", "H$_2$ (kmol h$^{-1}$)", "Fuel (MJ kmol$^{-1}$)", "Rel.\\ life cost", "$T_{wo,max}$ (K)", "CH$_4$ slip (\\%)"]
    return _table(rows, hdr, cap, "tab:regimes", colspec="lrrrrrrrrrr")


def table_steam() -> str:
    s = json.loads(cfg.steam_credit_json.read_text()); rows = []
    for a in s["alphas"]:
        p = s["per_alpha"][str(float(a))]
        for lab in ("min_life", "knee", "min_fuel"):
            if lab in p:
                r = p[lab]; rows.append([f"{a:.2f}", lab.replace("_", " "), f(r["steam_to_carbon"], 2), f(r["load"], 2), f(r["firing"], 3), f(r["excess_air"], 1), f(r["life_rel_base"], 3), f"{r['fuel_total_rel_base_pct']:+.1f}", f(r["T_wo_max_K"], 0)])
    cap = ("Steam-credit sweep of the iso-production (H$_2$ = base $\\pm$1\\,\\%) life--fuel front: $\\alpha=1$ counts firebox fuel only, $\\alpha=0$ charges steam raising and feed preheat in full. "
           f"Physics-verified. Source: \\texttt{{{_tex(str(cfg.steam_credit_json.relative_to(C.ROOT)))}}}.")
    return _table(rows, ["$\\alpha$", "Regime", "S/C", "Load", "Firing", "Excess air (\\%)", "Rel.\\ life cost", "Total fuel vs base (\\%)", "$T_{wo,max}$ (K)"], cap, "tab:steam", colspec="llrrrrrrr")


def table_uq() -> str:
    S = json.loads(cfg.uq_summary.read_text()); rows = []
    labels = {"S1_base": "S1 base", "RQ3_knee": "RQ3 knee", "RQ3_min_life_iso_H2": "RQ3 min-life iso-H$_2$", "S5_y4_holdslip_end": "S5 year 4, hold slip", "S6_SC2.5": "S/C 2.5", "S6_SC3.5": "S/C 3.5", "load_0.70": "load 0.70", "load_0.85": "load 0.85"}
    for n, lab in labels.items():
        i = S["intervals"][n]; t = i["T_wo_max_K"]; l = i["log10_t_r"]; r = i["life_rate_per_kmol_H2_rel_base"]
        rows.append([lab, f"{t['median']:.0f} [{t['p05']:.0f}, {t['p95']:.0f}]", f"{l['median']:.2f} [{l['p05']:.2f}, {l['p95']:.2f}]", f"{r['median']:.3f} [{r['p05']:.3f}, {r['p95']:.3f}]"])
    cap = (f"RQ4 Monte Carlo (scrambled Sobol design, $N={S['N']}$, all uncertainty groups) through the physics model: median and 5--95\\,\\% interval of the hot-spot temperature, of $\\log_{{10}}t_r$ and of the paired "
           f"life consumption per kmol H$_2$ relative to the base. Source: \\texttt{{{_tex(str(cfg.uq_summary.relative_to(C.ROOT)))}}}.")
    t1 = _table(rows, ["Operating point", "$T_{wo,max}$ (K)", "$\\log_{10}t_r$ (h)", "Rel.\\ life cost"], cap, "tab:uq", colspec="llll")
    vs = S["variance_shares_S1_base"]; sob = S.get("sobol_on_surrogate_S1_base", {}); rows2 = []
    for gname in ("kinetics", "heat_transfer", "creep", "measurement"):
        rows2.append([gname.replace("_", " "), f"{vs['log10_t_r']['shares'][gname]:.3f}", f"{sob.get('log10_t_r', {}).get('group_ST', {}).get(gname, float('nan')):.3f}", f"{vs['T_wo_max_K']['shares'][gname]:.3f}", f"{sob.get('T_wo_max_K', {}).get('group_ST', {}).get(gname, float('nan')):.3f}"])
    rows2.append(["sum / interaction residual", f"{vs['log10_t_r']['sum_of_shares']:.3f} / {vs['log10_t_r']['interaction_residual']:+.3f}", "", f"{vs['T_wo_max_K']['sum_of_shares']:.3f} / {vs['T_wo_max_K']['interaction_residual']:+.3f}", ""])
    cap2 = ("Variance decomposition at the S1 base point: one-group-at-a-time share of the joint variance ($N=1024$ per group) and group total-effect Sobol index from a Gaussian-process surrogate of the joint sample. "
            f"Source: \\texttt{{{_tex(str(cfg.uq_summary.relative_to(C.ROOT)))}}}.")
    t2 = _table(rows2, ["Group", "share $\\log_{10}t_r$", "$S_T$ $\\log_{10}t_r$", "share $T_{wo,max}$", "$S_T$ $T_{wo,max}$"], cap2, "tab:var", colspec="lrrrr")
    rf = S["robustness_fractions"]; rows3 = []
    labels3 = {"S5y4_life_gt_2x_S1": "ageing year 4 costs $>2\\times$ S1 life per kmol", "knee_life_lt_0.5_base": "knee regime $<0.5\\times$ base", "SC3.5_less_life_per_kmol_than_SC2.5": "S/C 3.5 cheaper in life than S/C 2.5", "S2_life_per_kmol_le_S1": "S2 daily $\\le$ S1", "S3_life_per_kmol_le_S1": "S3 renewable $\\le$ S1"}
    for k, lab in labels3.items():
        rows3.append([lab, f"{rf['nominal_L_q'][k]:.3f}", f"{rf['load_dependent_L_q'][k]:.3f}"])
    cap3 = ("Robustness of the RQ2/RQ3 conclusions: fraction of paired Monte Carlo samples supporting each statement, with the calibrated flame length $L_q$ and with the structural variant $L_q\\,\\mathrm{load}^{0.5}$. "
            f"Source: \\texttt{{{_tex(str(cfg.uq_summary.relative_to(C.ROOT)))}}}.")
    t3 = _table(rows3, ["Statement", "nominal $L_q$", "$L_q \\propto \\mathrm{load}^{0.5}$"], cap3, "tab:robust", colspec="lrr")
    return t1 + "\n" + t2 + "\n" + t3


def table_surrogate() -> str:
    M = json.loads(cfg.surrogate_metrics.read_text()); rows = []
    units = {"T_wo_max_K": "K", "T_out_K": "K", "P_out_bar": "bar", "CH4_slip_dry_pct": "\\%", "H2_net_kmol_h": "kmol h$^{-1}$", "duty_W": "W", "log10_t_r_hot": "--", "z_frac_T_wo_max": "--"}
    for t, rt in M["targets"].items():
        b = rt["models"][rt["best"]]; c = M["conformal"][t]["levels"]["0.9"]
        rows.append([t.replace("_", "\\_"), units[t], rt["best"].replace("_", "\\_"), g(b["test"]["r2"], 4), g(b["test"]["rmse"], 3), g(b["cv5"]["rmse"], 3), f"{c['picp']:.3f}", g(c["mpiw"], 3)])
    cap = (f"Best surrogate per target ($N_{{train}}$ = {M['n_train']} LHS, $N_{{test}}$ = {M['n_test']} Saltelli points), test $R^2$ and RMSE, 5-fold CV RMSE, and empirical coverage (PICP) and mean width (MPIW) of the nominal 90\\,\\% CV+ conformal interval. "
           f"Source: \\texttt{{{_tex(str(cfg.surrogate_metrics.relative_to(C.ROOT)))}}}.")
    return _table(rows, ["Target", "Unit", "Model", "$R^2$", "RMSE", "CV RMSE", "PICP$_{90}$", "MPIW$_{90}$"], cap, "tab:surr", colspec="lllrrrrr")


def numbers() -> str:
    """LaTeX macros for numbers quoted in the text (all read from result files)."""
    S = json.loads(cfg.uq_summary.read_text()); reg = pd.read_csv(cfg.regimes_csv, index_col=0); m = json.loads(cfg.pareto_meta.read_text()); sc = pd.read_csv(cfg.scenarios_summary).set_index("scenario")
    scm = json.loads(cfg.scenarios_meta.read_text()); sob = json.loads(cfg.sobol_json.read_text()); M = json.loads(cfg.surrogate_metrics.read_text()); st = json.loads(cfg.steam_credit_json.read_text())
    L = yaml.safe_load((LV / "latham_fit.yaml").read_text()); X = yaml.safe_load((LV / "xf1989_fit.yaml").read_text()); lhs = json.loads((C.DATA / "lhs_runs" / f"{cfg.lhs_name}_meta.json").read_text())
    rt = json.loads(cfg.runtimes_json.read_text())
    cal = L["calibration_alpha_top_fixed"]["best_fit"]; camp = scm["S5_campaigns"]
    v = X["verification_matched_alpha"]; rec = v["runs"][v["recommended_configuration"]]
    mac = {
        "nLHS": lhs["n_samples"], "nSaltelli": sob["n_runs"], "nMC": S["N"], "speedup": f"{M['speed']['speedup_factor']:.0f}", "surrRMSETwo": f"{M['targets']['T_wo_max_K']['models'][M['targets']['T_wo_max_K']['best']]['test']['rmse']:.2f}",
        "surrRMSElogtr": f"{M['targets']['log10_t_r_hot']['models'][M['targets']['log10_t_r_hot']['best']]['test']['rmse']:.4f}", "PICPtwo": f"{M['conformal']['T_wo_max_K']['levels']['0.9']['picp']:.3f}",
        "calFgt": f"{cal['params']['F_gt']:.3f}", "calLq": f"{cal['params']['L_q']:.3f}", "calFhtg": f"{cal['params']['f_htg']:.2f}", "calChi": f"{cal['chi2']:.1f}", "calRedChi": f"{cal['covariance']['reduced_chi2']:.2f}",
        "alphaIplant": f"{L['calibration']['best_fit']['alpha_i_plant_A_midheight_W_m2K']:.0f}", "alphaIxf": f"{v['alpha_i_back_calculated']['median_W_m2K']:.0f}", "xfRMSET": f"{rec['metrics']['T_gas']['rmse']:.1f}", "xfRMSEinc": f"{rec['x_CH4_increment']['rmse']:.3f}",
        "SobolFiringTwo": f"{sob['outputs']['T_wo_max_K']['ST']['specific_firing_factor']:.2f}", "SobolSCTwo": f"{sob['outputs']['T_wo_max_K']['ST']['steam_to_carbon']:.2f}", "SobolEAtwo": f"{sob['outputs']['T_wo_max_K']['ST']['excess_air']:.2f}",
        "SobolPlogtr": f"{sob['outputs']['log10_t_r_hot']['ST']['inlet_P']:.2f}", "SobolFiringSlip": f"{sob['outputs']['CH4_slip_dry_pct']['ST']['specific_firing_factor']:.2f}",
        "SxDaily": f"{sc.loc['S2_daily', 'life_consumption_per_kmol_H2_rel_S1']:.3f}", "SxRenew": f"{sc.loc['S3_renewable', 'life_consumption_per_kmol_H2_rel_S1']:.3f}", "SxSClow": f"{sc.loc['S6_SC2.5', 'life_consumption_per_kmol_H2_rel_S1']:.2f}", "SxSChigh": f"{sc.loc['S6_SC3.5', 'life_consumption_per_kmol_H2_rel_S1']:.2f}",
        "SxAgeSlip": f"{sc.loc['S5_hold_CH4_slip_ageing_y4', 'life_consumption_per_kmol_H2_rel_S1']:.2f}", "SxAgeTout": f"{sc.loc['S5_hold_T_out_ageing_y4', 'life_consumption_per_kmol_H2_rel_S1']:.2f}", "SxAgeSlipYone": f"{sc.loc['S5_hold_CH4_slip_ageing_y1', 'life_consumption_per_kmol_H2_rel_S1']:.2f}",
        "TwoAgeStart": f"{camp['hold_CH4_slip_ageing'][0]['T_wo_start']:.0f}", "TwoAgeEnd": f"{camp['hold_CH4_slip_ageing'][-1]['T_wo_end']:.0f}", "SxAvgErr": f"{sc.loc['S3_renewable', 'avg_condition_error']:.3f}",
        "evMild": f"{scm['S4_events']['S4a_mild_overfire']['ratio_event_to_base']:.1f}", "evSevere": f"{scm['S4_events']['S4b_severe_overfire']['ratio_event_to_base']:.1f}", "evHot": f"{scm['S4_events']['S4c_hot_band']['ratio_event_to_base']:.1f}", "evHotShare": f"{scm['S4_events']['S4c_hot_band']['share_of_annual_damage_pct']:.1f}",
        "nPareto": m["n_pareto"], "nFeasible": m["n_feasible_phys"], "nDominating": m["base_domination"]["n_dominating"], "kneeHtwo": f"{100*(reg.loc['knee','H2_net_kmol_h_phys']/reg.loc['base','H2_net_kmol_h_phys']-1):+.0f}", "kneeLife": f"{reg.loc['knee','life_rate_per_kmol_H2_rel_base_phys']:.2f}",
        "kneeFuel": f"{100*(reg.loc['knee','fuel_MJ_per_kmol_H2_phys']/reg.loc['base','fuel_MJ_per_kmol_H2_phys']-1):+.0f}", "isoLife": f"{reg.loc['min_life_iso_H2','life_rate_per_kmol_H2_rel_base_phys']:.3f}", "isoFuel": f"{100*(reg.loc['min_life_iso_H2','fuel_MJ_per_kmol_H2_phys']/reg.loc['base','fuel_MJ_per_kmol_H2_phys']-1):+.0f}",
        "isoTwo": f"{reg.loc['min_life_iso_H2','T_wo_max_K_phys']:.0f}", "baseTwo": f"{reg.loc['base','T_wo_max_K_phys']:.0f}", "maxHtwoLife": f"{reg.loc['max_H2','life_rate_per_kmol_H2_rel_base_phys']:.1f}", "maxHtwoGain": f"{100*(reg.loc['max_H2','H2_net_kmol_h_phys']/reg.loc['base','H2_net_kmol_h_phys']-1):+.0f}",
        "steamMinLifeZero": f"{st['per_alpha']['0.0']['min_life']['life_rel_base']:.3f}", "steamMinLifeZeroFuel": f"{st['per_alpha']['0.0']['min_life']['fuel_total_rel_base_pct']:+.1f}", "steamMinFuelZeroSC": f"{st['per_alpha']['0.0']['min_fuel']['steam_to_carbon']:.2f}", "steamMinFuelOneSC": f"{st['per_alpha']['1.0']['min_fuel']['steam_to_carbon']:.2f}",
        "uqTwoWidth": f"{S['intervals']['S1_base']['T_wo_max_K']['width_90']:.0f}", "uqLogtrLo": f"{S['intervals']['S1_base']['log10_t_r']['p05']:.2f}", "uqLogtrHi": f"{S['intervals']['S1_base']['log10_t_r']['p95']:.2f}", "uqKneeLo": f"{S['intervals']['RQ3_knee']['life_rate_per_kmol_H2_rel_base']['p05']:.2f}", "uqKneeHi": f"{S['intervals']['RQ3_knee']['life_rate_per_kmol_H2_rel_base']['p95']:.2f}",
        "uqAgeLo": f"{S['intervals']['S5_y4_holdslip_end']['life_rate_per_kmol_H2_rel_base']['p05']:.2f}", "uqAgeHi": f"{S['intervals']['S5_y4_holdslip_end']['life_rate_per_kmol_H2_rel_base']['p95']:.2f}",
        "shareCreep": f"{100*S['variance_shares_S1_base']['log10_t_r']['shares']['creep']:.0f}", "shareKinTwo": f"{100*S['variance_shares_S1_base']['T_wo_max_K']['shares']['kinetics']:.0f}", "shareHT": f"{100*S['variance_shares_S1_base']['log10_t_r']['shares']['heat_transfer']:.0f}",
        "robAge": f"{S['robustness_fractions']['nominal_L_q']['S5y4_life_gt_2x_S1']:.2f}", "robKnee": f"{S['robustness_fractions']['nominal_L_q']['knee_life_lt_0.5_base']:.2f}", "robSC": f"{S['robustness_fractions']['nominal_L_q']['SC3.5_less_life_per_kmol_than_SC2.5']:.2f}",
        "robFlexNom": f"{S['robustness_fractions']['nominal_L_q']['S3_life_per_kmol_le_S1']:.2f}", "robFlexVar": f"{S['robustness_fractions']['load_dependent_L_q']['S3_life_per_kmol_le_S1']:.3f}", "lqvarLoadLow": f"{S['intervals_lq_variant']['load_0.70']['life_rate_per_kmol_H2_rel_base']['median']:.2f}",
        "pipelineMin": f"{sum(v['wall_time_s'] for v in rt.values())/60:.0f}",
    }
    lines = ["% auto-generated by rdt.paper_tables.numbers() from v2 result files; do not edit"]
    for k, v in mac.items():
        lines.append(f"\\newcommand{{\\{k}}}{{{v}}}")
    return "\n".join(lines) + "\n"


def make_all() -> Dict[str, Path]:
    TAB.mkdir(exist_ok=True); out = {}
    for name, fn in (("tab_xf_verification", table_xf_verification), ("tab_latham_cases", table_latham_cases), ("tab_calibration", table_calibration), ("tab_sobol", table_sobol), ("tab_scenarios", table_scenarios),
                     ("tab_regimes", table_regimes), ("tab_steam_credit", table_steam), ("tab_uq", table_uq), ("tab_surrogate", table_surrogate)):
        p = TAB / f"{name}.tex"; p.write_text(fn()); out[name] = p
    NUMBERS.write_text(numbers()); out["numbers"] = NUMBERS
    return out
