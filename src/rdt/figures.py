"""Journal-style regeneration of all paper figures from the results of a :class:`rdt.config.RunConfig`.

``make_all(cfg)`` archives the existing PNGs in ``paper/figures/v1/`` (once) and writes ``figNN_*.pdf`` and ``.png``
using :mod:`rdt.plotting`. Figures 00-06 are physics/validation figures (recomputed, independent of the campaign
version); figures 07-15 read the ``*_<tag>`` result files. No in-figure titles: descriptions belong in captions.
"""

from __future__ import annotations

import json
import shutil
from dataclasses import replace
from pathlib import Path
from typing import Dict, List

import joblib
import numpy as np
import pandas as pd

from rdt import config as C
from rdt import plotting as P

FIG_DIR = C.ROOT / "paper" / "figures"
V1_DIR = FIG_DIR / "v1"


def archive_v1() -> None:
    V1_DIR.mkdir(exist_ok=True)
    for p in FIG_DIR.glob("fig*.png"):
        if not (V1_DIR / p.name).exists():
            shutil.copy2(p, V1_DIR / p.name)


# ---------------------------------------------------------------------------
# 00-06: physics and validation figures (version independent)
# ---------------------------------------------------------------------------
def fig00_equilibrium(stem):
    import cantera as ct
    gas = ct.Solution("gri30.yaml"); T = np.array([700.0, 750.0, 800.0, 850.0, 900.0]); conv = []
    for t in T:
        gas.TPX = t + 273.15, 25e5, {"CH4": 1.0, "H2O": 3.0}; mw0 = gas.mean_molecular_weight; gas.equilibrate("TP")
        conv.append(1.0 - gas.X[gas.species_index("CH4")] * 4.0 * mw0 / gas.mean_molecular_weight)
    with P.style():
        fig, ax = P.figure(); ax.plot(T, 100 * np.array(conv), marker="o"); P.tidy(ax, "Temperature (°C)", "Equilibrium CH$_4$ conversion (%)")
        return P.save(fig, stem)


def fig01_kinetics(stem):
    from rdt import kinetics as kin
    prm = kin.load_params(); T = np.linspace(573.0, 1123.0, 60); k1, k2, k3 = prm.arrhenius(T)
    feed = {"CH4": 1.0, "H2O": 3.0, "H2": 0.1}; tot = sum(feed.values()); p = {s: 25.0 * feed.get(s, 0.0) / tot for s in kin.SPECIES}
    Tr = np.linspace(700.0, 1100.0, 41); r1, r2, r3 = kin.rates(Tr, p)
    with P.style():
        fig, ax = P.figure(P.DOUBLE, 6.5 * P.CM, 1, 2)
        ax[0].semilogy(1000 / T, k1, label="$k_1$"); ax[0].semilogy(1000 / T, k2, label="$k_2$"); ax[0].semilogy(1000 / T, k3, label="$k_3$")
        ax[0].semilogy(1000 / np.array([773.0, 798.0, 823.0]), [0.2088, 0.5254, 2.069], "ko", label="$k_1$, Table 4"); P.tidy(ax[0], "1000/T (K$^{-1}$)", "Rate coefficient (kmol kg$^{-1}$ h$^{-1}$ bar$^n$)", True)
        ax[1].plot(Tr, r1, label="$r_1$"); ax[1].plot(Tr, r2, label="$r_2$"); ax[1].plot(Tr, r3, label="$r_3$"); P.tidy(ax[1], "T (K)", "Rate (kmol kg$^{-1}$ h$^{-1}$)", True)
        P.panel_label(ax[0], "a"); P.panel_label(ax[1], "b"); fig.tight_layout(); return P.save(fig, stem)


def fig02_reactor_profiles(stem):
    from rdt import reactor1d as r1
    tube = r1.tube_from_xu_froment(); bed = r1.bed_from_xu_froment(); feed = r1.feed_from_xu_froment()
    res = r1.simulate(tube, bed, feed, r1.WallBC.linear(1000.0, 1150.0, tube.L_heated))
    with P.style():
        fig, ax = P.figure(P.DOUBLE, 5.5 * P.CM, 1, 3)
        ax[0].plot(res.z, res.conversion_CH4, label="CH$_4$ conversion"); ax[0].plot(res.z, res.yield_CO2, label="CO$_2$ yield"); P.tidy(ax[0], "z (m)", "Conversion (–)", True)
        ax[1].plot(res.z, res.T_wall_outer, "k", label="outer wall"); ax[1].plot(res.z, res.T_wall_inner, label="inner wall"); ax[1].plot(res.z, res.T, label="process gas"); P.tidy(ax[1], "z (m)", "T (K)", True)
        ax[2].plot(res.z, res.P, "k"); P.tidy(ax[2], "z (m)", "p (bar)")
        for a, l in zip(ax, "abc"): P.panel_label(a, l)
        fig.tight_layout(); return P.save(fig, stem)


def fig03_validation(stem):
    import yaml
    from rdt import validate_xf1989 as vx
    tidy = vx.load_tidy(); cv = vx.curves(tidy); wall = vx.wall_from_digitised(tidy)
    d = yaml.safe_load((C.DATA / "literature_validation" / "xf1989_fit.yaml").read_text())["verification_matched_alpha"]
    alpha = d["alpha_i_back_calculated"]["median_W_m2K"]
    rec = vx.run_matched(alpha, wall, latham_inlet=True); lg = vx.run_case("leva_grummer", "latham", wall)
    with P.style():
        fig, ax = P.figure(P.DOUBLE, 11 * P.CM, 2, 2); ax = ax.ravel(); mk = dict(ls="none", marker="o", ms=2, alpha=0.5); sub = 8
        for c, col in (("x_CH4", P.PALETTE[0]), ("x_CO2", P.PALETTE[1])):
            z, v = cv[c]; ax[0].plot(z[::sub], v[::sub], color=col, **mk); ax[0].plot(rec.z, vx.model_curve(rec, c, rec.z), color=col, label=f"{c.replace('x_', '')} matched $\\alpha_i$"); ax[0].plot(lg.z, vx.model_curve(lg, c, lg.z), color=col, ls=":", label=f"{c.replace('x_', '')} Leva–Grummer")
        P.tidy(ax[0], "z (m)", "Conversion (–)", True)
        z, v = cv["T_wall_outer"]; ax[1].plot(z[::sub], v[::sub], color="k", **mk)
        for c, col, lab in (("T_wall_inner", P.PALETTE[1], "inner wall"), ("T_gas", P.PALETTE[2], "process gas")):
            z, v = cv[c]; ax[1].plot(z[::sub], v[::sub], color=col, **mk); ax[1].plot(rec.z, vx.model_curve(rec, c, rec.z), color=col, label=lab); ax[1].plot(lg.z, vx.model_curve(lg, c, lg.z), color=col, ls=":")
        P.tidy(ax[1], "z (m)", "T (K)", True)
        z, v = cv["p_t"]; ax[2].plot(z[::sub], v[::sub], color="k", **mk); ax[2].plot(rec.z, rec.P, "k", label="matched $\\alpha_i$"); ax[2].plot(lg.z, lg.P, "k:", label="Leva–Grummer"); P.tidy(ax[2], "z (m)", "p (bar)", True)
        for res_, ls, lab in ((rec, "-", "matched"), (lg, ":", "Leva–Grummer")):
            z, v = cv["T_gas"]; ax[3].plot(z, vx.model_curve(res_, "T_gas", z) - v, color=P.PALETTE[2], ls=ls, label=f"gas, {lab}")
            z, v = cv["T_wall_inner"]; ax[3].plot(z, vx.model_curve(res_, "T_wall_inner", z) - v, color=P.PALETTE[1], ls=ls, label=f"inner wall, {lab}")
        ax[3].axhline(0, color="k", lw=0.5); P.tidy(ax[3], "z (m)", "Model − digitised (K)", True)
        for a, l in zip(ax, "abcd"): P.panel_label(a, l)
        fig.tight_layout(); return P.save(fig, stem)


def _latham_profiles(stem, params, cases, parity=False):
    from rdt import latham_cases as lc
    df = lc.load_cases(); n = len(cases)
    with P.style():
        if parity:
            fig, axes = P.figure(P.DOUBLE, 11 * P.CM, 2, 3); axes_p = axes[0]; axp = axes[1]
        else:
            fig, axes = P.figure(P.DOUBLE, 11 * P.CM, 2, 2); axes_p = axes.ravel()
        pred_m = {}
        for ax, case in zip(axes_p, cases):
            row = df[df.case == case].iloc[0]; res = lc.run_case(row, params); c = lc.compare(row, res); zf = res.z_frac
            ax.plot(zf, res.T_fg, color=P.PALETTE[1], label="furnace gas"); ax.plot(zf, res.T_wo, "k", label="outer wall"); ax.plot(zf, res.T_wi, color=P.PALETTE[4], label="inner wall"); ax.plot(zf, res.tube.T, color=P.PALETTE[0], label="process gas")
            ax.errorbar([row.TWT_upper_frac, row.TWT_lower_frac], [row.TWT_upper_K, row.TWT_lower_K], yerr=1.96 * row.sigma_TWT_K, fmt="ks", ms=3, capsize=2, label="measured wall T")
            ax.errorbar([1.0], [row.out_T_K], yerr=1.96 * row.sigma_out_T_K, fmt="o", color=P.PALETTE[0], ms=3, capsize=2); ax.errorbar([1.0], [row.flue_out_T_K], yerr=1.96 * row.sigma_flue_out_T_K, fmt="o", color=P.PALETTE[1], ms=3, capsize=2)
            ax.text(0.02, 0.97, f"{case.replace('_', ' ')} ({int(row.plant_rate_pct)} %)", transform=ax.transAxes, va="top", fontsize=7); P.tidy(ax, "z/L (–)", "T (K)")
            pred_m[case] = c
        axes_p[0].legend(loc="lower right")
        if parity:
            for a in axp: a.axis("off")
            axp = fig.add_subplot(2, 1, 2)
            for q, col, lab in (("out_T_K", P.PALETTE[0], "outlet T"), ("flue_out_T_K", P.PALETTE[1], "flue-gas outlet T"), ("TWT_upper_K", "k", "wall T, upper"), ("TWT_lower_K", "grey", "wall T, lower")):
                axp.plot([pred_m[c][q]["measured"] for c in cases], [pred_m[c][q]["predicted"] for c in cases], "o", color=col, label=lab)
            lo, hi = 1050, 1320; axp.plot([lo, hi], [lo, hi], "k-", lw=0.5); axp.plot([lo, hi], [lo + 10, hi + 10], "k:", lw=0.4); axp.plot([lo, hi], [lo - 10, hi - 10], "k:", lw=0.4)
            P.tidy(axp, "Measured T (K)", "Predicted T (K)", True); axp.set_aspect("equal")
        fig.tight_layout(); return P.save(fig, stem)


def fig04_latham_baseline(stem):
    from rdt import latham_cases as lc
    return _latham_profiles(stem, lc.BASELINE, ["Plant_A", "Plant_B", "Plant_C1", "Plant_C2"])


def fig05_latham_calibrated(stem):
    from rdt import operate as op
    return _latham_profiles(stem, op.calibrated_params(), ["Plant_A", "Plant_B", "Plant_C2"], parity=True)


def fig06_creep_profile(stem):
    from rdt import creep, latham_cases as lc, operate as op
    p = op.calibrated_params(); df = lc.load_cases(); row = df[df.case == "Plant_A"].iloc[0]; res = lc.run_case(row, p); curve = creep.LarsonMillerCurve.from_config()
    T_mid = 0.5 * (res.T_wo + res.T_wi); lo = creep.life_along_tube(res.z, res.T_wo, res.tube.P, p.tube_od_m, p.wall_thickness_m, curve); lm = creep.life_along_tube(res.z, T_mid, res.tube.P, p.tube_od_m, p.wall_thickness_m, curve)
    with P.style():
        fig, ax = P.figure(P.SINGLE, 7 * P.CM); zf = res.z_frac
        ax.plot(zf, res.T_wo, "k", label="outer wall"); ax.plot(zf, T_mid, color=P.PALETTE[1], label="mid-wall"); P.tidy(ax, "z/L (–)", "T (K)")
        ax2 = ax.twinx(); ax2.plot(zf, lo["log10_t_r"], color=P.PALETTE[0], ls="--", label="log$_{10}$ $t_r$, outer"); ax2.plot(zf, lm["log10_t_r"], color=P.PALETTE[0], ls=":", label="log$_{10}$ $t_r$, mid"); ax2.set_ylabel("log$_{10}$ $t_r$ (h)"); ax2.spines["right"].set_visible(True)
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1 + h2, l1 + l2, loc="lower right")
        fig.tight_layout(); return P.save(fig, stem)


# ---------------------------------------------------------------------------
# 07-15: campaign figures (tagged results)
# ---------------------------------------------------------------------------
def fig07_maps(stem, cfg):
    from rdt import operate as op
    g = pd.read_csv(C.DATA / "lhs_runs" / f"{cfg.grid_name}.csv.gz"); g = g[g.converged == True]; base = op.base_case()
    sc = np.sort(g.steam_to_carbon.unique()); ld = np.sort(g.feed_per_tube_fraction.unique())
    piv = lambda col: g.pivot(index="steam_to_carbon", columns="feed_per_tube_fraction", values=col).reindex(index=sc, columns=ld).to_numpy()
    LD, SC = np.meshgrid(ld, sc)
    with P.style():
        fig, ax = P.figure(P.DOUBLE, 5.5 * P.CM, 1, 3)
        for a, (Z, lab, cmap) in zip(ax, ((piv("T_wo_max_K"), "Max. outer-wall T (K)", "magma"), (np.log10(piv("life_consumption_ratio_to_base")), "log$_{10}$ relative life consumption (–)", "viridis"), (piv("H2_net_kmol_h"), "H$_2$ per tube (kmol h$^{-1}$)", "cividis"))):
            cs = a.contourf(LD, SC, Z, levels=14, cmap=cmap); fig.colorbar(cs, ax=a, label=lab, pad=0.02); a.plot([1.0], [base.steam_to_carbon], "w*", ms=8, mec="k"); P.tidy(a, "Load (–)", "Steam-to-carbon (–)")
        cs2 = ax[2].contour(LD, SC, piv("specific_firing_factor"), levels=[0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15], colors="w", linewidths=0.5); ax[2].clabel(cs2, fmt="%.2f", fontsize=6)
        for a, l in zip(ax, "abc"): P.panel_label(a, l)
        fig.tight_layout(); return P.save(fig, stem)


def fig08_sobol(stem, cfg):
    from rdt import operate as op
    sob = json.loads(cfg.sobol_json.read_text()); names = list(op.INPUT_NAMES); x = np.arange(len(names)); w = 0.38
    short = {"steam_to_carbon": "S/C", "feed_per_tube_fraction": "load", "specific_firing_factor": "firing", "inlet_T": "$T_{in}$", "inlet_P": "$p_{in}$", "catalyst_activity": "activity", "excess_air": "excess air"}
    with P.style():
        fig, ax = P.figure(P.DOUBLE, 5.5 * P.CM, 1, 2)
        for a, y in zip(ax, ("T_wo_max_K", "log10_t_r_hot")):
            s = sob["outputs"][y]
            a.bar(x - w / 2, [s["S1"][k] for k in names], w, yerr=[s["S1_conf"][k] for k in names], capsize=1.5, label="$S_1$"); a.bar(x + w / 2, [s["ST"][k] for k in names], w, yerr=[s["ST_conf"][k] for k in names], capsize=1.5, label="$S_T$")
            a.set_xticks(x); a.set_xticklabels([short[k] for k in names], rotation=40, ha="right"); P.tidy(a, None, "Sobol index (–)", True)
        P.panel_label(ax[0], "a"); P.panel_label(ax[1], "b"); fig.tight_layout(); return P.save(fig, stem)


def fig09_parity(stem, cfg):
    from rdt import surrogate as s
    M = json.loads(cfg.surrogate_metrics.read_text()); te = pd.read_csv(cfg.saltelli_file); te = te[te.converged == True].reset_index(drop=True); te = s.add_physics_features(te)
    with P.style():
        fig, ax = P.figure(P.DOUBLE, 5.5 * P.CM, 1, 3)
        for a, (tgt, unit) in zip(ax, (("T_wo_max_K", "K"), ("CH4_slip_dry_pct", "%"), ("log10_t_r_hot", "–"))):
            b = joblib.load(cfg.models_dir / f"{tgt}__conformal.joblib"); X = te[b["features"]].to_numpy(float); y = te[tgt].to_numpy(float); yp, iv = b["conformal"].predict_interval(X)
            sub = np.argsort(y)[:: max(1, len(y) // 300)]
            a.errorbar(y[sub], yp[sub], yerr=[yp[sub] - iv[sub, 0, 0], iv[sub, 1, 0] - yp[sub]], fmt="o", ms=2, alpha=0.5, elinewidth=0.4)
            lim = [y.min(), y.max()]; a.plot(lim, lim, "k-", lw=0.5); P.tidy(a, f"Physics model ({unit})", f"Surrogate ({unit})"); a.set_aspect("equal")
        for a, l in zip(ax, "abc"): P.panel_label(a, l)
        fig.tight_layout(); return P.save(fig, stem)


def fig10_timeseries(stem, cfg):
    h = pd.read_csv(cfg.scenario_hourly("S3_renewable")); wk = slice(24 * 7 * 20, 24 * 7 * 21); t = np.arange(len(h))[wk] / 24.0
    with P.style():
        fig, ax = P.figure(P.SINGLE, 12 * P.CM, 4, 1, sharex=True)
        for a, (col, lab, c) in zip(ax, (("feed_per_tube_fraction", "Load (–)", 0), ("specific_firing_factor", "Firing factor (–)", 1), ("T_wo_max_K", "Max. wall T (K)", 7), ("dD", "Damage rate (10$^{-9}$ h$^{-1}$)", 3))):
            a.plot(t, (1e9 if col == "dD" else 1) * h[col][wk], color=P.PALETTE[c]); P.tidy(a, None, lab)
        ax[-1].set_xlabel("Day of year"); fig.tight_layout(h_pad=0.3); return P.save(fig, stem)


def fig11_scenario_summary(stem, cfg):
    summ = pd.read_csv(cfg.scenarios_summary).set_index("scenario")
    order = [o for o in ["S1_steady", "S2_daily", "S3_renewable", "S4a_mild_overfire", "S4b_severe_overfire", "S4c_hot_band", "S6_SC2.5", "S6_SC3.5", "S5_hold_CH4_slip_ageing_y1", "S5_hold_CH4_slip_ageing_y4", "S5_hold_T_out_ageing_y1", "S5_hold_T_out_ageing_y4"] if o in summ.index]
    labels = {"S1_steady": "S1 steady", "S2_daily": "S2 daily", "S3_renewable": "S3 renewable", "S4a_mild_overfire": "S4a mild", "S4b_severe_overfire": "S4b severe", "S4c_hot_band": "S4c hot band", "S6_SC2.5": "S6 S/C 2.5", "S6_SC3.5": "S6 S/C 3.5",
              "S5_hold_CH4_slip_ageing_y1": "S5 slip y1", "S5_hold_CH4_slip_ageing_y4": "S5 slip y4", "S5_hold_T_out_ageing_y1": "S5 $T_{out}$ y1", "S5_hold_T_out_ageing_y4": "S5 $T_{out}$ y4"}
    D1, H1 = summ.loc["S1_steady", "annual_damage_D"], summ.loc["S1_steady", "annual_H2_kmol_per_tube"]
    rel = summ.loc[order, "life_consumption_per_kmol_H2_rel_S1"]; rel_avg = (summ.loc[order, "D_avg_condition"] / summ.loc[order, "annual_H2_kmol_per_tube"]) / (D1 / H1)
    with P.style():
        fig, ax = P.figure(P.DOUBLE, 6 * P.CM); x = np.arange(len(order))
        ax.bar(x, rel, color=P.PALETTE[0], label="hourly accumulation"); ax.plot(x, rel_avg, "D", color="k", ms=3, label="annual-mean inputs"); ax.axhline(1.0, color="grey", lw=0.5)
        ax.set_xticks(x); ax.set_xticklabels([labels[o] for o in order], rotation=40, ha="right"); ax.set_yscale("log"); P.tidy(ax, None, "Life consumption per kmol H$_2$, rel. S1 (–)", True)
        fig.tight_layout(); return P.save(fig, stem)


def fig12_pareto(stem, cfg):
    par = pd.read_csv(cfg.pareto_csv); reg = pd.read_csv(cfg.regimes_csv, index_col=0); feas = par[par.feasible_phys]; infeas = par[~par.feasible_phys]
    marks = {"base": ("*", "k", 9), "min_life_iso_H2": ("P", P.PALETTE[2], 6), "min_fuel_overall": ("v", P.PALETTE[1], 5), "knee": ("D", P.PALETTE[3], 5), "max_H2": ("^", P.PALETTE[4], 5)}
    def panel(a, xk, yk, xl, yl):
        a.scatter(feas[xk], feas[yk], s=6, color=P.PALETTE[0], alpha=0.7)
        if len(infeas): a.scatter(infeas[xk], infeas[yk], s=8, facecolors="none", edgecolors="grey")
        for n, (m, c, s_) in marks.items():
            if n in reg.index: a.plot(reg.loc[n, xk], reg.loc[n, yk], marker=m, color=c, ms=s_, ls="none", mec="k", mew=0.4, label=n.replace("_", " "))
        P.tidy(a, xl, yl)
    with P.style():
        fig, ax = P.figure(P.DOUBLE, 5.5 * P.CM, 1, 3)
        panel(ax[0], "fuel_MJ_per_kmol_H2_phys", "life_rate_per_kmol_H2_rel_base_phys", "Fuel (MJ kmol$^{-1}$ H$_2$)", "Life consumption per kmol H$_2$ (rel.)"); ax[0].set_yscale("log")
        panel(ax[1], "H2_net_kmol_h_phys", "life_rate_per_kmol_H2_rel_base_phys", "H$_2$ per tube (kmol h$^{-1}$)", "Life consumption per kmol H$_2$ (rel.)"); ax[1].set_yscale("log")
        panel(ax[2], "H2_net_kmol_h_phys", "fuel_MJ_per_kmol_H2_phys", "H$_2$ per tube (kmol h$^{-1}$)", "Fuel (MJ kmol$^{-1}$ H$_2$)"); ax[2].legend(loc="upper left")
        for a, l in zip(ax, "abc"): P.panel_label(a, l)
        fig.tight_layout(); return P.save(fig, stem)


def fig13_steam_credit(stem, cfg):
    sc_ = json.loads(cfg.steam_credit_json.read_text()); rows = []
    for a in sc_["alphas"]:
        r = sc_["per_alpha"][str(float(a))]
        for lab in ("min_life", "knee", "min_fuel"):
            if lab in r: rows.append({"alpha": float(a), "regime": lab, **r[lab]})
    tab = pd.DataFrame(rows)
    with P.style():
        fig, ax = P.figure(P.DOUBLE, 5.5 * P.CM, 1, 2)
        for lab, mk in (("min_life", "P"), ("knee", "D"), ("min_fuel", "v")):
            t = tab[tab.regime == lab]; ax[0].plot(t.alpha, t.steam_to_carbon, marker=mk, label=lab.replace("_", " ")); ax[1].plot(t.alpha, t.life_rel_base, marker=mk, label=lab.replace("_", " "))
        P.tidy(ax[0], "Steam credit $\\alpha$ (–)", "Optimal steam-to-carbon (–)", True); ax[0].set_ylim(2.4, 4.1)
        ax[1].set_yscale("log"); P.tidy(ax[1], "Steam credit $\\alpha$ (–)", "Life consumption per kmol H$_2$ (rel.)", True)
        P.panel_label(ax[0], "a"); P.panel_label(ax[1], "b"); fig.tight_layout(); return P.save(fig, stem)


def fig14_uq(stem, cfg):
    names = {"S1_base": "S1 base", "RQ3_knee": "knee", "RQ3_min_life_iso_H2": "min-life iso-H$_2$", "S5_y4_holdslip_end": "S5 year 4"}
    with P.style():
        fig, ax = P.figure(P.DOUBLE, 5.5 * P.CM, 1, 2)
        for n, lab in names.items():
            d = pd.read_csv(cfg.uq_dir / f"{cfg.mc_tag}_{n}.csv.gz"); d = d[d.converged == True]
            ax[0].hist(d.T_wo_max_K, bins=50, histtype="step", density=True, label=lab); ax[1].hist(d.log10_t_r, bins=50, histtype="step", density=True, label=lab)
        ax[0].axvline(1193, color="k", ls=":", lw=0.6); P.tidy(ax[0], "Max. outer-wall T (K)", "Density (K$^{-1}$)", True); P.tidy(ax[1], "log$_{10}$ $t_r$ at hot spot (h)", "Density (–)", True)
        P.panel_label(ax[0], "a"); P.panel_label(ax[1], "b"); fig.tight_layout(); return P.save(fig, stem)


def fig15_variance(stem, cfg):
    S = json.loads(cfg.uq_summary.read_text()); vs = S["variance_shares_S1_base"]; cols = ["log10_t_r", "T_wo_max_K"]; labels = ["log$_{10}$ $t_r$", "$T_{wo,max}$"]
    with P.style():
        fig, ax = P.figure(P.SINGLE, 6.5 * P.CM); bottom = np.zeros(2)
        for g in ("kinetics", "heat_transfer", "creep", "measurement"):
            v = np.array([vs[c]["shares"][g] for c in cols]); ax.bar(labels, v, bottom=bottom, label=g.replace("_", " ")); bottom += v
        ax.axhline(1.0, color="k", lw=0.5); P.tidy(ax, None, "Share of variance (–)", True); fig.tight_layout(); return P.save(fig, stem)


PHYSICS_FIGS = {"fig00_equilibrium_conversion": fig00_equilibrium, "fig01_kinetics": fig01_kinetics, "fig02_reactor_profiles": fig02_reactor_profiles, "fig03_validation_xf1989": fig03_validation,
                "fig04_latham_baseline": fig04_latham_baseline, "fig05_latham_calibrated": fig05_latham_calibrated, "fig06_creep_life_profile": fig06_creep_profile}
CAMPAIGN_FIGS = {"fig07_operating_maps": fig07_maps, "fig08_sobol": fig08_sobol, "fig09_surrogate_parity": fig09_parity, "fig10_scenarios_timeseries": fig10_timeseries, "fig11_scenarios_summary": fig11_scenario_summary,
                 "fig12_pareto": fig12_pareto, "fig13_steam_credit": fig13_steam_credit, "fig14_uq_distributions": fig14_uq, "fig15_variance_shares": fig15_variance}


def make_all(cfg: C.RunConfig = C.V2, only=None, log=print) -> List[Path]:
    archive_v1(); made = []
    for name, fn in {**PHYSICS_FIGS, **CAMPAIGN_FIGS}.items():
        if only and name not in only: continue
        try:
            paths = fn(FIG_DIR / name) if name in PHYSICS_FIGS else fn(FIG_DIR / name, cfg)
            made.extend(paths); log(f"{name}: ok")
        except Exception as e:  # noqa: BLE001
            log(f"{name}: FAILED {type(e).__name__}: {e}")
    return made
