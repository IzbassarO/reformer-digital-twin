"""Change report v1 -> v2 (headline numbers) and the data/README.md inventory."""

from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from rdt import config as C

ROOT, DATA = C.ROOT, C.DATA
EXPLAIN = {
    "sobol": "excess-air range widened from 5-20 % to 5-25 %, which changes the sampled space and therefore every variance share",
    "scenarios": "scenario base excess air is now the calibrated 21.6 % instead of 20 %, changing the flue-gas flow and hot-spot temperature of every scenario",
    "pareto": "decision bound on excess air raised to 25 % and base point moved to 21.6 % excess air, which shifts the base objectives that the regimes are compared with",
    "uq": "operating points inherit the v2 scenario/pareto base (21.6 % excess air) and the Monte Carlo uses N = 2048 instead of 4096 samples",
    "robustness": "paired fractions inherit the v2 operating points and the smaller Monte Carlo sample (N = 2048)",
    "steam": "base point moved to 21.6 % excess air and the excess-air decision bound rose to 25 %, changing the base fuel and the feasible front",
    "surrogate": "training and test campaigns resampled over the widened excess-air range",
}


#: item-specific explanations (substring match on the item name, first hit wins) override the family default
ITEM_EXPLAIN = [
    ("ST excess_air", "excess air now spans 5-25 % instead of 5-20 %: a wider input range carries more variance, so its total-effect index rises while all others fall slightly"),
    ("min_life_iso_H2 life", "the minimum-life regime now uses the excess-air headroom above 20 % (leaner, cooler flame), lowering its hot spot by 13 K; the base denominator also lengthened by ~22 %"),
    ("min_life_iso_H2 T_wo", "the minimum-life regime now uses excess air above 20 %, which dilutes the flame and lowers the hot-spot temperature"),
    ("knee life", "the knee point is selected by distance to the utopia point and moved along the front; its own hot spot is 5 K hotter and the base denominator lengthened by ~22 %"),
    ("max_H2 life", "the regime itself is unchanged (same hot spot); the ratio rises because the v2 base at 21.6 % excess air is 3.3 K cooler and lives ~22 % longer"),
    ("min_fuel_overall life", "the regime itself is unchanged (same hot spot); the ratio rises because the v2 base at 21.6 % excess air is 3.3 K cooler and lives ~22 % longer"),
    ("members dominating", "the v2 base point (21.6 % excess air, cooler hot spot) is a better reference, so fewer Pareto members dominate it"),
    ("RQ3_knee life per kmol", "inherits the v2 knee regime (moved along the front) and the ~22 % longer base life"),
    ("RQ3_knee log10", "the v2 knee regime sits at a 5 K hotter hot spot"),
    ("RQ3_min_life_iso_H2", "the v2 minimum-life regime uses excess air above 20 % and runs 13 K cooler at the hot spot"),
    ("knee_life_lt_0.5_base", "the v2 knee regime consumes 0.34 of base life instead of 0.23 (median), so fewer paired samples fall below the 0.5 threshold"),
    ("min_fuel S/C", "the minimum-fuel end of the iso-production front moved because the base fuel changed (21.6 % excess air) and the front is flat in S/C near its fuel optimum"),
    ("min_fuel life", "the minimum-fuel end of the iso-production front is a different member (flat fuel objective in S/C), and the base life denominator lengthened by ~22 %"),
    ("min_life life", "the v2 minimum-life points use excess air above 20 %, cooling the hot spot; the base denominator lengthened by ~22 %"),
    ("best test RMSE", "same 2000 training points spread over a 25 % wider excess-air range (lower point density); the surrogate remains far inside the test thresholds (3 K, 0.05)"),
]


def _num(x, nd=3):
    try:
        return f"{float(x):.{nd}g}"
    except Exception:  # noqa: BLE001
        return str(x)


def _flag(v1, v2, abs_tol: Optional[float] = None) -> Tuple[str, float]:
    try:
        v1f, v2f = float(v1), float(v2)
    except Exception:  # noqa: BLE001
        return ("", float("nan"))
    if abs_tol is not None:
        d = v2f - v1f; return ("**changed**" if abs(d) > abs_tol else "", d)
    rel = (v2f - v1f) / abs(v1f) if v1f != 0 else float("inf") if v2f != 0 else 0.0
    return ("**changed**" if abs(rel) > 0.05 else "", rel)


def collect(cfg1: C.RunConfig = C.V1, cfg2: C.RunConfig = C.V2) -> List[Dict[str, object]]:
    rows = []
    def add(family, item, v1, v2, abs_tol=None, unit=""):
        flag, d = _flag(v1, v2, abs_tol)
        expl = ""
        if flag:
            expl = next((e for key, e in ITEM_EXPLAIN if key in item), EXPLAIN[family])
        rows.append({"family": family, "item": item, "unit": unit, "v1": _num(v1), "v2": _num(v2), "change": (f"{d:+.2f} {unit}" if abs_tol is not None else f"{100*d:+.1f} %") if d == d else "", "flag": flag,
                     "explanation": expl})
    # Sobol
    s1, s2 = json.loads(cfg1.sobol_json.read_text()), json.loads(cfg2.sobol_json.read_text())
    for y in ("T_wo_max_K", "log10_t_r_hot"):
        for k in s1["outputs"][y]["ST"]:
            add("sobol", f"{y} ST {k}", s1["outputs"][y]["ST"][k], s2["outputs"][y]["ST"][k], abs_tol=0.05, unit="")
    # scenarios
    a, b = pd.read_csv(cfg1.scenarios_summary).set_index("scenario"), pd.read_csv(cfg2.scenarios_summary).set_index("scenario")
    for sc in ("S2_daily", "S3_renewable", "S4a_mild_overfire", "S4b_severe_overfire", "S4c_hot_band", "S6_SC2.5", "S6_SC3.5", "S5_hold_CH4_slip_ageing_y4", "S5_hold_T_out_ageing_y4"):
        if sc in a.index and sc in b.index:
            add("scenarios", f"{sc} life per kmol H2 rel. S1", a.loc[sc, "life_consumption_per_kmol_H2_rel_S1"], b.loc[sc, "life_consumption_per_kmol_H2_rel_S1"])
    add("scenarios", "S1 mean T_wo,max", a.loc["S1_steady", "T_wo_max_mean_K"], b.loc["S1_steady", "T_wo_max_mean_K"], abs_tol=5.0, unit="K")
    add("scenarios", "S3 avg-condition error", a.loc["S3_renewable", "avg_condition_error"], b.loc["S3_renewable", "avg_condition_error"])
    # pareto regimes
    r1, r2 = pd.read_csv(cfg1.regimes_csv, index_col=0), pd.read_csv(cfg2.regimes_csv, index_col=0)
    for reg in ("base", "min_life_iso_H2", "knee", "max_H2", "min_fuel_overall"):
        if reg in r1.index and reg in r2.index:
            add("pareto", f"{reg} H2 (kmol/h)", r1.loc[reg, "H2_net_kmol_h_phys"], r2.loc[reg, "H2_net_kmol_h_phys"])
            add("pareto", f"{reg} fuel (MJ/kmol)", r1.loc[reg, "fuel_MJ_per_kmol_H2_phys"], r2.loc[reg, "fuel_MJ_per_kmol_H2_phys"])
            add("pareto", f"{reg} life rel. base", r1.loc[reg, "life_rate_per_kmol_H2_rel_base_phys"], r2.loc[reg, "life_rate_per_kmol_H2_rel_base_phys"])
            add("pareto", f"{reg} T_wo,max", r1.loc[reg, "T_wo_max_K_phys"], r2.loc[reg, "T_wo_max_K_phys"], abs_tol=5.0, unit="K")
            add("pareto", f"{reg} S/C", r1.loc[reg, "steam_to_carbon"], r2.loc[reg, "steam_to_carbon"])
    m1, m2 = json.loads(cfg1.pareto_meta.read_text()), json.loads(cfg2.pareto_meta.read_text())
    add("pareto", "members dominating the base", m1["base_domination"]["n_dominating"], m2["base_domination"]["n_dominating"])
    # UQ intervals
    u1, u2 = json.loads(cfg1.uq_summary.read_text()), json.loads(cfg2.uq_summary.read_text())
    for pt in ("S1_base", "RQ3_knee", "RQ3_min_life_iso_H2", "S5_y4_holdslip_end"):
        i1, i2 = u1["intervals"][pt], u2["intervals"][pt]
        add("uq", f"{pt} T_wo,max median", i1["T_wo_max_K"]["median"], i2["T_wo_max_K"]["median"], abs_tol=5.0, unit="K")
        add("uq", f"{pt} T_wo,max 90 % width", i1["T_wo_max_K"]["width_90"], i2["T_wo_max_K"]["width_90"], abs_tol=5.0, unit="K")
        add("uq", f"{pt} log10 t_r median", i1["log10_t_r"]["median"], i2["log10_t_r"]["median"], abs_tol=0.1, unit="dec")
        add("uq", f"{pt} life per kmol rel. base median", i1["life_rate_per_kmol_H2_rel_base"]["median"], i2["life_rate_per_kmol_H2_rel_base"]["median"])
    for c in ("log10_t_r", "T_wo_max_K"):
        for g in ("kinetics", "heat_transfer", "creep", "measurement"):
            add("uq", f"variance share {c} {g}", u1["variance_shares_S1_base"][c]["shares"][g], u2["variance_shares_S1_base"][c]["shares"][g], abs_tol=0.05)
    for tag in ("nominal_L_q", "load_dependent_L_q"):
        for k, v in u1["robustness_fractions"][tag].items():
            add("robustness", f"{tag} {k}", v, u2["robustness_fractions"][tag].get(k, float("nan")), abs_tol=0.05)
    # steam credit
    c1, c2 = json.loads(cfg1.steam_credit_json.read_text()), json.loads(cfg2.steam_credit_json.read_text())
    for al in c1["alphas"]:
        for lab in ("min_life", "min_fuel"):
            p1, p2 = c1["per_alpha"][str(float(al))], c2["per_alpha"][str(float(al))]
            if lab in p1 and lab in p2:
                add("steam", f"alpha {al} {lab} S/C", p1[lab]["steam_to_carbon"], p2[lab]["steam_to_carbon"])
                add("steam", f"alpha {al} {lab} life rel. base", p1[lab]["life_rel_base"], p2[lab]["life_rel_base"])
    # surrogate
    g1, g2 = json.loads(cfg1.surrogate_metrics.read_text()), json.loads(cfg2.surrogate_metrics.read_text())
    for t in ("T_wo_max_K", "CH4_slip_dry_pct", "log10_t_r_hot"):
        add("surrogate", f"best test RMSE {t}", g1["targets"][t]["models"][g1["targets"][t]["best"]]["test"]["rmse"], g2["targets"][t]["models"][g2["targets"][t]["best"]]["test"]["rmse"])
    return rows


def write_diff_report(path: Path = DATA / "design" / "v1_vs_v2_diff.md") -> pd.DataFrame:
    rows = collect(); df = pd.DataFrame(rows)
    rt = json.loads(C.V2.runtimes_json.read_text()) if C.V2.runtimes_json.exists() else {}
    lines = ["# Change report: v1 -> v2 headline numbers", "",
             "v2 = consolidated release: excess-air range 5-25 % (v1: 5-20 %), base case = calibrated Plant A state (21.6 % excess air; v1 scenarios/optimisation used 20 %), Monte Carlo N = 2048 (v1: 4096). Everything else unchanged. Flag: relative change > 5 % (temperatures: > 5 K; indices, shares and fractions: > 0.05 absolute).", "",
             f"Generated {time.strftime('%Y-%m-%d %H:%M')} at commit {git_hash()[:8]}. Runtimes per v2 stage: " + ", ".join(f"{k} {v['wall_time_s']:.0f} s" for k, v in rt.items()), "",
             "| family | item | v1 | v2 | change | flag | explanation |", "|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['family']} | {r['item']} | {r['v1']} | {r['v2']} | {r['change']} | {r['flag']} | {r['explanation']} |")
    n_flag = int((df.flag != "").sum())
    lines += ["", f"{n_flag} of {len(df)} headline numbers changed beyond the threshold."]
    path.write_text("\n".join(lines) + "\n")
    return df


def git_hash() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


PRODUCERS = [
    ("data/literature_validation/digitized/", "manual WebPlotDigitizer export + validate_xf1989 tidy step"), ("data/literature_validation/xf1989_fit.yaml", "notebooks/03_validation_xf1989.ipynb (rdt.validate_xf1989)"),
    ("data/literature_validation/latham2008_plant_cases.csv", "scratch build script from thesis Appendix H (see README section)"), ("data/literature_validation/latham_fit.yaml", "notebooks/04-06 (rdt.latham_cases, rdt.latham_calibration, rdt.creep)"),
    ("data/kinetics/", "hand transcription of Xu & Froment 1989 Table 5-7"), ("data/creep_derived/", "rdt.creep_ingest_datasheet (Schmidt+Clemens Centralloy data sheets), rdt.creep (Yeh digitisation, NIMS ingestion template/config)"),
    ("data/design/parameter_ranges", "hand-written operating-space definition"), ("data/design/sobol_", "rdt.pipeline stage design (rdt.operate)"), ("data/design/surrogate_metrics_", "rdt.pipeline stage surrogate (rdt.surrogate)"),
    ("data/design/v1_vs_v2_diff.md", "rdt.report.write_diff_report"), ("data/design/operating_space.yaml", "hand-written"),
    ("data/lhs_runs/", "rdt.pipeline stage design (rdt.operate: LHS, closed-loop grid, Saltelli)"), ("data/scenarios/", "rdt.pipeline stage scenarios (rdt.scenarios)"),
    ("data/optimization/", "rdt.pipeline stage pareto (rdt.optimize, rdt.steam_credit)"), ("data/uq/", "rdt.pipeline stage uq (rdt.uq)"), ("data/pipeline_", "rdt.pipeline"),
]


def _producer(rel: str) -> str:
    for prefix, prod in PRODUCERS:
        if rel.startswith(prefix):
            return prod
    return "—"


def _commit(rel: str) -> str:
    try:
        out = subprocess.check_output(["git", "log", "-1", "--format=%h", "--", rel], cwd=ROOT, text=True).strip()
        return out or "uncommitted"
    except Exception:  # noqa: BLE001
        return "?"


def write_data_readme(path: Path = DATA / "README.md", big_mb: float = 5.0) -> pd.DataFrame:
    rows = []
    for p in sorted(DATA.rglob("*")):
        if p.is_dir() or p.name in (".gitkeep",) or p.suffix in (".log",) or p.name == "README.md":
            continue
        rel = str(p.relative_to(ROOT)); size = p.stat().st_size / 1e6
        rows.append({"file": rel, "size_MB": round(size, 3), "commit": _commit(rel), "producer": _producer(rel), "large": size > big_mb})
    df = pd.DataFrame(rows)
    zen = [r["file"] for r in rows if r["large"] or "/hourly_" in r["file"] or "/mc_v" in r["file"] or "sobol_saltelli" in r["file"] or "lhs_plantA" in r["file"]]
    lines = ["# Data inventory", "", f"Generated {time.strftime('%Y-%m-%d %H:%M')} by `rdt.report.write_data_readme` at commit {git_hash()[:8]}. Sizes in MB; commit = last commit touching the file (`uncommitted` = not yet in git).", "",
             "| file | size (MB) | commit | produced by |", "|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['file']} | {r['size_MB']:.3f}{' **(>5 MB)**' if r['large'] else ''} | {r['commit']} | {r['producer']} |")
    lines += ["", f"## Files larger than {big_mb:.0f} MB", ""] + [f"- {r['file']} ({r['size_MB']:.1f} MB)" for r in rows if r["large"]]
    lines += ["", "## Proposed Zenodo deposit at publication (list only; nothing deleted)", "",
              "Campaign outputs that are reproducible from `python -m rdt.pipeline --stage all` and are too large or too numerous for the repository:", ""]
    lines += [f"- {f}" for f in zen] + ["- models/ (trained surrogates, ~400 MB, git-ignored; regenerate with the surrogate stage)", "- literature/ (copyrighted PDFs: never deposited)"]
    path.write_text("\n".join(lines) + "\n")
    return df
