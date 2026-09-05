"""Latham (2008) plant cases: build tube, bed, feed and flue gas from the plant-case CSV and run the coupled model.

Data: ``data/literature_validation/latham2008_plant_cases.csv`` (thesis Appendix H, Tables 29-31).
Baseline parameters: Latham et al. (2011) Table 6 best-fit values (L_q 0.48, alpha_top 0.182, f_htg 1.68,
voidage 0.607, f_prx 0.20, f_ctube 0.38), eta 0.1 with 0.05 in the top section, plus geometry from the
thesis (OD 0.146 m, heated length 12.5 m, lambda_tube 29.6 W/(m K), d_p 5.40 mm, rho_bed 1100 kg/m3) and an
ASSUMED wall thickness of 15 mm (ID 0.116 m) and initial exchange factor F_gt = 0.35.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd

from rdt import furnace as fu
from rdt import reactor1d as r1

ROOT = Path(__file__).resolve().parents[2]
CASES_CSV = ROOT / "data" / "literature_validation" / "latham2008_plant_cases.csv"
FIT_YAML = ROOT / "data" / "literature_validation" / "latham_fit.yaml"


@dataclass(frozen=True)
class LathamParams:
    L_q: float = 0.48
    alpha_top: float = 0.182
    f_htg: float = 1.68
    voidage: float = 0.607
    activity: float = 0.20          # Latham's f_prx
    eta: float = 0.1
    eta_top: float = 0.05
    f_ctube: float = 0.38
    F_gt: float = 0.35              # initial guess, to be fitted
    f_loss: float = 0.02
    tube_od_m: float = 0.146
    wall_thickness_m: float = 0.015  # ASSUMED (not given by Latham)
    heated_length_m: float = 12.5
    lambda_tube: float = 29.6
    d_p_m: float = 5.40e-3
    rho_bed: float = 1100.0
    heat_transfer: str = "leva_grummer"
    n_sections: int = 15

    @property
    def tube_id_m(self) -> float:
        return self.tube_od_m - 2.0 * self.wall_thickness_m


BASELINE = LathamParams()


def load_cases(path: Path = CASES_CSV) -> pd.DataFrame:
    return pd.read_csv(path)


def build_case(row: pd.Series, p: LathamParams = BASELINE):
    """Return (tube, bed, feed, flue, geom) for one CSV row."""
    tube = r1.TubeGeometry(d_i=p.tube_id_m, d_o=p.tube_od_m, L_heated=p.heated_length_m, lambda_tube=p.lambda_tube)
    z_top = p.heated_length_m / p.n_sections

    def eta(z):
        e = p.eta_top if z < z_top else p.eta
        return (e, e, e)

    bed = r1.CatalystBed(rho_bed=p.rho_bed, voidage=p.voidage, d_p=p.d_p_m, eta=eta, activity=p.activity)
    n = float(row.feed_per_tube_kmol_h)
    F = {"CH4": n * row.feed_x_C1_molpct / 100, "C2H6": n * row.feed_x_C2_molpct / 100,
         "C3H8": n * row.feed_x_C3_molpct / 100, "C4H10": n * row.feed_x_C4_molpct / 100,
         "C5H12": n * row.feed_x_C5_molpct / 100, "C6H14": n * row.feed_x_C6plus_molpct / 100,
         "H2": n * row.feed_x_H2_molpct / 100, "CO": n * row.feed_x_CO_molpct / 100,
         "N2": n * row.feed_x_N2_molpct / 100, "CO2": n * row.feed_x_CO2_molpct / 100,
         "H2O": n * row.feed_x_H2O_molpct / 100}
    feed = r1.Feed(T_in=float(row.T_in_K), P_in=float(row.P_in_Pa) / 1e5, F=F, inlet_higher_alkanes="latham")
    fuel = {c.replace("fuel_x_", ""): float(row[c]) for c in row.index if c.startswith("fuel_x_")}
    off = {c.replace("offgas_x_", ""): float(row[c]) for c in row.index if c.startswith("offgas_x_")}
    air = {"N2": float(row.air_x_N2), "CO2": float(row.air_x_CO2), "H2O": float(row.air_x_H2O),
           "O2": float(row.air_x_O2_by_difference)}
    flue = fu.FlueGas.from_combustion([(float(row.fuel_gas_kmol_h), fuel, float(row.fuel_gas_T_K)),
                                       (float(row.offgas_kmol_h), off, float(row.offgas_T_K)),
                                       (float(row.air_kmol_h), air, float(row.air_T_K))])
    geom = fu.FurnaceGeometry(N_tubes=int(row.n_tubes), L=p.heated_length_m,
                              A_flue_free=fu.LATHAM_FURNACE.A_flue_free, D_h=fu.LATHAM_FURNACE.D_h,
                              description=fu.LATHAM_FURNACE.description)
    return tube, bed, feed, flue, geom


def run_case(row: pd.Series, p: LathamParams = BASELINE, n_out: int = 201) -> fu.CoupledResult:
    tube, bed, feed, flue, geom = build_case(row, p)
    return fu.simulate_coupled(tube, bed, feed, geom, flue, fu.HeatRelease(p.L_q, p.alpha_top, p.f_loss, p.n_sections),
                               fu.FurnaceParams(F_gt=p.F_gt, f_ctube=p.f_ctube), f_htg=p.f_htg,
                               heat_transfer=p.heat_transfer, n_out=n_out)


def compare(row: pd.Series, res: fu.CoupledResult) -> Dict[str, Dict[str, float]]:
    """Predicted vs measured outlet quantities and tube-wall temperatures for one case."""
    tr = res.tube
    zu, zl = float(row.TWT_upper_z_m), float(row.TWT_lower_z_m)
    pred = {"out_T_K": float(tr.T[-1]), "out_P_Pa": float(tr.P[-1] * 1e5),
            "out_x_H2_wet_molpct": float(100 * tr.X["H2"][-1]), "out_x_CH4_wet_molpct": float(100 * tr.X["CH4"][-1]),
            "out_x_CO_wet_molpct": float(100 * tr.X["CO"][-1]), "out_x_CO2_wet_molpct": float(100 * tr.X["CO2"][-1]),
            "flue_out_T_K": float(res.T_fg[-1]),
            "TWT_upper_K": float(np.interp(zu, res.z, res.T_wo)), "TWT_lower_K": float(np.interp(zl, res.z, res.T_wo))}
    meas = {k: float(row[k]) for k in pred}
    out = {k: {"predicted": pred[k], "measured": meas[k], "difference": pred[k] - meas[k]} for k in pred}
    zf, tpk = res.twt_peak()
    out["TWT_peak"] = {"z_frac": zf, "z_m": zf * res.z[-1], "T_wo_K": tpk}
    out["T_fg_max_K"] = {"value": float(res.T_fg.max()), "z_m": float(res.z[int(np.argmax(res.T_fg))])}
    out["conversion_CH4"] = {"predicted": float(tr.conversion_CH4[-1])}
    out["energy"] = dict(res.energy)
    out["runtime_s"] = res.wall_time_s
    return out


def run_all(p: LathamParams = BASELINE, cases: Optional[pd.DataFrame] = None) -> Dict[str, Dict[str, object]]:
    df = load_cases() if cases is None else cases
    out = {}
    for _, row in df.iterrows():
        res = run_case(row, p)
        out[str(row.case)] = {"comparison": compare(row, res), "result": res, "plant_rate_pct": int(row.plant_rate_pct)}
    return out


def summary_table(runs: Dict[str, Dict[str, object]]) -> pd.DataFrame:
    rows = []
    for case, d in runs.items():
        c = d["comparison"]
        r = {"case": case, "rate_pct": d["plant_rate_pct"]}
        for k in ("out_T_K", "out_P_Pa", "out_x_H2_wet_molpct", "out_x_CH4_wet_molpct", "flue_out_T_K", "TWT_upper_K", "TWT_lower_K"):
            r[f"{k}_pred"] = c[k]["predicted"]; r[f"{k}_meas"] = c[k]["measured"]; r[f"{k}_diff"] = c[k]["difference"]
        r["TWT_peak_frac"] = c["TWT_peak"]["z_frac"]; r["TWT_peak_K"] = c["TWT_peak"]["T_wo_K"]; r["T_fg_max_K"] = c["T_fg_max_K"]["value"]
        r["energy_closure"] = c["energy"]["closure_rel_error"]
        rows.append(r)
    return pd.DataFrame(rows).set_index("case")
