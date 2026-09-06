"""Steam-credit parameter for the RQ3 fuel objective (RQ4 item 5).

``fuel_total = firebox fuel + (1 - alpha_credit) * (steam-raising + feed-preheat heat)`` per kmol H2, with
``alpha_credit`` in {0, 0.25, 0.5, 0.75, 1}: alpha = 1 is the firebox-only objective of RQ3, alpha = 0 charges all
process heat. For each alpha a two-objective NSGA-II (min life-consumption per kmol H2, min fuel_total per kmol H2)
is run on the GP surrogates under the iso-production constraint |H2 - H2_base| <= 1 % and the RQ3 constraints, and
the front is verified with the physics model.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
from pymoo.optimize import minimize
from pymoo.termination import get_termination

from rdt import optimize as ox
from rdt import scenarios as sc

OUT_DIR = ox.OUT_DIR
ALPHAS = (0.0, 0.25, 0.5, 0.75, 1.0)


class IsoProductionProblem(Problem):
    def __init__(self, tw: sc.TwinSurrogate, ref: ox.BaseRef, alpha: float, band: float = 0.01):
        super().__init__(n_var=5, n_obj=2, n_ieq_constr=5, xl=ox.XL.copy(), xu=ox.XU.copy())
        self.inner = ox.ReformerProblemTotalHeat(tw, ref)
        self.tw, self.ref, self.alpha, self.band = tw, ref, alpha, band

    def evaluate_table(self, X: np.ndarray) -> pd.DataFrame:
        t = self.inner.evaluate_table(X)   # has process_heat_MJ_h and Q_comb_W
        firebox = t.Q_comb_W * 3600.0 / 1e6 / self.tw.base.row.n_tubes
        t["fuel_total_MJ_per_kmol_H2"] = (firebox + (1.0 - self.alpha) * t.process_heat_MJ_h) / t.H2_net_kmol_h
        return t

    def _evaluate(self, X, out, *args, **kwargs):
        t = self.evaluate_table(X)
        out["F"] = np.column_stack([t.life_rate_per_kmol_H2_rel_base, t.fuel_total_MJ_per_kmol_H2])
        dev = np.abs(t.H2_net_kmol_h / self.ref.H2 - 1.0) - self.band
        out["G"] = np.column_stack([t.CH4_slip_dry_pct - self.ref.slip, t.T_wo_max_K - ox.T_WO_LIMIT_K, t.T_out_K - ox.T_OUT_LIMIT_K, dev, dev])


def run_steam_credit(alphas=ALPHAS, pop: int = 120, gens: int = 150, seed: int = 0, save: bool = True) -> Dict[str, object]:
    tw = sc.TwinSurrogate(); ref = ox.base_reference(tw)
    base_tab = ox.ReformerProblemTotalHeat(tw, ref).evaluate_table(np.array([[ref.inputs[k] for k in ox.DECISION]]))
    base_ph = float(base_tab.process_heat_MJ_h.iloc[0]); base_firebox = ref.fuel
    results = {"alphas": list(alphas), "base": {"firebox_MJ_per_kmol_H2": base_firebox, "process_heat_MJ_per_kmol_H2": base_ph / ref.H2}, "per_alpha": {}}
    fronts = []
    for a in alphas:
        prob = IsoProductionProblem(tw, ref, a)
        t0 = time.perf_counter()
        res = minimize(prob, NSGA2(pop_size=pop, eliminate_duplicates=True), get_termination("n_gen", gens), seed=seed, verbose=False)
        dt = time.perf_counter() - t0
        if res.X is None:
            results["per_alpha"][str(a)] = {"n_front": 0, "note": "no feasible solution"}; continue
        front = prob.evaluate_table(np.atleast_2d(res.X))
        ver = ox.physics_verify(front, tw, ref)
        firebox_ph = ver.Q_comb_W * 3600.0 / 1e6 / tw.base.row.n_tubes
        ver["fuel_total_MJ_per_kmol_H2_phys"] = (firebox_ph + (1.0 - a) * ver.process_heat_MJ_h) / ver.H2_net_kmol_h_phys
        ver["alpha_credit"] = a
        feas = ver[ver.feasible_phys & (np.abs(ver.H2_net_kmol_h_phys / ref.H2 - 1.0) <= 0.015)]
        fronts.append(ver)
        base_total = base_firebox + (1.0 - a) * base_ph / ref.H2
        rec = {"n_front": int(len(front)), "n_feasible_phys": int(len(feas)), "wall_time_s": dt, "base_fuel_total_MJ_per_kmol_H2": base_total}
        if len(feas):
            i = feas.life_rate_per_kmol_H2_rel_base_phys.idxmin(); j = feas.fuel_total_MJ_per_kmol_H2_phys.idxmin()
            k = ox.knee_point(feas, cols=("life_rate_per_kmol_H2_rel_base_phys", "fuel_total_MJ_per_kmol_H2_phys"), signs=(1, 1))
            for lab, idx in (("min_life", i), ("min_fuel", j), ("knee", k)):
                r = feas.loc[idx]
                rec[lab] = {"steam_to_carbon": float(r.steam_to_carbon), "load": float(r.feed_per_tube_fraction), "firing": float(r.specific_firing_factor), "inlet_T": float(r.inlet_T),
                            "excess_air": float(r.excess_air), "life_rel_base": float(r.life_rate_per_kmol_H2_rel_base_phys),
                            "fuel_total_rel_base_pct": float(100 * (r.fuel_total_MJ_per_kmol_H2_phys / base_total - 1)), "T_wo_max_K": float(r.T_wo_max_K_phys), "H2_rel_base_pct": float(100 * (r.H2_net_kmol_h_phys / ref.H2 - 1))}
            rec["SC_range_front"] = [float(feas.steam_to_carbon.min()), float(feas.steam_to_carbon.max())]
        results["per_alpha"][str(a)] = rec
    if save:
        pd.concat(fronts).to_csv(FRONTS_CSV, index=False)
        json.dump(results, open(RESULT_JSON, "w"), indent=2, default=float)
    return results


FRONTS_CSV = OUT_DIR / "steam_credit_fronts_v1.csv"; RESULT_JSON = OUT_DIR / "steam_credit_v1.json"


def use_config(cfg) -> None:
    global FRONTS_CSV, RESULT_JSON
    FRONTS_CSV = Path(cfg.steam_credit_fronts); RESULT_JSON = Path(cfg.steam_credit_json)
