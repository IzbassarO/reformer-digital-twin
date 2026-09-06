"""RQ3: Pareto-optimal operating regimes of the Plant A twin (NSGA-II on GP surrogates, physics-verified).

Decision variables: load fraction, steam-to-carbon, inlet temperature, excess air, firing factor (inlet pressure
and catalyst activity fixed at base). Objectives: maximise H2 per tube, minimise fuel per kmol H2
(Q_comb / H2 rate), minimise life-consumption rate per kmol H2 (1/t_r / H2, Yeh placeholder curve, relative to the
base case). Constraints: dry CH4 slip <= base, T_wo,max <= 1193 K (920 degC, design-limit ASSUMPTION for HP-Nb
tubes), outlet T <= 1125 K. The final Pareto set is re-evaluated with the physics model and members that violate
a constraint under physics are dropped.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
from pymoo.optimize import minimize
from pymoo.termination import get_termination

from rdt import creep
from rdt import operate as op
from rdt import scenarios as sc

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "data" / "optimization"
DECISION = ("feed_per_tube_fraction", "steam_to_carbon", "inlet_T", "excess_air", "specific_firing_factor")
XL = np.array([0.60, 2.5, 800.0, 5.0, 0.85]); XU = np.array([1.10, 4.0, 900.0, 20.0, 1.15])
T_WO_LIMIT_K = 1193.0     # 920 degC: design-limit assumption for HP-Nb reformer tubes
T_OUT_LIMIT_K = 1125.0
OBJECTIVES = ("H2_net_kmol_h", "fuel_MJ_per_kmol_H2", "life_rate_per_kmol_H2_rel_base")
EXCESS_AIR_BASE = sc.EXCESS_AIR_SCEN     # 20 %, inside the decision bounds (Plant A derived value 21.6 % is not)


@dataclass
class BaseRef:
    inputs: Dict[str, float]
    H2: float
    fuel: float
    life_rate_per_H2: float
    slip: float
    T_wo_max: float
    T_out: float
    Q_comb_W: float
    row: Dict[str, float]


def _fuel_MJ_per_kmol(Q_comb_W: np.ndarray, H2: np.ndarray, n_tubes: int) -> np.ndarray:
    return Q_comb_W * 3600.0 / 1e6 / n_tubes / H2      # MJ per hour per tube / (kmol H2 per hour per tube)


def base_reference(tw: sc.TwinSurrogate) -> BaseRef:
    """Physics evaluation of the base point (Plant A inputs, excess air 20 %, firing 1.0)."""
    x = tw.base.inputs(); x["excess_air"] = EXCESS_AIR_BASE
    r = op.run_case(x, tw.base, tw.curve)
    fuel = float(_fuel_MJ_per_kmol(np.array([r["Q_comb_W"]]), np.array([r["H2_net_kmol_h"]]), tw.base.row.n_tubes)[0])
    return BaseRef(inputs=x, H2=r["H2_net_kmol_h"], fuel=fuel, life_rate_per_H2=r["life_consumption_rate_per_h"] / r["H2_net_kmol_h"],
                   slip=r["CH4_slip_dry_pct"], T_wo_max=r["T_wo_max_K"], T_out=r["T_out_K"], Q_comb_W=r["Q_comb_W"], row=r)


class ReformerProblem(Problem):
    def __init__(self, tw: sc.TwinSurrogate, ref: BaseRef, T_shift_K: float = 0.0):
        super().__init__(n_var=5, n_obj=3, n_ieq_constr=3, xl=XL.copy(), xu=XU.copy())
        self.tw, self.ref, self.T_shift = tw, ref, T_shift_K
        self.n_evals = 0

    def table(self, X: np.ndarray) -> pd.DataFrame:
        df = pd.DataFrame(X, columns=DECISION)
        df["inlet_P"] = self.ref.inputs["inlet_P"]; df["catalyst_activity"] = self.ref.inputs["catalyst_activity"]
        df["steam_to_carbon"] = np.round(df.steam_to_carbon, 4)   # caching of the equilibrium feature
        return df

    def evaluate_table(self, X: np.ndarray) -> pd.DataFrame:
        df = self.table(X)
        p = self.tw.predict(df, physics_fallback=False)
        Q = df.specific_firing_factor.to_numpy() * df.feed_per_tube_fraction.to_numpy() * self.tw.base.Q_comb_W
        H2 = p["H2_net_kmol_h"].to_numpy(); T_wo = p["T_wo_max_K"].to_numpy(); sig = p["sigma_hot_MPa"].to_numpy()
        t_r = self.tw.curve.time_to_rupture(T_wo + self.T_shift, sig)
        out = df.copy()
        out["H2_net_kmol_h"] = H2; out["Q_comb_W"] = Q; out["fuel_MJ_per_kmol_H2"] = _fuel_MJ_per_kmol(Q, H2, self.tw.base.row.n_tubes)
        out["T_wo_max_K"] = T_wo; out["sigma_hot_MPa"] = sig; out["log10_t_r_hot"] = np.log10(t_r)
        out["life_rate_per_kmol_H2_rel_base"] = (1.0 / t_r / H2) / self.ref.life_rate_per_H2
        out["CH4_slip_dry_pct"] = p["CH4_slip_dry_pct"].to_numpy(); out["T_out_K"] = p["T_out_K"].to_numpy()
        return out

    def _evaluate(self, X, out, *args, **kwargs):
        self.n_evals += len(X)
        t = self.evaluate_table(X)
        out["F"] = np.column_stack([-t.H2_net_kmol_h, t.fuel_MJ_per_kmol_H2, t.life_rate_per_kmol_H2_rel_base])
        out["G"] = np.column_stack([t.CH4_slip_dry_pct - self.ref.slip, t.T_wo_max_K - T_WO_LIMIT_K, t.T_out_K - T_OUT_LIMIT_K])


def run_nsga2(problem: ReformerProblem, pop: int = 200, gens: int = 300, seed: int = 0):
    algo = NSGA2(pop_size=pop, eliminate_duplicates=True)
    t0 = time.perf_counter()
    res = minimize(problem, algo, get_termination("n_gen", gens), seed=seed, verbose=False, save_history=False)
    return res, time.perf_counter() - t0


def physics_verify(df: pd.DataFrame, tw: sc.TwinSurrogate, ref: BaseRef, T_shift_K: float = 0.0) -> pd.DataFrame:
    """Re-evaluate rows with the physics model; add *_phys columns, discrepancies and constraint flags."""
    from joblib import Parallel, delayed
    X = df[list(op.INPUT_NAMES)].to_dict(orient="records")
    rows = Parallel(n_jobs=-1)(delayed(_phys)(x) for x in X)
    ph = pd.DataFrame(rows, index=df.index)
    out = df.copy()
    out["H2_net_kmol_h_phys"] = ph.H2_net_kmol_h; out["T_wo_max_K_phys"] = ph.T_wo_max_K; out["T_out_K_phys"] = ph.T_out_K
    out["CH4_slip_dry_pct_phys"] = ph.CH4_slip_dry_pct; out["sigma_hot_MPa_phys"] = ph.sigma_hot_MPa; out["converged_phys"] = ph.converged
    Q = out.Q_comb_W.to_numpy()
    out["fuel_MJ_per_kmol_H2_phys"] = _fuel_MJ_per_kmol(Q, ph.H2_net_kmol_h.to_numpy(), tw.base.row.n_tubes)
    t_r = tw.curve.time_to_rupture(ph.T_wo_max_K.to_numpy() + T_shift_K, ph.sigma_hot_MPa.to_numpy())
    out["life_rate_per_kmol_H2_rel_base_phys"] = (1.0 / t_r / ph.H2_net_kmol_h.to_numpy()) / ref.life_rate_per_H2
    out["feasible_phys"] = (ph.converged == True) & (ph.CH4_slip_dry_pct <= ref.slip + 1e-6) & (ph.T_wo_max_K <= T_WO_LIMIT_K + 1e-6) & (ph.T_out_K <= T_OUT_LIMIT_K + 1e-6)  # noqa: E712
    for o in OBJECTIVES:
        out[f"disc_{o}"] = out[o] - out[f"{o}_phys"]
    return out


def _phys(x):
    r = op.run_case(x, op.base_case(), creep.LarsonMillerCurve.from_yaml())
    keys = ("H2_net_kmol_h", "T_wo_max_K", "T_out_K", "CH4_slip_dry_pct", "sigma_hot_MPa", "converged")
    return {k: r.get(k, np.nan) for k in keys}


# ---------------------------------------------------------------------------
# Analyses
# ---------------------------------------------------------------------------
def knee_point(df: pd.DataFrame, cols=("H2_net_kmol_h", "fuel_MJ_per_kmol_H2", "life_rate_per_kmol_H2_rel_base"), signs=(-1, 1, 1)) -> int:
    """Index of the member closest to the utopia point after min-max normalisation (all objectives minimised)."""
    F = np.column_stack([s * df[c].to_numpy() for c, s in zip(cols, signs)])
    Fn = (F - F.min(0)) / np.where(F.max(0) - F.min(0) > 0, F.max(0) - F.min(0), 1.0)
    return int(df.index[np.argmin(np.linalg.norm(Fn, axis=1))])


def iso_production(df: pd.DataFrame, ref: BaseRef, band: float = 0.01) -> Dict[str, object]:
    m = df[np.abs(df.H2_net_kmol_h_phys / ref.H2 - 1.0) <= band]
    if m.empty:
        return {"n_in_band": 0}
    i_life = m.life_rate_per_kmol_H2_rel_base_phys.idxmin(); i_fuel = m.fuel_MJ_per_kmol_H2_phys.idxmin()
    return {"n_in_band": int(len(m)), "min_life_idx": int(i_life), "min_fuel_idx": int(i_fuel),
            "min_life": {"life_rel_base": float(m.loc[i_life, "life_rate_per_kmol_H2_rel_base_phys"]), "fuel_penalty_pct": float(100 * (m.loc[i_life, "fuel_MJ_per_kmol_H2_phys"] / ref.fuel - 1)),
                         "H2_rel_base_pct": float(100 * (m.loc[i_life, "H2_net_kmol_h_phys"] / ref.H2 - 1))},
            "min_fuel": {"fuel_rel_base_pct": float(100 * (m.loc[i_fuel, "fuel_MJ_per_kmol_H2_phys"] / ref.fuel - 1)), "life_rel_base": float(m.loc[i_fuel, "life_rate_per_kmol_H2_rel_base_phys"]),
                         "H2_rel_base_pct": float(100 * (m.loc[i_fuel, "H2_net_kmol_h_phys"] / ref.H2 - 1))}}


def base_domination(df: pd.DataFrame, ref: BaseRef) -> Dict[str, object]:
    """Pareto members that dominate the base in all three (physics-verified) objectives, and the margins."""
    dom = df[(df.H2_net_kmol_h_phys >= ref.H2) & (df.fuel_MJ_per_kmol_H2_phys <= ref.fuel) & (df.life_rate_per_kmol_H2_rel_base_phys <= 1.0)]
    out = {"n_dominating": int(len(dom)), "n_pareto": int(len(df)), "base_is_dominated": bool(len(dom) > 0)}
    if len(dom):
        out["max_H2_gain_pct"] = float(100 * (dom.H2_net_kmol_h_phys.max() / ref.H2 - 1))
        out["max_fuel_saving_pct"] = float(100 * (1 - dom.fuel_MJ_per_kmol_H2_phys.min() / ref.fuel))
        out["max_life_saving_pct"] = float(100 * (1 - dom.life_rate_per_kmol_H2_rel_base_phys.min()))
    return out


def regimes_table(df: pd.DataFrame, ref: BaseRef, iso: Dict[str, object]) -> pd.DataFrame:
    rows = {}
    base_row = {**ref.inputs, "H2_net_kmol_h_phys": ref.H2, "fuel_MJ_per_kmol_H2_phys": ref.fuel, "life_rate_per_kmol_H2_rel_base_phys": 1.0,
                "CH4_slip_dry_pct_phys": ref.slip, "T_wo_max_K_phys": ref.T_wo_max, "T_out_K_phys": ref.T_out, "Q_comb_W": ref.Q_comb_W}
    rows["base"] = base_row
    if iso.get("n_in_band"):
        rows["min_life_iso_H2"] = df.loc[iso["min_life_idx"]].to_dict(); rows["min_fuel_iso_H2"] = df.loc[iso["min_fuel_idx"]].to_dict()
    rows["min_life_overall"] = df.loc[df.life_rate_per_kmol_H2_rel_base_phys.idxmin()].to_dict()
    rows["min_fuel_overall"] = df.loc[df.fuel_MJ_per_kmol_H2_phys.idxmin()].to_dict()
    rows["knee"] = df.loc[knee_point(df, cols=("H2_net_kmol_h_phys", "fuel_MJ_per_kmol_H2_phys", "life_rate_per_kmol_H2_rel_base_phys"))].to_dict()
    rows["max_H2"] = df.loc[df.H2_net_kmol_h_phys.idxmax()].to_dict()
    cols = list(op.INPUT_NAMES) + ["Q_comb_W", "H2_net_kmol_h_phys", "fuel_MJ_per_kmol_H2_phys", "life_rate_per_kmol_H2_rel_base_phys", "CH4_slip_dry_pct_phys", "T_wo_max_K_phys", "T_out_K_phys"]
    t = pd.DataFrame(rows).T[cols]
    t["H2_rel_base_pct"] = 100 * (t.H2_net_kmol_h_phys / ref.H2 - 1); t["fuel_rel_base_pct"] = 100 * (t.fuel_MJ_per_kmol_H2_phys / ref.fuel - 1)
    return t.astype(float)


def run_all(pop: int = 200, gens: int = 300, seed: int = 0, save: bool = True) -> Dict[str, object]:
    tw = sc.TwinSurrogate(); ref = base_reference(tw)
    prob = ReformerProblem(tw, ref)
    res, dt = run_nsga2(prob, pop, gens, seed)
    pareto = prob.evaluate_table(res.X); pareto["rank"] = 0
    ver = physics_verify(pareto, tw, ref)
    disc = {o: float(np.nanmax(np.abs(ver[f"disc_{o}"]))) for o in OBJECTIVES}
    feas = ver[ver.feasible_phys].copy()
    iso = iso_production(feas, ref); dom = base_domination(feas, ref); reg = regimes_table(feas, ref, iso)
    # robustness: creep-curve shift +/-20 K
    rob = {}
    for shift in (-20.0, 20.0):
        p2 = ReformerProblem(tw, ref, T_shift_K=shift); r2, dt2 = run_nsga2(p2, pop, gens, seed)
        X2 = r2.X; X1 = res.X
        # normalised nearest-neighbour distance between decision sets (both directions)
        n1 = (X1 - XL) / (XU - XL); n2 = (X2 - XL) / (XU - XL)
        d12 = np.array([np.min(np.linalg.norm(n2 - a, axis=1)) for a in n1]); d21 = np.array([np.min(np.linalg.norm(n1 - b, axis=1)) for b in n2])
        rob[f"shift_{shift:+.0f}K"] = {"n_pareto": int(len(X2)), "mean_nn_distance_norm": float(0.5 * (d12.mean() + d21.mean())), "max_nn_distance_norm": float(max(d12.max(), d21.max())),
                                       "decision_ranges": {k: [float(X2[:, i].min()), float(X2[:, i].max())] for i, k in enumerate(DECISION)}, "wall_time_s": dt2}
    out = {"pareto_verified": ver, "feasible": feas, "regimes": reg, "iso_production": iso, "base_domination": dom, "ref": ref,
           "discrepancy_max": disc, "n_pareto": int(len(pareto)), "n_feasible_phys": int(len(feas)), "wall_time_s": dt, "n_evals": prob.n_evals,
           "robustness": rob, "settings": {"pop": pop, "gens": gens, "seed": seed, "T_wo_limit_K": T_WO_LIMIT_K, "T_out_limit_K": T_OUT_LIMIT_K, "slip_limit_pct": ref.slip, "excess_air_base": EXCESS_AIR_BASE}}
    if save:
        OUT_DIR.mkdir(exist_ok=True)
        ver.to_csv(PARETO_CSV, index=False); reg.to_csv(REGIMES_CSV)
        meta = {k: v for k, v in out.items() if k not in ("pareto_verified", "feasible", "regimes", "ref")}
        meta["base_reference"] = {"inputs": ref.inputs, "H2": ref.H2, "fuel_MJ_per_kmol_H2": ref.fuel, "slip": ref.slip, "T_wo_max": ref.T_wo_max, "T_out": ref.T_out, "Q_comb_W": ref.Q_comb_W}
        meta["decision_ranges_pareto"] = {k: [float(pareto[k].min()), float(pareto[k].max())] for k in DECISION}
        meta["tag"] = TAG; meta["decision_bounds"] = {k: [float(XL[i]), float(XU[i])] for i, k in enumerate(DECISION)}
        json.dump(meta, open(PARETO_META, "w"), indent=2, default=float)
    return out


# ---------------------------------------------------------------------------
# Variant: fuel objective including steam raising and feed preheat (documented sensitivity)
# ---------------------------------------------------------------------------
DH_VAP_298_MJ_PER_KMOL = 44.0   # enthalpy of vaporisation of water at 25 degC


def process_heat_per_tube_MJ_h(df: pd.DataFrame, tw: sc.TwinSurrogate) -> np.ndarray:
    """Steam-raising (liquid water at 298 K -> vapour at T_in) plus dry-feed preheat (298 K -> T_in) per tube [MJ/h].

    Gas-phase enthalpies from gri30; C4-C6 counted as C3H8. Pressure effect on steam enthalpy neglected.
    Not part of the specified fuel objective; used only for the sensitivity variant."""
    import cantera as ct
    gas = ct.Solution("gri30.yaml")
    b = tw.base
    dry = {("C3H8" if s in ("C4H10", "C5H12", "C6H14") else s): v for s, v in b.dry_feed_kmol_h.items() if v > 0}
    agg: Dict[str, float] = {}
    for s, v in dry.items():
        agg[s] = agg.get(s, 0.0) + v
    out = np.empty(len(df))
    cache: Dict[tuple, float] = {}
    for i, (sc_, load, T) in enumerate(zip(df.steam_to_carbon, df.feed_per_tube_fraction, df.inlet_T)):
        key = (round(float(sc_), 4), round(float(T), 2))
        if key not in cache:
            gas.TPX = 298.15, ct.one_atm, agg; h0 = gas.enthalpy_mole; gas.TP = float(T), ct.one_atm; h1 = gas.enthalpy_mole
            q_dry = sum(agg.values()) * (h1 - h0) / 1e6                       # MJ/h per unit load
            gas.TPX = 298.15, ct.one_atm, {"H2O": 1.0}; hw0 = gas.enthalpy_mole; gas.TP = float(T), ct.one_atm; hw1 = gas.enthalpy_mole
            n_steam = float(sc_) * b.carbon_kmol_h
            q_steam = n_steam * ((hw1 - hw0) / 1e6 + DH_VAP_298_MJ_PER_KMOL)
            cache[key] = q_dry + q_steam
        out[i] = cache[key] * float(load)
    return out


class ReformerProblemTotalHeat(ReformerProblem):
    """Same problem but the fuel objective is (Q_comb + steam raising + feed preheat) per kmol H2."""

    def evaluate_table(self, X: np.ndarray) -> pd.DataFrame:
        t = super().evaluate_table(X)
        t["process_heat_MJ_h"] = process_heat_per_tube_MJ_h(t, self.tw)
        t["fuel_MJ_per_kmol_H2"] = (t.Q_comb_W * 3600.0 / 1e6 / self.tw.base.row.n_tubes + t.process_heat_MJ_h) / t.H2_net_kmol_h
        return t


def run_variant_total_heat(pop: int = 200, gens: int = 300, seed: int = 0, save: bool = True) -> Dict[str, object]:
    tw = sc.TwinSurrogate(); ref = base_reference(tw)
    prob = ReformerProblemTotalHeat(tw, ref)
    base_t = prob.evaluate_table(np.array([[ref.inputs[k] for k in DECISION]]))
    ref_total_fuel = float(base_t.fuel_MJ_per_kmol_H2.iloc[0])
    res, dt = run_nsga2(prob, pop, gens, seed)
    pareto = prob.evaluate_table(res.X)
    ver = physics_verify(pareto, tw, ref)
    ver["fuel_MJ_per_kmol_H2_phys"] = (ver.Q_comb_W * 3600.0 / 1e6 / tw.base.row.n_tubes + ver.process_heat_MJ_h) / ver.H2_net_kmol_h_phys
    feas = ver[ver.feasible_phys].copy()
    m = feas[np.abs(feas.H2_net_kmol_h_phys / ref.H2 - 1.0) <= 0.01]
    iso = {}
    if len(m):
        i = m.life_rate_per_kmol_H2_rel_base_phys.idxmin(); j = m.fuel_MJ_per_kmol_H2_phys.idxmin()
        iso = {"min_life": {"life_rel_base": float(m.loc[i, "life_rate_per_kmol_H2_rel_base_phys"]), "total_fuel_rel_base_pct": float(100 * (m.loc[i, "fuel_MJ_per_kmol_H2_phys"] / ref_total_fuel - 1)),
                            "steam_to_carbon": float(m.loc[i, "steam_to_carbon"]), "firing": float(m.loc[i, "specific_firing_factor"]), "inlet_T": float(m.loc[i, "inlet_T"]), "excess_air": float(m.loc[i, "excess_air"]), "load": float(m.loc[i, "feed_per_tube_fraction"])},
               "min_fuel": {"total_fuel_rel_base_pct": float(100 * (m.loc[j, "fuel_MJ_per_kmol_H2_phys"] / ref_total_fuel - 1)), "life_rel_base": float(m.loc[j, "life_rate_per_kmol_H2_rel_base_phys"]),
                            "steam_to_carbon": float(m.loc[j, "steam_to_carbon"]), "firing": float(m.loc[j, "specific_firing_factor"]), "inlet_T": float(m.loc[j, "inlet_T"]), "excess_air": float(m.loc[j, "excess_air"]), "load": float(m.loc[j, "feed_per_tube_fraction"])}}
    out = {"pareto_verified": ver, "feasible": feas, "ref_total_fuel_MJ_per_kmol_H2": ref_total_fuel, "iso_production": iso, "wall_time_s": dt,
           "decision_ranges": {k: [float(feas[k].min()), float(feas[k].max())] for k in DECISION}, "n_feasible": int(len(feas))}
    if save:
        ver.to_csv(PARETO_ALT_CSV, index=False)
        json.dump({k: v for k, v in out.items() if k not in ("pareto_verified", "feasible")}, open(PARETO_ALT_META, "w"), indent=2, default=float)
    return out


PARETO_CSV = OUT_DIR / "pareto_v1.csv"; REGIMES_CSV = OUT_DIR / "regimes_v1.csv"; PARETO_META = OUT_DIR / "pareto_v1_meta.json"
PARETO_ALT_CSV = OUT_DIR / "pareto_v1_alt_totalheat.csv"; PARETO_ALT_META = OUT_DIR / "pareto_v1_alt_totalheat_meta.json"; TAG = "v1"


def use_config(cfg) -> None:
    """Point the module at a rdt.config.RunConfig: base excess air, decision bounds from the ranges file, output names."""
    global EXCESS_AIR_BASE, XL, XU, PARETO_CSV, REGIMES_CSV, PARETO_META, PARETO_ALT_CSV, PARETO_ALT_META, TAG
    EXCESS_AIR_BASE = float(cfg.excess_air_base)
    r = pd.read_csv(cfg.ranges_csv).set_index("parameter")
    XL = np.array([r.loc[k, "min"] for k in DECISION], float); XU = np.array([r.loc[k, "max"] for k in DECISION], float)
    PARETO_CSV, REGIMES_CSV, PARETO_META = Path(cfg.pareto_csv), Path(cfg.regimes_csv), Path(cfg.pareto_meta)
    PARETO_ALT_CSV, PARETO_ALT_META, TAG = Path(cfg.pareto_alt_csv), Path(cfg.pareto_alt_meta), cfg.tag
