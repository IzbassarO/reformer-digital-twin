"""Operating-space evaluation of the calibrated digital twin (Plant A geometry and feed composition).

Inputs (see ``data/design/parameter_ranges.csv`` / ``operating_space.yaml``):

* ``steam_to_carbon`` [-]: steam changed only, hydrocarbon/inert dry feed fixed in composition;
* ``feed_per_tube_fraction`` [-]: dry (steam-free) feed per tube relative to Plant A;
* ``specific_firing_factor`` [-]: Q_comb per unit dry feed relative to Plant A, i.e.
  ``Q_comb = factor * feed_fraction * Q_comb,A`` (fuel and PSA off-gas streams scaled together);
* ``inlet_T`` [K], ``inlet_P`` [bar], ``catalyst_activity`` [-] (Latham's f_prx);
* ``excess_air`` [%]: combustion air set so that O2 supplied = (1 + EA/100) x stoichiometric; the flue-gas
  flow and composition follow from complete combustion.

Fixed: wall 15 mm, calibrated furnace/tube parameters (``latham_fit.yaml``, ``calibration_alpha_top_fixed``).
Life metrics use the Larson-Miller master curve selected by the run configuration (``cfg.creep_alloy``,
see :func:`rdt.creep.available_alloys`) and are reported as the consumption rate ``1/t_r`` relative to the
Plant A base case.
"""

from __future__ import annotations

import json
import subprocess
import time
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence

import numpy as np
import pandas as pd
import yaml
from scipy.optimize import brentq
from scipy.stats import qmc

from rdt import creep
from rdt import furnace as fu
from rdt import latham_cases as lc
from rdt import reactor1d as r1

ROOT = Path(__file__).resolve().parents[2]
RANGES_CSV = ROOT / "data" / "design" / "parameter_ranges.csv"
SPACE_YAML = ROOT / "data" / "design" / "operating_space.yaml"
LHS_DIR = ROOT / "data" / "lhs_runs"
INPUT_NAMES = ("steam_to_carbon", "feed_per_tube_fraction", "specific_firing_factor", "inlet_T", "inlet_P",
               "catalyst_activity", "excess_air")
BASE_CASE = "Plant_A"
TARGET_T_OUT_K = 1105.55   # Plant A measured outlet process-gas temperature (thesis Table 31)


# ---------------------------------------------------------------------------
# Ranges and base case
# ---------------------------------------------------------------------------
def load_ranges(path: Optional[Path] = None) -> pd.DataFrame:
    df = pd.read_csv(path or RANGES_CSV).set_index("parameter")
    missing = set(INPUT_NAMES) - set(df.index)
    if missing:
        raise ValueError(f"parameter_ranges.csv lacks {sorted(missing)}")
    return df.loc[list(INPUT_NAMES)]


def bounds_arrays(ranges: Optional[pd.DataFrame] = None):
    r = load_ranges() if ranges is None else ranges
    return r["min"].to_numpy(float), r["max"].to_numpy(float)


@dataclass(frozen=True)
class BaseCase:
    """Plant A calibrated state used as the reference (thesis Appendix H, Table 29-31)."""

    params: lc.LathamParams
    row: pd.Series
    dry_feed_kmol_h: Dict[str, float]      # per tube, steam excluded
    carbon_kmol_h: float                   # carbon atoms per tube
    steam_to_carbon: float
    T_in_K: float
    P_in_bar: float
    Q_comb_W: float
    excess_air_pct: float
    fuel_streams: tuple                    # ((n, X, T), (n, X, T)) fuel gas, off-gas
    air_X: Dict[str, float]
    air_T_K: float
    o2_stoich_kmol_h: float

    def inputs(self) -> Dict[str, float]:
        return {"steam_to_carbon": self.steam_to_carbon, "feed_per_tube_fraction": 1.0, "specific_firing_factor": 1.0,
                "inlet_T": self.T_in_K, "inlet_P": self.P_in_bar, "catalyst_activity": self.params.activity,
                "excess_air": self.excess_air_pct}


def calibrated_params(key: str = "calibration_alpha_top_fixed") -> lc.LathamParams:
    doc = yaml.safe_load(lc.FIT_YAML.read_text())
    bf = doc[key]["best_fit"]
    if "full_parameters" in bf:
        return lc.LathamParams(**bf["full_parameters"])
    from rdt import latham_calibration as cal
    return replace(cal.BASE, alpha_top=0.182, **bf["params"])


def base_case(case: str = BASE_CASE) -> BaseCase:
    df = lc.load_cases(); row = df[df.case == case].iloc[0]
    p = calibrated_params()
    n = float(row.feed_per_tube_kmol_h)
    dry = {"CH4": n * row.feed_x_C1_molpct / 100, "C2H6": n * row.feed_x_C2_molpct / 100, "C3H8": n * row.feed_x_C3_molpct / 100,
           "C4H10": n * row.feed_x_C4_molpct / 100, "C5H12": n * row.feed_x_C5_molpct / 100, "C6H14": n * row.feed_x_C6plus_molpct / 100,
           "H2": n * row.feed_x_H2_molpct / 100, "CO": n * row.feed_x_CO_molpct / 100, "N2": n * row.feed_x_N2_molpct / 100,
           "CO2": n * row.feed_x_CO2_molpct / 100}
    carbon = sum(k * dry[s] for s, k in {"CH4": 1, "C2H6": 2, "C3H8": 3, "C4H10": 4, "C5H12": 5, "C6H14": 6}.items())
    steam = n * row.feed_x_H2O_molpct / 100
    fuel = {c.replace("fuel_x_", ""): float(row[c]) for c in row.index if c.startswith("fuel_x_")}
    off = {c.replace("offgas_x_", ""): float(row[c]) for c in row.index if c.startswith("offgas_x_")}
    air = {"N2": float(row.air_x_N2), "CO2": float(row.air_x_CO2), "H2O": float(row.air_x_H2O), "O2": float(row.air_x_O2_by_difference)}
    streams = ((float(row.fuel_gas_kmol_h), fuel, float(row.fuel_gas_T_K)), (float(row.offgas_kmol_h), off, float(row.offgas_T_K)))
    o2_st = sum(nn * x * fu.O2_STOICH.get(sp, 0.0) for nn, X, _ in streams for sp, x in X.items())
    flue = fu.FlueGas.from_combustion([*streams, (float(row.air_kmol_h), air, float(row.air_T_K))])
    return BaseCase(params=p, row=row, dry_feed_kmol_h=dry, carbon_kmol_h=carbon, steam_to_carbon=steam / carbon,
                    T_in_K=float(row.T_in_K), P_in_bar=float(row.P_in_Pa) / 1e5, Q_comb_W=flue.Q_comb_W,
                    excess_air_pct=flue.excess_air_pct, fuel_streams=streams, air_X=air, air_T_K=float(row.air_T_K),
                    o2_stoich_kmol_h=o2_st)


# ---------------------------------------------------------------------------
# One run
# ---------------------------------------------------------------------------
def build_inputs(base: BaseCase, x: Mapping[str, float]):
    """Feed, flue gas and parameters for an input dict ``x`` (keys :data:`INPUT_NAMES`)."""
    load = float(x["feed_per_tube_fraction"])
    F = {s: v * load for s, v in base.dry_feed_kmol_h.items()}
    F["H2O"] = float(x["steam_to_carbon"]) * base.carbon_kmol_h * load
    feed = r1.Feed(T_in=float(x["inlet_T"]), P_in=float(x["inlet_P"]), F=F, inlet_higher_alkanes="latham")
    scale = float(x["specific_firing_factor"]) * load
    streams = [(n * scale, X, T) for n, X, T in base.fuel_streams]
    n_air = (1.0 + float(x["excess_air"]) / 100.0) * base.o2_stoich_kmol_h * scale / base.air_X["O2"]
    flue = fu.FlueGas.from_combustion([*streams, (n_air, base.air_X, base.air_T_K)])
    p = replace(base.params, activity=float(x["catalyst_activity"]))
    return feed, flue, p


def run_case(x: Mapping[str, float], base: Optional[BaseCase] = None,
             curve: Optional[creep.LarsonMillerCurve] = None, base_rate: Optional[float] = None) -> Dict[str, float]:
    """Run the coupled model for one input point and return a flat dict of inputs and outputs."""
    base = base or base_case()
    curve = curve or creep.active_curve()
    out = {k: float(x[k]) for k in INPUT_NAMES}
    t0 = time.perf_counter()
    try:
        feed, flue, p = build_inputs(base, x)
        tube, bed, _, _, geom = lc.build_case(base.row, p)
        res = fu.simulate_coupled(tube, bed, feed, geom, flue, fu.HeatRelease(p.L_q, p.alpha_top, p.f_loss, p.n_sections, p.inlet_mode),
                                  fu.FurnaceParams(F_gt=p.F_gt, f_ctube=p.f_ctube), f_htg=p.f_htg, heat_transfer=p.heat_transfer)
        if not res.success:
            raise RuntimeError(res.message)
        tr = res.tube; dry = tr.dry_mole_fractions
        k = int(np.argmax(res.T_wo)); T_mid = 0.5 * (res.T_wo + res.T_wi)
        sigma = creep.hoop_stress(float(tr.P[k]), p.tube_od_m, p.wall_thickness_m)
        t_r = float(curve.time_to_rupture(float(res.T_wo[k]), sigma))
        rate = 1.0 / t_r
        out.update({
            "converged": True, "Q_comb_W": flue.Q_comb_W, "flue_gas_kmol_h": flue.n_kmol_h, "T_fg_in_K": res.extras["T_fg_0plus_K"],
            "feed_per_tube_kmol_h": float(sum(feed.F.values())), "dry_feed_kmol_h": float(sum(feed.F.values()) - feed.F["H2O"]),
            "H2_out_kmol_h": float(tr.F["H2"][-1]), "H2_net_kmol_h": float(tr.F["H2"][-1] - feed.F.get("H2", 0.0)),
            "dry_H2_pct": 100 * float(dry["H2"][-1]), "CH4_slip_dry_pct": 100 * float(dry["CH4"][-1]),
            "CH4_conversion": float(tr.conversion_CH4[-1]), "T_out_K": float(tr.T[-1]), "P_out_bar": float(tr.P[-1]),
            "T_wo_max_K": float(res.T_wo[k]), "z_frac_T_wo_max": float(res.z[k] / res.z[-1]), "T_mid_max_K": float(T_mid.max()),
            "T_fg_out_K": float(res.T_fg[-1]), "T_fg_max_K": float(res.T_fg.max()),
            "duty_W": res.energy["tube_duty_W"], "thermal_efficiency": res.energy["tube_duty_W"] / res.energy["Q_eff_W"],
            "sigma_hot_MPa": sigma, "log10_t_r_hot": float(np.log10(t_r)), "life_consumption_rate_per_h": rate,
            "life_consumption_ratio_to_base": (rate / base_rate) if base_rate else np.nan,
            "energy_closure": res.energy["closure_rel_error"], "runtime_s": time.perf_counter() - t0, "error": "",
        })
    except Exception as e:  # noqa: BLE001 - record and continue
        out.update({"converged": False, "error": f"{type(e).__name__}: {e}", "runtime_s": time.perf_counter() - t0})
    return out


def base_run(base: Optional[BaseCase] = None, curve=None) -> Dict[str, float]:
    base = base or base_case()
    r = run_case(base.inputs(), base, curve)
    r["life_consumption_ratio_to_base"] = 1.0
    return r


# ---------------------------------------------------------------------------
# Closed loop: firing for a target outlet temperature
# ---------------------------------------------------------------------------
def solve_firing_for_outlet_T(x: Mapping[str, float], T_target: float = TARGET_T_OUT_K, base: Optional[BaseCase] = None,
                              curve=None, base_rate: Optional[float] = None, bracket=(0.6, 1.6), xtol: float = 1e-4) -> Dict[str, float]:
    """Adjust ``specific_firing_factor`` (brentq) so that the outlet process-gas temperature equals ``T_target``."""
    base = base or base_case()
    curve = curve or creep.active_curve()
    cache = {}

    def g(f):
        xx = dict(x); xx["specific_firing_factor"] = float(f)
        r = run_case(xx, base, curve, base_rate)
        cache[float(f)] = r
        if not r["converged"]:
            raise RuntimeError(r["error"])
        return r["T_out_K"] - T_target

    try:
        f_star = brentq(g, bracket[0], bracket[1], xtol=xtol)
        r = cache.get(float(f_star)) or run_case({**x, "specific_firing_factor": f_star}, base, curve, base_rate)
        r["closed_loop"] = True; r["T_target_K"] = T_target; r["firing_solver_evals"] = len(cache)
        return r
    except Exception as e:  # noqa: BLE001
        out = {k: float(x[k]) for k in INPUT_NAMES}
        out.update({"converged": False, "closed_loop": True, "T_target_K": T_target, "error": f"{type(e).__name__}: {e}"})
        return out


# ---------------------------------------------------------------------------
# Campaigns
# ---------------------------------------------------------------------------
def lhs_samples(n: int, seed: int = 0, ranges: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    lo, up = bounds_arrays(ranges)
    X = qmc.scale(qmc.LatinHypercube(d=len(INPUT_NAMES), seed=seed).random(n), lo, up)
    return pd.DataFrame(X, columns=INPUT_NAMES)


def _worker(x: Dict[str, float], base_rate: float) -> Dict[str, float]:
    return run_case(x, base_case(), creep.active_curve(), base_rate)


def run_batch(X: pd.DataFrame, base_rate: float, n_jobs: int = -1, closed_loop: bool = False,
              T_target: float = TARGET_T_OUT_K) -> pd.DataFrame:
    from joblib import Parallel, delayed
    recs = X.to_dict(orient="records")
    if closed_loop:
        rows = Parallel(n_jobs=n_jobs)(delayed(_worker_closed)(x, base_rate, T_target) for x in recs)
    else:
        rows = Parallel(n_jobs=n_jobs)(delayed(_worker)(x, base_rate) for x in recs)
    return pd.DataFrame(rows)


def _worker_closed(x, base_rate, T_target):
    return solve_firing_for_outlet_T(x, T_target, base_case(), creep.active_curve(), base_rate)


def git_hash() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def run_lhs_campaign(n: int = 2000, seed: int = 0, n_jobs: int = -1, name: str = "lhs_plantA_v1", ranges_csv: Optional[Path] = None) -> Dict[str, object]:
    """LHS open-loop campaign -> data/lhs_runs/<name>.csv.gz (+ metadata JSON, failed-runs CSV)."""
    LHS_DIR.mkdir(parents=True, exist_ok=True)
    ranges = load_ranges(ranges_csv); base = base_case(); b = base_run(base)
    X = lhs_samples(n, seed, ranges)
    t0 = time.perf_counter()
    df = run_batch(X, b["life_consumption_rate_per_h"], n_jobs)
    dt = time.perf_counter() - t0
    ok = df[df.converged == True]  # noqa: E712
    bad = df[df.converged != True]  # noqa: E712
    ok.to_csv(LHS_DIR / f"{name}.csv.gz", index=False, compression="gzip")
    bad[list(INPUT_NAMES) + ["error"]].to_csv(LHS_DIR / f"{name}_failed.csv", index=False)
    meta = {"name": name, "n_samples": int(n), "seed": seed, "sampler": "scipy.stats.qmc.LatinHypercube", "git_commit": git_hash(),
            "created": time.strftime("%Y-%m-%d %H:%M:%S"), "wall_time_s": dt, "n_jobs": n_jobs, "n_converged": int(len(ok)),
            "n_failed": int(len(bad)), "convergence_rate": float(len(ok) / n),
            "ranges_file": str(ranges_csv or RANGES_CSV), "ranges": {k: {"min": float(ranges.loc[k, "min"]), "max": float(ranges.loc[k, "max"]), "unit": str(ranges.loc[k, "unit"])} for k in INPUT_NAMES},
            "base_case": {**base.inputs(), "Q_comb_W": base.Q_comb_W, "life_consumption_rate_per_h": b["life_consumption_rate_per_h"],
                          "T_wo_max_K": b["T_wo_max_K"], "T_out_K": b["T_out_K"]},
            "life_curve": creep.active_curve().source or creep.alloy_key(), "outputs": [c for c in ok.columns if c not in INPUT_NAMES]}
    (LHS_DIR / f"{name}_meta.json").write_text(json.dumps(meta, indent=2))
    return {"data": ok, "failed": bad, "meta": meta}


def run_closed_loop_grid(n_sc: int = 15, n_load: int = 15, n_jobs: int = -1, name: str = "grid_sc_load_v1",
                         T_target: float = TARGET_T_OUT_K, ranges_csv: Optional[Path] = None) -> pd.DataFrame:
    """15 x 15 closed-loop grid over steam-to-carbon x load, other inputs at base."""
    ranges = load_ranges(ranges_csv); base = base_case(); b = base_run(base)
    sc = np.linspace(ranges.loc["steam_to_carbon", "min"], ranges.loc["steam_to_carbon", "max"], n_sc)
    ld = np.linspace(ranges.loc["feed_per_tube_fraction", "min"], ranges.loc["feed_per_tube_fraction", "max"], n_load)
    recs = []
    for s in sc:
        for l in ld:
            x = base.inputs(); x["steam_to_carbon"] = float(s); x["feed_per_tube_fraction"] = float(l); recs.append(x)
    t0 = time.perf_counter()
    df = run_batch(pd.DataFrame(recs), b["life_consumption_rate_per_h"], n_jobs, closed_loop=True, T_target=T_target)
    df.attrs["wall_time_s"] = time.perf_counter() - t0
    LHS_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(LHS_DIR / f"{name}.csv.gz", index=False, compression="gzip")
    return df


def use_config(cfg) -> None:
    """Point the module at a rdt.config.RunConfig (ranges file)."""
    global RANGES_CSV
    RANGES_CSV = Path(cfg.ranges_csv)
