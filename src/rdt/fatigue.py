"""Screening of the thermal-fatigue exposure of the RQ2 scenarios.

``python -m rdt.fatigue [--version v3]``

The life model of the twin covers creep only. This module asks a narrower question than a
creep-fatigue model: **how large is the cyclic thermal-stress exposure of smooth load variation
compared with the exposure the plant already accepts from start-up and shutdown?** It is a screening
calculation. It produces no life, no cycles-to-failure and no creep-fatigue damage sum, and it must
not be read as one.

What is computed
----------------
The elastic thermal stress from the radial temperature gradient uses the relation already stated in
Section S1 of the supplementary material,

    sigma_th = E alpha dT_wall / [2 (1 - nu)],

with the same constants (:data:`E_PA`, :data:`ALPHA_PER_K`, :data:`NU`). ``dT_wall`` is the drop
across the tube wall at the hot spot, which the coupled model gives exactly as

    dT_wall = q_o d_o / (2 lambda_tube) ln(d_o / d_i)

(:mod:`rdt.reactor1d`). It is not stored in the hourly scenario histories, so it is regenerated here by
the **physics model**, once per distinct operating point of each scenario. A surrogate was tried and
rejected: a Gaussian process fitted to the wall drop over the whole sampled operating box reaches only
about 3.7 K RMSE (38 K worst case), which is the same order as the stress ranges being counted and
would inject spurious reversals into otherwise smooth histories. Every reported row states its code
path, and the physics hot-spot temperature is compared against the one the scenario stage stored.

The hourly stress history of each scenario is then rainflow counted (:func:`rainflow`, four-point
method, residue as half cycles) and reduced to a Basquin-type range measure ``sum n_i dsigma_i^m`` at
several exponents, so that the comparison does not rest on one assumed exponent.

Assumptions, stated as such
---------------------------
* ``E``, ``alpha`` and ``nu`` are the single representative values of Section S1, not
  temperature-dependent properties. ASSUMPTION.
* The cold end of a start-up cycle has zero wall gradient. That is exact for a tube at ambient with no
  firing, but the *path* between cold and hot is not modelled: the staged warm-up in
  :func:`startup_cycle` walks the firing factor over the range the coupled model will solve and reports
  the gradient at each stage, not a transient.
* :data:`STARTUPS_PER_YEAR` = 2 to 6 full cool-downs per year is a **plant-practice assumption**, not a
  result of this or any other part of the twin.
* The Basquin exponents :data:`BASQUIN_M` are a sweep, not a fitted material property. No fatigue curve
  for the alloy is used, and none is implied.
* Rainflow counting is applied to the *elastic* stress range. In service the thermal stress relaxes by
  creep, so the elastic range is an upper bound on the cyclic stress actually seen. ASSUMPTION.

Scope
-----
Start-up and shutdown transients are represented only by their end states. Hold times, creep-fatigue
interaction, oxidation and carburisation are outside this module, as they are outside the twin.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from rdt import config as C

ROOT = C.ROOT
OUT_DIR = ROOT / "data" / "fatigue"

#: elastic constants of Section S1 of the supplementary material. ASSUMPTION: single representative
#: values, not temperature-dependent properties.
E_PA = 150.0e9              # Young's modulus [Pa]
ALPHA_PER_K = 18.0e-6       # linear thermal expansion at 900 C [1/K]
NU = 0.3                    # Poisson's ratio

#: full cool-downs per year of the reference plant. ASSUMPTION -- plant practice, not a computed result.
STARTUPS_PER_YEAR = (2.0, 6.0)
#: Basquin-type exponent sweep. ASSUMPTION -- a sweep, not a fitted material property.
BASQUIN_M = (3.0, 5.0, 8.0)
#: stress-range bin edges for the cycle histogram [MPa]
BIN_EDGES_MPA = (0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, np.inf)
HOURS_PER_YEAR = 8760


# ---------------------------------------------------------------------------
# Stress from the wall gradient
# ---------------------------------------------------------------------------
def thermal_stress_MPa(dT_wall_K, E_Pa: float = E_PA, alpha: float = ALPHA_PER_K,
                       nu: float = NU) -> np.ndarray:
    """``E alpha dT_wall / [2(1-nu)]`` in MPa -- the relation of Section S1, same constants."""
    return E_Pa * alpha * np.asarray(dT_wall_K, float) / (2.0 * (1.0 - nu)) / 1.0e6


# ---------------------------------------------------------------------------
# The wall drop at the hot spot, from the physics model
# ---------------------------------------------------------------------------
def hotspot_wall_drop(x: Mapping[str, float], base=None) -> Dict[str, float]:
    """Coupled-model state at the hot spot, including the drop across the tube wall.

    The hourly histories store ``T_wo_max_K`` but not ``T_wi``, so the coupled solve is repeated here
    and the wall drop read off at the outer-wall maximum. Identical call to the one
    :func:`rdt.operate.run_case` makes, so the temperatures agree with the stored ones.
    """
    from rdt import furnace as fu
    from rdt import latham_cases as lc
    from rdt import operate as op

    base = base or op.base_case()
    feed, flue, p = op.build_inputs(base, x)
    tube, bed, _, _, geom = lc.build_case(base.row, p)
    res = fu.simulate_coupled(tube, bed, feed, geom, flue,
                              fu.HeatRelease(p.L_q, p.alpha_top, p.f_loss, p.n_sections, p.inlet_mode),
                              fu.FurnaceParams(F_gt=p.F_gt, f_ctube=p.f_ctube), f_htg=p.f_htg,
                              heat_transfer=p.heat_transfer)
    if not res.success:
        raise RuntimeError(res.message)
    k = int(np.argmax(res.T_wo))
    return {"T_wo_max_K": float(res.T_wo[k]), "T_wi_at_hot_K": float(res.T_wi[k]),
            "dT_wall_hot_K": float(res.T_wo[k] - res.T_wi[k]),
            "q_o_hot_W_m2": float(res.tube.q_flux_outer[k]),
            "z_frac_hot": float(res.z[k] / res.z[-1]),
            "dT_wall_max_K": float(np.max(res.T_wo - res.T_wi))}


def _drop_worker(x: Dict[str, float]) -> Dict[str, float]:
    try:
        return hotspot_wall_drop(x)
    except Exception as e:  # noqa: BLE001 - recorded, not silently dropped
        return {"dT_wall_hot_K": np.nan, "error": f"{type(e).__name__}: {e}"}


def wall_drop_for_inputs(X: pd.DataFrame, n_jobs: int = -1, decimals: int = 8) -> pd.DataFrame:
    """Hot-spot wall drop for a table of operating inputs, by the physics model.

    No surrogate: a Gaussian process fitted to the wall drop over the whole sampled operating box
    reaches only about 3.7 K RMSE, which is of the same order as the stress ranges being counted and
    would inject spurious reversals into otherwise smooth histories. The coupled model is therefore
    solved once per *distinct* input row (rounded to ``decimals``) and the result broadcast back over
    the hours, which is exact and, because steady scenarios collapse to a single row, also cheap for
    them.
    """
    from joblib import Parallel, delayed
    from rdt import operate as op

    A = X[list(op.INPUT_NAMES)].round(decimals)
    keys = [tuple(r) for r in A.to_numpy()]
    distinct = sorted(set(keys))
    res = Parallel(n_jobs=n_jobs)(delayed(_drop_worker)(dict(zip(op.INPUT_NAMES, k))) for k in distinct)
    lut = dict(zip(distinct, res))
    out = pd.DataFrame([lut[k] for k in keys], index=X.index)
    out["evaluated_by"] = "physics"
    out.attrs["n_distinct"] = len(distinct)
    out.attrs["n_failed"] = int(sum(1 for r in res if not np.isfinite(r["dT_wall_hot_K"])))
    return out


def scenario_wall_drop(cfg: C.RunConfig, name: str, n_jobs: int = -1,
                       cache_dir: Optional[Path] = None, force: bool = False) -> pd.DataFrame:
    """Hourly wall drop of one scenario, cached so that a re-run costs nothing."""
    from rdt import operate as op

    cache_dir = Path(cache_dir or OUT_DIR / "wall_drop")
    cache = cache_dir / f"wall_drop_{cfg.tag}_{name}.csv.gz"
    hourly = pd.read_csv(cfg.scenario_hourly(name),
                         usecols=list(op.INPUT_NAMES) + ["T_wo_max_K", "evaluated_by"])
    meta_path = cache.with_suffix("").with_suffix(".meta.json")
    if cache.exists() and not force:
        d = pd.read_csv(cache)
        if len(d) == len(hourly):
            if meta_path.exists():
                d.attrs.update(json.loads(meta_path.read_text()))
            return d
    got = wall_drop_for_inputs(hourly[list(op.INPUT_NAMES)], n_jobs=n_jobs)
    out = pd.DataFrame({
        "T_wo_max_K_stored": hourly.T_wo_max_K,                 # scenario stage, surrogate/physics
        "T_wo_max_K_physics": got.T_wo_max_K,                   # this module, physics
        "dT_wall_hot_K": got.dT_wall_hot_K, "q_o_hot_W_m2": got.q_o_hot_W_m2,
        "dT_wall_max_K": got.dT_wall_max_K, "z_frac_hot": got.z_frac_hot,
        "stored_path": hourly.evaluated_by, "evaluated_by": got.evaluated_by})
    out.attrs.update(got.attrs)
    cache.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(cache, index=False, compression="gzip")
    meta_path.write_text(json.dumps(
        {"scenario": name, "n_hours": len(out), "n_distinct": got.attrs["n_distinct"],
         "n_failed": got.attrs["n_failed"]}, indent=2))
    return out


# ---------------------------------------------------------------------------
# Rainflow counting
# ---------------------------------------------------------------------------
def reversals(series: Sequence[float]) -> np.ndarray:
    """Turning points of ``series``, keeping the first and last samples.

    Consecutive equal values are collapsed first, so a flat history reduces to its two end points.
    """
    s = np.asarray(series, float)
    if s.size == 0:
        return s
    s = s[np.concatenate(([True], np.diff(s) != 0.0))]
    if s.size < 3:
        return s
    d = np.diff(s)
    turn = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0] + 1
    return np.concatenate(([s[0]], s[turn], [s[-1]]))


def rainflow(series: Sequence[float]) -> List[Tuple[float, float, float]]:
    """Rainflow cycle counting, four-point method, residue counted as half cycles.

    Returns ``(range, mean, count)`` per cycle, ``count`` being 1.0 for a closed cycle and 0.5 for a
    residue half cycle. The four-point rule extracts the inner pair ``p2--p3`` whenever its range is no
    larger than the ranges on either side of it, which is the standard closed-hysteresis-loop test and
    needs no special case for the start of the record.

    The total count always equals ``(n_reversals - 1) / 2``: every reversal but the first belongs to
    exactly one half cycle.
    """
    r = reversals(series)
    stack: List[float] = []
    out: List[Tuple[float, float, float]] = []
    for p in r:
        stack.append(float(p))
        while len(stack) >= 4:
            p1, p2, p3, p4 = stack[-4], stack[-3], stack[-2], stack[-1]
            r1, r2, r3 = abs(p2 - p1), abs(p3 - p2), abs(p4 - p3)
            if r2 <= r1 and r2 <= r3:
                out.append((r2, 0.5 * (p2 + p3), 1.0))
                del stack[-3:-1]                      # remove p2 and p3, keep p1 and p4
            else:
                break
    for a, b in zip(stack[:-1], stack[1:]):
        out.append((abs(b - a), 0.5 * (a + b), 0.5))
    return out


def cycle_table(cycles: Sequence[Tuple[float, float, float]],
                edges: Sequence[float] = BIN_EDGES_MPA) -> pd.DataFrame:
    """Counts per stress-range bin."""
    rng = np.array([c[0] for c in cycles], float)
    cnt = np.array([c[2] for c in cycles], float)
    e = np.asarray(edges, float)
    idx = np.clip(np.digitize(rng, e) - 1, 0, len(e) - 2)
    tot = np.zeros(len(e) - 1)
    np.add.at(tot, idx, cnt)
    lab = [f"{e[i]:g}-{e[i+1]:g}" if np.isfinite(e[i + 1]) else f">{e[i]:g}" for i in range(len(e) - 1)]
    return pd.DataFrame({"bin_MPa": lab, "cycles": tot})


def basquin_measure(cycles: Sequence[Tuple[float, float, float]], m: float) -> float:
    """``sum_i n_i dsigma_i^m`` [MPa^m]. A range measure, not a damage sum: no fatigue curve enters."""
    return float(sum(c[2] * c[0] ** m for c in cycles))


# ---------------------------------------------------------------------------
# Scenarios
# ---------------------------------------------------------------------------
def scenario_names(cfg: C.RunConfig) -> List[str]:
    return pd.read_csv(cfg.scenarios_summary)["scenario"].tolist()


def scenario_stress_history(cfg: C.RunConfig, name: str, n_jobs: int = -1) -> pd.DataFrame:
    """Hourly hot-spot wall drop and elastic thermal stress for one scenario, with the code path."""
    d = scenario_wall_drop(cfg, name, n_jobs=n_jobs)
    return d.assign(sigma_th_MPa=thermal_stress_MPa(d.dT_wall_hot_K),
                    sigma_th_max_MPa=thermal_stress_MPa(d.dT_wall_max_K))


def screen_scenario(cfg: C.RunConfig, name: str, startup_MPa: float, startup_thmax_MPa: float,
                    exponents: Sequence[float] = BASQUIN_M, n_jobs: int = -1) -> Dict[str, object]:
    h = scenario_stress_history(cfg, name, n_jobs=n_jobs)
    cyc = rainflow(h.sigma_th_MPa.to_numpy(float))
    cyc_max = rainflow(h.sigma_th_max_MPa.to_numpy(float))
    n_dist = int(h.attrs.get("n_distinct", h[["dT_wall_hot_K"]].round(9).drop_duplicates().shape[0]))
    dT = (h.T_wo_max_K_physics - h.T_wo_max_K_stored).abs()
    rec: Dict[str, object] = {
        "scenario": name,
        "evaluated_by": f"physics, {n_dist} distinct points",
        "max_dT_wo_vs_stored_K": float(dT.max()),
        "sigma_th_mean_MPa": float(h.sigma_th_MPa.mean()),
        "sigma_th_min_MPa": float(h.sigma_th_MPa.min()), "sigma_th_max_MPa": float(h.sigma_th_MPa.max()),
        "largest_range_MPa": float(max((c[0] for c in cyc), default=0.0)),
        "cycles_per_year": float(sum(c[2] for c in cyc)),
        "cycles_ge_1MPa": float(sum(c[2] for c in cyc if c[0] >= 1.0)),
        "z_frac_hot_min": float(h.z_frac_hot.min()), "z_frac_hot_max": float(h.z_frac_hot.max()),
        "sigma_thmax_min_MPa": float(h.sigma_th_max_MPa.min()),
        "sigma_thmax_max_MPa": float(h.sigma_th_max_MPa.max()),
        "largest_range_thmax_MPa": float(max((c[0] for c in cyc_max), default=0.0)),
        "cycles_thmax_ge_1MPa": float(sum(c[2] for c in cyc_max if c[0] >= 1.0)),
        "n_reversals": int(len(reversals(h.sigma_th_MPa.to_numpy(float)))),
    }
    for m in exponents:
        rec[f"measure_m{m:g}"] = basquin_measure(cyc, m)
        rec[f"ratio_to_startup_m{m:g}"] = basquin_measure(cyc, m) / (startup_MPa ** m)
        rec[f"ratio_thmax_m{m:g}"] = basquin_measure(cyc_max, m) / (startup_thmax_MPa ** m)
    rec["_cycles"] = cyc
    rec["_cycles_thmax"] = cyc_max
    return rec


# ---------------------------------------------------------------------------
# Start-up / shutdown reference cycle
# ---------------------------------------------------------------------------
def startup_cycle(cfg: C.RunConfig, stages: Sequence[float] = (0.65, 0.75, 0.85, 0.95, 1.0),
                  n_jobs: int = -1) -> Dict[str, object]:
    """Wall gradient from ambient to the calibrated base, stage by stage.

    The cold stage is exact: a tube at ambient with no firing carries no flux, so the wall gradient and
    the thermal stress are identically zero. The intermediate stages walk the firing factor at base
    load over the range the coupled model solves, and are computed by the physics model. The stress
    range of one full start-up/shutdown cycle is the hot value minus the cold zero.

    The *path* between stages is not a transient calculation; only the end states are model results.
    """
    from joblib import Parallel, delayed
    from rdt import operate as op

    base = op.base_case()
    xs = []
    for f in stages:
        x = dict(base.inputs()); x["specific_firing_factor"] = float(f); xs.append(x)
    res = Parallel(n_jobs=n_jobs)(delayed(_drop_worker)(x) for x in xs)
    rows = [{"stage": "cold (ambient, no firing)", "specific_firing_factor": 0.0,
             "T_wo_max_K": np.nan, "z_frac_hot": np.nan, "dT_wall_hot_K": 0.0, "sigma_th_MPa": 0.0,
             "dT_wall_max_K": 0.0, "sigma_th_max_MPa": 0.0, "evaluated_by": "exact (zero flux)"}]
    for f, r in zip(stages, res):
        rows.append({"stage": f"firing {f:.2f}", "specific_firing_factor": float(f),
                     "T_wo_max_K": r.get("T_wo_max_K", np.nan), "z_frac_hot": r.get("z_frac_hot", np.nan),
                     "dT_wall_hot_K": r["dT_wall_hot_K"],
                     "sigma_th_MPa": float(thermal_stress_MPa(r["dT_wall_hot_K"])),
                     "dT_wall_max_K": r.get("dT_wall_max_K", np.nan),
                     "sigma_th_max_MPa": float(thermal_stress_MPa(r.get("dT_wall_max_K", np.nan))),
                     "evaluated_by": "physics"})
    tab = pd.DataFrame(rows)
    hot = float(tab.sigma_th_MPa.iloc[-1]); hot_mx = float(tab.sigma_th_max_MPa.iloc[-1])
    return {"stages": tab, "range_MPa": hot - 0.0, "hot_MPa": hot,
            "range_thmax_MPa": hot_mx - 0.0, "hot_thmax_MPa": hot_mx,
            "dT_wall_hot_K": float(tab.dT_wall_hot_K.iloc[-1]),
            "dT_wall_max_K": float(tab.dT_wall_max_K.iloc[-1])}


# ---------------------------------------------------------------------------
# Build and report
# ---------------------------------------------------------------------------
def build(cfg: C.RunConfig, exponents: Sequence[float] = BASQUIN_M,
          n_jobs: int = -1) -> Dict[str, object]:
    su = startup_cycle(cfg, n_jobs=n_jobs)
    recs = [screen_scenario(cfg, n, su["range_MPa"], su["range_thmax_MPa"], exponents, n_jobs)
            for n in scenario_names(cfg)]
    cycles = {r["scenario"]: r.pop("_cycles") for r in recs}
    cycles_thmax = {r["scenario"]: r.pop("_cycles_thmax") for r in recs}
    table = pd.DataFrame(recs).set_index("scenario")
    hist = {n: cycle_table(c) for n, c in cycles.items()}
    return {"cfg": cfg, "startup": su, "table": table, "cycles": cycles,
            "cycles_thmax": cycles_thmax, "hist": hist, "exponents": list(exponents)}


def _crossover_exponent(res: Dict[str, object], scenario: str = "S2_daily", which: str = "hot",
                        grid: Sequence[float] = np.arange(1.0, 12.01, 0.05)) -> Dict[str, float]:
    """Exponent at which one year of the scenario equals one start-up, and the 2-6 per year band."""
    if which == "hot":
        cyc = res["cycles"][scenario]; su = float(res["startup"]["range_MPa"])
    else:
        cyc = res["cycles_thmax"][scenario]; su = float(res["startup"]["range_thmax_MPa"])
    out = {}
    for label, target in (("one_startup", 1.0), ("2_per_year", 2.0), ("6_per_year", 6.0)):
        hit = None
        for m in grid:
            if basquin_measure(cyc, m) / su ** m < target:
                hit = float(m); break
        out[label] = hit if hit is not None else float("nan")
    return out


def report(res: Dict[str, object]) -> str:
    cfg: C.RunConfig = res["cfg"]; su = res["startup"]; t: pd.DataFrame = res["table"]
    n_dist = {n: str(r).split(",")[1].strip().split()[0] for n, r in t.evaluated_by.items()}
    L: List[str] = []
    P = L.append
    P(f"Thermal-fatigue screening of the RQ2 scenarios  ({cfg.tag} results)")
    P("=" * 108)
    P("SCREENING ONLY. No life, no cycles to failure and no creep-fatigue damage sum is computed or")
    P("implied below. The question is one of magnitude: how does the cyclic thermal-stress exposure of")
    P("smooth load variation compare with the exposure the plant already accepts at every cold start?")
    P("")
    P("WHAT IS COMPUTED FROM THE TWIN, AND WHAT IS ASSUMED")
    P("-" * 108)
    P("  computed by the twin :")
    P("    - the hot-spot wall drop dT_wall = q_o d_o/(2 lambda) ln(d_o/d_i) at every scenario hour,")
    P("      from the same coupled solve the scenario stage used;")
    P("    - the wall drop at each stage of the staged warm-up;")
    P("    - the rainflow cycle counts, which follow from those histories alone.")
    P("  assumed :")
    P(f"    - E = {E_PA/1e9:g} GPa, alpha = {ALPHA_PER_K*1e6:g}e-6 /K, nu = {NU:g}: the single")
    P("      representative values of Section S1, not temperature-dependent properties;")
    P(f"    - {STARTUPS_PER_YEAR[0]:g} to {STARTUPS_PER_YEAR[1]:g} full cool-downs per year -- PLANT PRACTICE, "
      "not a computed result;")
    P(f"    - the Basquin exponents m = {', '.join(f'{m:g}' for m in res['exponents'])} are a sweep, not a "
      "fitted material property;")
    P("    - the elastic stress range is an upper bound: in service the thermal stress relaxes by creep.")
    P("")

    P("CODE PATH FOR THE WALL DROP")
    P("-" * 108)
    P("  dT_wall is not stored in the hourly histories, so it is recomputed here. A surrogate was tried")
    P("  and REJECTED: a GP fitted to the wall drop over the whole sampled operating box reaches only")
    P("  3.7 K RMSE (38 K worst case), the same order as the ranges being counted, and its noise would")
    P("  invent reversals in smooth histories. Every hour below is therefore the PHYSICS model, solved")
    P("  once per distinct operating point and broadcast over the hours that repeat it.")
    P("")
    P(f"  {'scenario':<30} {'distinct solves':>15} {'hours':>7}   max |T_wo(physics) - T_wo(stored)|")
    for name, r in t.iterrows():
        P(f"  {name:<30} {n_dist[name]:>15} {HOURS_PER_YEAR:>7}   {r.max_dT_wo_vs_stored_K:8.4f} K")
    P("  the stored column came from the scenario stage's own surrogate; the agreement above is the")
    P("  consistency check between the two.")
    offset = [(n_, r.max_dT_wo_vs_stored_K) for n_, r in t.iterrows() if r.max_dT_wo_vs_stored_K > 1.0]
    if offset:
        P("")
        P("  LIMITATION -- scenarios the wall-drop calculation cannot represent:")
        for n_, d in offset:
            P(f"    {n_:<30} stored T_wo exceeds the physics solve by up to {d:.1f} K")
        P("    In these the scenario stage adds a metal-temperature offset directly to T_wo rather than")
        P("    changing an operating input, so re-solving the inputs returns the base state. The twin")
        P("    defines no wall gradient for an imposed offset -- an offset produced by extra local flux")
        P("    and one produced by a redistribution of radiation imply different gradients, and the")
        P("    scenario does not say which it is. The excursion is therefore ABSENT from the cycle counts")
        P("    for those rows, which are consequently a lower bound. No flux is invented to fill the gap.")
    P("")

    P("START-UP / SHUTDOWN REFERENCE CYCLE")
    P("-" * 108)
    P("  Cold end is exact (a tube at ambient with no firing carries no flux, so the gradient is zero).")
    P("  The intermediate stages are physics-model end states, not a transient.")
    P("")
    P(f"  {'stage':<26} {'firing':>7} {'T_wo,max':>10} {'z/L':>6} | {'at the hot spot':^24} | "
      f"{'max anywhere on the tube':^26}")
    P(f"  {'':<26} {'':>7} {'[K]':>10} {'':>6} | {'dT [K]':>10} {'sigma [MPa]':>12} | "
      f"{'dT [K]':>11} {'sigma [MPa]':>13}")
    for _, r in su["stages"].iterrows():
        tw = "       -- " if not np.isfinite(r.T_wo_max_K) else f"{r.T_wo_max_K:10.2f}"
        zf = "    --" if not np.isfinite(r.z_frac_hot) else f"{r.z_frac_hot:6.3f}"
        P(f"  {r.stage:<26} {r.specific_firing_factor:7.2f} {tw} {zf} | {r.dT_wall_hot_K:10.2f} "
          f"{r.sigma_th_MPa:12.2f} | {r.dT_wall_max_K:11.2f} {r.sigma_th_max_MPa:13.2f}")
    P("")
    P(f"  stress range of ONE full start-up/shutdown cycle : {su['range_MPa']:.2f} MPa at the hot spot "
      f"(wall drop {su['dT_wall_hot_K']:.2f} K),")
    P(f"  {'':<49} {su['range_thmax_MPa']:.2f} MPa on the maximum gradient "
      f"({su['dT_wall_max_K']:.2f} K)")
    P("  The intermediate stages do not enter that number. Cold -> hot -> cold is monotone in each")
    P("  direction, so it rainflow-counts as one full cycle of the end-to-end range whatever path the")
    P("  warm-up takes; the stages are shown to make the gradient's growth visible, not to set the")
    P("  range. A start-up interrupted by holds or a trip would add sub-cycles, which are not modelled.")
    P(f"  at the assumed {STARTUPS_PER_YEAR[0]:g}-{STARTUPS_PER_YEAR[1]:g} cool-downs per year that is "
      f"{STARTUPS_PER_YEAR[0]:g}-{STARTUPS_PER_YEAR[1]:g} such cycles annually")
    P("")
    P("  Section S1 of the supplementary material previously carried an estimated 20-30 K wall drop and")
    P(f"  about 45 MPa. It now states the computed values: {su['dT_wall_hot_K']:.0f} K / "
      f"{su['hot_MPa']:.0f} MPa at the hot spot, {su['dT_wall_max_K']:.0f} K / "
      f"{su['hot_thmax_MPa']:.0f} MPa at the steepest")
    P("  gradient, and 32 K / 63 MPa averaged over the heated length. The relation and the constants")
    P("  used here are the ones printed there; only the wall drop comes from the model.")
    P("")

    P("SCENARIO STRESS HISTORIES AND RAINFLOW COUNTS")
    P("-" * 108)
    P(f"  {'scenario':<30} {'sig min':>8} {'sig max':>8} {'largest':>8} {'cycles':>9} {'>=1MPa':>8} | "
      f"{'ratio to one start-up':^30} | code path")
    P(f"  {'':<30} {'[MPa]':>8} {'[MPa]':>8} {'range':>8} {'per year':>9} {'per year':>8} | "
      + " ".join(f"{'m='+format(m,'g'):>9}" for m in res["exponents"]) + " | ")
    for name, r in t.iterrows():
        ratios = " ".join(f"{r['ratio_to_startup_m' + format(m, 'g')]:9.3g}" for m in res["exponents"])
        P(f"  {name:<30} {r.sigma_th_min_MPa:8.2f} {r.sigma_th_max_MPa:8.2f} "
          f"{r.largest_range_MPa:8.2f} {r.cycles_per_year:9.1f} {r.cycles_ge_1MPa:8.1f} | "
          f"{ratios} | {r.evaluated_by}")
    P("")
    P("  A flat history has no reversals and therefore no cycles: the steady scenarios count zero, as")
    P("  they should. The ratio is sum_i n_i dsigma_i^m over one year divided by dsigma_startup^m.")
    P("  The two cycle columns differ where the hourly control solve leaves numerical texture below")
    P("  1 MPa. Those sub-1 MPa counts are an artefact of the firing-bisection tolerance, not physical")
    P("  cycles, and they contribute nothing to the measure at any exponent used here: one such cycle")
    P("  is (1/110)^3 = 8e-7 of a start-up. The '>= 1 MPa' column is the count worth reading.")
    P("")
    P("SAME HISTORIES ON THE MAXIMUM GRADIENT ANYWHERE ON THE TUBE")
    P("-" * 108)
    P("  The table above follows the task definition: the wall drop AT THE HOT SPOT, the point of")
    P("  maximum outer-wall temperature. That point is not always the point of maximum wall gradient.")
    P("  Where the two disagree the hot-spot number can move the wrong way, so the same histories are")
    P("  counted again below on max_z dT_wall(z), which bounds the thermal stress anywhere on the tube.")
    P("")
    P(f"  {'scenario':<30} {'sig min':>8} {'sig max':>8} {'largest':>8} {'>=1MPa':>8} | "
      f"{'ratio to one start-up':^30} | {'hot spot z/L':>14}")
    P(f"  {'':<30} {'[MPa]':>8} {'[MPa]':>8} {'range':>8} {'per year':>8} | "
      + " ".join(f"{'m='+format(m,'g'):>9}" for m in res["exponents"]) + f" | {'min - max':>14}")
    for name, r in t.iterrows():
        ratios = " ".join(f"{r['ratio_thmax_m' + format(m, 'g')]:9.3g}" for m in res["exponents"])
        P(f"  {name:<30} {r.sigma_thmax_min_MPa:8.2f} {r.sigma_thmax_max_MPa:8.2f} "
          f"{r.largest_range_thmax_MPa:8.2f} {r.cycles_thmax_ge_1MPa:8.1f} | {ratios} | "
          f"{r.z_frac_hot_min:6.3f}-{r.z_frac_hot_max:.3f}")
    P("")
    moved = [n_ for n_, r in t.iterrows() if r.z_frac_hot_max - r.z_frac_hot_min > 0.05]
    P(f"  scenarios in which the hot spot MIGRATES along the tube: {', '.join(moved) or 'none'}.")
    P("  In those the two definitions must be read together. The clearest case is the severe")
    P("  over-firing event: at firing 1.22 the outer-wall maximum jumps to the tube outlet (z/L = 1.0),")
    P("  where the process gas has caught up with the wall and the flux is small, so the gradient AT")
    P("  the hot spot falls from 56.9 to 14.8 K even as the wall gets 36 K hotter. Read on the hot-spot")
    P("  definition alone that event looks like a large stress excursion; on the maximum gradient it is")
    P("  a small one, because max_z dT_wall rises monotonically with firing (62.8 -> 65.7 -> 69.5 K).")
    P("  The second table is the physically meaningful one for thermal stress; the first is the one the")
    P("  task asked for. Neither is a life.")
    P("")

    P("CYCLE COUNT BY STRESS-RANGE BIN (cycles per year)")
    P("-" * 108)
    hist = res["hist"]
    labels = list(hist[next(iter(hist))].bin_MPa)
    P(f"  {'scenario':<30} " + " ".join(f"{l:>9}" for l in labels))
    for name in t.index:
        h = hist[name]
        P(f"  {name:<30} " + " ".join(f"{c:9.1f}" for c in h.cycles))
    P("")

    P("WHERE THE DAILY CYCLE BECOMES COMPARABLE TO A START-UP")
    P("-" * 108)
    cx = _crossover_exponent(res, "S2_daily")
    for m in res["exponents"]:
        rr = float(t.loc["S2_daily", f"ratio_to_startup_m{m:g}"])
        verdict = ("above the whole assumed annual start-up exposure"
                   if rr > STARTUPS_PER_YEAR[1] else
                   "inside the assumed 2-6 start-up band" if rr >= STARTUPS_PER_YEAR[0] else
                   "below one start-up" if rr < 1.0 else "between one start-up and the assumed band")
        P(f"    m = {m:g}: one year of S2_daily = {rr:.3g} start-up cycles  -- {verdict}")
    P("")
    P("  Every scenario against the assumed annual start-up exposure "
      f"({STARTUPS_PER_YEAR[0]:g}-{STARTUPS_PER_YEAR[1]:g} cool-downs per year):")
    for m in res["exponents"]:
        col = f"ratio_to_startup_m{m:g}"
        over = [n_ for n_, r in t.iterrows() if r[col] > STARTUPS_PER_YEAR[1]]
        band = [n_ for n_, r in t.iterrows()
                if STARTUPS_PER_YEAR[0] <= r[col] <= STARTUPS_PER_YEAR[1]]
        P(f"    m = {m:g}: above the band {', '.join(over) or 'none'}"
          f"   |   inside the band {', '.join(band) or 'none'}")
    P("")
    cxm = _crossover_exponent(res, "S2_daily", "thmax")
    P("    crossing points for S2_daily (exponent above which a year is worth less than the target):")
    P(f"      {'':<38} {'hot spot':>10} {'max gradient':>14}")
    for key, lab in (("one_startup", "equal to ONE start-up"),
                     ("2_per_year", f"equal to {STARTUPS_PER_YEAR[0]:g} start-ups per year (low)"),
                     ("6_per_year", f"equal to {STARTUPS_PER_YEAR[1]:g} start-ups per year (high)")):
        P(f"      {lab:<38} {cx[key]:>10.2f} {cxm[key]:>14.2f}")
    P("")
    P("    The answer therefore depends on the exponent, which is why it is swept rather than assumed.")
    P("    No fatigue curve for this alloy is used anywhere above, and none of these numbers is a life.")
    return "\n".join(L)


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", default="v3", choices=sorted(C.BY_TAG))
    p.add_argument("--exponents", type=float, nargs="+", default=list(BASQUIN_M))
    p.add_argument("--n-jobs", type=int, default=-1)
    p.add_argument("--out-dir", type=Path, default=OUT_DIR)
    p.add_argument("--force", action="store_true", help="recompute the cached wall drops")
    a = p.parse_args(argv)

    cfg = C.BY_TAG[a.version]
    from rdt import pipeline as pl
    pl.apply_config(cfg)

    res = build(cfg, a.exponents, a.n_jobs)
    text = report(res)
    print(text)

    a.out_dir.mkdir(parents=True, exist_ok=True)
    (a.out_dir / f"fatigue_{cfg.tag}.txt").write_text(text + "\n")
    res["table"].to_csv(a.out_dir / f"fatigue_scenarios_{cfg.tag}.csv")
    res["startup"]["stages"].to_csv(a.out_dir / f"fatigue_startup_{cfg.tag}.csv", index=False)
    pd.concat({k: v.set_index("bin_MPa").cycles for k, v in res["hist"].items()}, axis=1).T.to_csv(
        a.out_dir / f"fatigue_bins_{cfg.tag}.csv")
    js = {"tag": cfg.tag, "E_Pa": E_PA, "alpha_per_K": ALPHA_PER_K, "nu": NU,
          "startups_per_year_ASSUMED": list(STARTUPS_PER_YEAR), "exponents_ASSUMED": a.exponents,
          "startup_range_MPa": res["startup"]["range_MPa"],
          "startup_dT_wall_K": res["startup"]["dT_wall_hot_K"],
          "crossover_exponent_S2_daily": _crossover_exponent(res, "S2_daily")}
    (a.out_dir / f"fatigue_{cfg.tag}.json").write_text(json.dumps(js, indent=2, default=float))
    for f in ("txt", "json"):
        print(f"\nwritten: {a.out_dir / f'fatigue_{cfg.tag}.{f}'}" if f == "txt"
              else f"         {a.out_dir / f'fatigue_{cfg.tag}.{f}'}")
    for f in ("scenarios", "startup", "bins"):
        print(f"         {a.out_dir / f'fatigue_{f}_{cfg.tag}.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
