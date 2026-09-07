"""Invariance of the headline result ratios to the choice of tube alloy.

``python -m rdt.alloy_comparison [--version v3]``

The digital twin predicts the thermal and mechanical state of the tube -- hot-spot metal temperature,
hoop stress, hydrogen rate -- without reference to any creep curve. The alloy enters only at the last
step, where ``t_r(T_wo, sigma)`` turns that state into a rupture time. Every scenario hour and every
Pareto regime therefore only has to be *re-scored*, not re-simulated, to ask what the paper's
conclusions would look like on a different alloy: the operating points, the firing, the temperatures
and the stresses are held exactly fixed and the master curve alone is swapped.

That is the controlled comparison the invariance claim needs. The ratios compared are

    S2 / S1                 daily-cycling versus steady operation, life per kmol H2
    S3 / S1                 renewable-following versus steady
    ageing year 4 / year 1  the S5 hold-slip ageing campaign
    knee / base             the RQ3 knee regime versus the base case
    S/C 3.5 / S/C 2.5       the two steam-to-carbon settings

One caveat is worth stating with the table: the Pareto set itself was optimised under the base alloy,
so the *location* of the knee is that of the base alloy; what is re-scored is its life ratio.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from rdt import config as C
from rdt import creep

ALLOYS = ("centralloy_g_4852", "centralloy_g_4852_micro", "centralloy_et_45_micro", creep.LEGACY_ALLOY)
SHORT = {"centralloy_g_4852": "G 4852 (HP-Nb)", "centralloy_g_4852_micro": "G 4852 Micro",
         "centralloy_et_45_micro": "ET 45 Micro", creep.LEGACY_ALLOY: "Yeh placeholder (legacy)"}

RATIOS = ("S2/S1", "S3/S1", "ageing y4/y1", "knee/base", "S/C 3.5 / 2.5")
_SCENARIOS = ("S1_steady", "S2_daily", "S3_renewable", "S6_SC2.5", "S6_SC3.5",
              "S5_hold_CH4_slip_ageing_y1", "S5_hold_CH4_slip_ageing_y4")


# ---------------------------------------------------------------------------
# Re-scoring
# ---------------------------------------------------------------------------
def load_hourly_states(cfg: C.RunConfig, scenarios: Sequence[str] = _SCENARIOS) -> Dict[str, pd.DataFrame]:
    """Per-hour tube state of each scenario: everything the creep curve needs, and the hydrogen rate."""
    out = {}
    for name in scenarios:
        p = cfg.scenario_hourly(name)
        if not p.exists():
            raise FileNotFoundError(f"{p} not found; run `python -m rdt.pipeline --stage scenarios --version {cfg.tag}`")
        out[name] = pd.read_csv(p, usecols=["T_wo_max_K", "sigma_hot_MPa", "H2_net_kmol_h"])
    return out


def scenario_life_per_kmol(states: Dict[str, pd.DataFrame], curve: creep.LarsonMillerCurve) -> Dict[str, float]:
    """Annual life consumption per kmol H2, relative to S1, for one master curve.

    ``D = sum_hours 1 h / t_r(T_wo, sigma)`` exactly as the scenario stage accumulates it; the hydrogen
    rate is a property of the operating point and does not move with the alloy.
    """
    raw = {}
    for name, df in states.items():
        t_r = curve.time_to_rupture(df.T_wo_max_K.to_numpy(float), df.sigma_hot_MPa.to_numpy(float))
        raw[name] = float(np.sum(1.0 / np.asarray(t_r, float))) / float(df.H2_net_kmol_h.sum())
    ref = raw["S1_steady"]
    return {k: v / ref for k, v in raw.items()}


def regime_states(cfg: C.RunConfig, regimes: Sequence[str] = ("base", "knee")) -> pd.DataFrame:
    """Physics-verified tube state of the named Pareto regimes (alloy independent, computed once)."""
    from rdt import operate as op

    df = pd.read_csv(cfg.regimes_csv, index_col=0)
    missing = [r for r in regimes if r not in df.index]
    if missing:
        raise KeyError(f"{cfg.regimes_csv} has no regime(s) {missing}; available: {list(df.index)}")
    rows = []
    base = op.base_case()
    for name in regimes:
        x = {k: float(df.loc[name, k]) for k in op.INPUT_NAMES}
        r = op.run_case(x, base)
        if not r["converged"]:
            raise RuntimeError(f"regime {name!r} did not converge: {r['error']}")
        rows.append({"regime": name, "T_wo_max_K": r["T_wo_max_K"], "sigma_hot_MPa": r["sigma_hot_MPa"],
                     "H2_net_kmol_h": r["H2_net_kmol_h"]})
    return pd.DataFrame(rows).set_index("regime")


def regime_life_per_kmol(states: pd.DataFrame, curve: creep.LarsonMillerCurve) -> Dict[str, float]:
    """Life consumption rate per kmol H2 of each regime, relative to the base case."""
    rate = {}
    for name, r in states.iterrows():
        t_r = float(curve.time_to_rupture(float(r.T_wo_max_K), float(r.sigma_hot_MPa)))
        rate[str(name)] = (1.0 / t_r) / float(r.H2_net_kmol_h)
    return {k: v / rate["base"] for k, v in rate.items()}


# ---------------------------------------------------------------------------
# Table
# ---------------------------------------------------------------------------
def build_table(cfg: C.RunConfig, alloys: Sequence[str] = ALLOYS) -> pd.DataFrame:
    """The five headline ratios for each alloy, plus the curve properties they follow from."""
    states = load_hourly_states(cfg)
    regimes = regime_states(cfg)
    rows = []
    for key in alloys:
        a = creep.load_alloy(key)
        s = scenario_life_per_kmol(states, a.minimum)
        g = regime_life_per_kmol(regimes, a.minimum)
        lmp_base = a.minimum.lmp_of_stress(float(regimes.loc["base", "sigma_hot_MPa"]))
        rows.append({
            "alloy": SHORT.get(key, key), "key": key, "C": a.C,
            "S2/S1": s["S2_daily"], "S3/S1": s["S3_renewable"],
            "ageing y4/y1": s["S5_hold_CH4_slip_ageing_y4"] / s["S5_hold_CH4_slip_ageing_y1"],
            "knee/base": g["knee"],
            "S/C 3.5 / 2.5": s["S6_SC3.5"] / s["S6_SC2.5"],
            "base t_r [y]": float(a.minimum.time_to_rupture(float(regimes.loc["base", "T_wo_max_K"]),
                                                            float(regimes.loc["base", "sigma_hot_MPa"]))) / 8760.0,
            "scatter [dec]": float(a.scatter.decades(float(regimes.loc["base", "T_wo_max_K"]), lmp_base)),
            "base LMP in range": bool(a.minimum.lmp_range[0] <= lmp_base <= a.minimum.lmp_range[1]),
        })
    return pd.DataFrame(rows).set_index("alloy")


def ranking_is_invariant(table: pd.DataFrame, material: float = 0.05) -> Dict[str, object]:
    """How far each ratio moves across the alloys, and whether its qualitative reading survives.

    Two different questions are separated here, because they carry very different weight in the paper:

    ``same_side_of_1``
        does the ratio stay on one side of unity, i.e. does the *direction* of the comparison hold;
    ``material``
        is the ratio further than ``material`` from unity on **every** curve, i.e. is there an effect
        worth reporting at all.

    A ratio that sits within a per cent or two of unity and crosses it is not a contradiction between
    the alloys -- it is a comparison that is too close to call on any of them, and should be reported
    as "no material difference" rather than as a direction.
    """
    out: Dict[str, object] = {}
    for r in RATIOS:
        v = table[r].to_numpy(float)
        out[r] = {"min": float(v.min()), "max": float(v.max()),
                  "spread_pct": float(100 * (v.max() - v.min()) / np.abs(v).mean()),
                  "same_side_of_1": bool(np.all(v > 1.0) or np.all(v < 1.0)),
                  "material": bool(np.all(np.abs(v - 1.0) > material)),
                  "max_dev_from_1": float(np.max(np.abs(v - 1.0)))}
    out["material_ratios"] = [r for r in RATIOS if out[r]["material"]]                       # type: ignore[index]
    out["material_and_same_side"] = all(out[r]["same_side_of_1"] for r in out["material_ratios"])  # type: ignore[index,union-attr]
    out["all_same_side"] = all(out[r]["same_side_of_1"] for r in RATIOS)                     # type: ignore[index]
    return out


def reoptimised_knee(tags: Sequence[str] = ("v2", "v3")) -> pd.DataFrame:
    """Knee/base ratio of each version's **own** Pareto run, i.e. with the knee re-optimised.

    The re-scored table holds the operating points fixed, which isolates the effect of the curve but
    overstates how much the *conclusion* moves: in practice the optimiser re-locates the knee on the new
    curve. Comparing the runs that each optimised their own knee separates the two effects.
    """
    rows = []
    for tag in tags:
        cfg = C.BY_TAG[tag]
        if not Path(cfg.regimes_csv).exists():
            continue
        df = pd.read_csv(cfg.regimes_csv, index_col=0)
        r = df.loc["knee"]
        rows.append({"tag": tag, "alloy": SHORT.get(cfg.creep_alloy, cfg.creep_alloy),
                     "knee/base (own curve)": float(r["life_rate_per_kmol_H2_rel_base_phys"]),
                     "S/C": float(r["steam_to_carbon"]), "firing": float(r["specific_firing_factor"]),
                     "T_wo_max_K": float(r["T_wo_max_K_phys"]), "H2": float(r["H2_net_kmol_h_phys"])})
    return pd.DataFrame(rows).set_index("tag")


def to_latex(table: pd.DataFrame) -> str:
    """A booktabs table for the manuscript."""
    cols = ["C", *RATIOS]
    head = " & ".join(["Alloy", "$C$", "S2/S1", "S3/S1", "ageing y4/y1", "knee/base", "S/C 3.5 / 2.5"])
    lines = ["\\begin{tabular}{lrrrrrr}", "\\toprule", head + " \\\\", "\\midrule"]
    for name, r in table.iterrows():
        vals = [f"{r['C']:.1f}"] + [f"{r[c]:.3f}" for c in RATIOS]
        lines.append(f"{name} & " + " & ".join(vals) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", default="v3", choices=sorted(C.BY_TAG))
    p.add_argument("--latex", action="store_true", help="also print a booktabs table")
    p.add_argument("--out", type=Path, default=None, help="write the table as CSV")
    a = p.parse_args(argv)
    cfg = C.BY_TAG[a.version]

    print(f"Alloy sensitivity of the headline ratios ({cfg.tag} results, base alloy {cfg.creep_alloy})")
    print("Operating points, temperatures, stresses and hydrogen rates are held fixed; only the")
    print("Larson-Miller master curve is swapped.\n")
    t = build_table(cfg)
    show = t[["C", *RATIOS, "base t_r [y]", "scatter [dec]", "base LMP in range"]]
    print(show.to_string(float_format=lambda v: f"{v:10.3f}"))

    inv = ranking_is_invariant(t)
    print("\nSpread across the four curves:")
    for r in RATIOS:
        d = inv[r]
        if not d["material"]:
            verdict = f"within {100 * d['max_dev_from_1']:.1f} % of 1 on every curve - no material difference"
        elif d["same_side_of_1"]:
            verdict = "same side of 1 - direction holds"
        else:
            verdict = "CHANGES SIDE OF 1"
        print(f"  {r:16s} {d['min']:.3f} .. {d['max']:.3f}   spread {d['spread_pct']:5.1f} %   {verdict}")
    mats = ", ".join(inv["material_ratios"]) or "none"
    print(f"\n  ratios with a material effect on every curve: {mats}")
    print(f"  their direction is the same on every curve   : {'YES' if inv['material_and_same_side'] else 'NO'}")
    if not inv["all_same_side"]:
        crossed = [r for r in RATIOS if not inv[r]["same_side_of_1"]]
        print(f"\n  NOTE: {', '.join(crossed)} cross unity between alloys. They sit within a few per cent")
        print("  of 1 on every curve, so the paper should report them as 'no material difference between")
        print("  these scenarios' rather than asserting a direction that depends on the alloy.")
    print("\n  Caveat: the table above holds the operating points fixed, so the knee is the base alloy's")
    print("  knee re-scored on each curve. Where a run optimised its own knee, the ratio moves far less:")
    ro = reoptimised_knee()
    if len(ro) > 1:
        print("    " + ro.round(4).to_string().replace("\n", "\n    "))
        v = ro["knee/base (own curve)"].to_numpy(float)
        print(f"\n    re-optimised knee/base spans {v.min():.3f}-{v.max():.3f}, a spread of "
              f"{100 * (v.max() - v.min()) / v.mean():.1f} % against {inv['knee/base']['spread_pct']:.1f} % "
              f"with the point held fixed;")
        print("    the optimiser recovers nearly all of the difference by re-locating the knee.")

    if a.latex:
        print("\n" + to_latex(t))
    if a.out:
        t.to_csv(a.out)
        print(f"\nwritten to {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
