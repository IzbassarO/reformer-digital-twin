"""Validation of the digitised manufacturer creep curves: ``python -m rdt.creep_validate``.

Three checks, all printed:

1. **100 000 h rupture strength.** The stress giving ``t_r = 1e5 h`` is computed from each fitted
   *minimum* (lower scatter band) curve as a function of temperature, and the temperature at which it
   equals the quoted 100 000 h design strength is reported and compared with a typical design tube-metal
   temperature of 900 C.
2. **Ingestion regression.** The synthetic-data recovery test of :func:`rdt.creep.ingest_nims` is re-run,
   confirming that the Larson-Miller constant and the scatter are still recovered within 10 %.
3. **Base-case extrapolation check.** The base operating point (1147 K, 12.9 MPa) is placed on each
   digitised curve and its LMP is checked to lie inside the range the data sheet actually covers.
"""

from __future__ import annotations

import argparse
from typing import Dict, List, Optional, Sequence

import numpy as np
from scipy.optimize import brentq

from rdt import creep

DESIGN_LIFE_H = 1.0e5
DESIGN_T_C = 900.0                  # typical design tube-metal temperature
BASE_T_K = 1147.0                   # base-case hot-spot tube-metal temperature
BASE_SIGMA_MPA = 12.9               # base-case hoop stress
ALLOYS = ("centralloy_g_4852", "centralloy_g_4852_micro", "centralloy_et_45_micro")

# Quoted 100 000 h rupture strengths to be reproduced by the fitted lower-scatter-band curves.
REFERENCE_100KH_MPA = {"centralloy_g_4852": 18.3, "centralloy_g_4852_micro": 21.2}

TEMPERATURES_C = (850.0, 875.0, 900.0, 925.0, 950.0, 975.0)


def _hr(title: str) -> str:
    return f"\n{title}\n" + "-" * len(title)


# ---------------------------------------------------------------------------
# Check 1
# ---------------------------------------------------------------------------
def rupture_strength_vs_temperature(key: str, t_r_h: float = DESIGN_LIFE_H,
                                    temps_C: Sequence[float] = TEMPERATURES_C) -> List[Dict[str, float]]:
    """Stress giving ``t_r_h`` hours on the minimum curve, at each temperature."""
    c = creep.load_alloy(key).minimum
    rows = []
    for T_C in temps_C:
        T = T_C + 273.15
        L = float(c.lmp(T, t_r_h))
        rows.append({"T_C": T_C, "lmp": L, "sigma_MPa": float(10 ** c.log10_stress(L)),
                     "in_range": bool(c.lmp_range[0] <= L <= c.lmp_range[1])})
    return rows


def temperature_for_100kh_strength(key: str, sigma_MPa: float, t_r_h: float = DESIGN_LIFE_H) -> Optional[float]:
    """Temperature [C] at which the minimum curve gives ``sigma_MPa`` after ``t_r_h`` hours."""
    c = creep.load_alloy(key).minimum
    L = c.lmp_of_stress(sigma_MPa)
    T = L / (c.scale * (c.C + np.log10(t_r_h)))     # invert LMP = scale T (C + log10 t_r)
    return float(T - 273.15)


def check_100kh(keys: Sequence[str] = ALLOYS) -> Dict[str, Dict[str, float]]:
    print(_hr("CHECK 1  100 000 h rupture strength from the fitted lower-scatter-band curves"))
    out: Dict[str, Dict[str, float]] = {}
    for key in keys:
        a = creep.load_alloy(key)
        print(f"\n{a.alloy}  (C = {a.C}, {a.source_file})")
        print("    T [C]      LMP    sigma(1e5 h) [MPa]   within digitised LMP range")
        for r in rupture_strength_vs_temperature(key):
            flag = "yes" if r["in_range"] else "NO (extrapolated)"
            print(f"    {r['T_C']:6.0f}   {r['lmp']:7.3f}   {r['sigma_MPa']:14.2f}   {flag:>22s}")
        ref = REFERENCE_100KH_MPA.get(key)
        if ref is None:
            print(f"    (no quoted 100 000 h strength for this alloy; nothing to match)")
            continue
        T_match = temperature_for_100kh_strength(key, ref)
        sigma_at_design = float(10 ** a.minimum.log10_stress(a.minimum.lmp(DESIGN_T_C + 273.15, DESIGN_LIFE_H)))
        out[key] = {"reference_MPa": ref, "T_match_C": T_match, "delta_vs_design_C": T_match - DESIGN_T_C,
                    "sigma_at_900C_MPa": sigma_at_design}
        print(f"    quoted 100 000 h strength {ref} MPa is matched at {T_match:.1f} C "
              f"({T_match - DESIGN_T_C:+.1f} C vs the {DESIGN_T_C:.0f} C design temperature)")
        print(f"    at {DESIGN_T_C:.0f} C the same curve gives {sigma_at_design:.2f} MPa for 100 000 h")
    if len(out) == 2:
        spread = abs(out[ALLOYS[0]]["T_match_C"] - out[ALLOYS[1]]["T_match_C"])
        print(f"\n  Both quoted values are matched within {spread:.1f} C of each other, i.e. they are a "
              f"consistent pair read at one temperature, not two independent numbers.")
    return out


# ---------------------------------------------------------------------------
# Check 2
# ---------------------------------------------------------------------------
def check_ingestion_regression(n_heats: int = 100, n_per_heat: int = 8, sigma_heat: float = 0.3,
                               tol: float = 0.10, tmp_dir: Optional[str] = None) -> Dict[str, float]:
    """Re-run the synthetic-data recovery test of :func:`rdt.creep.ingest_nims`."""
    import tempfile
    from pathlib import Path

    print(_hr("CHECK 2  ingestion regression: synthetic recovery of C and scatter"))
    ref = creep.LarsonMillerCurve.from_yaml()
    df = creep.synthetic_nims_dataset(n_heats=n_heats, n_per_heat=n_per_heat, sigma_heat=sigma_heat, noise=0.05, seed=0)
    with tempfile.TemporaryDirectory(dir=tmp_dir) as td:
        csv = Path(td) / "synthetic.csv"
        df.to_csv(csv, index=False)
        d = creep.ingest_nims(csv, out_dir=td)["SYNTH_XM"]["derived"]
    err_C = abs(d["C"] - ref.C) / ref.C
    err_s = abs(d["sigma_log10_heat"] - sigma_heat) / sigma_heat
    print(f"  generated {n_heats} heats x {n_per_heat} points from the {ref.curve_kind} curve of the "
          f"placeholder (C = {ref.C}, heat scatter {sigma_heat} decades)")
    print(f"  recovered C               = {d['C']:.3f}   (target {ref.C}, error {100 * err_C:.2f} %)   "
          f"{'PASS' if err_C < tol else 'FAIL'}")
    print(f"  recovered heat scatter    = {d['sigma_log10_heat']:.4f} decades   (target {sigma_heat}, "
          f"error {100 * err_s:.2f} %)   {'PASS' if err_s < tol else 'FAIL'}")
    print(f"  fitted degree {d['degree']}, grouped-CV MSE {d['cv_mse_log10_sigma']:.2e} in log10 sigma")
    return {"C": d["C"], "C_error": err_C, "sigma_log10_heat": d["sigma_log10_heat"], "scatter_error": err_s,
            "pass": bool(err_C < tol and err_s < tol)}


# ---------------------------------------------------------------------------
# Check 3
# ---------------------------------------------------------------------------
def check_base_case_in_range(T_K: float = BASE_T_K, sigma_MPa: float = BASE_SIGMA_MPA,
                             keys: Sequence[str] = ALLOYS) -> Dict[str, Dict[str, object]]:
    print(_hr(f"CHECK 3  base operating point ({T_K:.0f} K, {sigma_MPa} MPa) inside the digitised range"))
    print(f"  hoop stress {sigma_MPa} MPa at a hot-spot tube-metal temperature of {T_K:.1f} K "
          f"({T_K - 273.15:.1f} C)\n")
    out: Dict[str, Dict[str, object]] = {}
    for key in list(keys) + [creep.LEGACY_ALLOY]:
        a = creep.load_alloy(key)
        c = a.minimum
        L = c.lmp_of_stress(sigma_MPa)
        lo, hi = c.lmp_range
        t_r = float(c.time_to_rupture(T_K, sigma_MPa))
        inside = lo <= L <= hi
        margin = min(L - lo, hi - L)
        out[key] = {"alloy": a.alloy, "C": a.C, "lmp": L, "lmp_range": (lo, hi), "inside": inside,
                    "margin": margin, "t_r_h": t_r, "t_r_years": t_r / 8760.0,
                    "scatter_decades": float(a.scatter.decades(T_K, L))}
        print(f"  {a.alloy:26s} C = {a.C:5.2f}   LMP = {L:7.3f}   digitised {lo:.3f}-{hi:.3f}   "
              f"{'INSIDE' if inside else 'OUTSIDE - EXTRAPOLATED'} (margin {margin:+.3f})")
        print(f"  {'':26s} t_r = {t_r:.3e} h = {t_r / 8760.0:.0f} years   "
              f"scatter band {out[key]['scatter_decades']:.3f} decades")
    return out


# ---------------------------------------------------------------------------
def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--T-base", type=float, default=BASE_T_K)
    p.add_argument("--sigma-base", type=float, default=BASE_SIGMA_MPA)
    a = p.parse_args(argv)

    print("Validation of the Schmidt + Clemens Centralloy master curves")
    print(f"active alloy for the pipeline: {creep.alloy_key()}")
    r1 = check_100kh()
    r2 = check_ingestion_regression()
    r3 = check_base_case_in_range(a.T_base, a.sigma_base)

    print(_hr("SUMMARY"))
    ok1 = all(abs(v["delta_vs_design_C"]) < 60.0 for v in r1.values())
    ok3 = all(v["inside"] for v in r3.values())
    print(f"  1  100 000 h strengths matched near the design temperature : {'PASS' if ok1 else 'REVIEW'}")
    print(f"  2  ingestion recovers C and scatter within 10 %            : {'PASS' if r2['pass'] else 'FAIL'}")
    print(f"  3  base case inside every digitised LMP range              : {'PASS' if ok3 else 'FAIL'}")
    return 0 if (ok1 and r2["pass"] and ok3) else 1


if __name__ == "__main__":
    raise SystemExit(main())
