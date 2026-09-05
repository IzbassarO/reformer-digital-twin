"""Calibration of the furnace parameters on the Latham (2008) plant cases with cross-validation.

Fitted: F_gt, L_q, alpha_top, f_htg (bounds in :data:`BOUNDS`). Fixed: activity 0.20, voidage 0.607, eta
profile, f_ctube 0.38, wall 15 mm, well-mixed top zone. Objective: sum over cases and measured quantities of
((pred - meas)/s)^2 with the measurement uncertainties :data:`SIGMA` (tube-wall temperatures use the per-case
thesis standard deviation from the CSV, fallback 20 K). ``scipy.optimize.least_squares`` (trf, bounds),
multi-start from Latin-hypercube points; leave-one-case-out cross-validation; covariance from the Jacobian.
"""

from __future__ import annotations

from dataclasses import asdict, replace
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import least_squares
from scipy.stats import qmc

from rdt import latham_cases as lc

PARAM_NAMES = ("F_gt", "L_q", "alpha_top", "f_htg")
BOUNDS = {"F_gt": (0.2, 0.9), "L_q": (0.29, 0.60), "alpha_top": (0.05, 0.40), "f_htg": (0.8, 3.0)}
LOWER = np.array([BOUNDS[k][0] for k in PARAM_NAMES]); UPPER = np.array([BOUNDS[k][1] for k in PARAM_NAMES])
CAL_CASES = ("Plant_A", "Plant_B", "Plant_C2")   # C1 excluded (thesis table inconsistency, see README)
QUANTITIES = ("out_T_K", "out_P_Pa", "out_x_H2_wet_molpct", "out_x_CH4_wet_molpct", "out_x_CO_wet_molpct",
              "out_x_CO2_wet_molpct", "flue_out_T_K", "TWT_upper_K", "TWT_lower_K")
SIGMA = {"out_T_K": 3.0, "out_P_Pa": 0.1e5, "out_x_H2_wet_molpct": 0.3, "out_x_CH4_wet_molpct": 0.3,
         "out_x_CO_wet_molpct": 0.3, "out_x_CO2_wet_molpct": 0.3, "flue_out_T_K": 10.0}
TWT_SIGMA_FALLBACK = 20.0
BASE = replace(lc.BASELINE, inlet_mode="well_mixed_top")


def params_from_vector(x: Sequence[float], base: lc.LathamParams = BASE) -> lc.LathamParams:
    return replace(base, **{k: float(v) for k, v in zip(PARAM_NAMES, x)})


def case_rows(cases: Sequence[str] = CAL_CASES) -> List[pd.Series]:
    df = lc.load_cases()
    return [df[df.case == c].iloc[0] for c in cases]


def sigma_for(row: pd.Series, q: str) -> float:
    if q in ("TWT_upper_K", "TWT_lower_K"):
        s = float(row.get("sigma_TWT_K", np.nan))
        return s if np.isfinite(s) and s > 0 else TWT_SIGMA_FALLBACK
    return SIGMA[q]


def residuals(x: Sequence[float], rows: Sequence[pd.Series], base: lc.LathamParams = BASE,
              return_details: bool = False):
    p = params_from_vector(x, base)
    res_all = []; details = {}
    for row in rows:
        res = lc.run_case(row, p)
        c = lc.compare(row, res)
        r = [(c[q]["predicted"] - c[q]["measured"]) / sigma_for(row, q) for q in QUANTITIES]
        res_all.extend(r)
        details[str(row.case)] = {"comparison": c, "result": res, "scaled_residuals": dict(zip(QUANTITIES, r))}
    r = np.array(res_all)
    return (r, details) if return_details else r


def fit(rows: Sequence[pd.Series], x0: Sequence[float], base: lc.LathamParams = BASE) -> Dict[str, object]:
    n = [0]

    def fun(x):
        n[0] += 1
        return residuals(x, rows, base)

    sol = least_squares(fun, x0=np.asarray(x0, float), bounds=(LOWER, UPPER), method="trf",
                        x_scale=UPPER - LOWER, diff_step=2e-3, xtol=1e-6, ftol=1e-8, gtol=1e-8, max_nfev=400)
    return {"x": sol.x, "params": dict(zip(PARAM_NAMES, sol.x.tolist())), "cost": float(sol.cost),
            "chi2": float(2 * sol.cost), "n_obs": int(sol.fun.size), "jac": sol.jac, "fun": sol.fun,
            "status": int(sol.status), "message": sol.message, "n_evals": n[0], "x0": list(map(float, x0)),
            "at_bound": {k: bool(np.isclose(v, BOUNDS[k][0]) or np.isclose(v, BOUNDS[k][1])) for k, v in zip(PARAM_NAMES, sol.x)}}


def lhs_starts(n: int = 5, seed: int = 7) -> np.ndarray:
    sampler = qmc.LatinHypercube(d=len(PARAM_NAMES), seed=seed)
    return qmc.scale(sampler.random(n), LOWER, UPPER)


def multistart(rows: Sequence[pd.Series], n_starts: int = 5, seed: int = 7, base: lc.LathamParams = BASE,
               extra_starts: Optional[Sequence[Sequence[float]]] = None) -> Dict[str, object]:
    starts = list(lhs_starts(n_starts, seed))
    if extra_starts:
        starts = [np.asarray(s, float) for s in extra_starts] + starts
    fits = [fit(rows, x0, base) for x0 in starts]
    fits.sort(key=lambda f: f["cost"])
    # distinct local optima: parameter vectors differing by more than 2 % of the bound range
    optima = []
    for f in fits:
        if not any(np.all(np.abs(f["x"] - g["x"]) / (UPPER - LOWER) < 0.02) for g in optima):
            optima.append(f)
    return {"fits": fits, "best": fits[0], "distinct_optima": optima}


def covariance(f: Dict[str, object]) -> Dict[str, object]:
    """Covariance from the scaled-residual Jacobian: cov = s2 (J^T J)^-1, s2 = chi2/(n - p); 95 % t-intervals."""
    J = f["jac"]; n, p = J.shape
    dof = max(n - p, 1)
    s2 = f["chi2"] / dof
    JtJ = J.T @ J
    cov = s2 * np.linalg.pinv(JtJ)
    se = np.sqrt(np.diag(cov))
    t = stats.t.ppf(0.975, dof)
    corr = cov / np.outer(se, se)
    flags = [(PARAM_NAMES[i], PARAM_NAMES[j], float(corr[i, j])) for i in range(p) for j in range(i + 1, p) if abs(corr[i, j]) > 0.9]
    return {"cov": cov, "se": dict(zip(PARAM_NAMES, se.tolist())), "t_975": float(t), "dof": int(dof),
            "reduced_chi2": float(s2),
            "ci95": {k: [float(v - t * s), float(v + t * s)] for k, v, s in zip(PARAM_NAMES, f["x"], se)},
            "corr": corr, "corr_table": {a: {b: float(corr[i, j]) for j, b in enumerate(PARAM_NAMES)} for i, a in enumerate(PARAM_NAMES)},
            "high_correlation_pairs": flags, "condition_number_JtJ": float(np.linalg.cond(JtJ))}


def errors_table(x: Sequence[float], rows: Sequence[pd.Series], base: lc.LathamParams = BASE) -> Dict[str, Dict[str, float]]:
    """Per-case prediction errors (pred - meas) for all quantities at parameters x."""
    _, det = residuals(x, rows, base, return_details=True)
    return {c: {q: d["comparison"][q]["difference"] for q in QUANTITIES} | {"TWT_peak_frac": d["comparison"]["TWT_peak"]["z_frac"],
            "T_fg_0plus_K": d["result"].extras["T_fg_0plus_K"], "T_fg_max_K": d["comparison"]["T_fg_max_K"]["value"]}
            for c, d in det.items()}


def leave_one_out(rows: Sequence[pd.Series], x0: Sequence[float], base: lc.LathamParams = BASE) -> Dict[str, object]:
    folds = {}
    for k, held in enumerate(rows):
        train = [r for j, r in enumerate(rows) if j != k]
        f = fit(train, x0, base)
        err = errors_table(f["x"], [held], base)[str(held.case)]
        folds[str(held.case)] = {"trained_on": [str(r.case) for r in train], "params": f["params"], "cost_train": f["cost"],
                                 "held_out_errors": err, "at_bound": f["at_bound"]}
    return folds


def sensitivity_screen(row: pd.Series, p_fit: lc.LathamParams) -> pd.DataFrame:
    """One-at-a-time screening for one case: +/-20 % on the fitted parameters, wall 12/15/18 mm, activity 0.15/0.20/0.25."""
    def outputs(p):
        res = lc.run_case(row, p); c = lc.compare(row, res)
        dry = res.tube.dry_mole_fractions
        return {"TWT_upper_K": c["TWT_upper_K"]["predicted"], "TWT_lower_K": c["TWT_lower_K"]["predicted"],
                "out_T_K": c["out_T_K"]["predicted"], "dry_H2_pct": 100 * float(dry["H2"][-1])}
    ref = outputs(p_fit)
    rows_out = [{"perturbation": "reference", **ref}]
    for k in PARAM_NAMES:
        for fac in (0.8, 1.2):
            v = getattr(p_fit, k) * fac
            v = min(max(v, BOUNDS[k][0]), BOUNDS[k][1])
            o = outputs(replace(p_fit, **{k: v}))
            rows_out.append({"perturbation": f"{k} x{fac:.1f} (= {v:.3g})", **{q: o[q] - ref[q] for q in ref}})
    for w in (0.012, 0.015, 0.018):
        o = outputs(replace(p_fit, wall_thickness_m=w))
        rows_out.append({"perturbation": f"wall {1e3*w:.0f} mm", **{q: o[q] - ref[q] for q in ref}})
    for a in (0.15, 0.20, 0.25):
        o = outputs(replace(p_fit, activity=a))
        rows_out.append({"perturbation": f"activity {a:.2f}", **{q: o[q] - ref[q] for q in ref}})
    return pd.DataFrame(rows_out).set_index("perturbation")
