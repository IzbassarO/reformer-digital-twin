"""Validation of the 1-D tube model against Xu & Froment (1989, Part II) Fig. 3.

Digitised curves (WebPlotDigitizer) of Fig. 3 (p. 101) are stored in
``data/literature_validation/digitized/xf1989_fig3_tidy.csv`` with columns
``curve, z_m, value_raw, value, unit, note``. Curves: ``T_wall_outer``,
``T_wall_inner``, ``T_gas`` [K], ``x_CH4``, ``x_CO2`` [-] (left axis
``x/(1+zeta)``, corrected by ``1+zeta`` with ``zeta = CO2/CH4 = 0.056``) and
``p_t`` [bar].

The digitised outer-wall temperature over the heated length (0-11.12 m) is the
boundary condition; the remaining 11.12-12 m are treated as adiabatic, as Fig. 3
shows the wall temperatures collapsing onto the gas temperature there.

Calibration (two stages, bounded least squares):

* stage A: bed voidage and equivalent particle diameter fitted to ``p_t`` only;
* stage B: with stage-A values fixed, one multiplier ``m`` on all three
  effectiveness factors (``eta = 0.1 m``) fitted to ``x_CH4`` and ``T_gas``.

Objective: residuals ``(model - data)`` on the digitised z grid of each curve,
divided by the data range of that curve and by ``sqrt(n_curve)`` so that curves
with different units and numbers of points are weighted equally; the sum of
squares of these normalised residuals is minimised (``scipy.optimize.least_squares``,
``trf`` with bounds).
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Callable, Dict, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from rdt import reactor1d as r1

ROOT = Path(__file__).resolve().parents[2]
TIDY_PATH = ROOT / "data" / "literature_validation" / "digitized" / "xf1989_fig3_tidy.csv"
FIT_PATH = ROOT / "data" / "literature_validation" / "xf1989_fit.yaml"

CURVES = ("T_wall_outer", "T_wall_inner", "T_gas", "x_CH4", "x_CO2", "p_t")
COMPARED = ("x_CH4", "x_CO2", "T_gas", "T_wall_inner", "p_t")
UNITS = {"T_wall_outer": "K", "T_wall_inner": "K", "T_gas": "K", "x_CH4": "-", "x_CO2": "-", "p_t": "bar"}
ZETA = 0.056
L_TOTAL = 12.0
ETA_BASE = 0.1

#: Sanity ranges for the digitised data (Xu & Froment Table 2 / Fig. 3 reading).
SANITY = {
    "T_wall_outer_z0": (1000.0, 1040.0),
    "T_wall_outer_max": (1160.0, 1190.0),
    "p_t_z0": (28.5, 29.5),
    "p_t_z11.5": (25.0, 26.0),
    "x_CH4_z11.12": (0.55, 0.65),
}


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
def load_tidy(path: Path = TIDY_PATH) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = set(CURVES) - set(df.curve.unique())
    if missing:
        raise ValueError(f"tidy CSV lacks curves {sorted(missing)}")
    return df


def curves(tidy: pd.DataFrame) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    out = {}
    for c in CURVES:
        d = tidy[tidy.curve == c].sort_values("z_m")
        out[c] = (d.z_m.to_numpy(float), d.value.to_numpy(float))
    return out


def value_near(z: np.ndarray, v: np.ndarray, z0: float) -> float:
    return float(v[np.argmin(np.abs(z - z0))])


def sanity_checks(tidy: pd.DataFrame) -> Dict[str, Tuple[float, Tuple[float, float], bool]]:
    """Return {name: (value, (lo, hi), ok)} for the five sanity checks."""
    cv = curves(tidy)
    zo, To = cv["T_wall_outer"]
    zp, p = cv["p_t"]
    zx, x = cv["x_CH4"]
    vals = {
        "T_wall_outer_z0": value_near(zo, To, 0.0),
        "T_wall_outer_max": float(To.max()),
        "p_t_z0": value_near(zp, p, 0.0),
        "p_t_z11.5": value_near(zp, p, 11.5),
        "x_CH4_z11.12": value_near(zx, x, 11.12),
    }
    return {k: (v, SANITY[k], SANITY[k][0] <= v <= SANITY[k][1]) for k, v in vals.items()}


def wall_from_digitised(tidy: pd.DataFrame, L_heated: float = 11.12) -> r1.WallBC:
    """PCHIP wall BC from the digitised T_wall_outer over 0 <= z <= L_heated."""
    z, T = curves(tidy)["T_wall_outer"]
    m = z <= L_heated + 1e-9
    z, T = z[m], T[m]
    if z[0] > 0.0:  # hold the first digitised value back to z = 0
        z = np.concatenate([[0.0], z]); T = np.concatenate([[T[0]], T])
    return r1.WallBC.from_table(z, T, description="digitised Xu & Froment Fig. 3 T_wall_outer (PCHIP), adiabatic beyond 11.12 m")


# ---------------------------------------------------------------------------
# Model runs
# ---------------------------------------------------------------------------
def eta_spec(mult: float = 1.0, latham: bool = False, L: float = L_TOTAL):
    """Effectiveness factors: ``0.1*mult`` everywhere, or Latham's 0.05 in the first 10 % of ``L``."""
    if not latham:
        return (ETA_BASE * mult,) * 3
    z_cut = 0.1 * L

    def eta(z):
        e = (0.5 if z < z_cut else 1.0) * ETA_BASE * mult
        return (e, e, e)

    return eta


def run(voidage: float, d_p: float, wall: r1.WallBC, eta=None, L: float = L_TOTAL,
        n_out: int = 241, f_htg: float = 1.0) -> r1.Result:
    """Run the Xu & Froment case with the given bed parameters and wall BC (adiabatic beyond 11.12 m)."""
    tube = r1.tube_from_xu_froment()
    bed = r1.bed_from_xu_froment(voidage=voidage, eta=eta if eta is not None else (ETA_BASE,) * 3)
    bed = r1.CatalystBed(rho_bed=bed.rho_bed, voidage=voidage, d_p=d_p, eta=bed.eta, activity=bed.activity)
    feed = r1.feed_from_xu_froment()
    return r1.simulate(tube, bed, feed, wall, L=L, n_out=n_out, f_htg=f_htg, adiabatic_beyond_heated=True)


def model_curve(res: r1.Result, curve: str, z: np.ndarray) -> np.ndarray:
    src = {"x_CH4": res.conversion_CH4, "x_CO2": res.yield_CO2, "T_gas": res.T,
           "T_wall_inner": res.T_wall_inner, "T_wall_outer": res.T_wall_outer, "p_t": res.P}[curve]
    return np.interp(z, res.z, src)


def metrics(res: r1.Result, cv: Mapping[str, Tuple[np.ndarray, np.ndarray]],
            which: Sequence[str] = COMPARED) -> Dict[str, Dict[str, float]]:
    out = {}
    for c in which:
        z, v = cv[c]
        err = model_curve(res, c, z) - v
        out[c] = {"rmse": float(np.sqrt(np.mean(err**2))), "max_abs": float(np.max(np.abs(err))),
                  "mean_bias": float(np.mean(err)), "n": int(len(z)), "unit": UNITS[c]}
    return out


def outlet_summary(res: r1.Result) -> Dict[str, float]:
    o = res.outlet()
    # also the values at the end of the heated length
    k = int(np.argmin(np.abs(res.z - 11.12)))
    o.update({"conversion_CH4_at_11.12": float(res.conversion_CH4[k]), "T_K_at_11.12": float(res.T[k]),
              "P_bar_at_11.12": float(res.P[k]), "dP_bar_total": float(res.P[0] - res.P[-1]),
              "T_min_K": float(res.T.min()), "z_Tmin_m": float(res.z[int(np.argmin(res.T))]),
              "inlet_dip_K": float(res.T[0] - res.T.min())})
    return o


# ---------------------------------------------------------------------------
# Calibration
# ---------------------------------------------------------------------------
def _norm_residuals(res: r1.Result, cv, which) -> np.ndarray:
    parts = []
    for c in which:
        z, v = cv[c]
        rng = float(v.max() - v.min())
        parts.append((model_curve(res, c, z) - v) / rng / np.sqrt(len(z)))
    return np.concatenate(parts)


def calibrate_stage_a(cv, wall: r1.WallBC, x0=(0.6, 9.24e-3),
                      bounds=((0.35, 4e-3), (0.60, 12e-3))) -> Dict[str, object]:
    """Fit voidage and d_p [m] to p_t only."""
    n = [0]

    def fun(x):
        n[0] += 1
        res = run(x[0], x[1], wall)
        return _norm_residuals(res, cv, ("p_t",))

    t0 = time.perf_counter()
    sol = least_squares(fun, x0=np.array(x0), bounds=(np.array(bounds[0]), np.array(bounds[1])),
                        method="trf", x_scale=np.array([0.1, 1e-3]), diff_step=1e-3, xtol=1e-6, ftol=1e-8)
    return {"voidage": float(sol.x[0]), "d_p_m": float(sol.x[1]), "cost": float(sol.cost),
            "n_runs": n[0], "status": int(sol.status), "message": sol.message,
            "wall_time_s": time.perf_counter() - t0, "at_bound": bool(np.any(np.isclose(sol.x, bounds[0])) or np.any(np.isclose(sol.x, bounds[1])))}


def calibrate_stage_b(cv, wall: r1.WallBC, voidage: float, d_p: float, x0=1.0,
                      bounds=(0.2, 2.0), latham: bool = False) -> Dict[str, object]:
    """Fit the effectiveness multiplier m (eta = 0.1 m) to x_CH4 and T_gas jointly."""
    n = [0]

    def fun(x):
        n[0] += 1
        res = run(voidage, d_p, wall, eta=eta_spec(float(x[0]), latham=latham))
        return _norm_residuals(res, cv, ("x_CH4", "T_gas"))

    t0 = time.perf_counter()
    sol = least_squares(fun, x0=np.array([x0]), bounds=(np.array([bounds[0]]), np.array([bounds[1]])),
                        method="trf", diff_step=1e-3, xtol=1e-6, ftol=1e-8)
    return {"eta_multiplier": float(sol.x[0]), "eta_values": [ETA_BASE * float(sol.x[0])] * 3,
            "cost": float(sol.cost), "n_runs": n[0], "status": int(sol.status), "message": sol.message,
            "wall_time_s": time.perf_counter() - t0, "latham_inlet_variant": latham,
            "at_bound": bool(np.isclose(sol.x[0], bounds[0]) or np.isclose(sol.x[0], bounds[1]))}
