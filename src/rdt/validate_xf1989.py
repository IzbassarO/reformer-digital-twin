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


# ---------------------------------------------------------------------------
# Heat-transfer-correlation comparison (Xu & Froment Eqs. 11-12 vs Leva-Grummer)
# ---------------------------------------------------------------------------
#: Digitised T_gas starts at ~764 K at z = 0 although Table 2 gives T0 = 793.15 K, so the
#: first 0.3 m of the temperature curves are excluded from the comparison metrics.
Z_MIN_TEMPERATURE = 0.3
TEMPERATURE_CURVES = ("T_gas", "T_wall_inner", "T_wall_outer")
STAGE_A_BED = {"voidage": 0.526, "d_p_m": 9.44e-3}  # from the stage-A pressure fit (rounded)


def run_case(heat_transfer: str, inlet: str, wall: r1.WallBC, voidage: float = STAGE_A_BED["voidage"],
             d_p: float = STAGE_A_BED["d_p_m"], eta=None, L: float = L_TOTAL, n_out: int = 241,
             f_htg: float = 1.0) -> r1.Result:
    """Xu & Froment case with the natural gas split into CH4/C2+ and the chosen inlet rule and
    heat-transfer correlation; eta = 0.1 unless given; adiabatic beyond 11.12 m."""
    tube = r1.tube_from_xu_froment()
    bed0 = r1.bed_from_xu_froment(voidage=voidage, eta=eta if eta is not None else (ETA_BASE,) * 3)
    bed = r1.CatalystBed(rho_bed=bed0.rho_bed, voidage=voidage, d_p=d_p, eta=bed0.eta, activity=bed0.activity)
    feed = r1.feed_from_xu_froment(split_alkanes=True, inlet_higher_alkanes=inlet)
    return r1.simulate(tube, bed, feed, wall, L=L, n_out=n_out, f_htg=f_htg,
                       adiabatic_beyond_heated=True, heat_transfer=heat_transfer)


def metrics_zmin(res: r1.Result, cv, which: Sequence[str] = COMPARED,
                 z_min_temperature: float = Z_MIN_TEMPERATURE) -> Dict[str, Dict[str, float]]:
    """Like :func:`metrics` but temperature curves exclude z < ``z_min_temperature``."""
    out = {}
    for c in which:
        z, v = cv[c]
        if c in TEMPERATURE_CURVES:
            m = z >= z_min_temperature
            z, v = z[m], v[m]
        err = model_curve(res, c, z) - v
        out[c] = {"rmse": float(np.sqrt(np.mean(err**2))), "max_abs": float(np.max(np.abs(err))),
                  "mean_bias": float(np.mean(err)), "n": int(len(z)), "unit": UNITS[c],
                  "z_min": float(z[0])}
    return out


def inlet_film_dT(res: r1.Result) -> float:
    """Inner-wall minus gas temperature at z = 0 [K]."""
    return float(res.T_wall_inner[0] - res.T[0])


def four_combinations(cv, wall: r1.WallBC, **kw) -> Dict[str, Dict[str, object]]:
    """Run {leva_grummer, xu_froment} x {latham, xu_froment} and collect metrics and outlets."""
    out = {}
    for ht in ("leva_grummer", "xu_froment"):
        for inlet in ("latham", "xu_froment"):
            res = run_case(ht, inlet, wall, **kw)
            out[f"{ht}__inlet_{inlet}"] = {
                "heat_transfer": ht, "inlet_higher_alkanes": inlet,
                "metrics": metrics_zmin(res, cv), "outlet": outlet_summary(res),
                "inlet_conversion": float(res.conversion_CH4[0]),
                "inlet_film_dT_K": inlet_film_dT(res),
                "alpha_i_W_m2K": {"inlet": float(res.alpha_i[0]),
                                  "z_5m": float(np.interp(5.0, res.z, res.alpha_i)),
                                  "z_11m": float(np.interp(11.0, res.z, res.alpha_i))},
                "U_W_m2K": {"inlet": float(res.U[0]), "z_11m": float(np.interp(11.0, res.z, res.U))},
                "result": res,
            }
    return out


def alpha_i_from_digitised(cv, lambda_tube: float = r1.LAMBDA_TUBE_DEFAULT,
                           d_i: float = 0.1016, d_o: float = 0.1322,
                           z_grid: Sequence[float] = (0.5, 2.0, 4.0, 6.0, 8.0, 10.0, 11.0)) -> Dict[str, list]:
    """Back-calculate the bed-side coefficient implied by the digitised Fig. 3 curves.

    ``q_o = (T_wo - T_wi) 2 lambda_tube / (d_o ln(d_o/d_i))`` (wall conduction),
    ``q_i = q_o d_o / d_i``, ``alpha_i = q_i / (T_wi - T_gas)``. Depends on the assumed tube
    conductivity (not stated by Xu & Froment); ``lambda_tube`` default is Latham's 29.6 W/(m K).
    """
    zo, To = cv["T_wall_outer"]; zi, Ti = cv["T_wall_inner"]; zg, Tg = cv["T_gas"]
    rows = {"z_m": [], "T_wo": [], "T_wi": [], "T_gas": [], "q_o_kW_m2": [], "alpha_i_W_m2K": []}
    for z in z_grid:
        two, twi, tg = np.interp(z, zo, To), np.interp(z, zi, Ti), np.interp(z, zg, Tg)
        q_o = (two - twi) * 2.0 * lambda_tube / (d_o * np.log(d_o / d_i))
        q_i = q_o * d_o / d_i
        rows["z_m"].append(float(z)); rows["T_wo"].append(float(two)); rows["T_wi"].append(float(twi))
        rows["T_gas"].append(float(tg)); rows["q_o_kW_m2"].append(float(q_o / 1e3))
        rows["alpha_i_W_m2K"].append(float(q_i / (twi - tg)))
    return rows


# ---------------------------------------------------------------------------
# Verification with a matched (back-calculated) film coefficient
# ---------------------------------------------------------------------------
Z_INCREMENT_REF = 0.5   # reference point for the conversion increment x(z) - x(z_ref)


def alpha_i_profile_from_digitised(cv, lambda_tube: float = r1.LAMBDA_TUBE_DEFAULT,
                                   d_i: float = 0.1016, d_o: float = 0.1322,
                                   z_min: float = 0.5, z_max: float = 11.0, n: int = 106) -> Dict[str, object]:
    """Back-calculated bed-side coefficient alpha_i(z) on a dense grid over [z_min, z_max].

    Same method as :func:`alpha_i_from_digitised`: wall conduction from
    ``T_wall_outer - T_wall_inner`` with ``lambda_tube`` (29.6 W/(m K) assumed, Latham's value),
    then ``alpha_i = q_i / (T_wall_inner - T_gas)``. Returns the grid, the profile, its median and
    interquartile range, and a PCHIP :class:`~rdt.reactor1d.WallBC`-like callable held constant
    outside the grid.
    """
    from scipy.interpolate import PchipInterpolator

    zo, To = cv["T_wall_outer"]; zi, Ti = cv["T_wall_inner"]; zg, Tg = cv["T_gas"]
    z = np.linspace(z_min, z_max, n)
    two, twi, tg = np.interp(z, zo, To), np.interp(z, zi, Ti), np.interp(z, zg, Tg)
    q_o = (two - twi) * 2.0 * lambda_tube / (d_o * np.log(d_o / d_i))
    alpha = q_o * d_o / d_i / (twi - tg)
    q25, q50, q75 = np.percentile(alpha, [25, 50, 75])
    f = PchipInterpolator(z, alpha, extrapolate=False)
    lo, hi = z[0], z[-1]

    def profile(zz):
        return float(f(np.clip(zz, lo, hi)))

    return {"z": z, "alpha_i": alpha, "q_o_W_m2": q_o, "median": float(q50), "q25": float(q25),
            "q75": float(q75), "iqr": float(q75 - q25), "min": float(alpha.min()), "max": float(alpha.max()),
            "lambda_tube": lambda_tube, "profile": profile}


def run_matched(alpha_i, wall: r1.WallBC, latham_inlet: bool = False, inlet: str = "latham",
                voidage: float = STAGE_A_BED["voidage"], d_p: float = STAGE_A_BED["d_p_m"],
                L: float = L_TOTAL, n_out: int = 241) -> r1.Result:
    """Xu & Froment case with a prescribed bed-side coefficient (number or callable of z)."""
    tube = r1.tube_from_xu_froment()
    eta = eta_spec(1.0, latham=latham_inlet, L=L)
    bed0 = r1.bed_from_xu_froment(voidage=voidage, eta=eta)
    bed = r1.CatalystBed(rho_bed=bed0.rho_bed, voidage=voidage, d_p=d_p, eta=bed0.eta, activity=bed0.activity)
    feed = r1.feed_from_xu_froment(split_alkanes=True, inlet_higher_alkanes=inlet)
    return r1.simulate(tube, bed, feed, wall, L=L, n_out=n_out, adiabatic_beyond_heated=True,
                       heat_transfer="constant", alpha_i_const=alpha_i)


def increment_metric(res: r1.Result, cv, z_ref: float = Z_INCREMENT_REF) -> Dict[str, float]:
    """RMSE / bias of the conversion increment x_CH4(z) - x_CH4(z_ref) for z >= z_ref.

    Removes the dependence on how conversion is defined at z = 0 (inlet-definition ambiguity).
    """
    z, v = cv["x_CH4"]
    m = z >= z_ref
    data_inc = v[m] - np.interp(z_ref, z, v)
    model_inc = model_curve(res, "x_CH4", z[m]) - np.interp(z_ref, res.z, res.conversion_CH4)
    err = model_inc - data_inc
    return {"rmse": float(np.sqrt(np.mean(err**2))), "max_abs": float(np.max(np.abs(err))),
            "mean_bias": float(np.mean(err)), "n": int(m.sum()), "z_ref": z_ref, "unit": "-"}


def verification_runs(cv, wall: r1.WallBC) -> Dict[str, object]:
    """(a) constant alpha_i = median, (b) alpha_i(z) PCHIP profile; each with eta = 0.1 and with
    Latham's inlet eta profile. Returns per-run metrics, outlets and the recommended configuration."""
    prof = alpha_i_profile_from_digitised(cv)
    runs = {}
    for a_name, a_val in (("const_median", prof["median"]), ("profile_pchip", prof["profile"])):
        for e_name, lat in (("eta_0.1", False), ("eta_latham_inlet", True)):
            res = run_matched(a_val, wall, latham_inlet=lat)
            key = f"alpha_{a_name}__{e_name}"
            runs[key] = {"alpha_i": a_name, "eta": e_name, "metrics": metrics_zmin(res, cv),
                         "x_CH4_increment": increment_metric(res, cv), "outlet": outlet_summary(res),
                         "inlet_film_dT_K": inlet_film_dT(res), "result": res}
    # recommendation: smallest RMSE(T_gas) + 1000*RMSE(increment) among runs (equal weighting, K vs 1e-3 conversion)
    score = {k: d["metrics"]["T_gas"]["rmse"] + 1000.0 * d["x_CH4_increment"]["rmse"] for k, d in runs.items()}
    best = min(score, key=score.get)
    return {"alpha_profile": prof, "runs": runs, "recommended": best, "score": score}
