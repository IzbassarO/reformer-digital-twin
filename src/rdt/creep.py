"""Tube creep-life module: hoop stress, Larson-Miller master curve, Robinson life-fraction rule.

Data (Larson-Miller constant, master-curve rupture points) are loaded from YAML files in
``data/creep_derived/``; nothing alloy-specific is hard-coded here.

* :func:`hoop_stress` (mean-diameter, API 530 style), :func:`hoop_stress_thin`, :func:`hoop_stress_lame`.
* :class:`LarsonMillerCurve`: ``LMP = scale * T_K * (C + log10 t_r)``; master curve
  ``log10(sigma) = polynomial(LMP)`` fitted to tabulated (T, sigma, t_r) rupture points;
  :meth:`time_to_rupture`, :meth:`rupture_stress`, :meth:`lmp_of_stress`.
* :func:`robinson_damage` and :func:`remaining_life`: linear life-fraction rule over a piecewise-constant history.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import yaml

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "creep_derived"
DEFAULT_CURVE_FILE = DATA_DIR / "yeh2021_manaurite_xm.yaml"


# ---------------------------------------------------------------------------
# Hoop stress
# ---------------------------------------------------------------------------
def hoop_stress(P_bar: float, d_o: float, t_wall: float) -> float:
    """Mean-diameter hoop stress [MPa] (API 530 style): sigma = P (d_o - t) / (2 t).

    ``P_bar`` gauge pressure in bar, ``d_o`` outer diameter and ``t_wall`` wall thickness in m.
    """
    P = P_bar * 0.1  # MPa
    return P * (d_o - t_wall) / (2.0 * t_wall)


def hoop_stress_thin(P_bar: float, d_o: float, t_wall: float) -> float:
    """Thin-wall hoop stress on the outer diameter [MPa]: sigma = P d_o / (2 t) (Yeh 2021 Eq. 16)."""
    return P_bar * 0.1 * d_o / (2.0 * t_wall)


def hoop_stress_lame(P_bar: float, d_o: float, t_wall: float, at: str = "inner") -> float:
    """Lame thick-wall hoop stress [MPa] for internal pressure P and zero external pressure.

    ``at="inner"``: sigma = P (r_o^2 + r_i^2)/(r_o^2 - r_i^2) (maximum); ``at="outer"``:
    sigma = 2 P r_i^2/(r_o^2 - r_i^2); ``at="mean"``: average over the wall thickness.
    """
    P = P_bar * 0.1
    r_o = d_o / 2.0
    r_i = r_o - t_wall
    k = P * r_i**2 / (r_o**2 - r_i**2)
    if at == "inner":
        return k * (1.0 + r_o**2 / r_i**2)
    if at == "outer":
        return 2.0 * k
    if at == "mean":
        r = np.linspace(r_i, r_o, 201)
        return float(np.trapezoid(k * (1.0 + r_o**2 / r**2), r) / t_wall)
    raise ValueError("at must be 'inner', 'outer' or 'mean'")


# ---------------------------------------------------------------------------
# Larson-Miller
# ---------------------------------------------------------------------------
@dataclass
class LarsonMillerCurve:
    """Larson-Miller master curve ``log10(sigma_MPa) = poly(LMP)``, ``LMP = scale * T_K * (C + log10 t_r_h)``.

    ``points`` are rupture points (T_K, sigma_MPa, t_r_h) used for the polynomial fit; the fit
    range in LMP is stored so that extrapolation can be flagged.
    """

    C: float
    scale: float = 1.0e-3
    degree: int = 3
    points: List[Tuple[float, float, float]] = field(default_factory=list)
    coeffs: Optional[np.ndarray] = None          # highest power first (np.polyval)
    lmp_range: Tuple[float, float] = (np.nan, np.nan)
    alloy: str = ""
    source: str = ""
    curve_kind: str = ""

    # -- construction --------------------------------------------------------
    @classmethod
    def from_yaml(cls, path: Union[str, Path] = DEFAULT_CURVE_FILE, curve: Optional[str] = None,
                  degree: Optional[int] = None) -> "LarsonMillerCurve":
        """Load constant, scale and rupture points from a data file (see data/creep_derived/*.yaml)."""
        d = yaml.safe_load(Path(path).read_text())
        if "derived" in d:   # derived NIMS fit written by ingest_nims: coefficients only, no raw points
            dv = d["derived"]
            obj = cls(C=float(dv["C"]), scale=float(dv.get("scale", 1e-3)), degree=int(dv["degree"]), points=[], alloy=d.get("alloy", ""),
                      source=d.get("source", ""), curve_kind="nims_derived")
            obj.coeffs = np.array(dv["coefficients_high_to_low"], float); obj.lmp_range = (float(dv["lmp_range"][0]), float(dv["lmp_range"][1]))
            obj.sigma_log10_heat = float(dv.get("sigma_log10_heat", float("nan")))
            return obj
        lm = d["larson_miller"]
        curve = curve or lm.get("default_curve")
        cdef = lm["curves"][curve]
        pts = [(float(p["T_K"]), float(p["sigma_MPa"]), float(p["t_r_h"])) for p in cdef["points"]]
        obj = cls(C=float(lm["C"]), scale=float(lm.get("scale", 1e-3)), degree=int(degree or cdef.get("fit_degree", 3)),
                  points=pts, alloy=d.get("alloy", ""), source=d.get("source", ""), curve_kind=curve)
        obj.fit()
        return obj

    @classmethod
    def from_config(cls, config: Union[str, Path] = DATA_DIR / "creep_config.yaml") -> "LarsonMillerCurve":
        """Load the curve named by ``active:`` in ``data/creep_derived/creep_config.yaml`` (one-line switch).

        Superseded for pipeline runs by :func:`active_curve`, which follows ``RunConfig.creep_alloy``;
        kept for interactive use of the YAML curve files.
        """
        cfg = yaml.safe_load(Path(config).read_text())
        return cls.from_yaml(Path(config).parent / cfg["active"])

    def lmp(self, T_K, t_r_h):
        return self.scale * np.asarray(T_K, float) * (self.C + np.log10(np.asarray(t_r_h, float)))

    def fit(self) -> "LarsonMillerCurve":
        T = np.array([p[0] for p in self.points]); s = np.array([p[1] for p in self.points]); t = np.array([p[2] for p in self.points])
        L = self.lmp(T, t)
        self.coeffs = np.polyfit(L, np.log10(s), self.degree)
        self.lmp_range = (float(L.min()), float(L.max()))
        return self

    # -- evaluation ----------------------------------------------------------
    def log10_stress(self, lmp):
        return np.polyval(self.coeffs, np.asarray(lmp, float))

    def rupture_stress(self, T_K, t_r_h):
        """Stress [MPa] giving rupture after ``t_r_h`` hours at ``T_K``."""
        return 10.0 ** self.log10_stress(self.lmp(T_K, t_r_h))

    def lmp_of_stress(self, sigma_MPa: float) -> float:
        """Invert the master curve: LMP at which the rupture stress equals ``sigma_MPa``.

        The polynomial is inverted by root finding; the real root inside (or nearest to) the fit
        range is returned. The master curve is monotone decreasing in LMP over the fit range.
        """
        c = self.coeffs.copy(); c[-1] -= np.log10(sigma_MPa)
        roots = np.roots(c)
        real = np.real(roots[np.abs(np.imag(roots)) < 1e-8])
        if real.size == 0:
            raise ValueError("no real LMP root for this stress")
        lo, hi = self.lmp_range
        inside = real[(real >= lo - 1.0) & (real <= hi + 1.0)]
        if inside.size:
            return float(inside[np.argmin(np.abs(inside - 0.5 * (lo + hi)))])
        return float(real[np.argmin(np.minimum(np.abs(real - lo), np.abs(real - hi)))])

    def time_to_rupture(self, T_K, sigma_MPa):
        """Rupture time [h] at temperature ``T_K`` and stress ``sigma_MPa`` (scalar or arrays)."""
        T = np.asarray(T_K, float); s = np.asarray(sigma_MPa, float)
        Lm = np.vectorize(self.lmp_of_stress)(s) if s.ndim else self.lmp_of_stress(float(s))
        return 10.0 ** (np.asarray(Lm) / (self.scale * T) - self.C)

    def temperature_for_life(self, sigma_MPa: float, t_r_h: float) -> float:
        """Temperature [K] at which the rupture time equals ``t_r_h`` for stress ``sigma_MPa``."""
        return self.lmp_of_stress(sigma_MPa) / (self.scale * (self.C + np.log10(t_r_h)))

    def life_halving_dT(self, T_K: float, sigma_MPa: float) -> float:
        """Temperature increase [K] that halves the rupture time at fixed stress.

        From ``LMP = scale T (C + log10 t)`` at fixed LMP: ``dT = T log10(2) / (C + log10 t_r)``.
        """
        t = float(self.time_to_rupture(T_K, sigma_MPa))
        return T_K * np.log10(2.0) / (self.C + np.log10(t))

    def in_range(self, T_K, sigma_MPa) -> np.ndarray:
        L = np.vectorize(self.lmp_of_stress)(np.asarray(sigma_MPa, float))
        return (L >= self.lmp_range[0]) & (L <= self.lmp_range[1])


# ---------------------------------------------------------------------------
# Robinson life fraction
# ---------------------------------------------------------------------------
def robinson_damage(history: Sequence[Tuple[float, float, float]], curve: LarsonMillerCurve) -> Dict[str, object]:
    """Linear life-fraction damage over a piecewise-constant history of ``(dt_h, T_K, sigma_MPa)`` segments.

    ``D = sum_i dt_i / t_r(T_i, sigma_i)``; rupture is predicted at D = 1.
    """
    fractions = []
    for dt, T, s in history:
        tr = float(curve.time_to_rupture(T, s))
        fractions.append(dt / tr)
    D = float(np.sum(fractions))
    return {"damage": D, "fractions": fractions, "consumed_life_fraction_pct": 100.0 * D}


def remaining_life(history: Sequence[Tuple[float, float, float]], T_now: float, sigma_now: float,
                   curve: LarsonMillerCurve) -> Dict[str, float]:
    """Remaining time [h] under the current condition after the accumulated damage of ``history``."""
    D = robinson_damage(history, curve)["damage"]
    tr_now = float(curve.time_to_rupture(T_now, sigma_now))
    return {"damage_so_far": D, "t_r_current_h": tr_now, "remaining_h": max(0.0, (1.0 - D) * tr_now),
            "remaining_years": max(0.0, (1.0 - D) * tr_now) / 8760.0}


# ---------------------------------------------------------------------------
# Profile helpers (for the coupled model output)
# ---------------------------------------------------------------------------
def life_along_tube(z: np.ndarray, T_K: np.ndarray, P_bar: np.ndarray, d_o: float, t_wall: float,
                    curve: LarsonMillerCurve, stress: str = "mean_diameter") -> Dict[str, object]:
    """Local hoop stress and time-to-rupture along a tube from temperature and pressure profiles."""
    f = {"mean_diameter": hoop_stress, "thin": hoop_stress_thin,
         "lame_inner": lambda P, d, t: hoop_stress_lame(P, d, t, "inner")}[stress]
    sigma = np.array([f(float(p), d_o, t_wall) for p in P_bar])
    tr = np.array([float(curve.time_to_rupture(float(T), float(s))) for T, s in zip(T_K, sigma)])
    k = int(np.argmin(tr))
    return {"sigma_MPa": sigma, "t_r_h": tr, "log10_t_r": np.log10(tr), "t_r_min_h": float(tr[k]),
            "t_r_min_years": float(tr[k] / 8760.0), "z_min_m": float(z[k]), "z_frac_min": float(z[k] / z[-1]),
            "T_at_min_K": float(T_K[k]), "sigma_at_min_MPa": float(sigma[k]),
            "T_mean_K": float(np.trapezoid(T_K, z) / (z[-1] - z[0])), "sigma_mean_MPa": float(np.trapezoid(sigma, z) / (z[-1] - z[0]))}


# ---------------------------------------------------------------------------
# NIMS ingestion: raw transcription -> derived master curve (only derived quantities are stored)
# ---------------------------------------------------------------------------
TEMPLATE_COLUMNS = ["heat_id", "alloy", "T_C", "sigma_MPa", "t_r_h", "elongation_pct", "RoA_pct", "source_sheet", "page", "note"]


def _lm_fit_cv(T_K: np.ndarray, sigma: np.ndarray, t_r: np.ndarray, heats: np.ndarray, C_grid: np.ndarray,
               degrees=(1, 2, 3, 4), n_folds: int = 5, scale: float = 1e-3, seed: int = 0) -> Dict[str, object]:
    """Grid over the Larson-Miller constant and polynomial degree with grouped K-fold CV on log10(sigma)."""
    from sklearn.model_selection import GroupKFold
    y = np.log10(sigma); best = None
    gkf = GroupKFold(n_splits=min(n_folds, len(np.unique(heats))))
    for C in C_grid:
        L = scale * T_K * (C + np.log10(t_r))
        for deg in degrees:
            err = []
            for tr_idx, te_idx in gkf.split(L, y, heats):
                coef = np.polyfit(L[tr_idx], y[tr_idx], deg)
                err.append(np.mean((np.polyval(coef, L[te_idx]) - y[te_idx]) ** 2))
            cv = float(np.mean(err))
            if best is None or cv < best["cv_mse"]:
                best = {"C": float(C), "degree": int(deg), "cv_mse": cv}
    return best


def ingest_nims(csv_path: Union[str, Path], out_dir: Union[str, Path] = DATA_DIR, C_grid: Optional[np.ndarray] = None,
                degrees=(1, 2, 3), scale: float = 1e-3, min_points: int = 8, n_iter: int = 4) -> Dict[str, Dict[str, object]]:
    """Fit a Larson-Miller master curve per alloy from a transcription CSV (see nims_transcription_template.csv).

    C is optimised on a grid (default 15-30, step 0.25) and the polynomial degree chosen by grouped K-fold cross-
    validation (groups = heats), alternating ``n_iter`` times with per-heat offsets in log10 t_r so that the
    heat-to-heat scatter does not bias C (heat-normalised master curve). The scatter ``sigma_log10_heat`` is the
    standard deviation of the heat offsets (residuals converted from log10 sigma through the local slope of the
    master curve). Only derived quantities are written to ``<out_dir>/<alloy>_derived.yaml``; raw data are not copied.
    """
    df = pd.read_csv(csv_path)
    missing = set(TEMPLATE_COLUMNS[:5]) - set(df.columns)
    if missing:
        raise ValueError(f"transcription CSV lacks columns {sorted(missing)}")
    C_grid = np.arange(15.0, 30.01, 0.25) if C_grid is None else np.asarray(C_grid, float)
    out = {}
    for alloy, g in df.groupby("alloy"):
        g = g.dropna(subset=["T_C", "sigma_MPa", "t_r_h"])
        if len(g) < min_points:
            continue
        T = g.T_C.to_numpy(float) + 273.15; sig = g.sigma_MPa.to_numpy(float); tr = g.t_r_h.to_numpy(float); heats = g.heat_id.astype(str).to_numpy()
        # heat-normalised master-curve fit: alternate between (C, degree, polynomial) on heat-adjusted rupture
        # times and per-heat offsets in log10 t_r (random-effect style), so that heat scatter does not bias C
        offsets = pd.Series(0.0, index=np.unique(heats))
        for _ in range(n_iter):
            tr_adj = tr / 10 ** offsets.loc[heats].to_numpy()
            best = _lm_fit_cv(T, sig, tr_adj, heats, C_grid, degrees, scale=scale)
            L = scale * T * (best["C"] + np.log10(tr_adj)); coef = np.polyfit(L, np.log10(sig), best["degree"])
            # residuals in log10 t_r of the raw data: d(log10 t_r) = d(log10 sigma) / slope / (scale T)
            L_raw = scale * T * (best["C"] + np.log10(tr)); slope = np.polyval(np.polyder(coef), L_raw)
            res_logtr = -(np.log10(sig) - np.polyval(coef, L_raw)) / slope / (scale * T)   # slope < 0: positive offset = longer life
            offsets = pd.Series(res_logtr).groupby(heats).mean()
        heat_means = offsets
        derived = {"alloy": str(alloy), "source": "derived from NIMS transcription " + str(Path(csv_path).name) + "; raw data not redistributed",
                   "n_points": int(len(g)), "n_heats": int(len(heat_means)), "sheets": sorted(set(g.get("source_sheet", pd.Series(dtype=str)).dropna().astype(str))),
                   "derived": {"C": best["C"], "scale": scale, "degree": best["degree"], "coefficients_high_to_low": [float(c) for c in coef],
                               "lmp_range": [float(L.min()), float(L.max())], "cv_mse_log10_sigma": best["cv_mse"],
                               "sigma_log10_heat": float(heat_means.std(ddof=1)) if len(heat_means) > 1 else float("nan"),
                               "sigma_log10_total": float(np.std(res_logtr, ddof=1)),
                               "T_range_C": [float(g.T_C.min()), float(g.T_C.max())], "sigma_range_MPa": [float(sig.min()), float(sig.max())],
                               "t_r_range_h": [float(tr.min()), float(tr.max())]}}
        path = Path(out_dir) / f"{str(alloy).replace(' ', '_')}_derived.yaml"
        yaml.safe_dump(derived, open(path, "w"), sort_keys=False)
        derived["file"] = str(path); out[str(alloy)] = derived
    return out


def synthetic_nims_dataset(n_heats: int = 40, n_per_heat: int = 10, sigma_heat: float = 0.3, noise: float = 0.05,
                           seed: int = 0, alloy: str = "SYNTH_XM") -> pd.DataFrame:
    """Synthetic transcription generated from the Yeh placeholder curve with log-normal heat scatter (for tests)."""
    rng = np.random.default_rng(seed); curve = LarsonMillerCurve.from_yaml(); rows = []
    for h in range(n_heats):
        off = rng.normal(0.0, sigma_heat)
        for _ in range(n_per_heat):
            T_C = rng.uniform(850.0, 1050.0); logt = rng.uniform(2.0, 5.0)
            sig = float(curve.rupture_stress(T_C + 273.15, 10 ** logt))
            logt_obs = logt + off + rng.normal(0.0, noise)
            rows.append({"heat_id": f"H{h:03d}", "alloy": alloy, "T_C": T_C, "sigma_MPa": sig, "t_r_h": 10 ** logt_obs, "elongation_pct": np.nan, "RoA_pct": np.nan,
                         "source_sheet": "synthetic", "page": "", "note": "generated"})
    return pd.DataFrame(rows, columns=TEMPLATE_COLUMNS)


# ---------------------------------------------------------------------------
# Manufacturer data-sheet master curves (Schmidt + Clemens Centralloy)
# ---------------------------------------------------------------------------
# Curves digitised by rdt.creep_ingest_datasheet from the "Parametric stress rupture strength" chart
# of each data sheet; see data/creep_derived/sources.yaml. Each alloy carries its own Larson-Miller
# constant, read off its own sheet (G 4852: 18.6, G 4852 Micro: 22.9, ET 45 Micro: 19.3).
SOURCES_YAML = DATA_DIR / "sources.yaml"
LEGACY_ALLOY = "legacy_yeh_manaurite_xm"     # the v1/v2 placeholder, kept so those results reproduce
DEFAULT_ALLOY = "centralloy_g_4852"          # base alloy from v3 onwards (HP-Nb)
ALLOY_ENV = "RDT_CREEP_ALLOY"                # honoured by joblib worker processes, which re-import
LEGACY_SCATTER_DECADES = 0.3                 # the assumption the data-sheet scatter band replaces


def datasheet_C(alloy: str, sources: Union[str, Path] = SOURCES_YAML) -> float:
    """Larson-Miller constant printed on the data sheet of ``alloy`` (exact name match)."""
    d = yaml.safe_load(Path(sources).read_text())
    for s in d["sources"]:
        if s["alloy"] == alloy:
            return float(s["larson_miller_C"])
    raise KeyError(f"{alloy!r} not in {sources}; known: {[s['alloy'] for s in d['sources']]}")


@dataclass
class ScatterModel:
    """Rupture-time scatter taken from the width of the printed scatter band.

    The data sheets draw an *Average* curve and a *Lower Scatter Band* (95 % confidence). Measured
    **vertically** the two differ in stress, which is not directly usable; measured **horizontally**, at
    constant stress, they differ in Larson-Miller parameter, and that gap converts straight into rupture
    time. From ``LMP = scale T (C + log10 t_r)`` at fixed stress and temperature,

        delta_log10_tr = (LMP_avg - LMP_min) / (scale * T)

    so the gap in LMP is a temperature-independent property of the pair of curves, and the scatter in
    decades of life follows by dividing by ``scale * T``. ``delta_lmp`` is tabulated against the LMP of
    the average curve; ``assumed`` marks the legacy constant-0.3-decade placeholder.
    """

    alloy: str
    scale: float
    lmp: np.ndarray                     # LMP on the average curve
    delta_lmp: np.ndarray               # LMP_avg - LMP_min at the same stress
    stress_MPa: np.ndarray              # the constant-stress levels the gap was measured at
    assumed: bool = False
    note: str = ""

    def delta_lmp_at(self, lmp):
        """Horizontal gap [LMP units] at the given LMP, clamped to the digitised range."""
        return np.interp(np.asarray(lmp, float), self.lmp, self.delta_lmp)

    def decades(self, T_K, lmp=None):
        """Scatter in decades of rupture time at temperature ``T_K`` (and optionally a given LMP)."""
        g = self.delta_lmp.mean() if lmp is None else self.delta_lmp_at(lmp)
        return g / (self.scale * np.asarray(T_K, float))

    def summary(self, T_K: float) -> Dict[str, float]:
        """Range of the scatter in decades over the digitised LMP range at temperature ``T_K``."""
        d = self.delta_lmp / (self.scale * T_K)
        return {"T_K": float(T_K), "min_decades": float(d.min()), "median_decades": float(np.median(d)),
                "max_decades": float(d.max()), "mean_decades": float(d.mean()),
                "at_lmp_min": float(d[0]), "at_lmp_max": float(d[-1]),
                "lmp_range": (float(self.lmp[0]), float(self.lmp[-1])),
                "delta_lmp_min": float(self.delta_lmp.min()), "delta_lmp_max": float(self.delta_lmp.max()),
                "assumed": self.assumed}

    def describe(self, T_K: float, n: int = 6) -> str:
        s = self.summary(T_K)
        head = (f"scatter from the data-sheet band, {self.alloy}: {s['min_decades']:.3f}-{s['max_decades']:.3f} "
                f"decades of t_r at {T_K:.0f} K (median {s['median_decades']:.3f})")
        rows = [f"    LMP {l:6.2f}   dLMP {g:6.4f}   {g / (self.scale * T_K):6.3f} decades"
                for l, g in zip(np.linspace(self.lmp[0], self.lmp[-1], n),
                                self.delta_lmp_at(np.linspace(self.lmp[0], self.lmp[-1], n)))]
        return "\n".join([head, *rows])


@dataclass
class AlloyCurves:
    """The pair of master curves of one alloy plus the scatter model derived from their separation."""

    alloy: str
    C: float
    scale: float
    average: LarsonMillerCurve
    minimum: LarsonMillerCurve
    scatter: ScatterModel
    source_file: str = ""
    method: str = ""
    extracted_on: str = ""

    @property
    def lmp_range(self) -> Tuple[float, float]:
        """LMP range covered by both digitised curves (outside it the master curve is extrapolated)."""
        return (max(self.average.lmp_range[0], self.minimum.lmp_range[0]),
                min(self.average.lmp_range[1], self.minimum.lmp_range[1]))

    def curve(self, kind: str = "minimum") -> LarsonMillerCurve:
        if kind not in ("average", "minimum"):
            raise ValueError("kind must be 'average' or 'minimum'")
        return self.average if kind == "average" else self.minimum


def _fit_from_lmp(lmp: np.ndarray, sigma: np.ndarray, C: float, scale: float, degree: int,
                  alloy: str, source: str, kind: str) -> LarsonMillerCurve:
    """Fit ``log10(sigma) = poly(LMP)`` directly to a digitised curve (no (T, sigma, t_r) triples)."""
    obj = LarsonMillerCurve(C=C, scale=scale, degree=degree, points=[], alloy=alloy, source=source, curve_kind=kind)
    obj.coeffs = np.polyfit(lmp, np.log10(sigma), degree)
    obj.lmp_range = (float(lmp.min()), float(lmp.max()))
    return obj


def alloy_curves_from_csv(path: Union[str, Path], C: Optional[float] = None, degree: int = 3,
                          n_scatter: int = 200, sources: Union[str, Path] = SOURCES_YAML) -> AlloyCurves:
    """Build both master curves and the scatter model from a ``*_rupture_curve.csv`` written by
    :mod:`rdt.creep_ingest_datasheet`.

    ``log10(sigma) = P3(LMP)`` is fitted separately to the average and to the minimum (lower scatter
    band) curve, both with the alloy's own Larson-Miller constant taken from the data sheet.
    """
    df = pd.read_csv(path)
    alloy = str(df.alloy.iloc[0])
    C = float(datasheet_C(alloy, sources) if C is None else C)
    scale = 1e-3
    src = f"{df.source_file.iloc[0]} p.{df.page.iloc[0]} ({df.method.iloc[0]}, {df.extracted_on.iloc[0]})"
    fits = {}
    for kind in ("average", "minimum"):
        g = df[df.curve == kind].sort_values("lmp")
        if g.empty:
            raise ValueError(f"{path}: no {kind!r} curve")
        fits[kind] = _fit_from_lmp(g.lmp.to_numpy(float), g.stress_mpa.to_numpy(float), C, scale, degree,
                                   alloy, src, kind)

    # Horizontal gap at constant stress, over the stress range the two curves share.
    a, m = fits["average"], fits["minimum"]
    lo = max(10 ** a.log10_stress(a.lmp_range[1]), 10 ** m.log10_stress(m.lmp_range[1]))
    hi = min(10 ** a.log10_stress(a.lmp_range[0]), 10 ** m.log10_stress(m.lmp_range[0]))
    sigma = np.logspace(np.log10(lo), np.log10(hi), n_scatter)
    lmp_a = np.array([a.lmp_of_stress(float(s)) for s in sigma])
    lmp_m = np.array([m.lmp_of_stress(float(s)) for s in sigma])
    order = np.argsort(lmp_a)
    scatter = ScatterModel(alloy=alloy, scale=scale, lmp=lmp_a[order], delta_lmp=(lmp_a - lmp_m)[order],
                           stress_MPa=sigma[order], assumed=False,
                           note=f"horizontal separation of the Average and Lower Scatter Band curves of {src}")
    return AlloyCurves(alloy=alloy, C=C, scale=scale, average=a, minimum=m, scatter=scatter,
                       source_file=str(df.source_file.iloc[0]), method=str(df.method.iloc[0]),
                       extracted_on=str(df.extracted_on.iloc[0]))


def _legacy_alloy_curves() -> AlloyCurves:
    """The Yeh (2021) Manaurite XM placeholder, with its assumed constant 0.3-decade scatter.

    Kept selectable so that the v1 and v2 result files remain reproducible.
    """
    a = LarsonMillerCurve.from_yaml(DEFAULT_CURVE_FILE, curve="average")
    m = LarsonMillerCurve.from_yaml(DEFAULT_CURVE_FILE, curve="minimum")
    lmp = np.linspace(*m.lmp_range, 200)
    scatter = ScatterModel(alloy=a.alloy, scale=m.scale, lmp=lmp,
                           delta_lmp=np.full_like(lmp, LEGACY_SCATTER_DECADES * m.scale * 1150.0),
                           stress_MPa=10 ** m.log10_stress(lmp), assumed=True,
                           note=f"PLACEHOLDER: assumed constant {LEGACY_SCATTER_DECADES} decades of heat-to-heat "
                                f"scatter (v1/v2), not measured; expressed at 1150 K")
    return AlloyCurves(alloy=a.alloy, C=m.C, scale=m.scale, average=a, minimum=m, scatter=scatter,
                       source_file=DEFAULT_CURVE_FILE.name, method="figure digitisation", extracted_on="")


def available_alloys() -> Dict[str, Path]:
    """Selectable alloy keys mapped to the file they load from."""
    out = {LEGACY_ALLOY: DEFAULT_CURVE_FILE}
    for p in sorted(DATA_DIR.glob("*_rupture_curve.csv")):
        out[p.name[: -len("_rupture_curve.csv")]] = p
    return out


@lru_cache(maxsize=None)
def load_alloy(key: str) -> AlloyCurves:
    """Load an alloy by key (see :func:`available_alloys`); cached per process."""
    if key == LEGACY_ALLOY:
        return _legacy_alloy_curves()
    path = DATA_DIR / f"{key}_rupture_curve.csv"
    if not path.exists():
        raise KeyError(f"unknown alloy {key!r}; available: {sorted(available_alloys())}")
    return alloy_curves_from_csv(path)


_ACTIVE_ALLOY = DEFAULT_ALLOY


def alloy_key() -> str:
    """Key of the alloy currently selected, honouring ``RDT_CREEP_ALLOY``.

    The environment variable is what carries the selection into joblib worker processes, which
    re-import this module and would otherwise fall back to :data:`DEFAULT_ALLOY`.
    """
    return os.environ.get(ALLOY_ENV) or _ACTIVE_ALLOY


def use_config(cfg) -> None:
    """Point the module at a :class:`rdt.config.RunConfig` (selects the alloy master curve)."""
    global _ACTIVE_ALLOY
    _ACTIVE_ALLOY = str(getattr(cfg, "creep_alloy", DEFAULT_ALLOY))
    os.environ[ALLOY_ENV] = _ACTIVE_ALLOY


def active_alloy() -> AlloyCurves:
    """Both curves plus the scatter model of the currently selected alloy."""
    return load_alloy(alloy_key())


def active_curve(kind: str = "minimum") -> LarsonMillerCurve:
    """The master curve used for life calculations (the lower scatter band by default)."""
    return active_alloy().curve(kind)


def curve_for(key: str, kind: str = "minimum") -> LarsonMillerCurve:
    """Master curve of a named alloy, for code that must not rely on the ambient selection.

    Parallel workers take this rather than :func:`active_curve`: joblib reuses an existing process pool,
    and a reused worker keeps the environment it was spawned with, so a selection made after the pool
    started would not reach it. Callers capture :func:`alloy_key` in the parent and pass it down.
    """
    return load_alloy(key).curve(kind)
