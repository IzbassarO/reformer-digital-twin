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

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
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
        lm = d["larson_miller"]
        curve = curve or lm.get("default_curve")
        cdef = lm["curves"][curve]
        pts = [(float(p["T_K"]), float(p["sigma_MPa"]), float(p["t_r_h"])) for p in cdef["points"]]
        obj = cls(C=float(lm["C"]), scale=float(lm.get("scale", 1e-3)), degree=int(degree or cdef.get("fit_degree", 3)),
                  points=pts, alloy=d.get("alloy", ""), source=d.get("source", ""), curve_kind=curve)
        obj.fit()
        return obj

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
