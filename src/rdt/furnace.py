"""One-dimensional "long-furnace" model of a top-fired, co-current box reformer, coupled to the tube model.

Plug-flow furnace gas flows downward, co-current with the process gas (Latham 2008/2011 geometry:
top-fired box, 7 rows x 48 tubes, 8 rows x 12 burners), without the Hottel zone method. Radiation
between furnace gas and tubes is lumped into one effective exchange factor ``F_gt``.

Per unit length of furnace height ``z`` [m] (0 at the roof):

    n_fg cp_fg dT_fg/dz = r(z) - q(z) pi d_o N_tubes            (furnace-gas energy balance)
    q = sigma F_gt (T_fg^4 - T_wo^4) + h_conv (T_fg - T_wo)     (gas-to-tube flux per OUTER tube area)
    q = U (T_wo - T_gas) d_i/d_o                                 (tube-side flux, inner-area U from reactor1d)

with the heat-release density ``r(z)`` [W/m] a downward-opening parabola on ``0 <= z <= L_q L`` that is
zero at ``z = L_q L``, releases the fraction ``alpha_top`` of the heat in the top ``1/n_sections``
of the height (Latham used 15 vertical sections) and integrates to ``(1 - f_loss) Q_comb`` (Latham
2008 Sec. 3.2, Eqs. 36-42, made continuous). Flue-gas properties from Cantera (``gri30.yaml``).
The outer-wall temperature is solved at every ``z`` from the flux balance by a bracketed root
(brentq) inside the ODE right-hand side. The tube physics (kinetics, bed heat transfer, wall
conduction, pressure drop) are reused from :mod:`rdt.reactor1d` through :func:`reactor1d.tube_derivatives`.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from rdt import reactor1d as r1

SIGMA_SB = 5.670374419e-8  # W/(m2 K4)

#: Lower heating values at 298.15 K (water vapour product) for alkanes absent from gri30 [J/mol]
LHV_TABULATED_J_MOL = {"C4H10": 2657.3e3, "C5H12": 3244.9e3, "C6H14": 3855.1e3}
O2_STOICH = {"CH4": 2.0, "C2H6": 3.5, "C3H8": 5.0, "C4H10": 6.5, "C5H12": 8.0, "C6H14": 9.5, "H2": 0.5, "CO": 0.5}
C_ATOMS = {"CH4": 1, "C2H6": 2, "C3H8": 3, "C4H10": 4, "C5H12": 5, "C6H14": 6, "CO": 1, "CO2": 1}
H_ATOMS = {"CH4": 4, "C2H6": 6, "C3H8": 8, "C4H10": 10, "C5H12": 12, "C6H14": 14, "H2": 2, "H2O": 2}
#: cp surrogate for species absent from gri30
CP_SURROGATE = {"C4H10": "C3H8", "C5H12": "C3H8", "C6H14": "C3H8"}


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class FurnaceGeometry:
    """Box-furnace geometry.

    N_tubes : tubes in the furnace
    L : heated height [m] (tube exposed length)
    A_flue_free : free cross-section for the downward flue-gas flow [m2] (plan area minus tubes)
    D_h : hydraulic diameter of that flow passage [m]
    """

    N_tubes: int
    L: float
    A_flue_free: float
    D_h: float
    description: str = ""


#: Latham furnace: 7 x 48 tubes of 0.146 m OD; plan dimensions are not given in the thesis. ASSUMPTIONS:
#: tube pitch 0.30 m along a row (48 tubes -> 14.4 m), row-to-row pitch 2.0 m (7 rows -> 14 m) give a plan
#: area of 201.6 m2; minus 336 tubes (5.6 m2) -> 196 m2 free area; wetted perimeter = tubes (154 m) + walls
#: (56.8 m) -> D_h = 4 A / P = 3.7 m. Only the weak convective term depends on these numbers.
LATHAM_FURNACE = FurnaceGeometry(N_tubes=336, L=12.5, A_flue_free=196.0, D_h=3.7,
                                 description="Latham 2008 top-fired box; plan dimensions assumed (0.30 m tube pitch, 2.0 m row pitch)")


@dataclass(frozen=True)
class HeatRelease:
    """Latham's two-parameter heat-release shape.

    L_q : heat-release length as a fraction of the height (Latham 2011 Table 6: 0.48)
    alpha_top : fraction of the released heat in the top 1/n_sections of the height (0.182)
    f_loss : refractory loss fraction of Q_comb (design value 0.02, Latham 2008 p. 96)
    n_sections : number of vertical sections defining the "top section" (15)
    """

    L_q: float = 0.48
    alpha_top: float = 0.182
    f_loss: float = 0.02
    n_sections: int = 15

    def coefficients(self, L: float) -> Tuple[float, float, float, float]:
        """(a, b, c, z_q) of ``r(z)/Q_eff = a z^2 + b z + c`` on ``[0, z_q]``, ``z_q = L_q L``.

        Conditions: r(z_q) = 0; integral over [0, z_q] = 1; integral over [0, L/n_sections] = alpha_top.
        Raises ``ValueError`` if the parabola does not open downward or becomes negative on [0, z_q].
        """
        zq = self.L_q * L
        h = L / self.n_sections
        A = np.array([[zq**2, zq, 1.0],
                      [zq**3 / 3.0, zq**2 / 2.0, zq],
                      [h**3 / 3.0, h**2 / 2.0, h]])
        a, b, c = np.linalg.solve(A, np.array([0.0, 1.0, self.alpha_top]))
        zz = np.linspace(0.0, zq, 401)
        if a >= 0.0:
            raise ValueError(f"heat-release parabola opens upward (a = {a:.3g}); adjust L_q/alpha_top")
        if np.min(a * zz**2 + b * zz + c) < -1e-9:
            raise ValueError("heat-release parabola becomes negative on [0, z_q]; adjust L_q/alpha_top")
        return float(a), float(b), float(c), float(zq)

    def density(self, z, Q_comb: float, L: float):
        """Heat-release density r(z) [W/m]; integrates to (1 - f_loss) Q_comb."""
        a, b, c, zq = self.coefficients(L)
        z = np.asarray(z, dtype=float)
        r = (a * z**2 + b * z + c) * (1.0 - self.f_loss) * Q_comb
        return np.where((z >= 0.0) & (z <= zq), r, 0.0)


@dataclass(frozen=True)
class FlueGas:
    """Furnace gas: molar flow [kmol/h], mole fractions, inlet (mixed) temperature [K], combustion heat [W]."""

    n_kmol_h: float
    X: Mapping[str, float]
    T_in_K: float
    Q_comb_W: float
    excess_air_pct: Optional[float] = None
    notes: str = ""

    @classmethod
    def from_composition(cls, n_kmol_h: float, X: Mapping[str, float], T_in_K: float, Q_comb_W: float,
                         excess_air_pct: Optional[float] = None) -> "FlueGas":
        """Flue gas given directly (e.g. Latham 2011 Table 8: flow, heat of combustion, excess air)."""
        tot = sum(X.values())
        return cls(n_kmol_h, {k: v / tot for k, v in X.items()}, T_in_K, Q_comb_W, excess_air_pct, "given")

    @classmethod
    def from_combustion(cls, streams: Sequence[Tuple[float, Mapping[str, float], float]]) -> "FlueGas":
        """Complete combustion in excess air of any number of streams ``(n_kmol_h, X, T_K)``.

        Products CO2, H2O, N2, O2 from element balances; ``Q_comb`` = sum of lower heating values at
        298.15 K (gri30 enthalpies of formation for CH4, C2H6, C3H8, H2, CO; tabulated for C4-C6);
        ``T_in`` = adiabatic-mixing temperature of the unburnt streams (Latham combusts isothermally
        at the mixed inlet temperature and distributes the heat over the zones).
        """
        ct, gas = _ct()
        lhv = _lhv_table()
        Q = 0.0; o2_st = 0.0; C = 0.0; H = 0.0; N2 = 0.0; O2 = 0.0
        H_tot = 0.0; n_tot = 0.0; xmix: Dict[str, float] = {}
        for n, X, T in streams:
            for sp, x in X.items():
                if x <= 0.0:
                    continue
                if sp in lhv:
                    Q += n * x * lhv[sp] * 1e3 / 3600.0     # kmol/h * J/mol * 1e3 mol/kmol / 3600 -> W
                o2_st += n * x * O2_STOICH.get(sp, 0.0)
                C += n * x * C_ATOMS.get(sp, 0.0)
                H += n * x * H_ATOMS.get(sp, 0.0)
                if sp == "N2":
                    N2 += n * x
                if sp == "O2":
                    O2 += n * x
                key = CP_SURROGATE.get(sp, sp)
                xmix[key] = xmix.get(key, 0.0) + n * x
            gas.TPX = T, ct.one_atm, {CP_SURROGATE.get(sp, sp): x for sp, x in X.items() if x > 0.0}
            H_tot += n * gas.enthalpy_mole
            n_tot += n
        if o2_st <= 0.0:
            raise ValueError("no combustible species in the furnace streams")
        prod = {"CO2": C, "H2O": H / 2.0, "N2": N2, "O2": O2 - o2_st}
        if prod["O2"] < 0.0:
            raise ValueError("sub-stoichiometric air: complete combustion not possible")
        n_fg = sum(prod.values())
        Xmix = {k: v / n_tot for k, v in xmix.items()}
        h_target = H_tot / n_tot

        def g(T):
            gas.TPX = T, ct.one_atm, Xmix
            return gas.enthalpy_mole - h_target

        T_mix = brentq(g, 200.0, 1500.0)
        return cls(n_fg, {k: v / n_fg for k, v in prod.items()}, T_mix, Q,
                   100.0 * (O2 / o2_st - 1.0), "complete combustion of given streams (LHV basis)")


@dataclass(frozen=True)
class FurnaceParams:
    """F_gt : effective gas-to-tube radiative exchange factor [-] (initial guess 0.35, to be fitted)
    f_ctube : multiplier on the furnace-side convective coefficient (Latham 2011 Table 6: 0.38)"""

    F_gt: float = 0.35
    f_ctube: float = 1.0


@dataclass
class CoupledResult:
    tube: r1.Result
    z: np.ndarray
    T_fg: np.ndarray            # furnace gas [K]
    T_wo: np.ndarray            # outer wall [K]
    T_wi: np.ndarray            # inner wall [K]
    q_outer: np.ndarray         # W/m2 (outer tube area)
    release: np.ndarray         # W/m
    h_conv: np.ndarray          # W/(m2 K)
    energy: Dict[str, float]
    success: bool
    message: str
    n_rhs_evals: int
    wall_time_s: float
    extras: Dict[str, object] = field(default_factory=dict)

    @property
    def z_frac(self) -> np.ndarray:
        return self.z / self.z[-1]

    def twt_peak(self) -> Tuple[float, float]:
        """(z_frac, T_wo) at the outer-wall temperature maximum."""
        k = int(np.argmax(self.T_wo))
        return float(self.z_frac[k]), float(self.T_wo[k])


# ---------------------------------------------------------------------------
# Cantera helpers
# ---------------------------------------------------------------------------
@lru_cache(maxsize=1)
def _ct():
    import cantera as ct

    return ct, ct.Solution("gri30.yaml")


@lru_cache(maxsize=1)
def _lhv_table() -> Dict[str, float]:
    """LHV [J/mol] at 298.15 K from gri30 enthalpies (gas-phase water) plus tabulated C4-C6."""
    ct, gas = _ct()

    def hf(sp):
        gas.TPX = 298.15, ct.one_atm, {sp: 1.0}
        return gas.enthalpy_mole / 1e3  # J/mol

    out = {}
    for sp, nC, nH in (("CH4", 1, 4), ("C2H6", 2, 6), ("C3H8", 3, 8), ("H2", 0, 2)):
        out[sp] = -(nC * hf("CO2") + nH / 2 * hf("H2O") - hf(sp))
    out["CO"] = -(hf("CO2") - hf("CO"))
    out.update(LHV_TABULATED_J_MOL)
    return out


def flue_props(T: float, X: Mapping[str, float], P: float = 101325.0) -> Dict[str, float]:
    """cp [J/(kmol K)], h [J/kmol], rho [kg/m3], mu [Pa s], lam [W/(m K)], Pr, M [kg/kmol] of the flue gas."""
    ct, gas = _ct()
    gas.TPX = float(T), P, dict(X)
    return {"cp": gas.cp_mole, "h": gas.enthalpy_mole, "rho": gas.density, "mu": gas.viscosity,
            "lam": gas.thermal_conductivity, "Pr": gas.cp_mass * gas.viscosity / gas.thermal_conductivity,
            "M": gas.mean_molecular_weight}


# ---------------------------------------------------------------------------
# Furnace-side closures
# ---------------------------------------------------------------------------
def h_conv_dittus_boelter(props: Mapping[str, float], n_kmol_h: float, geom: FurnaceGeometry,
                          f_ctube: float = 1.0) -> float:
    """Furnace-gas-to-surface convective coefficient [W/(m2 K)]: Nu = 0.023 Re^0.8 Pr^0.4 on the flue-gas
    flow passage (Latham 2008 Eq. 12 / Latham 2011 [24]), times ``f_ctube``."""
    m_dot = n_kmol_h / 3600.0 * props["M"]          # kg/s
    G = m_dot / geom.A_flue_free                    # kg/(m2 s)
    Re = G * geom.D_h / props["mu"]
    Nu = 0.023 * Re**0.8 * props["Pr"]**0.4
    return f_ctube * Nu * props["lam"] / geom.D_h


def furnace_flux(T_fg: float, T_wo: float, h_conv: float, F_gt: float) -> float:
    """Gas-to-tube flux per outer tube area [W/m2]: sigma F_gt (T_fg^4 - T_wo^4) + h_conv (T_fg - T_wo)."""
    return SIGMA_SB * F_gt * (T_fg**4 - T_wo**4) + h_conv * (T_fg - T_wo)


def solve_wall_temperature(T_fg: float, T_gas: float, U_inner: float, tube: r1.TubeGeometry,
                           h_conv: float, F_gt: float) -> float:
    """Outer-wall temperature from ``q_furnace(T_fg, T_wo) = U (T_wo - T_gas) d_i/d_o`` (bracketed root)."""
    ratio = tube.d_i / tube.d_o

    def g(T_wo):
        return furnace_flux(T_fg, T_wo, h_conv, F_gt) - U_inner * (T_wo - T_gas) * ratio

    lo, hi = (T_gas, T_fg) if T_fg >= T_gas else (T_fg, T_gas)
    if hi - lo < 1e-9:
        return float(T_gas)
    return float(brentq(g, lo, hi, xtol=1e-6))


# ---------------------------------------------------------------------------
# Coupled simulation
# ---------------------------------------------------------------------------
def simulate_coupled(tube: r1.TubeGeometry, bed: r1.CatalystBed, feed: r1.Feed, geom: FurnaceGeometry,
                     flue: FlueGas, release: HeatRelease = HeatRelease(), fparams: FurnaceParams = FurnaceParams(),
                     f_htg: float = 1.0, heat_transfer: str = "leva_grummer", alpha_i_const=None,
                     n_out: int = 201, method: str = "LSODA", rtol: float = 1e-6) -> CoupledResult:
    """Integrate tube and furnace together from the roof (z = 0) to the floor (z = L = geom.L).

    State ``[F_i (6), T_gas, P, T_fg]``; at each z the outer-wall temperature is solved from the flux
    balance. ``feed`` is per tube; ``flue`` is for the whole furnace (``geom.N_tubes`` tubes).
    """
    L = geom.L
    release.coefficients(L)  # validate the parabola early
    nS = len(r1.SPECIES)
    F0_map = feed.state_flows()
    F0 = np.array([F0_map[sp] for sp in r1.SPECIES])
    y0 = np.concatenate([F0, [feed.T_in, feed.P_in, flue.T_in_K]])
    n_fg_s = flue.n_kmol_h / 3600.0  # kmol/s
    n_eval = [0]

    def wall_temperature(z, F, T, P, T_fg):
        loc0 = r1._local(z, T, P, F, tube, bed, None, f_htg, False, False, heat_transfer, alpha_i_const, T_wo_value=T)
        fp = flue_props(T_fg, flue.X)
        h = h_conv_dittus_boelter(fp, flue.n_kmol_h, geom, fparams.f_ctube)
        T_wo = solve_wall_temperature(T_fg, T, loc0["U"], tube, h, fparams.F_gt)
        return T_wo, h, fp

    def rhs(z, y):
        n_eval[0] += 1
        F, T, P, T_fg = y[:nS], y[nS], y[nS + 1], y[nS + 2]
        T_wo, h, fp = wall_temperature(z, F, T, P, T_fg)
        dy_tube, loc = r1.tube_derivatives(z, F, T, P, tube, bed, None, f_htg, False, False, heat_transfer,
                                           alpha_i_const, T_wo_value=T_wo)
        q_o = loc["q_o"]                                         # W/m2 outer
        r_z = float(release.density(z, flue.Q_comb_W, L))        # W/m
        dT_fg = (r_z - q_o * math.pi * tube.d_o * geom.N_tubes) / (n_fg_s * fp["cp"])
        return np.concatenate([dy_tube, [dT_fg]])

    atol = np.concatenate([r1.default_atol(F0), [1e-6]])
    z_eval = np.linspace(0.0, L, n_out)
    t0 = time.perf_counter()
    sol = solve_ivp(rhs, (0.0, L), y0, method=method, t_eval=z_eval, rtol=rtol, atol=atol)
    wall_time = time.perf_counter() - t0

    z = sol.t
    Y = sol.y
    T_fg = Y[nS + 2]
    n = len(z)
    T_wo = np.zeros(n); h_arr = np.zeros(n)
    for k in range(n):
        T_wo[k], h_arr[k], _ = wall_temperature(z[k], Y[:nS, k], Y[nS, k], Y[nS + 1, k], T_fg[k])
    wall = r1.WallBC.from_table(z, T_wo, description="coupled furnace solution")
    tube_res = r1._build_result(z, Y[:nS + 2], tube, bed, feed, wall, f_htg, False, False, heat_transfer,
                                alpha_i_const, bool(sol.success), str(sol.message), n_eval[0], wall_time,
                                extras={"coupled": True})
    rel = release.density(z, flue.Q_comb_W, L)
    duty = float(np.trapezoid(tube_res.q_flux_outer * math.pi * tube.d_o * geom.N_tubes, z))
    h_in = flue_props(flue.T_in_K, flue.X)["h"]; h_out = flue_props(float(T_fg[-1]), flue.X)["h"]
    dH_fg = n_fg_s * (h_out - h_in)
    Q_eff = (1.0 - release.f_loss) * flue.Q_comb_W
    energy = {"Q_comb_W": flue.Q_comb_W, "Q_eff_W": Q_eff, "flue_enthalpy_rise_W": dH_fg, "tube_duty_W": duty,
              "release_integral_W": float(np.trapezoid(rel, z)), "closure_rel_error": (dH_fg + duty - Q_eff) / Q_eff,
              "duty_fraction_of_Q_comb": duty / flue.Q_comb_W}
    return CoupledResult(tube=tube_res, z=z, T_fg=T_fg, T_wo=T_wo, T_wi=tube_res.T_wall_inner,
                         q_outer=tube_res.q_flux_outer, release=rel, h_conv=h_arr, energy=energy,
                         success=bool(sol.success), message=str(sol.message), n_rhs_evals=n_eval[0],
                         wall_time_s=wall_time,
                         extras={"F_gt": fparams.F_gt, "f_ctube": fparams.f_ctube, "L_q": release.L_q,
                                 "alpha_top": release.alpha_top, "f_loss": release.f_loss, "T_fg_in_K": flue.T_in_K,
                                 "n_fg_kmol_h": flue.n_kmol_h, "excess_air_pct": flue.excess_air_pct})
