"""One-dimensional pseudo-homogeneous plug-flow model of a single reformer tube.

Steady-state model of a catalyst-filled steam-methane-reformer tube with a
*prescribed* outer-wall temperature profile (no furnace model yet). Reaction
rates come from :mod:`rdt.kinetics` (Xu & Froment 1989); thermodynamic and
transport properties of the process gas come from Cantera (``gri30.yaml``) so
that heats of reaction, heat capacities and the equilibrium constants used by
the kinetics are mutually consistent.

State along the axial coordinate ``z`` [m]:

* molar flows ``F_i`` [kmol/h] of CH4, H2O, H2, CO, CO2, N2,
* gas temperature ``T`` [K],
* pressure ``P`` [bar].

Balances (per tube)::

    dF_i/dz = A_cs * rho_bed * sum_j nu_ij * eta_j * r_j                     (species)
    sum_i F_i Cp_i dT/dz = pi d_i U (T_wo(z) - T) - A_cs rho_bed sum_j eta_j r_j dH_j(T)   (energy)
    -dP/dz = f * rho_g * v_s^2 / d_p,  f = (1-phi)/phi^3 * (1.75 + 150 (1-phi)/Re_p)      (Ergun)

with ``1/U = 1/alpha_i + d_i/(2 lambda_tube) * ln(d_o/d_i)`` (Xu & Froment Part II,
Eq. 11; ``U`` referred to the inner tube surface), ``alpha_i`` from the
Leva-Grummer correlation as written in Latham et al. (2011) Eq. (20) with a
multiplier ``f_htg``, and the inner-wall temperature from the local flux
``T_wi = T_wo - q''_o * d_o/(2 lambda_tube) * ln(d_o/d_i)``.

Higher alkanes (C2H6, C3H8, C4H10, C5H12) in the feed are converted at the inlet
by Latham's overall rule ``C_k H_2k+2 + (k-1)/3 H2O -> (k-1)/3 CO + (2k+1)/3 CH4``
(Latham et al. 2011, Eq. 16). The small enthalpy of that step (about 3.7 kJ/mol
for ethane, Latham 2011 p. 1578) is neglected: the inlet temperature is kept.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Callable, Dict, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import yaml
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from rdt import kinetics as kin

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SPECIES: Tuple[str, ...] = ("CH4", "H2O", "H2", "CO", "CO2", "N2")
REACTING: Tuple[str, ...] = kin.SPECIES  # CH4, H2O, H2, CO, CO2
HIGHER_ALKANES: Dict[str, int] = {"C2H6": 2, "C3H8": 3, "C4H10": 4, "C5H12": 5}

ATOMS: Dict[str, Dict[str, int]] = {
    **kin.ATOMS,
    "N2": {"C": 0, "H": 0, "O": 0},
}

R_GAS = 8314.46261815324  # J/(kmol K)
BAR = 1.0e5  # Pa
DEFAULT_XF_PATH = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "literature_validation"
    / "xu_froment_1989_industrial.yaml"
)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class TubeGeometry:
    """Reformer tube geometry.

    d_i, d_o : inner / outer diameter [m]
    L_heated : heated length [m]
    lambda_tube : tube-wall thermal conductivity [W/(m K)]
    """

    d_i: float
    d_o: float
    L_heated: float
    lambda_tube: float

    @property
    def A_cs(self) -> float:
        """Inner cross-sectional area [m2]."""
        return math.pi * self.d_i**2 / 4.0

    @property
    def wall_resistance(self) -> float:
        """Conduction resistance of the wall referred to the inner surface, d_i/(2 lambda) ln(d_o/d_i) [m2 K/W]."""
        return self.d_i / (2.0 * self.lambda_tube) * math.log(self.d_o / self.d_i)


EtaSpec = Union[Sequence[float], Callable[[float], Sequence[float]]]


@dataclass(frozen=True)
class CatalystBed:
    """Catalyst bed properties.

    rho_bed : bulk catalyst density [kg_cat / m3 bed]
    voidage : bed void fraction phi [-]
    d_p : equivalent particle diameter [m]
    eta : effectiveness factors (eta1, eta2, eta3) or a callable ``eta(z)`` returning them.
          Default 0.1 for all three (Latham et al. 2011, p. 1578).
    activity : multiplier on all intrinsic rates (Latham's f_prx)
    """

    rho_bed: float
    voidage: float
    d_p: float
    eta: EtaSpec = (0.1, 0.1, 0.1)
    activity: float = 1.0

    def eta_at(self, z: float) -> np.ndarray:
        e = self.eta(z) if callable(self.eta) else self.eta
        return np.asarray(e, dtype=float)


@dataclass(frozen=True)
class Feed:
    """Inlet condition of one tube.

    T_in : temperature [K]
    P_in : pressure [bar]
    F : molar flows [kmol/h] keyed by species; may include C2H6, C3H8, C4H10, C5H12,
        which are converted to CH4/CO at the inlet by :func:`convert_higher_alkanes`.
    """

    T_in: float
    P_in: float
    F: Mapping[str, float]

    def state_flows(self) -> Dict[str, float]:
        """Molar flows of the six state species after higher-alkane conversion."""
        return convert_higher_alkanes(self.F)


class WallBC:
    """Prescribed outer-wall temperature ``T_wo(z)`` [K].

    Construct with :meth:`constant`, :meth:`linear` or :meth:`from_callable`.
    """

    def __init__(self, func: Callable[[float], float], description: str = "custom"):
        self._f = func
        self.description = description

    def __call__(self, z):
        return np.asarray(self._f(np.asarray(z, dtype=float)), dtype=float)

    @classmethod
    def constant(cls, T: float) -> "WallBC":
        return cls(lambda z: np.full_like(np.asarray(z, dtype=float), float(T)), f"constant {T} K")

    @classmethod
    def linear(cls, T_in: float, T_out: float, L: float) -> "WallBC":
        """Linear profile from ``T_in`` at z=0 to ``T_out`` at z=L (extrapolated beyond L)."""
        return cls(lambda z: T_in + (T_out - T_in) * np.asarray(z, dtype=float) / L,
                   f"linear {T_in} K -> {T_out} K over {L} m")

    @classmethod
    def from_callable(cls, f: Callable[[float], float], description: str = "callable") -> "WallBC":
        return cls(f, description)


@dataclass
class Result:
    """Axial profiles of a :func:`simulate` run. Arrays are indexed by ``z``."""

    z: np.ndarray                       # [m]
    T: np.ndarray                       # gas temperature [K]
    P: np.ndarray                       # [bar]
    F: Dict[str, np.ndarray]            # molar flows [kmol/h]
    X: Dict[str, np.ndarray]            # mole fractions (wet)
    conversion_CH4: np.ndarray          # (F_CH4,in - F_CH4)/F_CH4,in
    yield_CO2: np.ndarray               # (F_CO2 - F_CO2,in)/F_CH4,in
    q_flux_outer: np.ndarray            # heat flux at the OUTER tube surface [W/m2]
    q_flux_inner: np.ndarray            # heat flux at the inner tube surface [W/m2]
    T_wall_outer: np.ndarray            # [K]
    T_wall_inner: np.ndarray            # [K]
    dT_approach_I: np.ndarray           # T - T_eq(Q_I) [K]; nan where Q_I <= 0
    U: np.ndarray                       # overall coefficient, inner-area basis [W/(m2 K)]
    alpha_i: np.ndarray                 # bed-side coefficient [W/(m2 K)]
    rates: Dict[str, np.ndarray]        # effective rates eta_j r_j [kmol/(kg_cat h)] keyed r1, r2, r3
    F_in: Dict[str, float]
    success: bool
    message: str
    n_rhs_evals: int
    wall_time_s: float
    extras: Dict[str, object] = field(default_factory=dict)

    @property
    def dry_mole_fractions(self) -> Dict[str, np.ndarray]:
        dry = {s: self.X[s] for s in SPECIES if s != "H2O"}
        tot = sum(dry.values())
        return {s: v / tot for s, v in dry.items()}

    def outlet(self) -> Dict[str, float]:
        """Convenience summary of the last axial point."""
        dry = self.dry_mole_fractions
        return {
            "z_m": float(self.z[-1]),
            "T_K": float(self.T[-1]),
            "P_bar": float(self.P[-1]),
            "conversion_CH4": float(self.conversion_CH4[-1]),
            "yield_CO2": float(self.yield_CO2[-1]),
            "X_H2_dry": float(dry["H2"][-1]),
            "T_wall_outer_K": float(self.T_wall_outer[-1]),
            "T_wall_inner_K": float(self.T_wall_inner[-1]),
            "dT_approach_I_K": float(self.dT_approach_I[-1]),
        }


# ---------------------------------------------------------------------------
# Feed helpers
# ---------------------------------------------------------------------------
def convert_higher_alkanes(F: Mapping[str, float]) -> Dict[str, float]:
    """Apply Latham's inlet conversion rule for C2-C5 alkanes (Latham 2011, Eq. 16).

    ``C_k H_2k+2 + (k-1)/3 H2O -> (k-1)/3 CO + (2k+1)/3 CH4`` for each higher alkane.
    The rule conserves C, H and O and produces no net hydrogen. Its small
    enthalpy is neglected (see module docstring). Returns flows [kmol/h] of the
    six state species.
    """
    out = {s: float(F.get(s, 0.0)) for s in SPECIES}
    unknown = set(F) - set(SPECIES) - set(HIGHER_ALKANES)
    if unknown:
        raise KeyError(f"unsupported feed species {sorted(unknown)}")
    for sp, k in HIGHER_ALKANES.items():
        n = float(F.get(sp, 0.0))
        if n <= 0.0:
            continue
        out["H2O"] -= n * (k - 1) / 3.0
        out["CO"] += n * (k - 1) / 3.0
        out["CH4"] += n * (2 * k + 1) / 3.0
    if out["H2O"] < 0.0:
        raise ValueError("not enough steam to convert the higher alkanes at the inlet")
    return out


def ring_equivalent_diameter(d_pe: float, d_pi: float, H: float) -> float:
    """Equivalent diameter of a hollow cylinder: sphere with the same surface-to-volume ratio, 6 V/S [m].

    Brauer (1957), cited by Xu & Froment Part II p. 100 for the ring equivalent
    diameter, is not reproduced here; this is the plain Sauter-type definition
    including the inner hole surface. Flagged as a project choice.
    """
    V = math.pi / 4.0 * (d_pe**2 - d_pi**2) * H
    S = math.pi * (d_pe + d_pi) * H + 2.0 * math.pi / 4.0 * (d_pe**2 - d_pi**2)
    return 6.0 * V / S


def feed_from_xu_froment(path: Union[str, Path] = DEFAULT_XF_PATH) -> Feed:
    """Inlet of one tube for the Xu & Froment (1989) Part II Table 2 case.

    Uses the stated molar quantities: 5.168 kmol/h *equivalent* CH4 (treated as
    already containing the higher alkanes, so no C2+ is added) and the ratios
    H2O/CH4, CO2/CH4, H2/CH4, N2/CH4. T_in and P_in from the same table.
    """
    d = _load_yaml(path)
    n_ch4 = float(d["equivalent_CH4_feed_kmol_per_h"])
    r = d["molar_ratios_to_CH4"]
    F = {
        "CH4": n_ch4,
        "H2O": n_ch4 * float(r["H2O_over_CH4"]),
        "CO2": n_ch4 * float(r["CO2_over_CH4"]),
        "H2": n_ch4 * float(r["H2_over_CH4"]),
        "N2": n_ch4 * float(r["N2_over_CH4"]),
        "CO": 0.0,
    }
    return Feed(T_in=float(d["inlet"]["T0_K"]), P_in=float(d["inlet"]["p_t0_bar"]), F=F)


#: Tube-wall conductivity used when the source gives none: 106 500 J/(m h K) =
#: 29.6 W/(m K), Latham (2008) thesis Table 17 p. 93, citing Davis (2000) for a
#: reformer-tube stainless alloy. ASSUMPTION for the Xu & Froment tube (material not stated).
LAMBDA_TUBE_DEFAULT = 106500.0 / 3600.0

#: Bed voidage assumed for the Xu & Froment ring catalyst (not given in the paper;
#: "calculated according to Reichelt and Blasz (1971)" without a value). 0.6 is a
#: typical value for ring packings with the ring hole counted as bed void and is
#: close to Latham's fitted 0.607 for quadralobes. ASSUMPTION.
XF_VOIDAGE_ASSUMED = 0.6


def tube_from_xu_froment(path: Union[str, Path] = DEFAULT_XF_PATH,
                         lambda_tube: float = LAMBDA_TUBE_DEFAULT) -> TubeGeometry:
    """Tube geometry of the Xu & Froment Part II Table 2 case (ID 0.1016 m, OD 0.1322 m, heated 11.12 m)."""
    t = _load_yaml(path)["tube"]
    return TubeGeometry(d_i=float(t["internal_diameter_m"]), d_o=float(t["external_diameter_m"]),
                        L_heated=float(t["heated_length_m"]), lambda_tube=lambda_tube)


def bed_from_xu_froment(path: Union[str, Path] = DEFAULT_XF_PATH,
                        voidage: float = XF_VOIDAGE_ASSUMED,
                        eta: EtaSpec = (0.1, 0.1, 0.1), activity: float = 1.0) -> CatalystBed:
    """Catalyst bed for the Xu & Froment case.

    ``rho_bed = rho_s * (1 - voidage)`` with ``rho_s = 2355.2 kg/m3`` from Table 2
    interpreted as the density of the ring material, and the ring hole counted in
    the bed voidage (default :data:`XF_VOIDAGE_ASSUMED`). The equivalent particle
    diameter is :func:`ring_equivalent_diameter` of the Table 2 ring. Both the
    voidage and the equivalent-diameter definition are project assumptions.
    """
    c = _load_yaml(path)["catalyst_ring"]
    rho_s = float(c["rho_s_kg_per_m3"])
    d_p = ring_equivalent_diameter(float(c["d_pe_m"]), float(c["d_pi_m"]), float(c["H_m"]))
    return CatalystBed(rho_bed=rho_s * (1.0 - voidage), voidage=voidage, d_p=d_p, eta=eta, activity=activity)


@lru_cache(maxsize=4)
def _load_yaml(path: Union[str, Path]):
    with Path(path).open() as f:
        return yaml.safe_load(f)


# ---------------------------------------------------------------------------
# Gas properties (Cantera)
# ---------------------------------------------------------------------------
@lru_cache(maxsize=1)
def _gas():
    import cantera as ct

    gas = ct.Solution("gri30.yaml")  # mixture-averaged transport
    idx = {s: gas.species_index(s) for s in SPECIES}
    mw = {s: gas.molecular_weights[idx[s]] for s in SPECIES}  # kg/kmol
    return gas, idx, mw


def _gas_props(T: float, P_bar: float, X: Mapping[str, float]) -> Dict[str, object]:
    """Density, viscosity, conductivity, species Cp and reaction enthalpies at (T, P, X)."""
    gas, idx, _ = _gas()
    gas.TPX = float(T), float(P_bar) * BAR, {s: max(float(X[s]), 0.0) for s in SPECIES}
    cp = gas.partial_molar_cp  # J/(kmol K)  (= pure-species cp for an ideal gas)
    h = gas.partial_molar_enthalpies  # J/kmol
    dH = np.array([
        sum(kin.STOICHIOMETRY[s][r] * h[idx[s]] for s in REACTING) for r in range(3)
    ])  # J/kmol per reaction
    return {
        "rho": gas.density,                       # kg/m3
        "mu": gas.viscosity,                      # Pa s
        "lam": gas.thermal_conductivity,          # W/(m K)
        "cp": np.array([cp[idx[s]] for s in SPECIES]),
        "dH": dH,
    }


# ---------------------------------------------------------------------------
# Local closures
# ---------------------------------------------------------------------------
def leva_grummer_alpha(lam_g: float, mu: float, G_s: float, d_p: float, d_i: float,
                       f_htg: float = 1.0) -> float:
    """Bed-to-wall heat transfer coefficient [W/(m2 K)], Leva & Grummer (1948) as in Latham 2011 Eq. (20).

    ``alpha = f_htg * 0.813 * (lam_g / d_i) * exp(-6 d_p / d_i) * Re_p^0.9``,
    ``Re_p = d_p G_s / mu`` (Latham writes lam/(2 r_in) exp(-3 d_p / r_in), identical).
    """
    Re = d_p * G_s / mu
    return f_htg * 0.813 * lam_g / d_i * math.exp(-6.0 * d_p / d_i) * Re**0.9


def ergun_dPdz(rho: float, mu: float, v_s: float, d_p: float, phi: float) -> float:
    """Pressure gradient -dP/dz [Pa/m] from the Ergun equation (Latham 2011 Eq. 14 form)."""
    G_s = rho * v_s
    Re = d_p * G_s / mu
    f = (1.0 - phi) / phi**3 * (1.75 + 150.0 * (1.0 - phi) / Re)
    return f * rho * v_s**2 / d_p


def _local(z: float, T: float, P: float, F: np.ndarray, tube: TubeGeometry, bed: CatalystBed,
           wall: WallBC, f_htg: float, adiabatic: bool) -> Dict[str, object]:
    """Everything needed for the RHS and for post-processing at one axial point."""
    _, _, mw = _gas()
    F = np.maximum(F, 0.0)
    F_tot = F.sum()
    X = {s: F[i] / F_tot for i, s in enumerate(SPECIES)}
    props = _gas_props(T, P, X)

    m_dot = sum(F[i] * mw[s] for i, s in enumerate(SPECIES)) / 3600.0  # kg/s
    G_s = m_dot / tube.A_cs                                              # kg/(m2 s)
    v_s = G_s / props["rho"]                                             # m/s

    p = {s: P * X[s] for s in REACTING}                                  # bar
    r = np.array(kin.rates(T, p, activity=bed.activity))                 # kmol/(kg_cat h)
    eta = bed.eta_at(z)
    r_eff = eta * r

    T_wo = float(wall(z))
    if adiabatic:
        alpha_i = U = q_i = q_o = 0.0
        T_wi = T_wo
    else:
        alpha_i = leva_grummer_alpha(props["lam"], props["mu"], G_s, bed.d_p, tube.d_i, f_htg)
        U = 1.0 / (1.0 / alpha_i + tube.wall_resistance)
        q_i = U * (T_wo - T)                       # W/m2 at inner surface
        q_o = q_i * tube.d_i / tube.d_o            # W/m2 at outer surface
        T_wi = T_wo - q_o * tube.d_o / (2.0 * tube.lambda_tube) * math.log(tube.d_o / tube.d_i)

    return {"X": X, "props": props, "G_s": G_s, "v_s": v_s, "r_eff": r_eff, "r": r,
            "alpha_i": alpha_i, "U": U, "q_i": q_i, "q_o": q_o, "T_wo": T_wo, "T_wi": T_wi}


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------
def simulate(tube: TubeGeometry, bed: CatalystBed, feed: Feed, wall: WallBC,
             f_htg: float = 1.0, adiabatic: bool = False, L: Optional[float] = None,
             n_out: int = 201, method: str = "LSODA", rtol: float = 1e-7,
             atol: Optional[np.ndarray] = None) -> Result:
    """Integrate the tube model from z = 0 to ``L`` (default ``tube.L_heated``)."""
    L = tube.L_heated if L is None else float(L)
    F0_map = feed.state_flows()
    F0 = np.array([F0_map[s] for s in SPECIES])
    y0 = np.concatenate([F0, [feed.T_in, feed.P_in]])
    nS = len(SPECIES)
    n_eval = [0]

    def rhs(z, y):
        n_eval[0] += 1
        F, T, P = y[:nS], y[nS], y[nS + 1]
        loc = _local(z, T, P, F, tube, bed, wall, f_htg, adiabatic)
        pr = loc["props"]
        r_eff = loc["r_eff"]
        dF = np.zeros(nS)
        for i, s in enumerate(REACTING):
            dF[i] = tube.A_cs * bed.rho_bed * sum(nu * rj for nu, rj in zip(kin.STOICHIOMETRY[s], r_eff))
        # energy: J/(h m)
        q_lin = 3600.0 * math.pi * tube.d_i * loc["q_i"]
        q_rxn = tube.A_cs * bed.rho_bed * float(np.dot(r_eff, pr["dH"]))
        cp_flow = float(np.dot(np.maximum(F, 0.0), pr["cp"]))  # J/(h K)
        dT = (q_lin - q_rxn) / cp_flow
        dP = -ergun_dPdz(pr["rho"], pr["mu"], loc["v_s"], bed.d_p, bed.voidage) / BAR
        return np.concatenate([dF, [dT, dP]])

    if atol is None:
        atol = np.concatenate([np.full(nS, 1e-9 * max(F0.sum(), 1e-12)), [1e-6, 1e-8]])
    z_eval = np.linspace(0.0, L, n_out)
    t0 = time.perf_counter()
    sol = solve_ivp(rhs, (0.0, L), y0, method=method, t_eval=z_eval, rtol=rtol, atol=atol)
    wall_time = time.perf_counter() - t0

    z = sol.t
    Y = sol.y
    F = {s: Y[i] for i, s in enumerate(SPECIES)}
    T = Y[nS]
    P = Y[nS + 1]

    n = len(z)
    X = {s: np.zeros(n) for s in SPECIES}
    q_o = np.zeros(n); q_i = np.zeros(n); T_wo = np.zeros(n); T_wi = np.zeros(n)
    U = np.zeros(n); alpha = np.zeros(n); r1 = np.zeros(n); r2 = np.zeros(n); r3 = np.zeros(n)
    dT_app = np.full(n, np.nan)
    for k in range(n):
        Fk = Y[:nS, k]
        loc = _local(z[k], T[k], P[k], Fk, tube, bed, wall, f_htg, adiabatic)
        for s in SPECIES:
            X[s][k] = loc["X"][s]
        q_o[k], q_i[k], T_wo[k], T_wi[k] = loc["q_o"], loc["q_i"], loc["T_wo"], loc["T_wi"]
        U[k], alpha[k] = loc["U"], loc["alpha_i"]
        r1[k], r2[k], r3[k] = loc["r_eff"]
        dT_app[k] = approach_to_equilibrium_I(T[k], {s: P[k] * loc["X"][s] for s in REACTING})

    F_ch4_in = F0_map["CH4"]
    return Result(
        z=z, T=T, P=P, F=F, X=X,
        conversion_CH4=(F_ch4_in - F["CH4"]) / F_ch4_in,
        yield_CO2=(F["CO2"] - F0_map["CO2"]) / F_ch4_in,
        q_flux_outer=q_o, q_flux_inner=q_i, T_wall_outer=T_wo, T_wall_inner=T_wi,
        dT_approach_I=dT_app, U=U, alpha_i=alpha, rates={"r1": r1, "r2": r2, "r3": r3},
        F_in=F0_map, success=bool(sol.success), message=str(sol.message),
        n_rhs_evals=n_eval[0], wall_time_s=wall_time,
        extras={"wall_bc": wall.description, "method": method, "f_htg": f_htg, "adiabatic": adiabatic},
    )


def approach_to_equilibrium_I(T: float, p: Mapping[str, float],
                              T_bounds: Tuple[float, float] = (400.0, 2000.0)) -> float:
    """Approach-to-equilibrium temperature for reaction I, ``T - T_eq`` [K].

    ``T_eq`` solves ``K1(T_eq) = p_H2^3 p_CO / (p_CH4 p_H2O)`` with Cantera K1.
    Positive when the gas is hotter than the temperature at which its
    composition would be at equilibrium (Xu & Froment Part II, Fig. 5 "delta T").
    Returns nan if the reaction quotient is not positive (e.g. no CO yet).
    """
    Q = p["H2"] ** 3 * p["CO"] / (p["CH4"] * p["H2O"]) if p["CH4"] > 0 and p["H2O"] > 0 else 0.0
    if not np.isfinite(Q) or Q <= 0.0:
        return float("nan")
    lnQ = math.log(Q)

    def g(Te):
        return math.log(kin.equilibrium_constants_cantera(Te)[0]) - lnQ

    lo, hi = T_bounds
    if g(lo) * g(hi) > 0:
        return float("nan")
    return T - brentq(g, lo, hi, xtol=1e-6)


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------
def element_flows(F: Mapping[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Atomic flows [kmol atoms/h] of C, H, O (and N) from species molar flows."""
    out = {}
    for el in ("C", "H", "O"):
        out[el] = sum(ATOMS[s][el] * np.asarray(F[s]) for s in SPECIES)
    out["N"] = 2.0 * np.asarray(F["N2"])
    return out
