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
by one of two rules (``Feed.inlet_higher_alkanes``):

* ``"latham"``: ``C_k H_2k+2 + (k-1)/3 H2O -> (k-1)/3 CO + (2k+1)/3 CH4`` (Latham
  et al. 2011, Eq. 16; no net hydrogen);
* ``"xu_froment"``: ``C_k H_2k+2 + k H2O -> k CO + (2k+1) H2`` (complete steam
  reforming of the C2+ fraction, so the equivalent-CH4 conversion starts at the
  C2+ carbon fraction of the feed, about 0.09 for the Xu & Froment Table 2 gas).

The chemical engineer must confirm which definition Xu & Froment (1989) used
for the conversion plotted in Part II Fig. 3. The small enthalpy of the inlet
step is neglected in both cases: the inlet temperature is kept.

Bed-side heat transfer (``simulate(..., heat_transfer=...)``):

* ``"leva_grummer"`` (default): Leva & Grummer (1948) as in Latham 2011 Eq. (20);
* ``"xu_froment"``: Xu & Froment Part II Eqs. (11)-(12) with the De Wasch &
  Froment / Yagi-Kunii dynamic terms and the Kunii & Smith (1960) static bed
  conductivity, see :func:`xu_froment_alpha_i`.
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
    lambda_s : pellet thermal conductivity [W/(m K)]; default 1.0 is a literature-range value
          for porous ceramic (MgAl2O4 / alumina) supported Ni catalyst (typ. 0.5-2). ASSUMPTION.
    emissivity : pellet surface emissivity; default 0.8 for oxide ceramics at 800-1100 K
          (typ. 0.7-0.9). ASSUMPTION. Both only enter the Kunii-Smith static conductivity.
    """

    rho_bed: float
    voidage: float
    d_p: float
    eta: EtaSpec = (0.1, 0.1, 0.1)
    activity: float = 1.0
    lambda_s: float = 1.0      # solid (pellet) thermal conductivity [W/(m K)], used by Kunii-Smith only
    emissivity: float = 0.8    # pellet surface emissivity [-], used by Kunii-Smith only

    def eta_at(self, z: float) -> np.ndarray:
        e = self.eta(z) if callable(self.eta) else self.eta
        return np.asarray(e, dtype=float)


@dataclass(frozen=True)
class Feed:
    """Inlet condition of one tube.

    T_in : temperature [K]
    P_in : pressure [bar]
    F : molar flows [kmol/h] keyed by species; may include C2H6, C3H8, C4H10, C5H12,
        which are converted at the inlet by :func:`convert_higher_alkanes`.
    inlet_higher_alkanes : "latham" (C2+ -> CO + CH4, Latham 2011 Eq. 16) or "xu_froment"
        (C2+ -> CO + H2, complete reforming). Conversion is always reported relative to the
        carbon-equivalent CH4 feed, so the two rules give different inlet conversions
        (about 0.02 and 0.09 for the Xu & Froment gas). The chemical engineer must confirm
        which definition Xu & Froment used.
    """

    T_in: float
    P_in: float
    F: Mapping[str, float]
    inlet_higher_alkanes: str = "latham"   # "latham" or "xu_froment", see convert_higher_alkanes

    def state_flows(self) -> Dict[str, float]:
        """Molar flows of the six state species after higher-alkane conversion."""
        return convert_higher_alkanes(self.F, rule=self.inlet_higher_alkanes)

    @property
    def F_CH4_equivalent(self) -> float:
        """Carbon-equivalent methane feed [kmol/h]: CH4 + sum_k k * F(C_k H_2k+2). Conversion basis."""
        return float(self.F.get("CH4", 0.0)) + sum(k * float(self.F.get(sp, 0.0)) for sp, k in HIGHER_ALKANES.items())


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

    @classmethod
    def from_table(cls, z: Sequence[float], T: Sequence[float], description: str = "table (PCHIP)") -> "WallBC":
        """Shape-preserving (PCHIP) interpolation of tabulated ``T_wo(z)``, held constant beyond the table ends."""
        from scipy.interpolate import PchipInterpolator

        z = np.asarray(z, dtype=float); T = np.asarray(T, dtype=float)
        order = np.argsort(z); z, T = z[order], T[order]
        keep = np.concatenate([[True], np.diff(z) > 0]); z, T = z[keep], T[keep]
        f = PchipInterpolator(z, T, extrapolate=False)
        lo, hi = z[0], z[-1]

        def g(zz):
            zz = np.asarray(zz, dtype=float)
            return f(np.clip(zz, lo, hi))

        return cls(g, description)


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
def convert_higher_alkanes(F: Mapping[str, float], rule: str = "latham") -> Dict[str, float]:
    """Convert C2-C5 alkanes at the inlet and return flows [kmol/h] of the six state species.

    ``rule="latham"``: ``C_k H_2k+2 + (k-1)/3 H2O -> (k-1)/3 CO + (2k+1)/3 CH4``
    (Latham et al. 2011, Eq. 16; conserves C, H, O; no net hydrogen).
    ``rule="xu_froment"``: ``C_k H_2k+2 + k H2O -> k CO + (2k+1) H2`` (complete
    steam reforming of the C2+ fraction to CO and H2, so that the carbon-
    equivalent CH4 conversion at z = 0 equals the C2+ carbon fraction).
    The enthalpy of either step is neglected (see module docstring).
    """
    if rule not in ("latham", "xu_froment"):
        raise ValueError(f"unknown inlet rule {rule!r}; use 'latham' or 'xu_froment'")
    out = {s: float(F.get(s, 0.0)) for s in SPECIES}
    unknown = set(F) - set(SPECIES) - set(HIGHER_ALKANES)
    if unknown:
        raise KeyError(f"unsupported feed species {sorted(unknown)}")
    for sp, k in HIGHER_ALKANES.items():
        n = float(F.get(sp, 0.0))
        if n <= 0.0:
            continue
        if rule == "latham":
            out["H2O"] -= n * (k - 1) / 3.0
            out["CO"] += n * (k - 1) / 3.0
            out["CH4"] += n * (2 * k + 1) / 3.0
        else:
            out["H2O"] -= n * k
            out["CO"] += n * k
            out["H2"] += n * (2 * k + 1)
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


def feed_from_xu_froment(path: Union[str, Path] = DEFAULT_XF_PATH,
                         split_alkanes: bool = False,
                         inlet_higher_alkanes: str = "latham") -> Feed:
    """Inlet of one tube for the Xu & Froment (1989) Part II Table 2 case.

    Uses the stated molar quantities: 5.168 kmol/h *equivalent* CH4 and the ratios
    H2O/CH4, CO2/CH4, H2/CH4, N2/CH4; T_in and P_in from the same table.

    ``split_alkanes=False`` (default, previous behaviour): all 5.168 kmol/h are
    fed as CH4. ``split_alkanes=True``: the equivalent CH4 is distributed over
    CH4, C2H6, C3H8, C4H10, C5H12 in the carbon proportions of the Table 2
    natural-gas analysis (CH4 81.5, C2H6 2.8, C3H8 0.4, C4H10 0.1, C5H12 0.2 vol%),
    so that the chosen ``inlet_higher_alkanes`` rule sets the inlet conversion.
    The carbon-equivalent CH4 feed (conversion basis) is 5.168 kmol/h either way.
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
    if split_alkanes:
        ng = d["natural_gas_composition_vol_pct"]
        carbon = {"CH4": 1, "C2H6": 2, "C3H8": 3, "C4H10": 4, "C5H12": 5}
        total_c = sum(k * float(ng[sp]) for sp, k in carbon.items())
        for sp, k in carbon.items():
            F[sp] = n_ch4 * float(ng[sp]) / total_c   # moles of species per hour
    return Feed(T_in=float(d["inlet"]["T0_K"]), P_in=float(d["inlet"]["p_t0_bar"]), F=F,
                inlet_higher_alkanes=inlet_higher_alkanes)


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
        "cp_mass": gas.cp_mass,                   # J/(kg K)
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


SIGMA_SB = 5.670374419e-8  # W/(m2 K4)


def kunii_smith_lambda_er0(lam_g: float, lam_s: float, eps: float, T: float, emissivity: float,
                           d_p: float, beta: float = 1.0, gamma: float = 2.0 / 3.0) -> float:
    """Static effective radial bed conductivity, Kunii & Smith (1960) [W/(m K)].

    ``lam_er0/lam_g = eps (1 + beta h_rv d_p/lam_g) + beta (1-eps) / [ 1/(1/phi + h_rs d_p/lam_g) + gamma lam_g/lam_s ]``
    (Kunii & Smith 1960, Eq. 19; Froment & Bischoff 1979, Eq. 11.7.1-2), with
    ``h_rv = 4 sigma T^3 / (1 + eps/(2(1-eps)) (1-p)/p)`` (void-to-void radiation),
    ``h_rs = 4 sigma T^3 p/(2-p)`` (solid-to-solid radiation), ``p`` the emissivity,
    ``beta = 1`` and ``gamma = 2/3``. ``phi`` interpolates between the Kunii-Smith
    dense-packing (eps = 0.26) and loose-packing (eps = 0.476) values ``phi_2``,
    ``phi_1`` evaluated from their analytical expression with ``kappa = lam_s/lam_g``,
    ``sin^2 theta = 1/n``, ``n = 4 sqrt(3)`` (dense) and ``n = 1.5`` (loose); for
    eps > 0.476 (this ring bed) the loose value ``phi_1`` is used.
    In the original kcal units ``4 sigma`` is the familiar 0.1952 (T/100)^3.
    """
    kappa = lam_s / lam_g

    def phi_i(n: float) -> float:
        s2 = 1.0 / n
        c = math.sqrt(1.0 - s2)
        a = (kappa - 1.0) / kappa
        return 0.5 * a * a * s2 / (math.log(kappa - (kappa - 1.0) * c) - a * (1.0 - c)) - 2.0 / (3.0 * kappa)

    phi1, phi2 = phi_i(1.5), phi_i(4.0 * math.sqrt(3.0))
    if eps <= 0.26:
        phi = phi2
    elif eps >= 0.476:
        phi = phi1
    else:
        phi = phi2 + (phi1 - phi2) * (eps - 0.26) / 0.216
    p = emissivity
    h_rv = 4.0 * SIGMA_SB * T**3 / (1.0 + eps / (2.0 * (1.0 - eps)) * (1.0 - p) / p)
    h_rs = 4.0 * SIGMA_SB * T**3 * p / (2.0 - p)
    return lam_g * (eps * (1.0 + beta * h_rv * d_p / lam_g)
                    + beta * (1.0 - eps) / (1.0 / (1.0 / phi + h_rs * d_p / lam_g) + gamma / kappa))


def xu_froment_alpha_i(lam_g: float, mu: float, cp_mass: float, G_s: float, d_p: float, d_ti: float,
                       lam_er0: float) -> Dict[str, float]:
    """Bed-side coefficient alpha_i [W/(m2 K)] of Xu & Froment (1989) Part II, Eq. (12) and the
    two correlations below it (p. 100), in SI.

    ``Re = d_p G_s / mu``, ``Pr = cp mu / lam_g`` (superficial mass velocity, particle diameter);
    ``alpha_w0 = 8.694 lam_er0 / d_ti^(4/3)`` (De Wasch & Froment 1972; the constant carries
    m^(1/3) and is unchanged between kJ/(m h K) and W/(m K) because alpha and lambda scale alike);
    ``alpha_w = alpha_w0 + 0.444 Re Pr lam_g / d_p``;
    ``lam_er = lam_er0 + 0.14 lam_g Re Pr`` (Yagi & Kunii form);
    ``alpha_i = 8 lam_er alpha_w / (8 lam_er + alpha_w d_ti)``  (Eq. 12, i.e. 1/alpha_i = 1/alpha_w + d_ti/(8 lam_er)).
    Returns a dict with alpha_i and the intermediate quantities.
    """
    Re = d_p * G_s / mu
    Pr = cp_mass * mu / lam_g
    alpha_w0 = 8.694 * lam_er0 / d_ti ** (4.0 / 3.0)
    alpha_w = alpha_w0 + 0.444 * Re * Pr * lam_g / d_p
    lam_er = lam_er0 + 0.14 * lam_g * Re * Pr
    alpha_i = 8.0 * lam_er * alpha_w / (8.0 * lam_er + alpha_w * d_ti)
    return {"alpha_i": alpha_i, "alpha_w0": alpha_w0, "alpha_w": alpha_w, "lam_er": lam_er,
            "lam_er0": lam_er0, "Re": Re, "Pr": Pr}


HEAT_TRANSFER_OPTIONS = ("leva_grummer", "xu_froment", "constant")


def _local(z: float, T: float, P: float, F: np.ndarray, tube: TubeGeometry, bed: CatalystBed,
           wall: WallBC, f_htg: float, adiabatic: bool,
           adiabatic_beyond_heated: bool = False, heat_transfer: str = "leva_grummer",
           alpha_i_const: Union[float, Callable[[float], float], None] = None) -> Dict[str, object]:
    """Everything needed for the RHS and for post-processing at one axial point.

    If ``adiabatic_beyond_heated`` is true, no wall heat is exchanged for
    ``z > tube.L_heated`` (unheated tube tail); the wall temperatures are then
    reported equal to the gas temperature.
    """
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
    if adiabatic_beyond_heated and z > tube.L_heated:
        alpha_i = U = q_i = q_o = 0.0
        T_wo = T_wi = float(T)
    elif adiabatic:
        alpha_i = U = q_i = q_o = 0.0
        T_wi = T_wo
    else:
        if heat_transfer == "leva_grummer":
            alpha_i = leva_grummer_alpha(props["lam"], props["mu"], G_s, bed.d_p, tube.d_i, f_htg)
        elif heat_transfer == "xu_froment":
            lam0 = kunii_smith_lambda_er0(props["lam"], bed.lambda_s, bed.voidage, T, bed.emissivity, bed.d_p)
            alpha_i = f_htg * xu_froment_alpha_i(props["lam"], props["mu"], props["cp_mass"], G_s,
                                                 bed.d_p, tube.d_i, lam0)["alpha_i"]
        elif heat_transfer == "constant":
            if alpha_i_const is None:
                raise ValueError("heat_transfer='constant' requires alpha_i_const (float or callable of z)")
            alpha_i = f_htg * float(alpha_i_const(z) if callable(alpha_i_const) else alpha_i_const)
        else:
            raise ValueError(f"unknown heat_transfer {heat_transfer!r}; use {HEAT_TRANSFER_OPTIONS}")
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
             atol: Optional[np.ndarray] = None, adiabatic_beyond_heated: bool = False,
             heat_transfer: str = "leva_grummer",
             alpha_i_const: Union[float, Callable[[float], float], None] = None) -> Result:
    """Integrate the tube model from z = 0 to ``L`` (default ``tube.L_heated``).

    With ``adiabatic_beyond_heated=True`` the section ``tube.L_heated < z <= L`` exchanges
    no heat with the wall (unheated tail, as in Xu & Froment Part II Fig. 3 for 11.12-12 m).
    ``heat_transfer`` selects the bed-side coefficient: ``"leva_grummer"`` (Latham 2011
    Eq. 20, scaled by ``f_htg``) or ``"xu_froment"`` (Part II Eqs. 11-12 with Kunii-Smith
    static conductivity; ``f_htg`` also multiplies it, default 1) or ``"constant"``, which uses
    ``alpha_i_const`` [W/(m2 K)] directly, either a number or a callable ``alpha_i(z)`` (e.g. a
    coefficient profile back-calculated from measured wall and gas temperatures). CH4 conversion
    is reported relative to the carbon-equivalent CH4 feed (``feed.F_CH4_equivalent``).
    """
    if heat_transfer not in HEAT_TRANSFER_OPTIONS:
        raise ValueError(f"unknown heat_transfer {heat_transfer!r}; use {HEAT_TRANSFER_OPTIONS}")
    L = tube.L_heated if L is None else float(L)
    F0_map = feed.state_flows()
    F0 = np.array([F0_map[s] for s in SPECIES])
    y0 = np.concatenate([F0, [feed.T_in, feed.P_in]])
    nS = len(SPECIES)
    n_eval = [0]

    def rhs(z, y):
        n_eval[0] += 1
        F, T, P = y[:nS], y[nS], y[nS + 1]
        loc = _local(z, T, P, F, tube, bed, wall, f_htg, adiabatic, adiabatic_beyond_heated, heat_transfer, alpha_i_const)
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
        loc = _local(z[k], T[k], P[k], Fk, tube, bed, wall, f_htg, adiabatic, adiabatic_beyond_heated, heat_transfer, alpha_i_const)
        for s in SPECIES:
            X[s][k] = loc["X"][s]
        q_o[k], q_i[k], T_wo[k], T_wi[k] = loc["q_o"], loc["q_i"], loc["T_wo"], loc["T_wi"]
        U[k], alpha[k] = loc["U"], loc["alpha_i"]
        r1[k], r2[k], r3[k] = loc["r_eff"]
        dT_app[k] = approach_to_equilibrium_I(T[k], {s: P[k] * loc["X"][s] for s in REACTING})

    F_ch4_in = feed.F_CH4_equivalent
    return Result(
        z=z, T=T, P=P, F=F, X=X,
        conversion_CH4=(F_ch4_in - F["CH4"]) / F_ch4_in,
        yield_CO2=(F["CO2"] - F0_map["CO2"]) / F_ch4_in,
        q_flux_outer=q_o, q_flux_inner=q_i, T_wall_outer=T_wo, T_wall_inner=T_wi,
        dT_approach_I=dT_app, U=U, alpha_i=alpha, rates={"r1": r1, "r2": r2, "r3": r3},
        F_in=F0_map, success=bool(sol.success), message=str(sol.message),
        n_rhs_evals=n_eval[0], wall_time_s=wall_time,
        extras={"wall_bc": wall.description, "method": method, "f_htg": f_htg, "adiabatic": adiabatic,
                "adiabatic_beyond_heated": adiabatic_beyond_heated, "heat_transfer": heat_transfer,
                "inlet_higher_alkanes": feed.inlet_higher_alkanes, "F_CH4_equivalent": F_ch4_in},
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
