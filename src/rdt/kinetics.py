"""Xu & Froment (1989) intrinsic kinetics for methane steam reforming.

Implements the three-reaction Langmuir-Hinshelwood model of Xu and Froment
(AIChE J. 35 (1989) 88-96, Part I) on a Ni/MgAl2O4 catalyst:

    I    CH4 + H2O   <=> CO  + 3 H2
    II   CO  + H2O   <=> CO2 + H2
    III  CH4 + 2 H2O <=> CO2 + 4 H2

All numerical parameters are read from ``data/kinetics/xu_froment_1989.yaml``
(transcribed from Tables 5-7 of the paper); nothing that exists in that file is
hard-coded here. The module is pure NumPy except for
:func:`equilibrium_constants_cantera`, which uses Cantera standard-state Gibbs
energies (GRI-Mech 3.0 thermodynamics) to evaluate K1, K2, K3.

Units follow the paper: partial pressures in bar, rates in kmol/(kg_cat h),
k1 and k3 in kmol bar^0.5/(kg_cat h), k2 in kmol/(kg_cat h bar), adsorption
constants in bar^-1 (K_H2O dimensionless), energies in kJ/mol.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Dict, Mapping, Tuple, Union

import numpy as np
import yaml

ArrayLike = Union[float, np.ndarray]

DEFAULT_PARAMS_PATH = (
    Path(__file__).resolve().parents[2] / "data" / "kinetics" / "xu_froment_1989.yaml"
)

SPECIES = ("CH4", "H2O", "H2", "CO", "CO2")

#: Stoichiometric coefficients nu[species][reaction] for reactions I, II, III
#: (products positive, reactants negative).
STOICHIOMETRY: Dict[str, Tuple[int, int, int]] = {
    "CH4": (-1, 0, -1),
    "H2O": (-1, -1, -2),
    "H2": (3, 1, 4),
    "CO": (1, -1, 0),
    "CO2": (0, 1, 1),
}

#: Atom counts per species for the element-balance check.
ATOMS: Dict[str, Dict[str, int]] = {
    "CH4": {"C": 1, "H": 4, "O": 0},
    "H2O": {"C": 0, "H": 2, "O": 1},
    "H2": {"C": 0, "H": 2, "O": 0},
    "CO": {"C": 1, "H": 0, "O": 1},
    "CO2": {"C": 1, "H": 0, "O": 2},
}

#: Change in moles of gas per reaction (needed to convert K_p between pressure units).
DELTA_NU = (2, 0, 2)

#: Floor applied to the hydrogen partial pressure [bar]. The rate expressions
#: divide by p_H2^2.5, p_H2 and p_H2^3.5 (and DEN contains K_H2O p_H2O / p_H2),
#: so they diverge as p_H2 -> 0. Xu & Froment note (p. 95) that this is an
#: artefact of assuming H2 adsorption equilibrium and is harmless in practice
#: because industrial feeds always contain some hydrogen. The floor keeps the
#: expressions finite for a hydrogen-free inlet guess in a solver.
P_H2_FLOOR = 1.0e-6


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class XuFromentParams:
    """Xu & Froment (1989) parameter set at the reference activity.

    ``k_ref[i]`` and ``K_ref[j]`` are the values at the reference temperatures
    ``T_ref_k[i]`` / ``T_ref_K[j]`` (Table 5); ``E[i]`` and ``dH[j]`` are the
    activation energies and adsorption enthalpies in kJ/mol (Table 5);
    ``fresh_factor`` multiplies k1, k2, k3 for fresh catalyst (p. 94).
    ``ci95`` holds the 95 % confidence limits (UL, LL) of every parameter for
    later uncertainty quantification.
    """

    k_ref: Dict[int, float]
    E: Dict[int, float]
    T_ref_k: Dict[int, float]
    K_ref: Dict[str, float]
    dH: Dict[str, float]
    T_ref_K: Dict[str, float]
    fresh_factor: float
    R: float  # kJ/(mol K)
    ci95: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    source: str = ""

    # -- temperature dependence (Eqs. 6 and 7 of Part I) ---------------------
    def arrhenius(self, T: ArrayLike) -> Tuple[ArrayLike, ArrayLike, ArrayLike]:
        """Rate coefficients k1, k2, k3 at temperature ``T`` [K].

        Eq. (6): ``k_i = k_i,Tr * exp(-E_i/R * (1/T - 1/Tr))``.
        """
        T = np.asarray(T, dtype=float)
        return tuple(
            self.k_ref[i] * np.exp(-self.E[i] / self.R * (1.0 / T - 1.0 / self.T_ref_k[i]))
            for i in (1, 2, 3)
        )

    def adsorption(self, T: ArrayLike) -> Tuple[ArrayLike, ArrayLike, ArrayLike, ArrayLike]:
        """Adsorption constants K_CO, K_H2, K_CH4, K_H2O at ``T`` [K].

        Eq. (7): ``K_j = K_j,Tr * exp(-dH_j/R * (1/T - 1/Tr))``.
        """
        T = np.asarray(T, dtype=float)
        return tuple(
            self.K_ref[j] * np.exp(-self.dH[j] / self.R * (1.0 / T - 1.0 / self.T_ref_K[j]))
            for j in ("CO", "H2", "CH4", "H2O")
        )


def load_params(path: Union[str, Path] = DEFAULT_PARAMS_PATH) -> XuFromentParams:
    """Load :class:`XuFromentParams` from the project YAML file."""
    with Path(path).open() as f:
        d = yaml.safe_load(f)

    t5 = d["table5_reference_temperature_values"]
    t5e = d["table5_activation_energies_and_adsorption_enthalpies"]
    rep = d["reparameterised_forms"]
    tref = rep["T_ref_K"]

    k_ref = {i: float(t5[f"k{i}_{tref[f'k{i}']}"]["value"]) for i in (1, 2, 3)}
    E = {i: float(t5e[f"E{i}"]["value"]) for i in (1, 2, 3)}
    T_ref_k = {i: float(tref[f"k{i}"]) for i in (1, 2, 3)}

    K_ref = {j: float(t5[f"K_{j}_{tref[f'K_{j}']}"]["value"]) for j in ("CO", "H2", "CH4", "H2O")}
    dH = {j: float(t5e[f"dH_{j}"]["value"]) for j in ("CO", "H2", "CH4", "H2O")}
    T_ref_K = {j: float(tref[f"K_{j}"]) for j in ("CO", "H2", "CH4", "H2O")}

    ci95: Dict[str, Tuple[float, float]] = {}
    for block in (t5, t5e):
        for name, p in block.items():
            if isinstance(p, dict) and "UL" in p and "LL" in p:
                ci95[name] = (float(p["UL"]), float(p["LL"]))

    return XuFromentParams(
        k_ref=k_ref,
        E=E,
        T_ref_k=T_ref_k,
        K_ref=K_ref,
        dH=dH,
        T_ref_K=T_ref_K,
        fresh_factor=float(d["activity_factors"]["fresh_catalyst"]["value"]),
        R=float(rep["R_kJ_per_mol_K"]),
        ci95=ci95,
        source=str(d["metadata"].get("doi", "")),
    )


@lru_cache(maxsize=1)
def default_params() -> XuFromentParams:
    """Cached parameters from :data:`DEFAULT_PARAMS_PATH`."""
    return load_params(DEFAULT_PARAMS_PATH)


# ---------------------------------------------------------------------------
# Equilibrium constants
# ---------------------------------------------------------------------------
@lru_cache(maxsize=1)
def _gri30():
    import cantera as ct  # local import: Cantera is used only here

    gas = ct.Solution("gri30.yaml")
    idx = {s: gas.species_index(s) for s in SPECIES}
    return ct, gas, idx


def _equilibrium_constants_cantera_scalar(T: float) -> Tuple[float, float, float]:
    ct, gas, idx = _gri30()
    # Standard state at the reference pressure (1 atm) so that the standard-state
    # Gibbs energies carry no pressure term.
    gas.TP = float(T), ct.one_atm
    g_rt = gas.standard_gibbs_RT  # dimensionless g_i^0 / (R T)
    p_ref_bar = ct.one_atm / 1.0e5  # 1 atm expressed in bar (1.01325)
    out = []
    for r, dnu in enumerate(DELTA_NU):
        dg_rt = sum(STOICHIOMETRY[s][r] * g_rt[idx[s]] for s in SPECIES)
        # K_p in atm^dnu, then converted to bar^dnu: p[bar] = p[atm] * 1.01325
        out.append(float(np.exp(-dg_rt) * p_ref_bar**dnu))
    return out[0], out[1], out[2]


def equilibrium_constants_cantera(T: ArrayLike) -> Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """K1 [bar^2], K2 [-], K3 [bar^2] from Cantera (gri30.yaml) standard-state Gibbs energies.

    ``K_i = exp(-sum_j nu_ij g_j^0 / RT) * (p_ref / 1 bar)^{delta nu_i}`` with
    ``p_ref = 1 atm = 1.01325 bar``. Accepts a scalar or array of temperatures [K].
    """
    T_arr = np.atleast_1d(np.asarray(T, dtype=float))
    vals = np.array([_equilibrium_constants_cantera_scalar(t) for t in T_arr])
    K1, K2, K3 = vals[:, 0], vals[:, 1], vals[:, 2]
    if np.ndim(T) == 0:
        return float(K1[0]), float(K2[0]), float(K3[0])
    return K1, K2, K3


def equilibrium_constants_empirical(T: ArrayLike) -> Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """Empirical K1 [bar^2], K2 [-], K3 = K1*K2 [bar^2].

    ``K1 = exp(-26830/T + 30.114)`` and ``K2 = exp(4400/T - 4.036)`` are the
    correlations tabulated by Twigg (Catalyst Handbook, 2nd ed., 1989) and used
    for reformer simulation by Elnashaie and co-workers (e.g. Elnashaie &
    Elshishini, *Modelling, Simulation and Optimization of Industrial Fixed Bed
    Catalytic Reactors*, 1993). Provided as an independent cross-check on the
    Cantera values.
    """
    T = np.asarray(T, dtype=float)
    K1 = np.exp(-26830.0 / T + 30.114)
    K2 = np.exp(4400.0 / T - 4.036)
    return K1, K2, K1 * K2


def equilibrium_constants(T: ArrayLike, method: str = "cantera"):
    """Dispatch to :func:`equilibrium_constants_cantera` or ``_empirical``."""
    if method == "cantera":
        return equilibrium_constants_cantera(T)
    if method == "empirical":
        return equilibrium_constants_empirical(T)
    raise ValueError(f"unknown method {method!r}; use 'cantera' or 'empirical'")


# ---------------------------------------------------------------------------
# Rates
# ---------------------------------------------------------------------------
def _partial_pressures(p: Mapping[str, float]) -> Dict[str, np.ndarray]:
    missing = [s for s in SPECIES if s not in p]
    if missing:
        raise KeyError(f"partial pressures missing for {missing}")
    pp = {s: np.asarray(p[s], dtype=float) for s in SPECIES}
    pp["H2"] = np.maximum(pp["H2"], P_H2_FLOOR)
    return pp


def rates(
    T: ArrayLike,
    p: Mapping[str, float],
    activity: float = 1.0,
    fresh: bool = False,
    params: XuFromentParams | None = None,
    keq: Union[str, Tuple[ArrayLike, ArrayLike, ArrayLike]] = "cantera",
) -> Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """Intrinsic rates r1, r2, r3 [kmol/(kg_cat h)] at ``T`` [K] and partial pressures ``p`` [bar].

    Eq. (3) of Part I::

        r1 = k1 / p_H2^2.5 * (p_CH4 p_H2O   - p_H2^3 p_CO  / K1) / DEN^2
        r2 = k2 / p_H2     * (p_CO  p_H2O   - p_H2   p_CO2 / K2) / DEN^2
        r3 = k3 / p_H2^3.5 * (p_CH4 p_H2O^2 - p_H2^4 p_CO2 / K3) / DEN^2
        DEN = 1 + K_CO p_CO + K_H2 p_H2 + K_CH4 p_CH4 + K_H2O p_H2O / p_H2

    ``p`` must contain CH4, H2O, H2, CO, CO2. ``p_H2`` is floored at
    :data:`P_H2_FLOOR` (see its documentation). The result is multiplied by
    ``activity`` (e.g. Latham's f_prx) and, if ``fresh`` is true, by the
    fresh-catalyst factor 2.246 (Part I, p. 94). ``keq`` selects the source of
    K1, K2, K3 (``"cantera"``, ``"empirical"`` or a precomputed tuple).
    """
    prm = params or default_params()
    T = np.asarray(T, dtype=float)
    pp = _partial_pressures(p)

    k1, k2, k3 = prm.arrhenius(T)
    K_CO, K_H2, K_CH4, K_H2O = prm.adsorption(T)
    K1, K2, K3 = equilibrium_constants(T, keq) if isinstance(keq, str) else keq

    pH2, pCH4, pH2O, pCO, pCO2 = pp["H2"], pp["CH4"], pp["H2O"], pp["CO"], pp["CO2"]
    den = 1.0 + K_CO * pCO + K_H2 * pH2 + K_CH4 * pCH4 + K_H2O * pH2O / pH2

    r1 = k1 / pH2**2.5 * (pCH4 * pH2O - pH2**3 * pCO / K1) / den**2
    r2 = k2 / pH2 * (pCO * pH2O - pH2 * pCO2 / K2) / den**2
    r3 = k3 / pH2**3.5 * (pCH4 * pH2O**2 - pH2**4 * pCO2 / K3) / den**2

    scale = float(activity) * (prm.fresh_factor if fresh else 1.0)
    return r1 * scale, r2 * scale, r3 * scale


def species_rates(
    T: ArrayLike,
    p: Mapping[str, float],
    activity: float = 1.0,
    fresh: bool = False,
    params: XuFromentParams | None = None,
    keq: Union[str, Tuple[ArrayLike, ArrayLike, ArrayLike]] = "cantera",
) -> Dict[str, ArrayLike]:
    """Net production rates [kmol/(kg_cat h)] of CH4, H2O, H2, CO, CO2.

    ``R_k = sum_i nu_ik r_i`` with the stoichiometry of reactions I-III
    (:data:`STOICHIOMETRY`); consumption is negative.
    """
    r = rates(T, p, activity=activity, fresh=fresh, params=params, keq=keq)
    return {s: sum(nu * ri for nu, ri in zip(STOICHIOMETRY[s], r)) for s in SPECIES}


def element_balance(R: Mapping[str, ArrayLike]) -> Dict[str, ArrayLike]:
    """Net atom production rates of C, H and O implied by species rates ``R``.

    All three must be zero for a closed reaction network.
    """
    return {
        el: sum(ATOMS[s][el] * np.asarray(R[s]) for s in SPECIES) for el in ("C", "H", "O")
    }


def assert_element_balance(R: Mapping[str, ArrayLike], atol: float = 1e-10) -> None:
    """Raise ``AssertionError`` if C, H or O are not conserved to ``atol`` (relative to the largest rate)."""
    scale = max(float(np.max(np.abs(np.asarray(R[s])))) for s in SPECIES) or 1.0
    for el, val in element_balance(R).items():
        if not np.all(np.abs(val) <= atol * scale):
            raise AssertionError(f"{el} balance violated: {val} (tol {atol * scale:g})")
