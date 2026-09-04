"""Tests for rdt.kinetics (Xu & Froment 1989)."""

import cantera as ct
import numpy as np
import pytest
import yaml

from rdt import kinetics as kin

FEED = {"CH4": 1.0, "H2O": 3.0, "H2": 0.1}  # molar proportions, un-reacted feed


def partial_pressures(composition, P_bar):
    tot = sum(composition.values())
    p = {s: P_bar * composition.get(s, 0.0) / tot for s in kin.SPECIES}
    return p


@pytest.fixture(scope="module")
def params():
    return kin.load_params()


@pytest.fixture(scope="module")
def raw_yaml():
    with kin.DEFAULT_PARAMS_PATH.open() as f:
        return yaml.safe_load(f)


def test_k_at_reference_temperatures_match_yaml(params, raw_yaml):
    t5 = raw_yaml["table5_reference_temperature_values"]
    for i, key in ((1, "k1_648"), (2, "k2_648"), (3, "k3_648")):
        k = params.arrhenius(params.T_ref_k[i])[i - 1]
        assert float(k) == t5[key]["value"]


def test_K_at_reference_temperatures_match_yaml(params, raw_yaml):
    t5 = raw_yaml["table5_reference_temperature_values"]
    K = dict(zip(("CO", "H2", "CH4", "H2O"), [None] * 4))
    for j, key in (("CO", "K_CO_648"), ("H2", "K_H2_648"), ("CH4", "K_CH4_823"), ("H2O", "K_H2O_823")):
        idx = ("CO", "H2", "CH4", "H2O").index(j)
        assert float(params.adsorption(params.T_ref_K[j])[idx]) == t5[key]["value"]


def test_k1_at_823K_close_to_table4(params):
    """Table 4 (Part I, p. 93) per-temperature estimate of k1 at 823 K is 2.069."""
    k1 = float(params.arrhenius(823.0)[0])
    assert abs(k1 / 2.069 - 1.0) < 0.30


def test_cantera_vs_empirical_equilibrium_constants():
    T = np.linspace(800.0, 1100.0, 13)
    K1c, K2c, K3c = kin.equilibrium_constants_cantera(T)
    K1e, K2e, K3e = kin.equilibrium_constants_empirical(T)
    assert np.all(np.abs(K1c / K1e - 1.0) < 0.10)
    assert np.all(np.abs(K2c / K2e - 1.0) < 0.10)
    assert np.all(np.abs(K3c / K3e - 1.0) < 0.10)


def test_rates_vanish_at_cantera_equilibrium():
    T, P_bar = 1000.0, 25.0
    gas = ct.Solution("gri30.yaml")
    gas.TPX = T, P_bar * 1e5, {"CH4": 1.0, "H2O": 3.0}
    gas.equilibrate("TP")
    X = gas.mole_fraction_dict()
    p_eq = {s: P_bar * X.get(s, 0.0) for s in kin.SPECIES}
    r_eq = kin.rates(T, p_eq)
    r_feed = kin.rates(T, partial_pressures(FEED, P_bar))
    # r2 of the un-reacted feed is identically zero (no CO, no CO2), so the
    # characteristic rate scale is the largest feed rate (r3 here).
    scale = max(abs(float(r)) for r in r_feed)
    assert scale > 0
    for req in r_eq:
        assert abs(float(req)) < 1e-3 * scale, (req, scale)


def test_unreacted_feed_rates_positive():
    r1, r2, r3 = kin.rates(900.0, partial_pressures(FEED, 25.0))
    assert r1 > 0
    assert r3 > 0


def test_species_rates_conserve_elements():
    for T in (800.0, 900.0, 1000.0):
        p = partial_pressures({"CH4": 1.0, "H2O": 3.0, "H2": 0.5, "CO": 0.1, "CO2": 0.2}, 25.0)
        R = kin.species_rates(T, p)
        kin.assert_element_balance(R, atol=1e-10)
        bal = kin.element_balance(R)
        scale = max(abs(float(R[s])) for s in kin.SPECIES)
        for el in ("C", "H", "O"):
            assert abs(float(bal[el])) <= 1e-10 * scale


def test_fresh_factor_scales_all_rates(params):
    p = partial_pressures({"CH4": 1.0, "H2O": 3.0, "H2": 0.5, "CO": 0.1, "CO2": 0.2}, 25.0)
    r_ref = kin.rates(900.0, p)
    r_fresh = kin.rates(900.0, p, fresh=True)
    assert params.fresh_factor == 2.246
    for a, b in zip(r_ref, r_fresh):
        assert b == pytest.approx(2.246 * a, rel=1e-12)


def test_activity_scales_rates():
    p = partial_pressures({"CH4": 1.0, "H2O": 3.0, "H2": 0.5, "CO": 0.1, "CO2": 0.2}, 25.0)
    r_ref = kin.rates(900.0, p)
    r_act = kin.rates(900.0, p, activity=0.2)
    for a, b in zip(r_ref, r_act):
        assert b == pytest.approx(0.2 * a, rel=1e-12)


def test_hydrogen_floor_keeps_rates_finite():
    p = partial_pressures({"CH4": 1.0, "H2O": 3.0}, 25.0)  # no hydrogen
    r = kin.rates(900.0, p)
    assert all(np.isfinite(r))
