"""Tests for rdt.reactor1d (1-D pseudo-homogeneous tube with prescribed wall temperature)."""

import time

import cantera as ct
import numpy as np
import pytest

from rdt import reactor1d as r1


@pytest.fixture(scope="module")
def xf_case():
    tube = r1.tube_from_xu_froment()
    bed = r1.bed_from_xu_froment()
    feed = r1.feed_from_xu_froment()
    return tube, bed, feed


@pytest.fixture(scope="module")
def xf_result(xf_case):
    tube, bed, feed = xf_case
    wall = r1.WallBC.linear(1000.0, 1150.0, tube.L_heated)
    t0 = time.perf_counter()
    res = r1.simulate(tube, bed, feed, wall)
    res.extras["test_wall_time_s"] = time.perf_counter() - t0
    return res


def test_feed_from_yaml(xf_case):
    _, _, feed = xf_case
    assert feed.T_in == pytest.approx(793.15)
    assert feed.P_in == pytest.approx(29.0)
    assert feed.F["CH4"] == pytest.approx(5.168)
    assert feed.F["H2O"] / feed.F["CH4"] == pytest.approx(3.358)


def test_higher_alkane_rule_conserves_elements():
    F = {"CH4": 1.0, "H2O": 3.0, "C2H6": 0.1, "C3H8": 0.05, "C4H10": 0.01, "C5H12": 0.005}
    out = r1.convert_higher_alkanes(F)
    def atoms(d):
        c = d.get("CH4", 0) + 2 * d.get("C2H6", 0) + 3 * d.get("C3H8", 0) + 4 * d.get("C4H10", 0) + 5 * d.get("C5H12", 0) + d.get("CO", 0) + d.get("CO2", 0)
        h = 4 * d.get("CH4", 0) + 6 * d.get("C2H6", 0) + 8 * d.get("C3H8", 0) + 10 * d.get("C4H10", 0) + 12 * d.get("C5H12", 0) + 2 * d.get("H2O", 0) + 2 * d.get("H2", 0)
        o = d.get("H2O", 0) + d.get("CO", 0) + 2 * d.get("CO2", 0)
        return c, h, o
    assert atoms(F) == pytest.approx(atoms(out), rel=1e-12)
    assert out["H2"] == 0.0


def test_element_conservation_along_z(xf_result):
    ef = r1.element_flows(xf_result.F)
    for el in ("C", "H", "O", "N"):
        assert np.max(np.abs(ef[el] / ef[el][0] - 1.0)) < 1e-8, el


def test_adiabatic_tube_cools_and_converts_little(xf_case):
    tube, bed, feed = xf_case
    res = r1.simulate(tube, bed, feed, r1.WallBC.constant(feed.T_in), adiabatic=True)
    assert res.success
    assert np.all(np.diff(res.T) <= 1e-9)          # monotone decrease
    assert res.T[-1] < res.T[0]
    assert res.conversion_CH4[-1] < 0.1
    assert np.all(res.q_flux_outer == 0.0)


def test_hot_wall_long_tube_reaches_equilibrium(xf_case):
    tube, bed, feed = xf_case
    res = r1.simulate(tube, bed, feed, r1.WallBC.constant(1200.0), L=2 * tube.L_heated)
    assert res.success
    T_out, P_out = res.T[-1], res.P[-1]
    gas = ct.Solution("gri30.yaml")
    gas.TPX = T_out, P_out * 1e5, {s: res.X[s][-1] for s in r1.SPECIES}
    gas.equilibrate("TP")
    X_eq = gas.mole_fraction_dict()
    for s in r1.SPECIES:
        assert abs(res.X[s][-1] - X_eq.get(s, 0.0)) < 0.03, (s, res.X[s][-1], X_eq.get(s, 0.0))
    assert abs(res.dT_approach_I[-1]) < 30.0


def test_pressure_drop(xf_result):
    assert np.all(np.diff(xf_result.P) < 0.0)
    dP = xf_result.P[0] - xf_result.P[-1]
    assert 0.5 < dP < 5.0


def test_xu_froment_case_outlet(xf_result):
    assert xf_result.success
    assert xf_result.extras["test_wall_time_s"] < 10.0
    out = xf_result.outlet()
    assert 0.5 < out["conversion_CH4"] < 0.9
    assert 1000.0 < out["T_K"] < 1150.0
    assert out["T_wall_inner_K"] < out["T_wall_outer_K"]
    assert out["T_wall_inner_K"] > out["T_K"]


def test_wall_bc_variants():
    c = r1.WallBC.constant(1100.0); assert float(c(3.0)) == 1100.0
    l = r1.WallBC.linear(1000.0, 1200.0, 10.0); assert float(l(5.0)) == pytest.approx(1100.0)
    f = r1.WallBC.from_callable(lambda z: 1000.0 + 10.0 * z); assert float(f(2.0)) == pytest.approx(1020.0)
