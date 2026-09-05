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


# ---------------------------------------------------------------------------
# Xu & Froment (Part II) heat-transfer chain and inlet definitions
# ---------------------------------------------------------------------------
def test_xu_froment_alpha_i_recomputed_by_hand():
    """Recompute Part II Eq. 12 and its correlations by hand, going through kJ/(m h K)."""
    lam_g, mu, cp, G_s, d_p, d_ti, lam_er0 = 0.10, 3.0e-5, 2500.0, 15.0, 0.010, 0.10, 1.0
    out = r1.xu_froment_alpha_i(lam_g, mu, cp, G_s, d_p, d_ti, lam_er0)
    Re = d_p * G_s / mu                      # 5000
    Pr = cp * mu / lam_g                     # 0.75
    assert Re == pytest.approx(5000.0) and Pr == pytest.approx(0.75)
    # original units: lambda in kJ/(m h K) = W/(m K) * 3.6, alpha in kJ/(m2 h K)
    lam_g_kj, lam_er0_kj = lam_g * 3.6, lam_er0 * 3.6
    alpha_w0_kj = 8.694 * lam_er0_kj / d_ti ** (4.0 / 3.0)            # 674.3 kJ/(m2 h K)
    alpha_w_kj = alpha_w0_kj + 0.444 * Re * Pr * lam_g_kj / d_p        # + 59940
    lam_er_kj = lam_er0_kj + 0.14 * lam_g_kj * Re * Pr                 # 3.6 + 189 = 192.6
    alpha_i_kj = 8 * lam_er_kj * alpha_w_kj / (8 * lam_er_kj + alpha_w_kj * d_ti)
    assert out["alpha_w0"] == pytest.approx(alpha_w0_kj / 3.6, rel=1e-12)
    assert out["alpha_w0"] == pytest.approx(187.31, rel=1e-3)
    assert out["lam_er"] == pytest.approx(lam_er_kj / 3.6, rel=1e-12)
    assert out["lam_er"] == pytest.approx(53.5, rel=1e-6)
    assert out["alpha_i"] == pytest.approx(alpha_i_kj / 3.6, rel=1e-12)
    assert out["alpha_i"] == pytest.approx(3412.5, rel=1e-3)


def test_kunii_smith_phi_and_limits():
    """phi_1, phi_2 at kappa = 10 match hand evaluation of the Kunii-Smith expression; static
    conductivity exceeds the pure-gas void contribution and grows with temperature (radiation)."""
    import math
    kappa = 10.0
    def phi_i(n):
        s2 = 1 / n; c = math.sqrt(1 - s2); a = (kappa - 1) / kappa
        return 0.5 * a * a * s2 / (math.log(kappa - (kappa - 1) * c) - a * (1 - c)) - 2 / (3 * kappa)
    assert phi_i(1.5) == pytest.approx(0.160, abs=0.002)
    assert phi_i(4 * math.sqrt(3)) == pytest.approx(0.0635, abs=0.001)
    lam_g, lam_s, eps, d_p = 0.1, 1.0, 0.526, 9.44e-3
    l800 = r1.kunii_smith_lambda_er0(lam_g, lam_s, eps, 800.0, 0.8, d_p)
    l1100 = r1.kunii_smith_lambda_er0(lam_g, lam_s, eps, 1100.0, 0.8, d_p)
    assert l800 > eps * lam_g and l1100 > l800


def test_inlet_definitions_xu_froment_feed():
    f_none = r1.feed_from_xu_froment()
    f_lat = r1.feed_from_xu_froment(split_alkanes=True, inlet_higher_alkanes="latham")
    f_xf = r1.feed_from_xu_froment(split_alkanes=True, inlet_higher_alkanes="xu_froment")
    for f in (f_none, f_lat, f_xf):
        assert f.F_CH4_equivalent == pytest.approx(5.168)
    def inlet_conv(f):
        return (f.F_CH4_equivalent - f.state_flows()["CH4"]) / f.F_CH4_equivalent
    assert inlet_conv(f_none) == pytest.approx(0.0)
    assert inlet_conv(f_lat) == pytest.approx(0.0175, abs=5e-4)
    assert inlet_conv(f_xf) == pytest.approx(0.0914, abs=5e-4)   # C2+ carbon fraction of the Table 2 gas
    # both rules conserve C, H, O
    for f in (f_lat, f_xf):
        raw = f.F; st = f.state_flows()
        C_raw = sum(k * raw.get(s, 0) for s, k in {"CH4": 1, "C2H6": 2, "C3H8": 3, "C4H10": 4, "C5H12": 5, "CO": 1, "CO2": 1}.items())
        H_raw = sum(k * raw.get(s, 0) for s, k in {"CH4": 4, "C2H6": 6, "C3H8": 8, "C4H10": 10, "C5H12": 12, "H2O": 2, "H2": 2}.items())
        ef = r1.element_flows({s: np.array([st[s]]) for s in r1.SPECIES})
        assert float(ef["C"][0]) == pytest.approx(C_raw, rel=1e-12)
        assert float(ef["H"][0]) == pytest.approx(H_raw, rel=1e-12)


def test_heat_transfer_option_runs(xf_case):
    tube, bed, feed = xf_case
    wall = r1.WallBC.linear(1000.0, 1150.0, tube.L_heated)
    res = r1.simulate(tube, bed, feed, wall, heat_transfer="xu_froment")
    assert res.success and res.extras["heat_transfer"] == "xu_froment"
    assert np.all(res.alpha_i > 0)
    with pytest.raises(ValueError):
        r1.simulate(tube, bed, feed, wall, heat_transfer="nonsense")
