"""Tests for rdt.furnace (long-furnace model coupled to the tube model)."""

import time

import numpy as np
import pytest

from rdt import furnace as fu
from rdt import latham_cases as lc
from rdt import reactor1d as r1


@pytest.fixture(scope="module")
def plant_a():
    row = lc.load_cases().iloc[0]
    t0 = time.perf_counter()
    res = lc.run_case(row)
    res.extras["test_wall_time_s"] = time.perf_counter() - t0
    return row, res


def test_heat_release_parabola_normalisation():
    L = 12.5
    rel = fu.HeatRelease(L_q=0.48, alpha_top=0.182, f_loss=0.02, n_sections=15)
    a, b, c, zq = rel.coefficients(L)
    assert zq == pytest.approx(0.48 * L)
    assert a < 0.0
    z = np.linspace(0.0, L, 20001)
    r = rel.density(z, 1.0e6, L)
    assert np.trapezoid(r, z) == pytest.approx(0.98e6, rel=1e-4)
    top = z <= L / 15
    assert np.trapezoid(r[top], z[top]) == pytest.approx(0.182 * 0.98e6, rel=1e-3)
    assert np.all(r >= -1e-9) and np.all(r[z > zq] == 0.0)
    assert rel.density(zq, 1.0e6, L) == pytest.approx(0.0, abs=1e-6)


def test_heat_release_rejects_bad_shape():
    with pytest.raises(ValueError):
        fu.HeatRelease(L_q=0.48, alpha_top=0.9).coefficients(12.5)


def test_flue_gas_from_combustion_element_balance():
    # 1 kmol/h CH4 + 12 kmol/h air-like mixture (21 % O2, 79 % N2): 2 O2 stoichiometric, 26 % excess
    fg = fu.FlueGas.from_combustion([(1.0, {"CH4": 1.0}, 300.0), (12.0, {"O2": 0.21, "N2": 0.79}, 300.0)])
    n = fg.n_kmol_h
    assert n * fg.X["CO2"] == pytest.approx(1.0)
    assert n * fg.X["H2O"] == pytest.approx(2.0)
    assert n * fg.X["O2"] == pytest.approx(12 * 0.21 - 2.0)
    assert n * fg.X["N2"] == pytest.approx(12 * 0.79)
    assert fg.excess_air_pct == pytest.approx(26.0)
    assert fg.Q_comb_W == pytest.approx(802.6e3 * 1e3 / 3600.0, rel=2e-3)   # LHV of CH4 802.6 kJ/mol
    assert fg.T_in_K == pytest.approx(300.0, abs=0.5)


def test_flue_gas_from_composition():
    fg = fu.FlueGas.from_composition(12350.0, {"CO2": 0.08, "H2O": 0.17, "N2": 0.73, "O2": 0.02}, 600.0, 202.9e6, 9.2)
    assert sum(fg.X.values()) == pytest.approx(1.0)
    assert fg.Q_comb_W == 202.9e6


def test_energy_conservation(plant_a):
    _, res = plant_a
    e = res.energy
    assert abs(e["closure_rel_error"]) < 0.005
    assert e["release_integral_W"] == pytest.approx(e["Q_eff_W"], rel=1e-3)
    assert 0.3 < e["duty_fraction_of_Q_comb"] < 0.7


def test_furnace_gas_cools_after_release_zone(plant_a):
    _, res = plant_a
    zq = res.extras["L_q"] * res.z[-1]
    after = res.z > zq
    assert np.all(np.diff(res.T_fg[after]) < 0.0)


def test_temperature_ordering(plant_a):
    """T_gas < T_wi < T_wo < T_fg. The plug-flow furnace gas enters at the mixed fuel/air temperature
    (about 525 K, colder than the 885 K process gas) and is heated by the release parabola within the
    top section, so the ordering is checked from the end of the top section (z >= L/15) onward."""
    _, res = plant_a
    z_top = res.z[-1] / 15
    m = res.z >= z_top
    assert np.all(res.tube.T[m] < res.T_wi[m])
    assert np.all(res.T_wi[m] < res.T_wo[m])
    assert np.all(res.T_wo[m] < res.T_fg[m])
    # above the top section the furnace gas can be colder than the tube: heat then flows outward and the
    # ordering is reversed there (T_fg < T_wo < T_wi < T_gas); check that it is at least consistent
    n_rev = np.sum(res.T_fg[~m] < res.tube.T[~m])
    assert n_rev >= 1 and np.all(res.T_wo[~m][res.T_fg[~m] < res.tube.T[~m]] < res.tube.T[~m][res.T_fg[~m] < res.tube.T[~m]])


def test_coupled_run_time_and_success(plant_a):
    _, res = plant_a
    assert res.success
    assert res.extras["test_wall_time_s"] < 30.0
    assert res.tube.extras["coupled"] is True


def test_wall_temperature_root_brackets():
    tube = r1.TubeGeometry(d_i=0.116, d_o=0.146, L_heated=12.5, lambda_tube=29.6)
    T_wo = fu.solve_wall_temperature(1400.0, 900.0, 800.0, tube, 5.0, 0.35)
    assert 900.0 < T_wo < 1400.0
    q = fu.furnace_flux(1400.0, T_wo, 5.0, 0.35)
    assert q == pytest.approx(800.0 * (T_wo - 900.0) * tube.d_i / tube.d_o, rel=1e-6)
    # colder furnace gas than process gas: wall between the two
    T_wo2 = fu.solve_wall_temperature(600.0, 900.0, 800.0, tube, 5.0, 0.35)
    assert 600.0 < T_wo2 < 900.0
