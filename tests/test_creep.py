"""Tests for rdt.creep."""

import numpy as np
import pytest

from rdt import creep


@pytest.fixture(scope="module")
def curve():
    return creep.LarsonMillerCurve.from_yaml()


def test_hoop_stress_variants():
    s = creep.hoop_stress(28.0, 0.146, 0.015)
    assert s == pytest.approx(12.2, abs=0.3)                       # mean-diameter, API 530 style
    assert creep.hoop_stress_thin(28.0, 0.146, 0.015) == pytest.approx(13.63, abs=0.05)
    lame_in = creep.hoop_stress_lame(28.0, 0.146, 0.015, "inner"); lame_out = creep.hoop_stress_lame(28.0, 0.146, 0.015, "outer")
    assert lame_in > s > lame_out
    # wall-averaged Lame hoop stress equals P r_i / t by force balance (10.83 MPa here), which is
    # 11 % below the API 530 mean-diameter value P (d_o - t)/(2 t); the two are different conventions
    assert creep.hoop_stress_lame(28.0, 0.146, 0.015, "mean") == pytest.approx(2.8 * 0.058 / 0.015, rel=1e-3)


def test_lmp_round_trip(curve):
    for T, s in ((1150.0, 10.0), (1200.0, 12.0), (1250.0, 6.0)):
        t = float(curve.time_to_rupture(T, s))
        s_back = float(curve.rupture_stress(T, t))
        assert abs(s_back / s - 1.0) < 1e-6
        T_back = curve.temperature_for_life(s, t)
        assert abs(T_back / T - 1.0) < 1e-6


def test_master_curve_monotone_and_in_range(curve):
    L = np.linspace(*curve.lmp_range, 50)
    assert np.all(np.diff(curve.log10_stress(L)) < 0)
    assert curve.C == pytest.approx(22.96)


def test_yeh_design_anchor(curve):
    """Yeh (2021) p. 10-11: 925 degC, 11.9 MPa -> about 5e6 h on the manufacturer curve."""
    t = float(curve.time_to_rupture(925.0 + 273.15, 11.9))
    assert 2.5e6 < t < 1.0e7, t


def test_life_halving_temperature(curve):
    dT = curve.life_halving_dT(1150.0, 12.2)
    assert 8.0 < dT < 25.0, dT
    t1 = float(curve.time_to_rupture(1150.0, 12.2)); t2 = float(curve.time_to_rupture(1150.0 + dT, 12.2))
    assert t2 / t1 == pytest.approx(0.5, rel=0.03)


def test_robinson_two_segments(curve):
    hist = [(20000.0, 1150.0, 12.0), (5000.0, 1180.0, 12.0)]
    d = creep.robinson_damage(hist, curve)
    f1 = 20000.0 / float(curve.time_to_rupture(1150.0, 12.0)); f2 = 5000.0 / float(curve.time_to_rupture(1180.0, 12.0))
    assert d["damage"] == pytest.approx(f1 + f2, rel=1e-12)
    rem = creep.remaining_life(hist, 1150.0, 12.0, curve)
    assert rem["remaining_h"] == pytest.approx((1 - f1 - f2) * float(curve.time_to_rupture(1150.0, 12.0)), rel=1e-12)
