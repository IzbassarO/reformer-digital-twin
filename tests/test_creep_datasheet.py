"""Master curves and scatter model built from the digitised manufacturer data sheets (rdt.creep)."""

import numpy as np
import pandas as pd
import pytest

from rdt import config as C
from rdt import creep

ALLOYS = ("centralloy_g_4852", "centralloy_g_4852_micro", "centralloy_et_45_micro")
SHEET_C = {"centralloy_g_4852": 18.6, "centralloy_g_4852_micro": 22.9, "centralloy_et_45_micro": 19.3}


@pytest.fixture(autouse=True)
def _restore_active_alloy(monkeypatch):
    """use_config writes a process-wide env var; keep it out of the other tests."""
    monkeypatch.delenv(creep.ALLOY_ENV, raising=False)
    monkeypatch.setattr(creep, "_ACTIVE_ALLOY", creep.DEFAULT_ALLOY)


def test_available_alloys_lists_the_three_sheets_and_the_legacy_placeholder():
    keys = set(creep.available_alloys())
    assert set(ALLOYS) <= keys and creep.LEGACY_ALLOY in keys


@pytest.mark.parametrize("key", ALLOYS)
def test_each_alloy_uses_the_constant_from_its_own_sheet(key):
    """The three sheets print different constants; 22.9 belongs to G 4852 Micro alone."""
    a = creep.load_alloy(key)
    assert a.C == SHEET_C[key]
    assert a.average.C == a.C and a.minimum.C == a.C


@pytest.mark.parametrize("key", ALLOYS)
def test_cubic_master_curve_fits_the_digitised_points(key):
    a = creep.load_alloy(key)
    df = pd.read_csv(creep.DATA_DIR / f"{key}_rupture_curve.csv")
    for kind in ("average", "minimum"):
        g = df[df.curve == kind]
        c = a.curve(kind)
        assert c.degree == 3
        res = np.log10(g.stress_mpa.to_numpy()) - c.log10_stress(g.lmp.to_numpy())
        assert np.sqrt((res**2).mean()) < 0.01, f"{key}/{kind} RMS {np.sqrt((res ** 2).mean()):.4f}"
        # the master curve must fall monotonically over the digitised range
        L = np.linspace(*c.lmp_range, 200)
        assert np.all(np.diff(c.log10_stress(L)) < 0)


@pytest.mark.parametrize("key", ALLOYS)
def test_average_curve_is_above_the_minimum_curve(key):
    a = creep.load_alloy(key)
    L = np.linspace(*a.lmp_range, 100)
    assert np.all(a.average.log10_stress(L) > a.minimum.log10_stress(L))


@pytest.mark.parametrize("key", ALLOYS)
def test_scatter_is_the_horizontal_gap_and_scales_as_one_over_T(key):
    a = creep.load_alloy(key)
    s = a.scatter
    assert not s.assumed and np.all(s.delta_lmp > 0)
    # the LMP gap is a property of the curve pair alone; decades of life scale as 1/T
    assert np.allclose(s.decades(1000.0), s.decades(2000.0) * 2.0)
    # and it really is the horizontal separation at constant stress
    sigma = float(s.stress_MPa[len(s.stress_MPa) // 2])
    gap = a.average.lmp_of_stress(sigma) - a.minimum.lmp_of_stress(sigma)
    assert s.delta_lmp_at(a.average.lmp_of_stress(sigma)) == pytest.approx(gap, rel=1e-6)


@pytest.mark.parametrize("key", ALLOYS)
def test_scatter_decades_are_physically_plausible(key):
    """A factor of about two in rupture life, i.e. a few tenths of a decade."""
    s = creep.load_alloy(key).scatter.summary(1147.0)
    assert 0.1 < s["min_decades"] and s["max_decades"] < 1.0
    assert s["min_decades"] <= s["median_decades"] <= s["max_decades"]


def test_scatter_gap_equals_the_life_ratio_of_the_two_curves():
    """delta_log10_tr must reproduce the actual ratio of rupture times at the same T and stress."""
    a = creep.load_alloy("centralloy_g_4852")
    T = 1147.0
    for sigma in (10.0, 12.9, 20.0):
        t_avg = float(a.average.time_to_rupture(T, sigma))
        t_min = float(a.minimum.time_to_rupture(T, sigma))
        want = np.log10(t_avg / t_min)
        got = a.scatter.decades(T, a.average.lmp_of_stress(sigma))
        assert got == pytest.approx(want, rel=1e-4), (sigma, got, want)


def test_legacy_placeholder_still_selectable_and_unchanged():
    """v1/v2 reproducibility: the Yeh placeholder keeps its constant and its assumed 0.3-decade scatter."""
    a = creep.load_alloy(creep.LEGACY_ALLOY)
    assert a.C == pytest.approx(22.96)
    assert a.scatter.assumed
    assert a.scatter.summary(1150.0)["median_decades"] == pytest.approx(creep.LEGACY_SCATTER_DECADES, rel=1e-9)
    # identical to what the v2 code path loaded
    ref = creep.LarsonMillerCurve.from_yaml()
    assert np.allclose(a.minimum.coeffs, ref.coeffs)


def test_use_config_selects_the_alloy_and_exports_it_for_worker_processes():
    creep.use_config(C.V2)
    assert creep.alloy_key() == creep.LEGACY_ALLOY
    assert creep.active_curve().C == pytest.approx(22.96)
    creep.use_config(C.V3)
    assert creep.alloy_key() == "centralloy_g_4852"
    assert creep.active_curve().C == 18.6
    # joblib workers re-import the module, so the selection must travel in the environment
    import os

    assert os.environ[creep.ALLOY_ENV] == "centralloy_g_4852"


def test_v1_and_v2_pin_the_legacy_curve_v3_uses_the_data_sheet():
    assert C.V1.creep_alloy == creep.LEGACY_ALLOY and C.V2.creep_alloy == creep.LEGACY_ALLOY
    assert C.V3.creep_alloy == "centralloy_g_4852"


@pytest.mark.parametrize("key", ALLOYS)
def test_round_trip_time_and_stress_on_a_data_sheet_curve(key):
    c = creep.load_alloy(key).minimum
    L = float(np.mean(c.lmp_range))
    sigma = float(10 ** c.log10_stress(L))
    for T in (1100.0, 1147.0, 1200.0):
        t = float(c.time_to_rupture(T, sigma))
        assert float(c.rupture_stress(T, t)) == pytest.approx(sigma, rel=1e-6)
        assert c.lmp(T, t) == pytest.approx(L, rel=1e-6)


def test_unknown_alloy_key_names_the_available_ones():
    with pytest.raises(KeyError, match="available"):
        creep.load_alloy("no_such_alloy")
