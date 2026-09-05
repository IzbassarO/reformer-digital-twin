"""Sanity checks on data/literature_validation/latham2008_plant_cases.csv (Latham 2008 thesis Appendix H)."""

import numpy as np
import pytest

from rdt import latham_cases as lc


@pytest.fixture(scope="module")
def df():
    return lc.load_cases()


def test_cases_present(df):
    assert list(df.case) == ["Plant_A", "Plant_B", "Plant_C1", "Plant_C2"]
    assert list(df.plant_rate_pct) == [99, 92, 96, 71]


def test_inlet_ranges(df):
    assert df.T_in_K.between(850, 900).all()
    assert (df.P_in_Pa / 1e5).between(25, 35).all()
    assert df.steam_to_carbon_derived.between(2.5, 3.5).all()
    assert df.feed_per_tube_kmol_h.between(10, 40).all()


def test_tube_wall_temperature_ranges(df):
    assert df.TWT_upper_K.between(1050, 1250).all()
    assert df.TWT_lower_K.between(1050, 1250).all()
    assert (df.TWT_lower_K > df.TWT_upper_K).all()
    assert np.allclose(df.TWT_upper_frac, 3.66 / 12.5) and np.allclose(df.TWT_lower_frac, 8.53 / 12.5)


def test_compositions_sum_to_100(df):
    feed_cols = [c for c in df.columns if c.startswith("feed_x_")]
    assert np.allclose(df[feed_cols].sum(axis=1), 100.0, atol=0.2)
    out_cols = [c for c in df.columns if c.startswith("out_x_")]
    assert np.allclose(df[out_cols].sum(axis=1), 100.0, atol=0.2)


def test_plant_b_wet_h2_anchor(df):
    assert df.loc[df.case == "Plant_B", "out_x_H2_wet_molpct"].item() == pytest.approx(46.3, abs=0.05)


def test_given_furnace_quantities_are_empty(df):
    for c in ("heat_of_combustion_W_given", "flue_gas_kmol_h_given", "excess_air_pct_given"):
        assert df[c].isna().all(), c
    assert (df.heat_of_combustion_LHV_W_derived / 1e6).between(100, 300).all()
    assert df.excess_air_pct_derived.between(5, 40).all()
