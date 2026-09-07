"""The three validation checks on the digitised manufacturer curves (rdt.creep_validate)."""

import numpy as np
import pytest

from rdt import creep
from rdt import creep_validate as cv


# --- check 1: 100 000 h rupture strength ------------------------------------
@pytest.mark.parametrize("key", cv.ALLOYS)
def test_rupture_strength_falls_with_temperature(key):
    rows = cv.rupture_strength_vs_temperature(key)
    sigma = [r["sigma_MPa"] for r in rows]
    assert all(b < a for a, b in zip(sigma, sigma[1:]))
    assert all(1.0 < s < 100.0 for s in sigma)


@pytest.mark.parametrize("key,ref", sorted(cv.REFERENCE_100KH_MPA.items()))
def test_quoted_100kh_strength_is_reproduced_by_the_fitted_curve(key, ref):
    """The temperature at which the fitted minimum curve gives the quoted 100 000 h strength."""
    T = cv.temperature_for_100kh_strength(key, ref)
    # the curve must actually pass through the quoted value at that temperature
    c = creep.load_alloy(key).minimum
    assert float(c.rupture_stress(T + 273.15, cv.DESIGN_LIFE_H)) == pytest.approx(ref, rel=1e-6)
    # and land in a physically sensible band around the design tube-metal temperature
    assert 880.0 < T < 980.0, T


def test_the_two_quoted_strengths_are_matched_at_the_same_temperature():
    """18.3 and 21.2 MPa are a consistent pair: they are matched within a few degrees of each other.

    This is the real content of the check -- it says the two data sheets agree with one another, and
    it locates the temperature the pair was quoted at (about 930 C, not 900 C).
    """
    temps = [cv.temperature_for_100kh_strength(k, v) for k, v in cv.REFERENCE_100KH_MPA.items()]
    assert max(temps) - min(temps) < 5.0, temps
    assert 920.0 < np.mean(temps) < 940.0, temps


def test_at_900C_the_curves_give_more_than_the_quoted_strengths():
    """A design temperature of 900 C is less severe than the point the quoted pair refers to."""
    for key, ref in cv.REFERENCE_100KH_MPA.items():
        c = creep.load_alloy(key).minimum
        assert float(c.rupture_stress(cv.DESIGN_T_C + 273.15, cv.DESIGN_LIFE_H)) > ref


# --- check 2: ingestion regression ------------------------------------------
def test_ingestion_still_recovers_C_and_scatter(tmp_path):
    r = cv.check_ingestion_regression(n_heats=40, n_per_heat=8, tmp_dir=str(tmp_path))
    assert r["C_error"] < 0.10 and r["scatter_error"] < 0.10
    assert r["pass"]


# --- check 3: base case inside the digitised range --------------------------
def test_base_case_is_interpolation_not_extrapolation():
    out = cv.check_base_case_in_range()
    assert set(out) == set(cv.ALLOYS) | {creep.LEGACY_ALLOY}
    for key, v in out.items():
        assert v["inside"], f"{key}: LMP {v['lmp']:.3f} outside {v['lmp_range']}"
        assert v["margin"] > 0.5, f"{key}: only {v['margin']:.3f} of LMP margin"


def test_base_case_lmp_is_consistent_with_the_stress_and_the_constant():
    for key in cv.ALLOYS:
        c = creep.load_alloy(key).minimum
        L = c.lmp_of_stress(cv.BASE_SIGMA_MPA)
        t_r = float(c.time_to_rupture(cv.BASE_T_K, cv.BASE_SIGMA_MPA))
        assert c.lmp(cv.BASE_T_K, t_r) == pytest.approx(L, rel=1e-8)
        assert float(c.rupture_stress(cv.BASE_T_K, t_r)) == pytest.approx(cv.BASE_SIGMA_MPA, rel=1e-6)


def test_main_reports_all_three_checks_passing(capsys):
    assert cv.main([]) == 0
    out = capsys.readouterr().out
    for marker in ("CHECK 1", "CHECK 2", "CHECK 3", "SUMMARY"):
        assert marker in out
    assert "FAIL" not in out.split("SUMMARY")[-1]
