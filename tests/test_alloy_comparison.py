"""Alloy sensitivity of the headline ratios (rdt.alloy_comparison)."""

import numpy as np
import pandas as pd
import pytest

from rdt import alloy_comparison as AC
from rdt import config as C
from rdt import creep

V2_SUMMARY = C.DATA / "scenarios" / "summary_v2.csv"
V2_REGIMES = C.DATA / "optimization" / "regimes_v2.csv"
needs_v2 = pytest.mark.skipif(not (V2_SUMMARY.exists() and V2_REGIMES.exists()), reason="v2 results not present")


@pytest.fixture(scope="module")
def v2_states():
    return AC.load_hourly_states(C.V2)


@needs_v2
def test_rescoring_reproduces_the_published_v2_scenario_ratios(v2_states):
    """The whole method rests on this: re-scoring the stored hourly tube state with the same curve
    must return the numbers the scenario stage itself wrote."""
    got = AC.scenario_life_per_kmol(v2_states, creep.load_alloy(creep.LEGACY_ALLOY).minimum)
    published = pd.read_csv(V2_SUMMARY, index_col=0)["life_consumption_per_kmol_H2_rel_S1"]
    for name, value in got.items():
        assert value == pytest.approx(float(published[name]), rel=1e-10), name


@needs_v2
def test_rescoring_reproduces_the_published_v2_knee_ratio():
    states = AC.regime_states(C.V2)
    got = AC.regime_life_per_kmol(states, creep.load_alloy(creep.LEGACY_ALLOY).minimum)
    published = pd.read_csv(V2_REGIMES, index_col=0)["life_rate_per_kmol_H2_rel_base_phys"]["knee"]
    assert got["base"] == pytest.approx(1.0)
    assert got["knee"] == pytest.approx(float(published), rel=1e-6)


@needs_v2
def test_a_different_alloy_actually_changes_the_scores(v2_states):
    """Guard against the comparison silently re-scoring everything with one curve."""
    legacy = AC.scenario_life_per_kmol(v2_states, creep.load_alloy(creep.LEGACY_ALLOY).minimum)
    g4852 = AC.scenario_life_per_kmol(v2_states, creep.load_alloy("centralloy_g_4852").minimum)
    assert legacy["S1_steady"] == g4852["S1_steady"] == 1.0        # both normalised to S1
    assert legacy["S6_SC2.5"] != pytest.approx(g4852["S6_SC2.5"], rel=1e-6)


@needs_v2
def test_table_covers_every_alloy_and_ratio():
    t = AC.build_table(C.V2)
    assert len(t) == len(AC.ALLOYS)
    assert set(AC.RATIOS) <= set(t.columns)
    assert np.isfinite(t[list(AC.RATIOS)].to_numpy(float)).all()
    assert t["base LMP in range"].all(), "base case extrapolates on some curve"
    assert (t["C"] > 15.0).all() and (t["C"] < 25.0).all()


@needs_v2
def test_ranking_invariance_summary_is_consistent_with_the_table():
    t = AC.build_table(C.V2)
    inv = AC.ranking_is_invariant(t)
    for r in AC.RATIOS:
        v = t[r].to_numpy(float)
        assert inv[r]["min"] == pytest.approx(v.min())
        assert inv[r]["max"] == pytest.approx(v.max())
        assert inv[r]["same_side_of_1"] == bool(np.all(v > 1.0) or np.all(v < 1.0))


def test_ranking_invariance_detects_a_sign_change():
    t = pd.DataFrame({r: [0.9, 1.1, 1.0, 1.0] for r in AC.RATIOS})
    inv = AC.ranking_is_invariant(t)
    assert not inv["all_same_side"]


@needs_v2
def test_latex_table_has_one_row_per_alloy():
    tex = AC.to_latex(AC.build_table(C.V2))
    assert tex.count("\\\\") == len(AC.ALLOYS) + 1          # header plus one row per alloy
    assert "\\toprule" in tex and "\\bottomrule" in tex


@needs_v2
def test_missing_regime_names_the_available_ones():
    with pytest.raises(KeyError, match="available"):
        AC.regime_states(C.V2, regimes=("base", "no_such_regime"))


def test_material_flag_separates_direction_from_magnitude():
    """A ratio hugging 1 that crosses it is 'no material difference', not a contradiction."""
    near = pd.DataFrame({r: [0.99, 1.01, 0.995, 1.005] for r in AC.RATIOS})
    inv = AC.ranking_is_invariant(near, material=0.05)
    assert not inv["all_same_side"]
    assert inv["material_ratios"] == []
    assert inv["material_and_same_side"]          # vacuously: nothing material disagrees

    big = pd.DataFrame({r: [2.0, 2.5, 3.0, 2.2] for r in AC.RATIOS})
    inv = AC.ranking_is_invariant(big, material=0.05)
    assert inv["material_ratios"] == list(AC.RATIOS)
    assert inv["material_and_same_side"]


@needs_v2
def test_material_effects_agree_in_direction_across_alloys():
    """The conclusions the paper actually rests on must not depend on the alloy."""
    inv = AC.ranking_is_invariant(AC.build_table(C.V2))
    assert inv["material_and_same_side"]
    assert {"ageing y4/y1", "knee/base"} <= set(inv["material_ratios"])


@needs_v2
def test_reoptimised_knee_reads_each_runs_own_pareto_result():
    ro = AC.reoptimised_knee(("v2",))
    assert list(ro.index) == ["v2"]
    published = pd.read_csv(V2_REGIMES, index_col=0)["life_rate_per_kmol_H2_rel_base_phys"]["knee"]
    assert ro["knee/base (own curve)"]["v2"] == pytest.approx(float(published))


def test_reoptimised_knee_skips_versions_without_results():
    assert len(AC.reoptimised_knee(("v1", "v2", "v3"))) <= 3
