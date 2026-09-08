"""Per-tube population extension of the life model (rdt.population).

Two things need guarding: that the convexity correction the module reports is the right number (checked
against the closed-form lognormal result, which the code does not use), and that a fixed seed gives a
fixed bundle, since every number in the analysis is quoted for one draw.
"""

import numpy as np
import pandas as pd
import pytest

from rdt import config as C
from rdt import creep
from rdt import population as pop

V3_REGIMES = C.DATA / "optimization" / "regimes_v3.csv"
V3_HOURLY = C.DATA / "scenarios" / "hourly_v3_S1_steady.csv.gz"
needs_v3 = pytest.mark.skipif(not (V3_REGIMES.exists() and V3_HOURLY.exists()), reason="v3 results not present")

#: the calibrated base hot spot: 1146.7 K at 12.87 MPa on the G 4852 lower scatter band
T_BASE, SIGMA_BASE = 1146.7, 12.87


@pytest.fixture(scope="module")
def curve():
    return creep.load_alloy("centralloy_g_4852").minimum


@pytest.fixture(scope="module")
def state():
    return pop.TubeState(label="test", T_wo_max_K=T_BASE, sigma_hot_MPa=SIGMA_BASE,
                         H2_net_kmol_h=14.8, evaluated_by="fixture")


# ---------------------------------------------------------------------------
# The analytic cross-check
# ---------------------------------------------------------------------------
def test_local_slope_matches_a_finite_difference_of_the_curve(curve):
    """``a = LMP/(scale T^2)`` is derived by hand; confirm it against the curve it claims to describe."""
    a = pop.local_slope(curve, T_BASE, SIGMA_BASE)
    h = 0.05
    lo = np.log10(float(curve.time_to_rupture(T_BASE - h, SIGMA_BASE)))
    hi = np.log10(float(curve.time_to_rupture(T_BASE + h, SIGMA_BASE)))
    assert a == pytest.approx(-(hi - lo) / (2 * h), rel=1e-8)
    assert 0.01 < a < 0.05, "slope outside any physically sensible range for a reformer alloy"


def test_analytic_mean_ratio_is_the_lognormal_mean():
    """exp(s^2/2) with s = a ln10 sigma, written out independently of the module."""
    a, sg = 0.0224, 18.5
    s = a * np.log(10.0) * sg
    assert pop.analytic_mean_ratio(a, sg) == pytest.approx(float(np.exp(0.5 * s * s)), rel=1e-12)
    assert pop.analytic_mean_ratio(a, 0.0) == pytest.approx(1.0), "no scatter, no convexity correction"


def test_sampled_population_mean_converges_to_the_lognormal_result_for_small_scatter(curve, state):
    """Over a narrow band the slope is effectively constant, so the sampled ratio must hit the
    closed-form value to within Monte Carlo error and nothing else."""
    sg, n = 4.0, 400_000
    a = pop.local_slope(curve, state.T_wo_max_K, state.sigma_hot_MPa)
    rate_avg = 1.0 / float(curve.time_to_rupture(state.T_wo_max_K, state.sigma_hot_MPa))
    r = pop.damage_rates(curve, state.T_wo_max_K, state.sigma_hot_MPa, pop.offsets(sg, n, seed=3))
    assert float(np.mean(r)) / rate_avg == pytest.approx(pop.analytic_mean_ratio(a, sg), rel=3e-3)


@pytest.mark.parametrize("sg", pop.SIGMA_GRID_K)
def test_validate_analytic_agrees_with_the_lognormal_result_at_the_reported_scatter(curve, state, sg):
    """At 15-22 K the curvature of ``a`` is no longer negligible. The sampled mean must still land
    close to the analytic value, and must land *below* it: ``a`` falls as 1/T^2, so the hot tail gains
    less than the constant-slope approximation assumes."""
    v = pop.validate_analytic(curve, state, sg, n_large=400_000, seed=3)
    assert v["analytic"] == pytest.approx(pop.analytic_mean_ratio(v["a_dec_per_K"], sg), rel=1e-12)
    assert v["sampled_large"] == pytest.approx(v["analytic"], rel=0.06)
    assert -8.0 < v["rel_err_large_pct"] < 0.0, "curvature term has the wrong sign or size"
    # the 336-tube draw is noisier than the large one, and by roughly the predicted amount
    assert abs(v["rel_err_N336_pct"]) < 2.0 * v["mc_se_N336_pct"] + abs(v["rel_err_large_pct"])
    assert v["mc_se_N336_pct"] == pytest.approx(100 * v["lognormal_cv"] / np.sqrt(pop.N_TUBES))


def test_population_mean_exceeds_the_average_tube_and_grows_with_scatter(curve, state):
    """The motivating inequality: convexity makes the mean damage rate exceed the damage at the mean."""
    dT = [pop.offsets(sg, 200_000, seed=1) for sg in (0.0, 10.0, 20.0)]
    ratios = []
    for d, sg in zip(dT, (0.0, 10.0, 20.0)):
        s = pop.population_stats(curve, state, d, sg)
        ratios.append(s["ratio_pop_mean_to_avg"])
    assert ratios[0] == pytest.approx(1.0), "zero scatter must leave the average tube alone"
    assert ratios[0] < ratios[1] < ratios[2]


# ---------------------------------------------------------------------------
# Sampling reproducibility
# ---------------------------------------------------------------------------
def test_offsets_are_identical_for_a_fixed_seed():
    a = pop.offsets(18.5, 336, seed=0)
    b = pop.offsets(18.5, 336, seed=0)
    assert np.array_equal(a, b), "same seed must give a bit-identical bundle"
    assert not np.array_equal(a, pop.offsets(18.5, 336, seed=1))
    assert not np.array_equal(a, pop.offsets(22.0, 336, seed=0))


def test_offsets_have_the_requested_shape_and_scale():
    d = pop.offsets(18.5, pop.N_TUBES, seed=0)
    assert d.shape == (336,) and pop.N_TUBES == pop.ROWS * pop.TUBES_PER_ROW == 336
    big = pop.offsets(18.5, 500_000, seed=0)
    assert float(np.std(big, ddof=1)) == pytest.approx(18.5, rel=1e-2)
    assert abs(float(np.mean(big))) < 0.2
    # the draw is deliberately NOT re-centred: one real bundle misses zero
    assert abs(float(np.mean(d))) > 1e-6


def test_offsets_reject_impossible_arguments():
    with pytest.raises(ValueError):
        pop.offsets(18.5, 0)
    with pytest.raises(ValueError):
        pop.offsets(-1.0, 336)


def test_the_whole_statistic_is_reproducible_under_a_fixed_seed(curve, state):
    """Every number the report quotes for one draw must come back the same on a re-run."""
    first = pop.population_stats(curve, state, pop.offsets(18.5, 336, seed=0), 18.5)
    second = pop.population_stats(curve, state, pop.offsets(18.5, 336, seed=0), 18.5)
    numeric = [k for k, v in first.items() if isinstance(v, float)]
    assert numeric, "nothing numeric to compare"
    for k in numeric:
        assert first[k] == second[k], k
    other = pop.population_stats(curve, state, pop.offsets(18.5, 336, seed=7), 18.5)
    assert other["mult_hottest"] != first["mult_hottest"]


# ---------------------------------------------------------------------------
# The damage-rate kernel
# ---------------------------------------------------------------------------
def test_zero_offset_reproduces_the_average_tube_exactly(curve, state):
    r = pop.damage_rates(curve, state.T_wo_max_K, state.sigma_hot_MPa, np.zeros(5))
    expected = 1.0 / float(curve.time_to_rupture(state.T_wo_max_K, state.sigma_hot_MPa))
    assert r.shape == (5,)
    assert np.allclose(r, expected, rtol=1e-12)


def test_damage_rates_broadcast_over_an_hourly_history(curve):
    T = np.array([1140.0, 1150.0, 1145.0]); s = np.full(3, SIGMA_BASE)
    dT = np.array([-10.0, 0.0, 10.0])
    r = pop.damage_rates(curve, T, s, dT)
    assert r.shape == (3, 3)
    for i in range(3):                      # each row is the single-point result of that hour
        assert np.allclose(r[i], pop.damage_rates(curve, float(T[i]), float(s[i]), dT), rtol=1e-12)
    assert (np.diff(r, axis=1) > 0).all(), "a hotter tube must consume life faster"


def test_damage_rates_reject_mismatched_histories(curve):
    with pytest.raises(ValueError):
        pop.damage_rates(curve, np.zeros(3) + 1140.0, np.full(2, SIGMA_BASE), np.zeros(4))


def test_percentiles_and_the_decile_are_ordered_as_claimed(curve, state):
    s = pop.population_stats(curve, state, pop.offsets(18.5, 336, seed=0), 18.5)
    assert s["mult_p1"] < s["mult_p5"] < 1.0 < s["ratio_pop_mean_to_avg"] < s["mult_p95"]
    assert s["mult_p95"] < s["mult_hottest_decile"] < s["mult_hottest"]
    assert s["n_hottest_decile"] == 34                      # 10 % of 336
    for k in ("hottest", "hottest_decile", "p95", "p5", "p1"):
        assert s[f"life_{k}_y"] == pytest.approx(s["life_avg_tube_y"] / s[f"mult_{k}"], rel=1e-12)


# ---------------------------------------------------------------------------
# Wiring to the stored results
# ---------------------------------------------------------------------------
@needs_v3
def test_the_steady_scenario_anchor_matches_its_stored_history():
    st = pop.steady_scenario_state(C.V3, "S1_steady")
    df = pop.load_hourly(C.V3, "S1_steady")
    assert st.T_wo_max_K == pytest.approx(float(df.T_wo_max_K.iloc[0]))
    assert "surrogate" in st.evaluated_by, "the code path must be reported, not assumed"


@needs_v3
def test_a_non_steady_scenario_is_refused_as_a_single_anchor():
    with pytest.raises(ValueError, match="not steady"):
        pop.steady_scenario_state(C.V3, "S2_daily")


@needs_v3
def test_scenario_population_reduces_to_the_published_damage_at_zero_scatter():
    """With no offsets the per-tube sum must be the annual damage the scenario stage itself wrote."""
    curve = creep.load_alloy(C.V3.creep_alloy).minimum
    got = pop.scenario_population(C.V3, "S1_steady", np.zeros(4), curve)
    published = pd.read_csv(C.V3.scenarios_summary, index_col=0)["annual_damage_D"]["S1_steady"]
    assert got["D_avg_tube"] == pytest.approx(float(published), rel=1e-9)
    assert got["D_pop_mean"] == pytest.approx(got["D_avg_tube"], rel=1e-12)
    assert got["ratio_pop_mean_to_avg"] == pytest.approx(1.0)


@needs_v3
def test_a_static_offset_is_held_over_the_year_not_resampled():
    """A tube keeps its position: for a steady scenario the annual damage of tube i must be exactly
    8760 h times its hourly rate."""
    curve = creep.load_alloy(C.V3.creep_alloy).minimum
    dT = pop.offsets(18.5, 8, seed=0)
    got = pop.scenario_population(C.V3, "S1_steady", dT, curve)
    st = pop.steady_scenario_state(C.V3, "S1_steady")
    hourly = pop.damage_rates(curve, st.T_wo_max_K, st.sigma_hot_MPa, dT)
    assert got["D_hottest"] == pytest.approx(8760.0 * float(np.max(hourly)), rel=1e-9)


def test_ranking_flags_only_rows_that_actually_move():
    t = pd.DataFrame({"cost_avg": [1.0, 2.0, 3.0], "cost_pop_mean": [1.0, 2.0, 3.0],
                      "cost_hottest_decile": [2.0, 1.0, 3.0]}, index=["a", "b", "c"])
    rk = pop.ranking(t)
    assert not rk.rank_changes_pop.any()
    assert list(rk.rank_changes_decile) == [True, True, False]
    assert pop.swapped_pairs(t, "cost_avg", "cost_pop_mean") == []
    swaps = pop.swapped_pairs(t, "cost_avg", "cost_hottest_decile")
    assert [(i, j) for i, j, _, _ in swaps] == [("a", "b")]
    assert swaps[0][2] == pytest.approx(100.0)              # 1.0 vs 2.0 on the average tube
