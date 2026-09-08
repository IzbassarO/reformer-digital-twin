"""Thermal-fatigue screening of the RQ2 scenarios (rdt.fatigue).

The rainflow counter is the part of this module that could be wrong without anyone noticing, so it is
tested against synthetic histories whose cycle content can be worked out by hand. Each test states the
hand count in its docstring; none of the expected values is taken from another implementation.
"""

import numpy as np
import pandas as pd
import pytest

from rdt import config as C
from rdt import fatigue as F


def counts(cycles):
    """(range, count) pairs, rounded and sorted, for comparison with a hand count."""
    agg = {}
    for rng, _mean, cnt in cycles:
        agg[round(rng, 9)] = agg.get(round(rng, 9), 0.0) + cnt
    return sorted(agg.items())


# ---------------------------------------------------------------------------
# Reversal extraction
# ---------------------------------------------------------------------------
def test_reversals_keep_the_end_points_and_drop_the_interior_of_a_ramp():
    assert list(F.reversals([0, 1, 2, 3])) == [0, 3]
    assert list(F.reversals([3, 2, 1, 0])) == [3, 0]
    assert list(F.reversals([0, 2, 1, 3])) == [0, 2, 1, 3]


def test_reversals_collapse_repeated_values():
    """A history that never moves has no reversals at all, and therefore no cycles."""
    assert list(F.reversals([7, 7, 7, 7])) == [7]
    assert F.rainflow([7.0] * 8760) == []
    assert list(F.reversals([0, 0, 5, 5, 0, 0])) == [0, 5, 0]


# ---------------------------------------------------------------------------
# Rainflow, against hand counts
# ---------------------------------------------------------------------------
def test_constant_amplitude_gives_one_cycle_per_up_down_pair():
    """0-5-0-5-0-5-0 rises and falls three times, so it is three full cycles of range 5.

    The four-point rule closes two of them and leaves the residue 0-5-0, which is two half cycles of
    range 5 -- one more full cycle. Three in total, all of range 5, mean 2.5.
    """
    c = F.rainflow([0, 5, 0, 5, 0, 5, 0])
    assert counts(c) == [(5.0, 3.0)]
    assert all(m == pytest.approx(2.5) for _, m, _ in c)


def test_nested_history_recovers_each_nested_excursion():
    """0-10-4-6-4-10-0 has three nested excursions, by inspection:

    the innermost 4-6-4 (range 2), the 10-4-10 that contains it (range 6), and the overall 0-10-0
    (range 10). The first two close; the outermost is left in the residue as two half cycles.
    """
    c = F.rainflow([0, 10, 4, 6, 4, 10, 0])
    assert counts(c) == [(2.0, 1.0), (6.0, 1.0), (10.0, 1.0)]
    ranges = {round(r, 9): m for r, m, _ in c}
    assert ranges[2.0] == pytest.approx(5.0)      # 4 and 6
    assert ranges[6.0] == pytest.approx(7.0)      # 4 and 10


def test_monotonic_history_is_a_single_half_cycle():
    """Nothing closes on a ramp: 0 to 3 is half of a cycle of range 3."""
    assert F.rainflow([0, 1, 2, 3]) == [(3.0, 1.5, 0.5)]


def test_astm_e1049_worked_example():
    """The history -2, 1, -3, 5, -1, 3, -4, 4, -2 of ASTM E1049-85.

    Worked by hand with the four-point rule: 1--(-1) closes first (range 2? no -- the inner pair that
    satisfies the rule is -1/3, range 4), then -3/5 and -4/4 remain open. The full hand count is
    ranges 3, 4, 4, 6, 8, 8, 9 with counts 0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, which is the published
    result for this series.
    """
    c = F.rainflow([-2, 1, -3, 5, -1, 3, -4, 4, -2])
    assert counts(c) == [(3.0, 0.5), (4.0, 1.5), (6.0, 0.5), (8.0, 1.0), (9.0, 0.5)]
    assert sum(x[2] for x in c) == pytest.approx(4.0)


@pytest.mark.parametrize("series", [
    [0, 5, 0, 5, 0, 5, 0],
    [0, 10, 4, 6, 4, 10, 0],
    [-2, 1, -3, 5, -1, 3, -4, 4, -2],
    [0, 1, 2, 3],
    list(np.sin(np.linspace(0, 40, 500)) * 7 + np.linspace(0, 3, 500)),
])
def test_every_reversal_but_the_first_belongs_to_exactly_one_half_cycle(series):
    """The counting invariant: total count = (number of reversals - 1) / 2, always."""
    n = len(F.reversals(series))
    assert sum(c[2] for c in F.rainflow(series)) == pytest.approx((n - 1) / 2)


def test_counting_is_invariant_to_offset_and_scales_with_amplitude():
    s = [0, 10, 4, 6, 4, 10, 0]
    base = counts(F.rainflow(s))
    assert counts(F.rainflow([x + 1000.0 for x in s])) == base
    doubled = counts(F.rainflow([2.0 * x for x in s]))
    assert [(r / 2, n) for r, n in doubled] == base


def test_a_cycle_is_never_wider_than_the_history():
    rng = np.random.default_rng(0)
    s = rng.normal(0, 1, 3000).cumsum()
    span = float(s.max() - s.min())
    for r, _m, _c in F.rainflow(s):
        assert r <= span + 1e-9


# ---------------------------------------------------------------------------
# Binning and the Basquin measure
# ---------------------------------------------------------------------------
def test_bins_conserve_the_total_count():
    c = F.rainflow([0, 10, 4, 6, 4, 10, 0])
    t = F.cycle_table(c)
    assert float(t.cycles.sum()) == pytest.approx(sum(x[2] for x in c))
    assert list(t.bin_MPa)[-1].startswith(">")


def test_basquin_measure_is_the_weighted_sum_of_powers():
    c = [(10.0, 0.0, 1.0), (5.0, 0.0, 0.5)]
    assert F.basquin_measure(c, 3) == pytest.approx(10.0 ** 3 + 0.5 * 5.0 ** 3)
    assert F.basquin_measure(c, 1) == pytest.approx(12.5)
    assert F.basquin_measure([], 5) == 0.0


def test_a_larger_exponent_favours_the_large_cycle():
    """The whole point of the exponent sweep: which side wins depends on m, so it must be swept.

    365 cycles of range 10 against one of range 40. They are equal when 365 = 4**m, i.e. m = 4.24,
    so the many small cycles win below that exponent and the single large one above it.
    """
    many_small = [(10.0, 0.0, 365.0)]
    one_big = [(40.0, 0.0, 1.0)]
    assert F.basquin_measure(many_small, 3) > F.basquin_measure(one_big, 3)
    assert F.basquin_measure(many_small, 8) < F.basquin_measure(one_big, 8)
    cross = np.log(365.0) / np.log(4.0)
    assert F.basquin_measure(many_small, cross) == pytest.approx(F.basquin_measure(one_big, cross))


# ---------------------------------------------------------------------------
# The stress relation
# ---------------------------------------------------------------------------
def test_thermal_stress_matches_the_relation_printed_in_section_s1():
    """E alpha dT / [2(1-nu)] with E = 150 GPa, alpha = 18e-6/K, nu = 0.3."""
    assert F.thermal_stress_MPa(25.0) == pytest.approx(150e9 * 18e-6 * 25.0 / 1.4 / 1e6)
    assert float(F.thermal_stress_MPa(0.0)) == 0.0
    assert F.thermal_stress_MPa(50.0) == pytest.approx(2 * F.thermal_stress_MPa(25.0))
    # the 20-30 K drop quoted in Section S1 gives the 45 MPa quoted there
    assert 38.0 < float(F.thermal_stress_MPa(20.0)) < 48.0
    assert 55.0 < float(F.thermal_stress_MPa(30.0)) < 65.0


def test_the_assumptions_are_declared_as_module_constants():
    """These are assumptions the report has to be able to state; they must not be buried."""
    assert F.STARTUPS_PER_YEAR == (2.0, 6.0)
    assert F.BASQUIN_M == (3.0, 5.0, 8.0)
    assert (F.E_PA, F.ALPHA_PER_K, F.NU) == (150.0e9, 18.0e-6, 0.3)


# ---------------------------------------------------------------------------
# Wiring to the twin
# ---------------------------------------------------------------------------
needs_v3 = pytest.mark.skipif(not (C.DATA / "scenarios" / "hourly_v3_S1_steady.csv.gz").exists(),
                              reason="v3 scenario results not present")


@needs_v3
def test_the_hot_spot_wall_drop_obeys_the_cylindrical_conduction_relation():
    """dT_wall must equal q_o d_o/(2 lambda) ln(d_o/d_i) -- the relation of rdt.reactor1d."""
    import math
    from rdt import latham_cases as lc
    from rdt import operate as op
    from rdt import pipeline as pl

    pl.apply_config(C.V3)
    base = op.base_case()
    got = F.hotspot_wall_drop(base.inputs(), base)
    tube = lc.build_case(base.row, op.calibrated_params())[0]
    expect = got["q_o_hot_W_m2"] * tube.d_o / (2 * tube.lambda_tube) * math.log(tube.d_o / tube.d_i)
    assert got["dT_wall_hot_K"] == pytest.approx(expect, rel=1e-9)
    assert got["dT_wall_hot_K"] == pytest.approx(got["T_wo_max_K"] - got["T_wi_at_hot_K"], rel=1e-12)
    assert 0.0 < got["dT_wall_hot_K"] < got["dT_wall_max_K"] + 1e-9


@needs_v3
def test_a_steady_scenario_has_a_flat_history_and_therefore_no_cycles():
    """S1 holds one operating point all year. If the pipeline invents cycles there, it is noise."""
    cached = F.OUT_DIR / "wall_drop" / f"wall_drop_v3_S1_steady.csv.gz"
    if not cached.exists():
        pytest.skip("run `python -m rdt.fatigue` first")
    h = pd.read_csv(cached)
    assert h.dT_wall_hot_K.nunique() == 1, "a steady scenario must map to a single distinct solve"
    assert F.rainflow(F.thermal_stress_MPa(h.dT_wall_hot_K.to_numpy(float))) == []
