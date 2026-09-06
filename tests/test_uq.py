"""RQ4 uncertainty propagation (rdt.uq)."""

import json

import numpy as np
import pandas as pd
import pytest

from rdt import uq

SUMMARY = uq.OUT_DIR / "uq_summary_v1.json"


@pytest.fixture(scope="module")
def summary():
    if not SUMMARY.exists():
        pytest.skip("Monte Carlo not run yet")
    return json.loads(SUMMARY.read_text())


def test_spec_and_sampling():
    spec = uq.build_spec(save=False)
    X = uq.sample(spec, 64)
    assert list(X.columns) == uq.PARAMS
    cs = spec["groups"]["creep"]
    assert X.wall_thickness_m.between(cs["wall_thickness_m"]["low"], cs["wall_thickness_m"]["high"]).all()
    assert abs(X.E1.mean() - spec["groups"]["kinetics"]["E1"]["mean"]) < 2.0
    Xk = uq.sample(spec, 16, active_groups=["creep"])
    nom = uq.nominal(spec)
    assert (Xk.E1 == nom["E1"]).all() and (Xk.F_gt == nom["F_gt"]).all() and Xk.wall_thickness_m.std() > 0


def test_paired_ratio_of_base_is_one():
    df = pd.read_csv(uq.OUT_DIR / "mc_v1_S1_base.csv.gz") if (uq.OUT_DIR / "mc_v1_S1_base.csv.gz").exists() else None
    if df is None:
        pytest.skip("Monte Carlo not run yet")
    iv = uq.intervals(df, df)
    assert iv["life_rate_rel_base"]["median"] == 1.0 and iv["life_rate_rel_base"]["p05"] == 1.0 and iv["life_rate_rel_base"]["p95"] == 1.0


def test_intervals_contain_nominal(summary):
    for n, d in summary["nominal_inside_90pct_interval"].items():
        for c, ok in d.items():
            assert ok, (n, c)


def test_variance_shares_close_to_joint(summary):
    vs = summary["variance_shares_S1_base"]
    for c, v in vs.items():
        assert abs(v["sum_of_shares"] - 1.0) <= 0.2, (c, v["sum_of_shares"], "interaction residual", v["interaction_residual"])


def test_robustness_fractions_computed(summary):
    rf = summary["robustness_fractions"]
    for tag in ("nominal_L_q", "load_dependent_L_q"):
        for k in ("S5y4_life_gt_2x_S1", "knee_life_lt_0.5_base", "SC3.5_less_life_per_kmol_than_SC2.5", "S2_life_per_kmol_le_S1", "S3_life_per_kmol_le_S1"):
            assert k in rf[tag] and 0.0 <= rf[tag][k] <= 1.0, (tag, k)
