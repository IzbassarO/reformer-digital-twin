"""RQ2 scenario machinery (rdt.scenarios)."""

import json

import numpy as np
import pandas as pd
import pytest

from rdt import scenarios as sc

SUMMARY = sc.OUT_DIR / "summary_v1.csv"


@pytest.fixture(scope="module")
def tw():
    if not sc.sg.METRICS_JSON.exists():
        pytest.skip("surrogates not trained")
    return sc.TwinSurrogate()


@pytest.fixture(scope="module")
def summary():
    if not SUMMARY.exists():
        pytest.skip("scenario campaign not run yet")
    return pd.read_csv(SUMMARY).set_index("scenario")


def test_s1_hourly_equals_closed_form(tw):
    ss = sc.steady_state_closed_form(tw)
    h = sc.make_scenarios()["S1_steady"]
    r = sc.evaluate(h, tw, verify_every=8760)
    assert abs(r["annual_damage"] / ss["D_year"] - 1.0) < 1e-6
    assert r["avg_condition_error"] == pytest.approx(1.0, abs=1e-6)


def test_ramp_limits_respected():
    S = sc.make_scenarios()
    for name in ("S2_daily", "S3_renewable"):
        d = np.abs(np.diff(S[name].load))
        assert d.max() <= 0.10 + 1e-9, (name, d.max())
        assert S[name].load.min() >= 0.60 - 1e-9 and S[name].load.max() <= 1.10 + 1e-9


def test_surrogate_matches_physics_on_sampled_hours(summary):
    for name, row in summary.iterrows():
        assert row.verify_max_abs_dT_wo_K < 1.0, (name, row.verify_max_abs_dT_wo_K)


def test_severe_overfiring_adds_damage(summary):
    assert summary.loc["S4b_severe_overfire", "annual_damage_D"] > summary.loc["S1_steady", "annual_damage_D"]
    assert summary.loc["S4a_mild_overfire", "annual_damage_D"] > summary.loc["S1_steady", "annual_damage_D"]
    assert summary.loc["S4c_hot_band", "annual_damage_D"] > summary.loc["S1_steady", "annual_damage_D"]


def test_hold_slip_raises_wall_temperature_over_campaign():
    meta_path = sc.OUT_DIR / "summary_v1_meta.json"
    if not meta_path.exists():
        pytest.skip("scenario campaign not run yet")
    camp = json.loads(meta_path.read_text())["S5_campaigns"]["hold_CH4_slip_ageing"]
    assert camp[-1]["T_wo_end"] > camp[0]["T_wo_start"] + 5.0
    assert camp[-1]["firing_end"] > camp[0]["firing_start"]
