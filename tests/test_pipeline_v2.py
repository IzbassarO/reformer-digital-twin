"""Consolidated v2 release: configuration, result files and base-case-inside-space checks."""

import json

import pandas as pd
import pytest

from rdt import config as C


def test_v2_ranges_include_calibrated_base():
    r = pd.read_csv(C.V2.ranges_csv).set_index("parameter")
    assert r.loc["excess_air", "max"] == 25.0 and r.loc["excess_air", "min"] == 5.0
    assert r.loc["excess_air", "min"] <= C.V2.excess_air_base <= r.loc["excess_air", "max"]
    v1 = pd.read_csv(C.V1.ranges_csv).set_index("parameter")
    for k in ("steam_to_carbon", "feed_per_tube_fraction", "specific_firing_factor", "inlet_T", "inlet_P", "catalyst_activity"):
        assert (r.loc[k, ["min", "max"]] == v1.loc[k, ["min", "max"]]).all(), k


def test_v2_results_exist_and_v1_untouched():
    for p in (C.V2.lhs_file, C.V2.saltelli_file, C.V2.sobol_json, C.V2.surrogate_metrics, C.V2.scenarios_summary, C.V2.pareto_csv, C.V2.regimes_csv,
              C.V2.steam_credit_json, C.V2.uq_summary, C.V2.runtimes_json):
        if not p.exists():
            pytest.skip(f"v2 pipeline not run ({p.name} missing)")
    for p in (C.V1.lhs_file, C.V1.sobol_json, C.V1.scenarios_summary, C.V1.pareto_csv, C.V1.uq_summary):
        assert p.exists(), p
    rt = json.loads(C.V2.runtimes_json.read_text())
    assert set(rt) >= {"design", "surrogate", "scenarios", "pareto", "uq"}
    assert all(v["wall_time_s"] > 0 for v in rt.values())


def test_v2_uq_sample_size_documented():
    if not C.V2.uq_summary.exists():
        pytest.skip("v2 UQ not run")
    S = json.loads(C.V2.uq_summary.read_text())
    assert S["N"] == C.V2.n_mc == 2048


def test_diff_report_exists():
    p = C.DATA / "design" / "v1_vs_v2_diff.md"
    if not p.exists():
        pytest.skip("diff report not generated")
    txt = p.read_text()
    assert "| family | item | v1 | v2 |" in txt and "changed" in txt
