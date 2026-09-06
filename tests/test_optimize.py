"""RQ3 Pareto optimisation results (rdt.optimize)."""

import json

import pandas as pd
import pytest

from rdt import optimize as ox

META = ox.OUT_DIR / "pareto_v1_meta.json"


@pytest.fixture(scope="module")
def meta():
    if not META.exists():
        pytest.skip("optimisation not run yet")
    return json.loads(META.read_text())


@pytest.fixture(scope="module")
def regimes():
    return pd.read_csv(ox.OUT_DIR / "regimes_v1.csv", index_col=0)


@pytest.fixture(scope="module")
def pareto():
    return pd.read_csv(ox.OUT_DIR / "pareto_v1.csv")


def test_base_domination_is_documented(meta):
    d = meta["base_domination"]
    assert {"n_dominating", "n_pareto", "base_is_dominated"} <= set(d)
    if d["base_is_dominated"]:
        assert all(k in d for k in ("max_H2_gain_pct", "max_fuel_saving_pct", "max_life_saving_pct"))


def test_regimes_feasible_under_physics(regimes, meta):
    s = meta["settings"]
    for name, r in regimes.iterrows():
        assert r.CH4_slip_dry_pct_phys <= s["slip_limit_pct"] + 1e-3, (name, r.CH4_slip_dry_pct_phys)
        assert r.T_wo_max_K_phys <= s["T_wo_limit_K"] + 1e-3, (name, r.T_wo_max_K_phys)
        assert r.T_out_K_phys <= s["T_out_limit_K"] + 1e-3, (name, r.T_out_K_phys)


def test_min_life_iso_production_beats_base(regimes):
    if "min_life_iso_H2" not in regimes.index:
        pytest.skip("no Pareto member within the +/-1 % production band")
    assert regimes.loc["min_life_iso_H2", "life_rate_per_kmol_H2_rel_base_phys"] < 1.0
    assert abs(regimes.loc["min_life_iso_H2", "H2_rel_base_pct"]) <= 1.0


def test_surrogate_physics_discrepancy_small(meta, pareto):
    assert meta["discrepancy_max"]["H2_net_kmol_h"] < 0.05
    assert meta["discrepancy_max"]["life_rate_per_kmol_H2_rel_base"] < 0.1
    assert pareto.feasible_phys.mean() > 0.95
