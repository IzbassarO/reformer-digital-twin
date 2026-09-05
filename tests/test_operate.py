"""Operating-space dataset and closed-loop helper (rdt.operate)."""

import json

import numpy as np
import pandas as pd
import pytest

from rdt import operate as op


@pytest.fixture(scope="module")
def base():
    return op.base_case()


@pytest.fixture(scope="module")
def base_out(base):
    return op.base_run(base)


def test_lhs_samples_inside_bounds():
    lo, up = op.bounds_arrays()
    X = op.lhs_samples(200, seed=0)
    assert np.all(X.to_numpy() >= lo) and np.all(X.to_numpy() <= up)
    assert list(X.columns) == list(op.INPUT_NAMES)


def test_campaign_convergence_rate():
    meta_path = op.LHS_DIR / "lhs_plantA_v1_meta.json"
    if not meta_path.exists():
        pytest.skip("LHS campaign not run yet")
    meta = json.loads(meta_path.read_text())
    assert meta["convergence_rate"] >= 0.98, meta["convergence_rate"]
    df = pd.read_csv(op.LHS_DIR / "lhs_plantA_v1.csv.gz")
    lo, up = op.bounds_arrays()
    assert np.all(df[list(op.INPUT_NAMES)].to_numpy() >= lo - 1e-9) and np.all(df[list(op.INPUT_NAMES)].to_numpy() <= up + 1e-9)
    assert (df.T_wo_max_K > df.T_out_K).all()


def test_more_firing_raises_wall_temperature(base, base_out):
    x = base.inputs(); x["specific_firing_factor"] = 1.10
    r = op.run_case(x, base, base_rate=base_out["life_consumption_rate_per_h"])
    assert r["converged"] and r["T_wo_max_K"] > base_out["T_wo_max_K"] + 5.0
    assert r["life_consumption_ratio_to_base"] > 1.5


def test_closed_loop_reproduces_base_within_1K(base, base_out):
    r = op.solve_firing_for_outlet_T(base.inputs(), op.TARGET_T_OUT_K, base, base_rate=base_out["life_consumption_rate_per_h"])
    assert r["converged"] and abs(r["T_out_K"] - op.TARGET_T_OUT_K) < 1.0
    assert 0.9 < r["specific_firing_factor"] < 1.1   # measured Plant A firing reproduced within 10 %


def test_closed_loop_higher_sc_lowers_ch4_slip(base, base_out):
    br = base_out["life_consumption_rate_per_h"]
    lo = op.solve_firing_for_outlet_T({**base.inputs(), "steam_to_carbon": 2.6}, base=base, base_rate=br)
    hi = op.solve_firing_for_outlet_T({**base.inputs(), "steam_to_carbon": 3.8}, base=base, base_rate=br)
    assert lo["converged"] and hi["converged"]
    assert hi["CH4_slip_dry_pct"] < lo["CH4_slip_dry_pct"]
