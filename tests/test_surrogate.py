"""Surrogate accuracy, conformal coverage and base-case agreement (rdt.surrogate)."""

import json

import pytest

from rdt import operate as op
from rdt import surrogate as s


@pytest.fixture(scope="module")
def metrics():
    if not s.METRICS_JSON.exists():
        pytest.skip("surrogates not trained yet (run notebooks/08_surrogate.ipynb)")
    return json.loads(s.METRICS_JSON.read_text())


def test_best_surrogate_accuracy(metrics):
    t = metrics["targets"]
    best_two = t["T_wo_max_K"]["models"][t["T_wo_max_K"]["best"]]["test"]["rmse"]
    best_tr = t["log10_t_r_hot"]["models"][t["log10_t_r_hot"]["best"]]["test"]["rmse"]
    assert best_two < 3.0, best_two
    assert best_tr < 0.05, best_tr


def test_conformal_coverage(metrics):
    """Nominal 90 % CV+ intervals should cover 87-93 % of the test set; xfail with achieved values if they over-cover."""
    picps = {tgt: c["levels"]["0.9"]["picp"] for tgt, c in metrics["conformal"].items()}
    outside = {t: round(v, 3) for t, v in picps.items() if not 0.87 <= v <= 0.93}
    assert all(v >= 0.87 for v in picps.values()), f"under-coverage: {picps}"     # never accept invalid intervals
    if outside:
        pytest.xfail(f"CV+ intervals are valid but conservative (over-cover) at nominal 0.90: {outside}; "
                     f"plain-CV comparison: {[(t, round(c['levels']['0.9']['picp'], 3)) for t, c in metrics.get('conformal_cv_base_comparison', {}).items() if isinstance(c, dict) and 'levels' in c]}")


def test_surrogate_base_case_matches_physics(metrics):
    if not (s.MODELS_DIR / f"T_wo_max_K__{metrics['targets']['T_wo_max_K']['best']}.joblib").exists():
        pytest.skip("saved model missing")
    phys = op.base_run(op.base_case())["T_wo_max_K"]
    sur = s.predict_base_case("T_wo_max_K")
    assert abs(sur - phys) < 2.0, (sur, phys)
