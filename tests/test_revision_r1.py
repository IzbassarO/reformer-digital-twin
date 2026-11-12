import numpy as np
import pandas as pd

from rdt import revision_r1 as r1


def test_nondominated():
    F = np.array([[0, 0, 0], [1, 1, 1], [0, 1, -1], [0, 0, 0]])
    assert r1.nondominated(F).tolist() == [True, False, True, True]


def test_flex_ratios_identity():
    d = pd.DataFrame({"converged": [True, True], "life_rate_per_h": [2.0, 3.0], "H2_net_kmol_h": [4.0, 5.0],
                      "life_rate_per_kmol_H2": [0.5, 0.6]})
    s2, s3, ok = r1.flex_ratios(d, d, d)
    assert ok.all() and np.allclose(s2, 1.0) and np.allclose(s3, 1.0)


def test_make_bins_conserves_hours_and_separates_events():
    from rdt import operate as op
    n = 100
    h = pd.DataFrame({k: np.ones(n) for k in op.INPUT_NAMES})
    h["dD"] = 1e-7; h.loc[:4, "dD"] = 1e-6; h["T_wo_add_K"] = 0.0; h.loc[:4, "T_wo_add_K"] = 40.0
    b = r1.make_bins(h, 6)
    assert len(b) == 2 and b.hours.sum() == n and sorted(b.T_wo_add_K) == [0.0, 40.0]
