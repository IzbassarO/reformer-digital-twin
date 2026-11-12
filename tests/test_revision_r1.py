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


def test_activity_profiles():
    a = r1.activity_profile("E0"); assert len(a) == 4 * 8760 and np.isclose(a[0], 0.30) and abs(a[-1] - 0.10) < 1e-4
    e3 = r1.activity_profile("E3"); i6 = int(0.5 * 8760)
    assert np.isclose(e3[i6 - 1], 0.18, atol=1e-3) and np.isclose(e3[i6], 0.18, atol=1e-3) and abs(e3[-1] - 0.10) < 1e-4
    assert len(r1.activity_profile("E4")) == 3 * 8760 and len(r1.activity_profile("E5")) == 5 * 8760


def test_c1_restored_feed_closes_balances():
    from rdt import latham_cases as lc
    c = lc.load_cases().set_index("case")
    b = r1.atom_balance(c.loc["Plant_C1"])
    assert abs(b["feed_from_C_balance_kmol_h"] - 7190) < 10
    b2 = r1.atom_balance(c.loc["Plant_C1"], b["feed_from_C_balance_kmol_h"])
    assert abs(b2["C_closure_out_over_in"] - 1) < 1e-12 and abs(b2["H_closure_out_over_in"] - 1) < 1e-3


def test_e_upper_line_and_weight():
    assert np.isclose(r1.e_upper(0.99), -12.8) and np.isclose(r1.e_upper(0.71), 18.2)
    assert np.isclose(r1.W_HOT, (0.682 - 0.43) / (0.682 - 0.293))


def test_student_t_scatter_matches_lower_band():
    from scipy import stats
    s = 0.181
    x = np.array([stats.norm.ppf(0.05, 0, s)])
    assert np.isclose(r1._scatter(x, s, "t")[0], -1.645 * s, rtol=1e-3)
    assert np.allclose(r1._scatter(x, s, "normal"), x)
