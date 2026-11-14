import numpy as np
import pandas as pd
import pytest

from rdt import creep
from rdt import reactor1d as r1
from rdt import revision_r1b as rb

CHAIN_CSV = rb.OUT / "chain_terms_midheight.csv"


@pytest.mark.parametrize("case, alpha_i", [("Xu-Froment 1989", 3358.1), ("Plant A", 1562.8)])
def test_printed_chain_pins_table_S17(case, alpha_i):
    """Xu & Froment (1989) Part II, p. 100: Eq. (12) with alpha_w = alpha_w0 + 0.444 Re Pr lam_g/d_p and
    lam_er = lam_er0 + 0.14 lam_g Re Pr, evaluated by hand from the stored mid-height inputs."""
    t = pd.read_csv(CHAIN_CSV).query("case == @case").set_index("term").value
    lam, mu, cp, G = t["lambda_g"], t["mu_g"], t["cp_g"], t["G (superficial mass velocity)"]
    d_p, d_ti, lam0 = t["d_p equivalent"], t["d_t inner"], t["lambda_er0 (static, Kunii-Smith)"]
    Re, Pr = d_p * G / mu, cp * mu / lam
    a_w = 8.694 * lam0 / d_ti ** (4 / 3) + 0.444 * Re * Pr * lam / d_p
    l_er = lam0 + 0.14 * lam * Re * Pr
    by_hand = 1.0 / (1.0 / a_w + d_ti / (8.0 * l_er))
    out = r1.xu_froment_alpha_i(lam, mu, cp, G, d_p, d_ti, lam0)["alpha_i"]
    assert out == pytest.approx(by_hand, rel=1e-12)
    assert out == pytest.approx(alpha_i, abs=0.1)
    # the form quoted from papers citing De Wasch & Froment is a different correlation
    assert rb.quoted_alternative_form(lam, Re, Pr, d_p, d_ti, lam0) < 0.8 * out


def test_rescore_on_base_curve_is_identity():
    key = "centralloy_g_4852"; cur = creep.load_alloy(key).minimum
    df = pd.DataFrame({"T_wo_max_K": [1120.0, 1150.0], "sigma_hot_MPa": [12.9, 12.9], "C_LM": [cur.C, cur.C],
                       "log10_tr_scatter": [0.1, -0.2], "H2_net_kmol_h": [14.0, 15.0]})
    out = rb.rescore_mc(df, key, key)
    ref = np.log10(cur.time_to_rupture(df.T_wo_max_K.to_numpy(), df.sigma_hot_MPa.to_numpy())) + df.log10_tr_scatter.to_numpy()
    assert np.allclose(out.log10_t_r, ref, atol=1e-12)
    assert np.allclose(out.life_rate_per_kmol_H2, 10.0 ** -ref / df.H2_net_kmol_h)
