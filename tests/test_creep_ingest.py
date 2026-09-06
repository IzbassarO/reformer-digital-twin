"""NIMS-ready creep ingestion: synthetic recovery of C and heat scatter, config switch."""

import numpy as np
import pytest

from rdt import creep


def test_template_columns_exist():
    import pandas as pd
    t = pd.read_csv(creep.DATA_DIR / "nims_transcription_template.csv")
    assert list(t.columns) == creep.TEMPLATE_COLUMNS and len(t) == 0


def test_from_config_selects_active_curve():
    c = creep.LarsonMillerCurve.from_config()
    assert c.C == pytest.approx(22.96)


def test_ingest_recovers_C_and_scatter(tmp_path):
    df = creep.synthetic_nims_dataset(n_heats=100, n_per_heat=8, sigma_heat=0.3, noise=0.05, seed=0)
    csv = tmp_path / "synthetic.csv"; df.to_csv(csv, index=False)
    out = creep.ingest_nims(csv, out_dir=tmp_path)
    d = out["SYNTH_XM"]["derived"]
    assert abs(d["C"] - 22.96) / 22.96 < 0.10, d["C"]
    assert abs(d["sigma_log10_heat"] - 0.3) / 0.3 < 0.10, d["sigma_log10_heat"]
    assert out["SYNTH_XM"]["n_heats"] == 100 and out["SYNTH_XM"]["n_points"] == 800
    # the derived file loads as a curve and reproduces the generating curve within the data range
    c_new = creep.LarsonMillerCurve.from_yaml(tmp_path / "SYNTH_XM_derived.yaml")
    c_ref = creep.LarsonMillerCurve.from_yaml()
    for T, s in ((1150.0, 12.0), (1200.0, 8.0)):
        t_new = float(c_new.time_to_rupture(T, s)); t_ref = float(c_ref.time_to_rupture(T, s))
        assert abs(np.log10(t_new) - np.log10(t_ref)) < 0.15, (T, s, t_new, t_ref)
    assert not (tmp_path / "SYNTH_XM_derived.yaml").read_text().count("heat_id")   # no raw data in the derived file
