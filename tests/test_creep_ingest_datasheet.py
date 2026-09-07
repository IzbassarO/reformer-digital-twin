"""Digitisation of the manufacturer stress-rupture charts (rdt.creep_ingest_datasheet).

The Schmidt + Clemens PDFs are copyright and are not redistributed, so the tests that need them are
skipped when ``data/creep_derived/datasheets/`` is empty. The unit tests of the geometry and the
schema tests of the committed CSVs always run.
"""

import numpy as np
import pandas as pd
import pytest

from rdt import creep_ingest_datasheet as ing

SPECS = ing.load_sources()
HAVE_PDFS = all(s.pdf_path.exists() for s in SPECS.values())
needs_pdfs = pytest.mark.skipif(not HAVE_PDFS, reason="Schmidt+Clemens data sheets not present (not redistributed)")


# --- geometry ---------------------------------------------------------------
def test_flatten_path_matches_analytic_bezier():
    p0, p1, p2, p3 = (0.0, 0.0), (1.0, 2.0), (3.0, 2.0), (4.0, 0.0)
    got = ing.flatten_path([("m", p0), ("c", p1, p2, p3)], n_per_segment=101)
    t = np.linspace(0, 1, 101)[:, None]
    P = np.array([p0, p1, p2, p3], float)
    want = ((1 - t) ** 3 * P[0] + 3 * (1 - t) ** 2 * t * P[1] + 3 * (1 - t) * t**2 * P[2] + t**3 * P[3])
    assert np.allclose(got, want, atol=1e-12)


def test_flatten_path_v_and_y_operators_imply_control_points():
    """``v`` reuses the current point as control 1; ``y`` reuses the endpoint as control 2."""
    v = ing.flatten_path([("m", (0.0, 0.0)), ("v", (1.0, 1.0), (2.0, 0.0))], n_per_segment=51)
    y = ing.flatten_path([("m", (0.0, 0.0)), ("y", (1.0, 1.0), (2.0, 0.0))], n_per_segment=51)
    ref_v = ing.flatten_path([("m", (0.0, 0.0)), ("c", (0.0, 0.0), (1.0, 1.0), (2.0, 0.0))], n_per_segment=51)
    ref_y = ing.flatten_path([("m", (0.0, 0.0)), ("c", (1.0, 1.0), (2.0, 0.0), (2.0, 0.0))], n_per_segment=51)
    assert np.allclose(v, ref_v) and np.allclose(y, ref_y)
    assert not np.allclose(v, y)


def test_axis_round_trip_linear_and_log():
    lin = ing.Axis("lmp", "top", "linear", ((25.0, 762.4), (32.0, 244.1)))
    log = ing.Axis("stress", "x", "log10", ((100.0, 181.7), (1.0, 555.6)))
    for ax, vals in ((lin, [25.0, 28.5, 32.0]), (log, [1.0, 10.0, 100.0])):
        assert np.allclose(ax.to_data(ax.to_page(vals)), vals)
    # the calibration ticks themselves must land exactly
    assert log.to_data(181.7) == pytest.approx(100.0)
    assert log.to_data(555.6) == pytest.approx(1.0)
    assert log.to_data(0.5 * (181.7 + 555.6)) == pytest.approx(10.0, rel=1e-9)


def test_axis_rejects_degenerate_calibration():
    with pytest.raises(ValueError):
        ing.Axis("lmp", "top", "linear", ((25.0, 300.0), (32.0, 300.0))).to_data(300.0)


# --- sources.yaml -----------------------------------------------------------
def test_sources_declare_the_three_alloys_with_their_own_constants():
    assert set(SPECS) == {"G4852", "G4852Micro", "ET45Micro"}
    # read off each sheet individually; they are NOT all 22.9
    assert SPECS["G4852"].C == 18.6
    assert SPECS["G4852Micro"].C == 22.9
    assert SPECS["ET45Micro"].C == 19.3
    for s in SPECS.values():
        assert s.chart_title == "Parametric stress rupture strength" and s.page == 6
        assert s.lmp.transform == "linear" and s.stress.transform == "log10"


# --- committed CSVs ---------------------------------------------------------
@pytest.mark.parametrize("key", sorted(SPECS))
def test_committed_curve_csv_schema_and_monotonicity(key):
    spec = SPECS[key]
    df = pd.read_csv(spec.out_csv)
    assert list(df.columns) == ["alloy", "curve", "lmp", "stress_mpa", "source_file", "page", "method", "extracted_on"]
    assert set(df.curve) == {"average", "minimum"}
    assert (df.source_file == spec.datasheet_file).all() and (df.page == spec.page).all()
    assert set(df.method) <= {"vector", "raster"}
    # the curves may run past the last labelled gridline, but never outside the drawn plot frame
    edges = spec.lmp.to_data(np.array([spec.frame["top"], spec.frame["bottom"]]))
    lo, hi = float(min(edges)), float(max(edges))
    for kind, g in df.groupby("curve"):
        assert len(g) >= 100
        assert g.lmp.is_monotonic_increasing
        # rupture stress must fall as the parameter rises, and stay inside the printed axis ranges
        assert g.stress_mpa.iloc[-1] < g.stress_mpa.iloc[0]
        assert lo - 0.05 <= g.lmp.min() and g.lmp.max() <= hi + 0.05
        assert 1.0 <= g.stress_mpa.min() and g.stress_mpa.max() <= 100.0


@pytest.mark.parametrize("key", sorted(SPECS))
def test_average_curve_lies_above_the_lower_scatter_band(key):
    df = pd.read_csv(SPECS[key].out_csv)
    a = df[df.curve == "average"]; m = df[df.curve == "minimum"]
    grid = np.linspace(max(a.lmp.min(), m.lmp.min()), min(a.lmp.max(), m.lmp.max()), 200)
    gap = np.interp(grid, a.lmp, np.log10(a.stress_mpa)) - np.interp(grid, m.lmp, np.log10(m.stress_mpa))
    assert gap.min() > 0.0, f"{key}: curves cross"


# --- extraction against the real PDFs ---------------------------------------
@needs_pdfs
@pytest.mark.parametrize("key", sorted(SPECS))
def test_vector_and_raster_extraction_agree(key):
    """The two independent extraction paths must land on the same curve.

    Vector reads the Bezier control points; raster traces a 600 dpi rendering pixel by pixel. They share
    no code beyond the axis calibration, so agreement validates the calibration, the CropBox offset of
    the rendering and both tracers at once.
    """
    spec = SPECS[key]
    v, r = ing.extract_vector(spec), ing.extract_raster(spec, dpi=600)
    for kind in ing.CURVE_KINDS:
        a, b = v[kind], r[kind]
        grid = np.linspace(max(a[:, 0].min(), b[:, 0].min()), min(a[:, 0].max(), b[:, 0].max()), 300)
        rel = 10 ** np.abs(np.interp(grid, a[:, 0], np.log10(a[:, 1])) - np.interp(grid, b[:, 0], np.log10(b[:, 1]))) - 1
        assert np.median(rel) < 0.005, f"{key}/{kind}: median disagreement {np.median(rel):.4f}"
        assert np.percentile(rel, 99) < 0.02, f"{key}/{kind}: p99 disagreement {np.percentile(rel, 99):.4f}"


@needs_pdfs
def test_digitise_reproduces_the_committed_csv(tmp_path, monkeypatch):
    monkeypatch.setattr(ing, "DATA_DIR", tmp_path)
    spec = SPECS["G4852"]
    fresh = ing.digitise(spec, method="vector", verbose=False)
    old = pd.read_csv(spec.out_csv)
    for kind in ing.CURVE_KINDS:
        f = fresh[fresh.curve == kind]; o = old[old.curve == kind]
        assert len(f) == len(o)
        assert np.allclose(f.lmp.to_numpy(), o.lmp.to_numpy())
        assert np.allclose(f.stress_mpa.to_numpy(), o.stress_mpa.to_numpy())


@needs_pdfs
def test_missing_datasheet_raises_a_pointed_error(tmp_path):
    import dataclasses

    spec = dataclasses.replace(SPECS["G4852"], datasheet_dir=tmp_path)
    with pytest.raises(FileNotFoundError, match="not redistributed"):
        ing.digitise(spec, verbose=False)
