"""Smoke test: manuscript tables and number macros are generated from the v2 result files."""

import re

import pytest

from rdt import config as C

pytest.importorskip("yaml")


def test_tables_and_numbers_generate(tmp_path, monkeypatch):
    from rdt import paper_tables as T

    monkeypatch.setattr(T, "TAB", tmp_path / "tables")
    monkeypatch.setattr(T, "NUMBERS", tmp_path / "numbers.tex")
    if not C.V2.uq_summary.exists():
        pytest.skip("v2 results not present")
    out = T.make_all()
    assert len(out) == 10
    for name, p in out.items():
        txt = p.read_text()
        assert p.exists() and txt.strip()
        if name != "numbers":
            assert "\\begin{table}" in txt and "\\bottomrule" in txt and "Source:" in txt
    macros = re.findall(r"\\newcommand\{\\(\w+)\}\{([^}]*)\}", out["numbers"].read_text())
    assert len(macros) > 60
    assert all(v and v != "nan" for _, v in macros)
