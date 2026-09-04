"""Consistency checks on the transcribed Xu & Froment (1989) kinetic parameters."""

from pathlib import Path

import pytest
import yaml

DATA = Path(__file__).resolve().parents[1] / "data" / "kinetics" / "xu_froment_1989.yaml"


@pytest.fixture(scope="module")
def kin():
    with DATA.open() as f:
        return yaml.safe_load(f)


def test_activation_energies(kin):
    e = kin["table5_activation_energies_and_adsorption_enthalpies"]
    assert e["E1"]["value"] == pytest.approx(240.1)
    assert e["E2"]["value"] == pytest.approx(67.13)
    assert e["E3"]["value"] == pytest.approx(243.9)


def test_pre_exponential_K_H2O_positive(kin):
    assert kin["table6_pre_exponential_factors"]["A_K_H2O"]["value"] > 0


def test_adsorption_enthalpy_signs(kin):
    e = kin["table5_activation_energies_and_adsorption_enthalpies"]
    assert e["dH_H2O"]["value"] > 0
    for key in ("dH_CO", "dH_H2", "dH_CH4"):
        assert e[key]["value"] < 0, key


def test_confidence_bounds_bracket_estimates(kin):
    for block in ("table5_reference_temperature_values",
                  "table5_activation_energies_and_adsorption_enthalpies"):
        for name, p in kin[block].items():
            if not isinstance(p, dict) or "UL" not in p:
                continue
            assert p["LL"] <= p["value"] <= p["UL"], name
