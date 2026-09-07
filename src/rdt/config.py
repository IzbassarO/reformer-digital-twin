"""Run configurations (v1 = original campaigns, v2 = consolidated release, v3 = manufacturer creep curves).

Every result-producing stage takes a :class:`RunConfig` and writes files tagged with ``cfg.tag`` next to the
files of the earlier versions, which are never overwritten.

``creep_alloy`` selects the Larson-Miller master curve (see :func:`rdt.creep.available_alloys`). v1 and v2
pin the legacy Yeh Manaurite XM placeholder so that their published results stay reproducible; v3 uses the
Schmidt + Clemens Centralloy G 4852 (HP-Nb) data-sheet curves digitised by
:mod:`rdt.creep_ingest_datasheet`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"


@dataclass(frozen=True)
class RunConfig:
    tag: str = "v1"
    ranges_csv: Path = DATA / "design" / "parameter_ranges.csv"
    excess_air_base: float = 20.0        # base excess air used by scenarios / optimisation (v1: top of the 5-20 % range)
    n_lhs: int = 2000
    lhs_seed: int = 0
    sobol_N: int = 256                   # Saltelli base sample (4096 runs at 7 inputs)
    grid_n: int = 15
    n_mc: int = 4096                     # RQ4 Monte Carlo samples
    n_mc_group: int = 1024
    nsga_pop: int = 200
    nsga_gens: int = 300
    creep_alloy: str = "centralloy_g_4852"   # Larson-Miller master curve; see rdt.creep.available_alloys()

    # ---- file names -------------------------------------------------------
    @property
    def lhs_name(self) -> str:
        return f"lhs_plantA_{self.tag}"

    @property
    def grid_name(self) -> str:
        return f"grid_sc_load_{self.tag}"

    @property
    def sobol_json(self) -> Path:
        return DATA / "design" / f"sobol_{self.tag}.json"

    @property
    def saltelli_file(self) -> Path:
        return DATA / "lhs_runs" / f"sobol_saltelli_{self.tag}.csv.gz"

    @property
    def lhs_file(self) -> Path:
        return DATA / "lhs_runs" / f"{self.lhs_name}.csv.gz"

    @property
    def surrogate_metrics(self) -> Path:
        return DATA / "design" / f"surrogate_metrics_{self.tag}.json"

    @property
    def models_dir(self) -> Path:
        return ROOT / "models" if self.tag == "v1" else ROOT / "models" / self.tag

    @property
    def scenarios_summary(self) -> Path:
        return DATA / "scenarios" / f"summary_{self.tag}.csv"

    @property
    def scenarios_meta(self) -> Path:
        return DATA / "scenarios" / f"summary_{self.tag}_meta.json"

    def scenario_hourly(self, name: str) -> Path:
        return DATA / "scenarios" / (f"hourly_{name}.csv.gz" if self.tag == "v1" else f"hourly_{self.tag}_{name}.csv.gz")

    @property
    def pareto_csv(self) -> Path:
        return DATA / "optimization" / f"pareto_{self.tag}.csv"

    @property
    def regimes_csv(self) -> Path:
        return DATA / "optimization" / f"regimes_{self.tag}.csv"

    @property
    def pareto_meta(self) -> Path:
        return DATA / "optimization" / f"pareto_{self.tag}_meta.json"

    @property
    def pareto_alt_csv(self) -> Path:
        return DATA / "optimization" / f"pareto_{self.tag}_alt_totalheat.csv"

    @property
    def pareto_alt_meta(self) -> Path:
        return DATA / "optimization" / f"pareto_{self.tag}_alt_totalheat_meta.json"

    @property
    def steam_credit_json(self) -> Path:
        return DATA / "optimization" / f"steam_credit_{self.tag}.json"

    @property
    def steam_credit_fronts(self) -> Path:
        return DATA / "optimization" / f"steam_credit_fronts_{self.tag}.csv"

    @property
    def uq_dir(self) -> Path:
        return DATA / "uq"

    @property
    def mc_tag(self) -> str:
        return f"mc_{self.tag}"

    @property
    def uq_summary(self) -> Path:
        return DATA / "uq" / f"uq_summary_{self.tag}.json"

    @property
    def uq_spec(self) -> Path:
        # v1/v2 wrote a single untagged spec; keep that name for them and tag it from v3 on so that
        # a later run can never overwrite the spec a published version was produced with.
        return DATA / "uq" / ("uncertainty_spec.yaml" if self.tag in ("v1", "v2") else f"uncertainty_spec_{self.tag}.yaml")

    @property
    def figures_dir(self) -> Path:
        return ROOT / "paper" / "figures" / ("" if self.tag in ("v1", "v2") else self.tag)

    @property
    def runtimes_json(self) -> Path:
        return DATA / f"pipeline_{self.tag}_runtimes.json"


LEGACY_CREEP_ALLOY = "legacy_yeh_manaurite_xm"

V1 = RunConfig(creep_alloy=LEGACY_CREEP_ALLOY)
V2 = RunConfig(tag="v2", ranges_csv=DATA / "design" / "parameter_ranges_v2.csv", excess_air_base=21.64128020288727,
               n_mc=2048, creep_alloy=LEGACY_CREEP_ALLOY)
V3 = RunConfig(tag="v3", ranges_csv=DATA / "design" / "parameter_ranges_v2.csv", excess_air_base=21.64128020288727,
               n_mc=2048, creep_alloy="centralloy_g_4852")
# v2 notes: excess-air range 5-25 % so that the calibrated Plant A base (21.6 %) is inside; base case = calibrated state;
# Monte Carlo N = 2048 (half of v1) to keep the UQ stage near 20 minutes on 10 cores (the v1 4096-sample run took 34 min).
# v3 notes: identical operating space and sample sizes to v2; the only change is the creep master curve, which moves
# from the Yeh placeholder to the Schmidt + Clemens Centralloy G 4852 data sheet (C = 18.6, own scatter band).

BY_TAG = {c.tag: c for c in (V1, V2, V3)}
