"""Run configurations (v1 = original campaigns, v2 = consolidated release with the corrected operating space).

Every result-producing stage takes a :class:`RunConfig` and writes files tagged with ``cfg.tag`` next to the
v1 files; v1 files are never overwritten.
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
    def runtimes_json(self) -> Path:
        return DATA / f"pipeline_{self.tag}_runtimes.json"


V1 = RunConfig()
V2 = RunConfig(tag="v2", ranges_csv=DATA / "design" / "parameter_ranges_v2.csv", excess_air_base=21.64128020288727, n_mc=2048)
# v2 notes: excess-air range 5-25 % so that the calibrated Plant A base (21.6 %) is inside; base case = calibrated state;
# Monte Carlo N = 2048 (half of v1) to keep the UQ stage near 20 minutes on 10 cores (the v1 4096-sample run took 34 min).
