# Reformer Digital Twin

[![DOI](https://zenodo.org/badge/1357536128.svg)](https://doi.org/10.5281/zenodo.22666191)

Code, input data and result files for the manuscript

> M. Sissenbay, I. Orynbassar, "A physics-informed digital twin for tube-life-aware operation of industrial
> steam methane reformers under uncertainty".

The DOI above (10.5281/zenodo.22666191) always resolves to the latest archived version; each release has its
own DOI on the same Zenodo page.

The package `rdt` couples a one-dimensional model of the catalyst tube (Xu–Froment kinetics, Cantera
properties) with a furnace-side radiation model and a Larson–Miller creep-life model, trains surrogates of the
coupled model, and uses them for annual operating scenarios, multi-objective optimization and uncertainty
propagation.

## What the repository contains

```
src/rdt/                       the package
  kinetics.py reactor1d.py furnace.py latham_cases.py      tube and furnace models
  latham_calibration.py validate_xf1989.py                 calibration on plant cases, verification
  creep.py creep_ingest_datasheet.py creep_validate.py     rupture curves
  operate.py surrogate.py scenarios.py optimize.py uq.py   campaigns of the pipeline
  steam_credit.py flame_sweep.py fatigue.py population.py alloy_comparison.py   further analyses
  revision_r1.py revision_r1b.py                           sensitivity and robustness analyses of the revision
  figures.py graphical_abstract.py paper_tables.py         figures and tables (see below)
  pipeline.py                                              entry point
data/
  kinetics/ literature_validation/   numbers transcribed from the literature, with source, page and table
  literature_validation/digitized/   curves digitized from Fig. 3 of Xu & Froment (1989, Part II)
  creep_derived/                     rupture curves digitized from the manufacturer data sheets and from Yeh (2021)
  design/ lhs_runs/ scenarios/ optimization/ uq/ fatigue/ population/   result files (tags v1, v2, v3)
  revision_r1/ revision_r1b/         result files of the revision analyses
notebooks/                     step-by-step notebooks 00–11 (outputs cleared); 03–06 produce the two fit files
paper/figures/v3/              figures, written by the figures stage and by rdt.population
paper/figures/                 graphical abstract (rdt.graphical_abstract) and the schematic
tests/                         test suite
```

The paper uses the result files tagged `v3` (Centralloy G 4852 rupture curve) and those in
`data/revision_r1/` and `data/revision_r1b/`. The files tagged `v1` and `v2` are earlier campaigns on the
Manaurite XM curve; they are kept because the comparison across rupture curves and several tests read them.

Not in the repository: the source articles, the thesis and the manufacturer data sheets (copyright of their
publishers; `data/creep_derived/sources.yaml` and the headers of the YAML files say where each number comes
from), and the trained surrogates (see below). `data/README.md` is an inventory of the data files;
the data are licensed CC BY 4.0 (`data/LICENSE.md`), the code MIT (`LICENSE`).

## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"        # Python >= 3.11; exact versions in requirements-lock.txt
make test                      # or: python -m pytest -q
```

`notebooks/00_sanity_check.ipynb` checks the Cantera installation: equilibrium of a CH4/H2O mixture at 850 °C,
25 bar and a steam-to-carbon ratio of 3 with GRI-Mech 3.0 should give 48.54 % H2, 8.77 % CO, 5.56 % CO2,
3.51 % CH4 and 33.62 % H2O (wet basis), a methane conversion of 80.34 %.

## Trained surrogates

The surrogate models are not stored in the repository (about 400 MB). The scenario stage, the parity figure
and parts of `rdt.revision_r1` and `rdt.revision_r1b` need them. Regenerate them from the stored design with

```bash
python -m rdt.pipeline --stage surrogate --version v3
```

This writes `models/v3/` and rewrites `data/design/surrogate_metrics_v3.json`. It took 561 s (about 9.5 minutes)
on a 10-core Apple M1 Pro (`data/pipeline_v3_runtimes.json`).

## Reproducing the results

All commands write next to the stored files of the same tag; seeds are fixed in the code.

```bash
python -m rdt.pipeline --stage all --version v3     # design, surrogate, scenarios, pareto, uq, figures
make all                                             # the same (VERSION defaults to v3)
python -m rdt.pipeline --stage figures --version v3  # only the figures, from the stored results
```

Wall time of the stages on the machine above: design 233 s, surrogate 561 s, scenarios 182 s, pareto 182 s,
uq 860 s, figures 9 s.

Analyses outside the pipeline, each from the stored v3 results:

```bash
python -m rdt.flame_sweep            # robustness fractions against the flame-length exponent
python -m rdt.fatigue                # thermal-fatigue screening of the scenarios
python -m rdt.population             # tube population
python -m rdt.alloy_comparison       # headline ratios on the four rupture curves
python -m rdt.revision_r1 --calc A B C D E F G H I J K L M    # data/revision_r1/
python -m rdt.revision_r1b           # data/revision_r1b/ (four rupture curves, heat-transfer chain, checks)
python -m rdt.graphical_abstract     # paper/figures/graphical_abstract.*
python -c "from rdt import paper_tables as T; T.make_all()"    # LaTeX tables in paper/tables/
```

## Figures and tables

| Figure file | Tracked at | Generated by |
|---|---|---|
| `fig00_equilibrium_conversion`, `fig01_kinetics`, `fig02_reactor_profiles`, `fig03_validation_xf1989`, `fig04_latham_baseline`, `fig05_latham_calibrated`, `fig06_creep_life_profile` | `paper/figures/v3/` | `python -m rdt.pipeline --stage figures --version v3` (`rdt.figures`, from the models and the input data) |
| `fig07_operating_maps`, `fig08_sobol`, `fig09_surrogate_parity`, `fig10_scenarios_timeseries`, `fig12_pareto`, `fig13_steam_credit`, `fig14_uq_distributions`, `fig15_variance_shares` | `paper/figures/v3/` | the same command (`rdt.figures`, from the stored v3 results; `fig09` needs the trained surrogates) |
| `fig11_scenarios_summary` (with Monte Carlo intervals) | `data/revision_r1/` | `python -m rdt.revision_r1 --calc D` (the version without intervals in `paper/figures/v3/` comes from the figures stage) |
| `fig16_population_damage` | `paper/figures/v3/` | `python -m rdt.population` |
| `figS_flame_length_sweep` | `data/revision_r1/flame_length_sweep.pdf` | `python -m rdt.revision_r1 --calc J` |
| `fig01_schematic` | `paper/figures/` | drawn by hand; no script |
| `graphical_abstract` | `paper/figures/` | `python -m rdt.graphical_abstract` |

`rdt.paper_tables` regenerates the numbers of main Tables 1–4 and of Tables S5 and S6 from the v3 result files
(Table 1 without the Plant C1 row, which is in `data/revision_r1/c1_heldout.csv`). Captions and layout of the
tables in the manuscript differ. The other tables of the supplement are printed from the result files named in
their captions and from the outputs of `rdt.revision_r1` and `rdt.revision_r1b`.

The calibration on the plant cases and the verification against Xu & Froment are in `rdt.latham_calibration`
and `rdt.validate_xf1989`; notebooks 03–06 run them and write `data/literature_validation/latham_fit.yaml`
and `xf1989_fit.yaml`. Re-digitizing the rupture curves (`python -m rdt.creep_ingest_datasheet --all`) needs
the three data sheets named in `data/creep_derived/sources.yaml`, which must be obtained from the manufacturer.
