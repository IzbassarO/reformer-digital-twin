# Reformer Digital Twin

A physics-informed digital twin of an industrial steam methane reformer.
It couples a 1D process model of the catalyst-filled reformer tubes with an
ML surrogate providing uncertainty quantification and a tube creep-life
module, enabling fast, uncertainty-aware prediction of reformer performance
and remaining tube life.

## Folder structure

```
.
├── src/rdt/                    # Python package
│   ├── kinetics.py             # reaction kinetics
│   ├── reactor1d.py            # 1D tube process model
│   ├── furnace.py              # furnace-side heat transfer
│   ├── tube_thermal.py         # tube wall thermal model
│   ├── creep.py                # creep-life module
│   └── surrogate.py            # ML surrogate with UQ
├── data/
│   ├── literature_validation/  # curated literature data for validation
│   ├── creep_derived/          # derived creep-property datasets
│   └── lhs_runs/               # Latin-hypercube sampling runs (large CSVs ignored)
├── notebooks/                  # exploratory notebooks
├── app/                        # interactive application
├── paper/figures/              # manuscript figures
└── tests/                      # test suite
```

## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Sanity check

Running `notebooks/00_sanity_check.ipynb` verifies the Cantera installation
by computing the SMR equilibrium composition (GRI-Mech 3.0) at 850 °C, 25 bar,
S/C = 3 and should reproduce the following output:

```
Equilibrium at 850 °C, 25 bar, S/C = 3 (GRI-Mech 3.0)
Species    Wet basis / %
H2                 48.54
CO                  8.77
CO2                 5.56
CH4                 3.51
H2O                33.62

Dry-basis H2 fraction : 73.13 %
CH4 conversion        : 80.34 %
```

## Manuscript

The journal manuscript (International Journal of Hydrogen Energy, `elsarticle` with `final,3p,times`) lives in
`paper/`: `main.tex` and `supplementary.tex`, both self-contained -- the numbers are inline, the bibliography is
embedded (`thebibliography`, no `.bib`), and the Elsevier highlights are a `highlights` block inside `main.tex`.
Figures are read from `paper/figures/`. Build with

```bash
make paper        # latexmk -pdf for main and supplementary
```

The earlier sectioned draft (`main.tex` + `sections/*.tex` + generated `tables/*.tex` and `numbers.tex`, 11 tables
in the main text, 41 open `\todo{}` markers) is preserved unchanged in `paper/archive/draft_step21/` and can still
be rebuilt with `make paper-draft`, which regenerates its tables and macros from the v2 result files through
`rdt.paper_tables`. A section-by-section comparison of the two is in
`paper/notes/reconcile_overleaf_vs_repo.md`.
