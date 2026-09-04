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
