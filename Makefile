PY ?= .venv/bin/python
VERSION ?= v2

.PHONY: all design surrogate scenarios pareto uq figures test install lock

install:
	$(PY) -m pip install -e ".[dev]"

lock:
	$(PY) -m pip freeze --exclude-editable > requirements-lock.txt

all:
	$(PY) -m rdt.pipeline --stage all --version $(VERSION)

design surrogate scenarios pareto uq figures:
	$(PY) -m rdt.pipeline --stage $@ --version $(VERSION)

test:
	$(PY) -m pytest -q
