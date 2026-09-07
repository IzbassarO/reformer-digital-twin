PY ?= .venv/bin/python
VERSION ?= v2

.PHONY: all design surrogate scenarios pareto uq figures test install lock paper paper-draft

# manuscript: the consolidated Overleaf pair in paper/ (numbers inline, bibliography embedded,
# highlights inside main.tex). The archived sectioned draft is rebuilt by `make paper-draft`.
paper:
	cd paper && latexmk -pdf -interaction=nonstopmode main.tex && latexmk -pdf -interaction=nonstopmode supplementary.tex && latexmk -c

paper-draft:
	$(PY) -c "from rdt import paper_tables as T; T.make_all()"
	cd paper/archive/draft_step21 && latexmk -pdf -interaction=nonstopmode main.tex && latexmk -pdf -interaction=nonstopmode supplementary.tex && latexmk -pdf -interaction=nonstopmode highlights.tex && latexmk -c

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
