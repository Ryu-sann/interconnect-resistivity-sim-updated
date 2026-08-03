# Cu interconnect size-effect study

This repository contains code, digitized data, generated figures, and a
technical report for a summer research project on size-dependent resistivity
in nanoscale copper interconnects.

The project implements the Fuchs-Sondheimer (FS) surface-scattering model and
the Mayadas-Shatzkes (MS) grain-boundary-scattering model from the original
equations, checks the implementation against published limits/tables/curves,
fits the Steinhogl Cu wire dataset, validates against the Yarimbiyik Cu film
dataset, and tests whether a grain-boundary reflection coefficient fitted from
one dataset transfers to the other.

## Main files

- `docs/report.md` - full technical report with figures embedded
- `docs/PARAMETERS.md` - parameter values, assumptions, and sources
- `docs/references.md` - source papers, DOI links, and data provenance notes
- `docs/AI_USE.md` - brief disclosure of AI-assisted work
- `src/` - model and analysis scripts
- `data/` - digitized or transcribed literature data used by the scripts
- `figures/` - generated figures used in the report
- `tests/test_models.py` - physics/regression tests for the model code

## Repository structure

```text
.
|-- README.md
|-- requirements.txt
|-- data/
|   |-- steinhogl2002_fig3_points.csv
|   |-- steinhogl2002_fig3_curves.csv
|   |-- yarimbiyik2009_table1_films.csv
|   `-- yarimbiyik2009_table2_grains.csv
|-- docs/
|   |-- report.md
|   |-- PARAMETERS.md
|   |-- references.md
|   |-- AI_USE.md
|   `-- archive/
|       `-- report_before_transfer_revision.md
|-- figures/
|   |-- benchmark_fig3.png
|   |-- fit_steinhogl.png
|   |-- validate_yarimbiyik.png
|   |-- validate_yarimbiyik_transfer.png
|   |-- sensitivity.png
|   `-- ru_case_study.png
|-- notebooks/
|   `-- walkthrough.ipynb
|-- src/
|   |-- models.py
|   |-- benchmark_fig3.py
|   |-- fit_steinhogl.py
|   |-- validate_yarimbiyik.py
|   |-- validate_yarimbiyik_transfer.py
|   |-- sensitivity.py
|   `-- ru_case_study.py
`-- tests/
    `-- test_models.py
```

## Setup

Use Python 3.10 or newer.

```bash
pip install -r requirements.txt
```

## Reproducing the analysis

Run commands from the repository root. The first script that needs the
Chambers wire interpolation table will create `data/cache/`; that directory
is generated locally and is not committed.

```bash
python src/benchmark_fig3.py
python src/fit_steinhogl.py
python src/validate_yarimbiyik.py
python src/validate_yarimbiyik_transfer.py
python src/sensitivity.py
python src/ru_case_study.py
```

Expected generated outputs:

- `figures/benchmark_fig3.png`
- `figures/fit_steinhogl.png`
- `figures/validate_yarimbiyik.png`
- `figures/validate_yarimbiyik_transfer.png`
- `figures/validate_yarimbiyik_transfer.csv`
- `figures/validate_yarimbiyik_transfer_output.txt`
- `figures/sensitivity.png`
- `figures/ru_case_study.png`

Run the test suite with:

```bash
python -m pytest tests/ -q
```

The tests check published Sondheimer table values, analytic limiting cases,
FS/MS consistency checks, and regression behavior for the combined model.

## Data provenance

The data files in `data/` were extracted from published experimental papers.
The Steinhogl wire data were digitized from Fig. 3. The Yarimbiyik film
resistivity and grain-size data were transcribed from Tables 1 and 2. The
original articles are cited in `docs/references.md`; copyrighted article PDFs
are not committed to this public repository.

## Notes for review

The core scientific claims are in `docs/report.md`. The most important review
targets are:

- the FS/MS implementations in `src/models.py`
- the Steinhogl fit in `src/fit_steinhogl.py`
- the Yarimbiyik no-refit validation in `src/validate_yarimbiyik.py`
- the cross-dataset transfer test in `src/validate_yarimbiyik_transfer.py`
- the model tests in `tests/test_models.py`

