# Source code

This directory contains the model implementation and the analysis scripts used
to generate the figures and numerical results in `docs/report.md`.

## Modules

- `models.py` - core physics code:
  - Sondheimer FS film integral
  - Chambers rectangular-wire surface-scattering integral
  - partial-specularity series for wires
  - Mayadas-Shatzkes grain-boundary model
  - additive and multiplicative FS+MS combination helpers
  - disk-cached interpolation wrapper for the p=0 wire calculation

## Analysis scripts

- `benchmark_fig3.py` - reproduces the Steinhogl Fig. 3 model curves from the
  paper parameter set and writes `figures/benchmark_fig3.png`.
- `fit_steinhogl.py` - fits p and R to the digitized Steinhogl wire data and
  writes `figures/fit_steinhogl.png`.
- `validate_yarimbiyik.py` - evaluates the Yarimbiyik Cu film data using
  literature parameters without refitting and writes
  `figures/validate_yarimbiyik.png`.
- `validate_yarimbiyik_transfer.py` - applies the Steinhogl-fitted R value to
  the Yarimbiyik film dataset and writes the transfer-test figure, CSV, and
  summary text files in `figures/`.
- `sensitivity.py` - performs one-at-a-time sensitivity sweeps for p, R, and
  grain-size scaling and writes `figures/sensitivity.png`.
- `ru_case_study.py` - exploratory Cu-with-liner versus barrierless-Ru
  comparison and writes `figures/ru_case_study.png`.

All scripts should be run from the repository root, for example:

```bash
python src/fit_steinhogl.py
```

The shared model tests are in `tests/test_models.py`.

