# Figure 3 reproduction package

This folder contains the frozen observation-level data and standalone linear mixed-effects model (LMM) analysis used for Figure 3.

## Files

- `figure_input_data.csv`: frozen NIRSIT-X/INVOS observations with source rows and inclusion flags.
- `fit_rso2_sjvo2_lmm.m`: standalone LMM estimation and plotting code.
- `figure3_results.csv`: regenerated full-precision results.
- `figure3_display.csv`: regenerated manuscript-style values.
- `figure3_reproduced.png`: regenerated two-panel LMM figure.
- `provenance_sha256.txt`: extraction date and archived source/workbook/input hashes.

## Run

From the repository root:

```matlab
run(fullfile('Figure3', 'fit_rso2_sjvo2_lmm.m'))
```

The script fits one random-intercept LMM per device and regenerates the result tables and two-panel figure.

## Provenance

The package was isolated from the archived internal `LMM_AI_IN.m` workflow. The original `Pig_AI.xlsx` and `Pig_IN.xlsx` workbooks were used to create `figure_input_data.csv`; neither workbook nor the upstream script is required or included in this flat public repository. Their SHA-256 values are retained in `provenance_sha256.txt`.

INVOS source row 4 is excluded, matching the archived source.

## Critical column-orientation finding

The archived workbook headers identify column 2 as NIRSIT-X/INVOS rSO2 and column 3 as SjvO2. The original source assigns column 2 to `Y`, column 3 to `z`, fits `Y ~ z`, and labels the plotted equation and axes as if `Y` were SjvO2 and `z` were rSO2.

This package intentionally preserves the numerical calculation:

```text
rSO2 ~ SjvO2 + (1 | Subject)
```

It also preserves the current source plot labels so the regenerated figure follows the archived implementation. Reversing the model direction is a distinct analysis and must be an explicit manuscript decision.

Marginal and conditional R-squared analyses and supporting materials are not included in this public release. They are available from the corresponding author upon reasonable request.

## Requirements

- MATLAB R2025b or a compatible release
- Statistics and Machine Learning Toolbox
