# Supplementary Figure 1 reproduction package

This folder contains the frozen observation-level data and standalone Bland–Altman analysis used for Supplementary Figure 1.

## Files

- `figure_input_data.csv`: frozen NIRSIT-X/INVOS observations with source rows and inclusion flags.
- `analyze_bland_altman_agreement.m`: standalone Bland–Altman, subject-cluster bootstrap, and plotting code.
- `supple_figure1_results.csv`: regenerated full-precision results.
- `supple_figure1_display.csv`: regenerated manuscript-style values.
- `supple_figure1_point_data.csv`: mean and difference for every plotted observation.
- `supple_figure1_reproduced.png`: regenerated two-panel figure.
- `provenance_sha256.txt`: extraction date and archived source/workbook/input hashes.

## Run

From the repository root:

```matlab
run(fullfile('SuppleFigure1', 'analyze_bland_altman_agreement.m'))
```

## Provenance

The package was isolated from the Bland–Altman portion of the archived internal `LMM_AI_IN.m` workflow. The original `Pig_AI.xlsx` and `Pig_IN.xlsx` workbooks were used to create `figure_input_data.csv`; neither workbook nor the upstream script is required or included in this flat public repository. Their SHA-256 values are retained in `provenance_sha256.txt`.

INVOS source row 4 is excluded, matching the archived source.

## Critical sign/label finding

The archived workbook headers identify column 2 as rSO2 and column 3 as SjvO2. The source calculates:

```text
difference = rSO2 - SjvO2
```

However, the archived plot label says `SjvO2 - rSO2`, which has the opposite sign. This package preserves the negative source values and current source label to reproduce the existing figure; `supple_figure1_point_data.csv` explicitly names the calculated difference.

Correcting the label or reversing the difference is a manuscript-level decision and reverses the signs of the bias and limits of agreement.

The script uses `rng(2026)`, 2,000 subject-cluster bootstrap samples per device, population SD for limits of agreement, and percentile 95% confidence intervals.

## Requirements

- MATLAB R2025b or a compatible release
- Statistics and Machine Learning Toolbox
