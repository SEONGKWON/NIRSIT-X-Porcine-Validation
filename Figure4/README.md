# Figure 4 reproduction package

This folder contains the frozen observation-level data and standalone four-quadrant relative-change analysis used for Figure 4.

## Files

- `figure_input_data.csv`: frozen NIRSIT-X/INVOS observations with source rows and inclusion flags.
- `analyze_four_quadrant_concordance.m`: pseudo-baseline, relative-change, concordance, Wilson CI, and plotting code.
- `figure4_results.csv`: regenerated concordance statistics.
- `figure4_display.csv`: regenerated one-decimal percentage values.
- `figure4_point_data.csv`: every plotted coordinate and pseudo-baseline flag.
- `figure4_reproduced.png`: regenerated two-panel four-quadrant plot.
- `provenance_sha256.txt`: extraction date and archived source/workbook/input hashes.

## Run

From the repository root:

```matlab
run(fullfile('Figure4', 'analyze_four_quadrant_concordance.m'))
```

## Provenance

The package was isolated from the archived internal `Four_axis_relative_3.m` workflow. The original `Pig_AI.xlsx` and `Pig_IN.xlsx` workbooks were used to create `figure_input_data.csv`; neither workbook nor the upstream script is required or included in this flat public repository. Their SHA-256 values are retained in `provenance_sha256.txt`.

## Exact source rules retained

- INVOS source row 4 is excluded.
- A subject's pseudo-baseline is the observation closest to that subject's bivariate median `(median(rSO2), median(SjvO2))`; it is not an experimental Baseline time point.
- Pseudo-baseline points are plotted but excluded from concordance.
- A point is included when `abs(delta SjvO2) > threshold OR abs(delta rSO2) > threshold`.
- Concordance is the fraction with a strictly positive coordinate product; zero products are non-concordant.
- Wilson 95% confidence intervals are reported.
- Plot axes are fixed from -80% to +80%.

This relative-change branch must not be confused with alternative absolute-change or interval-baseline implementations, which produce different concordance results.

## Requirements

- MATLAB R2025b or a compatible release
