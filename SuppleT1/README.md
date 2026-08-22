# Supplementary Table 1 reproduction package

This folder contains the extracted subject-level physiology stage values and standalone MATLAB code needed to reproduce Supplementary Table 1.

## Files

- `physiology_stage_values.csv`: subject-level MBP, SBP, DBP, SpO2, HR, and EtCO2 values at Baseline, T1, T2, and T3.
- `analyze_physiology_wilcoxon_holm.m`: paired Wilcoxon signed-rank tests, Holm–Bonferroni correction, Z values, and effect-size r.
- `supple_t1_results.csv`: full-precision output containing all four comparisons in each Holm family.
- `supple_t1_display.csv`: three-decimal output containing the three Baseline comparisons displayed in the supplementary table.

## Run

From the repository root:

```matlab
run(fullfile('SuppleT1', 'analyze_physiology_wilcoxon_holm.m'))
```


## Provenance

The stage values were extracted on August 21, 2026 using MATLAB R2025b from the archived internal source `Pig_Reulst_rev2.m`. That upstream script is not required and is not included in this flat public repository.

SHA-256 of the source script at extraction:

```text
00505FEA8BFE0A20EFDCC781F59925D7986B88184FA1766779BB32668A490EFD
```

T1, T2, and T3 correspond to `after_clamping`, `before_open`, and `after_open`.

## Important behavior preserved

1. Descriptive medians and IQRs are recalculated from complete pairs for each comparison.
2. Missing values are excluded pairwise.
3. Zero paired differences remain in descriptive statistics but are removed before the Wilcoxon test and effect-size calculation.
4. `N_complete` is used for paired descriptive values.
5. `N_nonzero` is the effective N for the test and for `r = abs(Z) / sqrt(N_nonzero)`.
6. The Holm family contains four comparisons per variable: Baseline–T1, Baseline–T2, Baseline–T3, and T1–T2.
7. The display table contains the first three comparisons while retaining adjusted p values calculated from all four.

Pig 11 has a T2 value but no Baseline value. It contributes to the manuscript Table 1 T2 summary but is excluded from every Baseline–T2 pair in Supplementary Table 1.

The calculation-derived three-decimal HR values are 135.133 and 134.167. Values such as 135.130 and 134.170 result from padding previously rounded two-decimal values rather than rounding the full-precision values.

## Requirements

- MATLAB R2025b or a compatible release
- Statistics and Machine Learning Toolbox
