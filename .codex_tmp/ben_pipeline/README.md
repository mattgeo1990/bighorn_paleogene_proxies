# IPL clumped-isotope corrections

Reproducible R workflow for Δ47–Δ48 corrections of Nu Dog clumped-isotope
analyses. The project is deliberately organized as an auditable scientific
workflow rather than a single opaque script.

## Start here

1. Open `IPL_clumped_corrections.Rproj` in RStudio.
2. Run `source("run_pipeline.R")`.
3. Open `reports/correction_workflow.html`.
4. Inspect `outputs/IPL_clumped_corrections.xlsx`.

The Excel output is generated from the original manual-correction workbook.
Its original `Manual Correction` sheet is preserved, while `data_cleaned`,
`Parameters`, `Calibration`, and `Full Correction` provide the full cleaned
dataset and a live, traceable Excel formula chain. Changing an accepted value
or correction parameter recalculates the downstream workbook columns.

The input snapshot is in `data/raw/`. Accepted reference values and model
choices are visible in `config/`. Update those files before treating the
results as publication-ready.

## Correction sequence

The active pipeline reproduces Ben's manual-workbook Δ47 procedure and applies
it to included analyses with IPL numbers of 4800 or greater:

1. Read the linearity slope, CDES transfer slope/intercept, and carbonate
   offset directly from Ben's controlling workbook cells.
2. Project every analysis to δ47 = 0 with Ben's stored slope.
3. Apply Ben's stored CDES transfer and carbonate offset.
4. Calculate residuals for 102-GC-AZ01, ETH-1 through ETH-4, and IAEA-C1.
5. Mark extreme within-material standard residuals as excluded using a 4 scaled-
   MAD threshold. The values remain in the audit trail, while
   `correction_qc_include` is false for downstream filtering.
7. Calculate Petersen temperatures with Ben's explicit +0.088 acid-
   fractionation term.
7. Re-run the same implementation on Ben's March-May subset and compare every
   overlapping reported result.

Outputs include the corrected run-level CSV, fitted parameters, a standard
outlier audit, the row-level overlap comparison, and compact validation metrics.

## Important status

The older experimental projection/transfer/drift scripts remain in `R/` for
provenance, but `run_pipeline.R` now calls `R/09_ben_correction.R`.
