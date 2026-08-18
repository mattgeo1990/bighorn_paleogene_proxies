# IPL Δ′17O analytical QC/QA workflow

This document defines the analytical quality-control workflow for the IPL triple-oxygen dataset. It is implemented in section **3. Clean IPL Δ′17O data** of [`scripts/01_build_CFB_soilcarb_dataset.R`](../scripts/01_build_CFB_soilcarb_dataset.R).

This is deliberately an **analytical** workflow. It asks whether an analysis is technically reliable enough to contribute to a sample summary. It does not decide whether the carbonate is primary, recrystallized, burial-altered, or otherwise useful for paleoenvironmental interpretation. Those questions remain in the separate diagenetic-screening workflow.

## QC decision structure

Each analysis is evaluated in four independent ways:

1. **Mismatch:** Are the mass-33 and mass-34 mismatch diagnostics acceptably small?
2. **Precision:** Are the reported isotope uncertainties reasonable and present?
3. **Value validity:** Is a finite corrected Δ′17O value available?
4. **Replicate context:** Does the analysis agree sufficiently with other analyses of the same horizon/sample, or does it require additional replication?

The output preserves the raw measurements and adds explicit audit fields:

| Field | Meaning |
|---|---|
| `D17O_qc_mismatch_max` | Larger of `X33_mismatch` and `X34_mismatch` |
| `D17O_qc_mismatch_flag` | `pass`, `review`, `exclude`, or `missing` |
| `D17O_qc_precision_flag` | `pass`, `review`, or `missing` based on reported errors |
| `D17O_qc_value_flag` | Whether the accepted corrected Δ′17O is finite |
| `D17O_qc_manual_decision` | Explicit retained/excluded/pending-replicate decision |
| `D17O_qc_decision` | Final analytical disposition used by the section |
| `D17O_qc_reason` | Human-readable reason for an explicit exception |

The analysis-level and summary QC tables are written to:

- `data/processed/CFB_IPL_D17O_analytical_QC_analysis_level.csv`
- `data/processed/CFB_IPL_D17O_analytical_QC_summary.csv`

## Step 1 — Start from corrected IPL values

Section 3 begins with `IPL_D17O_data`, the project table containing corrected carbonate Δ′17O, normalized δ′18O, mismatch diagnostics, and reported analytical errors. The QC does not recalculate the isotope correction; that is handled by the dedicated reactor-correction script and its audit products.

## Step 2 — Evaluate mismatch

The larger of `X33_mismatch` and `X34_mismatch` is used as the conservative mismatch diagnostic:

```r
D17O_qc_mismatch_max = pmax(X33_mismatch, X34_mismatch, na.rm = TRUE)
```

The project review thresholds are:

| Maximum mismatch | Interpretation |
|---:|---|
| `< 0.05` | Pass mismatch screen |
| `0.05–<0.13` | Review in context of precision, comments, and replicates |
| `≥ 0.13` | Hard review flag; normally exclude unless a documented contextual decision retains it |

These are project QC thresholds, not universal IPL instrument specifications. They make the current decisions explicit and reproducible.

## Step 3 — Evaluate reported precision

The script checks `d17O.err`, `d18O.err`, and `CAP17O.err` for missing values. It flags an analysis for precision review when `d17O.err` or `CAP17O.err` is at least 0.010. Precision review is not automatically equivalent to rejection: the value is considered alongside mismatch, run context, comments, and replicate behavior.

This distinction matters because an analysis can have elevated uncertainty but still be informative, while a low-mismatch analysis can still be problematic if the result is poorly constrained or lacks a valid corrected value.

## Step 4 — Verify a finite corrected value

An analysis cannot contribute to a numerical summary if `Dp17Ocarb_permeg_ACCEPTED` is missing or non-finite. Such rows receive `exclude_missing_value` in `D17O_qc_decision`.

## Step 5 — Review replicates without using geology as an analytical filter

The script counts accepted analyses by `MLA_horizon_id` and identifies singletons. Singletons are not automatically rejected; they are reported as less secure because repeatability cannot be assessed.

Replicate disagreement is treated as an analytical-confidence issue only when it is large enough to make the sample summary dependent on a particular analysis. The current explicit example is IPL 5780 from `PB-00-02-09L`, which is temporarily omitted pending additional replication. This is not a claim that one value is geologically altered; it is a statement that the current replicate set does not establish a stable sample value.

## Step 6 — Apply the current explicit decisions

The existing decisions are retained and named in the code so the locked dataset does not change silently:

| IPL analysis | Decision | Reason |
|---:|---|---|
| 5699 | Exclude | High mass-33/34 mismatch |
| 5841 | Exclude | High mismatch and poor analytical precision |
| 6332 | Exclude | Poor analytical precision |
| 6325 | Exclude | High mass-33/34 mismatch |
| 6330 | Retain pending contextual review | High mismatch flag, but retained in the current locked dataset |
| 6344 | Exclude | High mass-33/34 mismatch |
| 5780 | Omit pending replicates | Large within-horizon replicate disagreement |

All other analyses with finite accepted Δ′17O values are retained for the analytical summary stage, including unusual Δ′17O values that lack independent analytical evidence for failure.

## What is intentionally *not* done here

The analytical QC does not:

- reject a value because it is high or low in Δ′17O;
- reject a value because it plots away from a soil-water trend;
- label a carbonate as recrystallized or altered;
- use Δ47 temperature disagreement as an analytical IPL failure criterion;
- use thin-section fabric or EMPA phase maps;
- remove a statistical Δ′17O outlier automatically.

Those are geological or proxy-interpretation questions and belong in the separate diagenetic-history workflow.

## How to use the QC outputs

For each analysis, inspect the QC table in this order:

1. Confirm a finite accepted Δ′17O value.
2. Check `D17O_qc_mismatch_flag` and the two raw mismatch values.
3. Check `D17O_qc_precision_flag` and the reported errors.
4. Read `D17O_qc_manual_decision` and `D17O_qc_reason`.
5. Compare retained replicates within the same horizon/sample.

At the sample-summary stage, report the number of analyses retained and identify singleton or pending-replicate samples. Do not hide exclusions by overwriting the raw table.

## Reproducibility and development history

Matthew Allen and Codex developed this workflow by reviewing the existing IPL correction outputs, manual correction files, reactor-level reconciliation tables, and the current section-3 code. The goal was to preserve the prior locked decisions while converting implicit judgments into named fields and reproducible thresholds.

The implementation was deliberately placed in `01_build_CFB_soilcarb_dataset.R` because this is where the project dataset is assembled and summarized. The separate correction script remains responsible for reconstructing the reactor-level isotope correction itself. The separate diagenetic-screening script remains responsible for geological alteration history. This separation prevents a technically valid but geologically unusual analysis from being mislabeled as an analytical failure.

## Recommended interpretation labels

For downstream tables, the following labels are preferable to a single opaque pass/fail field:

- `retained_for_summary`: technically usable under the current analytical rules;
- `singleton`: retained, but without within-sample replication;
- `review`: elevated mismatch or precision requiring context;
- `omit_pending_replicates`: temporarily excluded because replicate agreement is inadequate;
- `exclude`: analytical evidence supports removal;
- `geological_screening_pending`: analytical QC passed, but primary-versus-altered status has not yet been evaluated.

The key principle is to keep these questions separate: **Was the analysis measured reliably? Are the replicates reproducible? Is the carbonate geologically primary?** This document addresses only the first two.
