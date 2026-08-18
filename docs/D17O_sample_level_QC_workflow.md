# IPL Δ′17O sample-level QC/QA workflow

This workflow is implemented in section **8.1 Summarize IPL Δ′17O data** of [`scripts/01_build_CFB_soilcarb_dataset.R`](../scripts/01_build_CFB_soilcarb_dataset.R). It follows the analysis-level QC workflow in [`D17O_analytical_QC_workflow.md`](D17O_analytical_QC_workflow.md).

The purpose is to decide how much confidence to assign to a horizon-level Δ′17O summary after individual analytical failures have been removed. It does not determine whether the carbonate is primary or diagenetically altered; that remains a separate geological screening question.

## Reporting unit

The project reports IPL triple-oxygen data by `MLA_horizon_id`. Each horizon receives one summary row containing its retained mean, replicate information, uncertainty, and sample-level QC status.

The analysis-level table is not overwritten. Excluded analyses remain visible in the raw and analysis-level QC products. The sample-level summary records both the number attempted and the number retained.

## Step 1 — Count attempted and retained analyses

For every horizon, the script records:

- `D17O_n_attempted`: all IPL analyses associated with the horizon before analysis-level exclusions;
- `D17O_n_analysis_excluded`: analyses removed by the section-3 analytical QC decisions;
- `D17O_n_retained`: analyses contributing to the horizon mean;
- `D17O_pending_replicates`: whether a deliberately omitted analysis is awaiting additional replication.

This prevents a singleton retained value from being mistaken for a horizon that was only ever analyzed once.

## Step 2 — Calculate the horizon summary

For retained analyses, the script calculates:

```text
mean Δ′17O
observed replicate SD
observed replicate SE = SD / sqrt(n)
replicate count
```

The same aggregation is retained for carbonate δ′17O and δ′18O because those values are needed downstream for soil-water reconstruction.

## Step 3 — Apply the minimum reproducibility uncertainty

The project uses a minimum Δ′17O reproducibility SD of **12 per meg**. This is a project analytical-reproducibility floor, not a claim that every analysis has exactly 12 per meg uncertainty.

For each horizon:

```text
SD_used = max(observed replicate SD, 12 per meg)
SE_used = SD_used / sqrt(n_retained)
95% half-width = 1.96 × SE_used
```

For a singleton, the observed SD is undefined, so the 12 per meg floor is used directly and the resulting SE is 12 per meg. This avoids assigning zero uncertainty to an unreplicated value.

The output fields are:

- `D17O_sd_used_permeg`
- `D17O_se_used_permeg`
- `IPL17O_sd_Dp17Ocarb_adj`
- `IPL17O_se_Dp17Ocarb_adj`
- `IPL17O_ci95_Dp17Ocarb_adj`

## Step 4 — Classify replicate support

The script uses a 24 per meg maximum replicate-range criterion. This is twice the 12 per meg project reproducibility floor and provides a transparent first-pass test for two or more retained analyses.

| Condition | `D17O_replicate_status` | Interpretation |
|---|---|---|
| No retained analyses | `no_usable_analysis` | No numerical horizon summary |
| One retained analysis | `singleton_provisional` | Usable but unreplicated |
| At least two retained analyses and range ≤24 per meg | `replicated_concordant` | Minimum replicated support for primary use |
| At least two retained analyses and range >24 per meg | `replicated_heterogeneous` | Replicates require review or additional analysis |
| An analysis was deliberately omitted pending replication | `pending_replicates` | Do not use in the locked primary summary |

The replicate range is:

```text
max(retained Δ′17O) − min(retained Δ′17O)
```

It is intentionally reported alongside SD because SD is unstable for a pair of analyses.

## Step 5 — Assign primary versus sensitivity use

The script adds `D17O_primary_use`:

| Status | Primary use |
|---|---|
| `replicated_concordant` | `primary` |
| `singleton_provisional` | `sensitivity_only` |
| `replicated_heterogeneous` | `sensitivity_only` |
| `pending_replicates` | `exclude` |
| `no_usable_analysis` | `exclude` |

This is a confidence classification, not a deletion rule. Singleton and heterogeneous horizons remain available for sensitivity analyses and transparent reporting.

## Why two concordant analyses are the minimum

Requiring three analyses for every horizon would remove a large amount of otherwise usable information. Requiring two concordant analyses provides the minimum direct test of within-horizon repeatability. Three or more analyses are stronger support, but the current implementation does not create a separate three-replicate category because the principal distinction is replicated versus unreplicated.

The 24 per meg range criterion is deliberately simple and auditable. It should be revisited if a larger set of standards or replicate materials establishes a better empirical repeatability distribution.

## Current treatment of heterogeneous replicates

Heterogeneous replicates are not automatically discarded. Their values and range are retained, but the horizon is marked `sensitivity_only`. Before promoting one to primary use, inspect:

1. analysis-level mismatch and precision flags;
2. comments and reactor/session information;
3. whether one analysis is an obvious technical failure;
4. petrographic or phase context, separately from analytical QC; and
5. whether additional replication resolves the disagreement.

The current pending-replicate example is PB-00-02-09L / IPL 5780. Its omission is an explicit replicate-confidence decision, not a diagenetic diagnosis.

## Outputs

The sample-level audit table is:

[`data/processed/CFB_IPL_D17O_sample_level_QC.csv`](../data/processed/CFB_IPL_D17O_sample_level_QC.csv)

The same sample-level fields are joined into the IPL horizon summary and therefore flow into downstream soil-water reconstruction inputs.

## Recommended downstream use

For the primary reconstruction, use horizons with:

```text
D17O_primary_use == "primary"
```

Then repeat the reconstruction with singleton and heterogeneous horizons included as sensitivity cases. Report how many horizons fall into each category and whether the reconstructed soil-water trend changes materially.

## Conceptual separation

This workflow answers:

> Are the retained analyses sufficiently replicated and internally consistent to summarize at the horizon level?

It does not answer:

> Is the carbonate primary, recrystallized, or altered?

That second question belongs to the separate diagenetic-history workflow.
