# Bighorn Basin Paleogene Proxies

This repository integrates Paleocene–Eocene proxy records from the Bighorn
Basin, with a primary CFB soil-carbonate workflow, section-specific age models,
regional reference datasets, uncertainty propagation, and publication-ready
stratigraphic- and age-domain figures.

## Project guide

The complete computational and interpretive documentation is here:

- [Data pipeline, chronostratigraphic framework, and interpretive guide](docs/BHB_DATA_PIPELINE_AND_INTERPRETIVE_FRAMEWORK.md)

It explains the data model, production scripts, age priors, section framework,
temperature and soil-water models, proxy interpretation, GCM comparisons,
uncertainty limitations, and recommended downstream analyses.

## Run the production pipeline

From the project root in R:

```r
source(here::here("scripts", "master_run_BHB_pipeline.R"))
```

The master script runs the numbered production scripts in dependency order and
writes processed data to `data/processed/` and figures to `figures/`.

