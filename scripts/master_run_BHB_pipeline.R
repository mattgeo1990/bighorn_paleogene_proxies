# master_run_BHB_pipeline.R
# Run full BHB multiproxy workflow in order

library(here)

message("Starting BHB multiproxy pipeline...")

scripts <- c(
  "01_build_BHB_dataset.R",
  "02_temperature_model.R",
  "03_analyze_d18Ocarb_data.R",
  "04_soilwater_d18O_reconstructions.R",
  "05_soilwater_D17O_reconstructions.R",
  "06_process_reference_datasets.R",
  "07_strat_domain_plots.R",
  "08_age_domain_plots.R"
)

for (script in scripts) {
  message("\n----------------------------------------")
  message("Running: ", script)
  message("----------------------------------------\n")
  
  source(here("scripts", script))
  
  message("\nFinished: ", script)
}

message("\nBHB multiproxy pipeline complete.")