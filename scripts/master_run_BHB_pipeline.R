# master_run_BHB_pipeline.R
# Run the production Bighorn Basin proxy workflow in dependency order.
#
# The workflow uses explicit CFB primary products and separately processed
# regional-reference records. Plotting scripts run only after all production
# data, chronology, soil-water, and regional-temperature products exist.

library(here)

message("Starting Bighorn Basin proxy-data pipeline...")

# Route plots without an explicit ggsave() target to a null device during the
# non-interactive production run. This prevents exploratory expressions in
# retained legacy blocks from creating an untracked Rplots.pdf; all declared
# pipeline figures are still written by their explicit ggsave() calls.
opened_pipeline_null_device <- FALSE
if (!interactive()) {
  grDevices::pdf(file = NULL)
  opened_pipeline_null_device <- TRUE
}

scripts <- c(
  "01_build_CFB_soilcarb_dataset.R",
  "02_analyze_CFB_carbonate_agreement.R",
  "03_screen_CFB_clumped_diagenesis.R",
  "04_model_CFB_temperatures.R",
  "05_reconstruct_CFB_soilwater.R",
  "06_process_reference_datasets.R",
  "07_build_and_apply_BHB_age_models.R",
  "08_model_BHB_temperatures.R",
  "09_plot_CFB_strat_domain.R",
  "10_plot_BHB_age_domain.R"
)

for (script in scripts) {
  # Patchwork determines annotation and guide dimensions from the active
  # graphics device. The null device is useful while data-processing scripts
  # retain exploratory plot expressions, but it can compress multipanel titles
  # during export. Close it before the explicit publication-plot stages; those
  # scripts write every intended figure through ggsave().
  if (
    script == "09_plot_CFB_strat_domain.R" &&
      opened_pipeline_null_device
  ) {
    grDevices::dev.off()
    opened_pipeline_null_device <- FALSE
  }

  message("\n----------------------------------------")
  message("Running: ", script)
  message("----------------------------------------\n")
  
  source(here("scripts", script))
  
  message("\nFinished: ", script)
}

if (opened_pipeline_null_device && grDevices::dev.cur() > 1) {
  grDevices::dev.off()
}

message("\nBighorn Basin proxy-data pipeline complete.")
