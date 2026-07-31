source("R/utils.R")
source("R/09_ben_correction.R")

cfg <- read_simple_yaml(here("config", "pipeline.yml"))
pipeline <- run_ben_pipeline(
  input_path = here(cfg$input_file),
  ben_workbook = here(
    "data", "reference", "Clumped March - May 2026 Manual Correction.xlsx"
  ),
  output_dir = here("outputs"),
  minimum_IPL = 4800,
  outlier_threshold = 4
)

print(pipeline$validation$metrics)
message("Ben-equivalent pipeline complete. See outputs/.")
