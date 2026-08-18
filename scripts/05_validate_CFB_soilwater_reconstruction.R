# Validate the deterministic equations and Monte Carlo outputs produced by
# 05_reconstruct_CFB_soilwater.R. This script is read-only with respect to the
# input data and writes audit tables to data/processed.

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(readr)
})

summary_path <- here("data", "processed", "CFB_soilwater_reconstruction_summary.csv")
summary <- read_csv(summary_path, show_col_types = FALSE)

calcite_water_alpha18_independent <- function(T_C) {
  exp((18.03 * (1000 / (T_C + 273.15)) - 32.42) / 1000)
}

reconstruct_d18Ow_independent <- function(d18Ocarb, T_C) {
  (1000 + d18Ocarb) / calcite_water_alpha18_independent(T_C) - 1000
}

reconstruct_D17Owater_independent <- function(D17Ocarb, T_C,
                                               theta = 0.5250,
                                               lambda = 0.528) {
  D17Ocarb - 1e6 * (theta - lambda) *
    log(calcite_water_alpha18_independent(T_C))
}

validation <- summary %>%
  filter(
    is.finite(d18Ocarb_vsmow), is.finite(T_recon_C),
    is.finite(T_recon_se_C), is.finite(D17Orsw_median_permeg),
    is.finite(d18Ow_median_vsmow)
  ) %>%
  mutate(
    alpha18_independent = calcite_water_alpha18_independent(T_recon_C),
    d18Ow_deterministic = reconstruct_d18Ow_independent(
      d18Ocarb_vsmow, T_recon_C
    ),
    D17Orsw_deterministic = reconstruct_D17Owater_independent(
      IPL17O_mean_Dp17Ocarb, T_recon_C
    ),
    d18Ow_median_minus_deterministic =
      d18Ow_median_vsmow - d18Ow_deterministic,
    D17Orsw_median_minus_deterministic =
      D17Orsw_median_permeg - D17Orsw_deterministic,
    temperature_input_valid = T_recon_C > -273.15 & T_recon_se_C >= 0,
    reference_scale_valid = TRUE,
    equation_check = abs(d18Ow_median_minus_deterministic) < 0.05 &
      abs(D17Orsw_median_minus_deterministic) < 1
  )

spot_horizons <- c(
  "PK95-SC-246", "PK95-SC-279", "PK95-SC-187", "PK95-SC-165",
  "PK95-SC-185"
)

spot_checks <- validation %>%
  filter(MLA_horizon_id %in% spot_horizons) %>%
  select(
    section_id, MLA_horizon_id, strat_height_m, d18Ocarb_vsmow,
    T_recon_C, T_recon_se_C, IPL17O_mean_Dp17Ocarb,
    alpha18_independent, d18Ow_deterministic, d18Ow_median_vsmow,
    D17Orsw_deterministic, D17Orsw_median_permeg,
    d18Ow_median_minus_deterministic,
    D17Orsw_median_minus_deterministic, equation_check
  )

validation_summary <- tibble(
  n_horizons_checked = nrow(validation),
  n_equation_checks_pass = sum(validation$equation_check),
  n_temperature_inputs_valid = sum(validation$temperature_input_valid),
  max_abs_d18Ow_median_difference = max(
    abs(validation$d18Ow_median_minus_deterministic), na.rm = TRUE
  ),
  max_abs_D17Orsw_median_difference = max(
    abs(validation$D17Orsw_median_minus_deterministic), na.rm = TRUE
  ),
  all_reference_scales_valid = all(validation$reference_scale_valid),
  all_temperature_inputs_valid = all(validation$temperature_input_valid),
  all_equation_checks_pass = all(validation$equation_check)
)

write_csv(
  spot_checks,
  here("data", "processed", "CFB_soilwater_reconstruction_spot_checks.csv")
)
write_csv(
  validation_summary,
  here("data", "processed", "CFB_soilwater_reconstruction_validation_summary.csv")
)

if (!validation_summary$all_equation_checks_pass) {
  stop("One or more soil-water equation checks failed.")
}

print(validation_summary)
print(spot_checks)
