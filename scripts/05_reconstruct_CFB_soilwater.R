# 05_reconstruct_CFB_soilwater.R
# Purpose: Reconstruct CFB soil-water d18O and D17O from the integrated
#          carbonate-isotope dataset and the CFB temperature model.

#-- 1.) Setup ---------------------------------------------------------------
library(tidyverse)
library(here)

# ---- Settings ----
n_mc <- 10000
set.seed(123)

#-- 2.) Load CFB Carbonate and Temperature Data ----------------------------
CFB_soilcarb_with_temperature <- read_csv(
  here("data", "processed", "CFB_soilcarb_with_temperature.csv"),
  show_col_types = FALSE
)

if (any(CFB_soilcarb_with_temperature$section_id != "CFB")) {
  stop("Soil-water reconstruction input must contain only section_id == 'CFB'.")
}

# Helper: reconstruct D17Owater after the d18Owater table has been created.
reconstruct_CFB_D17Owater <- function(CFB_d18Ow_reconstruction) {

theta_carb <- 0.5250
lambda_ref <- 0.528

CFB_D17Owater_inputs <- CFB_d18Ow_reconstruction %>%
  filter(
    !is.na(IPL17O_mean_Dp17Ocarb),
    !is.na(IPL17O_mean_d17Ocarb),
    !is.na(IPL17O_mean_d18Ocarb),
    !is.na(T_recon_C),
    !is.na(T_recon_se_C)
  ) %>%
  mutate(
    D17Ocarb_se_used_permeg = pmax(
      case_when(
        !is.na(IPL17O_se_Dp17Ocarb_adj) & IPL17O_se_Dp17Ocarb_adj > 0 ~
          IPL17O_se_Dp17Ocarb_adj,
        !is.na(IPL17O_se_Dp17Ocarb) & IPL17O_se_Dp17Ocarb > 0 ~
          IPL17O_se_Dp17Ocarb,
        TRUE ~ 12
      ),
      12,
      na.rm = TRUE
    ),
    d17Ocarb_se_used = case_when(
      !is.na(IPL17O_se_dp17Ocarb) & IPL17O_se_dp17Ocarb > 0 ~
        IPL17O_se_dp17Ocarb,
      TRUE ~ 0.15
    ),
    d18Ocarb_se_used = case_when(
      !is.na(IPL17O_se_dp18Ocarb) & IPL17O_se_dp18Ocarb > 0 ~
        IPL17O_se_dp18Ocarb,
      TRUE ~ 0.15
    )
  )

# Prepare and reconstruct CFB soil-water D17O.
#
# IMPORTANT UNCERTAINTY TREATMENT
# Delta-prime-17O is a correlated quantity calculated from delta-prime-17O
# and delta-prime-18O. Independently resampling those two isotope ratios
# destroys their analytical covariance and creates spurious uncertainties of
# thousands of per meg. Instead, resample the measured carbonate
# Delta-prime-17O and propagate that quantity directly through fractionation.
# In logarithmic notation:
#
#   D17O_water = D17O_carb - 1e6 * (theta - lambda) * ln(alpha18)
#
# where D17O is in per meg. Carbonate d18O is therefore not an independent
# uncertainty term in this transformation; its contribution is already
# embodied in the measured carbonate D17O uncertainty.

compute_D17Owater <- function(
    D17O_carb_permeg,
    T_C,
    theta = theta_carb,
    lambda = lambda_ref) {
  alpha18 <- calcite_water_alpha18(T_C)

  D17O_carb_permeg -
    1e6 * (theta - lambda) * log(alpha18)
}

mc_D17Owater_one <- function(
    D17O_mean_permeg,
    D17O_se_permeg,
    T_mean,
    T_se,
    n = n_mc) {
  D17Owater_sim <- compute_D17Owater(
    D17O_carb_permeg = rnorm(n, D17O_mean_permeg, D17O_se_permeg),
    T_C = rnorm(n, T_mean, T_se)
  )

  tibble(
    D17Orsw_mean_permeg = mean(D17Owater_sim, na.rm = TRUE),
    D17Orsw_median_permeg = median(D17Owater_sim, na.rm = TRUE),
    D17Orsw_sd_permeg = sd(D17Owater_sim, na.rm = TRUE),
    D17Orsw_lower95_permeg = quantile(D17Owater_sim, 0.025, na.rm = TRUE),
    D17Orsw_upper95_permeg = quantile(D17Owater_sim, 0.975, na.rm = TRUE)
  )
}

CFB_D17Owater_reconstruction <- CFB_D17Owater_inputs %>%
  mutate(
    mc_D17Owater = pmap(
      list(
        IPL17O_mean_Dp17Ocarb,
        D17Ocarb_se_used_permeg,
        T_recon_C,
        T_recon_se_C
      ),
      mc_D17Owater_one
    )
  ) %>%
  unnest(mc_D17Owater) %>%
  mutate(
    D17O_uncertainty_method = paste(
      "Monte Carlo propagation of measured carbonate D17O SE and",
      "temperature SE; correlated d17O-d18O errors preserved by",
      "resampling D17O directly"
    )
  )

# Assemble and export the combined CFB soil-water summary.

CFB_soilwater_reconstruction_summary <- CFB_d18Ow_reconstruction %>%
  left_join(
    CFB_D17Owater_reconstruction %>%
      select(
        section_id,
        MLA_horizon_id,
        strat_height_m,
        D17Ocarb_se_used_permeg,
        D17O_uncertainty_method,
        D17Orsw_mean_permeg,
        D17Orsw_median_permeg,
        D17Orsw_sd_permeg,
        D17Orsw_lower95_permeg,
        D17Orsw_upper95_permeg
      ),
    by = c("section_id", "MLA_horizon_id", "strat_height_m")
  )

write_csv(
  CFB_D17Owater_inputs,
  here("data", "processed", "CFB_D17Owater_inputs.csv")
)

write_csv(
  CFB_D17Owater_reconstruction,
  here("data", "processed", "CFB_D17Owater_reconstruction.csv")
)

write_csv(
  CFB_soilwater_reconstruction_summary,
  here("data", "processed", "CFB_soilwater_reconstruction_summary.csv")
)

invisible(CFB_soilwater_reconstruction_summary)
}

# ---- Add placeholder columns if future U-M data are not present yet ----
if (!"IPL_NuDog_d18Ocarb_VSMOW" %in% names(CFB_soilcarb_with_temperature)) {
  CFB_soilcarb_with_temperature <- CFB_soilcarb_with_temperature %>%
    mutate(IPL_NuDog_d18Ocarb_VSMOW = NA_real_)
}

#-- 3.) Select Carbonate d18O and Temperature Inputs -----------------------
# Priority:
#   Snell samples: use Snell d18Ocarb
#   PB samples: use Bowen d18Ocarb
#   PK95-SC samples: use Koch; if missing, use IPL_NuDog; if missing, use IPL17O

CFB_soilwater_inputs <- CFB_soilcarb_with_temperature %>%
  mutate(
    is_snell = !is.na(Snell_sample_id) & !is.na(Snell_mean_d18Ocarb_vsmow),
    
    d18Ocarb_vsmow = case_when(
      is_snell ~ Snell_mean_d18Ocarb_vsmow,
      str_detect(MLA_horizon_id, "^PB") ~ Bowen_mean_d18Ocarb_vsmow,
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(Koch_mean_d18Ocarb_vsmow) ~ Koch_mean_d18Ocarb_vsmow,
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(IPL_NuDog_d18Ocarb_VSMOW) ~ IPL_NuDog_d18Ocarb_VSMOW,
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(IPL17O_mean_dp18Ocarb) ~ IPL17O_mean_dp18Ocarb,
      TRUE ~ NA_real_
    ),
    
    d18Ocarb_se_vsmow = case_when(
      is_snell ~ NA_real_,
      str_detect(MLA_horizon_id, "^PB") ~ Bowen_se_d18Ocarb_vsmow,
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(Koch_mean_d18Ocarb_vsmow) ~ Koch_se_d18Ocarb_vsmow,
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(IPL_NuDog_d18Ocarb_VSMOW) ~ NA_real_,
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(IPL17O_mean_dp18Ocarb) ~ IPL17O_se_dp18Ocarb,
      TRUE ~ NA_real_
    ),
    
    d18Ocarb_source = case_when(
      is_snell ~ "Snell",
      str_detect(MLA_horizon_id, "^PB") & !is.na(Bowen_mean_d18Ocarb_vsmow) ~ "Bowen",
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(Koch_mean_d18Ocarb_vsmow) ~ "Koch",
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(IPL_NuDog_d18Ocarb_VSMOW) ~ "IPL_NuDog",
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(IPL17O_mean_dp18Ocarb) ~ "IPL17O",
      TRUE ~ NA_character_
    ),
    
    d18Ocarb_missing_flag = case_when(
      str_detect(MLA_horizon_id, "^PB") & is.na(Bowen_mean_d18Ocarb_vsmow) ~ "PB missing Bowen d18Ocarb",
      str_detect(MLA_horizon_id, "^PK95-SC") &
        is.na(Koch_mean_d18Ocarb_vsmow) &
        is.na(IPL_NuDog_d18Ocarb_VSMOW) &
        is.na(IPL17O_mean_dp18Ocarb) ~ "PK95-SC missing Koch, IPL_NuDog, and IPL17O d18Ocarb",
      is.na(d18Ocarb_vsmow) ~ "No selected d18Ocarb",
      TRUE ~ NA_character_
    )
  )

# ---- Choose temperature source ----
# Prioritize measured U-M T47. Otherwise use modeled/interpolated temperature.

CFB_soilwater_inputs <- CFB_soilwater_inputs %>%
  mutate(
    T_recon_C = case_when(
      !is.na(IPLD47_mean_T47_C) ~ IPLD47_mean_T47_C,
      TRUE ~ T_model_C
    ),
    
    T_recon_se_C = case_when(
      !is.na(IPLD47_mean_T47_C) & !is.na(IPLD47_se_T47_C) ~ IPLD47_se_T47_C,
      !is.na(IPLD47_mean_T47_C) & is.na(IPLD47_se_T47_C) ~ T_model_se_C,
      TRUE ~ T_model_se_C
    ),
    
    T_recon_source = case_when(
      !is.na(IPLD47_mean_T47_C) ~ "U-M measured T47",
      !is.na(T_model_C) ~ "modeled/interpolated T",
      TRUE ~ NA_character_
    )
  )

# ---- Fill missing carbonate SE values ----
default_d18Ocarb_se <- 0.15

CFB_soilwater_inputs <- CFB_soilwater_inputs %>%
  mutate(
    d18Ocarb_se_used_vsmow = case_when(
      !is.na(d18Ocarb_se_vsmow) & d18Ocarb_se_vsmow > 0 ~ d18Ocarb_se_vsmow,
      !is.na(d18Ocarb_vsmow) ~ default_d18Ocarb_se,
      TRUE ~ NA_real_
    )
  )

# ---- Calcite-water fractionation: Kim & O'Neil style ----
calcite_water_alpha18 <- function(T_C) {
  T_K <- T_C + 273.15
  exp((18.03 * (1000 / T_K) - 32.42) / 1000)
}

#-- 4.) Reconstruct CFB Soil-Water d18O ------------------------------------
reconstruct_d18Ow <- function(d18Ocarb_vsmow, T_C) {
  alpha18 <- calcite_water_alpha18(T_C)
  ((1000 + d18Ocarb_vsmow) / alpha18) - 1000
}

# ---- Monte Carlo reconstruction for one row ----
mc_d18Ow_one <- function(d18Ocarb_mean, d18Ocarb_se, T_mean, T_se, n = 10000) {
  d18Ocarb_sim <- rnorm(n, mean = d18Ocarb_mean, sd = d18Ocarb_se)
  T_sim <- rnorm(n, mean = T_mean, sd = T_se)
  
  d18Ow_sim <- reconstruct_d18Ow(
    d18Ocarb_vsmow = d18Ocarb_sim,
    T_C = T_sim
  )
  
  tibble(
    d18Ow_mean_vsmow = mean(d18Ow_sim, na.rm = TRUE),
    d18Ow_median_vsmow = median(d18Ow_sim, na.rm = TRUE),
    d18Ow_sd_vsmow = sd(d18Ow_sim, na.rm = TRUE),
    d18Ow_lower95_vsmow = quantile(d18Ow_sim, 0.025, na.rm = TRUE),
    d18Ow_upper95_vsmow = quantile(d18Ow_sim, 0.975, na.rm = TRUE)
  )
}

# ---- Run reconstruction ----
CFB_d18Ow_reconstruction <- CFB_soilwater_inputs %>%
  filter(
    !is.na(d18Ocarb_vsmow),
    !is.na(d18Ocarb_se_used_vsmow),
    !is.na(T_recon_C),
    !is.na(T_recon_se_C)
  ) %>%
  mutate(
    mc = pmap(
      list(
        d18Ocarb_vsmow,
        d18Ocarb_se_used_vsmow,
        T_recon_C,
        T_recon_se_C
      ),
      ~ mc_d18Ow_one(
        d18Ocarb_mean = ..1,
        d18Ocarb_se = ..2,
        T_mean = ..3,
        T_se = ..4,
        n = n_mc
      )
    )
  ) %>%
  unnest(mc)

# ---- QC outputs ----
missing_d18Ocarb <- CFB_soilwater_inputs %>%
  filter(!is.na(d18Ocarb_missing_flag)) %>%
  select(
    MLA_horizon_id,
    strat_height_m,
    d18Ocarb_missing_flag
  )

CFB_soilwater_inputs %>%
  count(d18Ocarb_source)

missing_d18Ocarb

# Quick look plots ------

ggplot(
  CFB_d18Ow_reconstruction,
  aes(x = d18Ow_mean_vsmow, y = strat_height_m)
) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 1500,
    ymax = 1540,
    fill = "grey70",
    alpha = 0.25
  ) +
  geom_errorbarh(
    aes(
      xmin = d18Ow_lower95_vsmow,
      xmax = d18Ow_upper95_vsmow
    ),
    height = 0,
    alpha = 0.5
  ) +
  geom_point(
    aes(shape = d18Ocarb_source, color = T_recon_source),
    size = 2.5,
    alpha = 0.85
  ) +
  annotate(
    "text",
    x = Inf,
    y = 1520,
    label = "PETM",
    hjust = 1.1,
    fontface = "bold"
  ) +
  scale_y_continuous(
    breaks = seq(500, 2300, by = 100)
  ) +
  scale_x_continuous(
    breaks = seq(-12, 2, by = 2)
  ) +
  labs(
    x = expression(delta^18 * O[sw] ~ "(‰ VSMOW)"),
    y = "Stratigraphic height (m)",
    shape = expression(delta^18 * O[carb] ~ "source"),
    color = "Temperature source"
  ) +
  theme_classic(base_size = 18)

ggplot(
  CFB_d18Ow_reconstruction,
  aes(
    x = T_recon_C,
    y = d18Ow_mean_vsmow
  )
) +
  geom_point(
    size = 2,
    alpha = 0.7
  ) +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = expression(delta^18 * O[soil-water] ~ "(‰ VSMOW)")
  ) +
  theme_classic(base_size = 18)

d18Ow_measuredT47 <- CFB_d18Ow_reconstruction %>%
  filter(T_recon_source == "U-M measured T47")

ggplot(
  d18Ow_measuredT47,
  aes(
    x = T_recon_C,
    y = d18Ow_mean_vsmow
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = expression(delta^18 * O[soil-water] ~ "(‰ VSMOW)")
  ) +
  theme_classic(base_size = 18)


#-- 5.) Export CFB d18O-Water Products -------------------------------------
write_csv(
  CFB_soilwater_inputs,
  here("data", "processed", "CFB_soilwater_inputs.csv")
)

write_csv(
  CFB_d18Ow_reconstruction,
  here("data", "processed", "CFB_d18Owater_reconstruction.csv")
)

write_csv(
  missing_d18Ocarb,
  here("data", "processed", "CFB_missing_selected_d18Ocarb.csv")
)


# Temporary evaluation retained for provenance; not part of production.
if (FALSE) {

# TEMPORARY EVAL----------
# ---- Check which horizons fail the d18Ow reconstruction filter ----

recon_filter_check <- CFB_soilwater_inputs %>%
  mutate(
    has_d18Ocarb = !is.na(d18Ocarb_vsmow),
    has_d18Ocarb_se = !is.na(d18Ocarb_se_used_vsmow),
    has_T = !is.na(T_recon_C),
    has_T_se = !is.na(T_recon_se_C),
    
    passes_d18Ow_filter =
      has_d18Ocarb &
      has_d18Ocarb_se &
      has_T &
      has_T_se,
    
    missing_for_filter = pmap_chr(
      list(has_d18Ocarb, has_d18Ocarb_se, has_T, has_T_se),
      function(has_d18Ocarb, has_d18Ocarb_se, has_T, has_T_se) {
        missing <- c(
          if (!has_d18Ocarb) "d18Ocarb",
          if (!has_d18Ocarb_se) "d18Ocarb_se",
          if (!has_T) "T_recon_C",
          if (!has_T_se) "T_recon_se_C"
        )
        paste(missing, collapse = "; ")
      }
    )
  )

# Summary count
recon_filter_check %>%
  count(passes_d18Ow_filter, missing_for_filter)

# Horizon IDs that fail
horizons_failing_d18Ow_filter <- recon_filter_check %>%
  filter(!passes_d18Ow_filter) %>%
  select(
    MLA_horizon_id,
    strat_height_m,
    missing_for_filter,
    d18Ocarb_source,
    d18Ocarb_vsmow,
    d18Ocarb_se_used_vsmow,
    T_recon_source,
    T_recon_C,
    T_recon_se_C
  ) %>%
  arrange(strat_height_m)

horizons_failing_d18Ow_filter

}

#-- 6.) Reconstruct D17Owater and Assemble the Final CFB Summary ------------

CFB_soilwater_reconstruction_summary <- reconstruct_CFB_D17Owater(
  CFB_d18Ow_reconstruction
)
