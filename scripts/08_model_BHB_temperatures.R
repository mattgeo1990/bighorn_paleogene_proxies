# 08_model_BHB_temperatures.R
# Purpose: Build two distinct, non-extrapolating regional temperature models:
#          (1) mean annual air temperature (MAT) from Wing et al. LMA data and
#          (2) soil-carbonate formation temperature from all age-resolved BHB
#              Delta47 observations (CFB plus MCP).
#
# The two models are not combined. LMA estimates MAT, whereas pedogenic
# carbonate Delta47 primarily records warm-season soil temperature. Published
# GCM results are retained as climate-state benchmarks, not fitted through-time
# observations: their simulations represent Eocene equilibrium experiments at
# specified CO2 concentrations rather than dated horizons.

#-- 1.) Setup ---------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)
source(here("scripts", "helpers", "save_figure_variants.R"))

processed_dir <- here("data", "processed")
figure_dir <- here("figures", "temperature_models", "regional_BHB")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(47013)
n_simulations <- 2000L
t47_age_bin_width_ma <- 0.10
t47_spline_spar <- 0.45
lma_spline_spar <- 0.65
minimum_simulation_support <- 0.80

quantile_safe <- function(x, probability) {
  if (all(is.na(x))) return(NA_real_)
  unname(quantile(x, probability, na.rm = TRUE, names = FALSE))
}

predict_inside_data <- function(fit, xout, xmin, xmax) {
  prediction <- rep(NA_real_, length(xout))
  inside <- is.finite(xout) & xout >= xmin & xout <= xmax
  prediction[inside] <- predict(fit, x = xout[inside])$y
  prediction
}

#-- 2.) Assemble All Age-Resolved BHB Delta47 Observations -----------------
CFB_temperature_observations <- read_csv(
  here("data", "processed", "CFB_temperature_observations.csv"),
  show_col_types = FALSE
)
CFB_age_lookup <- read_csv(
  here(
    "data", "processed",
    "CFB_soilcarb_isotope_summary_age_calibrated.csv"
  ),
  show_col_types = FALSE
) %>%
  distinct(section_id, MLA_horizon_id, strat_height_m, Age_Ma)

BHB_regional_soilcarb <- read_csv(
  here(
    "data", "processed",
    "BHB_regional_soilcarb_reference_summary_age_calibrated.csv"
  ),
  show_col_types = FALSE
)

BHB_D47_temperature_observations <- bind_rows(
  CFB_temperature_observations %>%
    left_join(
      CFB_age_lookup,
      by = c("section_id", "MLA_horizon_id", "strat_height_m")
    ) %>%
    transmute(
      section_id,
      dataset = paste0(source, " CFB Delta47"),
      study_status = if_else(source == "U-M", "This study", "Published data"),
      MLA_horizon_id,
      Age_Ma,
      temperature_C = T_C,
      temperature_se_C = T_se_C,
      age_uncertainty_status = "not propagated"
    ),
  BHB_regional_soilcarb %>%
    transmute(
      section_id,
      dataset,
      study_status = "Published data",
      MLA_horizon_id,
      Age_Ma,
      temperature_C = T47_C,
      temperature_se_C = T47_se_C,
      age_uncertainty_status = "not propagated"
    )
) %>%
  filter(
    section_id %in% c("CFB", "MCP"),
    is.finite(Age_Ma),
    is.finite(temperature_C)
  )

valid_t47_se <- BHB_D47_temperature_observations$temperature_se_C %>%
  .[is.finite(.) & . > 0]
fallback_t47_se <- median(valid_t47_se)

BHB_D47_temperature_observations <- BHB_D47_temperature_observations %>%
  mutate(
    temperature_se_imputed =
      !is.finite(temperature_se_C) | temperature_se_C <= 0,
    temperature_se_C = if_else(
      temperature_se_imputed, fallback_t47_se, temperature_se_C
    ),
    inverse_variance_weight = 1 / temperature_se_C^2,
    age_bin_ma = round(Age_Ma / t47_age_bin_width_ma) *
      t47_age_bin_width_ma
  )

# Binning before smoothing prevents densely sampled horizons (especially the
# PETM interval) from controlling the regional curve simply because many more
# carbonates were measured there. Within each 0.10 Myr bin, analytical
# estimates are combined with inverse-variance weights. This is a fixed-effect
# mean conditional on reported analytical errors; it does not absorb age-model
# uncertainty, interlaboratory systematics, or ecological/seasonal differences.
BHB_D47_temperature_age_bins <- BHB_D47_temperature_observations %>%
  group_by(age_bin_ma) %>%
  summarise(
    Age_Ma = weighted.mean(
      Age_Ma, inverse_variance_weight, na.rm = TRUE
    ),
    temperature_C = weighted.mean(
      temperature_C, inverse_variance_weight, na.rm = TRUE
    ),
    temperature_se_C = sqrt(
      1 / sum(inverse_variance_weight, na.rm = TRUE)
    ),
    n_observations = n(),
    n_sections = n_distinct(section_id),
    sections = paste(sort(unique(section_id)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(Age_Ma)

if (nrow(BHB_D47_temperature_age_bins) < 4) {
  stop("At least four age bins are required for the BHB Delta47 model.")
}

t47_age_grid <- tibble(
  Age_Ma = seq(
    min(BHB_D47_temperature_age_bins$Age_Ma),
    max(BHB_D47_temperature_age_bins$Age_Ma),
    by = 0.01
  )
)

t47_central_fit <- smooth.spline(
  x = BHB_D47_temperature_age_bins$Age_Ma,
  y = BHB_D47_temperature_age_bins$temperature_C,
  w = 1 / BHB_D47_temperature_age_bins$temperature_se_C^2,
  spar = t47_spline_spar
)

t47_simulations <- map_dfr(seq_len(n_simulations), function(simulation_id) {
  simulated_temperature <- rnorm(
    nrow(BHB_D47_temperature_age_bins),
    mean = BHB_D47_temperature_age_bins$temperature_C,
    sd = BHB_D47_temperature_age_bins$temperature_se_C
  )
  simulated_fit <- smooth.spline(
    x = BHB_D47_temperature_age_bins$Age_Ma,
    y = simulated_temperature,
    w = 1 / BHB_D47_temperature_age_bins$temperature_se_C^2,
    spar = t47_spline_spar
  )
  tibble(
    simulation_id,
    Age_Ma = t47_age_grid$Age_Ma,
    temperature_C = predict_inside_data(
      simulated_fit,
      t47_age_grid$Age_Ma,
      min(BHB_D47_temperature_age_bins$Age_Ma),
      max(BHB_D47_temperature_age_bins$Age_Ma)
    )
  )
})

BHB_D47_temperature_model <- t47_age_grid %>%
  mutate(
    temperature_model_C = predict_inside_data(
      t47_central_fit,
      Age_Ma,
      min(BHB_D47_temperature_age_bins$Age_Ma),
      max(BHB_D47_temperature_age_bins$Age_Ma)
    )
  ) %>%
  left_join(
    t47_simulations %>%
      group_by(Age_Ma) %>%
      summarise(
        temperature_lower95_C = quantile_safe(temperature_C, 0.025),
        temperature_lower80_C = quantile_safe(temperature_C, 0.10),
        temperature_median_C = quantile_safe(temperature_C, 0.50),
        temperature_upper80_C = quantile_safe(temperature_C, 0.90),
        temperature_upper95_C = quantile_safe(temperature_C, 0.975),
        simulation_support_n = sum(is.finite(temperature_C)),
        .groups = "drop"
      ),
    by = "Age_Ma"
  ) %>%
  mutate(
    simulation_support_fraction = simulation_support_n / n_simulations,
    model_target = "soil carbonate formation temperature",
    proxy_basis = "all age-resolved CFB and MCP pedogenic carbonate Delta47",
    uncertainty_scope = paste(
      "Monte Carlo propagation of reported temperature SE at fixed ages;",
      "age-model and proxy-season uncertainties not propagated"
    ),
    extrapolated = FALSE
  )

#-- 3.) Build the Wing LMA Mean-Annual-Air-Temperature Model ---------------
Wing_LMA_temperature_observations <- read_csv(
  here("data", "processed", "WingEtAl2000_LMA_MAT_processed.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    sample,
    sections_with_samples,
    Age_Ma = published_age_model_2_Ma,
    age_duration_ma = duration_Myr,
    age_younger_ma = Age_Ma - age_duration_ma / 2,
    age_older_ma = Age_Ma + age_duration_ma / 2,
    MAT_C,
    MAT_se_C = MAT_error_C,
    age_basis = "Wing et al. aggregate sample interval; published age model 2"
  ) %>%
  filter(is.finite(Age_Ma), is.finite(MAT_C), is.finite(MAT_se_C)) %>%
  arrange(Age_Ma)

if (nrow(Wing_LMA_temperature_observations) < 4) {
  stop("At least four Wing LMA estimates are required for the MAT model.")
}

lma_age_grid <- tibble(
  Age_Ma = seq(
    min(Wing_LMA_temperature_observations$Age_Ma),
    max(Wing_LMA_temperature_observations$Age_Ma),
    by = 0.01
  )
)

lma_central_fit <- smooth.spline(
  x = Wing_LMA_temperature_observations$Age_Ma,
  y = Wing_LMA_temperature_observations$MAT_C,
  w = 1 / Wing_LMA_temperature_observations$MAT_se_C^2,
  spar = lma_spline_spar
)

# Each LMA value pools floras from a finite interval and commonly from several
# sections. The Monte Carlo ensemble therefore samples a uniform age anywhere
# within the published interval and a normal temperature distribution using
# the reported LMA error. This treats the published interval as an age range,
# not a Gaussian analytical error. It is deliberately conservative about
# temporal placement, but it cannot reconstruct within-bin climate structure.
lma_simulations <- map_dfr(seq_len(n_simulations), function(simulation_id) {
  simulated <- Wing_LMA_temperature_observations %>%
    mutate(
      sampled_age_ma = runif(n(), age_younger_ma, age_older_ma),
      sampled_MAT_C = rnorm(n(), MAT_C, MAT_se_C)
    ) %>%
    arrange(sampled_age_ma)

  simulated_fit <- smooth.spline(
    x = simulated$sampled_age_ma,
    y = simulated$sampled_MAT_C,
    w = 1 / simulated$MAT_se_C^2,
    spar = lma_spline_spar
  )

  tibble(
    simulation_id,
    Age_Ma = lma_age_grid$Age_Ma,
    MAT_C = predict_inside_data(
      simulated_fit,
      lma_age_grid$Age_Ma,
      min(simulated$sampled_age_ma),
      max(simulated$sampled_age_ma)
    )
  )
})

BHB_LMA_MAT_model <- lma_age_grid %>%
  mutate(
    MAT_model_C = predict_inside_data(
      lma_central_fit,
      Age_Ma,
      min(Wing_LMA_temperature_observations$Age_Ma),
      max(Wing_LMA_temperature_observations$Age_Ma)
    )
  ) %>%
  left_join(
    lma_simulations %>%
      group_by(Age_Ma) %>%
      summarise(
        MAT_lower95_C = quantile_safe(MAT_C, 0.025),
        MAT_lower80_C = quantile_safe(MAT_C, 0.10),
        MAT_median_C = quantile_safe(MAT_C, 0.50),
        MAT_upper80_C = quantile_safe(MAT_C, 0.90),
        MAT_upper95_C = quantile_safe(MAT_C, 0.975),
        simulation_support_n = sum(is.finite(MAT_C)),
        .groups = "drop"
      ),
    by = "Age_Ma"
  ) %>%
  mutate(
    simulation_support_fraction = simulation_support_n / n_simulations,
    retained_for_plot =
      simulation_support_fraction >= minimum_simulation_support,
    model_target = "mean annual air temperature",
    proxy_basis = "Wing et al. (2000) aggregate leaf-margin analysis",
    uncertainty_scope = paste(
      "Monte Carlo propagation of published MAT error and uniform sampling",
      "within each aggregate age interval"
    ),
    extrapolated = FALSE
  )

#-- 4.) Standardize Published Regional GCM Benchmarks ----------------------
BHB_regional_GCM_temperature_benchmarks <- read_csv(
  here("data", "raw", "BHB_regional_GCM_temperature_benchmarks.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    temperature_metric = factor(
      temperature_metric,
      levels = c("CMMT", "MAT", "WMMT", "MART")
    ),
    comparison_scope = if_else(
      region == "Bighorn Basin",
      "BHB-specific historical regional model",
      "nearby continental-interior regional model"
    )
  )

#-- 5.) Export Reproducible Model Products ---------------------------------
write_csv(
  BHB_D47_temperature_observations,
  file.path(processed_dir, "BHB_D47_temperature_observations.csv")
)
write_csv(
  BHB_D47_temperature_age_bins,
  file.path(processed_dir, "BHB_D47_temperature_age_bins.csv")
)
write_csv(
  BHB_D47_temperature_model,
  file.path(processed_dir, "BHB_D47_temperature_model.csv")
)
write_csv(
  Wing_LMA_temperature_observations,
  file.path(processed_dir, "BHB_LMA_MAT_observations.csv")
)
write_csv(
  BHB_LMA_MAT_model,
  file.path(processed_dir, "BHB_LMA_MAT_model.csv")
)
write_csv(
  BHB_regional_GCM_temperature_benchmarks,
  file.path(processed_dir, "BHB_regional_GCM_temperature_benchmarks.csv")
)

#-- 6.) Diagnostic Figures --------------------------------------------------
temperature_theme <- theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_BHB_LMA_MAT_model <- ggplot() +
  geom_ribbon(
    data = BHB_LMA_MAT_model %>% filter(retained_for_plot),
    aes(Age_Ma, ymin = MAT_lower95_C, ymax = MAT_upper95_C),
    fill = "#1B7837", alpha = 0.16
  ) +
  geom_ribbon(
    data = BHB_LMA_MAT_model %>% filter(retained_for_plot),
    aes(Age_Ma, ymin = MAT_lower80_C, ymax = MAT_upper80_C),
    fill = "#1B7837", alpha = 0.25
  ) +
  geom_line(
    data = BHB_LMA_MAT_model %>% filter(retained_for_plot),
    aes(Age_Ma, MAT_median_C), color = "#1B7837", linewidth = 0.9
  ) +
  geom_segment(
    data = Wing_LMA_temperature_observations,
    aes(
      x = age_younger_ma, xend = age_older_ma,
      y = MAT_C, yend = MAT_C
    ),
    color = "grey35"
  ) +
  geom_errorbar(
    data = Wing_LMA_temperature_observations,
    aes(Age_Ma, ymin = MAT_C - MAT_se_C, ymax = MAT_C + MAT_se_C),
    width = 0, color = "grey35"
  ) +
  geom_point(
    data = Wing_LMA_temperature_observations,
    aes(Age_Ma, MAT_C), shape = 21, fill = "white", size = 2.5
  ) +
  scale_x_reverse() +
  labs(
    x = "Age (Ma)", y = expression("MAT (" * degree * "C)"),
    title = "BHB mean annual air temperature",
    subtitle = "Wing leaf-margin analysis; bands propagate MAT and aggregate-age uncertainty"
  ) +
  temperature_theme

p_BHB_D47_temperature_model <- ggplot() +
  geom_ribbon(
    data = BHB_D47_temperature_model,
    aes(
      Age_Ma,
      ymin = temperature_lower95_C,
      ymax = temperature_upper95_C
    ),
    fill = "#B2182B", alpha = 0.16
  ) +
  geom_ribbon(
    data = BHB_D47_temperature_model,
    aes(
      Age_Ma,
      ymin = temperature_lower80_C,
      ymax = temperature_upper80_C
    ),
    fill = "#B2182B", alpha = 0.25
  ) +
  geom_line(
    data = BHB_D47_temperature_model,
    aes(Age_Ma, temperature_median_C),
    color = "#B2182B", linewidth = 0.9
  ) +
  geom_errorbar(
    data = BHB_D47_temperature_observations,
    aes(
      Age_Ma,
      ymin = temperature_C - temperature_se_C,
      ymax = temperature_C + temperature_se_C,
      color = study_status
    ),
    width = 0, linewidth = 0.3, alpha = 0.45
  ) +
  geom_point(
    data = BHB_D47_temperature_observations,
    aes(Age_Ma, temperature_C, color = study_status, shape = study_status),
    size = 2, alpha = 0.82
  ) +
  scale_color_manual(
    values = c("This study" = "#B2182B", "Published data" = "grey30")
  ) +
  scale_shape_manual(
    values = c("This study" = 21, "Published data" = 22)
  ) +
  scale_x_reverse() +
  labs(
    x = "Age (Ma)",
    y = expression(Delta[47] * " formation temperature (" * degree * "C)"),
    color = NULL, shape = NULL,
    title = "BHB soil-carbonate formation temperature",
    subtitle = paste(
      "U-M measurements from this study versus published CU and Caltech data;",
      "no extrapolation"
    )
  ) +
  temperature_theme

p_regional_GCM_benchmarks <- BHB_regional_GCM_temperature_benchmarks %>%
  ggplot(aes(temperature_metric, estimate_c, color = scenario, shape = region)) +
  geom_errorbar(
    aes(ymin = lower_c, ymax = upper_c),
    width = 0.12, linewidth = 0.6, na.rm = TRUE
  ) +
  geom_point(size = 3, position = position_dodge(width = 0.25)) +
  labs(
    x = NULL, y = expression("Temperature metric (" * degree * "C)"),
    color = "GCM experiment", shape = "Region",
    title = "Regional Eocene GCM benchmarks",
    subtitle = "Equilibrium climate states; points are not assigned to proxy ages"
  ) +
  temperature_theme

p_BHB_integrated_temperature_models <-
  (p_BHB_LMA_MAT_model | p_BHB_D47_temperature_model) /
  p_regional_GCM_benchmarks +
  plot_layout(heights = c(1.25, 0.85)) +
  plot_annotation(
    title = "Bighorn Basin regional temperature models and GCM context",
    subtitle = paste(
      "MAT and carbonate formation temperature are modeled separately;",
      "GCM values are comparison climate states"
    )
  )

save_figure_variants(
  p_BHB_integrated_temperature_models, figure_dir,
  "BHB_integrated_temperature_models_GCM", 12, 9,
  presentation_width = 12
)
