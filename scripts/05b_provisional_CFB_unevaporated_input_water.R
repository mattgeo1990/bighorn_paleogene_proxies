# Provisional sensitivity analysis: back-project reconstructed CFB soil water
# to an unevaporated input-water d18O distribution.
#
# This is intentionally separate from the production pipeline. The preferred
# scenario uses a broad soil-evaporation-slope prior informed by modern soil
# waters and Barnes-Allison model results. Two Passey-Ji lake-polynomial
# scenarios are retained only as sensitivity comparisons.

library(tidyverse)
library(here)

set.seed(20260728)
n_mc <- 20000L
n_saved_per_sample <- 2000L

input_file <- here(
  "data", "processed", "CFB_soilwater_reconstruction_summary.csv"
)
output_dir <- here("data", "processed")
figure_dir <- here("figures", "soilwater_backprojection_provisional")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

draw_truncated_normal <- function(n, mean, sd, lower, upper) {
  out <- rnorm(n, mean, sd)
  bad <- out < lower | out > upper
  while (any(bad)) {
    out[bad] <- rnorm(sum(bad), mean, sd)
    bad <- out < lower | out > upper
  }
  out
}

calcite_water_alpha18 <- function(T_C) {
  exp((18030 / (T_C + 273.15) - 32.42) / 1000)
}

lambda_passey_ji <- function(D17O_water_permil, offset = 0) {
  1.9856 * D17O_water_permil^3 +
    0.5730 * D17O_water_permil^2 +
    0.0601 * D17O_water_permil +
    0.5236 + offset
}

inputs <- read_csv(input_file, show_col_types = FALSE) %>%
  filter(
    IPLD47_n_T47 > 0,
    !is.na(T_measured_C),
    !is.na(IPL17O_mean_dp18Ocarb),
    !is.na(IPL17O_mean_Dp17Ocarb)
  ) %>%
  transmute(
    section_id,
    MLA_horizon_id,
    strat_height_m,
    alteration_likelihood,
    exclude_from_temp_model,
    n_clumped = IPLD47_n_T47,
    n_triple_oxygen = IPL17O_n_analyses,
    T47_C = T_measured_C,
    T47_se_C = coalesce(T_measured_se_C, T_recon_se_C, 2.1),
    dp18Ocarb_permil = IPL17O_mean_dp18Ocarb,
    dp18Ocarb_se_permil = pmax(
      coalesce(IPL17O_se_dp18Ocarb, 0.15), 0.05
    ),
    D17Ocarb_permeg = IPL17O_mean_Dp17Ocarb,
    D17Ocarb_se_permeg = pmax(
      coalesce(IPL17O_se_Dp17Ocarb_adj, IPL17O_se_Dp17Ocarb, 12),
      12
    )
  )

if (nrow(inputs) == 0) {
  stop("No horizons have paired triple-oxygen and measured clumped data.")
}

simulate_one <- function(row) {
  T_draw <- rnorm(n_mc, row$T47_C, row$T47_se_C)
  T_draw[T_draw <= -10] <- NA_real_

  # Preserve the measured carbonate Delta-prime-17O uncertainty directly.
  # The covariance between carbonate D17O and d-prime-18O is not available,
  # so these two measurement terms are provisionally sampled independently.
  dp18Ocarb_draw <- rnorm(
    n_mc, row$dp18Ocarb_permil, row$dp18Ocarb_se_permil
  )
  D17Ocarb_draw_permil <- rnorm(
    n_mc, row$D17Ocarb_permeg, row$D17Ocarb_se_permeg
  ) / 1000

  theta_draw <- rnorm(n_mc, 0.5250, 0.000276)
  alpha18 <- calcite_water_alpha18(T_draw)
  alpha18 <- alpha18 + rnorm(n_mc, 0, 0.0005)

  d18Ocarb_draw <- 1000 * (exp(dp18Ocarb_draw / 1000) - 1)
  d18Osoil_draw <- (1000 + d18Ocarb_draw) / alpha18 - 1000
  dp18Osoil_draw <- 1000 * log1p(d18Osoil_draw / 1000)
  D17Osoil_draw <- D17Ocarb_draw_permil +
    1000 * log(alpha18) * (0.528 - theta_draw)

  D17O_mwl_draw <- rnorm(n_mc, 0.032, 0.005)
  lambda_noise <- rnorm(n_mc, 0, 0.000521)

  scenarios <- list(
    soil_slope_prior = draw_truncated_normal(
      n_mc, mean = 0.5200, sd = 0.0015, lower = 0.518, upper = 0.524
    ),
    passey_ji_lake_published =
      lambda_passey_ji(D17Osoil_draw, offset = 0) + lambda_noise,
    ben_lake_adjusted =
      lambda_passey_ji(D17Osoil_draw, offset = 0.0007) + lambda_noise
  )

  map_dfr(names(scenarios), function(scenario_name) {
    lambda_draw <- scenarios[[scenario_name]]
    dp18Oinput_draw <- dp18Osoil_draw +
      (D17O_mwl_draw - D17Osoil_draw) / (lambda_draw - 0.528)

    valid <- is.finite(dp18Oinput_draw) &
      is.finite(dp18Osoil_draw) &
      lambda_draw < 0.528 &
      dp18Oinput_draw <= dp18Osoil_draw

    tibble(
      scenario = scenario_name,
      draw = seq_len(n_mc),
      dp18Osoil_permil = dp18Osoil_draw,
      D17Osoil_permeg = D17Osoil_draw * 1000,
      lambda_evap = lambda_draw,
      D17O_MWL_permeg = D17O_mwl_draw * 1000,
      dp18Oinput_permil = if_else(valid, dp18Oinput_draw, NA_real_),
      valid = valid
    )
  }) %>%
    mutate(
      section_id = row$section_id,
      MLA_horizon_id = row$MLA_horizon_id,
      strat_height_m = row$strat_height_m,
      alteration_likelihood = row$alteration_likelihood,
      exclude_from_temp_model = row$exclude_from_temp_model,
      .before = 1
    )
}

posterior_draws_full <- map_dfr(
  split(inputs, seq_len(nrow(inputs))),
  ~ simulate_one(.x)
)

posterior_summary <- posterior_draws_full %>%
  group_by(
    section_id, MLA_horizon_id, strat_height_m, alteration_likelihood,
    exclude_from_temp_model, scenario
  ) %>%
  summarise(
    n_trials = n(),
    n_valid = sum(valid),
    fraction_valid = mean(valid),
    dp18Osoil_mean_permil = mean(dp18Osoil_permil, na.rm = TRUE),
    D17Osoil_mean_permeg = mean(D17Osoil_permeg, na.rm = TRUE),
    lambda_mean = mean(lambda_evap[valid], na.rm = TRUE),
    d18Oinput_mean_vsmow = mean(
      1000 * (exp(dp18Oinput_permil[valid] / 1000) - 1), na.rm = TRUE
    ),
    d18Oinput_median_vsmow = median(
      1000 * (exp(dp18Oinput_permil[valid] / 1000) - 1), na.rm = TRUE
    ),
    d18Oinput_sd_vsmow = sd(
      1000 * (exp(dp18Oinput_permil[valid] / 1000) - 1), na.rm = TRUE
    ),
    d18Oinput_lower95_vsmow = quantile(
      1000 * (exp(dp18Oinput_permil[valid] / 1000) - 1),
      0.025, na.rm = TRUE
    ),
    d18Oinput_upper95_vsmow = quantile(
      1000 * (exp(dp18Oinput_permil[valid] / 1000) - 1),
      0.975, na.rm = TRUE
    ),
    .groups = "drop"
  )

# Save an equal number of valid draws per sample so pooled plots give every
# horizon equal weight, regardless of its valid-draw fraction.
posterior_draws_saved <- posterior_draws_full %>%
  filter(valid) %>%
  group_by(MLA_horizon_id, scenario) %>%
  group_modify(
    ~ slice_sample(
      .x,
      n = n_saved_per_sample,
      replace = nrow(.x) < n_saved_per_sample
    )
  ) %>%
  ungroup() %>%
  mutate(
    d18Oinput_vsmow =
      1000 * (exp(dp18Oinput_permil / 1000) - 1),
    d18Osoil_vsmow =
      1000 * (exp(dp18Osoil_permil / 1000) - 1)
  )

overall_summary <- posterior_draws_saved %>%
  group_by(scenario) %>%
  summarise(
    n_samples = n_distinct(MLA_horizon_id),
    d18Oinput_mean_vsmow = mean(d18Oinput_vsmow),
    d18Oinput_median_vsmow = median(d18Oinput_vsmow),
    d18Oinput_sd_vsmow = sd(d18Oinput_vsmow),
    d18Oinput_lower95_vsmow = quantile(d18Oinput_vsmow, 0.025),
    d18Oinput_upper95_vsmow = quantile(d18Oinput_vsmow, 0.975),
    .groups = "drop"
  )

scenario_labels <- c(
  soil_slope_prior = "Soil slope prior",
  passey_ji_lake_published = "Passey-Ji lake polynomial",
  ben_lake_adjusted = "Ben-adjusted lake polynomial"
)

plot_draws <- posterior_draws_saved %>%
  mutate(
    scenario_label = factor(
      recode(scenario, !!!scenario_labels),
      levels = unname(scenario_labels)
    )
  )

p_overall <- ggplot(
  plot_draws,
  aes(d18Oinput_vsmow, color = scenario_label, fill = scenario_label)
) +
  geom_density(alpha = 0.13, linewidth = 0.9, adjust = 1.05) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.4) +
  scale_color_manual(values = c("#1B6CA8", "#8E6BBE", "#C44E52")) +
  scale_fill_manual(values = c("#1B6CA8", "#8E6BBE", "#C44E52")) +
  coord_cartesian(xlim = c(-35, 0)) +
  labs(
    x = expression("Back-projected unevaporated input " * delta^18 *
                     "O (per mil VSMOW)"),
    y = "Probability density",
    color = NULL,
    fill = NULL,
    title = "Provisional distribution of unevaporated input-water d18O",
    subtitle = paste(
      "Equal sample weighting; only horizons with paired triple-oxygen",
      "and measured clumped-isotope data"
    ),
    caption = paste0(
      "Sensitivity analysis, not a production paleoprecipitation estimate.",
      "\n",
      "Soil prior: lambda ~ N(0.520, 0.0015), truncated to 0.518-0.524;",
      " MWL Delta-prime-17O = 32 +/- 5 per meg."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    plot.caption = element_text(hjust = 0, color = "grey30")
  )

p_samples <- posterior_summary %>%
  mutate(
    scenario_label = factor(
      recode(scenario, !!!scenario_labels),
      levels = unname(scenario_labels)
    ),
    MLA_horizon_id = fct_reorder(
      MLA_horizon_id, d18Oinput_median_vsmow, .fun = median
    )
  ) %>%
  ggplot(
    aes(
      d18Oinput_median_vsmow, MLA_horizon_id,
      xmin = d18Oinput_lower95_vsmow,
      xmax = d18Oinput_upper95_vsmow,
      color = scenario_label
    )
  ) +
  geom_errorbarh(
    position = position_dodge(width = 0.55),
    height = 0,
    linewidth = 0.55
  ) +
  geom_point(position = position_dodge(width = 0.55), size = 1.8) +
  scale_color_manual(values = c("#1B6CA8", "#8E6BBE", "#C44E52")) +
  labs(
    x = expression("Back-projected unevaporated input " * delta^18 *
                     "O (per mil VSMOW)"),
    y = "CFB horizon",
    color = NULL,
    title = "Sample-level posterior medians and 95% intervals",
    subtitle = "Intervals include isotope, temperature, fractionation, MWL, and slope uncertainty"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top")

p_validity <- posterior_summary %>%
  mutate(
    scenario_label = factor(
      recode(scenario, !!!scenario_labels),
      levels = unname(scenario_labels)
    )
  ) %>%
  ggplot(aes(fraction_valid, fill = scenario_label)) +
  geom_histogram(binwidth = 0.05, boundary = 0, color = "white") +
  facet_wrap(vars(scenario_label), ncol = 1) +
  scale_fill_manual(values = c("#1B6CA8", "#8E6BBE", "#C44E52")) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x = "Fraction of Monte Carlo draws yielding a physical back-projection",
    y = "Number of horizons",
    title = "Back-projection validity diagnostic"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

write_csv(
  posterior_summary,
  file.path(output_dir, "CFB_unevaporated_input_d18O_provisional_summary.csv")
)
write_csv(
  overall_summary,
  file.path(output_dir, "CFB_unevaporated_input_d18O_provisional_overall.csv")
)
write_csv(
  posterior_draws_saved,
  file.path(output_dir, "CFB_unevaporated_input_d18O_provisional_draws.csv.gz")
)

ggsave(
  file.path(figure_dir, "CFB_unevaporated_input_d18O_overall.png"),
  p_overall, width = 8.2, height = 5.6, dpi = 350
)
ggsave(
  file.path(figure_dir, "CFB_unevaporated_input_d18O_overall.pdf"),
  p_overall, width = 8.2, height = 5.6
)
ggsave(
  file.path(figure_dir, "CFB_unevaporated_input_d18O_by_sample.png"),
  p_samples, width = 8.2, height = max(6, 0.30 * nrow(inputs) + 2.1),
  dpi = 350, limitsize = FALSE
)
ggsave(
  file.path(figure_dir, "CFB_unevaporated_input_d18O_validity.png"),
  p_validity, width = 7.2, height = 7.8, dpi = 350
)

print(overall_summary)
print(
  posterior_summary %>%
    group_by(scenario) %>%
    summarise(
      n_samples = n(),
      median_fraction_valid = median(fraction_valid),
      min_fraction_valid = min(fraction_valid),
      .groups = "drop"
    )
)
