# 08b_synthesize_BHB_seasonal_climate.R
# Purpose: Preserve published Bighorn Basin and regional seasonal-temperature
#          evidence as reproducible probability ensembles, construct a clearly
#          labeled relevance-weighted synthesis, and test whether the resolved
#          BHB Delta47 curve aligns with 47 N summer-solstice insolation.
#
# These ensembles are not formal Bayesian posteriors. Published point values,
# reported ranges, and explicit structural-error assumptions are sampled by
# Monte Carlo. The integrated distribution is a relevance-weighted evidence
# mixture; source-specific draws remain available so the judgmental weighting
# can be changed without reconstructing the literature database.

#-- 1.) Setup ---------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)
source(here("scripts", "helpers", "save_figure_variants.R"))

set.seed(47023)
n_draws_per_constraint <- 20000L
n_integrated_draws <- 50000L
n_orbital_null <- 5000L

processed_dir <- here("data", "processed")
figure_dir <- here("figures", "temperature_models", "seasonal_synthesis")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

q_safe <- function(x, p) {
  if (!any(is.finite(x))) return(NA_real_)
  unname(quantile(x, p, na.rm = TRUE, names = FALSE))
}

#-- 2.) Read and Audit Literature Constraints -------------------------------
constraints <- read_csv(
  here("data", "raw", "BHB_seasonal_temperature_constraints.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    metric = factor(metric, levels = c("MAAT", "CMMT", "WMMT", "MART")),
    evidence_label = paste(study, scenario, sep = ": ")
  )

stopifnot(
  all(constraints$distribution %in% c("normal", "uniform")),
  all(is.finite(constraints$estimate_c)),
  all(constraints$primary_weight >= 0)
)

#-- 3.) Generate Source-Specific Probability Ensembles ----------------------
draw_constraint <- function(distribution, estimate, lower, upper, sd, n) {
  if (distribution == "uniform") {
    return(runif(n, min = lower, max = upper))
  }
  rnorm(n, mean = estimate, sd = sd)
}

constraint_draws <- pmap_dfr(
  constraints,
  function(constraint_id, study, year, doi, region, age_context,
           evidence_class, environment, scenario, metric, estimate_c,
           lower_c, upper_c, distribution, sd_c, primary_weight,
           include_primary, uncertainty_basis, notes, evidence_label) {
    tibble(
      constraint_id,
      draw_id = seq_len(n_draws_per_constraint),
      metric = as.character(metric),
      temperature_c = draw_constraint(
        distribution, estimate_c, lower_c, upper_c, sd_c,
        n_draws_per_constraint
      ),
      study, year, doi, region, age_context, evidence_class, environment,
      scenario, primary_weight, include_primary, uncertainty_basis, notes
    )
  }
)

summarize_draws <- function(data, ...) {
  data %>%
    group_by(...) %>%
    summarise(
      n_draws = sum(is.finite(temperature_c)),
      mean_c = mean(temperature_c, na.rm = TRUE),
      sd_c = sd(temperature_c, na.rm = TRUE),
      median_c = q_safe(temperature_c, 0.50),
      lower50_c = q_safe(temperature_c, 0.25),
      upper50_c = q_safe(temperature_c, 0.75),
      lower80_c = q_safe(temperature_c, 0.10),
      upper80_c = q_safe(temperature_c, 0.90),
      lower95_c = q_safe(temperature_c, 0.025),
      upper95_c = q_safe(temperature_c, 0.975),
      .groups = "drop"
    )
}

constraint_summary <- summarize_draws(
  constraint_draws, constraint_id, metric, study, year, doi,
  evidence_class, environment, scenario, uncertainty_basis
)

#-- 4.) Construct the Relevance-Weighted BHB Synthesis ----------------------
# Sampling weights express geographic and proxy relevance, not likelihoods.
# Direct BHB proxy evidence receives the greatest weight; plausible BHB GCM
# states and the dry-soil Green River analog receive lesser weight. The
# intentionally extreme Kiehl orbital-maximum experiment, historical Sewall
# simulation, and wetter-setting leaf MART estimates remain comparisons but
# are excluded from the central synthesis. Marginal metrics are synthesized
# separately, so users must not assume the four draws form internally coupled
# seasonal climates.
primary_constraints <- constraints %>%
  filter(include_primary, primary_weight > 0)

integrated_draws <- map_dfr(
  levels(constraints$metric),
  function(metric_name) {
    candidates <- primary_constraints %>%
      filter(as.character(metric) == metric_name)
    selected <- sample(
      candidates$constraint_id,
      size = n_integrated_draws,
      replace = TRUE,
      prob = candidates$primary_weight
    )
    selected_index <- match(selected, constraints$constraint_id)
    values <- map_dbl(selected_index, function(i) {
      draw_constraint(
        constraints$distribution[i],
        constraints$estimate_c[i],
        constraints$lower_c[i],
        constraints$upper_c[i],
        constraints$sd_c[i],
        1
      )
    })
    tibble(
      synthesis_id = "BHB_dry_interior_relevance_weighted_v1",
      metric = metric_name,
      draw_id = seq_len(n_integrated_draws),
      temperature_c = values,
      sampled_constraint_id = selected
    )
  }
)

integrated_summary <- summarize_draws(
  integrated_draws, synthesis_id, metric
) %>%
  mutate(
    interpretation = paste(
      "Literature-informed relevance-weighted marginal ensemble;",
      "not a formal posterior and not jointly constrained across metrics"
    )
  )

#-- 5.) Test Temperature-Insolation Alignment -------------------------------
temperature_model <- read_csv(
  here("data", "processed", "BHB_D47_temperature_model.csv"),
  show_col_types = FALSE
) %>%
  filter(is.finite(Age_Ma), is.finite(temperature_model_C)) %>%
  arrange(Age_Ma)

insolation <- read_csv(
  here("data", "processed", "BHB_ZB20a_summer_insolation_47N.csv"),
  show_col_types = FALSE
) %>%
  filter(is.finite(Age_Ma), is.finite(summer_solstice_insolation_w_m2)) %>%
  arrange(Age_Ma)

age_min <- max(min(temperature_model$Age_Ma), min(insolation$Age_Ma))
age_max <- min(max(temperature_model$Age_Ma), max(insolation$Age_Ma))
comparison_grid <- tibble(Age_Ma = seq(age_min, age_max, by = 0.001)) %>%
  mutate(
    temperature_c = approx(
      temperature_model$Age_Ma,
      temperature_model$temperature_model_C,
      xout = Age_Ma
    )$y,
    insolation_w_m2 = approx(
      insolation$Age_Ma,
      insolation$summer_solstice_insolation_w_m2,
      xout = Age_Ma
    )$y
  )

lag_results <- tibble(lag_kyr = seq(-100, 100, by = 1)) %>%
  mutate(
    pearson_r = map_dbl(lag_kyr, function(lag) {
      shifted <- approx(
        comparison_grid$Age_Ma,
        comparison_grid$insolation_w_m2,
        xout = comparison_grid$Age_Ma + lag / 1000,
        rule = 1
      )$y
      cor(comparison_grid$temperature_c, shifted, use = "complete.obs")
    }),
    spearman_rho = map_dbl(lag_kyr, function(lag) {
      shifted <- approx(
        comparison_grid$Age_Ma,
        comparison_grid$insolation_w_m2,
        xout = comparison_grid$Age_Ma + lag / 1000,
        rule = 1
      )$y
      cor(
        comparison_grid$temperature_c, shifted,
        use = "complete.obs", method = "spearman"
      )
    })
  )

# Circular shifts preserve each series' autocorrelation while breaking their
# absolute phase relationship. This is a more conservative null comparison
# than treating 1-kyr points as independent observations.
observed_r <- cor(
  comparison_grid$temperature_c,
  comparison_grid$insolation_w_m2,
  use = "complete.obs"
)
null_shifts <- sample(
  seq_len(nrow(comparison_grid) - 1),
  n_orbital_null,
  replace = TRUE
)
null_correlations <- tibble(
  simulation_id = seq_len(n_orbital_null),
  circular_shift_kyr = null_shifts,
  pearson_r = map_dbl(null_shifts, function(k) {
    shifted <- c(
      tail(comparison_grid$insolation_w_m2, k),
      head(comparison_grid$insolation_w_m2, -k)
    )
    cor(comparison_grid$temperature_c, shifted, use = "complete.obs")
  })
)

orbital_test_summary <- tibble(
  test = "zero-lag Pearson correlation against circular-shift null",
  observed = observed_r,
  null_mean = mean(null_correlations$pearson_r),
  null_sd = sd(null_correlations$pearson_r),
  two_sided_p = mean(
    abs(null_correlations$pearson_r) >= abs(observed_r)
  ),
  caution = paste(
    "Exploratory only: astronomical and section age-model phase uncertainty",
    "is not fully propagated; exact precession phase becomes less secure",
    "toward the older part of the record."
  )
)

# Age-jitter sensitivity quantifies how rapidly the zero-lag result degrades
# under plausible-but-not-calibrated age perturbations. These are scenarios,
# not estimates of the actual age uncertainty.
age_jitter_results <- crossing(
  age_jitter_sd_kyr = c(0, 5, 10, 25, 50),
  simulation_id = seq_len(1000)
) %>%
  mutate(
    pearson_r = map2_dbl(
      age_jitter_sd_kyr, simulation_id,
      function(jitter_sd, id) {
        jittered_age <- comparison_grid$Age_Ma +
          rnorm(nrow(comparison_grid), 0, jitter_sd / 1000)
        jittered_insolation <- approx(
          insolation$Age_Ma,
          insolation$summer_solstice_insolation_w_m2,
          xout = jittered_age,
          rule = 1
        )$y
        cor(
          comparison_grid$temperature_c,
          jittered_insolation,
          use = "complete.obs"
        )
      }
    )
  )

#-- 6.) Save Reusable Data Products -----------------------------------------
write_csv(
  constraints %>% select(-evidence_label),
  file.path(processed_dir, "BHB_seasonal_temperature_constraints.csv")
)
write_csv(
  constraint_draws,
  file.path(processed_dir, "BHB_seasonal_temperature_constraint_draws.csv.gz")
)
write_csv(
  constraint_summary,
  file.path(processed_dir, "BHB_seasonal_temperature_constraint_summary.csv")
)
write_csv(
  integrated_draws,
  file.path(processed_dir, "BHB_seasonal_temperature_integrated_draws.csv.gz")
)
write_csv(
  integrated_summary,
  file.path(processed_dir, "BHB_seasonal_temperature_integrated_summary.csv")
)
write_csv(
  lag_results,
  file.path(processed_dir, "BHB_temperature_insolation_lag_analysis.csv")
)
write_csv(
  null_correlations,
  file.path(processed_dir, "BHB_temperature_insolation_null_draws.csv.gz")
)
write_csv(
  orbital_test_summary,
  file.path(processed_dir, "BHB_temperature_insolation_test_summary.csv")
)
write_csv(
  age_jitter_results,
  file.path(processed_dir, "BHB_temperature_insolation_age_jitter.csv.gz")
)

#-- 7.) Plot Literature Distributions and Orbital Diagnostics ---------------
plot_draws <- constraint_draws %>%
  group_by(constraint_id) %>%
  slice_sample(n = 2500) %>%
  ungroup() %>%
  mutate(
    metric = factor(metric, levels = c("MAAT", "CMMT", "WMMT", "MART"))
  )

p_density <- ggplot(
  plot_draws,
  aes(temperature_c, after_stat(density), color = evidence_class)
) +
  geom_density(linewidth = 0.75, adjust = 1.1) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  labs(
    title = "Published and literature-informed seasonal-temperature evidence",
    subtitle = "Curves remain separated by source class; they are not interchangeable observations",
    x = "Temperature or temperature range (deg C)",
    y = "Probability density",
    color = "Evidence class"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

p_forest <- ggplot(
  constraint_summary,
  aes(
    x = median_c,
    y = reorder(paste(study, scenario, sep = " — "), median_c),
    color = evidence_class
  )
) +
  geom_errorbarh(
    aes(xmin = lower95_c, xmax = upper95_c),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(size = 2) +
  facet_wrap(~metric, scales = "free", ncol = 2) +
  labs(
    title = "Central estimates and sampled 95% intervals",
    subtitle = "Intervals combine reported ranges with explicitly documented structural-error assumptions",
    x = "Temperature or temperature range (deg C)",
    y = NULL,
    color = "Evidence class"
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom")

p_integrated <- ggplot(
  integrated_draws %>%
    mutate(metric = factor(metric, levels = c("MAAT", "CMMT", "WMMT", "MART"))) %>%
    slice_sample(n = 50000),
  aes(temperature_c, fill = metric)
) +
  geom_density(alpha = 0.55, color = NA) +
  facet_wrap(~metric, scales = "free_y", ncol = 4) +
  labs(
    title = "Relevance-weighted Bighorn Basin dry-interior synthesis",
    subtitle = "Marginal evidence mixtures, not formal posteriors or a jointly coupled seasonal climate",
    x = "Temperature or temperature range (deg C)",
    y = "Probability density"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

p_lag <- lag_results %>%
  pivot_longer(
    c(pearson_r, spearman_rho),
    names_to = "statistic", values_to = "correlation"
  ) %>%
  ggplot(aes(lag_kyr, correlation, color = statistic)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_line(linewidth = 0.8) +
  labs(
    title = "Lag sensitivity",
    x = "Insolation age shift (kyr)",
    y = "Correlation",
    color = NULL
  ) +
  theme_classic(base_size = 10)

p_null <- ggplot(null_correlations, aes(pearson_r)) +
  geom_histogram(bins = 45, fill = "grey75", color = "white") +
  geom_vline(xintercept = observed_r, color = "#b2182b", linewidth = 1) +
  labs(
    title = "Circular-shift null",
    subtitle = paste0(
      "Observed r = ", round(observed_r, 2),
      "; two-sided p = ", round(orbital_test_summary$two_sided_p, 3)
    ),
    x = "Pearson r",
    y = "Null simulations"
  ) +
  theme_classic(base_size = 10)

p_jitter <- ggplot(
  age_jitter_results,
  aes(factor(age_jitter_sd_kyr), pearson_r)
) +
  geom_boxplot(outlier.shape = NA, fill = "#9ecae1") +
  labs(
    title = "Age-jitter sensitivity",
    subtitle = "Scenario test; jitter values are not formal age errors",
    x = "Age-jitter SD (kyr)",
    y = "Zero-lag Pearson r"
  ) +
  theme_classic(base_size = 10)

save_figure_variants(
  p_density, figure_dir, "BHB_literature_temperature_distributions",
  11, 8, presentation_width = 6
)
save_figure_variants(
  p_forest, figure_dir, "BHB_literature_temperature_forest",
  12, 10, presentation_width = 6
)
save_figure_variants(
  p_integrated, figure_dir, "BHB_primary_temperature_ranges",
  12, 4.5, presentation_width = 12
)
save_figure_variants(
  p_lag + p_null + p_jitter, figure_dir,
  "BHB_orbital_temperature_alignment_diagnostics",
  13, 4.5, presentation_width = 12
)

message("Seasonal-climate synthesis outputs written to: ", processed_dir)
message("Seasonal-climate figures written to: ", figure_dir)
