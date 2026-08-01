# 16_plot_Harper_Tornillo_BHB_temperature_synthesis.R
# Purpose: Three-panel age-domain comparison of Harper et al. Pacific SST,
#          screened Kelson et al. Tornillo soil-carbonate temperatures, and
#          the Bighorn Basin multiproxy seasonal-temperature synthesis.
#
# Tornillo policy:
#   Only rows explicitly flagged used_in_tornillo_temperature_model by
#   06_process_reference_datasets.R enter the model or plot. This follows the
#   sample exclusions in Kelson et al. (2018) and preserves the project's
#   documented legacy paired-nodule exclusion.

#-- 1.) Setup ---------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)
source(here("scripts", "helpers", "save_figure_variants.R"))

set.seed(4716)
n_simulations <- 2000L
tornillo_age_bin_width_ma <- 0.05
tornillo_spline_spar <- 0.35
age_limits <- c(59.0, 52.5)
petm_old_ma <- 55.93
petm_young_ma <- 55.75

processed_dir <- here("data", "processed")
figure_dir <- here(
  "figures", "age_domain", "regional_temperature_synthesis"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

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

add_petm <- function(alpha = 0.08) {
  annotate(
    "rect", xmin = -Inf, xmax = Inf,
    ymin = petm_young_ma, ymax = petm_old_ma,
    fill = "#D73027", alpha = alpha
  )
}

shared_age_scale <- scale_y_reverse(
  limits = age_limits,
  breaks = seq(59, 52.5, by = -0.5),
  minor_breaks = seq(59, 52.5, by = -0.1),
  expand = expansion(mult = c(0.008, 0.012))
)

theme_panel <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5),
    axis.title = element_text(size = 10.5),
    axis.text = element_text(size = 8.5, color = "black"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8.5),
    plot.margin = margin(4, 5, 4, 4)
  )

#-- 2.) Harper Pacific SST -------------------------------------------------
Harper_SST <- read_csv(
  here("data", "processed", "Harper2024_CO2_SST_processed.csv"),
  show_col_types = FALSE
) %>%
  filter(
    between(Age_Ma, min(age_limits), max(age_limits)),
    is.finite(Harper2024_mean_SST_C)
  ) %>%
  arrange(desc(Age_Ma))

#-- 3.) Screen and Model Tornillo Temperatures -----------------------------
Tornillo_inventory <- read_csv(
  here("data", "processed", "Kelson_Tornillo_D47_processed.csv"),
  show_col_types = FALSE
)

required_screen_columns <- c(
  "used_in_tornillo_temperature_model",
  "kelson_exclusion_reason",
  "kelson_screening_basis"
)
if (!all(required_screen_columns %in% names(Tornillo_inventory))) {
  stop(
    "Kelson screening flags are missing. Run ",
    "06_process_reference_datasets.R before this stage."
  )
}

Tornillo_observations <- Tornillo_inventory %>%
  filter(
    used_in_tornillo_temperature_model,
    between(Age_Ma, min(age_limits), max(age_limits)),
    is.finite(T47_C)
  ) %>%
  mutate(
    T47_se_imputed = !is.finite(T47_se_C) | T47_se_C <= 0,
    T47_se_C = if_else(
      T47_se_imputed,
      median(
        T47_se_C[is.finite(T47_se_C) & T47_se_C > 0],
        na.rm = TRUE
      ),
      T47_se_C
    ),
    inverse_variance_weight = 1 / T47_se_C^2,
    age_bin_ma =
      round(Age_Ma / tornillo_age_bin_width_ma) *
      tornillo_age_bin_width_ma
  )

if (any(!Tornillo_observations$used_in_tornillo_temperature_model)) {
  stop("An excluded Tornillo observation entered the temperature model.")
}

Tornillo_age_bins <- Tornillo_observations %>%
  group_by(age_bin_ma) %>%
  summarise(
    Age_Ma = weighted.mean(Age_Ma, inverse_variance_weight),
    temperature_C = weighted.mean(T47_C, inverse_variance_weight),
    temperature_se_C = sqrt(1 / sum(inverse_variance_weight)),
    n_observations = n(),
    sample_ids = paste(sort(sample_id), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(Age_Ma)

if (nrow(Tornillo_age_bins) < 4) {
  stop("At least four screened Tornillo age bins are required.")
}

Tornillo_age_grid <- tibble(
  Age_Ma = seq(
    min(Tornillo_age_bins$Age_Ma),
    max(Tornillo_age_bins$Age_Ma),
    by = 0.01
  )
)

Tornillo_central_fit <- smooth.spline(
  x = Tornillo_age_bins$Age_Ma,
  y = Tornillo_age_bins$temperature_C,
  w = 1 / Tornillo_age_bins$temperature_se_C^2,
  spar = tornillo_spline_spar
)

Tornillo_simulations <- map_dfr(seq_len(n_simulations), function(i) {
  simulated_fit <- smooth.spline(
    x = Tornillo_age_bins$Age_Ma,
    y = rnorm(
      nrow(Tornillo_age_bins),
      Tornillo_age_bins$temperature_C,
      Tornillo_age_bins$temperature_se_C
    ),
    w = 1 / Tornillo_age_bins$temperature_se_C^2,
    spar = tornillo_spline_spar
  )
  tibble(
    simulation_id = i,
    Age_Ma = Tornillo_age_grid$Age_Ma,
    temperature_C = predict_inside_data(
      simulated_fit,
      Tornillo_age_grid$Age_Ma,
      min(Tornillo_age_bins$Age_Ma),
      max(Tornillo_age_bins$Age_Ma)
    )
  )
})

Tornillo_temperature_model <- Tornillo_age_grid %>%
  mutate(
    temperature_model_C = predict_inside_data(
      Tornillo_central_fit,
      Age_Ma,
      min(Tornillo_age_bins$Age_Ma),
      max(Tornillo_age_bins$Age_Ma)
    )
  ) %>%
  left_join(
    Tornillo_simulations %>%
      group_by(Age_Ma) %>%
      summarise(
        temperature_lower95_C = quantile_safe(temperature_C, 0.025),
        temperature_lower80_C = quantile_safe(temperature_C, 0.10),
        temperature_median_C = quantile_safe(temperature_C, 0.50),
        temperature_upper80_C = quantile_safe(temperature_C, 0.90),
        temperature_upper95_C = quantile_safe(temperature_C, 0.975),
        .groups = "drop"
      ),
    by = "Age_Ma"
  ) %>%
  mutate(
    model_target = "Tornillo soil-carbonate formation temperature",
    screening_rule = paste(
      "Kelson et al. (2018) paleoclimate exclusions plus legacy project",
      "paired-nodule exclusion; primary micrite only"
    ),
    uncertainty_scope = paste(
      "Monte Carlo propagation of reported temperature SE at fixed ages;",
      "age-model and seasonal uncertainties not propagated"
    ),
    spline_spar = tornillo_spline_spar,
    extrapolated = FALSE
  )

write_csv(
  Tornillo_observations,
  file.path(processed_dir, "Tornillo_temperature_model_observations.csv"),
  na = ""
)
write_csv(
  Tornillo_age_bins,
  file.path(processed_dir, "Tornillo_temperature_model_age_bins.csv"),
  na = ""
)
write_csv(
  Tornillo_temperature_model,
  file.path(processed_dir, "Tornillo_temperature_model.csv"),
  na = ""
)

#-- 4.) BHB Multiproxy Synthesis ------------------------------------------
BHB_D47 <- read_csv(
  here("data", "processed", "BHB_D47_temperature_observations.csv"),
  show_col_types = FALSE
) %>%
  filter(
    used_in_temperature_model,
    between(Age_Ma, min(age_limits), max(age_limits))
  ) %>%
  mutate(
    record = if_else(
      study_status == "This study", "This study D47", "Published D47"
    )
  )

Wing_LMA <- read_csv(
  here("data", "processed", "WingEtAl2000_LMA_MAT_processed.csv"),
  show_col_types = FALSE
) %>%
  filter(
    is.finite(published_age_model_2_Ma),
    is.finite(MAT_C),
    between(published_age_model_2_Ma, min(age_limits), max(age_limits))
  ) %>%
  transmute(
    Age_Ma = published_age_model_2_Ma,
    age_old_ma = Age_Ma + duration_Myr / 2,
    age_young_ma = Age_Ma - duration_Myr / 2,
    temperature_C = MAT_C,
    temperature_se_C = MAT_error_C,
    record = "Wing LMA MAT"
  )

Fricke_Wing <- read_csv(
  here("data", "processed", "FrickeWing2004_BHB_MAAT_processed.csv"),
  show_col_types = FALSE
) %>%
  filter(
    is.finite(Age_Ma),
    is.finite(temperature_C),
    between(Age_Ma, min(age_limits), max(age_limits))
  ) %>%
  transmute(
    Age_Ma,
    age_old_ma = age_older_ma,
    age_young_ma = age_younger_ma,
    temperature_C,
    temperature_se_C = NA_real_,
    record = "Fricke-Wing MAAT"
  )

seasonal <- read_csv(
  here(
    "data", "processed",
    "BHB_seasonal_temperature_integrated_summary.csv"
  ),
  show_col_types = FALSE
) %>%
  filter(metric %in% c("CMMT", "MAAT", "WMMT")) %>%
  mutate(metric = factor(metric, levels = c("CMMT", "MAAT", "WMMT")))

seasonal_colors <- c(
  "CMMT" = "#2166AC", "MAAT" = "#8C6BB1", "WMMT" = "#D73027"
)

seasonal_rect <- function(metric_name, lower_col, upper_col, alpha_value) {
  geom_rect(
    data = seasonal %>% filter(metric == metric_name),
    aes(
      xmin = .data[[lower_col]], xmax = .data[[upper_col]],
      ymin = min(age_limits), ymax = max(age_limits)
    ),
    fill = seasonal_colors[[metric_name]],
    alpha = alpha_value,
    inherit.aes = FALSE
  )
}

#-- 5.) Panels -------------------------------------------------------------
p_Harper <- ggplot(Harper_SST) +
  add_petm() +
  geom_ribbon(
    aes(
      xmin = Harper2024_SST_lower95_C,
      xmax = Harper2024_SST_upper95_C,
      y = Age_Ma
    ),
    fill = "#FDAE61", alpha = 0.30
  ) +
  geom_path(
    aes(Harper2024_mean_SST_C, Age_Ma),
    color = "#D73027", linewidth = 0.95
  ) +
  shared_age_scale +
  scale_x_continuous(
    limits = c(30, 45), breaks = seq(30, 45, 5),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  labs(
    tag = "A", title = "Pacific marine SST",
    subtitle = "Harper et al.; mean and 95% interval",
    x = expression("SST (" * degree * "C)"), y = "Age (Ma)"
  ) +
  theme_panel

# Horizontal version for the standalone Harper export. The vertical version
# above remains in the three-panel regional synthesis so all three panels
# retain the same age-axis orientation.
p_Harper_horizontal <- ggplot(Harper_SST) +
  annotate(
    "rect",
    xmin = petm_young_ma,
    xmax = petm_old_ma,
    ymin = -Inf,
    ymax = Inf,
    fill = "#D73027",
    alpha = 0.08
  ) +
  geom_ribbon(
    aes(
      x = Age_Ma,
      ymin = Harper2024_SST_lower95_C,
      ymax = Harper2024_SST_upper95_C
    ),
    fill = "#FDAE61",
    alpha = 0.30
  ) +
  geom_path(
    aes(Age_Ma, Harper2024_mean_SST_C),
    color = "#D73027",
    linewidth = 0.95
  ) +
  scale_x_reverse(
    limits = age_limits,
    breaks = seq(59, 52.5, by = -0.5),
    minor_breaks = seq(59, 52.5, by = -0.1),
    expand = expansion(mult = c(0.008, 0.012))
  ) +
  scale_y_continuous(
    limits = c(30, 45),
    breaks = seq(30, 45, by = 5),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  labs(
    tag = "A",
    title = "Pacific marine SST",
    subtitle = "Harper et al.; mean and 95% interval",
    x = "Age (Ma)",
    y = expression("SST (" * degree * "C)")
  ) +
  theme_panel

p_Tornillo <- ggplot() +
  add_petm() +
  geom_ribbon(
    data = Tornillo_temperature_model,
    aes(
      xmin = temperature_lower95_C,
      xmax = temperature_upper95_C,
      y = Age_Ma
    ),
    fill = "#E69F00", alpha = 0.15
  ) +
  geom_ribbon(
    data = Tornillo_temperature_model,
    aes(
      xmin = temperature_lower80_C,
      xmax = temperature_upper80_C,
      y = Age_Ma
    ),
    fill = "#E69F00", alpha = 0.25
  ) +
  geom_path(
    data = Tornillo_temperature_model,
    aes(temperature_model_C, Age_Ma),
    color = "#9A5B00", linewidth = 1
  ) +
  geom_errorbarh(
    data = Tornillo_observations,
    aes(
      xmin = T47_C - T47_se_C,
      xmax = T47_C + T47_se_C,
      y = Age_Ma
    ),
    height = 0, linewidth = 0.35, alpha = 0.50
  ) +
  geom_point(
    data = Tornillo_observations,
    aes(T47_C, Age_Ma),
    shape = 24, fill = "#E69F00", color = "black",
    size = 2.5, stroke = 0.55
  ) +
  shared_age_scale +
  scale_x_continuous(
    limits = c(10, 47), breaks = seq(10, 45, 5),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  labs(
    tag = "B", title = "Tornillo Basin",
    subtitle = "Screened micrite D47 model; 80% and 95% intervals",
    x = "T47 temperature (degrees C)",
    y = NULL
  ) +
  theme_panel +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

BHB_record_colors <- c(
  "This study D47" = "#B2182B",
  "Published D47" = "white",
  "Wing LMA MAT" = "#55A868",
  "Fricke-Wing MAAT" = "#8C6BB1"
)
BHB_record_shapes <- c(
  "This study D47" = 21,
  "Published D47" = 21,
  "Wing LMA MAT" = 22,
  "Fricke-Wing MAAT" = 23
)
BHB_record_levels <- c(
  "This study D47", "Published D47",
  "Wing LMA MAT", "Fricke-Wing MAAT"
)

p_BHB <- ggplot() +
  seasonal_rect("CMMT", "lower95_c", "upper95_c", 0.055) +
  seasonal_rect("CMMT", "lower80_c", "upper80_c", 0.085) +
  seasonal_rect("CMMT", "lower50_c", "upper50_c", 0.12) +
  seasonal_rect("MAAT", "lower95_c", "upper95_c", 0.055) +
  seasonal_rect("MAAT", "lower80_c", "upper80_c", 0.085) +
  seasonal_rect("MAAT", "lower50_c", "upper50_c", 0.12) +
  seasonal_rect("WMMT", "lower95_c", "upper95_c", 0.055) +
  seasonal_rect("WMMT", "lower80_c", "upper80_c", 0.085) +
  seasonal_rect("WMMT", "lower50_c", "upper50_c", 0.12) +
  geom_vline(
    data = seasonal,
    aes(xintercept = mean_c, color = metric),
    linewidth = 0.7, alpha = 0.85, show.legend = FALSE
  ) +
  add_petm() +
  geom_errorbarh(
    data = BHB_D47,
    aes(
      xmin = temperature_C - temperature_se_C,
      xmax = temperature_C + temperature_se_C,
      y = Age_Ma
    ),
    height = 0, linewidth = 0.32, color = "grey35",
    alpha = 0.52, na.rm = TRUE
  ) +
  geom_errorbarh(
    data = Wing_LMA,
    aes(
      xmin = temperature_C - temperature_se_C,
      xmax = temperature_C + temperature_se_C,
      y = Age_Ma
    ),
    height = 0, linewidth = 0.32, color = "#1B7837",
    alpha = 0.55, na.rm = TRUE
  ) +
  geom_errorbar(
    data = bind_rows(Wing_LMA, Fricke_Wing),
    aes(
      x = temperature_C,
      ymin = age_young_ma,
      ymax = age_old_ma
    ),
    width = 0, linewidth = 0.38, color = "grey35", alpha = 0.45
  ) +
  geom_point(
    data = BHB_D47,
    aes(
      temperature_C, Age_Ma,
      shape = record, fill = record
    ),
    color = "black", size = 2.2, stroke = 0.65
  ) +
  geom_point(
    data = bind_rows(Wing_LMA, Fricke_Wing),
    aes(
      temperature_C, Age_Ma,
      shape = record, fill = record
    ),
    color = "black", size = 2.4, stroke = 0.65
  ) +
  geom_text(
    data = seasonal,
    aes(x = mean_c, y = 52.62, label = metric, color = metric),
    size = 2.8, fontface = "bold", show.legend = FALSE
  ) +
  scale_shape_manual(
    values = BHB_record_shapes,
    breaks = BHB_record_levels,
    name = NULL
  ) +
  scale_fill_manual(
    values = BHB_record_colors,
    breaks = BHB_record_levels,
    name = NULL
  ) +
  scale_color_manual(values = seasonal_colors, guide = "none") +
  shared_age_scale +
  scale_x_continuous(
    limits = c(-2, 55), breaks = seq(0, 50, 10),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    tag = "C", title = "Bighorn Basin synthesis",
    subtitle = "Screened D47, other proxies, and seasonal ranges",
    x = expression("Temperature (" * degree * "C)"), y = NULL
  ) +
  guides(
    shape = guide_legend(nrow = 1, byrow = TRUE),
    fill = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme_panel +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(0.45, "cm"),
    legend.spacing.x = unit(0.08, "cm"),
    legend.box.margin = margin(-4, 0, -2, 0)
  )

p_regional_temperature_synthesis <-
  p_Harper + p_Tornillo + p_BHB +
  plot_layout(widths = c(0.92, 1.03, 1.35), guides = "collect") +
  plot_annotation(
    title = "Early Paleogene temperature records from ocean to land",
    caption = str_wrap(
      paste(
        "Pink band: PETM. Tornillo model uses only screened primary micrite;",
        "dark/light orange ribbons are 80%/95% Monte Carlo intervals.",
        "BHB clumped points pass the current d18Ocarb and T47 screen.",
        "BHB blue/purple/red bands are literature-informed CMMT/MAAT/WMMT",
        "50%, 80%, and 95% ranges. These archives record different seasons",
        "and environmental reservoirs."
      ),
      width = 190
    ),
    theme = theme(
      plot.title = element_text(size = 15, face = "bold"),
      plot.caption = element_text(size = 7.3, hjust = 0),
      plot.tag = element_text(size = 12, face = "bold"),
      plot.margin = margin(5, 7, 4, 7)
    )
  ) &
  theme(legend.position = "bottom")

save_figure_variants(
  plot = p_regional_temperature_synthesis,
  presentation_plot = p_regional_temperature_synthesis,
  base_dir = figure_dir,
  stem = "Harper_Tornillo_BHB_temperature_age_synthesis",
  manuscript_width = 10,
  manuscript_height = 5.8,
  presentation_width = 12,
  presentation_height = 6
)

save_figure_variants(
  plot = p_Harper_horizontal,
  presentation_plot = p_Harper_horizontal,
  base_dir = figure_dir,
  stem = "Harper_SST",
  manuscript_width = 10,
  manuscript_height = 5.8,
  presentation_width = 6,
  presentation_height = 3
)

message(
  "Saved Harper-Tornillo-BHB synthesis with ",
  nrow(Tornillo_observations), " screened Tornillo observations, ",
  nrow(Tornillo_age_bins), " Tornillo age bins, and ",
  nrow(BHB_D47), " screened BHB D47 observations."
)
