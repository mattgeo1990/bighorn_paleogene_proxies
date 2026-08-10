# 12_plot_CFB_longterm_temperature_global_context.R
# Purpose: Compare the broad CFB Delta47 temperature trajectory with a
#          continuous marine temperature record from approximately 59-54 Ma.
#
# Trend definition:
#   - The fitted environmental trend uses every pedogenic-micrite T47
#     observation in the "none" screening scenario (no diagenetic exclusions).
#   - Explicit spar and altered-carbonate measurements are plotted as context
#     but are not used to fit the environmental trajectory because they do not
#     represent equivalent soil-carbonate formation-temperature observations.
#   - A low-complexity REML GAM (k = 5) estimates the broad nonlinear mean.
#   - An inverse-variance-weighted linear regression shows the simplest
#     monotonic interpretation.
#   - Neither fit extrapolates beyond the observed CFB age range.
#
# Global context:
#   Barnet et al. (2019) ODP Site 1262 bottom-water temperatures are plotted as
#   a deep-ocean/global-climate context record. They are not a direct estimate
#   of Bighorn Basin air or soil temperature and should not be compared in
#   absolute magnitude with pedogenic-carbonate T47.

#-- 1.) Setup ---------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)
library(mgcv)

source(here("scripts", "helpers", "save_figure_variants.R"))

figure_dir <- here("figures", "age_domain", "regional_comparison")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

age_limits <- c(59, 54)
petm_age_limits <- c(55.93, 55.75)

add_petm_age_band <- function(alpha = 0.10) {
  annotate(
    "rect",
    xmin = petm_age_limits[1], xmax = petm_age_limits[2],
    ymin = -Inf, ymax = Inf,
    fill = "#D73027", alpha = alpha
  )
}

age_x_scale <- scale_x_reverse(
  limits = age_limits,
  breaks = seq(59, 54, by = -0.5),
  minor_breaks = seq(59, 54, by = -0.1),
  expand = expansion(mult = c(0.005, 0.01))
)

theme_comparison <- theme_classic(base_size = 18) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18, color = "black"),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 18),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    plot.margin = margin(4, 6, 4, 6)
  )

#-- 2.) Build the Unscreened Pedogenic T47 Dataset -------------------------
CFB_age_lookup <- read_csv(
  here(
    "data", "processed",
    "CFB_soilcarb_with_temperature_age_calibrated.csv"
  ),
  show_col_types = FALSE
) %>%
  select(MLA_horizon_id, Age_Ma) %>%
  filter(is.finite(Age_Ma)) %>%
  distinct(MLA_horizon_id, .keep_all = TRUE)

CFB_T47_all_pedogenic <- read_csv(
  here("data", "processed", "CFB_temperature_observations.csv"),
  show_col_types = FALSE
) %>%
  left_join(CFB_age_lookup, by = "MLA_horizon_id") %>%
  mutate(
    source = recode(source, IPL = "U-M", Snell = "Caltech"),
    dataset_status = if_else(source == "U-M", "This study", "Published"),
    fit_weight = 1 / pmax(T_se_C, 0.5)^2
  ) %>%
  filter(
    is.finite(Age_Ma),
    is.finite(T_C),
    is.finite(fit_weight),
    between(Age_Ma, age_limits[2], age_limits[1])
  )

if (nrow(CFB_T47_all_pedogenic) < 8) {
  stop("Too few age-resolved pedogenic T47 observations for broad trend fitting.")
}

CFB_T47_model_pedogenic <- CFB_T47_all_pedogenic %>%
  filter(used_in_primary_temperature_model)

if (nrow(CFB_T47_model_pedogenic) < 8) {
  stop("Too few screened pedogenic T47 observations for broad trend fitting.")
}

#-- 3.) Recover and Age Non-Primary Clumped Materials ----------------------
CFB_age_grid <- read_csv(
  here("data", "processed", "BHB_section_age_model_grid.csv"),
  show_col_types = FALSE
) %>%
  filter(section_id == "CFB") %>%
  arrange(strat_height_m)

height_to_age <- function(height_m) {
  approx(
    x = CFB_age_grid$strat_height_m,
    y = CFB_age_grid$Age_Ma,
    xout = height_m,
    rule = 1,
    ties = mean
  )$y
}

CFB_horizon_height_lookup <- read_csv(
  here("data", "processed", "CFB_soilcarb_isotope_summary.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    lookup_key = MLA_horizon_id %>%
      str_remove("^SNELL-") %>%
      str_remove("^PK95-"),
    roster_strat_height_m = strat_height_m
  ) %>%
  filter(is.finite(roster_strat_height_m)) %>%
  group_by(lookup_key) %>%
  summarise(
    roster_strat_height_m = median(roster_strat_height_m),
    .groups = "drop"
  )

CFB_nonprimary_T47 <- read_csv(
  here(
    "data", "processed",
    "CFB_d18O_T47_screening_observations.csv"
  ),
  show_col_types = FALSE
) %>%
  filter(carbonate_type %in% c("Spar", "Altered carbonate")) %>%
  mutate(
    lookup_key = MLA_horizon_id %>%
      str_remove("^SNELL-") %>%
      str_remove("^PK95-")
  ) %>%
  left_join(CFB_horizon_height_lookup, by = "lookup_key") %>%
  mutate(
    strat_height_m = coalesce(strat_height_m, roster_strat_height_m),
    Age_Ma = height_to_age(strat_height_m)
  ) %>%
  filter(
    is.finite(Age_Ma),
    is.finite(T47_C),
    between(Age_Ma, age_limits[2], age_limits[1])
  )

#-- 4.) Fit Broad Smooth and Linear Trends ---------------------------------
CFB_gam <- gam(
  T_C ~ s(Age_Ma, k = 5, bs = "tp"),
  data = CFB_T47_model_pedogenic,
  weights = fit_weight,
  method = "REML"
)

CFB_linear <- lm(
  T_C ~ Age_Ma,
  data = CFB_T47_model_pedogenic,
  weights = fit_weight
)

CFB_prediction_grid <- tibble(
  Age_Ma = seq(
    min(CFB_T47_model_pedogenic$Age_Ma),
    max(CFB_T47_model_pedogenic$Age_Ma),
    length.out = 500
  )
)

gam_prediction <- predict(
  CFB_gam,
  newdata = CFB_prediction_grid,
  se.fit = TRUE,
  type = "response"
)
linear_prediction <- predict(
  CFB_linear,
  newdata = CFB_prediction_grid,
  se.fit = TRUE
)

CFB_prediction_grid <- CFB_prediction_grid %>%
  mutate(
    gam_mean_C = as.numeric(gam_prediction$fit),
    gam_se_C = as.numeric(gam_prediction$se.fit),
    gam_lower95_C = gam_mean_C - 1.96 * gam_se_C,
    gam_upper95_C = gam_mean_C + 1.96 * gam_se_C,
    linear_mean_C = as.numeric(linear_prediction$fit),
    linear_se_C = as.numeric(linear_prediction$se.fit),
    linear_lower95_C = linear_mean_C - 1.96 * linear_se_C,
    linear_upper95_C = linear_mean_C + 1.96 * linear_se_C
  )

CFB_linear_summary <- broom::tidy(CFB_linear, conf.int = TRUE) %>%
  mutate(
    interpretation = case_when(
      term == "Age_Ma" & estimate < 0 ~
        "Negative slope means warming toward younger ages.",
      term == "Age_Ma" ~
        "Positive slope means cooling toward younger ages.",
      TRUE ~ "Regression intercept."
    )
  )

#-- 5.) Prepare the Barnet Marine Temperature Record -----------------------
Barnet_BWT <- read_csv(
  here(
    "data", "processed",
    "BarnetEtAl2019_ODP1262_benthic_isotopes.csv"
  ),
  show_col_types = FALSE
) %>%
  filter(
    is.finite(Age_Ma),
    is.finite(bottom_water_temperature_C),
    between(Age_Ma, age_limits[2], age_limits[1])
  ) %>%
  arrange(Age_Ma)

Barnet_BWT_50kyr <- Barnet_BWT %>%
  mutate(age_bin_ma = floor(Age_Ma / 0.05) * 0.05 + 0.025) %>%
  group_by(age_bin_ma) %>%
  summarise(
    Age_Ma = mean(Age_Ma),
    bottom_water_temperature_C = mean(bottom_water_temperature_C),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(Age_Ma)

#-- 6.) Plot the Aligned Continental and Marine Records --------------------
point_colors <- c("This study" = "#B2182B", "Published" = "grey48")
point_fills <- c("This study" = "#B2182B", "Published" = "white")
point_shapes <- c("This study" = 21, "Published" = 22)

CFB_T47_plot_points <- CFB_T47_model_pedogenic %>%
  transmute(
    Age_Ma,
    temperature_C = T_C,
    temperature_se_C = T_se_C,
    plot_group = if_else(
      dataset_status == "This study", "U-M micrite", "Published micrite"
    )
  ) %>%
  mutate(
    plot_group = factor(
      plot_group,
      levels = c(
        "U-M micrite", "Published micrite",
        "Altered carbonate", "Spar"
      )
    )
  )

plot_group_colors <- c(
  "U-M micrite" = "#B2182B",
  "Published micrite" = "grey48",
  "Altered carbonate" = "#B27A00",
  "Spar" = "#54278F"
)
plot_group_fills <- c(
  "U-M micrite" = "#B2182B",
  "Published micrite" = "white",
  "Altered carbonate" = "#F6C85F",
  "Spar" = "#BCA7D7"
)
plot_group_shapes <- c(
  "U-M micrite" = 21,
  "Published micrite" = 22,
  "Altered carbonate" = 23,
  "Spar" = 24
)

p_CFB_longterm_T47 <- ggplot() +
  add_petm_age_band() +
  geom_ribbon(
    data = CFB_prediction_grid,
    aes(
      x = Age_Ma,
      ymin = gam_lower95_C, ymax = gam_upper95_C
    ),
    fill = "#B2182B", alpha = 0.17
  ) +
  geom_line(
    data = CFB_prediction_grid,
    aes(Age_Ma, gam_mean_C),
    color = "#B2182B", linewidth = 1.15
  ) +
  geom_line(
    data = CFB_prediction_grid,
    aes(Age_Ma, linear_mean_C),
    color = "black", linewidth = 0.75, linetype = 2
  ) +
  geom_errorbar(
    data = CFB_T47_plot_points,
    aes(
      x = Age_Ma,
      ymin = temperature_C - temperature_se_C,
      ymax = temperature_C + temperature_se_C,
      color = plot_group
    ),
    width = 0, linewidth = 0.28, alpha = 0.42,
    na.rm = TRUE
  ) +
  geom_point(
    data = CFB_T47_plot_points,
    aes(
      Age_Ma, temperature_C,
      color = plot_group,
      fill = plot_group,
      shape = plot_group
    ),
    size = 2.15, stroke = 0.75, alpha = 0.9
  ) +
  scale_color_manual(values = plot_group_colors) +
  scale_fill_manual(values = plot_group_fills) +
  scale_shape_manual(values = plot_group_shapes) +
  age_x_scale +
  labs(
    x = NULL,
    y = expression("CFB " * Delta[47] * " (" * degree * "C)"),
    title = "A  CFB pedogenic-carbonate temperature",
    subtitle = paste0(
      "Only observations retained by the temperature-model screen; ",
      "red = REML smooth +/-95%; dashed = weighted linear fit"
    ),
    color = NULL, fill = NULL, shape = NULL
  ) +
  guides(
    color = "none",
    fill = "none",
    shape = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme_comparison +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_Barnet_global_context <- ggplot(Barnet_BWT) +
  add_petm_age_band() +
  geom_line(
    aes(Age_Ma, bottom_water_temperature_C),
    color = "grey68", linewidth = 0.24, alpha = 0.58
  ) +
  geom_line(
    data = Barnet_BWT_50kyr,
    aes(Age_Ma, bottom_water_temperature_C),
    color = "#2166AC", linewidth = 0.9
  ) +
  age_x_scale +
  labs(
    x = "Age (Ma; older to younger)",
    y = expression("ODP 1262 BWT (" * degree * "C)"),
    title = "B  South Atlantic deep-ocean temperature",
    subtitle = paste0(
      "Grey = observations; blue = 50-kyr means; ",
      "deep-ocean context, not a direct CFB analogue"
    )
  ) +
  theme_comparison +
  theme(legend.position = "none")

p_CFB_global_longterm_temperature <-
  p_CFB_longterm_T47 /
  p_Barnet_global_context +
  plot_layout(heights = c(1.10, 0.90)) +
  plot_annotation(
    title = "Continental and marine temperature trajectories, 59-54 Ma",
    subtitle = paste(
      "Shared age direction; independently estimated curves are not tuned"
    ),
    caption = paste(
      "CFB band is uncertainty in the fitted conditional mean, not total",
      "proxy, age-model, seasonal, or preservation uncertainty.",
      "Barnet et al. (2019), ODP Site 1262;",
      "PANGAEA 10.1594/PANGAEA.884585."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 18),
      plot.caption = element_text(size = 18, hjust = 0),
      plot.margin = margin(8, 7, 7, 7)
    )
  )

#-- 7.) Export Reproducible Models, Data, and Figures ----------------------
write_csv(
  CFB_prediction_grid,
  here(
    "data", "processed",
    "CFB_T47_longterm_age_trend_predictions.csv"
  )
)
write_csv(
  CFB_linear_summary,
  here(
    "data", "processed",
    "CFB_T47_longterm_weighted_linear_summary.csv"
  )
)
saveRDS(
  CFB_gam,
  here("data", "processed", "CFB_T47_longterm_gam.rds")
)
saveRDS(
  CFB_linear,
  here("data", "processed", "CFB_T47_longterm_weighted_linear.rds")
)

save_figure_variants(
  p_CFB_global_longterm_temperature,
  figure_dir,
  "CFB_T47_longterm_vs_Barnet2019_global_temperature",
  manuscript_width = 8.5,
  manuscript_height = 8,
  presentation_width = 12,
  presentation_height = 6
)

message(
  "Long-term comparison fit to ",
  nrow(CFB_T47_model_pedogenic),
  " screened pedogenic T47 observations; ",
  nrow(CFB_T47_all_pedogenic) - nrow(CFB_T47_model_pedogenic),
  " primary observations shown but excluded by the talk screen; ",
  nrow(CFB_nonprimary_T47),
  " non-primary observations shown but excluded from the fit."
)
