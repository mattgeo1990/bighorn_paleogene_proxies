# 11_plot_CFB_published_context_slide.R
# Purpose: Build a presentation-sized introduction figure for the published
#          CFB stratigraphic, carbon-isotope, and temperature records.
#
# This figure intentionally excludes all U-M observations. It combines:
#   1. the project CFB chronostratigraphic framework;
#   2. published carbonate d13C measurements (CU, Caltech/Snell, Koch, Bowen);
#   3. published soil-carbonate Delta47 temperatures (CU and Caltech/Snell);
#   4. Wing et al. (2000) leaf-margin MAT estimates; and
#   5. Fricke and Wing (2004) phosphate/LMA-based MAAT estimates.
#
# Wing and Fricke-Wing estimates are published in age or aggregate-age space,
# not at a single CFB meter level. Their central ages and age ranges are
# projected to CFB height by inverting the same section-specific GTS2020 age
# model used elsewhere in this project. Estimates outside the modeled CFB
# height range are not plotted; the code records this explicitly.
#
# The MAAT, CMMT, and WMMT contextual distributions come from the project's
# literature-informed, relevance-weighted Monte Carlo synthesis. They combine
# proxy and model constraints discussed in the associated source table. The
# nested bands are marginal 50%, 80%, and 95% intervals. They are not formal
# Bayesian posteriors, are not fitted to the observations in this figure, and
# are not assumed to vary with stratigraphic height.

#-- 1.) Setup ---------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)

source(here("scripts", "helpers", "CFB_chronostrat_panels.R"))
source(here("scripts", "helpers", "save_figure_variants.R"))

figure_dir <- here("figures", "strat_domain", "CFB")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

strat_limits <- c(500, 2300)
petm_limits_m <- c(1500, 1540)

shared_y <- scale_y_continuous(
  limits = strat_limits,
  breaks = seq(500, 2300, by = 200),
  minor_breaks = seq(500, 2300, by = 100),
  expand = expansion(mult = c(0.005, 0.01))
)

theme_slide_panel <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.height = unit(0.35, "cm"),
    legend.spacing.x = unit(0.10, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )

add_petm_strat <- function(alpha = 0.11) {
  annotate(
    "rect", xmin = -Inf, xmax = Inf,
    ymin = petm_limits_m[1], ymax = petm_limits_m[2],
    fill = "#D73027", alpha = alpha
  )
}

#-- 2.) Load Existing Pipeline Products ------------------------------------
CFB_isotopes <- read_csv(
  here("data", "processed", "CFB_soilcarb_isotope_summary.csv"),
  show_col_types = FALSE
)

CFB_T47 <- read_csv(
  here("data", "processed", "CFB_temperature_observations.csv"),
  show_col_types = FALSE
) %>%
  mutate(source = recode(source, IPL = "U-M", Snell = "Caltech")) %>%
  filter(source != "U-M")

Wing_LMA <- read_csv(
  here("data", "processed", "WingEtAl2000_LMA_MAT_processed.csv"),
  show_col_types = FALSE
)

Fricke_Wing <- read_csv(
  here("data", "processed", "FrickeWing2004_BHB_MAAT_processed.csv"),
  show_col_types = FALSE
)

seasonal_synthesis <- read_csv(
  here(
    "data", "processed",
    "BHB_seasonal_temperature_integrated_summary.csv"
  ),
  show_col_types = FALSE
) %>%
  filter(metric %in% c("CMMT", "MAAT", "WMMT")) %>%
  mutate(
    metric = factor(metric, levels = c("CMMT", "MAAT", "WMMT"))
  )

CFB_age_grid <- read_csv(
  here("data", "processed", "BHB_section_age_model_grid.csv"),
  show_col_types = FALSE
) %>%
  filter(section_id == "CFB") %>%
  filter(is.finite(strat_height_m), is.finite(Age_Ma)) %>%
  arrange(Age_Ma)

if (nrow(CFB_age_grid) < 2) {
  stop("The CFB age-model grid is unavailable or contains fewer than two rows.")
}

#-- 3.) Project Aggregate Ages onto the CFB Meter Scale --------------------
# The CFB model is monotonic, allowing a direct inverse interpolation.
# rule = 1 deliberately prohibits extrapolation beyond the modeled section.
age_to_CFB_height <- function(age_ma) {
  approx(
    x = CFB_age_grid$Age_Ma,
    y = CFB_age_grid$strat_height_m,
    xout = age_ma,
    rule = 1,
    ties = mean
  )$y
}

Wing_CFB <- Wing_LMA %>%
  transmute(
    record = "Wing leaf-margin MAT",
    temperature_C = MAT_C,
    temperature_se_C = MAT_error_C,
    Age_Ma = published_age_model_2_Ma,
    age_older_ma = Age_Ma + duration_Myr / 2,
    age_younger_ma = Age_Ma - duration_Myr / 2
  ) %>%
  mutate(
    strat_height_m = age_to_CFB_height(Age_Ma),
    strat_min_m = age_to_CFB_height(age_older_ma),
    strat_max_m = age_to_CFB_height(age_younger_ma)
  )

Fricke_Wing_CFB <- Fricke_Wing %>%
  transmute(
    record = "Fricke–Wing phosphate MAAT",
    temperature_C,
    temperature_se_C = NA_real_,
    Age_Ma,
    age_older_ma,
    age_younger_ma
  ) %>%
  mutate(
    strat_height_m = age_to_CFB_height(Age_Ma),
    strat_min_m = age_to_CFB_height(age_older_ma),
    strat_max_m = age_to_CFB_height(age_younger_ma)
  )

aggregate_temperature_CFB <- bind_rows(Wing_CFB, Fricke_Wing_CFB) %>%
  filter(
    is.finite(strat_height_m),
    between(strat_height_m, strat_limits[1], strat_limits[2])
  )

projection_audit <- bind_rows(Wing_CFB, Fricke_Wing_CFB) %>%
  mutate(
    plotted = is.finite(strat_height_m) &
      between(strat_height_m, strat_limits[1], strat_limits[2]),
    exclusion_reason = if_else(
      plotted,
      NA_character_,
      "Outside the modeled 500-2300 m CFB age-height range"
    )
  )

write_csv(
  projection_audit,
  here(
    "data", "processed",
    "CFB_published_temperature_strat_projection_audit.csv"
  )
)

#-- 4.) Published Carbon-Isotope Data --------------------------------------
published_d13C <- bind_rows(
  CFB_isotopes %>% transmute(
    strat_height_m, source = "CU",
    value = CU_mean_d13Ccarb_vpdb, se = NA_real_
  ),
  CFB_isotopes %>% transmute(
    strat_height_m, source = "Caltech (Snell)",
    value = Snell_mean_d13Ccarb_vpdb, se = NA_real_
  ),
  CFB_isotopes %>% transmute(
    strat_height_m, source = "Koch",
    value = Koch_mean_d13Ccarb_vpdb, se = Koch_se_d13Ccarb_vpdb
  ),
  CFB_isotopes %>% transmute(
    strat_height_m, source = "Bowen",
    value = Bowen_mean_d13Ccarb_vpdb, se = Bowen_se_d13Ccarb_vpdb
  )
) %>%
  filter(
    is.finite(value),
    between(strat_height_m, strat_limits[1], strat_limits[2])
  ) %>%
  mutate(
    source = factor(
      source, levels = c("Koch", "Bowen", "CU", "Caltech (Snell)")
    )
  )

d13C_height_bins <- published_d13C %>%
  mutate(height_bin_m = round(strat_height_m / 50) * 50) %>%
  group_by(height_bin_m) %>%
  summarise(value = median(value, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(value))

d13C_loess <- loess(
  value ~ height_bin_m,
  data = d13C_height_bins,
  span = 0.22,
  degree = 1,
  control = loess.control(surface = "direct")
)

d13C_trend <- tibble(
  strat_height_m = seq(
    min(d13C_height_bins$height_bin_m),
    max(d13C_height_bins$height_bin_m),
    by = 5
  )
) %>%
  mutate(
    value = predict(
      d13C_loess,
      newdata = data.frame(height_bin_m = strat_height_m)
    )
  ) %>%
  filter(is.finite(value))

p_d13C_published <- ggplot(published_d13C) +
  add_petm_strat() +
  geom_errorbarh(
    aes(
      xmin = value - se, xmax = value + se,
      y = strat_height_m
    ),
    height = 0, color = "grey65", linewidth = 0.25,
    alpha = 0.45, na.rm = TRUE
  ) +
  geom_point(
    aes(value, strat_height_m, shape = source),
    color = "grey42", fill = "white", size = 1.8, stroke = 0.65,
    alpha = 0.86
  ) +
  geom_path(
    data = d13C_trend,
    aes(value, strat_height_m),
    color = "black", linewidth = 0.85
  ) +
  scale_shape_manual(values = c(21, 22, 24, 23), drop = FALSE) +
  scale_x_continuous(
    limits = c(-17, -4),
    breaks = c(-16, -12, -8, -4),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  shared_y +
  labs(
    x = expression(delta^13 * C[carbonate] ~ "(per mil VPDB)"),
    y = "CFB stratigraphic height (m)",
    title = expression("Published " * delta^13 * C[carbonate]),
    shape = NULL
  ) +
  theme_slide_panel +
  theme(
    legend.position = "none"
  )

#-- 5.) Published Temperature Data and Seasonal Context -------------------
temperature_colors <- c(
  "Published T47" = "grey35",
  "Wing LMA" = "#1B7837",
  "Fricke-Wing MAAT" = "#762A83"
)
temperature_shapes <- c(
  "Published T47" = 21,
  "Wing LMA" = 22,
  "Fricke-Wing MAAT" = 23
)
seasonal_colors <- c(
  "CMMT" = "#2166AC",
  "MAAT" = "#8C6BB1",
  "WMMT" = "#D73027"
)

published_T47_plot <- CFB_T47 %>%
  transmute(
    record = "Published T47",
    temperature_C = T_C,
    temperature_se_C = T_se_C,
    strat_height_m,
    strat_min_m = strat_height_m,
    strat_max_m = strat_height_m
  )

temperature_observations <- bind_rows(
  published_T47_plot,
  aggregate_temperature_CFB %>%
    mutate(
      record = recode(
        record,
        "Wing leaf-margin MAT" = "Wing LMA",
        "Fricke–Wing phosphate MAAT" = "Fricke-Wing MAAT"
      )
    )
) %>%
  mutate(
    record = factor(record, levels = names(temperature_colors))
  )

p_temperature_context <- ggplot() +
  add_petm_strat(alpha = 0.08) +
  geom_rect(
    data = seasonal_synthesis,
    aes(
      xmin = lower95_c, xmax = upper95_c,
      ymin = strat_limits[1], ymax = strat_limits[2],
      fill = metric
    ),
    alpha = 0.055, inherit.aes = FALSE
  ) +
  geom_rect(
    data = seasonal_synthesis,
    aes(
      xmin = lower80_c, xmax = upper80_c,
      ymin = strat_limits[1], ymax = strat_limits[2],
      fill = metric
    ),
    alpha = 0.075, inherit.aes = FALSE
  ) +
  geom_rect(
    data = seasonal_synthesis,
    aes(
      xmin = lower50_c, xmax = upper50_c,
      ymin = strat_limits[1], ymax = strat_limits[2],
      fill = metric
    ),
    alpha = 0.11, inherit.aes = FALSE
  ) +
  geom_vline(
    data = seasonal_synthesis,
    aes(xintercept = mean_c, color = metric),
    linewidth = 0.75, alpha = 0.88,
    show.legend = FALSE
  ) +
  geom_errorbarh(
    data = temperature_observations,
    aes(
      xmin = temperature_C - temperature_se_C,
      xmax = temperature_C + temperature_se_C,
      y = strat_height_m,
      color = record
    ),
    height = 0, linewidth = 0.32, alpha = 0.62,
    na.rm = TRUE
  ) +
  geom_errorbar(
    data = temperature_observations %>%
      filter(strat_min_m != strat_max_m),
    aes(
      x = temperature_C, ymin = strat_min_m, ymax = strat_max_m,
      color = record
    ),
    width = 0, linewidth = 0.45, alpha = 0.60,
    na.rm = TRUE
  ) +
  geom_point(
    data = temperature_observations,
    aes(
      temperature_C, strat_height_m,
      shape = record, color = record
    ),
    fill = "white", size = 2.7, stroke = 0.9
  ) +
  scale_fill_manual(values = seasonal_colors, drop = FALSE) +
  scale_color_manual(
    values = c(seasonal_colors, temperature_colors),
    drop = FALSE
  ) +
  scale_shape_manual(values = temperature_shapes, drop = FALSE) +
  scale_x_continuous(
    limits = c(-2, 55),
    breaks = seq(0, 50, by = 10),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  shared_y +
  labs(
    x = expression("Temperature (" * degree * "C)"),
    y = NULL,
    title = "Published temperature records",
    shape = NULL,
    color = NULL,
    fill = "Seasonal synthesis"
  ) +
  guides(
    fill = guide_legend(
      order = 2,
      override.aes = list(alpha = 0.13)
    ),
    shape = guide_legend(order = 1, nrow = 1, byrow = TRUE),
    color = "none"
  ) +
  theme_slide_panel +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top",
    legend.box = "vertical",
    legend.margin = margin(0, 0, 1, 0)
  )

#-- 6.) Assemble the 9 x 6 Inch Introduction Slide Figure -----------------
p_chronostrat <- build_CFB_chronostrat_strat_panel(
  strat_limits,
  breaks = seq(500, 2300, by = 200)
) +
  theme(
    axis.text.x.top = element_text(size = 6.8, face = "bold"),
    plot.margin = margin(4, 1, 4, 2)
  )

provenance_note <- paste(
  "Published proxy data only; U-M observations excluded. Wing and Fricke-Wing",
  "aggregate ages are projected to CFB height using the section-specific",
  "GTS2020 age model (no extrapolation); Fricke-Wing MAAT is phosphate-derived.",
  "\nMAAT/CMMT/WMMT context:",
  "literature-informed relevance-weighted marginal Monte Carlo intervals",
  "(50/80/95%); not fitted posteriors."
)

p_CFB_published_context_slide <-
  p_chronostrat +
  p_d13C_published +
  p_temperature_context +
  plot_layout(widths = c(1.25, 1.25, 2.25)) +
  plot_annotation(
    title = "Clarks Fork Basin stratigraphy and published climate records",
    subtitle = "Late Paleocene–early Eocene; shared CFB composite meter scale",
    caption = provenance_note,
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 11),
      plot.caption = element_text(
        size = 7.1, hjust = 0, lineheight = 1.05, color = "grey25"
      ),
      plot.margin = margin(5, 6, 5, 6)
    )
  )

#-- 7.) Export and Audit ----------------------------------------------------
save_figure_variants(
  plot = p_CFB_published_context_slide,
  presentation_plot = p_CFB_published_context_slide,
  base_dir = figure_dir,
  stem = "CFB_published_strat_d13C_temperature_context",
  manuscript_width = 9,
  manuscript_height = 6,
  presentation_width = 9,
  presentation_height = 6
)

message(
  "Saved 9 x 6 inch CFB introduction figure; ",
  sum(projection_audit$plotted), " of ", nrow(projection_audit),
  " aggregate temperature estimates fall within the modeled CFB section."
)
