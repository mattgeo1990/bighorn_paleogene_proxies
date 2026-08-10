# 03b_plot_BHB_burial_thermal_context.R
# Purpose: Translate the burial-history constraints of Roberts et al. (2008)
#          into a reproducible thermal-depth context for evaluating CFB
#          soil-carbonate alteration.
#
# This is a sensitivity reconstruction, not a sample-specific burial model.
# Roberts et al. calibrated eight 1-D basin models using down-hole temperature,
# vitrinite reflectance, assumed paleosurface temperature, and constant heat
# flow. Their tables 2 and 3 report a thermal gradient and the modeled maximum
# burial depth and temperature of the base of the Fort Union Formation at each
# locality. Maximum burial was placed at approximately 10 Ma.
#
# For comparison with the CFB composite, this script:
#   1. converts the eight reported gradients to degrees C per km;
#   2. infers the maximum-burial surface-temperature intercept required to
#      reproduce each reported base-Fort-Union temperature;
#   3. approximates the CFB K-Pg datum as the base of the Fort Union; and
#   4. projects temperature upward through CFB stratigraphic height using each
#      locality as a separate scenario.
#
# Limitations:
#   - The Roberts wells are not co-located with the CFB outcrop.
#   - Gradients are linearized from basin-scale 1-D models.
#   - The calculation does not model transient conductive equilibration,
#     mineral-specific reordering kinetics, permeability, or fluid flow.
#   - Differences among curves represent spatial/model sensitivity, not a
#     formal probabilistic confidence interval.

#-- 1.) Setup ---------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)
library(cowplot)
source(here("scripts", "helpers", "save_figure_variants.R"))

figure_dir <- here("figures", "diagenesis_screening", "burial_thermal_context")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

feet_to_m <- 0.3048
fahrenheit_to_celsius <- function(x) (x - 32) * 5 / 9

#-- 2.) Load Published Burial-Thermal Constraints --------------------------
thermal_scenarios <- read_csv(
  here("data", "raw", "RobertsEtAl2008_BHB_thermal_scenarios.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    thermal_gradient_C_per_km =
      thermal_gradient_F_per_100ft * (5 / 9) /
      (100 * feet_to_m) * 1000,
    base_Fort_Union_max_depth_km =
      base_Fort_Union_max_depth_ft * feet_to_m / 1000,
    base_Fort_Union_max_temperature_C =
      fahrenheit_to_celsius(base_Fort_Union_max_temperature_F),
    inferred_max_burial_surface_temperature_C =
      base_Fort_Union_max_temperature_C -
      thermal_gradient_C_per_km * base_Fort_Union_max_depth_km,
    display_group = if_else(
      location == "McCulloch Peak",
      "McCulloch Peak model",
      "Other Roberts et al. locations"
    )
  )

if (any(!is.finite(thermal_scenarios$thermal_gradient_C_per_km))) {
  stop("Non-finite thermal gradient after unit conversion.")
}

#-- 3.) Construct Depth and CFB Stratigraphic-Height Profiles --------------
depth_profiles <- thermal_scenarios %>%
  mutate(
    depth_km = map(base_Fort_Union_max_depth_km, ~ seq(0, .x, length.out = 250))
  ) %>%
  unnest(depth_km) %>%
  mutate(
    temperature_C =
      inferred_max_burial_surface_temperature_C +
      thermal_gradient_C_per_km * depth_km
  )

# Common-depth envelope showing the full range spanned by the Roberts et al.
# BHB gradients and calibrated intercepts. This is a scenario range, not a
# confidence interval.
depth_envelope <- crossing(
  thermal_scenarios,
  depth_km = seq(
    0, max(thermal_scenarios$base_Fort_Union_max_depth_km),
    by = 0.01
  )
) %>%
  mutate(
    temperature_C =
      inferred_max_burial_surface_temperature_C +
      thermal_gradient_C_per_km * depth_km
  ) %>%
  group_by(depth_km) %>%
  summarise(
    minimum_C = min(temperature_C),
    maximum_C = max(temperature_C),
    .groups = "drop"
  )

CFB_strat_columns <- read_csv(
  here("data", "raw", "NorBHB_strat_columns.csv"),
  show_col_types = FALSE
)

CFB_formation_intervals <- CFB_strat_columns %>%
  filter(column == "Formation", label %in% c("Fort Union", "Willwood")) %>%
  transmute(
    formation = label,
    base_m,
    top_m
  )

CFB_strat_limits <- range(
  CFB_formation_intervals$base_m,
  CFB_formation_intervals$top_m,
  na.rm = TRUE
)

CFB_height_grid <- tibble(
  strat_height_m = seq(CFB_strat_limits[1], CFB_strat_limits[2], by = 5)
)

CFB_thermal_profiles <- crossing(
  thermal_scenarios,
  CFB_height_grid
) %>%
  mutate(
    max_burial_depth_km = pmax(
      base_Fort_Union_max_depth_km - strat_height_m / 1000,
      0
    ),
    estimated_max_burial_temperature_C =
      inferred_max_burial_surface_temperature_C +
      thermal_gradient_C_per_km * max_burial_depth_km
  )

# The envelope summarizes spatial sensitivity among the eight published
# locality models. It is intentionally described as a scenario range rather
# than a confidence interval.
CFB_thermal_envelope <- CFB_thermal_profiles %>%
  group_by(strat_height_m) %>%
  summarise(
    minimum_C = min(estimated_max_burial_temperature_C),
    lower_quartile_C = quantile(estimated_max_burial_temperature_C, 0.25),
    median_C = median(estimated_max_burial_temperature_C),
    upper_quartile_C = quantile(estimated_max_burial_temperature_C, 0.75),
    maximum_C = max(estimated_max_burial_temperature_C),
    .groups = "drop"
  )

# Export temperatures at reproducible reference horizons, including the
# Fort Union-Willwood contact used by the CFB chronostratigraphic panels.
reference_heights_m <- sort(unique(c(
  CFB_strat_limits,
  CFB_formation_intervals$base_m,
  CFB_formation_intervals$top_m,
  500, 700, 1000, 1480, 1500, 2000, 2300
)))

CFB_thermal_reference_horizons <- crossing(
  thermal_scenarios,
  strat_height_m = reference_heights_m
) %>%
  mutate(
    max_burial_depth_km = pmax(
      base_Fort_Union_max_depth_km - strat_height_m / 1000,
      0
    ),
    estimated_max_burial_temperature_C =
      inferred_max_burial_surface_temperature_C +
      thermal_gradient_C_per_km * max_burial_depth_km
  ) %>%
  select(
    location, burial_depth_class, heat_flow_mW_m2,
    thermal_gradient_C_per_km,
    inferred_max_burial_surface_temperature_C,
    strat_height_m, max_burial_depth_km,
    estimated_max_burial_temperature_C, source
  )

#-- 4.) Plot Published Geotherms --------------------------------------------
scenario_colors <- c(
  "McCulloch Peak model" = "#D95F02",
  "Other Roberts et al. locations" = "grey62"
)

p_basin_geotherms <- ggplot() +
  geom_ribbon(
    data = depth_envelope,
    aes(
      xmin = minimum_C, xmax = maximum_C,
      y = depth_km
    ),
    fill = "grey65", alpha = 0.28
  ) +
  geom_path(
    data = depth_profiles %>% filter(location == "McCulloch Peak"),
    aes(temperature_C, depth_km),
    color = "#D95F02", linewidth = 1.15
  ) +
  scale_y_reverse(expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    x = expression("Estimated maximum-burial temperature (" * degree * "C)"),
    y = "Burial depth (km)",
    title = "A  Roberts et al. (2008) burial-temperature scenarios",
    subtitle = "Orange = McCulloch Peak; grey = full BHB model range"
  ) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none")

#-- 5.) Project Scenarios Through the CFB Composite ------------------------
p_CFB_max_burial_temperature <- ggplot() +
  geom_rect(
    data = CFB_formation_intervals,
    aes(
      xmin = -Inf, xmax = Inf,
      ymin = base_m, ymax = top_m,
      fill = formation
    ),
    alpha = 0.10
  ) +
  geom_ribbon(
    data = CFB_thermal_envelope,
    aes(
      xmin = minimum_C, xmax = maximum_C,
      y = strat_height_m
    ),
    fill = "grey70", alpha = 0.22
  ) +
  geom_path(
    data = CFB_thermal_profiles %>% filter(location == "McCulloch Peak"),
    aes(estimated_max_burial_temperature_C, strat_height_m),
    color = "#D95F02", linewidth = 1.1
  ) +
  scale_fill_manual(
    values = c("Fort Union" = "#4575B4", "Willwood" = "#D73027")
  ) +
  scale_y_continuous(
    limits = CFB_strat_limits,
    breaks = seq(CFB_strat_limits[1], CFB_strat_limits[2], by = 250),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = expression("Estimated maximum-burial temperature (" * degree * "C)"),
    y = "CFB stratigraphic height above K-Pg (m)",
    fill = "Formation",
    title = "B  Thermal projection through the CFB composite",
    subtitle = "Orange = McCulloch Peak; grey = full Roberts et al. BHB range"
  ) +
  theme_classic(base_size = 18) +
  theme(legend.position = "top")

p_BHB_burial_thermal_context <-
  (p_basin_geotherms | p_CFB_max_burial_temperature) +
  plot_annotation(
    title = "Bighorn Basin burial-temperature context for carbonate alteration",
    subtitle = paste(
      "Sensitivity reconstruction from Roberts et al. (2008) tables 2 and 3;",
      "not a sample-specific kinetic or fluid-flow model"
    )
  )

#-- 6.) Compare CFB T47 with the Roberts Thermal Scenarios ----------------
# Delta47 temperatures estimate soil-carbonate formation temperature, whereas
# the Roberts projections estimate maximum burial temperature. Overlaying them
# in temperature-height space is therefore a preservation-screening comparison,
# not evidence that the two quantities should coincide. The separation between
# a sample's T47 and the thermal scenarios indicates the magnitude of later
# heating potentially experienced after carbonate formation.
# Use the complete screening inventory rather than the primary-micrite
# temperature-model input. This retains altered carbonate and spar. Two Snell
# spar rows lack height in the source summary; their positions are recovered
# from matching identifiers in the authoritative CFB horizon roster.
CFB_horizon_lookup <- read_csv(
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

CFB_T47_observations <- read_csv(
  here(
    "data", "processed",
    "CFB_d18O_T47_screening_observations.csv"
  ),
  show_col_types = FALSE
) %>%
  mutate(
    lookup_key = MLA_horizon_id %>%
      str_remove("^SNELL-") %>%
      str_remove("^PK95-")
  ) %>%
  left_join(CFB_horizon_lookup, by = "lookup_key") %>%
  mutate(
    original_strat_height_m = strat_height_m,
    strat_height_m = coalesce(strat_height_m, roster_strat_height_m),
    strat_position_source = case_when(
      is.finite(original_strat_height_m) ~ "screening inventory",
      is.finite(roster_strat_height_m) ~ "matched CFB horizon roster",
      TRUE ~ "unresolved"
    ),
    dataset_status = if_else(
      str_detect(source, regex("U-M|IPL", ignore_case = TRUE)),
      "This study", "Published"
    ),
    carbonate_type = factor(
      carbonate_type,
      levels = c("Pedogenic micrite", "Altered carbonate", "Spar")
    ),
    p_altered_preservation = if_else(
      is.finite(p_altered_preservation),
      p_altered_preservation,
      NA_real_
    )
  ) %>%
  filter(
    is.finite(T47_C),
    is.finite(strat_height_m),
    between(
      strat_height_m,
      CFB_strat_limits[1],
      CFB_strat_limits[2]
    ) 
  ) 
  

# Fit the burial-panel trend to every primary pedogenic-micrite observation,
# including measurements omitted by the deterministic production-model screen.
# Missing SE values receive the median reported primary-micrite SE so that no
# primary observation is silently dropped from the weighted fit.
CFB_T47_primary_for_linear_fit <- CFB_T47_observations %>%
  filter(
    carbonate_type == "Pedogenic micrite",
    is.finite(T47_C),
    is.finite(strat_height_m)
  ) %>%
  mutate(
    fit_se_C = if_else(
      is.finite(T47_se_C),
      pmax(T47_se_C, 0.5),
      median(T47_se_C[is.finite(T47_se_C)], na.rm = TRUE)
    ),
    fit_weight = 1 / fit_se_C^2
  )

CFB_T47_strat_linear_fit <- lm(
  T47_C ~ strat_height_m,
  data = CFB_T47_primary_for_linear_fit,
  weights = fit_weight
)

CFB_T47_strat_linear_grid <- tibble(
  strat_height_m = seq(
    min(CFB_T47_primary_for_linear_fit$strat_height_m),
    max(CFB_T47_primary_for_linear_fit$strat_height_m),
    length.out = 300
  )
)

CFB_T47_strat_linear_prediction <- predict(
  CFB_T47_strat_linear_fit,
  newdata = CFB_T47_strat_linear_grid,
  se.fit = TRUE
)

CFB_T47_strat_linear_grid <- CFB_T47_strat_linear_grid %>%
  mutate(
    T_linear_C = as.numeric(CFB_T47_strat_linear_prediction$fit),
    T_linear_se_C = as.numeric(CFB_T47_strat_linear_prediction$se.fit),
    T_linear_lower95_C = T_linear_C - 1.96 * T_linear_se_C,
    T_linear_upper95_C = T_linear_C + 1.96 * T_linear_se_C
  )

CFB_T47_plot_audit <- CFB_T47_observations %>%
  count(
    carbonate_type, dataset_status, strat_position_source,
    name = "n_plotted"
  )

p_CFB_T47_burial_comparison <- ggplot() +
  geom_rect(
    data = CFB_formation_intervals %>%
      filter(formation == "Fort Union"),
    aes(
      xmin = -Inf, xmax = Inf,
      ymin = base_m, ymax = top_m
    ),
    fill = "#4575B4", alpha = 0.045
  ) +
  geom_rect(
    data = CFB_formation_intervals %>%
      filter(formation == "Willwood"),
    aes(
      xmin = -Inf, xmax = Inf,
      ymin = base_m, ymax = top_m
    ),
    fill = "#D73027", alpha = 0.035
  ) +
  geom_hline(
    yintercept = 1480,
    color = "grey45", linewidth = 0.45, linetype = "dashed"
  ) +
  geom_ribbon(
    data = CFB_thermal_envelope,
    aes(
      xmin = minimum_C, xmax = maximum_C,
      y = strat_height_m
    ),
    fill = "grey65", alpha = 0.18
  ) +
  geom_path(
    data = CFB_thermal_profiles %>% filter(location == "McCulloch Peak"),
    aes(estimated_max_burial_temperature_C, strat_height_m),
    color = "#D95F02", linewidth = 1.05
  ) +
  geom_ribbon(
    data = CFB_T47_strat_linear_grid,
    aes(
      xmin = T_linear_lower95_C, xmax = T_linear_upper95_C,
      y = strat_height_m
    ),
    fill = "#2166AC", alpha = 0.14
  ) +
  geom_path(
    data = CFB_T47_strat_linear_grid,
    aes(T_linear_C, strat_height_m),
    color = "#2166AC", linewidth = 1.0
  ) +
  geom_errorbarh(
    data = CFB_T47_observations,
    aes(
      xmin = T47_C - T47_se_C, xmax = T47_C + T47_se_C,
      y = strat_height_m, color = dataset_status
    ),
    height = 0, linewidth = 0.30, alpha = 0.52,
    na.rm = TRUE
  ) +
  geom_point(
    data = CFB_T47_observations,
    aes(
      T47_C, strat_height_m,
      color = dataset_status, fill = p_altered_preservation,
      shape = carbonate_type
    ),
    size = 2.65, stroke = 0.85
  ) +
  geom_text(
    data = CFB_formation_intervals,
    aes(
      x = 141,
      y = (base_m + top_m) / 2,
      label = paste0(formation, " Fm")
    ),
    angle = 90, fontface = "bold", size = 18 / ggplot2::.pt,
    color = "grey30"
  ) +
  scale_fill_gradientn(
    colors = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"),
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::label_percent(accuracy = 1),
    name = "d18O trajectory\nP(altered)"
  ) +
  scale_color_manual(
    values = c(
      "This study" = "black",
      "Published" = "grey42"
    ),
    name = "Data provenance"
  ) +
  scale_shape_manual(
    values = c(
      "Pedogenic micrite" = 21,
      "Altered carbonate" = 22,
      "Spar" = 24
    )
  ) +
  scale_y_continuous(
    limits = CFB_strat_limits,
    breaks = seq(CFB_strat_limits[1], CFB_strat_limits[2], by = 250),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_x_continuous(
    limits = c(0, 145),
    breaks = c(0, 30, 60, 90, 120, 140),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = expression("Temperature (" * degree * "C)"),
    y = "CFB stratigraphic height above K-Pg (m)",
    shape = "CFB clumped material",
    title = "B  CFB T47 and maximum-burial temperature",
    subtitle = paste0(
      "Points = carbonate formation T; orange = McCulloch Peak;\n",
      "grey = Roberts BHB range; blue = all-primary fit; fill = alteration index"
    )
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(
      order = 2, nrow = 1, byrow = TRUE,
      override.aes = list(fill = "white", color = "black")
    ),
    fill = guide_colorbar(order = 3, barwidth = grid::unit(4.2, "cm"))
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 18, lineheight = 1.05),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.justification = "center",
    legend.margin = margin(1, 1, 1, 1),
    plot.margin = margin(5, 5, 5, 5)
  )

p_CFB_T47_Roberts_geothermal_comparison <- p_CFB_T47_burial_comparison

#-- 7.) Export Data and Figures --------------------------------------------
write_csv(
  thermal_scenarios,
  here("data", "processed", "RobertsEtAl2008_BHB_thermal_scenarios_converted.csv")
)
write_csv(
  CFB_thermal_profiles,
  here("data", "processed", "CFB_maximum_burial_temperature_scenarios.csv")
)
write_csv(
  CFB_thermal_envelope,
  here("data", "processed", "CFB_maximum_burial_temperature_envelope.csv")
)
write_csv(
  CFB_thermal_reference_horizons,
  here(
    "data", "processed",
    "CFB_maximum_burial_temperature_reference_horizons.csv"
  )
)
write_csv(
  CFB_T47_plot_audit,
  here(
    "data", "processed",
    "CFB_T47_Roberts_geothermal_plot_audit.csv"
  )
)
write_csv(
  CFB_T47_strat_linear_grid,
  here(
    "data", "processed",
    "CFB_T47_strat_height_linear_predictions.csv"
  )
)
write_csv(
  broom::tidy(CFB_T47_strat_linear_fit, conf.int = TRUE),
  here(
    "data", "processed",
    "CFB_T47_strat_height_linear_summary.csv"
  )
)

save_figure_variants(
  p_BHB_burial_thermal_context, figure_dir,
  "BHB_Roberts2008_burial_temperature_context", 13, 7.5,
  presentation_width = 12
)

save_figure_variants(
  p_CFB_T47_Roberts_geothermal_comparison, figure_dir,
  "CFB_T47_vs_Roberts2008_geothermal_gradients",
  manuscript_width = 6,
  manuscript_height = 6,
  presentation_width = 6,
  presentation_height = 7
)

print(thermal_scenarios)
print(
  CFB_thermal_reference_horizons %>%
    filter(strat_height_m %in% c(0, 1480, 2300)) %>%
    select(location, strat_height_m, estimated_max_burial_temperature_C)
)
