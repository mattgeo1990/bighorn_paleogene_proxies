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

p_basin_geotherms <- ggplot(
  depth_profiles,
  aes(
    temperature_C, depth_km,
    group = location, color = display_group
  )
) +
  geom_path(linewidth = 0.85, alpha = 0.82) +
  geom_point(
    data = thermal_scenarios,
    aes(
      base_Fort_Union_max_temperature_C,
      base_Fort_Union_max_depth_km,
      color = display_group
    ),
    inherit.aes = FALSE, size = 2.2
  ) +
  scale_y_reverse(expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = scenario_colors) +
  labs(
    x = expression("Estimated maximum-burial temperature (" * degree * "C)"),
    y = "Burial depth (km)",
    color = NULL,
    title = "A  Roberts et al. (2008) burial-temperature scenarios",
    subtitle = paste(
      "Eight calibrated basin models;",
      "points mark the modeled base of the Fort Union"
    )
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top")

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
  geom_ribbon(
    data = CFB_thermal_envelope,
    aes(
      xmin = lower_quartile_C, xmax = upper_quartile_C,
      y = strat_height_m
    ),
    fill = "grey45", alpha = 0.25
  ) +
  geom_path(
    data = CFB_thermal_profiles %>% filter(location == "McCulloch Peak"),
    aes(estimated_max_burial_temperature_C, strat_height_m),
    color = "#D95F02", linewidth = 1.1
  ) +
  geom_path(
    data = CFB_thermal_envelope,
    aes(median_C, strat_height_m),
    color = "black", linewidth = 0.9, linetype = 2
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
    subtitle = "Orange = McCulloch Peak; dashed = median; bands = scenario ranges"
  ) +
  theme_classic(base_size = 11) +
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

#-- 6.) Export Data and Figures --------------------------------------------
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

save_figure_variants(
  p_BHB_burial_thermal_context, figure_dir,
  "BHB_Roberts2008_burial_temperature_context", 13, 7.5,
  presentation_width = 12
)

print(thermal_scenarios)
print(
  CFB_thermal_reference_horizons %>%
    filter(strat_height_m %in% c(0, 1480, 2300)) %>%
    select(location, strat_height_m, estimated_max_burial_temperature_C)
)
