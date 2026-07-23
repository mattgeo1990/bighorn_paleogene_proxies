# 10_plot_BHB_age_domain.R
# Purpose: Generate publication-ready age-domain figures for the CFB record
#          and explicitly separated regional/global reference datasets.

#-- 1.) Setup ---------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)

cfb_age_figure_dir <- here("figures", "age_domain", "CFB")
regional_age_figure_dir <- here("figures", "age_domain", "regional_comparison")
dir.create(cfb_age_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(regional_age_figure_dir, recursive = TRUE, showWarnings = FALSE)

petm_age_old_ma <- 55.93
petm_age_young_ma <- 55.75
cfb_age_limits <- c(58.9, 54.0)
regional_age_limits <- c(59.0, 52.5)

source_colors <- c(
  "IPL" = "#000000", "CU" = "#0072B2", "Snell" = "#D55E00",
  "Koch" = "#009E73", "Bowen" = "#CC79A7"
)
source_shapes <- c(
  "IPL" = 21, "CU" = 22, "Snell" = 24, "Koch" = 23, "Bowen" = 25
)

add_petm_age <- function(alpha = 0.14) {
  annotate(
    "rect", xmin = -Inf, xmax = Inf,
    ymin = petm_age_young_ma, ymax = petm_age_old_ma,
    fill = "#E41A1C", alpha = alpha
  )
}

theme_age <- theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 11),
    plot.margin = margin(5, 5, 5, 5)
  )

age_scale <- function(limits) {
  scale_y_reverse(
    limits = limits,
    breaks = seq(ceiling(limits[2] * 2) / 2,
                 floor(limits[1] * 2) / 2, by = 0.5),
    minor_breaks = seq(ceiling(limits[2] * 10) / 10,
                       floor(limits[1] * 10) / 10, by = 0.1),
    expand = expansion(mult = c(0.01, 0.02))
  )
}

save_age_plot <- function(plot, directory, stem, width, height) {
  ggsave(
    file.path(directory, paste0(stem, ".png")),
    plot, width = width, height = height, dpi = 500
  )
  ggsave(
    file.path(directory, paste0(stem, ".pdf")),
    plot, width = width, height = height, device = grDevices::pdf,
    useDingbats = FALSE
  )
}

#-- 2.) Load Age-Calibrated Primary Products -------------------------------
CFB_isotopes_age <- read_csv(
  here("data", "processed", "CFB_soilcarb_isotope_summary_age_calibrated.csv"),
  show_col_types = FALSE
)
CFB_temperature_age <- read_csv(
  here("data", "processed", "CFB_soilcarb_with_temperature_age_calibrated.csv"),
  show_col_types = FALSE
)
CFB_soilwater_age <- read_csv(
  here(
    "data", "processed",
    "CFB_soilwater_reconstruction_summary_age_calibrated.csv"
  ),
  show_col_types = FALSE
)
CFB_temperature_obs <- read_csv(
  here("data", "processed", "CFB_temperature_observations.csv"),
  show_col_types = FALSE
)

age_lookup <- CFB_isotopes_age %>%
  distinct(section_id, MLA_horizon_id, strat_height_m, Age_Ma)

CFB_temperature_obs_age <- CFB_temperature_obs %>%
  left_join(
    age_lookup,
    by = c("section_id", "MLA_horizon_id", "strat_height_m")
  ) %>%
  filter(!is.na(Age_Ma))

#-- 3.) Reshape CFB Carbon-Isotope Observations ----------------------------
d13C_age <- bind_rows(
  CFB_isotopes_age %>% transmute(
    MLA_horizon_id, Age_Ma, source = "IPL",
    value = IPL_NuDog_d13Ccarb_VPDB, se = NA_real_
  ),
  CFB_isotopes_age %>% transmute(
    MLA_horizon_id, Age_Ma, source = "CU",
    value = CU_mean_d13Ccarb_vpdb, se = NA_real_
  ),
  CFB_isotopes_age %>% transmute(
    MLA_horizon_id, Age_Ma, source = "Snell",
    value = Snell_mean_d13Ccarb_vpdb, se = NA_real_
  ),
  CFB_isotopes_age %>% transmute(
    MLA_horizon_id, Age_Ma, source = "Koch",
    value = Koch_mean_d13Ccarb_vpdb, se = Koch_se_d13Ccarb_vpdb
  ),
  CFB_isotopes_age %>% transmute(
    MLA_horizon_id, Age_Ma, source = "Bowen",
    value = Bowen_mean_d13Ccarb_vpdb, se = Bowen_se_d13Ccarb_vpdb
  )
) %>%
  filter(!is.na(value), !is.na(Age_Ma)) %>%
  mutate(source = factor(source, levels = names(source_colors)))

#-- 4.) Construct CFB Age-Domain Panels ------------------------------------
p_temperature_age <- ggplot() +
  add_petm_age() +
  geom_ribbon(
    data = CFB_temperature_age %>% filter(!is.na(T_model_C)),
    aes(
      xmin = T_model_lower95_C, xmax = T_model_upper95_C,
      y = Age_Ma
    ),
    fill = "#B2182B", alpha = 0.18
  ) +
  geom_path(
    data = CFB_temperature_age %>% filter(!is.na(T_model_C)) %>%
      arrange(desc(Age_Ma)),
    aes(T_model_C, Age_Ma), color = "#B2182B", linewidth = 1
  ) +
  geom_errorbarh(
    data = CFB_temperature_obs_age,
    aes(
      xmin = T_C - T_se_C, xmax = T_C + T_se_C,
      y = Age_Ma, color = source
    ),
    height = 0, linewidth = 0.35, alpha = 0.55
  ) +
  geom_point(
    data = CFB_temperature_obs_age,
    aes(T_C, Age_Ma, color = source, shape = source),
    size = 2, fill = "white", stroke = 0.7
  ) +
  scale_color_manual(values = source_colors, drop = FALSE) +
  scale_shape_manual(values = source_shapes, drop = FALSE) +
  age_scale(cfb_age_limits) +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = "Age (Ma)", color = "Source", shape = "Source",
    title = "A  CFB temperature"
  ) +
  theme_age

p_d13C_age <- ggplot(d13C_age) +
  add_petm_age() +
  geom_errorbarh(
    aes(xmin = value - se, xmax = value + se, y = Age_Ma, color = source),
    height = 0, linewidth = 0.3, alpha = 0.35, na.rm = TRUE
  ) +
  geom_point(
    aes(value, Age_Ma, color = source, shape = source),
    size = 1.7, alpha = 0.78, fill = "white", stroke = 0.55
  ) +
  scale_color_manual(values = source_colors, drop = FALSE) +
  scale_shape_manual(values = source_shapes, drop = FALSE) +
  age_scale(cfb_age_limits) +
  labs(
    x = expression(delta^13 * C[carbonate] ~ "(per mil VPDB)"),
    y = "Age (Ma)", color = "Source", shape = "Source",
    title = expression("B  CFB " * delta^13 * C[carbonate])
  ) +
  theme_age

p_d18Owater_age <- CFB_soilwater_age %>%
  filter(!is.na(d18Ow_mean_vsmow), !is.na(Age_Ma)) %>%
  ggplot() +
  add_petm_age() +
  geom_errorbarh(
    aes(
      xmin = d18Ow_lower95_vsmow, xmax = d18Ow_upper95_vsmow,
      y = Age_Ma
    ),
    height = 0, linewidth = 0.4, alpha = 0.45, color = "#2166AC"
  ) +
  geom_point(
    aes(d18Ow_mean_vsmow, Age_Ma),
    shape = 21, fill = "white", color = "#2166AC", size = 2
  ) +
  age_scale(cfb_age_limits) +
  labs(
    x = expression(delta^18 * O[water] ~ "(per mil VSMOW)"),
    y = "Age (Ma)",
    title = expression("C  CFB " * delta^18 * O[water])
  ) +
  theme_age + theme(legend.position = "none")

p_D17Owater_age <- CFB_soilwater_age %>%
  filter(!is.na(D17Orsw_mean_permeg), !is.na(Age_Ma)) %>%
  ggplot() +
  add_petm_age() +
  geom_errorbarh(
    aes(
      xmin = D17Orsw_lower95_permeg, xmax = D17Orsw_upper95_permeg,
      y = Age_Ma
    ),
    height = 0, linewidth = 0.4, alpha = 0.45, color = "#762A83"
  ) +
  geom_point(
    aes(D17Orsw_mean_permeg, Age_Ma),
    shape = 21, fill = "white", color = "#762A83", size = 2
  ) +
  age_scale(cfb_age_limits) +
  labs(
    x = expression(Delta*"'"^17 * O[water] ~ "(per meg)"),
    y = "Age (Ma)",
    title = expression("D  CFB " * Delta*"'"^17 * O[water])
  ) +
  theme_age + theme(legend.position = "none")

# Chronology-support panel: this is a distance-to-control diagnostic, not a
# formal age uncertainty. Larger values indicate greater distance from the
# nearest absolute-age tie point along the CFB section.
p_age_control <- CFB_isotopes_age %>%
  distinct(Age_Ma, distance_to_nearest_prior_m, age_control_distance_index) %>%
  filter(!is.na(Age_Ma), !is.na(distance_to_nearest_prior_m)) %>%
  ggplot(aes(distance_to_nearest_prior_m, Age_Ma,
             color = age_control_distance_index)) +
  add_petm_age() +
  geom_path(linewidth = 0.8) +
  geom_point(size = 1.2) +
  scale_color_viridis_c(option = "C", direction = -1) +
  age_scale(cfb_age_limits) +
  labs(
    x = "Distance to nearest age prior (m)", y = "Age (Ma)",
    color = "Relative control distance",
    title = "E  Chronology support"
  ) +
  theme_age

remove_repeated_age <- function(plot) {
  plot + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
}

p_CFB_age_multiproxy <-
  (p_temperature_age + theme(legend.position = "none")) +
  remove_repeated_age(p_d13C_age) +
  remove_repeated_age(p_d18Owater_age) +
  remove_repeated_age(p_D17Owater_age) +
  remove_repeated_age(p_age_control) +
  plot_layout(widths = c(1.15, 1, 1, 1, 0.95), guides = "collect") &
  theme(legend.position = "top")

#-- 5.) Regional and Global Temperature Comparisons ------------------------
BHB_regional_reference <- read_csv(
  here(
    "data", "processed",
    "BHB_regional_soilcarb_reference_summary_age_calibrated.csv"
  ),
  show_col_types = FALSE
)
Kelson_Tornillo <- read_csv(
  here("data", "processed", "Kelson_Tornillo_D47_processed.csv"),
  show_col_types = FALSE
)
Harper2024 <- read_csv(
  here("data", "processed", "Harper2024_CO2_SST_processed.csv"),
  show_col_types = FALSE
)
Wing2000_MAT <- read_csv(
  here("data", "processed", "WingEtAl2000_LMA_MAT_processed.csv"),
  show_col_types = FALSE
)
FrickeWing2004_BHB_MAAT <- read_csv(
  here("data", "processed", "FrickeWing2004_BHB_MAAT_processed.csv"),
  show_col_types = FALSE
)
BHB_insolation_47N <- read_csv(
  here("data", "processed", "BHB_ZB20a_summer_insolation_47N.csv"),
  show_col_types = FALSE
)
BHB_LMA_MAT_model <- read_csv(
  here("data", "processed", "BHB_LMA_MAT_model.csv"),
  show_col_types = FALSE
)
BHB_LMA_MAT_observations <- read_csv(
  here("data", "processed", "BHB_LMA_MAT_observations.csv"),
  show_col_types = FALSE
)
BHB_D47_temperature_model <- read_csv(
  here("data", "processed", "BHB_D47_temperature_model.csv"),
  show_col_types = FALSE
)
BHB_D47_temperature_observations <- read_csv(
  here("data", "processed", "BHB_D47_temperature_observations.csv"),
  show_col_types = FALSE
)

# Standardize every BHB temperature estimate that can be placed on a numerical
# age axis. These records measure different aspects of climate: carbonate D47
# primarily reflects warm-season soil temperature, whereas leaf-margin and the
# Fricke-Wing estimates represent mean annual air temperature. They are shown
# together for temporal comparison, not treated as interchangeable estimates.
BHB_temperature_proxy_data <- bind_rows(
  CFB_temperature_obs_age %>%
    transmute(
      dataset = paste0(source, " CFB D47"),
      temperature_record = "CFB soil carbonate D47",
      proxy_type = "Warm-season soil temperature",
      section_id = "CFB",
      Age_Ma,
      age_younger_ma = Age_Ma,
      age_older_ma = Age_Ma,
      temperature_C = T_C,
      temperature_se_C = T_se_C
    ),
  BHB_regional_reference %>%
    filter(!is.na(T47_C), !is.na(Age_Ma)) %>%
    transmute(
      dataset,
      temperature_record = "MCP soil carbonate D47",
      proxy_type = "Warm-season soil temperature",
      section_id,
      Age_Ma,
      age_younger_ma = Age_Ma,
      age_older_ma = Age_Ma,
      temperature_C = T47_C,
      temperature_se_C = T47_se_C
    ),
  Wing2000_MAT %>%
    filter(!is.na(MAT_C), !is.na(published_age_model_2_Ma)) %>%
    transmute(
      dataset,
      temperature_record = "Wing leaf-margin MAT",
      proxy_type = "Mean annual air temperature",
      section_id = "AGGREGATE_BHB",
      Age_Ma = published_age_model_2_Ma,
      age_younger_ma = Age_Ma - duration_Myr / 2,
      age_older_ma = Age_Ma + duration_Myr / 2,
      temperature_C = MAT_C,
      temperature_se_C = MAT_error_C
    ),
  FrickeWing2004_BHB_MAAT %>%
    filter(!is.na(temperature_C), !is.na(Age_Ma)) %>%
    transmute(
      dataset,
      temperature_record = "Fricke and Wing MAAT",
      proxy_type = "Mean annual air temperature",
      section_id = "AGGREGATE_BHB",
      Age_Ma,
      age_younger_ma,
      age_older_ma,
      temperature_C,
      temperature_se_C = NA_real_
    )
) %>%
  filter(
    Age_Ma >= regional_age_limits[2],
    Age_Ma <= regional_age_limits[1]
  ) %>%
  mutate(
    temperature_record = factor(
      temperature_record,
      levels = c(
        "CFB soil carbonate D47",
        "MCP soil carbonate D47",
        "Wing leaf-margin MAT",
        "Fricke and Wing MAAT"
      )
    )
  )

temperature_record_colors <- c(
  "CFB soil carbonate D47" = "#B2182B",
  "MCP soil carbonate D47" = "#E66101",
  "Wing leaf-margin MAT" = "#1B7837",
  "Fricke and Wing MAAT" = "#762A83"
)
temperature_record_shapes <- c(
  "CFB soil carbonate D47" = 21,
  "MCP soil carbonate D47" = 24,
  "Wing leaf-margin MAT" = 22,
  "Fricke and Wing MAAT" = 23
)

p_regional_temperature <- ggplot() +
  add_petm_age() +
  geom_ribbon(
    data = CFB_temperature_age %>% filter(!is.na(T_model_C)),
    aes(
      xmin = T_model_lower95_C, xmax = T_model_upper95_C,
      y = Age_Ma
    ),
    fill = "grey70", alpha = 0.25
  ) +
  geom_path(
    data = CFB_temperature_age %>% filter(!is.na(T_model_C)) %>%
      arrange(desc(Age_Ma)),
    aes(T_model_C, Age_Ma), color = "black", linewidth = 1
  ) +
  geom_errorbarh(
    data = BHB_regional_reference %>% filter(!is.na(T47_C)),
    aes(xmin = T47_C - T47_se_C, xmax = T47_C + T47_se_C, y = Age_Ma),
    height = 0, color = "#D55E00", alpha = 0.55
  ) +
  geom_point(
    data = BHB_regional_reference %>% filter(!is.na(T47_C)),
    aes(T47_C, Age_Ma, shape = section_id),
    color = "#D55E00", fill = "white", size = 2.6
  ) +
  geom_errorbarh(
    data = Kelson_Tornillo %>%
      filter(Age_Ma >= regional_age_limits[2],
             Age_Ma <= regional_age_limits[1], !is.na(T47_C)),
    aes(xmin = T47_C - T47_se_C, xmax = T47_C + T47_se_C, y = Age_Ma),
    height = 0, color = "#0072B2", alpha = 0.45
  ) +
  geom_point(
    data = Kelson_Tornillo %>%
      filter(Age_Ma >= regional_age_limits[2],
             Age_Ma <= regional_age_limits[1], !is.na(T47_C)),
    aes(T47_C, Age_Ma),
    shape = 22, fill = "white", color = "#0072B2", size = 2
  ) +
  age_scale(regional_age_limits) +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = "Age (Ma)", shape = "BHB section",
    title = "Continental temperature",
    subtitle = "Black: CFB; orange: other BHB sections; blue: Tornillo Basin"
  ) +
  theme_age

p_Harper_SST <- Harper2024 %>%
  filter(Age_Ma >= regional_age_limits[2], Age_Ma <= regional_age_limits[1]) %>%
  ggplot() +
  add_petm_age() +
  geom_ribbon(
    aes(
      xmin = Harper2024_SST_lower95_C,
      xmax = Harper2024_SST_upper95_C,
      y = Age_Ma
    ),
    fill = "#FDAE61", alpha = 0.28
  ) +
  geom_path(
    aes(Harper2024_mean_SST_C, Age_Ma),
    color = "#D73027", linewidth = 0.9
  ) +
  age_scale(regional_age_limits) +
  labs(x = "Marine SST (degrees C)", y = "Age (Ma)",
       title = "Marine temperature") +
  theme_age + theme(legend.position = "none")

p_Harper_CO2 <- Harper2024 %>%
  filter(Age_Ma >= regional_age_limits[2], Age_Ma <= regional_age_limits[1]) %>%
  ggplot() +
  add_petm_age() +
  geom_ribbon(
    aes(
      xmin = Harper2024_CO2_lower95_ppm,
      xmax = Harper2024_CO2_upper95_ppm,
      y = Age_Ma
    ),
    fill = "#91CF60", alpha = 0.28
  ) +
  geom_path(
    aes(Harper2024_mean_CO2_ppm, Age_Ma),
    color = "#1A9850", linewidth = 0.9
  ) +
  age_scale(regional_age_limits) +
  labs(x = expression(CO[2] ~ "(ppm)"), y = "Age (Ma)",
       title = expression("Atmospheric " * CO[2])) +
  theme_age + theme(legend.position = "none")

#-- 6.) BHB Temperature and Astronomical-Forcing Comparison ----------------

# MAT and pedogenic-carbonate formation temperature are plotted in separate
# panels because they target different parts of the climate system. Both model
# curves stop at their observation limits; neither is extrapolated.
p_BHB_LMA_MAT_age <- ggplot() +
  add_petm_age() +
  geom_ribbon(
    data = BHB_LMA_MAT_model %>% filter(retained_for_plot),
    aes(
      xmin = MAT_lower95_C, xmax = MAT_upper95_C,
      y = Age_Ma
    ),
    fill = "#1B7837", alpha = 0.16
  ) +
  geom_ribbon(
    data = BHB_LMA_MAT_model %>% filter(retained_for_plot),
    aes(
      xmin = MAT_lower80_C, xmax = MAT_upper80_C,
      y = Age_Ma
    ),
    fill = "#1B7837", alpha = 0.24
  ) +
  geom_path(
    data = BHB_LMA_MAT_model %>%
      filter(retained_for_plot) %>%
      arrange(desc(Age_Ma)),
    aes(MAT_median_C, Age_Ma), color = "#1B7837", linewidth = 0.95
  ) +
  geom_errorbarh(
    data = BHB_LMA_MAT_observations,
    aes(
      xmin = MAT_C - MAT_se_C, xmax = MAT_C + MAT_se_C,
      y = Age_Ma
    ),
    height = 0, color = "#1B7837", alpha = 0.65
  ) +
  geom_errorbar(
    data = BHB_LMA_MAT_observations,
    aes(x = MAT_C, ymin = age_younger_ma, ymax = age_older_ma),
    width = 0, color = "#1B7837", alpha = 0.65
  ) +
  geom_point(
    data = BHB_LMA_MAT_observations,
    aes(MAT_C, Age_Ma),
    shape = 21, fill = "white", color = "#1B7837", size = 2.4
  ) +
  age_scale(regional_age_limits) +
  labs(
    x = expression("LMA MAT (" * degree * "C)"), y = "Age (Ma)",
    title = "A  Mean annual air temperature",
    subtitle = "Wing LMA; 80% and 95% Monte Carlo intervals"
  ) +
  theme_age + theme(legend.position = "none")

p_BHB_D47_temperature_age <- ggplot() +
  add_petm_age() +
  geom_ribbon(
    data = BHB_D47_temperature_model,
    aes(
      xmin = temperature_lower95_C,
      xmax = temperature_upper95_C,
      y = Age_Ma
    ),
    fill = "#B2182B", alpha = 0.16
  ) +
  geom_ribbon(
    data = BHB_D47_temperature_model,
    aes(
      xmin = temperature_lower80_C,
      xmax = temperature_upper80_C,
      y = Age_Ma
    ),
    fill = "#B2182B", alpha = 0.24
  ) +
  geom_path(
    data = BHB_D47_temperature_model %>% arrange(desc(Age_Ma)),
    aes(temperature_median_C, Age_Ma),
    color = "#B2182B", linewidth = 0.95
  ) +
  geom_errorbarh(
    data = BHB_D47_temperature_observations,
    aes(
      xmin = temperature_C - temperature_se_C,
      xmax = temperature_C + temperature_se_C,
      y = Age_Ma,
      color = section_id
    ),
    height = 0, linewidth = 0.3, alpha = 0.42
  ) +
  geom_point(
    data = BHB_D47_temperature_observations,
    aes(temperature_C, Age_Ma, color = section_id),
    size = 1.7, alpha = 0.72
  ) +
  scale_color_manual(values = c(CFB = "#B2182B", MCP = "#E66101")) +
  age_scale(regional_age_limits) +
  labs(
    x = expression(Delta[47] * " formation temperature (" * degree * "C)"),
    y = "Age (Ma)", color = "Section",
    title = "B  Soil-carbonate formation temperature",
    subtitle = "All age-resolved CFB and MCP data; 80% and 95% intervals"
  ) +
  theme_age

p_all_BHB_temperature <- ggplot() +
  add_petm_age() +
  geom_ribbon(
    data = CFB_temperature_age %>% filter(!is.na(T_model_C)),
    aes(
      xmin = T_model_lower95_C, xmax = T_model_upper95_C,
      y = Age_Ma
    ),
    fill = "grey70", alpha = 0.23
  ) +
  geom_path(
    data = CFB_temperature_age %>%
      filter(!is.na(T_model_C)) %>%
      arrange(desc(Age_Ma)),
    aes(T_model_C, Age_Ma),
    color = "black", linewidth = 0.95
  ) +
  geom_errorbarh(
    data = BHB_temperature_proxy_data,
    aes(
      xmin = temperature_C - temperature_se_C,
      xmax = temperature_C + temperature_se_C,
      y = Age_Ma,
      color = temperature_record
    ),
    height = 0, linewidth = 0.35, alpha = 0.55, na.rm = TRUE
  ) +
  geom_errorbar(
    data = BHB_temperature_proxy_data %>%
      filter(age_older_ma > age_younger_ma),
    aes(
      x = temperature_C,
      ymin = age_younger_ma,
      ymax = age_older_ma,
      color = temperature_record
    ),
    width = 0, linewidth = 0.45, alpha = 0.75
  ) +
  geom_point(
    data = BHB_temperature_proxy_data,
    aes(
      temperature_C, Age_Ma,
      color = temperature_record,
      shape = temperature_record
    ),
    fill = "white", size = 2.4, stroke = 0.75
  ) +
  scale_color_manual(values = temperature_record_colors) +
  scale_shape_manual(values = temperature_record_shapes) +
  guides(
    color = guide_legend(title = NULL, nrow = 1),
    shape = guide_legend(title = NULL, nrow = 1)
  ) +
  age_scale(regional_age_limits) +
  labs(
    x = expression("Temperature estimate (" * degree * "C)"),
    y = "Age (Ma)",
    color = "BHB temperature record",
    shape = "BHB temperature record",
    title = "A  Bighorn Basin temperature proxies",
    subtitle = paste(
      "Black curve: CFB D47 model; vertical bars show aggregate",
      "age ranges"
    )
  ) +
  theme_age

p_insolation_47N <- BHB_insolation_47N %>%
  filter(
    Age_Ma >= regional_age_limits[2],
    Age_Ma <= regional_age_limits[1]
  ) %>%
  ggplot(aes(summer_solstice_insolation_w_m2, Age_Ma)) +
  annotate(
    "rect", xmin = -Inf, xmax = Inf, ymin = 58, ymax = Inf,
    fill = "grey50", alpha = 0.10
  ) +
  add_petm_age(alpha = 0.10) +
  geom_path(color = "#2C7FB8", linewidth = 0.30, alpha = 0.85) +
  age_scale(regional_age_limits) +
  labs(
    x = expression("Summer-solstice insolation at 47" * degree * "N (W " * m^-2 * ")"),
    y = "Age (Ma)",
    title = "C  Modeled astronomical forcing",
    subtitle = "ZB20a(1,1); grey interval marks less-secure phase (>58 Ma)"
  ) +
  theme_age +
  theme(legend.position = "none")

p_BHB_temperature_insolation <-
  patchwork::wrap_plots(
    p_BHB_LMA_MAT_age,
    remove_repeated_age(
      p_BHB_D47_temperature_age + theme(legend.position = "top")
    ),
    remove_repeated_age(p_insolation_47N),
    nrow = 1,
    widths = c(1, 1.25, 1)
  ) +
  patchwork::plot_annotation(
    title = "Bighorn Basin temperature models and orbital-scale insolation",
    subtitle = paste(
      "MAT and carbonate formation temperature are modeled separately;",
      "the orbital solution is not tuned to either proxy record"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10)
    )
  )

p_regional_climate <-
  p_regional_temperature +
  remove_repeated_age(p_Harper_SST) +
  remove_repeated_age(p_Harper_CO2) +
  plot_layout(widths = c(1.2, 1, 1)) +
  plot_annotation(
    title = "Continental temperature, marine temperature, and atmospheric CO2",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

#-- 7.) Export Age-Domain Figures ------------------------------------------
save_age_plot(p_temperature_age, cfb_age_figure_dir,
              "CFB_temperature_age", 5.2, 7.5)
save_age_plot(p_d13C_age, cfb_age_figure_dir,
              "CFB_d13Ccarb_age", 5.2, 7.5)
save_age_plot(p_d18Owater_age, cfb_age_figure_dir,
              "CFB_d18Owater_age", 5.2, 7.5)
save_age_plot(p_D17Owater_age, cfb_age_figure_dir,
              "CFB_D17Owater_age", 5.2, 7.5)
save_age_plot(p_age_control, cfb_age_figure_dir,
              "CFB_chronology_support_age", 5.2, 7.5)
save_age_plot(p_CFB_age_multiproxy, cfb_age_figure_dir,
              "CFB_multiproxy_age_full", 16, 8)
save_age_plot(p_regional_temperature, regional_age_figure_dir,
              "BHB_regional_temperature_age", 6.2, 8)
save_age_plot(p_regional_climate, regional_age_figure_dir,
              "BHB_Harper_Tornillo_temperature_CO2_age", 12, 8)
save_age_plot(p_insolation_47N, regional_age_figure_dir,
              "BHB_ZB20a_summer_insolation_47N_age", 6.2, 8)
save_age_plot(p_BHB_LMA_MAT_age, regional_age_figure_dir,
              "BHB_LMA_MAT_model_age", 5.5, 8)
save_age_plot(p_BHB_D47_temperature_age, regional_age_figure_dir,
              "BHB_D47_formation_temperature_model_age", 6.2, 8)
save_age_plot(p_BHB_temperature_insolation, regional_age_figure_dir,
              "BHB_temperature_models_ZB20a_insolation_47N_age", 15.5, 8.5)

# Supporting table used in the combined figure. This makes explicit which
# observations, age ranges, proxy meanings, and uncertainty values were drawn.
write_csv(
  BHB_temperature_proxy_data,
  here("data", "processed", "BHB_temperature_proxy_plot_data.csv")
)
