# 10_plot_BHB_age_domain.R
# Purpose: Generate publication-ready age-domain figures for the CFB record
#          and explicitly separated regional/global reference datasets.

#-- 1.) Setup ---------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)
source(here("scripts", "helpers", "CFB_chronostrat_panels.R"))
source(here("scripts", "helpers", "save_figure_variants.R"))

cfb_age_figure_dir <- here("figures", "age_domain", "CFB")
regional_age_figure_dir <- here("figures", "age_domain", "regional_comparison")
presentation_figure_dir <- here("figures", "presentation")
dir.create(cfb_age_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(regional_age_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(presentation_figure_dir, recursive = TRUE, showWarnings = FALSE)

petm_age_old_ma <- 55.93
petm_age_young_ma <- 55.75
cfb_age_limits <- c(58.9, 54.0)
regional_age_limits <- c(59.0, 52.5)

source_colors <- c(
  "U-M" = "#000000", "CU" = "#0072B2", "Caltech" = "#D55E00",
  "Koch" = "#009E73", "Bowen" = "#CC79A7"
)
source_shapes <- c(
  "U-M" = 21, "CU" = 22, "Caltech" = 24, "Koch" = 23, "Bowen" = 25
)

t47_status_colors <- c(
  "This study" = "#B2182B",
  "Published data" = "grey25"
)
t47_status_shapes <- c("This study" = 21, "Published data" = 22)

add_petm_age <- function(alpha = 0.14) {
  annotate(
    "rect", xmin = -Inf, xmax = Inf,
    ymin = petm_age_young_ma, ymax = petm_age_old_ma,
    fill = "#E41A1C", alpha = alpha
  )
}

theme_age <- theme_classic(base_size = 18) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    plot.title = element_text(face = "bold", size = 18),
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

save_age_plot <- function(
    plot, directory, stem, width, height, presentation_width = NULL,
    presentation_plot = NULL
) {
  save_figure_variants(
    plot = plot,
    base_dir = directory,
    stem = stem,
    manuscript_width = width,
    manuscript_height = height,
    presentation_width = presentation_width,
    presentation_plot = presentation_plot
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
seasonal_synthesis_draws <- read_csv(
  here(
    "data", "processed",
    "BHB_seasonal_temperature_integrated_draws.csv.gz"
  ),
  show_col_types = FALSE
)
seasonal_synthesis_summary <- read_csv(
  here(
    "data", "processed",
    "BHB_seasonal_temperature_integrated_summary.csv"
  ),
  show_col_types = FALSE
)
seasonal_constraint_summary <- read_csv(
  here(
    "data", "processed",
    "BHB_seasonal_temperature_constraint_summary.csv"
  ),
  show_col_types = FALSE
)
Barnet2019_ODP1262 <- read_csv(
  here(
    "data", "processed",
    "BarnetEtAl2019_ODP1262_benthic_isotopes.csv"
  ),
  show_col_types = FALSE
)
Westerhold2018_ODP1209 <- read_csv(
  here(
    "data", "processed",
    "WesterholdEtAl2018_ODP1209_benthic_isotopes.csv"
  ),
  show_col_types = FALSE
)
BHB_regional_d13C_reference <- read_csv(
  here(
    "data", "processed",
    "BHB_regional_soilcarb_reference_summary_age_calibrated.csv"
  ),
  show_col_types = FALSE
)

age_lookup <- CFB_isotopes_age %>%
  distinct(section_id, MLA_horizon_id, strat_height_m, Age_Ma)

p_CFB_chronostrat_age <- build_CFB_chronostrat_age_panel(
  age_lookup = age_lookup,
  age_limits = cfb_age_limits
)

CFB_temperature_obs_age <- CFB_temperature_obs %>%
  filter(used_in_primary_temperature_model) %>%
  mutate(
    source = recode(source, IPL = "U-M", Snell = "Caltech"),
    study_status = if_else(source == "U-M", "This study", "Published data"),
    study_status = factor(
      study_status, levels = c("This study", "Published data")
    )
  ) %>%
  left_join(
    age_lookup,
    by = c("section_id", "MLA_horizon_id", "strat_height_m")
  ) %>%
  filter(!is.na(Age_Ma))

#-- 3.) Reshape CFB Carbon-Isotope Observations ----------------------------
d13C_age <- bind_rows(
  CFB_isotopes_age %>% transmute(
    MLA_horizon_id, Age_Ma, source = "U-M",
    value = IPL_NuDog_d13Ccarb_VPDB, se = NA_real_
  ),
  CFB_isotopes_age %>% transmute(
    MLA_horizon_id, Age_Ma, source = "CU",
    value = CU_mean_d13Ccarb_vpdb, se = NA_real_
  ),
  CFB_isotopes_age %>% transmute(
    MLA_horizon_id, Age_Ma, source = "Caltech",
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
  mutate(
    source = factor(source, levels = names(source_colors)),
    study_status = if_else(source == "U-M", "U-M / this study", "Published"),
    study_status = factor(
      study_status, levels = c("U-M / this study", "Published")
    )
  )

# Published BHB soil-carbonate carbon-isotope observations used for global
# comparison. U-M/IPL values are deliberately omitted: this panel provides an
# independent published terrestrial record alongside the marine reference.
BHB_d13C_published_age <- d13C_age %>%
  filter(source != "U-M") %>%
  transmute(
    section_id = "CFB",
    MLA_horizon_id,
    Age_Ma,
    source = as.character(source),
    d13Ccarb_vpdb = value,
    d13Ccarb_se = se
  ) %>%
  bind_rows(
    BHB_regional_d13C_reference %>%
      filter(!is.na(Age_Ma), !is.na(d13Ccarb_vpdb)) %>%
      transmute(
        section_id,
        MLA_horizon_id,
        Age_Ma,
        source = case_when(
          str_detect(dataset, regex("^Snell", ignore_case = TRUE)) ~
            "Caltech",
          str_detect(dataset, regex("^Koch", ignore_case = TRUE)) ~
            "Koch",
          TRUE ~ dataset
        ),
        d13Ccarb_vpdb,
        d13Ccarb_se = NA_real_
      )
  ) %>%
  distinct(
    section_id, MLA_horizon_id, Age_Ma, source, d13Ccarb_vpdb,
    .keep_all = TRUE
  ) %>%
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
      y = Age_Ma, color = study_status
    ),
    height = 0, linewidth = 0.35, alpha = 0.55
  ) +
  geom_point(
    data = CFB_temperature_obs_age,
    aes(
      T_C, Age_Ma, color = study_status,
      fill = p_altered_preservation, shape = study_status
    ),
    size = 2.5, stroke = 0.8
  ) +
  scale_color_manual(
    values = c("This study" = "#B2182B", "Published data" = "grey55"),
    drop = FALSE
  ) +
  scale_fill_gradientn(
    colors = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"),
    limits = c(0, 1), breaks = c(0, 0.5, 1),
    labels = scales::label_percent(accuracy = 1),
    name = "d18O trajectory\nP(altered)"
  ) +
  scale_shape_manual(values = t47_status_shapes, drop = FALSE) +
  age_scale(cfb_age_limits) +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = "Age (Ma)", color = NULL, shape = "Dataset",
    title = "A  CFB temperature"
  ) +
  guides(
    color = "none",
    shape = guide_legend(
      order = 1, override.aes = list(fill = "white", color = "black")
    ),
    fill = guide_colorbar(
      order = 2,
      barwidth = grid::unit(2.0, "cm"),
      barheight = grid::unit(0.35, "cm"),
      title.position = "top"
    )
  ) +
  theme_age +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.justification = "left"
  )

p_d13C_age <- ggplot(d13C_age) +
  add_petm_age() +
  geom_errorbarh(
    aes(xmin = value - se, xmax = value + se, y = Age_Ma, color = source),
    height = 0, linewidth = 0.3, alpha = 0.35, na.rm = TRUE
  ) +
  geom_point(
    aes(
      value, Age_Ma, color = study_status,
      fill = study_status, shape = source
    ),
    size = 2.1, alpha = 0.9, stroke = 0.7
  ) +
  scale_color_manual(
    values = c("U-M / this study" = "black", "Published" = "grey62"),
    drop = FALSE
  ) +
  scale_fill_manual(
    values = c("U-M / this study" = "black", "Published" = "white"),
    drop = FALSE
  ) +
  scale_shape_manual(values = source_shapes, drop = FALSE) +
  age_scale(cfb_age_limits) +
  labs(
    x = expression(delta^13 * C[carbonate] ~ "(per mil VPDB)"),
    y = "Age (Ma)", color = NULL, fill = "Dataset", shape = NULL,
    title = expression("B  CFB " * delta^13 * C[carbonate])
  ) +
  guides(
    color = "none",
    shape = "none",
    fill = guide_legend(override.aes = list(shape = 21))
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
  # Primary manuscript inference requires replicated, concordant D17O.
  # Provisional, heterogeneous, and pending horizons remain in the processed
  # audit tables but are not shown as primary observations.
  filter(
    D17O_primary_use == "primary",
    !is.na(D17Orsw_mean_permeg),
    !is.na(Age_Ma)
  ) %>%
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
  p_CFB_chronostrat_age +
  (p_temperature_age + theme(legend.position = "none")) +
  remove_repeated_age(p_d13C_age) +
  remove_repeated_age(p_d18Owater_age) +
  remove_repeated_age(p_D17Owater_age) +
  remove_repeated_age(p_age_control) +
  plot_layout(
    widths = c(0.92, 1.15, 1, 1, 1, 0.95),
    guides = "collect"
  ) &
  theme(legend.position = "top")

#-- 4b.) Global Marine-Terrestrial Carbon-Cycle Context -------------------
#
# Unlike the section-oriented panels above, this synthesis uses age on the
# horizontal axis, running from older at left to younger at right. That
# orientation makes temporal alignment among independent archives explicit
# and is preferable for event comparison, lead-lag interpretation, and
# conference-slide reading. The raw marine observations remain visible; a
# 50-kyr bin mean emphasizes the long-term trajectory without implying a
# fitted process model.

global_context_age_limits <- c(59.0, 52.5)

age_scale_horizontal <- scale_x_reverse(
  limits = global_context_age_limits,
  breaks = seq(59, 52.5, by = -0.5),
  minor_breaks = seq(59, 52.5, by = -0.1),
  expand = expansion(mult = c(0.005, 0.005))
)

add_petm_horizontal <- function(alpha = 0.11) {
  annotate(
    "rect",
    xmin = petm_age_old_ma,
    xmax = petm_age_young_ma,
    ymin = -Inf,
    ymax = Inf,
    fill = "#E41A1C",
    alpha = alpha
  )
}

Barnet2019_LPEE <- Barnet2019_ODP1262 %>%
  filter(
    Age_Ma <= global_context_age_limits[1],
    Age_Ma >= global_context_age_limits[2]
  )

Barnet2019_50kyr <- Barnet2019_LPEE %>%
  mutate(age_bin_ma = floor(Age_Ma / 0.05) * 0.05 + 0.025) %>%
  group_by(age_bin_ma) %>%
  summarise(
    Age_Ma = mean(Age_Ma, na.rm = TRUE),
    d13C_benthic_vpdb = mean(d13C_benthic_vpdb, na.rm = TRUE),
    bottom_water_temperature_C =
      mean(bottom_water_temperature_C, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(Age_Ma)

p_Barnet_BWT_horizontal <- ggplot(Barnet2019_LPEE) +
  add_petm_horizontal() +
  geom_line(
    aes(Age_Ma, bottom_water_temperature_C),
    color = "grey67", linewidth = 0.22, alpha = 0.65
  ) +
  geom_line(
    data = Barnet2019_50kyr,
    aes(Age_Ma, bottom_water_temperature_C),
    color = "#2166AC", linewidth = 0.8
  ) +
  age_scale_horizontal +
  labs(
    x = NULL,
    y = expression("Bottom-water " * T ~ "(" * degree * "C)"),
    title = "A  South Atlantic deep-ocean temperature",
    subtitle = "ODP Site 1262; grey = observations, blue = 50-kyr means"
  ) +
  theme_age +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_Barnet_d13C_horizontal <- ggplot(Barnet2019_LPEE) +
  add_petm_horizontal() +
  geom_line(
    aes(Age_Ma, d13C_benthic_vpdb),
    color = "grey67", linewidth = 0.22, alpha = 0.65
  ) +
  geom_line(
    data = Barnet2019_50kyr,
    aes(Age_Ma, d13C_benthic_vpdb),
    color = "#7F3B08", linewidth = 0.8
  ) +
  age_scale_horizontal +
  labs(
    x = NULL,
    y = expression(delta^13 * C[benthic] ~ "(per mil VPDB)"),
    title = "B  South Atlantic marine carbon cycle"
  ) +
  theme_age +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Standalone Atlantic-Pacific isotope comparison for a 6 x 3 inch slide.
# Site 1262 is the South Atlantic record presented by Barnet et al. (2019);
# Site 1209 is the Pacific comparison record used in that paper.
barnet_full_age_limits <- c(67.1, 52.35)

barnet_atlantic_pacific <- bind_rows(
  Barnet2019_ODP1262 %>%
    transmute(
      Age_Ma,
      basin = "Atlantic",
      d13C_benthic_vpdb,
      d18O_benthic_vpdb = d18O_corrected_vpdb
    ),
  Westerhold2018_ODP1209 %>%
    transmute(Age_Ma, basin, d13C_benthic_vpdb, d18O_benthic_vpdb)
) %>%
  mutate(basin = factor(basin, levels = c("Atlantic", "Pacific"))) %>%
  arrange(basin, Age_Ma)

barnet_full_age_scale <- scale_x_reverse(
  limits = barnet_full_age_limits,
  breaks = seq(68, 52, by = -2),
  minor_breaks = seq(68, 52, by = -0.5),
  expand = expansion(mult = c(0.005, 0.005))
)

barnet_basin_linetypes <- c("Atlantic" = "solid", "Pacific" = "22")
barnet_d13C_colors <- c("Atlantic" = "#F03B20", "Pacific" = "#A50056")
barnet_d18O_colors <- c("Atlantic" = "#29A9E0", "Pacific" = "#173F90")

barnet_full_panel_theme <- theme_classic(base_size = 18) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18, color = "black"),
    plot.title = element_text(size = 18, face = "bold"),
    # The 6 x 3 export needs explicit room for 18-point vertical axis titles.
    plot.margin = margin(1, 4, 1, 28),
    axis.title.y = element_text(size = 18, margin = margin(r = 4)),
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    legend.key.width = unit(0.45, "cm"),
    legend.key.height = unit(0.28, "cm"),
    legend.spacing.x = unit(0.08, "cm"),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.86), color = "grey55", linewidth = 0.25
    )
  )

p_Barnet_d13C_full_horizontal <-
  ggplot(
    barnet_atlantic_pacific,
    aes(Age_Ma, d13C_benthic_vpdb, color = basin, linetype = basin)
  ) +
  geom_line(linewidth = 0.38, alpha = 0.9, na.rm = TRUE) +
  barnet_full_age_scale +
  scale_color_manual(values = barnet_d13C_colors) +
  scale_linetype_manual(values = barnet_basin_linetypes) +
  labs(
    x = NULL,
    # Keep the 18-point isotope label short enough for the 6 x 3 panel.
    y = expression(delta^13 * C)
  ) +
  barnet_full_panel_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.88, 0.76),
    legend.direction = "vertical"
  )

p_Barnet_d18O_full_horizontal <-
  ggplot(
    barnet_atlantic_pacific,
    aes(Age_Ma, d18O_benthic_vpdb, color = basin, linetype = basin)
  ) +
  geom_line(linewidth = 0.38, alpha = 0.9, na.rm = TRUE) +
  barnet_full_age_scale +
  scale_y_reverse() +
  scale_color_manual(values = barnet_d18O_colors) +
  scale_linetype_manual(values = barnet_basin_linetypes) +
  labs(
    x = NULL,
    y = expression(delta^18 * O)
  ) +
  barnet_full_panel_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

barnet_epochs <- tribble(
  ~unit,          ~older_ma, ~younger_ma, ~fill,
  "Late\nCret.",     67.10,       66.00,  "#B8DE8A",
  "Paleocene",       66.00,       56.00,  "#F6E8A6",
  "Eocene",          56.00,       52.35,  "#F7B267"
)

barnet_stages <- tribble(
  ~unit,       ~older_ma, ~younger_ma, ~fill,
  "Maast.",         67.10,       66.00,  "#D8E7A8",
  "Danian",         66.00,       61.66,  "#B8DE8A",
  "Sel.",           61.66,       59.24,  "#8FD080",
  "Thanetian",      59.24,       56.00,  "#63C178",
  "Ypresian",       56.00,       52.35,  "#F4C96B"
)

p_Barnet_chronostrat_full_horizontal <-
  ggplot() +
  geom_rect(
    data = barnet_epochs,
    aes(
      xmin = older_ma,
      xmax = younger_ma,
      ymin = 1,
      ymax = 2,
      fill = fill
    ),
    color = "grey25",
    linewidth = 0.35
  ) +
  geom_rect(
    data = barnet_stages,
    aes(
      xmin = older_ma,
      xmax = younger_ma,
      ymin = 0,
      ymax = 1,
      fill = fill
    ),
    color = "grey25",
    linewidth = 0.35
  ) +
  geom_text(
    data = barnet_epochs %>% filter(unit != "Late\nCret."),
    aes(
      x = (older_ma + younger_ma) / 2,
      y = 1.5,
      label = unit
    ),
    size = 18 / ggplot2::.pt
  ) +
  geom_text(
    data = barnet_epochs %>% filter(unit == "Late\nCret."),
    aes(
      x = (older_ma + younger_ma) / 2,
      y = 1.5,
      label = unit
    ),
    size = 10 / ggplot2::.pt,
    lineheight = 0.85
  ) +
  geom_text(
    data = barnet_stages %>% filter(unit != "Maast."),
    aes(
      x = (older_ma + younger_ma) / 2,
      y = 0.5,
      label = unit
    ),
    size = 18 / ggplot2::.pt
  ) +
  geom_text(
    data = barnet_stages %>% filter(unit == "Maast."),
    aes(
      x = (older_ma + younger_ma) / 2,
      y = 0.5,
      label = unit
    ),
    # The 1.1 Myr Maastrichtian sliver is the constrained-strat exception.
    size = 10 / ggplot2::.pt
  ) +
  barnet_full_age_scale +
  scale_fill_identity() +
  scale_y_continuous(
    limits = c(0, 2),
    expand = c(0, 0)
  ) +
  labs(
    x = "Age (Ma)",
    y = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_text(
      size = 18,
      margin = margin(t = 2)
    ),
    axis.text.x = element_text(
      size = 18,
      color = "black"
    ),
    plot.margin = margin(0, 4, 1, 4)
  )

p_Barnet_stable_isotopes_full_horizontal_6x3 <-
  p_Barnet_d13C_full_horizontal /
  p_Barnet_d18O_full_horizontal /
  p_Barnet_chronostrat_full_horizontal +
  plot_layout(
    heights = c(1, 1, 0.75)
  )

p_BHB_d13C_horizontal <- ggplot(BHB_d13C_published_age) +
  add_petm_horizontal() +
  geom_errorbar(
    aes(
      x = Age_Ma,
      ymin = d13Ccarb_vpdb - d13Ccarb_se,
      ymax = d13Ccarb_vpdb + d13Ccarb_se,
      color = source
    ),
    width = 0,
    linewidth = 0.35,
    alpha = 0.45,
    na.rm = TRUE
  ) +
  geom_point(
    aes(Age_Ma, d13Ccarb_vpdb, color = source, shape = source),
    size = 1.8,
    alpha = 0.82,
    fill = "white",
    stroke = 0.6
  ) +
  age_scale_horizontal +
  scale_color_manual(values = source_colors, drop = TRUE) +
  scale_shape_manual(values = source_shapes, drop = TRUE) +
  labs(
    x = "Age (Ma; older to younger)",
    y = expression(delta^13 * C[soil~carbonate] ~ "(per mil VPDB)"),
    color = "Published source",
    shape = "Published source",
    title = expression("C  Bighorn Basin soil-carbonate " * delta^13 * C),
    subtitle = paste(
      "Published CFB and regional-reference values from CU, Caltech, Koch,",
      "and Bowen; U-M/IPL omitted"
    )
  ) +
  theme_age +
  theme(legend.position = "bottom")

p_global_Barnet_BHB_d13C <-
  p_Barnet_BWT_horizontal /
  p_Barnet_d13C_horizontal /
  p_BHB_d13C_horizontal +
  plot_layout(heights = c(0.9, 0.9, 1.15), guides = "collect") +
  plot_annotation(
    title = "Global marine climate and carbon-cycle context for the Bighorn Basin",
    subtitle = paste(
      "Marine record: Barnet et al. (2019), ODP Site 1262;",
      "terrestrial panel excludes U-M/IPL data"
    ),
    caption = paste(
      "Barnet et al. data: PANGAEA 10.1594/PANGAEA.884585.",
      "PETM band = 55.93-55.75 Ma."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 18),
      plot.caption = element_text(size = 18, hjust = 0)
    )
  ) &
  theme(legend.position = "bottom")

# Presentation version: this is deliberately composed independently rather
# than resizing the publication figure. It is designed for a 12 x 6 inch
# placement on a 16:9 slide. Secondary titles are removed, typography and
# symbols are scaled for that physical size, and margins are tightened. The
# analytical content, age limits, PETM interval, and observations are unchanged.

slide_panel_theme <- theme(
  axis.title = element_text(size = 18),
  axis.text = element_text(size = 18, color = "black"),
  plot.title = element_text(face = "bold", size = 18),
  plot.subtitle = element_text(size = 18),
  plot.margin = margin(2, 6, 2, 6),
  legend.title = element_text(size = 18),
  legend.text = element_text(size = 18)
)

p_Barnet_BWT_slide <- p_Barnet_BWT_horizontal +
  labs(
    title = expression(
      "A  South Atlantic bottom-water temperature (" * degree * "C)"
    ),
    subtitle = NULL,
    y = NULL
  ) +
  slide_panel_theme

p_Barnet_d13C_slide <- p_Barnet_d13C_horizontal +
  labs(
    title = expression(
      "B  South Atlantic benthic " * delta^13 * C ~ "(per mil VPDB)"
    ),
    y = NULL
  ) +
  slide_panel_theme

p_BHB_d13C_slide <- p_BHB_d13C_horizontal +
  labs(
    title = expression(
      "C  Bighorn Basin soil-carbonate " * delta^13 * C ~ "(per mil VPDB)"
    ),
    subtitle = NULL,
    y = NULL
  ) +
  slide_panel_theme +
  theme(
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-3, 0, 0, 0)
  )

p_global_Barnet_BHB_d13C_slide <-
  p_Barnet_BWT_slide /
  p_Barnet_d13C_slide /
  p_BHB_d13C_slide +
  plot_layout(heights = c(0.82, 0.82, 1.08), guides = "collect") +
  plot_annotation(
    title = "Global marine context and Bighorn Basin carbon-cycle change",
    subtitle = paste(
      "Barnet et al. (2019), ODP Site 1262;",
      "grey = observations, colored marine curves = 50-kyr means;",
      "U-M/IPL omitted"
    ),
    caption = paste(
      "Marine data: PANGAEA 10.1594/PANGAEA.884585.",
      "Pink band: PETM, 55.93-55.75 Ma."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 20),
      plot.subtitle = element_text(size = 18),
      plot.caption = element_text(size = 18, hjust = 0),
      plot.margin = margin(8, 10, 5, 10)
    )
  ) &
  theme(legend.position = "bottom")

# Minimal-title, 12 x 6 inch version for direct placement on a presentation
# slide. Interpretive framing belongs in the slide title or spoken narrative;
# the figure itself retains only the panel identities and compact provenance.
p_global_Barnet_BHB_d13C_slide_12x6 <-
  p_Barnet_BWT_slide /
  p_Barnet_d13C_slide /
  p_BHB_d13C_slide +
  plot_layout(heights = c(0.82, 0.82, 1.08), guides = "collect") +
  plot_annotation(
    caption = paste(
      "Barnet et al. (2019), ODP Site 1262;",
      "grey = observations; colored marine curves = 50-kyr means;",
      "pink = PETM; U-M/IPL omitted.",
      "Data: PANGAEA 10.1594/PANGAEA.884585."
    ),
    theme = theme(
      plot.caption = element_text(size = 18, hjust = 0),
      plot.margin = margin(4, 8, 3, 8)
    )
  ) &
  theme(legend.position = "bottom")

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
) %>%
  filter(used_in_temperature_model)

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
      color = study_status
    ),
    height = 0, linewidth = 0.3, alpha = 0.42
  ) +
  geom_point(
    data = BHB_D47_temperature_observations,
    aes(
      temperature_C, Age_Ma,
      color = study_status, shape = study_status, fill = study_status
    ),
    size = 2.2, stroke = 0.75, alpha = 0.9
  ) +
  scale_color_manual(values = t47_status_colors) +
  scale_shape_manual(values = t47_status_shapes) +
  scale_fill_manual(
    values = c("This study" = "#B2182B", "Published data" = "white")
  ) +
  age_scale(regional_age_limits) +
  guides(
    color = "none",
    shape = "none",
    fill = guide_legend(title = "Dataset", override.aes = list(shape = 21))
  ) +
  labs(
    x = expression(Delta[47] * " formation temperature (" * degree * "C)"),
    y = "Age (Ma)", color = NULL, shape = NULL,
    title = "B  Soil-carbonate formation temperature",
    subtitle = "This study versus published data; 80% and 95% intervals"
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

#-- 6b.) Seasonal-Climate Context for BHB Delta47 --------------------------
# These are age-invariant literature-informed marginal distributions, not
# additional time-series observations. Nested vertical bands show the central
# 95%, 80%, and 50% ranges; the dark line is the median. CMMT, MAAT, and WMMT
# are absolute temperatures and can be compared with Delta47. MART is a
# temperature difference, so it is retained in the companion evidence panel
# rather than drawn as an absolute-temperature band.
seasonal_colors <- c(
  CMMT = "#2166AC",
  MAAT = "#8856A7",
  WMMT = "#D73027",
  MART = "#E08214"
)

seasonal_absolute_summary <- seasonal_synthesis_summary %>%
  filter(metric %in% c("CMMT", "MAAT", "WMMT")) %>%
  mutate(metric = factor(metric, levels = c("CMMT", "MAAT", "WMMT")))

add_seasonal_reference_bands <- function(age_limits) {
  list(
    geom_rect(
      data = seasonal_absolute_summary,
      aes(
        xmin = lower95_c, xmax = upper95_c,
        ymin = min(age_limits), ymax = max(age_limits),
        fill = metric
      ),
      alpha = 0.075, inherit.aes = FALSE
    ),
    geom_rect(
      data = seasonal_absolute_summary,
      aes(
        xmin = lower80_c, xmax = upper80_c,
        ymin = min(age_limits), ymax = max(age_limits),
        fill = metric
      ),
      alpha = 0.10, inherit.aes = FALSE
    ),
    geom_rect(
      data = seasonal_absolute_summary,
      aes(
        xmin = lower50_c, xmax = upper50_c,
        ymin = min(age_limits), ymax = max(age_limits),
        fill = metric
      ),
      alpha = 0.14, inherit.aes = FALSE
    ),
    geom_vline(
      data = seasonal_absolute_summary,
      aes(xintercept = mean_c, color = metric),
      linewidth = 0.65, alpha = 0.85, show.legend = FALSE
    )
  )
}

p_BHB_D47_seasonal_context <- ggplot() +
  add_seasonal_reference_bands(regional_age_limits) +
  add_petm_age(alpha = 0.10) +
  geom_ribbon(
    data = BHB_D47_temperature_model,
    aes(
      xmin = temperature_lower95_C,
      xmax = temperature_upper95_C,
      y = Age_Ma
    ),
    fill = "grey45", alpha = 0.16
  ) +
  geom_ribbon(
    data = BHB_D47_temperature_model,
    aes(
      xmin = temperature_lower80_C,
      xmax = temperature_upper80_C,
      y = Age_Ma
    ),
    fill = "grey35", alpha = 0.18
  ) +
  geom_path(
    data = BHB_D47_temperature_model %>% arrange(desc(Age_Ma)),
    aes(temperature_median_C, Age_Ma),
    color = "black", linewidth = 0.95
  ) +
  geom_errorbarh(
    data = BHB_D47_temperature_observations,
    aes(
      xmin = temperature_C - temperature_se_C,
      xmax = temperature_C + temperature_se_C,
      y = Age_Ma,
      color = study_status
    ),
    height = 0, linewidth = 0.3, alpha = 0.45
  ) +
  geom_point(
    data = BHB_D47_temperature_observations %>%
      filter(study_status == "Published data"),
    aes(temperature_C, Age_Ma),
    shape = 21, color = "grey50", fill = "white",
    size = 2, stroke = 0.7
  ) +
  geom_point(
    data = BHB_D47_temperature_observations %>%
      filter(study_status == "This study"),
    aes(temperature_C, Age_Ma),
    shape = 21, color = "#B2182B", fill = "#B2182B",
    size = 2.25, stroke = 0.7
  ) +
  scale_fill_manual(values = seasonal_colors, drop = FALSE) +
  scale_color_manual(
    values = c(
      seasonal_colors,
      t47_status_colors
    )
  ) +
  scale_shape_manual(values = t47_status_shapes) +
  age_scale(regional_age_limits) +
  guides(
    fill = guide_legend(
      title = "Literature synthesis",
      override.aes = list(alpha = 0.18)
    ),
    color = "none",
    shape = "none"
  ) +
  labs(
    x = expression(Delta[47] * " formation temperature (" * degree * "C)"),
    y = "Age (Ma)",
    title = "BHB soil-carbonate temperatures in seasonal-climate context",
    subtitle = "Lines = means; nested bands = central 50%, 80%, and 95%."
  ) +
  theme_age

# Single-panel synthesis of all BHB temperature proxy observations. The
# soil-carbonate Delta47 model is retained, but published and new Delta47
# observations share one symbol because their study provenance is not the
# comparison of interest here. Leaf-margin and Fricke-Wing estimates are shown
# as observations only; no LMA MAT model is drawn.
p_BHB_all_temperature_seasonal_context <- ggplot() +
  add_seasonal_reference_bands(regional_age_limits) +
  add_petm_age(alpha = 0.10) +
  geom_ribbon(
    data = BHB_D47_temperature_model,
    aes(
      xmin = temperature_lower95_C,
      xmax = temperature_upper95_C,
      y = Age_Ma
    ),
    fill = "grey45", alpha = 0.16
  ) +
  geom_ribbon(
    data = BHB_D47_temperature_model,
    aes(
      xmin = temperature_lower80_C,
      xmax = temperature_upper80_C,
      y = Age_Ma
    ),
    fill = "grey35", alpha = 0.18
  ) +
  geom_path(
    data = BHB_D47_temperature_model %>% arrange(desc(Age_Ma)),
    aes(temperature_median_C, Age_Ma),
    color = "black", linewidth = 0.95
  ) +
  geom_errorbarh(
    data = BHB_D47_temperature_observations,
    aes(
      xmin = temperature_C - temperature_se_C,
      xmax = temperature_C + temperature_se_C,
      y = Age_Ma
    ),
    height = 0, linewidth = 0.32, color = "#B2182B",
    alpha = 0.48
  ) +
  geom_point(
    data = BHB_D47_temperature_observations,
    aes(temperature_C, Age_Ma, shape = "Soil-carbonate D47"),
    color = "#B2182B", fill = "#B2182B",
    size = 2.35, stroke = 0.7
  ) +
  geom_errorbarh(
    data = BHB_LMA_MAT_observations,
    aes(
      xmin = MAT_C - MAT_se_C,
      xmax = MAT_C + MAT_se_C,
      y = Age_Ma
    ),
    height = 0, linewidth = 0.32, color = "#1B7837",
    alpha = 0.55
  ) +
  geom_errorbar(
    data = BHB_LMA_MAT_observations,
    aes(x = MAT_C, ymin = age_younger_ma, ymax = age_older_ma),
    width = 0, linewidth = 0.38, color = "#1B7837",
    alpha = 0.55
  ) +
  geom_point(
    data = BHB_LMA_MAT_observations,
    aes(MAT_C, Age_Ma, shape = "Wing leaf-margin MAT"),
    color = "#1B7837", fill = "white",
    size = 2.35, stroke = 0.75
  ) +
  geom_errorbar(
    data = FrickeWing2004_BHB_MAAT,
    aes(
      x = temperature_C,
      ymin = age_younger_ma,
      ymax = age_older_ma
    ),
    width = 0, linewidth = 0.38, color = "#762A83",
    alpha = 0.55
  ) +
  geom_point(
    data = FrickeWing2004_BHB_MAAT,
    aes(
      temperature_C, Age_Ma,
      shape = "Fricke-Wing phosphate/LMA MAAT"
    ),
    color = "#762A83", fill = "white",
    size = 2.35, stroke = 0.75
  ) +
  scale_fill_manual(values = seasonal_colors, drop = FALSE) +
  scale_color_manual(values = seasonal_colors) +
  scale_shape_manual(
    values = c(
      "Soil-carbonate D47" = 21,
      "Wing leaf-margin MAT" = 22,
      "Fricke-Wing phosphate/LMA MAAT" = 23
    ),
    breaks = c(
      "Soil-carbonate D47",
      "Wing leaf-margin MAT",
      "Fricke-Wing phosphate/LMA MAAT"
    )
  ) +
  age_scale(regional_age_limits) +
  guides(
    fill = guide_legend(
      title = "Literature seasonal synthesis",
      override.aes = list(alpha = 0.18),
      order = 2
    ),
    color = "none",
    shape = guide_legend(
      title = "Temperature proxy", order = 1,
      nrow = 1, byrow = TRUE
    )
  ) +
  labs(
    x = expression("Temperature (" * degree * "C)"),
    y = "Age (Ma)",
    title = "Bighorn Basin temperature proxies and seasonal context",
    subtitle = paste(
      "Black curve and grey bands: soil-carbonate D47 model;",
      "LMA observations are not modeled"
    )
  ) +
  theme_age +
  theme(
    legend.position = "top",
    legend.box = "vertical",
    legend.justification = "left",
    plot.margin = margin(9, 7, 5, 7)
  )

# Empirical phase distributions propagate each reported observation SE and
# give every horizon equal prior weight. They show the sampled T47 population,
# not an estimate of regional air temperature or a time-weighted climate mean.
set.seed(47024)
phase_draws <- BHB_D47_temperature_observations %>%
  mutate(
    petm_phase = case_when(
      Age_Ma > petm_age_old_ma ~ "Pre-PETM",
      Age_Ma >= petm_age_young_ma ~ "PETM",
      TRUE ~ "Post-PETM"
    ),
    petm_phase = factor(
      petm_phase,
      levels = c("Post-PETM", "PETM", "Pre-PETM")
    )
  ) %>%
  select(section_id, MLA_horizon_id, petm_phase, temperature_C,
         temperature_se_C) %>%
  mutate(draws = map2(
    temperature_C,
    temperature_se_C,
    ~ rnorm(2000, mean = .x, sd = .y)
  )) %>%
  unnest_longer(draws, values_to = "temperature_draw_c")

phase_summary <- phase_draws %>%
  group_by(petm_phase) %>%
  summarise(
    mean_c = mean(temperature_draw_c),
    median_c = median(temperature_draw_c),
    lower50_c = quantile(temperature_draw_c, 0.25),
    upper50_c = quantile(temperature_draw_c, 0.75),
    lower80_c = quantile(temperature_draw_c, 0.10),
    upper80_c = quantile(temperature_draw_c, 0.90),
    lower95_c = quantile(temperature_draw_c, 0.025),
    upper95_c = quantile(temperature_draw_c, 0.975),
    n_horizons = n_distinct(paste(section_id, MLA_horizon_id)),
    .groups = "drop"
  )

p_T47_petm_phase_context <- ggplot(
  phase_draws %>% group_by(petm_phase) %>% slice_sample(n = 12000),
  aes(temperature_draw_c, petm_phase)
) +
  geom_rect(
    data = seasonal_absolute_summary,
    aes(
      xmin = lower95_c, xmax = upper95_c,
      ymin = -Inf, ymax = Inf, fill = metric
    ),
    alpha = 0.065, inherit.aes = FALSE
  ) +
  geom_rect(
    data = seasonal_absolute_summary,
    aes(
      xmin = lower80_c, xmax = upper80_c,
      ymin = -Inf, ymax = Inf, fill = metric
    ),
    alpha = 0.085, inherit.aes = FALSE
  ) +
  geom_rect(
    data = seasonal_absolute_summary,
    aes(
      xmin = lower50_c, xmax = upper50_c,
      ymin = -Inf, ymax = Inf, fill = metric
    ),
    alpha = 0.11, inherit.aes = FALSE
  ) +
  geom_vline(
    data = seasonal_absolute_summary,
    aes(xintercept = mean_c, color = metric),
    linewidth = 0.7, show.legend = FALSE
  ) +
  geom_violin(
    fill = "grey70", color = "grey25",
    alpha = 0.62, scale = "width", trim = FALSE
  ) +
  geom_segment(
    data = phase_summary,
    aes(
      x = lower95_c, xend = upper95_c,
      y = petm_phase, yend = petm_phase
    ),
    linewidth = 1.1, color = "grey55", inherit.aes = FALSE
  ) +
  geom_segment(
    data = phase_summary,
    aes(
      x = lower80_c, xend = upper80_c,
      y = petm_phase, yend = petm_phase
    ),
    linewidth = 3.0, color = "grey30", inherit.aes = FALSE
  ) +
  geom_segment(
    data = phase_summary,
    aes(
      x = lower50_c, xend = upper50_c,
      y = petm_phase, yend = petm_phase
    ),
    linewidth = 6.0, color = "black", inherit.aes = FALSE
  ) +
  geom_point(
    data = phase_summary,
    aes(mean_c, petm_phase),
    shape = 21, fill = "white", size = 2.5,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = seasonal_colors) +
  scale_color_manual(values = seasonal_colors) +
  labs(
    x = expression(Delta[47] * " formation temperature (" * degree * "C)"),
    y = NULL,
    fill = "Seasonal synthesis",
    title = "BHB T47 before, during, and after the PETM",
    subtitle = paste(
      "Violin = propagated T47; white point = mean;",
      "bars = central 50%, 80%, and 95%."
    )
  ) +
  theme_classic(base_size = 18) +
  theme(legend.position = "top")

p_seasonal_evidence_forest <- seasonal_constraint_summary %>%
  mutate(
    metric = factor(metric, levels = c("MAAT", "CMMT", "WMMT", "MART")),
    evidence_label = paste(study, scenario, sep = " — ")
  ) %>%
  ggplot(
    aes(
      median_c,
      reorder(evidence_label, median_c),
      color = evidence_class
    )
  ) +
  geom_errorbarh(
    aes(xmin = lower95_c, xmax = upper95_c),
    height = 0, linewidth = 0.45
  ) +
  geom_point(size = 1.7) +
  facet_wrap(~metric, scales = "free", ncol = 2) +
  labs(
    x = "Temperature or temperature range (degrees C)",
    y = NULL,
    color = "Evidence",
    title = "Published proxy and model constraints",
    subtitle = "MART is a range; the other metrics are absolute temperatures"
  ) +
  theme_classic(base_size = 18) +
  theme(legend.position = "bottom")

p_T47_petm_phase_with_evidence <-
  p_T47_petm_phase_context +
  p_seasonal_evidence_forest +
  plot_layout(widths = c(1.0, 1.35))

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
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 18)
    )
  )

p_regional_climate <-
  p_regional_temperature +
  remove_repeated_age(p_Harper_SST) +
  remove_repeated_age(p_Harper_CO2) +
  plot_layout(widths = c(1.2, 1, 1)) +
  plot_annotation(
    title = "Continental temperature, marine temperature, and atmospheric CO2",
    theme = theme(plot.title = element_text(face = "bold", size = 18))
  )

#-- 7.) Build Horizontal-Age and PETM-Focus Variants ----------------------
# These variants preserve the same observations, uncertainty, and fitted
# models as the vertical-age figures. Only the coordinate orientation and
# displayed age window change. Full-range plots place older ages at left and
# younger ages at right; PETM-focus plots use the common 56.5--55.5 Ma window.
horizontal_age_plot <- function(plot, limits = cfb_age_limits) {
  plot +
    scale_y_reverse(
      limits = limits,
      breaks = seq(ceiling(limits[2] * 2) / 2,
                   floor(limits[1] * 2) / 2, by = 0.5),
      minor_breaks = seq(ceiling(limits[2] * 10) / 10,
                         floor(limits[1] * 10) / 10, by = 0.1),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    coord_flip() +
    labs(x = NULL, y = "Age (Ma)") +
    theme(
      legend.position = "top",
      plot.title.position = "plot"
    )
}

petm_age_limits <- c(56.5, 55.5)

p_temperature_age_horizontal <- horizontal_age_plot(p_temperature_age)
p_d13C_age_horizontal <- horizontal_age_plot(p_d13C_age)
p_d18Owater_age_horizontal <- horizontal_age_plot(p_d18Owater_age)
p_D17Owater_age_horizontal <- horizontal_age_plot(p_D17Owater_age)

p_temperature_age_PETM_horizontal <-
  horizontal_age_plot(p_temperature_age, petm_age_limits)
p_d13C_age_PETM_horizontal <-
  horizontal_age_plot(p_d13C_age, petm_age_limits)
p_d18Owater_age_PETM_horizontal <-
  horizontal_age_plot(p_d18Owater_age, petm_age_limits)
p_D17Owater_age_PETM_horizontal <-
  horizontal_age_plot(p_D17Owater_age, petm_age_limits)

# Standalone vertical-age counterpart to the horizontal PETM temperature plot.
# This retains temperature on x and age on y for direct comparison with the
# stratigraphic-style age figures used elsewhere in the project.
p_temperature_age_PETM_vertical <-
  p_temperature_age +
  age_scale(petm_age_limits) +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = "Age (Ma)",
    title = "CFB temperature across the PETM"
  )

# Regional single-panel records use the same age-on-x convention so they can
# be aligned directly with the CFB proxy panels and marine reference records.
p_insolation_47N_horizontal <-
  horizontal_age_plot(p_insolation_47N, regional_age_limits)
p_BHB_LMA_MAT_age_horizontal <-
  horizontal_age_plot(p_BHB_LMA_MAT_age, regional_age_limits)
p_BHB_D47_temperature_age_horizontal <-
  horizontal_age_plot(p_BHB_D47_temperature_age, regional_age_limits)
p_BHB_D47_seasonal_context_horizontal <-
  horizontal_age_plot(p_BHB_D47_seasonal_context, regional_age_limits)
p_BHB_D47_temperature_age_PETM_horizontal <-
  horizontal_age_plot(p_BHB_D47_temperature_age, petm_age_limits)
p_BHB_D47_seasonal_context_PETM_horizontal <-
  horizontal_age_plot(p_BHB_D47_seasonal_context, petm_age_limits)

p_CFB_age_multiproxy_horizontal <- wrap_plots(
  p_temperature_age_horizontal,
  p_d13C_age_horizontal,
  p_d18Owater_age_horizontal,
  p_D17Owater_age_horizontal,
  ncol = 2,
  guides = "collect"
) &
  theme(legend.position = "top")

p_CFB_age_multiproxy_PETM_horizontal <- wrap_plots(
  p_temperature_age_PETM_horizontal,
  p_d13C_age_PETM_horizontal,
  p_d18Owater_age_PETM_horizontal,
  p_D17Owater_age_PETM_horizontal,
  ncol = 2,
  guides = "collect"
) &
  theme(legend.position = "top")

#-- 8.) Export Age-Domain Figures ------------------------------------------
with_age_chronostrat <- function(plot) {
  p_CFB_chronostrat_age + plot +
    plot_layout(widths = c(0.92, 2.1))
}

save_age_plot(with_age_chronostrat(p_temperature_age), cfb_age_figure_dir,
              "CFB_temperature_age", 7.5, 7.5)
save_age_plot(with_age_chronostrat(p_d13C_age), cfb_age_figure_dir,
              "CFB_d13Ccarb_age", 7.5, 7.5)
save_age_plot(with_age_chronostrat(p_d18Owater_age), cfb_age_figure_dir,
              "CFB_d18Owater_age", 7.5, 7.5)
save_age_plot(with_age_chronostrat(p_D17Owater_age), cfb_age_figure_dir,
              "CFB_D17Owater_age", 7.5, 7.5)
save_age_plot(with_age_chronostrat(p_age_control), cfb_age_figure_dir,
              "CFB_chronology_support_age", 7.5, 7.5)
save_age_plot(p_CFB_age_multiproxy, cfb_age_figure_dir,
              "CFB_multiproxy_age_full", 18, 8)
save_age_plot(p_temperature_age_horizontal, cfb_age_figure_dir,
              "CFB_temperature_age_horizontal", 7.5, 5.5,
              presentation_width = 6)
save_age_plot(p_d13C_age_horizontal, cfb_age_figure_dir,
              "CFB_d13Ccarb_age_horizontal", 7.5, 5.5,
              presentation_width = 6)
save_age_plot(p_d18Owater_age_horizontal, cfb_age_figure_dir,
              "CFB_d18Owater_age_horizontal", 7.5, 5.5,
              presentation_width = 6)
save_age_plot(p_D17Owater_age_horizontal, cfb_age_figure_dir,
              "CFB_D17Owater_age_horizontal", 7.5, 5.5,
              presentation_width = 6)
save_age_plot(p_CFB_age_multiproxy_horizontal, cfb_age_figure_dir,
              "CFB_multiproxy_age_horizontal", 12, 8,
              presentation_width = 12,
              presentation_plot =
                add_presentation_theme(p_CFB_age_multiproxy_horizontal) &
                theme(legend.position = "none"))
save_age_plot(p_temperature_age_PETM_horizontal, cfb_age_figure_dir,
              "CFB_temperature_age_PETM_56.5-55.5Ma", 7.5, 5.5,
              presentation_width = 6)
save_age_plot(p_temperature_age_PETM_vertical, cfb_age_figure_dir,
              "CFB_temperature_age_PETM_56.5-55.5Ma_vertical", 6, 7.5,
              presentation_width = 6)
save_age_plot(p_d13C_age_PETM_horizontal, cfb_age_figure_dir,
              "CFB_d13Ccarb_age_PETM_56.5-55.5Ma", 7.5, 5.5,
              presentation_width = 6)
save_age_plot(p_d18Owater_age_PETM_horizontal, cfb_age_figure_dir,
              "CFB_d18Owater_age_PETM_56.5-55.5Ma", 7.5, 5.5,
              presentation_width = 6)
save_age_plot(p_D17Owater_age_PETM_horizontal, cfb_age_figure_dir,
              "CFB_D17Owater_age_PETM_56.5-55.5Ma", 7.5, 5.5,
              presentation_width = 6)
save_age_plot(p_CFB_age_multiproxy_PETM_horizontal, cfb_age_figure_dir,
              "CFB_multiproxy_age_PETM_56.5-55.5Ma", 12, 8,
              presentation_width = 12,
              presentation_plot =
                add_presentation_theme(
                  p_CFB_age_multiproxy_PETM_horizontal
                ) & theme(legend.position = "none"))
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
save_age_plot(p_BHB_D47_seasonal_context, regional_age_figure_dir,
              "BHB_D47_seasonal_synthesis_context_age", 7.4, 8.5)
save_age_plot(
  p_BHB_all_temperature_seasonal_context,
  regional_age_figure_dir,
  "BHB_all_temperature_proxies_seasonal_context_age",
  8.2, 8.5,
  presentation_width = 12
)
save_age_plot(p_insolation_47N_horizontal, regional_age_figure_dir,
              "BHB_ZB20a_summer_insolation_47N_age_horizontal", 7.5, 5.5,
              presentation_width = 6)
save_age_plot(p_BHB_LMA_MAT_age_horizontal, regional_age_figure_dir,
              "BHB_LMA_MAT_model_age_horizontal", 7.5, 5.5,
              presentation_width = 6)
save_age_plot(p_BHB_D47_temperature_age_horizontal, regional_age_figure_dir,
              "BHB_D47_formation_temperature_model_age_horizontal", 7.5, 5.5,
              presentation_width = 6)
save_age_plot(
  p_BHB_D47_seasonal_context_horizontal, regional_age_figure_dir,
  "BHB_D47_seasonal_synthesis_context_age_horizontal", 8.5, 6,
  presentation_width = 6
)
save_age_plot(
  p_BHB_D47_temperature_age_PETM_horizontal, regional_age_figure_dir,
  "BHB_D47_formation_temperature_PETM_56.5-55.5Ma", 7.5, 5.5,
  presentation_width = 6
)
save_age_plot(
  p_BHB_D47_seasonal_context_PETM_horizontal, regional_age_figure_dir,
  "BHB_D47_seasonal_synthesis_PETM_56.5-55.5Ma", 8.5, 6,
  presentation_width = 6
)
save_age_plot(p_T47_petm_phase_context, regional_age_figure_dir,
              "BHB_D47_pre_during_post_PETM_seasonal_context", 9.5, 5.8)
save_age_plot(p_T47_petm_phase_with_evidence, regional_age_figure_dir,
              "BHB_D47_PETM_phases_with_seasonal_evidence", 16, 7.5)
save_age_plot(p_BHB_temperature_insolation, regional_age_figure_dir,
              "BHB_temperature_models_ZB20a_insolation_47N_age", 15.5, 8.5)
save_figure_variants(
  plot = p_global_Barnet_BHB_d13C,
  presentation_plot = p_global_Barnet_BHB_d13C_slide_12x6,
  base_dir = regional_age_figure_dir,
  stem = "Barnet2019_marine_BHB_d13C_age",
  manuscript_width = 11,
  manuscript_height = 9.5,
  presentation_width = 12,
  presentation_height = 6
)
save_figure_variants(
  plot = p_Barnet_stable_isotopes_full_horizontal_6x3,
  presentation_plot = p_Barnet_stable_isotopes_full_horizontal_6x3,
  base_dir = regional_age_figure_dir,
  stem = "Barnet2019_Atlantic_Pacific_d13C_d18O_horizontal_6x3",
  manuscript_width = 6,
  manuscript_height = 4,
  presentation_width = 7,
  presentation_height = 4
)

# Supporting table used in the combined figure. This makes explicit which
# observations, age ranges, proxy meanings, and uncertainty values were drawn.
write_csv(
  BHB_temperature_proxy_data,
  here("data", "processed", "BHB_temperature_proxy_plot_data.csv")
)
write_csv(
  phase_summary,
  here("data", "processed", "BHB_D47_PETM_phase_distribution_summary.csv")
)
write_csv(
  BHB_d13C_published_age,
  here("data", "processed", "BHB_d13C_published_plot_data.csv")
)
