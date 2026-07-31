# 09_plot_CFB_strat_domain.R
# Purpose: Generate reproducible, publication-ready stratigraphic-domain
#          figures from the authoritative CFB soil-carbonate products.

#-- 1.) Setup ---------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)
source(here("scripts", "helpers", "CFB_chronostrat_panels.R"))
source(here("scripts", "helpers", "save_figure_variants.R"))

strat_figure_dir <- here("figures", "strat_domain", "CFB")
dir.create(strat_figure_dir, recursive = TRUE, showWarnings = FALSE)

# Centralized display settings make interval revisions explicit.
petm_strat_min_m <- 1500
petm_strat_max_m <- 1540

source_colors <- c(
  "U-M" = "#000000",
  "CU" = "#0072B2",
  "Caltech" = "#D55E00",
  "Koch" = "#009E73",
  "Bowen" = "#CC79A7"
)

source_shapes <- c(
  "U-M" = 21, "CU" = 22, "Caltech" = 24, "Koch" = 23, "Bowen" = 25
)

t47_status_colors <- c(
  "This study" = "#B2182B",
  "Published data" = "grey25"
)
t47_status_shapes <- c("This study" = 21, "Published data" = 22)

add_petm_strat <- function(alpha = 0.14) {
  annotate(
    "rect", xmin = -Inf, xmax = Inf,
    ymin = petm_strat_min_m, ymax = petm_strat_max_m,
    fill = "#E41A1C", alpha = alpha
  )
}

theme_strat <- theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 11),
    plot.margin = margin(5, 5, 5, 5)
  )

save_plot_pair <- function(
    plot, stem, width, height, presentation_width = NULL
) {
  save_figure_variants(
    plot = plot,
    base_dir = strat_figure_dir,
    stem = stem,
    manuscript_width = width,
    manuscript_height = height,
    presentation_width = presentation_width
  )
}

#-- 2.) Load Authoritative Pipeline Products -------------------------------
CFB_isotopes <- read_csv(
  here("data", "processed", "CFB_soilcarb_isotope_summary.csv"),
  show_col_types = FALSE
)
CFB_temperature <- read_csv(
  here("data", "processed", "CFB_temperature_model.csv"),
  show_col_types = FALSE
)
CFB_temperature_obs <- read_csv(
  here("data", "processed", "CFB_temperature_observations.csv"),
  show_col_types = FALSE
)
CFB_soilwater <- read_csv(
  here("data", "processed", "CFB_soilwater_reconstruction_summary.csv"),
  show_col_types = FALSE
)

CFB_temperature_obs <- CFB_temperature_obs %>%
  filter(used_in_primary_temperature_model) %>%
  mutate(
    source = recode(source, IPL = "U-M", Snell = "Caltech"),
    study_status = if_else(source == "U-M", "This study", "Published data"),
    study_status = factor(
      study_status, levels = c("This study", "Published data")
    )
  )

if (any(CFB_isotopes$section_id != "CFB")) {
  stop("Stratigraphic plotting input contains a non-CFB section.")
}

strat_limits <- range(CFB_isotopes$strat_height_m, na.rm = TRUE)
strat_limits <- c(
  floor(strat_limits[1] / 100) * 100,
  ceiling(strat_limits[2] / 100) * 100
)

shared_strat_scale <- scale_y_continuous(
  limits = strat_limits,
  breaks = seq(strat_limits[1], strat_limits[2], by = 200),
  minor_breaks = seq(strat_limits[1], strat_limits[2], by = 100),
  expand = expansion(mult = c(0.01, 0.02))
)

p_CFB_chronostrat_full <- build_CFB_chronostrat_strat_panel(strat_limits)

#-- 3.) Reshape Source-Specific Carbon-Isotope Observations ----------------
d13C_observations <- bind_rows(
  CFB_isotopes %>% transmute(
    MLA_horizon_id, strat_height_m, source = "U-M",
    value = IPL_NuDog_d13Ccarb_VPDB, se = NA_real_
  ),
  CFB_isotopes %>% transmute(
    MLA_horizon_id, strat_height_m, source = "CU",
    value = CU_mean_d13Ccarb_vpdb, se = NA_real_
  ),
  CFB_isotopes %>% transmute(
    MLA_horizon_id, strat_height_m, source = "Caltech",
    value = Snell_mean_d13Ccarb_vpdb, se = NA_real_
  ),
  CFB_isotopes %>% transmute(
    MLA_horizon_id, strat_height_m, source = "Koch",
    value = Koch_mean_d13Ccarb_vpdb, se = Koch_se_d13Ccarb_vpdb
  ),
  CFB_isotopes %>% transmute(
    MLA_horizon_id, strat_height_m, source = "Bowen",
    value = Bowen_mean_d13Ccarb_vpdb, se = Bowen_se_d13Ccarb_vpdb
  )
) %>%
  filter(!is.na(value), !is.na(strat_height_m)) %>%
  mutate(
    source = factor(source, levels = names(source_colors)),
    study_status = if_else(source == "U-M", "U-M / this study", "Published"),
    study_status = factor(
      study_status, levels = c("U-M / this study", "Published")
    )
  )

#-- 4.) Temperature in Stratigraphic Space ---------------------------------
p_temperature_strat <- ggplot() +
  add_petm_strat() +
  geom_ribbon(
    data = CFB_temperature %>% filter(!is.na(T_model_C)),
    aes(
      xmin = T_model_lower95_C, xmax = T_model_upper95_C,
      y = strat_height_m
    ),
    fill = "#B2182B", alpha = 0.18
  ) +
  geom_path(
    data = CFB_temperature %>% filter(!is.na(T_model_C)),
    aes(T_model_C, strat_height_m),
    color = "#B2182B", linewidth = 1
  ) +
  geom_errorbarh(
    data = CFB_temperature_obs,
    aes(
      xmin = T_C - T_se_C, xmax = T_C + T_se_C,
      y = strat_height_m, color = study_status
    ),
    height = 0, linewidth = 0.35, alpha = 0.55
  ) +
  geom_point(
    data = CFB_temperature_obs,
    aes(
      T_C, strat_height_m, color = study_status,
      fill = p_altered_preservation, shape = study_status
    ),
    size = 2.4, stroke = 0.8
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
  shared_strat_scale +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = "CFB stratigraphic height (m)",
    color = NULL, shape = "Dataset",
    title = "A  Temperature"
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
  theme_strat +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.justification = "left"
  )

#-- 5.) Carbon Isotopes in Stratigraphic Space -----------------------------
p_d13C_strat <- ggplot(d13C_observations) +
  add_petm_strat() +
  geom_errorbarh(
    aes(
      xmin = value - se, xmax = value + se,
      y = strat_height_m, color = source
    ),
    height = 0, linewidth = 0.3, alpha = 0.35, na.rm = TRUE
  ) +
  geom_point(
    aes(
      value, strat_height_m, color = study_status,
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
  shared_strat_scale +
  labs(
    x = expression(delta^13 * C[carbonate] ~ "(per mil VPDB)"),
    y = "CFB stratigraphic height (m)",
    color = NULL, fill = "Dataset", shape = NULL,
    title = expression("B  " * delta^13 * C[carbonate])
  ) +
  guides(
    color = "none",
    shape = "none",
    fill = guide_legend(override.aes = list(shape = 21))
  ) +
  theme_strat

#-- 6.) Reconstructed Soil Water in Stratigraphic Space -------------------
d18Owater_obs <- CFB_soilwater %>%
  filter(!is.na(d18Ow_mean_vsmow), !is.na(strat_height_m))

p_d18Owater_strat <- ggplot(d18Owater_obs) +
  add_petm_strat() +
  geom_errorbarh(
    aes(
      xmin = d18Ow_lower95_vsmow, xmax = d18Ow_upper95_vsmow,
      y = strat_height_m
    ),
    height = 0, linewidth = 0.4, alpha = 0.45, color = "#2166AC"
  ) +
  geom_point(
    aes(d18Ow_mean_vsmow, strat_height_m),
    shape = 21, fill = "white", color = "#2166AC", size = 2
  ) +
  shared_strat_scale +
  labs(
    x = expression(delta^18 * O[water] ~ "(per mil VSMOW)"),
    y = "CFB stratigraphic height (m)",
    title = expression("C  Reconstructed " * delta^18 * O[water])
  ) +
  theme_strat +
  theme(legend.position = "none")

D17Owater_obs <- CFB_soilwater %>%
  filter(!is.na(D17Orsw_mean_permeg), !is.na(strat_height_m))

p_D17Owater_strat <- ggplot(D17Owater_obs) +
  add_petm_strat() +
  geom_errorbarh(
    aes(
      xmin = D17Orsw_lower95_permeg, xmax = D17Orsw_upper95_permeg,
      y = strat_height_m
    ),
    height = 0, linewidth = 0.4, alpha = 0.45, color = "#762A83"
  ) +
  geom_point(
    aes(D17Orsw_mean_permeg, strat_height_m),
    shape = 21, fill = "white", color = "#762A83", size = 2
  ) +
  shared_strat_scale +
  labs(
    x = expression(Delta*"'"^17 * O[water] ~ "(per meg)"),
    y = "CFB stratigraphic height (m)",
    title = expression("D  Reconstructed " * Delta*"'"^17 * O[water])
  ) +
  theme_strat +
  theme(legend.position = "none")

#-- 7.) Assemble Full-Section and PETM-Focused Panels ----------------------
remove_repeated_y <- function(plot) {
  plot + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
}

p_CFB_strat_full <-
  p_CFB_chronostrat_full +
  (p_temperature_strat + theme(legend.position = "none")) +
  remove_repeated_y(p_d13C_strat) +
  remove_repeated_y(p_d18Owater_strat) +
  remove_repeated_y(p_D17Owater_strat) +
  plot_layout(widths = c(0.88, 1.12, 1, 1, 1), guides = "collect") &
  theme(legend.position = "top")

petm_zoom_limits <- c(1450, 1600)
p_CFB_chronostrat_petm <- build_CFB_chronostrat_strat_panel(
  petm_zoom_limits,
  breaks = seq(1450, 1600, by = 25)
)
p_CFB_strat_petm <-
  p_CFB_chronostrat_petm +
  (p_temperature_strat + coord_cartesian(ylim = petm_zoom_limits) +
     theme(legend.position = "none")) +
  remove_repeated_y(p_d13C_strat + coord_cartesian(ylim = petm_zoom_limits)) +
  remove_repeated_y(p_d18Owater_strat + coord_cartesian(ylim = petm_zoom_limits)) +
  remove_repeated_y(p_D17Owater_strat + coord_cartesian(ylim = petm_zoom_limits)) +
  plot_layout(widths = c(0.88, 1.12, 1, 1, 1), guides = "collect") &
  theme(legend.position = "top")

#-- 8.) Export Figures ------------------------------------------------------
with_chronostrat <- function(plot, framework = p_CFB_chronostrat_full) {
  framework + plot + plot_layout(widths = c(0.88, 2.1))
}

save_plot_pair(
  with_chronostrat(p_temperature_strat),
  "CFB_temperature_strat", 7.4, 7.5, presentation_width = 6
)
save_plot_pair(
  with_chronostrat(p_d13C_strat),
  "CFB_d13Ccarb_strat", 7.4, 7.5, presentation_width = 6
)
save_plot_pair(
  with_chronostrat(p_d18Owater_strat),
  "CFB_d18Owater_strat", 7.4, 7.5, presentation_width = 6
)
save_plot_pair(
  with_chronostrat(p_D17Owater_strat),
  "CFB_D17Owater_strat", 7.4, 7.5, presentation_width = 6
)
save_plot_pair(
  p_CFB_strat_full, "CFB_multiproxy_strat_full", 16.5, 8,
  presentation_width = 12
)
save_plot_pair(
  p_CFB_strat_petm, "CFB_multiproxy_strat_PETM_zoom", 16.5, 6,
  presentation_width = 12
)
