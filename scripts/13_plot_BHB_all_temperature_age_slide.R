# 13_plot_BHB_all_temperature_age_slide.R
# Purpose: Presentation-ready age-domain synthesis of all Bighorn Basin
#          temperature proxy observations, including the new IPL Delta47 data.
#
# Design:
#   - Preserves the seasonal-context visual language of the published-only
#     stratigraphic figure.
#   - Uses age rather than CFB height.
#   - Shows only the mammal-biozone column.
#   - IPL observations are filled by the d18O-only alteration index.
#   - Published proxy observations remain open symbols.
#   - Does not overwrite the published-only figure.

library(tidyverse)
library(here)
library(patchwork)

source(here("scripts", "helpers", "save_figure_variants.R"))

figure_dir <- here("figures", "age_domain", "regional_comparison")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

age_limits <- c(59.0, 52.5)
petm_old_ma <- 55.93
petm_young_ma <- 55.75

age_scale <- scale_y_reverse(
  limits = age_limits,
  breaks = seq(59, 52.5, by = -0.5),
  minor_breaks = seq(59, 52.5, by = -0.1),
  expand = expansion(mult = c(0.008, 0.012))
)

theme_slide <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 10.5, color = "black"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9.5),
    legend.key.height = unit(0.38, "cm"),
    legend.spacing.x = unit(0.12, "cm"),
    plot.margin = margin(5, 6, 5, 5)
  )

#-- Load and age-calibrate CFB Delta47 observations ------------------------
cfb_temperature <- read_csv(
  here("data", "processed", "CFB_temperature_observations.csv"),
  show_col_types = FALSE
)

cfb_age_lookup <- read_csv(
  here(
    "data", "processed",
    "CFB_soilcarb_with_temperature_age_calibrated.csv"
  ),
  show_col_types = FALSE
) %>%
  distinct(section_id, MLA_horizon_id, strat_height_m, Age_Ma) %>%
  filter(is.finite(Age_Ma))

cfb_temperature_age <- cfb_temperature %>%
  left_join(
    cfb_age_lookup,
    by = c("section_id", "MLA_horizon_id", "strat_height_m")
  ) %>%
  filter(
    is.finite(Age_Ma), is.finite(T_C),
    used_in_primary_temperature_model
  ) %>%
  mutate(
    dataset = if_else(source == "U-M", "IPL (this study)", "Published D47"),
    record = if_else(source == "U-M", "IPL D47", "Published D47")
  )

# Published Snell Delta47 data from McCullough Peaks extend the record beyond
# CFB. CFB Snell values are already represented in cfb_temperature_age.
mcp_snell <- read_csv(
  here(
    "data", "processed",
    "SnellEtAl2013_soilcarb_summary_age_calibrated.csv"
  ),
  show_col_types = FALSE
) %>%
  filter(
    section_id == "MCP",
    is.finite(Age_Ma),
    is.finite(Snell_mean_T47_C)
  ) %>%
  transmute(
    Age_Ma,
    T_C = Snell_mean_T47_C,
    T_se_C = Snell_se_T47_C,
    dataset = "Published D47",
    record = "Published D47"
  )

wing_lma <- read_csv(
  here("data", "processed", "WingEtAl2000_LMA_MAT_processed.csv"),
  show_col_types = FALSE
) %>%
  filter(is.finite(published_age_model_2_Ma), is.finite(MAT_C)) %>%
  transmute(
    Age_Ma = published_age_model_2_Ma,
    age_old_ma = Age_Ma + duration_Myr / 2,
    age_young_ma = Age_Ma - duration_Myr / 2,
    T_C = MAT_C,
    T_se_C = MAT_error_C,
    record = "Wing LMA"
  )

fricke_wing <- read_csv(
  here("data", "processed", "FrickeWing2004_BHB_MAAT_processed.csv"),
  show_col_types = FALSE
) %>%
  filter(is.finite(Age_Ma), is.finite(temperature_C)) %>%
  transmute(
    Age_Ma,
    age_old_ma = age_older_ma,
    age_young_ma = age_younger_ma,
    T_C = temperature_C,
    T_se_C = NA_real_,
    record = "Fricke-Wing MAAT"
  )

published_d47 <- bind_rows(
  cfb_temperature_age %>%
    filter(dataset == "Published D47") %>%
    select(Age_Ma, T_C, T_se_C, dataset, record),
  mcp_snell
) %>%
  filter(between(Age_Ma, min(age_limits), max(age_limits)))

ipl_d47 <- cfb_temperature_age %>%
  filter(dataset == "IPL (this study)") %>%
  filter(between(Age_Ma, min(age_limits), max(age_limits)))

#-- Seasonal temperature context ------------------------------------------
# Prefer the pipeline product when available. The embedded values are the
# versioned v1 summary currently tracked by the project and make this figure
# reproducible even when generated products have been cleaned.
seasonal_path <- here(
  "data", "processed", "BHB_seasonal_temperature_integrated_summary.csv"
)

if (file.exists(seasonal_path)) {
  seasonal <- read_csv(seasonal_path, show_col_types = FALSE)
} else {
  seasonal <- tribble(
    ~metric, ~mean_c, ~lower50_c, ~upper50_c, ~lower80_c, ~upper80_c,
    ~lower95_c, ~upper95_c,
    "CMMT", 5.97534, 4.11351, 8.05527, 1.90598, 9.83079, -1.15733, 11.92443,
    "MAAT", 19.40890, 17.42765, 21.42392, 15.46471, 23.09650, 13.50832, 25.01724,
    "WMMT", 34.73502, 30.83933, 38.62994, 28.34929, 41.96826, 25.81680, 44.87044
  )
}

seasonal <- seasonal %>%
  filter(metric %in% c("CMMT", "MAAT", "WMMT")) %>%
  mutate(metric = factor(metric, levels = c("CMMT", "MAAT", "WMMT")))

seasonal_colors <- c(
  "CMMT" = "#2166AC",
  "MAAT" = "#8C6BB1",
  "WMMT" = "#D73027"
)
record_colors <- c(
  "Published D47" = "grey35",
  "Wing LMA" = "#1B7837",
  "Fricke-Wing MAAT" = "#762A83"
)
record_shapes <- c(
  "Published D47" = 21,
  "Wing LMA" = 22,
  "Fricke-Wing MAAT" = 23
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

#-- Convert the CFB mammal zones to age ------------------------------------
meter_age <- cfb_age_lookup %>%
  distinct(strat_height_m, Age_Ma) %>%
  arrange(strat_height_m)

meter_to_age <- function(x) {
  approx(
    x = meter_age$strat_height_m,
    y = meter_age$Age_Ma,
    xout = x,
    rule = 2,
    ties = mean
  )$y
}

#biozones <- read_csv(
#  here("data", "raw", "NorBHB_strat_columns.csv"),
#  show_col_types = FALSE
#) %>%
#  filter(column == "Biozone") %>%
#  transmute(
#    label,
#    age_a = meter_to_age(base_m),
#    age_b = meter_to_age(top_m),
#    ymin = pmax(pmin(age_a, age_b), min(age_limits)),
#    ymax = pmin(pmax(age_a, age_b), max(age_limits)),
#    ymid = (ymin + ymax) / 2
#  ) %>%
#  filter(ymax > ymin)

#p_biozones <- ggplot(biozones) +
# geom_rect(
#   aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax),
#    fill = "#F7F7F7", color = "black", linewidth = 0.35
#  ) +
#  geom_text(
#    data = biozones %>% filter((ymax - ymin) >= 0.075),
#    aes(x = 0.5, y = ymid, label = label),
#    size = 3.2, lineheight = 0.9
#  ) +
#  scale_x_continuous(
#    limits = c(0, 1), breaks = 0.5, labels = "Mammal\nbiozone",
#    position = "top", expand = c(0, 0)
#  ) +
#  age_scale +
#  labs(x = NULL, y = NULL) +
#  coord_cartesian(clip = "off") +
#  theme_classic(base_size = 10) +
#  theme(
#    axis.text.y = element_blank(),
#    axis.ticks.y = element_blank(),
#    axis.line = element_blank(),
#    axis.text.x.top = element_text(size = 9, face = "bold"),
#    axis.ticks.x = element_blank(),
#    plot.margin = margin(5, 1, 5, 2)
# )

#-- Main temperature panel -------------------------------------------------
p_temperature <- ggplot() +
  seasonal_rect("CMMT", "lower95_c", "upper95_c", 0.055) +
  seasonal_rect("CMMT", "lower80_c", "upper80_c", 0.075) +
  seasonal_rect("CMMT", "lower50_c", "upper50_c", 0.11) +
  seasonal_rect("MAAT", "lower95_c", "upper95_c", 0.055) +
  seasonal_rect("MAAT", "lower80_c", "upper80_c", 0.075) +
  seasonal_rect("MAAT", "lower50_c", "upper50_c", 0.11) +
  seasonal_rect("WMMT", "lower95_c", "upper95_c", 0.055) +
  seasonal_rect("WMMT", "lower80_c", "upper80_c", 0.075) +
  seasonal_rect("WMMT", "lower50_c", "upper50_c", 0.11) +
  geom_vline(
    data = seasonal,
    aes(xintercept = mean_c, color = metric),
    linewidth = 0.75, alpha = 0.88, show.legend = FALSE
  ) +
  annotate(
    "rect", xmin = -Inf, xmax = Inf,
    ymin = petm_young_ma, ymax = petm_old_ma,
    fill = "#D73027", alpha = 0.075
  ) +
  geom_errorbarh(
    data = published_d47,
    aes(xmin = T_C - T_se_C, xmax = T_C + T_se_C, y = Age_Ma),
    height = 0, color = record_colors[["Published D47"]],
    linewidth = 0.34, alpha = 0.58, na.rm = TRUE
  ) +
  geom_errorbarh(
    data = wing_lma,
    aes(xmin = T_C - T_se_C, xmax = T_C + T_se_C, y = Age_Ma),
    height = 0, color = record_colors[["Wing LMA"]],
    linewidth = 0.34, alpha = 0.58, na.rm = TRUE
  ) +
  geom_errorbar(
    data = bind_rows(wing_lma, fricke_wing),
    aes(x = T_C, ymin = age_young_ma, ymax = age_old_ma, color = record),
    width = 0, linewidth = 0.45, alpha = 0.55
  ) +
  geom_errorbarh(
    data = ipl_d47,
    aes(xmin = T_C - T_se_C, xmax = T_C + T_se_C, y = Age_Ma),
    height = 0, color = "grey25", linewidth = 0.36, alpha = 0.55
  ) +
  geom_point(
    data = published_d47,
    aes(T_C, Age_Ma, shape = record, color = record),
    fill = "white", size = 2.7, stroke = 0.9
  ) +
  geom_point(
    data = bind_rows(wing_lma, fricke_wing),
    aes(T_C, Age_Ma, shape = record, color = record),
    fill = "white", size = 2.7, stroke = 0.9
  ) +
  geom_point(
    data = ipl_d47,
    aes(T_C, Age_Ma, fill = p_altered_preservation),
    shape = 21, color = "black", size = 3.0, stroke = 0.8
  ) +
  scale_color_manual(
    values = c(seasonal_colors, record_colors),
    breaks = names(record_colors), drop = FALSE
  ) +
  scale_shape_manual(values = record_shapes, drop = FALSE) +
  scale_fill_gradientn(
    colors = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"),
    limits = c(0, 1), breaks = c(0, 0.5, 1),
    labels = scales::label_percent(accuracy = 1),
    name = "IPL d18O trajectory\nP(altered)"
  ) +
  scale_x_continuous(
    limits = c(-2, 62), breaks = seq(0, 60, by = 10),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  age_scale +
  labs(
    x = expression("Temperature (" * degree * "C)"),
    y = "Age (Ma)",
    color = NULL, shape = NULL, fill = "Seasonal synthesis"
  ) +
  guides(
    color = "none",
    shape = guide_legend(
      order = 1, nrow = 1, byrow = TRUE,
      override.aes = list(fill = "white")
    ),
    fill = guide_colorbar(
      order = 2,
      barwidth = unit(2.4, "cm"),
      barheight = unit(0.35, "cm"),
      title.position = "top"
    )
  ) +
  theme_slide +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.margin = margin(0, 0, 1, 0)
  )

p_final <- p_temperature +
  plot_annotation(
    title = "Bighorn Basin temperature proxies",
    caption = paste0(
      "Filled circles: new IPL soil-carbonate \u039447; fill denotes estimated ",
      "d18O-only alteration index. Only model-screened \u039447 observations are shown.\n",
      "Published records: CU and Snell et al. soil-carbonate \u039447,",
      " Wing et al. leaf-margin MAT, and Fricke-Wing phosphate/LMA MAAT. ",
      "Blue/purple/red bands: CMMT/MAAT/WMMT 50/80/95% intervals; ",
      "pink band marks the PETM."
    ),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.caption = element_text(
        size = 7.5, hjust = 0, lineheight = 1.05, color = "grey25"
      ),
      plot.margin = margin(5, 7, 5, 7)
    )
  )

save_figure_variants(
  plot = p_final,
  presentation_plot = p_final,
  base_dir = figure_dir,
  stem = "BHB_all_temperature_proxies_age_IPL_alteration",
  manuscript_width = 9,
  manuscript_height = 6,
  presentation_width = 9,
  presentation_height = 6
)

message(
  "Saved age-domain BHB temperature synthesis: ",
  nrow(ipl_d47), " IPL D47; ",
  nrow(published_d47), " published D47; ",
  nrow(wing_lma), " Wing LMA; ",
  nrow(fricke_wing), " Fricke-Wing MAAT observations."
)
