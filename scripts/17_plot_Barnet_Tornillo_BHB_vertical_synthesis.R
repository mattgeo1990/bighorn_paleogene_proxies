# 17_plot_Barnet_Tornillo_BHB_vertical_synthesis.R
# Purpose: Age-aligned marine, Tornillo, and BHB temperature comparison for
#          presentation. The marine temperature axis is calibrated directly
#          from the Barnet et al. (2019) published Site 1262 BWT values.

library(tidyverse)
library(here)
library(patchwork)

source(here("scripts", "helpers", "save_figure_variants.R"))

age_limits <- c(59, 52.5)
petm_old_ma <- 55.93
petm_young_ma <- 55.75
figure_dir <- here("figures", "age_domain", "regional_temperature_synthesis")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

age_scale_left <- scale_y_reverse(
  limits = age_limits,
  breaks = seq(59, 52.5, by = -0.5),
  minor_breaks = seq(59, 52.5, by = -0.1),
  expand = expansion(mult = c(0.004, 0.006))
)
age_scale_blank <- scale_y_reverse(
  limits = age_limits,
  breaks = seq(59, 52.5, by = -0.5),
  minor_breaks = seq(59, 52.5, by = -0.1),
  labels = NULL,
  name = NULL,
  expand = expansion(mult = c(0.004, 0.006))
)

add_petm <- function(alpha = 0.075) {
  annotate(
    "rect", xmin = -Inf, xmax = Inf,
    ymin = petm_young_ma, ymax = petm_old_ma,
    fill = "#D73027", alpha = alpha
  )
}

theme_panel <- theme_classic(base_size = 18) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18, color = "black"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.margin = margin(2, 5, 2, 5),
    legend.position = "none"
  )

#-- Marine records and Barnet temperature conversion ----------------------
atlantic <- read_csv(
  here("data", "processed", "BarnetEtAl2019_ODP1262_benthic_isotopes.csv"),
  show_col_types = FALSE
) %>%
  filter(between(Age_Ma, min(age_limits), max(age_limits))) %>%
  transmute(
    Age_Ma, basin = "Atlantic",
    d18O_benthic_vpdb = d18O_corrected_vpdb,
    bottom_water_temperature_C
  )

pacific <- read_csv(
  here("data", "processed", "WesterholdEtAl2018_ODP1209_benthic_isotopes.csv"),
  show_col_types = FALSE
) %>%
  filter(between(Age_Ma, min(age_limits), max(age_limits))) %>%
  transmute(
    Age_Ma, basin = "Pacific", d18O_benthic_vpdb,
    bottom_water_temperature_C = NA_real_
  )

marine <- bind_rows(atlantic, pacific) %>%
  mutate(basin = factor(basin, levels = c("Atlantic", "Pacific"))) %>%
  arrange(basin, Age_Ma)

barnet_conversion <- lm(
  bottom_water_temperature_C ~ d18O_benthic_vpdb,
  data = atlantic %>%
    filter(is.finite(bottom_water_temperature_C), is.finite(d18O_benthic_vpdb))
)
barnet_intercept <- unname(coef(barnet_conversion)[1])
barnet_slope <- unname(coef(barnet_conversion)[2])

marine_colors <- c("Atlantic" = "#29A9E0", "Pacific" = "#173F90")
marine_linetypes <- c("Atlantic" = "solid", "Pacific" = "22")

#-- Tornillo and BHB observations -----------------------------------------
tornillo <- read_csv(
  here("data", "processed", "Tornillo_temperature_model_observations.csv"),
  show_col_types = FALSE
) %>%
  filter(between(Age_Ma, min(age_limits), max(age_limits)))

tornillo_model <- read_csv(
  here("data", "processed", "Tornillo_temperature_model.csv"),
  show_col_types = FALSE
) %>%
  filter(between(Age_Ma, min(age_limits), max(age_limits)))

bhb_d47 <- read_csv(
  here("data", "processed", "BHB_D47_temperature_observations.csv"),
  show_col_types = FALSE
) %>%
  filter(
    used_in_temperature_model,
    between(Age_Ma, min(age_limits), max(age_limits))
  )

bhb_model <- read_csv(
  here("data", "processed", "BHB_D47_temperature_model.csv"),
  show_col_types = FALSE
) %>%
  filter(between(Age_Ma, min(age_limits), max(age_limits)))

wing_lma <- read_csv(
  here("data", "processed", "WingEtAl2000_LMA_MAT_processed.csv"),
  show_col_types = FALSE
) %>%
  filter(
    is.finite(published_age_model_2_Ma), is.finite(MAT_C),
    between(published_age_model_2_Ma, min(age_limits), max(age_limits))
  ) %>%
  transmute(
    Age_Ma = published_age_model_2_Ma,
    age_old_ma = Age_Ma + duration_Myr / 2,
    age_young_ma = Age_Ma - duration_Myr / 2,
    temperature_C = MAT_C, temperature_se_C = MAT_error_C
  )

phosphate <- read_csv(
  here("data", "processed", "FrickeWing2004_BHB_MAAT_processed.csv"),
  show_col_types = FALSE
) %>%
  filter(
    is.finite(Age_Ma), is.finite(temperature_C),
    between(Age_Ma, min(age_limits), max(age_limits))
  ) %>%
  transmute(
    Age_Ma, age_old_ma = age_older_ma, age_young_ma = age_younger_ma,
    temperature_C
  )

seasonal <- read_csv(
  here("data", "processed", "BHB_seasonal_temperature_integrated_summary.csv"),
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
    fill = seasonal_colors[[metric_name]], alpha = alpha_value,
    inherit.aes = FALSE
  )
}

#-- Shared age/epoch strip -------------------------------------------------
epochs <- tribble(
  ~epoch, ~older_ma, ~younger_ma, ~fill,
  "Paleocene", 59.0, 56.0, "#F6E8A6",
  "Eocene", 56.0, 52.5, "#F7B267"
)

p_age_epoch <- ggplot(epochs) +
  geom_rect(
    aes(xmin = 0, xmax = 1, ymin = younger_ma, ymax = older_ma, fill = fill),
    color = "black", linewidth = 0.35
  ) +
  geom_text(
    aes(x = 0.5, y = (older_ma + younger_ma) / 2, label = epoch),
    angle = 90, size = 18 / ggplot2::.pt
  ) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  age_scale_left +
  labs(x = NULL, y = "Age (Ma)") +
  theme_classic(base_size = 18) +
  theme(
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(2, 1, 2, 5)
  )

#-- Data panels ------------------------------------------------------------
p_marine <- ggplot(
  marine,
  aes(d18O_benthic_vpdb, Age_Ma, color = basin, linetype = basin)
) +
  add_petm() +
  geom_path(linewidth = 0.48, alpha = 0.92, na.rm = TRUE) +
  scale_color_manual(values = marine_colors) +
  scale_linetype_manual(values = marine_linetypes) +
  scale_x_reverse(
    limits = c(1.25, -1.55), breaks = c(1, 0, -1),
    name = expression(delta^18 * O[benthic] ~ "(per mil VPDB)"),
    sec.axis = sec_axis(
      transform = ~ barnet_intercept + barnet_slope * ., 
      name = expression("Ice-free bottom-water T (" * degree * "C)"),
      breaks = seq(5, 20, by = 5)
    )
  ) +
  age_scale_blank +
  labs(title = "Marine reference") +
  theme_panel +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x.top = element_text(size = 18),
    axis.text.x.top = element_text(size = 18)
  )

p_tornillo <- ggplot() +
  add_petm() +
  geom_ribbon(
    data = tornillo_model,
    aes(xmin = temperature_lower95_C, xmax = temperature_upper95_C, y = Age_Ma),
    fill = "#B2182B", alpha = 0.14
  ) +
  geom_ribbon(
    data = tornillo_model,
    aes(xmin = temperature_lower80_C, xmax = temperature_upper80_C, y = Age_Ma),
    fill = "#B2182B", alpha = 0.23
  ) +
  geom_path(
    data = tornillo_model, aes(temperature_model_C, Age_Ma),
    color = "#B2182B", linewidth = 0.9
  ) +
  geom_errorbarh(
    data = tornillo,
    aes(xmin = T47_C - T47_se_C, xmax = T47_C + T47_se_C, y = Age_Ma),
    height = 0, color = "grey35", linewidth = 0.35, alpha = 0.55
  ) +
  geom_errorbar(
    data = tornillo,
    aes(x = T47_C, ymin = Age_Ma - 0.2, ymax = Age_Ma + 0.2),
    width = 0, color = "grey35", linewidth = 0.35, alpha = 0.55
  ) +
  geom_point(
    data = tornillo, aes(T47_C, Age_Ma),
    shape = 21, fill = "white", color = "black", size = 2.7, stroke = 0.8
  ) +
  scale_x_continuous(limits = c(10, 50), breaks = seq(10, 50, by = 10)) +
  age_scale_blank +
  labs(title = "Tornillo Basin", x = expression(T[47] * " soil carbonate (" * degree * "C)")) +
  theme_panel +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

alteration_colors <- c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B")

p_bhb <- ggplot() +
  seasonal_rect("CMMT", "lower95_c", "upper95_c", 0.045) +
  seasonal_rect("CMMT", "lower80_c", "upper80_c", 0.070) +
  seasonal_rect("CMMT", "lower50_c", "upper50_c", 0.105) +
  seasonal_rect("MAAT", "lower95_c", "upper95_c", 0.045) +
  seasonal_rect("MAAT", "lower80_c", "upper80_c", 0.070) +
  seasonal_rect("MAAT", "lower50_c", "upper50_c", 0.105) +
  seasonal_rect("WMMT", "lower95_c", "upper95_c", 0.045) +
  seasonal_rect("WMMT", "lower80_c", "upper80_c", 0.070) +
  seasonal_rect("WMMT", "lower50_c", "upper50_c", 0.105) +
  geom_vline(
    data = seasonal, aes(xintercept = mean_c, color = metric),
    linewidth = 0.55, alpha = 0.70, show.legend = FALSE
  ) +
  add_petm() +
  geom_ribbon(
    data = bhb_model,
    aes(xmin = temperature_lower95_C, xmax = temperature_upper95_C, y = Age_Ma),
    fill = "#B2182B", alpha = 0.14
  ) +
  geom_ribbon(
    data = bhb_model,
    aes(xmin = temperature_lower80_C, xmax = temperature_upper80_C, y = Age_Ma),
    fill = "#B2182B", alpha = 0.23
  ) +
  geom_path(
    data = bhb_model, aes(temperature_model_C, Age_Ma),
    color = "#B2182B", linewidth = 0.9
  ) +
  geom_errorbarh(
    data = bhb_d47,
    aes(
      xmin = temperature_C - temperature_se_C,
      xmax = temperature_C + temperature_se_C, y = Age_Ma
    ),
    height = 0, color = "grey35", linewidth = 0.34, alpha = 0.55
  ) +
  geom_errorbar(
    data = bhb_d47,
    aes(x = temperature_C, ymin = Age_Ma - 0.2, ymax = Age_Ma + 0.2),
    width = 0, color = "grey35", linewidth = 0.34, alpha = 0.55
  ) +
  geom_point(
    data = bhb_d47,
    aes(temperature_C, Age_Ma, fill = p_altered_preservation),
    shape = 21, color = "black", size = 2.8, stroke = 0.75
  ) +
  geom_errorbarh(
    data = wing_lma,
    aes(
      xmin = temperature_C - temperature_se_C,
      xmax = temperature_C + temperature_se_C, y = Age_Ma
    ),
    height = 0, color = "#1B7837", linewidth = 0.35, alpha = 0.6
  ) +
  geom_errorbar(
    data = wing_lma,
    aes(x = temperature_C, ymin = age_young_ma, ymax = age_old_ma),
    width = 0, color = "#1B7837", linewidth = 0.4, alpha = 0.55
  ) +
  geom_point(
    data = wing_lma, aes(temperature_C, Age_Ma),
    shape = 22, fill = "white", color = "#1B7837", size = 2.8, stroke = 0.85
  ) +
  geom_errorbar(
    data = phosphate,
    aes(x = temperature_C, ymin = age_young_ma, ymax = age_old_ma),
    width = 0, color = "#762A83", linewidth = 0.4, alpha = 0.55
  ) +
  geom_point(
    data = phosphate, aes(temperature_C, Age_Ma),
    shape = 23, fill = "white", color = "#762A83", size = 2.8, stroke = 0.85
  ) +
  scale_fill_gradientn(
    colors = alteration_colors, limits = c(0, 1),
    breaks = c(0, 0.5, 1), labels = scales::label_percent(accuracy = 1),
    name = "BHB T47 alteration probability"
  ) +
  scale_color_manual(values = seasonal_colors, guide = "none") +
  scale_x_continuous(breaks = seq(0, 60, by = 10)) +
  age_scale_blank +
  coord_cartesian(xlim = c(0, 60)) +
  labs(title = "Bighorn Basin synthesis", x = expression("Temperature (" * degree * "C)")) +
  guides(
    fill = guide_colorbar(
      direction = "horizontal", barwidth = unit(5.0, "cm"),
      barheight = unit(0.35, "cm"), title.position = "top"
    )
  ) +
  theme_panel +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    legend.position = "none"
  )

#-- Legends ---------------------------------------------------------------
legend_proxy <- tibble(
  x = 0.03, y = c(0.62, 0.36, 0.10),
  label = c(
    "T47 soil carbonate", "Vertebrate phosphate d18O",
    "Leaf-Margin Analysis"
  ),
  shape = factor(c("T47", "Phosphate", "LMA"), levels = c("T47", "Phosphate", "LMA"))
)

p_legend <- ggplot() +
  geom_segment(
    data = tibble(x = c(0.03, 0.18), xend = c(0.09, 0.24), y = 0.92),
    aes(x, y, xend = xend, yend = y),
    color = c(marine_colors[["Atlantic"]], marine_colors[["Pacific"]]),
    linetype = c("solid", "22"), linewidth = 0.8
  ) +
  annotate("text", x = c(0.10, 0.25), y = 0.92,
           label = c("Atlantic", "Pacific"), hjust = 0,
           size = 18 / ggplot2::.pt) +
  geom_point(
    data = legend_proxy, aes(x, y, shape = shape),
    size = 3.2, fill = "white", color = "black", stroke = 0.8
  ) +
  geom_text(
    data = legend_proxy, aes(x + 0.035, y, label = label),
    hjust = 0, size = 18 / ggplot2::.pt
  ) +
  annotate(
    "text", x = 0.35, y = 0.36, hjust = 0,
    label = "Tornillo T47 = open circles\nBHB T47 fill = alteration probability",
    size = 18 / ggplot2::.pt, lineheight = 1.05
  ) +
  geom_tile(
    data = tibble(
      x = seq(0.68, 0.96, length.out = 101), y = 0.21,
      probability = seq(0, 1, length.out = 101)
    ),
    aes(x, y, fill = probability), width = 0.003, height = 0.16
  ) +
  annotate(
    "text", x = 0.68, y = 0.43,
    label = "BHB T47 alteration probability", hjust = 0,
    size = 18 / ggplot2::.pt
  ) +
  annotate(
    "text", x = c(0.68, 0.82, 0.96), y = 0.03,
    label = c("0%", "50%", "100%"), size = 18 / ggplot2::.pt
  ) +
  annotate(
    "text", x = 0.68, y = 0.92,
    label = "BHB seasonal synthesis", hjust = 0,
    size = 18 / ggplot2::.pt
  ) +
  annotate(
    "rect",
    xmin = c(0.68, 0.78, 0.88), xmax = c(0.72, 0.82, 0.92),
    ymin = 0.66, ymax = 0.80,
    fill = unname(seasonal_colors), alpha = 0.30
  ) +
  annotate(
    "text", x = c(0.73, 0.83, 0.93), y = 0.73,
    label = c("CMMT", "MAAT", "WMMT"), hjust = 0,
    size = 18 / ggplot2::.pt
  ) +
  scale_shape_manual(values = c("T47" = 21, "Phosphate" = 23, "LMA" = 22)) +
  scale_fill_gradientn(colors = alteration_colors, limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() + theme(legend.position = "none")

main_panels <- p_age_epoch + p_marine + p_tornillo + p_bhb +
  plot_layout(widths = c(0.55, 2.15, 3.0, 4.3))

final_plot <- main_panels / p_legend +
  plot_layout(
    heights = c(4.35, 1.70)
  ) +
  plot_annotation(
    title = "Marine and continental temperature records, 59-52.5 Ma",
    caption = str_wrap(paste0(
      "Marine curves: Barnet et al. (2019) Atlantic Site 1262 and its Pacific Site 1209 comparison. ",
      "Top marine axis uses the Barnet published Site 1262 BWT-d18O relation (T = ",
      round(barnet_intercept, 2), " ", if_else(barnet_slope < 0, "-", "+"), " ",
      round(abs(barnet_slope), 2), " x d18O). Terrestrial T47 age bars show an ",
      "assumed +/-0.2 Myr plotting uncertainty; formal age uncertainty is not propagated. ",
      "Pink band marks the PETM."
    ), width = 145),
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.caption = element_text(size = 18, hjust = 0),
      plot.margin = margin(4, 7, 4, 7)
    )
  )

save_figure_variants(
  plot = final_plot,
  presentation_plot = final_plot,
  base_dir = figure_dir,
  stem = "Barnet_Tornillo_BHB_vertical_temperature_synthesis_59_52p5Ma",
  manuscript_width = 12,
  manuscript_height = 8,
  presentation_width = 12,
  presentation_height = 8
)

message(
  "Saved 12 x 8 synthesis with 18-point annotation space; Barnet conversion slope = ",
  round(barnet_slope, 4), " C per mil; R2 = ",
  round(summary(barnet_conversion)$r.squared, 5), "."
)
