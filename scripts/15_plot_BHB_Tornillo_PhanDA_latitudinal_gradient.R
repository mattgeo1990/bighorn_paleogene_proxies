# 15_plot_BHB_Tornillo_PhanDA_latitudinal_gradient.R
# Purpose: Compare screened Bighorn Basin and Tornillo clumped-isotope
#          formation temperatures with PhanDA stage-level latitudinal surface-
#          temperature gradients.
#
# Important:
#   - BHB observations are drawn from the production temperature-observation
#     table and ONLY rows passing the current d18Ocarb and temperature screens
#     are plotted (used_in_temperature_model == TRUE).
#   - BHB fill is the current d18Ocarb-based alteration index. This index is
#     retained for visualization even though the plotted observations have
#     already passed the hard screen.
#   - Tornillo points are primary micrite observations from Kelson et al. and
#     are not assigned the BHB-specific alteration index.
#   - PhanDA curves describe stage-level surface air temperature, whereas the
#     clumped-isotope points are soil-carbonate formation temperatures. Their
#     juxtaposition is a proxy/model context comparison, not an equivalence.

#-- 1.) Setup ---------------------------------------------------------------
library(tidyverse)
library(here)
source(here("scripts", "helpers", "save_figure_variants.R"))

BHB_paleolat_deg_n <- 48
Tornillo_paleolat_deg_n <- 33

figure_dir <- here(
  "figures", "regional_comparison", "PhanDA_latitudinal_gradient"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

stage_from_age <- function(age_ma) {
  case_when(
    age_ma > 56.00 & age_ma <= 59.24 ~ "Thanetian",
    age_ma >= 47.83 & age_ma <= 56.00 ~ "Ypresian",
    TRUE ~ NA_character_
  )
}

stage_levels <- c("Thanetian", "Ypresian")

#-- 2.) Load PhanDA Latitudinal Gradients ----------------------------------
PhanDA_ltg <- bind_rows(
  read_csv(
    here("data", "raw", "PhanDA_LTG_Thanetian.csv"),
    show_col_types = FALSE
  ) %>% mutate(stage_group = "Thanetian"),
  read_csv(
    here("data", "raw", "PhanDA_LTG_Ypresian.csv"),
    show_col_types = FALSE
  ) %>% mutate(stage_group = "Ypresian")
) %>%
  filter(Latitude >= 0, Latitude <= 90) %>%
  mutate(stage_group = factor(stage_group, levels = stage_levels))

#-- 3.) Prepare Screened BHB and Tornillo Observations ---------------------
BHB_points <- read_csv(
  here("data", "processed", "BHB_D47_temperature_observations.csv"),
  show_col_types = FALSE
) %>%
  filter(
    used_in_temperature_model,
    !is.na(Age_Ma),
    !is.na(temperature_C),
    !is.na(p_altered_preservation)
  ) %>%
  transmute(
    region = "Bighorn Basin",
    sample_id = MLA_horizon_id,
    dataset,
    study_status,
    Age_Ma,
    stage_group = stage_from_age(Age_Ma),
    paleolatitude_deg_n = BHB_paleolat_deg_n,
    temperature_C,
    temperature_se_C,
    p_altered_preservation,
    passed_BHB_screen = TRUE
  ) %>%
  filter(!is.na(stage_group))

Tornillo_points <- read_csv(
  here(
    "data", "processed",
    "Kelson_Tornillo_micrite_strat_age_filtered.csv"
  ),
  show_col_types = FALSE
) %>%
  filter(!is.na(Age_Ma), !is.na(T47_C)) %>%
  transmute(
    region = "Tornillo Basin",
    sample_id,
    dataset = "Kelson et al. Tornillo micrite",
    study_status = "Published data",
    Age_Ma,
    stage_group = stage_from_age(Age_Ma),
    paleolatitude_deg_n = Tornillo_paleolat_deg_n,
    temperature_C = T47_C,
    temperature_se_C = T47_se_C,
    p_altered_preservation = NA_real_,
    passed_BHB_screen = NA
  ) %>%
  filter(!is.na(stage_group))

if (nrow(BHB_points) == 0) {
  stop("No screened BHB observations are available for the PhanDA comparison.")
}
if (any(!BHB_points$passed_BHB_screen)) {
  stop("An unscreened BHB observation entered the PhanDA comparison.")
}

comparison_points <- bind_rows(BHB_points, Tornillo_points) %>%
  mutate(stage_group = factor(stage_group, levels = stage_levels)) %>%
  arrange(stage_group, region, desc(Age_Ma), sample_id)

# Calculate display-only latitude offsets once so each point remains aligned
# with its uncertainty bar. The unjittered paleolatitude is retained in the
# audit file and continues to define the scientific location.
set.seed(47033)
comparison_points <- comparison_points %>%
  group_by(stage_group, region) %>%
  mutate(
    display_paleolatitude_deg_n =
      paleolatitude_deg_n + runif(n(), min = -0.65, max = 0.65)
  ) %>%
  ungroup()

write_csv(
  comparison_points,
  here(
    "data", "processed",
    "BHB_Tornillo_PhanDA_latitudinal_gradient_plot_data.csv"
  ),
  na = ""
)

#-- 4.) Plot ---------------------------------------------------------------
BHB_plot_points <- comparison_points %>% filter(region == "Bighorn Basin")
Tornillo_plot_points <- comparison_points %>% filter(region == "Tornillo Basin")

p_PhanDA_clumped_comparison <- ggplot() +
  geom_ribbon(
    data = PhanDA_ltg,
    aes(x = Latitude, ymin = LTG_05, ymax = LTG_95),
    fill = "grey70", alpha = 0.28
  ) +
  geom_ribbon(
    data = PhanDA_ltg,
    aes(x = Latitude, ymin = LTG_16, ymax = LTG_84),
    fill = "grey50", alpha = 0.34
  ) +
  geom_line(
    data = PhanDA_ltg,
    aes(x = Latitude, y = LTG_50),
    color = "grey20", linewidth = 0.9
  ) +
  geom_vline(
    xintercept = c(Tornillo_paleolat_deg_n, BHB_paleolat_deg_n),
    color = "grey55", linewidth = 0.35, linetype = "dotted"
  ) +
  geom_errorbar(
    data = BHB_plot_points,
    aes(
      x = display_paleolatitude_deg_n,
      ymin = temperature_C - temperature_se_C,
      ymax = temperature_C + temperature_se_C
    ),
    width = 0, linewidth = 0.35,
    color = "grey20", alpha = 0.65, na.rm = TRUE
  ) +
  geom_point(
    data = BHB_plot_points,
    aes(
      x = display_paleolatitude_deg_n,
      y = temperature_C,
      fill = p_altered_preservation
    ),
    shape = 21, size = 2.8,
    stroke = 0.55, color = "black"
  ) +
  geom_errorbar(
    data = Tornillo_plot_points,
    aes(
      x = display_paleolatitude_deg_n,
      ymin = temperature_C - temperature_se_C,
      ymax = temperature_C + temperature_se_C
    ),
    width = 0, linewidth = 0.35,
    color = "#7A4A00", alpha = 0.7, na.rm = TRUE
  ) +
  geom_point(
    data = Tornillo_plot_points,
    aes(x = display_paleolatitude_deg_n, y = temperature_C),
    shape = 24, size = 3,
    stroke = 0.6, fill = "#E69F00", color = "black"
  ) +
  annotate(
    "text", x = Tornillo_paleolat_deg_n, y = Inf,
    label = "Tornillo", vjust = 1.35, size = 3.2, fontface = "bold"
  ) +
  annotate(
    "text", x = BHB_paleolat_deg_n, y = Inf,
    label = "BHB", vjust = 1.35, size = 3.2, fontface = "bold"
  ) +
  facet_wrap(~stage_group, nrow = 1) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0.5, limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::label_percent(accuracy = 1),
    name = "Modeled BHB\nalteration probability"
  ) +
  scale_x_continuous(
    breaks = seq(20, 70, by = 10),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = seq(0, 50, by = 10),
    expand = expansion(mult = c(0.01, 0.055))
  ) +
  coord_cartesian(xlim = c(20, 70), ylim = c(-5, 55)) +
  labs(
    title = "BHB and Tornillo clumped temperatures in PhanDA context",
    subtitle = paste(
      "Only BHB observations passing the current d18Ocarb >= 20 per mil",
      "and T47 <= 50 deg C screen are shown"
    ),
    x = "Paleolatitude (\u00b0N)",
    y = "Temperature (\u00b0C)",
    caption = str_wrap(
      paste0(
        "PhanDA stage-level surface-air temperature: median line, ",
        "16\u201384% dark ribbon, 5\u201395% light ribbon. ",
        "Circles: BHB soil-carbonate T47; triangles: Kelson et al. Tornillo ",
        "micrite T47. Error bars are \u00b11 SE. Formation temperatures and ",
        "surface-air temperatures need not represent the same season or reservoir."
      ),
      width = 155
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5),
    plot.caption = element_text(size = 7.5, hjust = 0),
    panel.spacing.x = unit(1.1, "lines"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9)
  )

save_figure_variants(
  plot = p_PhanDA_clumped_comparison,
  base_dir = figure_dir,
  stem = "BHB_Tornillo_D47_PhanDA_latitudinal_gradient",
  manuscript_width = 8.2,
  manuscript_height = 4.8,
  presentation_width = 12,
  presentation_height = 6
)

message(
  "Saved PhanDA latitudinal-gradient comparison with ",
  nrow(BHB_points), " screened BHB and ",
  nrow(Tornillo_points), " Tornillo observations."
)
