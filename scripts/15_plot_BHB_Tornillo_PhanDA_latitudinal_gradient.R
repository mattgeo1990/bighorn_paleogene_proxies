# 15_plot_BHB_Tornillo_PhanDA_latitudinal_gradient.R
# Purpose: Compare screened Bighorn Basin and Tornillo clumped-isotope
#          formation temperatures and published North American continental
#          temperature estimates with PhanDA stage-level latitudinal surface-
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
#   - Published continental estimates are plotted only when their reported age
#     can be assigned to one stage. West et al. high-latitude records reported
#     only as "late Paleocene to early Eocene" remain in the audit table but are
#     intentionally excluded from both stage-specific panels.

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

#-- 4.) Prepare Published North American Continental Constraints -----------
published_continental_source <- read_csv(
  here(
    "data", "raw",
    "published_North_America_Paleocene_Eocene_continental_temperature.csv"
  ),
  show_col_types = FALSE
)

Wing_LMA_points <- read_csv(
  here("data", "processed", "WingEtAl2000_LMA_MAT_processed.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    study = "Wing et al. (2000)",
    locality = "Bighorn Basin",
    interval_reported = sample,
    stage_assignment = stage_from_age(published_age_model_2_Ma),
    include_stage_plot = !is.na(stage_assignment),
    latitude_type = "paleolatitude",
    latitude_deg_n = BHB_paleolat_deg_n,
    temperature_target = "MAT",
    temperature_C = MAT_C,
    temperature_lower_C = MAT_C - MAT_error_C,
    temperature_upper_C = MAT_C + MAT_error_C,
    proxy_group = "Paleobotanical MAT",
    topographic_context = "basin",
    notes = paste0(
      "Leaf-margin MAT; published aggregate age model 2; ",
      "temperature interval is published +/- error"
    )
  )

published_continental_audit <- bind_rows(
  published_continental_source,
  Wing_LMA_points
) %>%
  mutate(
    stage_group = factor(stage_assignment, levels = stage_levels),
    plot_record = case_when(
      study == "Wing et al. (2000)" ~ "Wing et al. (2000) LMA",
      study == "West et al. (2020)" ~ "West et al. (2020) ensemble MAT",
      study == "Fricke and Wing (2004)" &
        proxy_group == "Vertebrate phosphate" ~
        "Fricke & Wing (2004) phosphate MAT",
      TRUE ~ "Fricke & Wing (2004) paleobotanical MAT"
    )
  )

published_continental_points <- published_continental_audit %>%
  filter(include_stage_plot, !is.na(stage_group)) %>%
  group_by(stage_group, plot_record, latitude_deg_n) %>%
  mutate(
    display_latitude_deg_n = latitude_deg_n +
      seq(-0.42, 0.42, length.out = n())
  ) %>%
  ungroup()

write_csv(
  published_continental_audit,
  here(
    "data", "processed",
    "published_North_America_stage_temperature_plot_audit.csv"
  ),
  na = ""
)

#-- 5.) Plot ---------------------------------------------------------------
BHB_plot_points <- comparison_points %>% filter(region == "Bighorn Basin")
Tornillo_plot_points <- comparison_points %>% filter(region == "Tornillo Basin")

published_record_shapes <- c(
  "Wing et al. (2000) LMA" = 0,
  "Fricke & Wing (2004) phosphate MAT" = 5,
  "Fricke & Wing (2004) paleobotanical MAT" = 1,
  "West et al. (2020) ensemble MAT" = 2
)

published_record_colors <- c(
  "Wing et al. (2000) LMA" = "#1B7837",
  "Fricke & Wing (2004) phosphate MAT" = "#762A83",
  "Fricke & Wing (2004) paleobotanical MAT" = "#A6611A",
  "West et al. (2020) ensemble MAT" = "#008C95"
)

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
    data = published_continental_points,
    aes(
      x = display_latitude_deg_n,
      ymin = temperature_lower_C,
      ymax = temperature_upper_C,
      color = plot_record
    ),
    width = 0, linewidth = 0.3, alpha = 0.42, na.rm = TRUE
  ) +
  geom_point(
    data = published_continental_points,
    aes(
      x = display_latitude_deg_n,
      y = temperature_C,
      color = plot_record,
      shape = plot_record
    ),
    size = 2.45, stroke = 0.72, fill = "white"
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
    label = "Tornillo", vjust = 1.35,
    size = 18 / ggplot2::.pt, fontface = "bold"
  ) +
  annotate(
    "text", x = BHB_paleolat_deg_n, y = Inf,
    label = "BHB", vjust = 1.35,
    size = 18 / ggplot2::.pt, fontface = "bold"
  ) +
  facet_wrap(~stage_group, nrow = 1) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0.5, limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::label_percent(accuracy = 1),
    name = "Modeled BHB\nalteration probability"
  ) +
  scale_color_manual(
    values = published_record_colors,
    breaks = names(published_record_colors),
    name = "Published North American continental estimates"
  ) +
  scale_shape_manual(
    values = published_record_shapes,
    breaks = names(published_record_shapes),
    name = "Published North American continental estimates"
  ) +
  scale_x_continuous(
    breaks = seq(20, 80, by = 10),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = seq(0, 50, by = 10),
    expand = expansion(mult = c(0.01, 0.055))
  ) +
  coord_cartesian(xlim = c(20, 82), ylim = c(-5, 55)) +
  labs(
    title = "North American continental temperatures in PhanDA context",
    subtitle = str_wrap(
      paste(
        "Stage-resolved published constraints are sparse in the Thanetian and",
        "substantially more numerous in the Ypresian"
      ),
      width = 95
    ),
    x = "Paleolatitude (\u00b0N)",
    y = "Temperature (\u00b0C)",
    caption = str_wrap(
      paste0(
        "PhanDA stage-level surface-air temperature: median line, ",
        "16\u201384% dark ribbon, 5\u201395% light ribbon. ",
        "Filled circles: screened BHB soil-carbonate T47; filled triangles: ",
        "Kelson et al. Tornillo micrite T47; open symbols: published continental ",
        "MAT estimates. Error definitions differ by source and are retained in ",
        "the audit table. West et al. high-latitude sites reported only as late ",
        "Paleocene to early Eocene are not assigned to either stage. West ",
        "latitudes are modern geographic coordinates; other points use published ",
        "paleolatitudes. Formation temperature and MAT are not equivalent targets."
      ),
      width = 110
    )
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 18),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 18),
    plot.caption = element_text(size = 18, hjust = 0),
    panel.spacing.x = unit(1.1, "lines"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  ) +
  guides(
    color = guide_legend(
      order = 1, ncol = 2, byrow = TRUE, title.position = "top"
    ),
    shape = guide_legend(
      order = 1, ncol = 2, byrow = TRUE, title.position = "top"
    ),
    fill = guide_colorbar(
      order = 2, direction = "horizontal", title.position = "top",
      barwidth = unit(4, "in")
    )
  )

# Pass an explicit presentation plot so the shared export helper does not
# replace the talk-sized typography with its smaller generic defaults.
p_PhanDA_clumped_comparison_presentation <-
  p_PhanDA_clumped_comparison +
  labs(title = NULL, subtitle = NULL, caption = NULL) +
  theme(
    strip.text = element_text(face = "bold", size = 20),
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 16),
    plot.caption = element_text(size = 12, hjust = 0),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20, color = "black"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title.position = "top",
    legend.justification = "center",
    legend.key.height = unit(1.15, "lines")
  )

save_figure_variants(
  plot = p_PhanDA_clumped_comparison,
  base_dir = figure_dir,
  stem = "BHB_Tornillo_D47_PhanDA_latitudinal_gradient",
  manuscript_width = 8.2,
  manuscript_height = 5.8,
  presentation_plot = p_PhanDA_clumped_comparison_presentation,
  presentation_width = 12,
  presentation_height = 8
)

message(
  "Saved PhanDA latitudinal-gradient comparison with ",
  nrow(BHB_points), " screened BHB and ",
  nrow(Tornillo_points), " Tornillo observations, plus ",
  nrow(published_continental_points),
  " stage-resolved published continental estimates."
)
