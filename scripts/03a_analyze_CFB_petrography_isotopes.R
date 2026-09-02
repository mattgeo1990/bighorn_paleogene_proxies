# 03a_analyze_CFB_petrography_isotopes.R
# Purpose: Join independently scored thin-section petrography to the
#          authoritative CFB isotope summary and generate reproducible
#          petrography-isotope diagnostic products.
#
# Petrographic observations are treated as independent descriptive evidence.
# This script does not revise alteration-screening flags or infer timing from
# crystal size alone.

#-- 1.) Setup ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(patchwork)
  library(ggrepel)
})
source(here("scripts", "helpers", "save_figure_variants.R"))

input_petrography <- here("data", "raw", "CFB_petrographic_scoring.csv")
input_isotopes <- here(
  "data", "processed", "CFB_soilcarb_isotope_summary.csv"
)
output_dir <- here("figures", "diagenetic_screening", "petrography")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_petrography_columns <- c(
  "SampleID", "dominant_fabric", "micrite_abundance",
  "microspar_abundance", "spar_abundance", "matrix_heterogeneity",
  "matrix_framework", "microfracture_density", "texture_gradients",
  "spar_distribution", "Notes"
)

#-- 2.) Read and validate inputs -------------------------------------------
petrography <- read_csv(input_petrography, show_col_types = FALSE)
isotopes <- read_csv(input_isotopes, show_col_types = FALSE)

missing_columns <- setdiff(required_petrography_columns, names(petrography))
if (length(missing_columns) > 0) {
  stop(
    "Petrography input is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}
if (nrow(petrography) != 18L || anyDuplicated(petrography$SampleID)) {
  stop("Expected 18 unique thin-section SampleID values.")
}
if (anyDuplicated(isotopes$MLA_horizon_id)) {
  stop("Isotope summary contains duplicated MLA_horizon_id values.")
}

# The microscope/order-form identifier uses 80n1, whereas the analytical
# summary uses 80-n1. Keep both identifiers in the joined product.
sample_id_aliases <- tribble(
  ~SampleID, ~MLA_horizon_id,
  "PK95-SC-80n1", "PK95-SC-80-n1"
)

petrography_isotopes <- petrography %>%
  left_join(sample_id_aliases, by = "SampleID") %>%
  mutate(MLA_horizon_id = coalesce(MLA_horizon_id, SampleID)) %>%
  left_join(
    isotopes %>%
      select(
        section_id, MLA_horizon_id, strat_height_m,
        d18Ocarb_vpdb = IPL_NuDog_d18Ocarb_VPDB,
        d18Ocarb_vsmow = IPL_NuDog_d18Ocarb_VSMOW,
        T47_C = IPLD47_mean_T47_C,
        T47_sd_C = IPLD47_sd_T47_C,
        T47_se_C = IPLD47_se_T47_C,
        T47_n = IPLD47_n_T47,
        D17O_primary_use
      ),
    by = "MLA_horizon_id"
  ) %>%
  mutate(
    dominant_fabric = factor(
      dominant_fabric,
      levels = c("micrite", "microspar", "mixed", "spar", "siliciclastics")
    ),
    matrix_heterogeneity = factor(
      matrix_heterogeneity, levels = c("low", "moderate", "high")
    ),
    texture_gradients = factor(texture_gradients, levels = c("no", "yes")),
    isotope_status = case_when(
      is.finite(d18Ocarb_vpdb) & is.finite(T47_C) ~ "Paired d18Ocarb and T47",
      is.finite(d18Ocarb_vpdb) ~ "d18Ocarb only",
      is.finite(T47_C) ~ "T47 only",
      TRUE ~ "No paired IPL values"
    ),
    carbonate_only = dominant_fabric != "siliciclastics"
  ) %>%
  arrange(strat_height_m, SampleID)

unmatched_samples <- petrography_isotopes %>%
  filter(is.na(section_id)) %>%
  pull(SampleID)
if (length(unmatched_samples) > 0) {
  stop(
    "Thin-section samples did not match the isotope horizon roster: ",
    paste(unmatched_samples, collapse = ", ")
  )
}

write_csv(
  petrography_isotopes,
  here("data", "processed", "CFB_petrography_isotope_comparison.csv")
)

#-- 3.) Reproducible descriptive summaries --------------------------------
summarize_category <- function(data, variable) {
  finite_summary <- function(x, fun) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else fun(x)
  }

  data %>%
    mutate(group_level = as.character(.data[[variable]])) %>%
    group_by(group_level) %>%
    summarize(
      n_thin_sections = n(),
      n_d18Ocarb = sum(is.finite(d18Ocarb_vpdb)),
      mean_d18Ocarb_vpdb = finite_summary(d18Ocarb_vpdb, mean),
      min_d18Ocarb_vpdb = finite_summary(d18Ocarb_vpdb, min),
      max_d18Ocarb_vpdb = finite_summary(d18Ocarb_vpdb, max),
      n_T47 = sum(is.finite(T47_C)),
      mean_T47_C = finite_summary(T47_C, mean),
      min_T47_C = finite_summary(T47_C, min),
      max_T47_C = finite_summary(T47_C, max),
      .groups = "drop"
    ) %>%
    mutate(
      group_variable = variable,
      .before = 1
    )
}

summary_variables <- c(
  "dominant_fabric", "matrix_heterogeneity", "matrix_framework",
  "microfracture_density", "texture_gradients", "spar_abundance",
  "spar_distribution"
)

category_summary <- map_dfr(
  summary_variables,
  ~ summarize_category(filter(petrography_isotopes, carbonate_only), .x)
)
write_csv(
  category_summary,
  here("data", "processed", "CFB_petrography_isotope_category_summary.csv")
)

paired_statistics <- bind_rows(
  petrography_isotopes %>%
    filter(is.finite(d18Ocarb_vpdb), is.finite(T47_C)) %>%
    summarize(
      subset = "All paired thin sections", n = n(),
      pearson_r = cor(d18Ocarb_vpdb, T47_C, method = "pearson"),
      spearman_rho = cor(d18Ocarb_vpdb, T47_C, method = "spearman")
    ),
  petrography_isotopes %>%
    filter(carbonate_only, is.finite(d18Ocarb_vpdb), is.finite(T47_C)) %>%
    summarize(
      subset = "Carbonate thin sections (SC-242 excluded)", n = n(),
      pearson_r = cor(d18Ocarb_vpdb, T47_C, method = "pearson"),
      spearman_rho = cor(d18Ocarb_vpdb, T47_C, method = "spearman")
    )
)
write_csv(
  paired_statistics,
  here("data", "processed", "CFB_petrography_isotope_statistics.csv")
)

#-- 4.) Shared plot settings -----------------------------------------------
fabric_shapes <- c(
  micrite = 21, microspar = 22, mixed = 23, spar = 24,
  siliciclastics = 25
)
heterogeneity_fills <- c(
  low = "#56B4E9", moderate = "#E69F00", high = "#D55E00"
)
fabric_labels <- c(
  micrite = "Micrite", microspar = "Microspar", mixed = "Mixed",
  spar = "Spar", siliciclastics = "Siliciclastics"
)
heterogeneity_labels <- c(
  low = "Low", moderate = "Moderate", high = "High"
)

theme_petro <- theme_classic(base_size = 18) +
  theme(
    legend.position = "top",
    legend.box = "vertical",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    plot.caption = element_text(size = 10, color = "grey30", hjust = 0),
    plot.margin = margin(7, 9, 7, 7)
  )

paired_plot_data <- petrography_isotopes %>%
  filter(is.finite(d18Ocarb_vpdb), is.finite(T47_C))
carbonate_plot_data <- paired_plot_data %>% filter(carbonate_only)

make_d18O_T47_plot <- function(data, title, subtitle, show_fit = FALSE) {
  p <- ggplot(data, aes(T47_C, d18Ocarb_vpdb)) +
    geom_errorbarh(
      aes(xmin = T47_C - T47_se_C, xmax = T47_C + T47_se_C),
      height = 0, linewidth = 0.4, alpha = 0.5, na.rm = TRUE
    ) +
    geom_point(
      aes(shape = dominant_fabric, fill = matrix_heterogeneity),
      size = 3.5, color = "black", stroke = 0.8
    ) +
    ggrepel::geom_text_repel(
      aes(label = str_remove(SampleID, "PK95-")),
      size = 3.2, seed = 123, min.segment.length = 0,
      box.padding = 0.35, max.overlaps = Inf, show.legend = FALSE
    ) +
    scale_shape_manual(
      values = fabric_shapes, limits = names(fabric_shapes),
      labels = fabric_labels, drop = FALSE
    ) +
    scale_fill_manual(
      values = heterogeneity_fills, limits = names(heterogeneity_fills),
      labels = heterogeneity_labels, drop = FALSE
    ) +
    labs(
      x = expression(Delta[47] * " temperature (" * degree * "C)"),
      y = expression(delta^18 * O[carbonate] ~ "(per mil VPDB)"),
      shape = "Dominant fabric", fill = "Matrix heterogeneity",
      title = title, subtitle = subtitle
    ) +
    guides(
      shape = guide_legend(nrow = 2, byrow = TRUE),
      fill = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    theme_petro

  if (show_fit) {
    p <- p + geom_smooth(
      method = "lm", formula = y ~ x, se = TRUE,
      color = "grey30", fill = "grey75", linewidth = 0.7,
      alpha = 0.25
    )
  }
  p
}

p_d18O_T47_all <- make_d18O_T47_plot(
  paired_plot_data,
  expression("Thin-section petrography in " * T[47] * "-" * delta^18 * O * " space"),
  "All paired samples; SC-242 is petrographically siliciclastic"
) + theme(
  legend.position = "bottom", legend.box = "vertical",
  legend.title = element_text(size = 11),
  legend.text = element_text(size = 10)
)

carbonate_r <- paired_statistics %>%
  filter(str_detect(subset, "Carbonate")) %>%
  pull(pearson_r)
p_d18O_T47_carbonate <- make_d18O_T47_plot(
  carbonate_plot_data,
  expression("Carbonate-only sensitivity: " * T[47] * "-" * delta^18 * O),
  paste0("SC-242 excluded; descriptive Pearson r = ", sprintf("%.2f", carbonate_r)),
  show_fit = TRUE
) + theme(legend.position = "none")

p_d18O_T47_sensitivity <- p_d18O_T47_all + p_d18O_T47_carbonate +
  plot_layout(guides = "keep") +
  plot_annotation(
    caption = paste(
      "Horizontal bars are T47 standard errors.",
      "Petrographic categories are descriptive and do not establish alteration timing."
    ),
    theme = theme(plot.caption = element_text(size = 10, hjust = 0))
  )

save_figure_variants(
  p_d18O_T47_sensitivity, output_dir,
  "CFB_petrography_d18O_T47_sensitivity",
  manuscript_width = 13, manuscript_height = 6.7,
  presentation_width = 12
)

#-- 5.) Category comparisons -----------------------------------------------
make_category_plot <- function(data, category, response, y_lab, title) {
  ggplot(data, aes(x = .data[[category]], y = .data[[response]])) +
    geom_errorbar(
      aes(
        ymin = if (response == "T47_C") T47_C - T47_se_C else .data[[response]],
        ymax = if (response == "T47_C") T47_C + T47_se_C else .data[[response]]
      ),
      width = 0.08, linewidth = 0.35, alpha = 0.4, na.rm = TRUE
    ) +
    geom_point(
      aes(shape = dominant_fabric, fill = matrix_heterogeneity),
      position = position_jitter(width = 0.10, height = 0, seed = 123),
      size = 3.1, color = "black", stroke = 0.7, na.rm = TRUE
    ) +
    stat_summary(
      fun = median, geom = "crossbar", width = 0.48,
      linewidth = 0.65, color = "#0072B2", na.rm = TRUE
    ) +
    scale_shape_manual(
      values = fabric_shapes, limits = names(fabric_shapes),
      breaks = c("micrite", "microspar", "mixed", "spar"),
      labels = fabric_labels[c("micrite", "microspar", "mixed", "spar")],
      drop = FALSE
    ) +
    scale_fill_manual(
      values = heterogeneity_fills, limits = names(heterogeneity_fills),
      labels = heterogeneity_labels, drop = FALSE
    ) +
    labs(x = str_to_sentence(str_replace_all(category, "_", " ")),
         y = y_lab, title = title, shape = "Dominant fabric",
         fill = "Matrix heterogeneity") +
    theme_petro +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

p_T47_fabric <- make_category_plot(
  filter(petrography_isotopes, carbonate_only),
  "dominant_fabric", "T47_C",
  expression(Delta[47] * " temperature (" * degree * "C)"),
  "T47 by dominant fabric"
)
p_T47_heterogeneity <- make_category_plot(
  filter(petrography_isotopes, carbonate_only),
  "matrix_heterogeneity", "T47_C",
  expression(Delta[47] * " temperature (" * degree * "C)"),
  "T47 by matrix heterogeneity"
)
p_d18O_fabric <- make_category_plot(
  filter(petrography_isotopes, carbonate_only),
  "dominant_fabric", "d18Ocarb_vpdb",
  expression(delta^18 * O[carbonate] ~ "(per mil VPDB)"),
  expression(delta^18 * O[carbonate] * " by dominant fabric")
)
p_d18O_heterogeneity <- make_category_plot(
  filter(petrography_isotopes, carbonate_only),
  "matrix_heterogeneity", "d18Ocarb_vpdb",
  expression(delta^18 * O[carbonate] ~ "(per mil VPDB)"),
  expression(delta^18 * O[carbonate] * " by matrix heterogeneity")
)

p_category_comparison <-
  (p_T47_fabric + p_T47_heterogeneity) /
  (p_d18O_fabric + p_d18O_heterogeneity) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "CFB thin-section petrography and isotope distributions",
    subtitle = paste(
      "Blue crossbars are category medians; points retain individual-sample",
      "identity and uncertainty"
    ),
    caption = paste(
      "SC-242 is excluded from carbonate-category summaries.",
      "SC-22 and SC-54 lack paired IPL values."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 13),
      plot.caption = element_text(size = 10, hjust = 0)
    )
  ) &
  theme(legend.position = "bottom", legend.box = "vertical")

save_figure_variants(
  p_category_comparison, output_dir,
  "CFB_petrography_isotope_category_comparison",
  manuscript_width = 12, manuscript_height = 10,
  presentation_width = 12
)

#-- 6.) Stratigraphic multipanel view --------------------------------------
strat_data <- petrography_isotopes %>% filter(is.finite(strat_height_m))
strat_range <- range(strat_data$strat_height_m, na.rm = TRUE)
shared_strat_scale <- scale_y_continuous(
  limits = strat_range,
  breaks = seq(ceiling(strat_range[1] / 200) * 200,
               floor(strat_range[2] / 200) * 200, by = 200),
  expand = expansion(mult = c(0.02, 0.02))
)

p_T47_strat <- ggplot(strat_data, aes(T47_C, strat_height_m)) +
  geom_errorbarh(
    aes(xmin = T47_C - T47_se_C, xmax = T47_C + T47_se_C),
    height = 0, linewidth = 0.35, alpha = 0.45, na.rm = TRUE
  ) +
  geom_point(
    aes(shape = dominant_fabric, fill = matrix_heterogeneity),
    size = 2.8, color = "black", stroke = 0.7, na.rm = TRUE
  ) +
  shared_strat_scale +
  scale_shape_manual(
    values = fabric_shapes, labels = fabric_labels, drop = FALSE
  ) +
  scale_fill_manual(
    values = heterogeneity_fills, labels = heterogeneity_labels, drop = FALSE
  ) +
  labs(
    x = expression(T[47] ~ "(" * degree * "C)"),
    y = "CFB stratigraphic height (m)", title = expression(T[47]),
    shape = "Dominant fabric", fill = "Matrix heterogeneity"
  ) + theme_petro +
  theme(legend.position = "bottom", legend.box = "vertical")

p_d18O_strat <- ggplot(strat_data, aes(d18Ocarb_vpdb, strat_height_m)) +
  geom_point(
    aes(shape = dominant_fabric, fill = matrix_heterogeneity),
    size = 2.8, color = "black", stroke = 0.7, na.rm = TRUE
  ) +
  shared_strat_scale +
  scale_shape_manual(
    values = fabric_shapes, labels = fabric_labels, drop = FALSE
  ) +
  scale_fill_manual(
    values = heterogeneity_fills, labels = heterogeneity_labels, drop = FALSE
  ) +
  labs(
    x = expression(delta^18 * O[carb] ~ "(per mil VPDB)"),
    y = NULL, title = expression(delta^18 * O[carbonate]),
    shape = "Dominant fabric", fill = "Matrix heterogeneity"
  ) + theme_petro + theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    legend.position = "none"
  )

p_fabric_strat <- ggplot(strat_data, aes(dominant_fabric, strat_height_m)) +
  geom_point(
    aes(shape = dominant_fabric, fill = matrix_heterogeneity),
    size = 3, color = "black", stroke = 0.7
  ) +
  shared_strat_scale +
  scale_shape_manual(
    values = fabric_shapes, labels = fabric_labels, drop = FALSE
  ) +
  scale_fill_manual(
    values = heterogeneity_fills, labels = heterogeneity_labels, drop = FALSE
  ) +
  labs(
    x = "Dominant fabric", y = NULL, title = "Fabric",
    shape = "Dominant fabric", fill = "Matrix heterogeneity"
  ) +
  theme_petro +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    legend.position = "none"
  )

p_gradient_strat <- ggplot(strat_data, aes(texture_gradients, strat_height_m)) +
  geom_point(
    aes(fill = matrix_heterogeneity), shape = 21,
    size = 3, color = "black", stroke = 0.7
  ) +
  shared_strat_scale +
  scale_fill_manual(
    values = heterogeneity_fills, labels = heterogeneity_labels, drop = FALSE
  ) +
  labs(
    x = "Texture gradients", y = NULL, title = "Gradients",
    fill = "Matrix heterogeneity"
  ) +
  theme_petro +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    legend.position = "none"
  )

p_stratigraphy <-
  p_T47_strat + p_d18O_strat + p_fabric_strat + p_gradient_strat +
  guide_area() +
  plot_layout(
    design = "ABCD\nEEEE",
    widths = c(1.1, 1.1, 1.25, 0.8), heights = c(1, 0.42),
    guides = "collect"
  ) +
  plot_annotation(
    title = "CFB thin-section petrography and isotopes in stratigraphic space",
    caption = paste(
      "Missing isotope values remain explicit.",
      "Fill shows matrix heterogeneity; shape shows dominant fabric."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 18),
      plot.caption = element_text(size = 10, hjust = 0)
    )
  )

save_figure_variants(
  p_stratigraphy, output_dir,
  "CFB_petrography_isotope_stratigraphy",
  manuscript_width = 13, manuscript_height = 9,
  presentation_width = 12
)

message(
  "Wrote petrography-isotope tables and figures for ",
  nrow(petrography_isotopes), " thin sections to: ", output_dir
)
