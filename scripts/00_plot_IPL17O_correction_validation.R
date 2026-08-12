# 00_plot_IPL17O_correction_validation.R
#
# Recreate the validation figures used by docs/IPL17O_correction_procedure.md.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(here)
  library(readr)
})

input_path <- here("data", "processed", "IPL17O_automated_reconciliation.csv")
figure_dir <- here("docs", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

reconciliation <- read_csv(input_path, show_col_types = FALSE) %>%
  mutate(
    reactor = factor(reactor_id),
    reference_short = case_when(
      reactor_id %in% c(31, 32) ~ "IPL R pipeline",
      reactor_id == 33 ~ "manual Excel",
      reactor_id == 34 ~ "BP-validated manual Excel",
      TRUE ~ "prior compiled result"
    )
  )

reactor_colors <- c(
  `31` = "#1b9e77", `32` = "#d95f02", `33` = "#7570b3",
  `34` = "#e7298a", `35` = "#66a61e", `36` = "#1f78b4"
)

theme_validation <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey30"),
    legend.position = "bottom"
  )

p_agreement <- ggplot(
  reconciliation,
  aes(x = legacy_standardized, y = new_automated, color = reactor)
) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.6, color = "grey35") +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_color_manual(values = reactor_colors) +
  coord_equal() +
  labs(
    title = "Automated corrections reproduce prior accepted results",
    subtitle = "The diagonal is exact agreement; all 91 comparisons are within 5 per meg",
    x = "Prior accepted Delta-prime-17O (per meg)",
    y = "New automated Delta-prime-17O (per meg)",
    color = "Reactor"
  ) +
  theme_validation

ggsave(
  file.path(figure_dir, "IPL17O_new_vs_prior.png"),
  p_agreement, width = 9, height = 5.7, dpi = 300, bg = "white"
)

p_residual <- ggplot(
  reconciliation,
  aes(x = reactor, y = new_minus_legacy, color = reactor)
) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -5, ymax = 5,
           fill = "#d9ead3", alpha = 0.45) +
  geom_hline(yintercept = 0, linewidth = 0.6, color = "grey30") +
  geom_hline(yintercept = c(-5, 5), linetype = "dashed", color = "#38761d") +
  geom_jitter(width = 0.14, height = 0, size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4,
               color = "black") +
  scale_color_manual(values = reactor_colors, guide = "none") +
  labs(
    title = "Differences from prior corrections by reactor",
    subtitle = "Black diamonds are reactor means; the green band is the +/-5 per meg validation criterion",
    x = "Reactor",
    y = "New minus prior accepted Delta-prime-17O (per meg)"
  ) +
  theme_validation

ggsave(
  file.path(figure_dir, "IPL17O_residuals_by_reactor.png"),
  p_residual, width = 8.2, height = 5.2, dpi = 300, bg = "white"
)

p_cache <- ggplot(
  reconciliation,
  aes(x = reactor, y = workbook_O2_difference, color = reactor)
) +
  geom_hline(yintercept = 0, linewidth = 0.6, color = "grey30") +
  geom_jitter(width = 0.14, height = 0, size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4,
               color = "black") +
  scale_color_manual(values = reactor_colors, guide = "none") +
  labs(
    title = "Fresh O2 reduction versus cached workbook formulas",
    subtitle = "Differences identify cached/external-link behavior; fresh raw-data calculations were retained",
    x = "Reactor",
    y = "Fresh minus cached Delta-prime-17O of O2 (per meg)"
  ) +
  theme_validation

ggsave(
  file.path(figure_dir, "IPL17O_fresh_vs_cached_O2.png"),
  p_cache, width = 8.2, height = 5.2, dpi = 300, bg = "white"
)

message("Wrote validation figures to: ", figure_dir)
