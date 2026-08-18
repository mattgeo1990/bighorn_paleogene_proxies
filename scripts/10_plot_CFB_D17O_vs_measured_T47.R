# Plot reconstructed CFB soil-water Delta-prime-17O against measured T47.
# This is distinct from the related T_recon_C plot in 10_other_plots.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(here)
  library(readr)
})

input_file <- here(
  "data", "processed", "CFB_soilwater_reconstruction_summary.csv"
)
output_dir <- here("figures", "CFB_D17O_vs_measured_T47")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plot_data <- read_csv(input_file, show_col_types = FALSE) %>%
  filter(
    is.finite(T_measured_C),
    is.finite(D17Orsw_median_permeg),
    is.finite(T_measured_se_C),
    is.finite(D17Orsw_lower95_permeg),
    is.finite(D17Orsw_upper95_permeg)
  ) %>%
  mutate(
    T_measured_se_C = pmax(T_measured_se_C, 0),
    T_lower95_C = T_measured_C - 1.96 * T_measured_se_C,
    T_upper95_C = T_measured_C + 1.96 * T_measured_se_C
  )

if ("D17O_replicate_status" %in% names(plot_data)) {
  plot_data <- plot_data %>%
    mutate(sample_support = D17O_replicate_status)
} else {
  plot_data <- plot_data %>%
    mutate(sample_support = "not classified")
}

write_csv(
  plot_data,
  here("data", "processed", "CFB_D17O_vs_measured_T47_plot_data.csv")
)

p <- ggplot(
  plot_data,
  aes(x = T_measured_C, y = D17Orsw_median_permeg)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_errorbar(
    aes(ymin = D17Orsw_lower95_permeg, ymax = D17Orsw_upper95_permeg),
    width = 0, color = "grey45", alpha = 0.7
  ) +
  geom_errorbarh(
    aes(xmin = T_lower95_C, xmax = T_upper95_C),
    height = 0, color = "grey45", alpha = 0.7
  ) +
  geom_point(aes(shape = sample_support), size = 2.7, color = "#145A86") +
  geom_text(
    aes(label = MLA_horizon_id),
    nudge_y = 1.5, size = 2.6, check_overlap = TRUE
  ) +
  labs(
    title = "Reconstructed soil-water Delta-prime-17O versus measured T47",
    subtitle = "CFB horizons with paired measured clumped-isotope temperature and reconstructed soil water",
    x = "Measured T47 (degrees C)",
    y = "Reconstructed soil-water Delta-prime-17O (per meg)",
    shape = "Sample-level support",
    caption = "Points = medians; vertical and horizontal bars = approximate 95% intervals. This plot uses measured T47, not T_recon_C."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    plot.caption = element_text(size = 8.5, color = "grey30", hjust = 0),
    legend.position = "bottom"
  )

ggsave(
  file.path(output_dir, "CFB_D17O_vs_measured_T47.png"),
  p, width = 8.5, height = 6.2, dpi = 300, bg = "white"
)
ggsave(
  file.path(output_dir, "CFB_D17O_vs_measured_T47.pdf"),
  p, width = 8.5, height = 6.2, bg = "white"
)

message("Wrote measured-T47 plot to: ", output_dir)
