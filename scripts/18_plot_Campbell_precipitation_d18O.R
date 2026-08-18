# Campbell et al. (2024) modeled precipitation isotope annual cycle
#
# Extracts monthly precipitation d18O for the four orbital configurations at
# each of the 3x and 6x PI-CO2 experiments, summarizes annual means, and plots
# orbital envelopes. The four configurations are sensitivity end members, not
# a chronological sequence.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(here)
  library(readr)
})

input_file <- here(
  "data", "raw", "campbell2024_polecat_monthly_precipitation.csv"
)
output_dir <- here("data", "processed")
figure_dir <- here("figures", "Campbell2024_precipitation_d18O")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

monthly <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    co2_multiple = as.integer(co2_multiple),
    month_name = factor(month_name, levels = month.abb),
    experiment_label = case_when(
      co2_multiple == 3 ~ "3x PI CO2",
      co2_multiple == 6 ~ "6x PI CO2",
      TRUE ~ paste0(co2_multiple, "x PI CO2")
    )
  ) %>%
  arrange(co2_multiple, experiment, month)

# Preserve the extracted scenario-level monthly table.
write_csv(
  monthly,
  file.path(output_dir, "Campbell2024_Polecat_monthly_precipitation_d18O.csv")
)

orbital_envelope <- monthly %>%
  group_by(co2_multiple, experiment_label, month, month_name) %>%
  summarise(
    d18O_mean = mean(d18O_precip_permil_vsmow, na.rm = TRUE),
    d18O_min = min(d18O_precip_permil_vsmow, na.rm = TRUE),
    d18O_max = max(d18O_precip_permil_vsmow, na.rm = TRUE),
    d18O_range = d18O_max - d18O_min,
    precip_mean_mm = mean(precip_mm_month, na.rm = TRUE),
    .groups = "drop"
  )

# Annual precipitation-weighted isotope means. These are the appropriate annual
# means for precipitation, rather than the unweighted mean of monthly isotope
# values.
annual_by_experiment <- monthly %>%
  group_by(co2_multiple, experiment_label, experiment) %>%
  summarise(
    d18O_annual_precip_weighted = weighted.mean(
      d18O_precip_permil_vsmow, precip_mm_month, na.rm = TRUE
    ),
    d18O_annual_unweighted = mean(d18O_precip_permil_vsmow, na.rm = TRUE),
    annual_precip_mm = sum(precip_mm_month, na.rm = TRUE),
    .groups = "drop"
  )

annual_envelope <- annual_by_experiment %>%
  group_by(co2_multiple, experiment_label) %>%
  summarise(
    d18O_mean = mean(d18O_annual_precip_weighted),
    d18O_min = min(d18O_annual_precip_weighted),
    d18O_max = max(d18O_annual_precip_weighted),
    d18O_range = d18O_max - d18O_min,
    annual_precip_mean_mm = mean(annual_precip_mm),
    .groups = "drop"
  )

write_csv(
  orbital_envelope,
  file.path(output_dir, "Campbell2024_Polecat_monthly_precipitation_d18O_orbital_envelope.csv")
)
write_csv(
  annual_by_experiment,
  file.path(output_dir, "Campbell2024_Polecat_annual_precipitation_d18O_by_experiment.csv")
)
write_csv(
  annual_envelope,
  file.path(output_dir, "Campbell2024_Polecat_annual_precipitation_d18O_orbital_envelope.csv")
)

month_labels <- setNames(month.abb, seq_along(month.abb))

p_d18O <- ggplot(
  orbital_envelope,
  aes(x = month, y = d18O_mean)
) +
  geom_ribbon(
    aes(ymin = d18O_min, ymax = d18O_max),
    fill = "#78A9D1", alpha = 0.32, color = NA
  ) +
  geom_line(color = "#145A86", linewidth = 1.05) +
  geom_point(color = "#0B405E", fill = "white", shape = 21, size = 2.3) +
  facet_wrap(vars(experiment_label), ncol = 1) +
  scale_x_continuous(
    breaks = 1:12, labels = month_labels,
    expand = expansion(mult = c(0.015, 0.015))
  ) +
  labs(
    title = "Modeled Bighorn Basin precipitation d18O annual cycle",
    subtitle = "Campbell et al. (2024) Polecat Bench grid cells; ribbons span four orbital configurations",
    x = NULL,
    y = "d18O precip (per mil VSMOW)",
    caption = "Line = mean among orbital configurations; ribbon = orbital range. The scenarios are sensitivity end members, not a chronological sequence."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    plot.caption = element_text(size = 8.5, color = "grey30", hjust = 0),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

ggsave(
  file.path(figure_dir, "Campbell2024_Polecat_precipitation_d18O_annual_cycle.png"),
  p_d18O, width = 7.2, height = 7.0, dpi = 300, bg = "white"
)
ggsave(
  file.path(figure_dir, "Campbell2024_Polecat_precipitation_d18O_annual_cycle.pdf"),
  p_d18O, width = 7.2, height = 7.0, bg = "white"
)

print(annual_envelope)
message("Wrote Campbell precipitation d18O tables to: ", output_dir)
message("Wrote plot to: ", figure_dir)
