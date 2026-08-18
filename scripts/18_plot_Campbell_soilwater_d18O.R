# Campbell et al. (2024) modeled Polecat/Bighorn Basin soil-water d18O cycle.
# Four orbital configurations are shown as an envelope within each CO2 case.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(here)
  library(readr)
})

input_file <- here("data", "raw", "campbell2024_polecat_monthly_soilwater.csv")
output_dir <- here("data", "processed")
figure_dir <- here("figures", "Campbell2024_soilwater_d18O")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

monthly <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    co2_multiple = if_else(startsWith(experiment, "3x"), 3L, 6L),
    experiment_label = if_else(co2_multiple == 3, "3x PI CO2", "6x PI CO2"),
    month_name = factor(month.abb[month], levels = month.abb)
  ) %>%
  arrange(co2_multiple, experiment, month)

envelope <- monthly %>%
  group_by(co2_multiple, experiment_label, month, month_name) %>%
  summarise(
    d18O_mean = mean(d18O_soil_water_permil_vsmow),
    d18O_min = min(d18O_soil_water_permil_vsmow),
    d18O_max = max(d18O_soil_water_permil_vsmow),
    d18O_range = d18O_max - d18O_min,
    .groups = "drop"
  )

annual <- monthly %>%
  group_by(co2_multiple, experiment_label, experiment) %>%
  summarise(
    d18O_annual_mean = mean(d18O_soil_water_permil_vsmow),
    d18O_min = min(d18O_soil_water_permil_vsmow),
    d18O_max = max(d18O_soil_water_permil_vsmow),
    d18O_range = d18O_max - d18O_min,
    .groups = "drop"
  )

annual_envelope <- annual %>%
  group_by(co2_multiple, experiment_label) %>%
  summarise(
    d18O_mean = mean(d18O_annual_mean),
    d18O_min = min(d18O_annual_mean),
    d18O_max = max(d18O_annual_mean),
    d18O_range = d18O_max - d18O_min,
    .groups = "drop"
  )

write_csv(monthly, file.path(output_dir, "Campbell2024_Polecat_monthly_soilwater_d18O.csv"))
write_csv(envelope, file.path(output_dir, "Campbell2024_Polecat_monthly_soilwater_d18O_orbital_envelope.csv"))
write_csv(annual, file.path(output_dir, "Campbell2024_Polecat_annual_soilwater_d18O_by_experiment.csv"))
write_csv(annual_envelope, file.path(output_dir, "Campbell2024_Polecat_annual_soilwater_d18O_orbital_envelope.csv"))

p <- ggplot(envelope, aes(month, d18O_mean)) +
  geom_ribbon(aes(ymin = d18O_min, ymax = d18O_max), fill = "#78A9D1", alpha = 0.32) +
  geom_line(color = "#145A86", linewidth = 1.05) +
  geom_point(color = "#0B405E", fill = "white", shape = 21, size = 2.3) +
  facet_wrap(vars(experiment_label), ncol = 1) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    title = "Modeled Bighorn Basin soil-water d18O annual cycle",
    subtitle = "Campbell et al. (2024) Polecat Bench grid cells; ribbons span four orbital configurations",
    x = NULL, y = "d18O soil water (per mil VSMOW)",
    caption = "Line = mean among orbital configurations; ribbon = orbital range. Soil-water isotope composition is model output, not precipitation alone."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    plot.caption = element_text(size = 8.5, color = "grey30", hjust = 0),
    strip.background = element_blank(), strip.text = element_text(face = "bold", hjust = 0),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(figure_dir, "Campbell2024_Polecat_soilwater_d18O_annual_cycle.png"), p, width = 7.2, height = 7, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "Campbell2024_Polecat_soilwater_d18O_annual_cycle.pdf"), p, width = 7.2, height = 7, bg = "white")

print(annual_envelope)
