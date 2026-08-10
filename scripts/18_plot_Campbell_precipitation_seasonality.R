# 18_plot_Campbell_precipitation_seasonality.R
# Purpose: Create a compact presentation schematic of modeled monthly
#          precipitation seasonality near Polecat Bench using the Campbell et
#          al. (2024) iCESM sensitivity experiments.
#
# Important: the 3x and 6x PI CO2 panels are forcing experiments on the same
# approximately 55 Ma paleogeography. They are not Thanetian and Ypresian time
# slices and must not be interpreted as a chronological reconstruction.

library(tidyverse)
library(here)

input_file <- here(
  "data", "raw", "campbell2024_polecat_monthly_precipitation.csv"
)
output_dir <- here("figures", "presentation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

monthly <- read_csv(input_file, show_col_types = FALSE) %>%
  group_by(co2_multiple, month, month_name) %>%
  summarise(
    precip_mean_mm = mean(precip_mm_month),
    precip_min_mm = min(precip_mm_month),
    precip_max_mm = max(precip_mm_month),
    .groups = "drop"
  ) %>%
  mutate(
    forcing = factor(
      co2_multiple,
      levels = c(3, 6),
      labels = c(
        "3x PI CO2 - summer-weighted",
        "6x PI CO2 - spring + autumn maxima"
      )
    )
  )

winter_labels <- tibble(
  forcing = factor(levels(monthly$forcing), levels = levels(monthly$forcing)),
  x = 1.55,
  y = 11,
  label = "Dry winter"
)

p_precip_seasonality <- ggplot(
  monthly,
  aes(x = month, y = precip_mean_mm)
) +
  annotate(
    "rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
    fill = "grey50", alpha = 0.09
  ) +
  annotate(
    "rect", xmin = 11.5, xmax = 12.5, ymin = -Inf, ymax = Inf,
    fill = "grey50", alpha = 0.09
  ) +
  geom_ribbon(
    aes(ymin = precip_min_mm, ymax = precip_max_mm),
    fill = "#9ECAE1", alpha = 0.48, color = NA
  ) +
  geom_line(color = "#17689A", linewidth = 1.15) +
  geom_point(
    shape = 21, size = 2.6, stroke = 0.65,
    color = "#0B405E", fill = "white"
  ) +
  geom_text(
    data = winter_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 10.5 / ggplot2::.pt,
    color = "grey35"
  ) +
  facet_wrap(vars(forcing), ncol = 1, strip.position = "top") +
  scale_x_continuous(
    breaks = seq(1, 11, by = 2),
    labels = month.abb[seq(1, 11, by = 2)],
    expand = expansion(mult = c(0.015, 0.015))
  ) +
  scale_y_continuous(
    breaks = seq(0, 180, 30),
    limits = c(0, 180),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Modeled precipitation seasonality",
    subtitle = "Polecat Bench region | ~55 Ma paleogeography",
    x = NULL,
    y = "Monthly precipitation (mm)",
    caption = str_wrap(
      paste0(
        "Line = mean; ribbon = range among four orbital configurations. ",
        "Sensitivity experiments, not a temporal sequence. ",
        "Campbell et al. (2024)."
      ),
      width = 80
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16.5, hjust = 0),
    plot.subtitle = element_text(size = 11.5, color = "grey30", hjust = 0),
    plot.caption = element_text(size = 8.5, color = "grey30", hjust = 0),
    axis.title.y = element_text(size = 14, margin = margin(r = 8)),
    axis.text = element_text(size = 11.5, color = "black"),
    axis.text.x = element_text(size = 11),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold", size = 13.5,
      hjust = 0, margin = margin(b = 4)
    ),
    panel.spacing.y = unit(1.05, "lines"),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 12, 8, 10)
  )

svg_file <- file.path(
  output_dir,
  "Campbell2024_Polecat_precipitation_seasonality_6x6.svg"
)
png_file <- file.path(
  output_dir,
  "Campbell2024_Polecat_precipitation_seasonality_6x6.png"
)

ggsave(
  svg_file, p_precip_seasonality,
  width = 6, height = 6, units = "in", device = svglite::svglite
)
ggsave(
  png_file, p_precip_seasonality,
  width = 6, height = 6, units = "in", dpi = 300, bg = "white"
)

message("Saved: ", svg_file)
