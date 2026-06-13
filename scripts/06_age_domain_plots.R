# ---- Load packages ----
library(tidyverse)
library(here)
library(patchwork)

# ---- Load processed BHB and reference datasets ----
BHB_temp <- read_csv(
  here("data", "processed", "BHB_multiproxy_with_temperature.csv")
)

Harper2024_CO2_SST <- read_csv(
  here("data", "processed", "Harper2024_CO2_SST_processed.csv")
)

# ---- Prepare BHB temperature data ----
# Use measured IPL T47 where available; otherwise use modeled/interpolated T.
# If you already created T_recon_C elsewhere, use that instead.

BHB_temp_plot <- BHB_temp %>%
  mutate(
    T_plot_C = case_when(
      !is.na(IPLD47_mean_T47_C) ~ IPLD47_mean_T47_C,
      TRUE ~ T_model_C
    ),
    T_plot_lower95_C = case_when(
      !is.na(IPLD47_mean_T47_C) & !is.na(IPLD47_se_T47_C) ~
        IPLD47_mean_T47_C - 1.96 * IPLD47_se_T47_C,
      TRUE ~ T_model_lower95_C
    ),
    T_plot_upper95_C = case_when(
      !is.na(IPLD47_mean_T47_C) & !is.na(IPLD47_se_T47_C) ~
        IPLD47_mean_T47_C + 1.96 * IPLD47_se_T47_C,
      TRUE ~ T_model_upper95_C
    ),
    T_source = case_when(
      !is.na(IPLD47_mean_T47_C) ~ "Measured IPL Δ47",
      !is.na(T_model_C) ~ "Modeled/interpolated",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(Age_Ma),
    !is.na(T_plot_C)
  )

# ---- Restrict Harper record to age range near BHB data ----
age_min <- min(BHB_temp_plot$Age_Ma, na.rm = TRUE)
age_max <- max(BHB_temp_plot$Age_Ma, na.rm = TRUE)

Harper_plot <- Harper2024_CO2_SST %>%
  filter(
    Age_Ma >= age_min,
    Age_Ma <= age_max
  )

# ---- Shared age-axis settings ----
age_breaks <- seq(
  floor(age_min * 10) / 10,
  ceiling(age_max * 10) / 10,
  by = 0.1
)
# ---- Panel A: BHB 3-point spline temperature model ----
p_BHB_temp_age <- ggplot(
  BHB_temp_plot,
  aes(x = T_model_C, y = Age_Ma)
) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 55.93, ymax = 56.05,
    fill = "grey70", alpha = 0.25
  ) +
  geom_ribbon(
    aes(
      xmin = T_model_lower95_C,
      xmax = T_model_upper95_C
    ),
    alpha = 0.25
  ) +
  geom_path(linewidth = 1) +
  geom_point(
    data = BHB_temp_plot %>% filter(!is.na(T_measured_C)),
    aes(x = T_measured_C, y = Age_Ma),
    size = 2,
    alpha = 0.75
  ) +
  scale_y_reverse(
    limits = c(age_max, age_min),
    breaks = age_breaks
  ) +
  labs(
    title = expression("BHB carbonate " * Delta[47] * " temperatures"),
    x = expression(T[47] ~ "(" * degree * "C)"),
    y = "Age (Ma)"
  ) +
  theme_classic()

# ---- Panel B: Harper Pacific SST ----
p_Harper_SST_age <- ggplot(
  Harper_plot,
  aes(x = Harper2024_mean_SST_C, y = Age_Ma)
) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 55.93, ymax = 56.05,
    fill = "grey70", alpha = 0.25
  ) +
  geom_ribbon(
    aes(
      xmin = Harper2024_SST_lower95_C,
      xmax = Harper2024_SST_upper95_C
    ),
    alpha = 0.25
  ) +
  geom_path(linewidth = 1) +
  scale_y_reverse(
    limits = c(age_max, age_min),
    breaks = age_breaks
  ) +
  labs(
    title = "Harper et al. (2024) Pacific SST",
    x = expression("SST (" * degree * "C)"),
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )
# ---- Combine panels on common age axis ----
p_temp_comparison_age <- p_BHB_temp_age + p_Harper_SST_age +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Bighorn Basin carbonate temperatures and Pacific SST on a common age axis"
  )

print(p_temp_comparison_age)

# ---- Save figure ----
ggsave(
  here("figures", "BHB_T47_Harper_SST_common_age_axis.png"),
  p_temp_comparison_age,
  width = 8,
  height = 6,
  dpi = 300
)
