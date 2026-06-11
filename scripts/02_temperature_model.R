# 02_temperature_model.R
# Purpose: Build a 3-point running-average temperature-depth model
#          and predict/interpolate temperature at every BHB horizon.

# ---- Load packages ----
library(tidyverse)
library(here)
library(zoo)

# ---- Load processed multiproxy summary ----
BHB_multiproxy_summary <- read_csv(
  here("data", "processed", "BHB_multiproxy_summary.csv")
)

# ---- Build temperature observation table ----
# IPL, CU, and Snell are kept as separate observations.
# CU reports 2SE, so convert to SE.

temp_obs <- BHB_multiproxy_summary %>%
  select(
    MLA_horizon_id,
    strat_height_m,
    IPLD47_mean_T47_C,
    IPLD47_se_T47_C,
    CU_mean_T47_C,
    CU_2se_T47_C,
    Snell_mean_T47_C,
    Snell_se_T47_C
  ) %>%
  pivot_longer(
    cols = c(IPLD47_mean_T47_C, CU_mean_T47_C, Snell_mean_T47_C),
    names_to = "source",
    values_to = "T_C"
  ) %>%
  mutate(
    source = case_when(
      source == "IPLD47_mean_T47_C" ~ "IPL",
      source == "CU_mean_T47_C"     ~ "CU",
      source == "Snell_mean_T47_C"  ~ "Snell"
    ),
    T_se_C = case_when(
      source == "IPL"   ~ IPLD47_se_T47_C,
      source == "CU"    ~ CU_2se_T47_C / 2,
      source == "Snell" ~ Snell_se_T47_C
    )
  ) %>%
  select(
    MLA_horizon_id,
    strat_height_m,
    source,
    T_C,
    T_se_C
  ) %>%
  filter(!is.na(T_C), !is.na(strat_height_m)) %>%
  mutate(
    # Fill missing/zero SE so weighting does not break.
    T_se_C = if_else(
      is.na(T_se_C) | T_se_C == 0,
      median(T_se_C, na.rm = TRUE),
      T_se_C
    ),
    weight = 1 / T_se_C^2
  )

# ---- Quick QC ----
temp_obs %>% count(source)
temp_obs %>% arrange(strat_height_m)

# ---- Collapse to one measured temperature per horizon ----
# If multiple datasets exist at one horizon, use inverse-variance weighting.

temp_horizon <- temp_obs %>%
  group_by(MLA_horizon_id, strat_height_m) %>%
  summarise(
    T_measured_C = weighted.mean(T_C, w = weight, na.rm = TRUE),
    T_measured_se_C = sqrt(1 / sum(weight, na.rm = TRUE)),
    n_T_obs = n(),
    T_sources = paste(sort(unique(source)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(strat_height_m)

# ---- 3-point running-average temperature model ----
# This smooths local noise but preserves sharper structure better than a GAM.

temp_horizon <- temp_horizon %>%
  mutate(
    T_model_anchor_C = rollapply(
      T_measured_C,
      width = 3,
      FUN = mean,
      fill = NA,
      align = "center",
      partial = TRUE
    )
  )

# ---- Helper function for local model uncertainty ----
# Estimate uncertainty from the 3 nearest measured horizons.
# This is simple and local: it does not represent full model uncertainty.

get_local_T_se <- function(x, temp_horizon, k = 3) {
  nearest <- temp_horizon %>%
    mutate(distance_m = abs(strat_height_m - x)) %>%
    arrange(distance_m) %>%
    slice_head(n = k)
  
  sqrt(mean(nearest$T_measured_se_C^2, na.rm = TRUE))
}

# ---- Predict/interpolate temperature at every BHB horizon ----
prediction_grid <- BHB_multiproxy_summary %>%
  distinct(MLA_horizon_id, strat_height_m) %>%
  filter(!is.na(strat_height_m)) %>%
  arrange(strat_height_m)

# Smooth spline through 3-point running-average anchors.
# Lower spar = follows data more closely; higher spar = smoother.
temp_spline <- smooth.spline(
  x = temp_horizon$strat_height_m,
  y = temp_horizon$T_model_anchor_C,
  spar = 0.35
)

temp_interp <- predict(
  temp_spline,
  x = prediction_grid$strat_height_m
)

BHB_temperature_model <- prediction_grid %>%
  mutate(
    T_model_C = as.numeric(temp_interp$y),
    
    T_model_se_C = map_dbl(
      strat_height_m,
      get_local_T_se,
      temp_horizon = temp_horizon,
      k = 3
    ),
    
    T_model_lower95_C = T_model_C - 1.96 * T_model_se_C,
    T_model_upper95_C = T_model_C + 1.96 * T_model_se_C
  ) %>%
  left_join(
    temp_horizon,
    by = c("MLA_horizon_id", "strat_height_m")
  ) %>%
  mutate(
    has_measured_T47 = !is.na(T_measured_C),
    n_T_obs = replace_na(n_T_obs, 0)
  )

# ---- Join modeled temperature back to full multiproxy summary ----
BHB_multiproxy_with_temperature <- BHB_multiproxy_summary %>%
  left_join(
    BHB_temperature_model,
    by = c("MLA_horizon_id", "strat_height_m")
  )

# ---- Diagnostic plot ----

ggplot() +
  geom_point(
    data = temp_obs,
    aes(x = T_C, y = strat_height_m, shape = source),
    size = 2
  ) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = 53,
    ymin = 1500,
    ymax = 1540,
    fill = "red",
    alpha = .6
  ) +
  annotate(
    "text",
    x = Inf,
    y = 1520,
    label = "PETM",
    hjust = 1.1,
    fontface = "bold",
    size = 4
  ) +
  geom_errorbarh(
    data = temp_obs,
    aes(
      xmin = T_C - T_se_C,
      xmax = T_C + T_se_C,
      y = strat_height_m
    ),
    height = 0
  ) +
  geom_point(
    data = temp_horizon,
    aes(x = T_measured_C, y = strat_height_m),
    size = 2.5
  ) +
  geom_ribbon(
    data = BHB_temperature_model,
    aes(
      xmin = T_model_lower95_C,
      xmax = T_model_upper95_C,
      y = strat_height_m
    ),
    fill = "firebrick",
    alpha = 0.2
  ) +
  geom_path(
    data = BHB_temperature_model,
    aes(x = T_model_C, y = strat_height_m),
    color = "firebrick",
    linewidth = 1.2
  ) +
  labs(
    x = expression(Delta[47] * "-derived temperature (" * degree * "C)"),
    y = "Stratigraphic height (m)",
    shape = "Data source"
  ) +
  scale_x_continuous(
    breaks = seq(20, 60, by = 5)
  ) +
  scale_y_continuous(
    breaks = seq(500, 2300, by = 100)
  ) +
  theme_classic()

# ---- Save outputs ----
write_csv(
  temp_obs,
  here("data", "processed", "BHB_temperature_observations.csv")
)

write_csv(
  temp_horizon,
  here("data", "processed", "BHB_temperature_horizon_summary.csv")
)

write_csv(
  BHB_temperature_model,
  here("data", "processed", "BHB_temperature_model.csv")
)

write_csv(
  BHB_multiproxy_with_temperature,
  here("data", "processed", "BHB_multiproxy_with_temperature.csv")
)