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
# ---- Temperature screening scheme ----

screening_scheme <- "exclude_all_suspect"
# options:
# "none"
# "exclude_heavy"
# "exclude_heavy_moderate"
# "exclude_all_suspect"

temp_screening_flags <- read_csv(
  here("data", "processed", "temperature_screening_flags.csv")
)

needed_screening_cols <- c(
  "exclude_heavy",
  "exclude_heavy_moderate",
  "exclude_all_suspect"
)

for (col in needed_screening_cols) {
  if (!col %in% names(temp_screening_flags)) {
    temp_screening_flags[[col]] <- FALSE
  }
}

BHB_multiproxy_summary <- BHB_multiproxy_summary %>%
  select(-any_of(c(
    needed_screening_cols,
    "exclude_from_temp_model"
  ))) %>%
  left_join(
    temp_screening_flags,
    by = "MLA_horizon_id"
  ) %>%
  mutate(
    across(
      any_of(needed_screening_cols),
      ~ replace_na(.x, FALSE)
    ),
    exclude_from_temp_model = case_when(
      screening_scheme == "none" ~ FALSE,
      screening_scheme == "exclude_heavy" ~ exclude_heavy,
      screening_scheme == "exclude_heavy_moderate" ~ exclude_heavy_moderate,
      screening_scheme == "exclude_all_suspect" ~ exclude_all_suspect,
      TRUE ~ FALSE
    )
  )
# ---- Build temperature observation table ----
# IPL, CU, and Snell are kept as separate observations.
# CU reports 2SE, so convert to SE.

temp_obs <- BHB_multiproxy_summary %>%
  filter(!exclude_from_temp_model) %>%
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
library(ggplot2)
library(dplyr)
library(here)

p_temp_teaser <- ggplot() +
  
  # PETM interval behind everything
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 1500, ymax = 1540,
    fill = "grey",
    alpha = 0.5
  ) +
  
  annotate(
    "text",
    x = 75,
    y = 1520,
    label = "PETM",
    hjust = 1,
    fontface = "bold",
    size = 4
  ) +
  
  # model uncertainty ribbon
  geom_ribbon(
    data = BHB_temperature_model,
    aes(
      xmin = T_model_lower95_C,
      xmax = T_model_upper95_C,
      y = strat_height_m
    ),
    fill = "firebrick",
    alpha = 0.18
  ) +
  
  # model line
  geom_path(
    data = BHB_temperature_model,
    aes(x = T_model_C, y = strat_height_m),
    color = "firebrick",
    linewidth = 1.2
  ) +
  
  # source-level SE bars
  geom_errorbarh(
    data = temp_obs,
    aes(
      xmin = T_C - T_se_C,
      xmax = T_C + T_se_C,
      y = strat_height_m,
      color = source
    ),
    height = 0,
    linewidth = 0.45,
    alpha = 0.65
  ) +
  
  # source-level observations
  geom_point(
    data = temp_obs,
    aes(
      x = T_C,
      y = strat_height_m,
      shape = source,
      fill = source,
      color = source
    ),
    size = 1.5,
    stroke = 0.8,
    alpha = 0.95
  ) +
  
  scale_shape_manual(
    values = c(
      "IPL" = 21,
      "CU" = 22,
      "Snell" = 24
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "IPL" = "black",
      "CU" = "white",
      "Snell" = "gray65"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "IPL" = "black",
      "CU" = "gray35",
      "Snell" = "gray35"
    )
  ) +
  
  scale_x_continuous(
    breaks = seq(10, 75, by = 5),
    limits = c(10, 75),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  
  scale_y_continuous(
    breaks = seq(500, 2300, by = 100),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  
  labs(
    x = expression(Delta[47] * "-derived temperature (" * degree * "C)"),
    y = "Stratigraphic height (m)",
    shape = "Data source",
    fill = "Data source",
    color = "Data source"
  ) +
  
  guides(
    color = guide_legend(override.aes = list(linewidth = 0, size = 3)),
    fill = guide_legend(override.aes = list(size = 3)),
    shape = guide_legend(override.aes = list(size = 3))
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = c(0.82, 0.18),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

p_temp_teaser


# PETM-focused
p_temp_teaser +
  coord_cartesian(
    ylim = c(1450, 1560)
  )

# LPEE focus
p_temp_teaser +
  coord_cartesian(
    ylim = c(1400, 1700)
  )

# Upper section
p_temp_teaser +
  coord_cartesian(
    ylim = c(1500, 2300)
  )

# Full stratigraphic section
ggsave(
  here("figures", "BHB_D47_temperature_full.png"),
  p_temp_teaser,
  width = 6,
  height = 7.2,
  dpi = 600
)

# PETM-focused
ggsave(
  here("figures", "BHB_D47_temperature_PETM.png"),
  p_temp_teaser +
    coord_cartesian(
      ylim = c(1450, 1560)
    ),
  width = 6,
  height = 5,
  dpi = 600
)

# LPEE focus
ggsave(
  here("figures", "BHB_D47_temperature_LPEE.png"),
  p_temp_teaser +
    coord_cartesian(
      ylim = c(1400, 1700)
    ),
  width = 6,
  height = 6,
  dpi = 600
)

# Upper section
ggsave(
  here("figures", "BHB_D47_temperature_Wasatchian.png"),
  p_temp_teaser +
    coord_cartesian(
      ylim = c(1500, 2300)
    ),
  width = 6,
  height = 6,
  dpi = 600
)
      
      
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

# interactive plot -----
library(plotly)
library(here)
library(htmlwidgets)

p_temp_interactive <- plot_ly() %>%
  add_trace(
    data = BHB_temperature_model,
    x = ~T_model_C,
    y = ~strat_height_m,
    type = "scatter",
    mode = "lines",
    line = list(width = 3),
    name = "Temperature model",
    text = ~paste0(
      "Horizon: ", MLA_horizon_id,
      "<br>Model T: ", round(T_model_C, 1), " °C",
      "<br>Strat height: ", strat_height_m, " m"
    ),
    hoverinfo = "text"
  ) %>%
  add_trace(
    data = temp_obs,
    x = ~T_C,
    y = ~strat_height_m,
    type = "scatter",
    mode = "markers",
    color = ~source,
    symbol = ~source,
    name = ~source,
    error_x = list(
      type = "data",
      array = ~T_se_C,
      visible = TRUE
    ),
    marker = list(size = 8),
    text = ~paste0(
      "Horizon: ", MLA_horizon_id,
      "<br>Source: ", source,
      "<br>T: ", round(T_C, 1), " °C",
      "<br>SE: ", round(T_se_C, 2), " °C",
      "<br>Strat height: ", strat_height_m, " m"
    ),
    hoverinfo = "text"
  ) %>%
  layout(
    shapes = list(
      list(
        type = "rect",
        xref = "paper",
        x0 = 0,
        x1 = 1,
        y0 = 1500,
        y1 = 1540,
        fillcolor = "red",
        opacity = 0.2,
        line = list(width = 0)
      )
    ),
    annotations = list(
      list(
        x = 1,
        y = 1520,
        xref = "paper",
        yref = "y",
        text = "PETM",
        showarrow = FALSE,
        xanchor = "right",
        font = list(size = 14)
      )
    ),
    xaxis = list(
      title = "Δ47-derived temperature (°C)"
    ),
    yaxis = list(
      title = "Stratigraphic height (m)"
    )
  )

p_temp_interactive

saveWidget(
  p_temp_interactive,
  here("figures", "interactive_temperature_model.html"),
  selfcontained = TRUE
)

# ---- Interactive age-space temperature plot ----

library(plotly)
library(here)
library(htmlwidgets)
library(dplyr)

# Use Matthew age model for all horizons
BHB_temp_age_plot <- BHB_multiproxy_with_temperature %>%
  mutate(
    plot_Age_Ma = Age_Ma,
    age_source = "Matthew age model"
  ) %>%
  filter(!is.na(plot_Age_Ma), !is.na(T_model_C))

temp_obs_age <- temp_obs %>%
  left_join(
    BHB_temp_age_plot %>%
      select(MLA_horizon_id, strat_height_m, plot_Age_Ma, age_source),
    by = c("MLA_horizon_id", "strat_height_m")
  ) %>%
  filter(!is.na(plot_Age_Ma), !is.na(T_C))

p_temp_age_interactive <- plot_ly() %>%
  add_trace(
    data = BHB_temp_age_plot,
    x = ~T_model_C,
    y = ~plot_Age_Ma,
    type = "scatter",
    mode = "lines",
    line = list(width = 3),
    name = "Temperature model",
    text = ~paste0(
      "Horizon: ", MLA_horizon_id,
      "<br>Age: ", round(plot_Age_Ma, 3), " Ma",
      "<br>Age source: ", age_source,
      "<br>Model T: ", round(T_model_C, 1), " °C",
      "<br>Strat height: ", strat_height_m, " m"
    ),
    hoverinfo = "text"
  ) %>%
  add_trace(
    data = temp_obs_age,
    x = ~T_C,
    y = ~plot_Age_Ma,
    type = "scatter",
    mode = "markers",
    color = ~source,
    symbol = ~source,
    name = ~source,
    error_x = list(
      type = "data",
      array = ~T_se_C,
      visible = TRUE
    ),
    marker = list(size = 8),
    text = ~paste0(
      "Horizon: ", MLA_horizon_id,
      "<br>Source: ", source,
      "<br>Age: ", round(plot_Age_Ma, 3), " Ma",
      "<br>Age source: ", age_source,
      "<br>T: ", round(T_C, 1), " °C",
      "<br>SE: ", round(T_se_C, 2), " °C",
      "<br>Strat height: ", strat_height_m, " m"
    ),
    hoverinfo = "text"
  ) %>%
  layout(
    xaxis = list(
      title = "Δ47-derived temperature (°C)"
    ),
    yaxis = list(
      title = "Age (Ma)",
      autorange = "reversed"
    )
  )

p_temp_age_interactive

saveWidget(
  p_temp_age_interactive,
  here("figures", "interactive_temperature_age_model.html"),
  selfcontained = TRUE
)