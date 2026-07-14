# 07_strat_domain_plots.R
# Purpose: Generate stratigraphic-domain plots for BHB multiproxy manuscript/meeting figures

library(tidyverse)
library(here)
library(patchwork)
library(plotly)

# ---- Shared plot settings ----

y_limits <- c(500, 2250)
y_breaks <- seq(500, 2250, by = 100)

petm_ymin <- 1500
petm_ymax <- 1540

add_petm <- function(fill = "red", alpha = 0.25) {
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = petm_ymin, ymax = petm_ymax,
    fill = fill,
    alpha = alpha
  )
}

shared_y_scale <- function() {
  scale_y_continuous(
    limits = y_limits,
    breaks = y_breaks,
    expand = expansion(mult = c(0.01, 0.02))
  )
}

theme_strat_panel <- theme_classic(base_size = 10) +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.position = "none",
    plot.margin = margin(4, 4, 4, 4)
  )

# ---- Load data ----

BHB_multiproxy_summary <- read_csv(here("data", "processed", "BHB_multiproxy_summary.csv"))
BHB_d18Ow <- read_csv(here("data", "processed", "BHB_d18Ow_reconstruction.csv"))
BHB_temperature_model <- read_csv(here("data", "processed", "BHB_temperature_model.csv"))
temp_obs <- read_csv(here("data", "processed", "BHB_temperature_observations.csv"))
BHB_multiproxy_final <- read_csv(
  here("data", "processed", "BHB_multiproxy_final.csv")
)

bowen_reps <- read_csv(
  here("data", "raw", "Bowen2001_IsotopeData.csv")
) %>%
  mutate(
    strat_height_m = round(strat_height_m, 1),
    Bowen_d13Ccarb_vpdb = as.numeric(d13C_VPDB),
    Bowen_d18Ocarb_vpdb = as.numeric(d18Ocarb_VPDB),
    Bowen_d18Ocarb_vsmow = as.numeric(d18Ocarb_VSMOW)
  ) %>%
  filter(
    !is.na(strat_height_m),
    !is.na(Bowen_d13Ccarb_vpdb),
    Bowen_d18Ocarb_vpdb > -11,
    !grepl("SPAR", MLA_sample_id)
  )

# ------PLOTS ------
# all d13C by strat -------
library(tidyverse)

d13c_strat <- BHB_multiproxy_summary %>%
  select(
    strat_height_m,
    Koch_mean_d13Ccarb_vpdb,
    Bowen_mean_d13Ccarb_vpdb,
    CU_mean_d13Ccarb_vpdb,
    Snell_mean_d13Ccarb_vpdb,
    IPL_NuDog_d13Ccarb_VPDB
  ) %>%
  pivot_longer(
    cols = -strat_height_m,
    names_to = "dataset",
    values_to = "d13C"
  ) %>%
  filter(!is.na(d13C))

ggplot(
  d13c_strat,
  aes(
    x = d13C,
    y = strat_height_m,
    color = dataset
  )
) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    x = expression(delta^13 * C[carb] ~ "(‰ VPDB)"),
    y = "Stratigraphic height (m)",
    color = "Dataset",
    title = expression(paste(delta^13, "C datasets by stratigraphic height"))
  ) +
  theme_classic()

# Standalone: Polecat Bench d13C Bowen et al. (2001) --------

# ---- Bowen 2001 primary-only replicate-level d13C strat plot ----
bowen_reps <- read_csv(
  here("data", "raw", "Bowen2001_IsotopeData.csv")
) %>%
  mutate(
    strat_height_m = round(strat_height_m, 1),
    Bowen_d13Ccarb_vpdb = as.numeric(d13C_VPDB),
    Bowen_d18Ocarb_vpdb = as.numeric(d18Ocarb_VPDB),
    Bowen_d18Ocarb_vsmow = as.numeric(d18Ocarb_VSMOW)
  ) %>%
  filter(
    !is.na(strat_height_m),
    !is.na(Bowen_d13Ccarb_vpdb),
    !is.na(Bowen_d18Ocarb_vpdb),
    Bowen_d18Ocarb_vpdb > -11,
    !grepl("SPAR", MLA_sample_id)
  )

p_bowen_d13C <- ggplot(
  bowen_reps,
  aes(
    x = strat_height_m,
    y = Bowen_d13Ccarb_vpdb
  )
) +
  geom_point(
    aes(
      tooltip = paste0(
        "Soil ID: ", Bowen_Soil_ID,
        "<br>Sample ID: ", Bowen_Sample_ID,
        "<br>MLA sample: ", MLA_sample_id,
        "<br>Horizon: ", MLA_horizon_id,
        "<br>Strat height: ", strat_height_m, " m",
        "<br>d13C: ", round(Bowen_d13Ccarb_vpdb, 2), " per mil VPDB",
        "<br>d18O: ", round(Bowen_d18Ocarb_vpdb, 2), " per mil VPDB"
      )
    ),
    size = 1.8,
    alpha = 0.65,
    color = "gray45"
  ) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    span = 0.06,
    se = FALSE,
    linewidth = 1.5,
    color = "black"
  ) +
  coord_flip() +
  scale_x_continuous(
    breaks = seq(
      floor(min(bowen_reps$strat_height_m, na.rm = TRUE)),
      ceiling(max(bowen_reps$strat_height_m, na.rm = TRUE)),
      by = 1
    )
  ) +
  labs(
    title = "Bowen et al. primary carbonate d13C",
    subtitle = "Replicate-level data; black curve = LOESS smooth through stratigraphic height",
    x = "Stratigraphic height above K-Pg (m)",
    y = "d13Ccarb (per mil VPDB)"
  ) +
  theme_classic(base_size = 14)

ggplotly(
  p_bowen_d13C,
  tooltip = "tooltip",
  height = 1400,
  width = 900
)


# ---- Temperature plot ----

temp_obs_all <- temp_obs %>%
  filter(!is.na(T_C), !is.na(strat_height_m))

p_T47_strat_clean <- ggplot() +
  geom_ribbon(
    data = BHB_temperature_model,
    aes(
      xmin = T_model_lower95_C,
      xmax = T_model_upper95_C,
      y = strat_height_m
    ),
    fill = "grey70",
    alpha = 0.45
  ) +
  add_petm() +
  geom_path(
    data = BHB_temperature_model,
    aes(x = T_model_C, y = strat_height_m),
    linewidth = 1
  ) +
  geom_errorbarh(
    data = temp_obs_all,
    aes(
      xmin = T_C - T_se_C,
      xmax = T_C + T_se_C,
      y = strat_height_m
    ),
    height = 0,
    linewidth = 0.35,
    alpha = 0.45
  ) +
  geom_point(
    data = temp_obs_all,
    aes(x = T_C, y = strat_height_m),
    size = 1.8,
    alpha = 0.75
  ) +
  scale_x_continuous(
    breaks = seq(10, 75, by = 10),
    limits = c(10, 75),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  shared_y_scale() +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = "Stratigraphic height (m)"
  ) +
  theme_strat_panel

# ---- d13C carbonate plot ----

d13C_obs_all <- BHB_multiproxy_summary %>%
  select(
    MLA_horizon_id,
    strat_height_m,
    Koch_mean_d13Ccarb_vpdb,
    Koch_se_d13Ccarb_vpdb,
    Bowen_mean_d13Ccarb_vpdb,
    Bowen_se_d13Ccarb_vpdb,
    CU_mean_d13Ccarb_vpdb,
    Snell_mean_d13Ccarb_vpdb
  ) %>%
  pivot_longer(
    cols = c(
      Koch_mean_d13Ccarb_vpdb,
      Bowen_mean_d13Ccarb_vpdb,
      CU_mean_d13Ccarb_vpdb,
      Snell_mean_d13Ccarb_vpdb
    ),
    names_to = "source",
    values_to = "d13Ccarb_vpdb"
  ) %>%
  mutate(
    d13C_se_vpdb = case_when(
      source == "Koch_mean_d13Ccarb_vpdb" ~ Koch_se_d13Ccarb_vpdb,
      source == "Bowen_mean_d13Ccarb_vpdb" ~ Bowen_se_d13Ccarb_vpdb,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(d13Ccarb_vpdb), !is.na(strat_height_m))

p_d13C_strat_clean <- ggplot() +
  add_petm() +
  geom_errorbarh(
    data = d13C_obs_all,
    aes(
      xmin = d13Ccarb_vpdb - d13C_se_vpdb,
      xmax = d13Ccarb_vpdb + d13C_se_vpdb,
      y = strat_height_m
    ),
    height = 0,
    linewidth = 0.35,
    alpha = 0.45,
    na.rm = TRUE
  ) +
  geom_point(
    data = d13C_obs_all,
    aes(x = d13Ccarb_vpdb, y = strat_height_m),
    size = 1.8,
    alpha = 0.75
  ) +
  scale_x_continuous(
    breaks = seq(-16, -4, by = 2),
    limits = c(-16, -4),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  shared_y_scale() +
  labs(
    x = expression(delta^13 * C[carb] ~ "(‰ VPDB)"),
    y = "Stratigraphic height (m)"
  ) +
  theme_strat_panel

# Compact rolling-mean d13Ccarb stratigraphic panel -----

library(dplyr)
library(ggplot2)
install.packages("slider")
library(slider)
library(svglite)

d13c_panel_width_in  <- 1
d13c_panel_height_in <- 3

# Plot limits
d13c_y_limits <- c(0, 2300)
d13c_y_breaks <- seq(0, 2300, by = 200)

# Rolling-window width in stratigraphic meters
rolling_window_m <- 100

# Calculate centered rolling mean by stratigraphic height
d13C_rolling <- d13C_obs_all %>%
  filter(
    !is.na(strat_height_m),
    !is.na(d13Ccarb_vpdb),
    strat_height_m >= d13c_y_limits[1],
    strat_height_m <= d13c_y_limits[2]
  ) %>%
  arrange(strat_height_m) %>%
  group_by(strat_height_m) %>%
  summarise(
    d13Ccarb_vpdb = mean(d13Ccarb_vpdb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    d13C_rolling_mean = slide_index_dbl(
      .x = d13Ccarb_vpdb,
      .i = strat_height_m,
      .f = ~ mean(.x, na.rm = TRUE),
      .before = rolling_window_m / 2,
      .after = rolling_window_m / 2,
      .complete = FALSE
    )
  )

# Horizontal limits based on both observations and rolling mean
d13c_x_limits <- range(
  c(
    d13C_obs_all$d13Ccarb_vpdb,
    d13C_rolling$d13C_rolling_mean
  ),
  na.rm = TRUE
)

d13c_x_padding <- diff(d13c_x_limits) * 0.04

d13c_x_limits <- c(
  d13c_x_limits[1] - d13c_x_padding,
  d13c_x_limits[2] + d13c_x_padding
)

p_d13C_compact <- ggplot() +
  
  # Raw observations shown lightly for context
  geom_point(
    data = d13C_obs_all %>%
      filter(
        !is.na(d13Ccarb_vpdb),
        !is.na(strat_height_m)
      ),
    aes(
      x = d13Ccarb_vpdb,
      y = strat_height_m
    ),
    size = 0.75,
    alpha = 0.20
  ) +
  
  # Rolling mean line
  geom_path(
    data = d13C_rolling,
    aes(
      x = d13C_rolling_mean,
      y = strat_height_m
    ),
    linewidth = 0.65,
    lineend = "round",
    na.rm = TRUE
  ) +
  
  scale_x_continuous(
    limits = d13c_x_limits,
    breaks = seq(-16, -4, by = 4),
    expand = expansion(mult = c(0, 0))
  ) +
  
  scale_y_continuous(
    limits = d13c_y_limits,
    breaks = d13c_y_breaks,
    expand = expansion(mult = c(0, 0))
  ) +
  
  coord_cartesian(
    xlim = d13c_x_limits,
    ylim = d13c_y_limits,
    expand = FALSE,
    clip = "on"
  ) +
  
  labs(
    x = expression(delta^13 * C[carb] ~ ("\u2030 VPDB")),
    y = NULL
  ) +
  
  theme_classic(
    base_family = "Arial",
    base_size = 6.5
  ) +
  
  theme(
    axis.title.x = element_text(
      family = "Arial",
      size = 6.5,
      margin = margin(t = 2)
    ),
    
    axis.text.x = element_text(
      family = "Arial",
      size = 6
    ),
    
    axis.ticks.x = element_line(
      linewidth = 0.25
    ),
    
    axis.ticks.length.x = grid::unit(
      1,
      "mm"
    ),
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    
    axis.line = element_line(
      linewidth = 0.28
    ),
    
    legend.position = "none",
    
    plot.margin = margin(
      t = 0,
      r = 2,
      b = 2,
      l = 0,
      unit = "pt"
    )
  )

p_d13C_compact

# -------------------------------------------------------------------
# Export
# -------------------------------------------------------------------

dir.create(
  here("figures", "strat_domain"),
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename = here(
    "figures",
    "strat_domain",
    "BHB_d13Ccarb_compact.png"
  ),
  plot = p_d13C_compact,
  width = d13c_panel_width_in,
  height = d13c_panel_height_in,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = here(
    "figures",
    "strat_domain",
    "BHB_d13Ccarb_compact.svg"
  ),
  plot = p_d13C_compact,
  width = d13c_panel_width_in,
  height = d13c_panel_height_in,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

# ---- d18O soil-water plot ----

d18Ow_obs <- BHB_d18Ow %>%
  filter(!is.na(d18Ow_mean_vsmow), !is.na(strat_height_m))

p_d18Ow_strat_clean <- ggplot() +
  add_petm() +
  geom_errorbarh(
    data = d18Ow_obs,
    aes(
      xmin = d18Ow_lower95_vsmow,
      xmax = d18Ow_upper95_vsmow,
      y = strat_height_m
    ),
    height = 0,
    linewidth = 0.30,
    alpha = 0.35,
    na.rm = TRUE
  ) +
  geom_point(
    data = d18Ow_obs,
    aes(x = d18Ow_mean_vsmow, y = strat_height_m),
    size = 1.6,
    alpha = 0.70
  ) +
  scale_x_continuous(
    breaks = seq(-10, 2, by = 2),
    limits = c(-10, 2),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  shared_y_scale() +
  labs(
    x = expression(delta^18 * O[soil-water] ~ "(‰ VSMOW)"),
    y = "Stratigraphic height (m)"
  ) +
  theme_strat_panel

# ---- D17O soil-water strat plot ----

D17Orsw_obs <- BHB_multiproxy_final %>%
  filter(
    !is.na(D17Orsw_mean_permeg),
    !is.na(strat_height_m)
  )

p_D17Orsw_strat_clean <- ggplot() +
  add_petm() +
  geom_errorbarh(
    data = D17Orsw_obs,
    aes(
      xmin = D17Orsw_mean_permeg - IPL17O_sd_Dp17Ocarb_adj,
      xmax = D17Orsw_mean_permeg + IPL17O_sd_Dp17Ocarb_adj,
      y = strat_height_m
    ),
    height = 0,
    linewidth = 0.30,
    alpha = 0.35,
    na.rm = TRUE
  ) +
  geom_point(
    data = D17Orsw_obs,
    aes(x = D17Orsw_mean_permeg, y = strat_height_m),
    size = 1.6,
    alpha = 0.70
  ) +
  scale_x_continuous(
    breaks = seq(-200, 200, by = 50),
    limits = c(-200, 200),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  shared_y_scale() +
  labs(
    x = expression(Delta*minute^17*O[soil-water] ~ "(per meg)"),
    y = "Stratigraphic height (m)"
  ) +
  theme_strat_panel

p_D17Orsw_strat_clean
# ---- Biozone side panel ----

biozones <- tribble(
  ~Zone,  ~ymin, ~ymax,
  "Wa-5", 2200, 2240,
  "Wa-4", 2020, 2200,
  "Wa-3b",1780, 2020,
  "Wa-3a",1750, 1780,
  "Wa-2", 1645, 1750,
  "Wa-1", 1543, 1645,
  "Wa-0", 1506, 1543,
  "Cf-3", 1335, 1506,
  "Cf-2", 1180, 1335,
  "Cf-1", 885,  1180,
  "Ti-6", 820,  885,
  "Ti-5b",655,  820,
  "Ti-5a",530,  655,
  "Ti-4", 415,  530,
  "Ti-3", 215,  415,
  "Ti-2", 155,  215
)

p_biozone_strat <- ggplot(biozones) +
  geom_rect(
    aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax),
    fill = "grey90",
    color = "grey45",
    linewidth = 0.25
  ) +
  geom_text(
    aes(x = 0.5, y = (ymin + ymax) / 2, label = Zone),
    size = 2.5
  ) +
  scale_y_continuous(
    limits = y_limits,
    breaks = y_breaks,
    expand = c(0, 0)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 10) +
  theme(
    plot.margin = margin(4, 2, 4, 4)
  )

# ---- Panel figures ----

p_T47_panel <- p_T47_strat_clean +
  labs(tag = "A")

p_d13C_panel <- p_d13C_strat_clean +
  labs(y = NULL, tag = "B") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_d18Ow_panel <- p_d18Ow_strat_clean +
  labs(y = NULL, tag = "C") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_strat_panel_T_d13C <- p_biozone_strat + p_T47_panel + p_d13C_panel +
  plot_layout(
    nrow = 1,
    ncol = 3,
    widths = c(0.25, 1, 1)
  ) +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(face = "bold", size = 12)
    )
  )

p_strat_panel_T_d13C_d18Ow <- p_biozone_strat + p_T47_panel + p_d13C_panel + p_d18Ow_panel +
  plot_layout(
    nrow = 1,
    ncol = 4,
    widths = c(0.25, 1, 1, 1)
  ) +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(face = "bold", size = 12)
    )
  )

p_D17Orsw_panel <- p_D17Orsw_strat_clean +
  labs(y = NULL, tag = "C") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_strat_panel_T_d18Ow_D17Orsw <- 
  p_biozone_strat + p_T47_panel + p_d18Ow_panel + p_D17Orsw_panel +
  plot_layout(
    nrow = 1,
    ncol = 4,
    widths = c(0.25, 1, 1, 1)
  ) +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(face = "bold", size = 12)
    )
  )


p_strat_panel_T_d13C
p_strat_panel_T_d13C_d18Ow
p_strat_panel_T_d18Ow_D17Orsw

# ---- Save outputs ----


dir.create(here("figures", "strat_domain"), recursive = TRUE, showWarnings = FALSE)


output_type <- "PPT"   # "PPT" or "MS"

if (output_type == "PPT") {
  base_size <- 16
  point_size <- 3.5
  line_width <- 1.3
} else {
  base_size <- 10
  point_size <- 1.8
  line_width <- 0.8
}


ggsave(
  here("figures", "strat_domain", paste0("BHB_strat_panel_T47_d13C_", output_type, ".png")),
  p_strat_panel_T_d13C,
  width = 11.5,
  height = 6,
  dpi = 600
)

ggsave(
  here("figures", "strat_domain", "BHB_strat_panel_T47_d13C_d18Ow.png"),
  p_strat_panel_T_d13C_d18Ow,
  width = 9.8,
  height = 6,
  dpi = 600
)

ggsave(
  here("figures", "strat_domain", "BHB_strat_panel_T47_d18Ow_D17Orsw.png"),
  p_strat_panel_T_d18Ow_D17Orsw,
  width = 12,
  height = 6,
  dpi = 600
)

# ---- P-E boundary zoom panels: 1200–1800 m ----

pe_y_limits <- c(1200, 1800)
pe_y_breaks <- seq(1200, 1800, by = 100)

shared_pe_y_scale <- function() {
  scale_y_continuous(
    limits = pe_y_limits,
    breaks = pe_y_breaks,
    expand = expansion(mult = c(0.01, 0.02))
  )
}

# Rebuild biozone panel for zoomed interval
p_biozone_strat_PE <- ggplot(biozones) +
  geom_rect(
    aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax),
    fill = "grey90",
    color = "grey45",
    linewidth = 0.25
  ) +
  geom_text(
    aes(x = 0.5, y = (ymin + ymax) / 2, label = Zone),
    size = 2.8
  ) +
  scale_y_continuous(
    limits = pe_y_limits,
    breaks = pe_y_breaks,
    expand = c(0, 0)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 10) +
  theme(
    plot.margin = margin(4, 2, 4, 4)
  )

# Apply zoom scale to each proxy panel
p_T47_panel_PE <- p_T47_strat_clean +
  shared_pe_y_scale() +
  labs(tag = "A")

p_d13C_panel_PE <- p_d13C_strat_clean +
  shared_pe_y_scale() +
  labs(y = NULL, tag = "B") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_d18Ow_panel_PE <- p_d18Ow_strat_clean +
  shared_pe_y_scale() +
  labs(y = NULL, tag = "C") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_D17Orsw_panel_PE <- p_D17Orsw_strat_clean +
  shared_pe_y_scale() +
  labs(y = NULL, tag = "D") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Panel: T47 + d13C
p_strat_panel_T_d13C_PE <- 
  p_biozone_strat_PE + p_T47_panel_PE + p_d13C_panel_PE +
  plot_layout(
    nrow = 1,
    ncol = 3,
    widths = c(0.25, 1, 1)
  ) +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(face = "bold", size = 12)
    )
  )

# Panel: T47 + d13C + d18Ow
p_strat_panel_T_d13C_d18Ow_PE <- 
  p_biozone_strat_PE + p_T47_panel_PE + p_d13C_panel_PE + p_d18Ow_panel_PE +
  plot_layout(
    nrow = 1,
    ncol = 4,
    widths = c(0.25, 1, 1, 1)
  ) +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(face = "bold", size = 12)
    )
  )

# Panel: T47 + d18Ow + D17Orsw
p_strat_panel_T_d18Ow_D17Orsw_PE <- 
  p_biozone_strat_PE + p_T47_panel_PE + p_d18Ow_panel_PE + p_D17Orsw_panel_PE +
  plot_layout(
    nrow = 1,
    ncol = 4,
    widths = c(0.25, 1, 1, 1)
  ) +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(face = "bold", size = 12)
    )
  )

p_strat_panel_T_d13C_PE
p_strat_panel_T_d13C_d18Ow_PE
p_strat_panel_T_d18Ow_D17Orsw_PE

ggsave(
  here("figures", "strat_domain", "BHB_strat_panel_T47_d13C_PE_zoom.png"),
  p_strat_panel_T_d13C_PE,
  width = 11.5,
  height = 6,
  dpi = 600
)

ggsave(
  here("figures", "strat_domain", "BHB_strat_panel_T47_d13C_d18Ow_PE_zoom.png"),
  p_strat_panel_T_d13C_d18Ow_PE,
  width = 12,
  height = 6,
  dpi = 600
)

ggsave(
  here("figures", "strat_domain", "BHB_strat_panel_T47_d18Ow_D17Orsw_PE_zoom.png"),
  p_strat_panel_T_d18Ow_D17Orsw_PE,
  width = 12,
  height = 6,
  dpi = 600
)

write_csv(
  biozones,
  here(
    "data",
    "processed",
    "BHB_biozones.csv"
  )
)
