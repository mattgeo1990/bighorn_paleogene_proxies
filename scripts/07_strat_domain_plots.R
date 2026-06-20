# 07_strat_domain_plots.R
# Purpose: Generate stratigraphic-domain plots for BHB multiproxy manuscript/meeting figures

library(tidyverse)
library(here)
library(patchwork)

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
    breaks = seq(-100, 50, by = 10),
    limits = c(-70, 10),
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
