# 08_age_domain_plots.R
# Purpose: Generate age-domain plots for BHB multiproxy manuscript/meeting figures

library(tidyverse)
library(here)
library(patchwork)
library(ggstance)

# ---- Shared age plot settings ----

age_limits <- c(59, 53.5)
age_breaks <- seq(53.5, 59, by = 0.1)

petm_age_old <- 55.935
petm_age_young <- 55.75

add_petm_age <- function(fill = "red", alpha = 0.25) {
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = petm_age_young, ymax = petm_age_old,
    fill = fill,
    alpha = alpha
  )
}

shared_age_scale <- function() {
  scale_y_reverse(
    limits = age_limits,
    breaks = age_breaks,
    expand = expansion(mult = c(0.01, 0.02))
  )
}

temp_x_scale <- function() {
  scale_x_continuous(
    breaks = seq(10, 75, by = 10),
    limits = c(10, 75),
    expand = expansion(mult = c(0.02, 0.04))
  )
}

theme_age_panel <- theme_classic(base_size = 10) +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.position = "none",
    plot.margin = margin(4, 4, 4, 4)
  )

# ---- Load data ----

BHB_multiproxy_summary <- read_csv(here("data", "processed", "BHB_multiproxy_summary.csv"))
BHB_d18Ow <- read_csv(here("data", "processed", "BHB_d18Ow_reconstruction.csv"))
BHB_multiproxy_final <- read_csv(here("data", "processed", "BHB_D17Orsw_reconstruction.csv"))
BHB_temperature_model <- read_csv(here("data", "processed", "BHB_temperature_model.csv"))
temp_obs <- read_csv(here("data", "processed", "BHB_temperature_observations.csv"))

Kelson_Tornillo_D47 <- read_csv(
  here("data", "processed", "Kelson_Tornillo_D47_strat_age_filtered.csv")
)


Harper2024_CO2_SST <- read_csv(
  here("data", "processed", "Harper2024_CO2_SST_processed.csv")
)

biozones <- read_csv(
  here(
    "data",
    "processed",
    "BHB_biozones.csv"
  )
)

age_priors <- read.csv(
  here("data", "raw", "PCB-CFB_age_priors.csv")
)

age_priors_clean <- age_priors %>%
  mutate(
    est_Depth_m_PCB_outcrop = as.numeric(est_Depth_m_PCB_outcrop),
    Age_Ma_best_estimate = as.numeric(Age_Ma_best_estimate)
  ) %>%
  filter(
    !is.na(est_Depth_m_PCB_outcrop),
    !is.na(Age_Ma_best_estimate)
  ) %>%
  arrange(est_Depth_m_PCB_outcrop)

biozones_age <- biozones %>%
  mutate(
    Age_ymin_Ma = approx(
      x = age_priors_clean$est_Depth_m_PCB_outcrop,
      y = age_priors_clean$Age_Ma_best_estimate,
      xout = ymin,
      rule = 1
    )$y,
    Age_ymax_Ma = approx(
      x = age_priors_clean$est_Depth_m_PCB_outcrop,
      y = age_priors_clean$Age_Ma_best_estimate,
      xout = ymax,
      rule = 1
    )$y
  )

biozones_age

# ---- Add age to temperature model and observations ----

age_lookup <- BHB_multiproxy_summary %>%
  select(MLA_horizon_id, strat_height_m, Age_Ma) %>%
  filter(!is.na(strat_height_m), !is.na(Age_Ma))

BHB_temperature_model_age <- BHB_temperature_model %>%
  left_join(age_lookup, by = c("MLA_horizon_id", "strat_height_m")) %>%
  filter(!is.na(Age_Ma))

temp_obs_age <- temp_obs %>%
  left_join(age_lookup, by = c("MLA_horizon_id", "strat_height_m")) %>%
  filter(!is.na(T_C), !is.na(Age_Ma))

# ---- Temperature age plot ----

p_T47_age_clean <- ggplot() +
  geom_ribbon(
    data = BHB_temperature_model_age,
    aes(
      xmin = T_model_lower95_C,
      xmax = T_model_upper95_C,
      y = Age_Ma
    ),
    fill = "grey70",
    alpha = 0.45
  ) +
  add_petm_age() +
  geom_path(
    data = BHB_temperature_model_age,
    aes(x = T_model_C, y = Age_Ma),
    linewidth = 1
  ) +
  geom_errorbarh(
    data = temp_obs_age,
    aes(
      xmin = T_C - T_se_C,
      xmax = T_C + T_se_C,
      y = Age_Ma
    ),
    height = 0,
    linewidth = 0.35,
    alpha = 0.45
  ) +
  geom_point(
    data = temp_obs_age,
    aes(x = T_C, y = Age_Ma),
    size = 1.8,
    alpha = 0.75
  ) +
  scale_x_continuous(
    breaks = seq(10, 75, by = 10),
    limits = c(10, 75),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  shared_age_scale() +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = "Age (Ma)"
  ) +
  theme_age_panel

library(tidyverse)
library(ggstance)


# ---- d13C carbonate age plot ----

d13C_obs_age <- BHB_multiproxy_summary %>%
  select(
    MLA_horizon_id,
    strat_height_m,
    Age_Ma,
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
  filter(!is.na(d13Ccarb_vpdb), !is.na(Age_Ma))

p_d13C_age_clean <- ggplot() +
  add_petm_age() +
  geom_errorbarh(
    data = d13C_obs_age,
    aes(
      xmin = d13Ccarb_vpdb - d13C_se_vpdb,
      xmax = d13Ccarb_vpdb + d13C_se_vpdb,
      y = Age_Ma
    ),
    height = 0,
    linewidth = 0.35,
    alpha = 0.45,
    na.rm = TRUE
  ) +
  geom_point(
    data = d13C_obs_age,
    aes(x = d13Ccarb_vpdb, y = Age_Ma),
    size = 1.8,
    alpha = 0.75
  ) +
  scale_x_continuous(
    breaks = seq(-16, -4, by = 2),
    limits = c(-16, -4),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  shared_age_scale() +
  labs(
    x = expression(delta^13 * C[carb] ~ "(‰ VPDB)"),
    y = "Age (Ma)"
  ) +
  theme_age_panel

# ---- d18O soil-water age plot ----

d18Ow_age <- BHB_d18Ow %>%
  filter(!is.na(d18Ow_mean_vsmow), !is.na(Age_Ma))

p_d18Ow_age_clean <- ggplot() +
  add_petm_age() +
  geom_errorbarh(
    data = d18Ow_age,
    aes(
      xmin = d18Ow_lower95_vsmow,
      xmax = d18Ow_upper95_vsmow,
      y = Age_Ma
    ),
    height = 0,
    linewidth = 0.30,
    alpha = 0.35,
    na.rm = TRUE
  ) +
  geom_point(
    data = d18Ow_age,
    aes(x = d18Ow_mean_vsmow, y = Age_Ma),
    size = 1.6,
    alpha = 0.70
  ) +
  scale_x_continuous(
    breaks = seq(-10, 2, by = 2),
    limits = c(-10, 2),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  shared_age_scale() +
  labs(
    x = expression(delta^18 * O[soil-water] ~ "(‰ VSMOW)"),
    y = "Age (Ma)"
  ) +
  theme_age_panel

# ---- D17O soil-water age plot ----

D17Orsw_age <- BHB_multiproxy_final %>%
  filter(!is.na(D17Orsw_mean_permeg), !is.na(Age_Ma))

p_D17Orsw_age_clean <- ggplot() +
  add_petm_age() +
  geom_errorbarh(
    data = D17Orsw_age,
    aes(
      xmin = D17Orsw_mean_permeg - IPL17O_sd_Dp17Ocarb_adj,
      xmax = D17Orsw_mean_permeg + IPL17O_sd_Dp17Ocarb_adj,
      y = Age_Ma
    ),
    height = 0,
    linewidth = 0.30,
    alpha = 0.35,
    na.rm = TRUE
  ) +
  geom_point(
    data = D17Orsw_age,
    aes(x = D17Orsw_mean_permeg, y = Age_Ma),
    size = 1.6,
    alpha = 0.70
  ) +
  scale_x_continuous(
    breaks = seq(-100, 50, by = 25),
    limits = c(-100, 50),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  shared_age_scale() +
  labs(
    x = expression(Delta*minute^17*O[soil-water] ~ "(per meg)"),
    y = "Age (Ma)"
  ) +
  theme_age_panel

# ---- Harper reference records ----
# ---- Smooth Harper reference records like BHB temperature model ----

Harper_plot <- Harper2024_CO2_SST %>%
  filter(
    Age_Ma >= min(age_limits),
    Age_Ma <= max(age_limits)
  ) %>%
  arrange(Age_Ma) %>%
  mutate(
    SST_anchor_C = zoo::rollapply(
      Harper2024_mean_SST_C,
      width = 3,
      FUN = mean,
      fill = NA,
      align = "center",
      partial = TRUE
    ),
    CO2_anchor_ppm = zoo::rollapply(
      Harper2024_mean_CO2_ppm,
      width = 3,
      FUN = mean,
      fill = NA,
      align = "center",
      partial = TRUE
    )
  )

Harper_SST_spline <- smooth.spline(
  x = Harper_plot$Age_Ma,
  y = Harper_plot$SST_anchor_C,
  spar = 0.35
)

Harper_CO2_spline <- smooth.spline(
  x = Harper_plot$Age_Ma,
  y = Harper_plot$CO2_anchor_ppm,
  spar = 0.35
)

Harper_plot <- Harper_plot %>%
  mutate(
    Harper2024_SST_smooth_C = as.numeric(
      predict(Harper_SST_spline, x = Age_Ma)$y
    ),
    Harper2024_CO2_smooth_ppm = as.numeric(
      predict(Harper_CO2_spline, x = Age_Ma)$y
    )
  )

p_Harper_SST_age <- ggplot(Harper_plot) +
  add_petm_age(fill = "red", alpha = 0.25) +
  geom_point(
    aes(x = Harper2024_mean_SST_C, y = Age_Ma),
    size = 1.0,
    alpha = 0.18
  ) +
  geom_ribbon(
    aes(
      xmin = Harper2024_SST_lower95_C,
      xmax = Harper2024_SST_upper95_C,
      y = Age_Ma
    ),
    fill = "grey70",
    alpha = 0.45
  ) +
  geom_path(
    aes(x = Harper2024_SST_smooth_C, y = Age_Ma),
    linewidth = 1
  ) +
  scale_x_continuous(
    breaks = seq(20, 40, by = 5),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  shared_age_scale() +
  labs(
    x = expression("Pacific SST (" * degree * "C)"),
    y = "Age (Ma)"
  ) +
  theme_age_panel

# Adjust CO2 column names here if your processed file uses different names
p_Harper_CO2_age <- ggplot(Harper_plot) +
  add_petm_age(fill = "red", alpha = 0.25) +
  geom_point(
    aes(x = Harper2024_mean_CO2_ppm, y = Age_Ma),
    size = 1.0,
    alpha = 0.18
  ) +
  geom_ribbon(
    aes(
      xmin = Harper2024_CO2_lower95_ppm,
      xmax = Harper2024_CO2_upper95_ppm,
      y = Age_Ma
    ),
    fill = "grey70",
    alpha = 0.45
  ) +
  geom_path(
    aes(x = Harper2024_CO2_smooth_ppm, y = Age_Ma),
    linewidth = 1
  ) +
  scale_x_continuous(
    breaks = seq(0, 3000, by = 500),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  shared_age_scale() +
  labs(
    x = expression(CO[2] ~ "(ppm)"),
    y = "Age (Ma)"
  ) +
  theme_age_panel

# ---- Kelson Tornillo primary D47 temperatures ----

# ---- Kelson Tornillo primary D47 temperatures ----

Kelson_Tornillo_primary_age <- Kelson_Tornillo_D47 %>%
  filter(
    petro_class == "Micrite",
    !is.na(Age_Ma),
    !is.na(T47_C)
  ) %>%
  arrange(Age_Ma) %>%
  mutate(
    T47_anchor_C = zoo::rollapply(
      T47_C,
      width = 3,
      FUN = mean,
      fill = NA,
      align = "center",
      partial = TRUE
    )
  )

Kelson_T47_spline <- smooth.spline(
  x = Kelson_Tornillo_primary_age$Age_Ma,
  y = Kelson_Tornillo_primary_age$T47_anchor_C,
  spar = 0.35
)

Kelson_Tornillo_primary_age <- Kelson_Tornillo_primary_age %>%
  mutate(
    T47_smooth_C = as.numeric(
      predict(Kelson_T47_spline, x = Age_Ma)$y
    )
  )

p_Kelson_Tornillo_T47_age <- ggplot() +
  add_petm_age(fill = "red", alpha = 0.25) +
  geom_errorbarh(
    data = Kelson_Tornillo_primary_age,
    aes(
      xmin = T47_C - T47_se_C,
      xmax = T47_C + T47_se_C,
      y = Age_Ma
    ),
    height = 0,
    linewidth = 0.35,
    alpha = 0.45,
    na.rm = TRUE
  ) +
  geom_point(
    data = Kelson_Tornillo_primary_age,
    aes(x = T47_C, y = Age_Ma),
    size = 1.8,
    alpha = 0.75
  ) +
  geom_path(
    data = Kelson_Tornillo_primary_age,
    aes(x = T47_smooth_C, y = Age_Ma),
    linewidth = 1,
    na.rm = TRUE
  ) +
  temp_x_scale() +
  shared_age_scale() +
  labs(
    x = expression("Tornillo " * Delta[47] * " temperature (" * degree * "C)"),
    y = "Age (Ma)"
  ) +
  theme_age_panel
# ---- Panel figures ----

p_T47_age_panel <- p_T47_age_clean +
  labs(tag = "A")

p_d13C_age_panel <- p_d13C_age_clean +
  labs(y = NULL, tag = "B") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_d18Ow_age_panel <- p_d18Ow_age_clean +
  labs(y = NULL, tag = "C") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_D17Orsw_age_panel <- p_D17Orsw_age_clean +
  labs(y = NULL, tag = "D") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_Harper_SST_panel <- p_Harper_SST_age +
  labs(y = NULL, tag = "E") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_Harper_CO2_panel <- p_Harper_CO2_age +
  labs(y = NULL, tag = "F") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_age_panel_BHB_only <- 
  p_T47_age_panel + p_d13C_age_panel + p_d18Ow_age_panel + p_D17Orsw_age_panel +
  plot_layout(nrow = 1, ncol = 4, widths = c(1, 1, 1, 1)) +
  plot_annotation(
    theme = theme(plot.tag = element_text(face = "bold", size = 12))
  )


p_age_panel_temp_comparison <- 
  p_T47_age_panel +  p_Harper_SST_panel + p_Harper_CO2_panel + 
  plot_layout(nrow = 1, ncol = 3, widths = c(1, 1, 1, 1, 1)) +
  plot_annotation(
    theme = theme(plot.tag = element_text(face = "bold", size = 12))
  )

# ---- Multi-panel temperature comparison: Harper SST, Tornillo, BHB ----

p_Harper_SST_panel <- p_Harper_SST_age +
  labs(y = "Age (Ma)", tag = "A")

p_Kelson_Tornillo_panel <- p_Kelson_Tornillo_T47_age +
  labs(y = NULL, tag = "B") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_BHB_T47_panel <- p_T47_age_clean +
  labs(y = NULL, tag = "C") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_age_panel_temp_comparison <- 
  p_Harper_SST_panel +
  p_Kelson_Tornillo_panel +
  p_BHB_T47_panel +
  plot_layout(nrow = 1, ncol = 3, widths = c(1, 1, 1)) +
  plot_annotation(
    theme = theme(plot.tag = element_text(face = "bold", size = 12))
  )

p_age_panel_temp_comparison

ggsave(
  here("figures", "age_domain", "Harper_Tornillo_BHB_temperature_age_panel.png"),
  p_age_panel_temp_comparison,
  width = 12,
  height = 6,
  dpi = 600
)
p_age_panel_BHB_only
p_age_panel_temp_comparison

#-- 17O through time ------

p_age_panel_17O <- 
  p_BHB_T47_panel +
  p_d18Ow_age_clean +
  p_D17Orsw_age_clean +
  plot_layout(nrow = 1, ncol = 3, widths = c(1, 1, 1)) +
  plot_annotation(
    theme = theme(plot.tag = element_text(face = "bold", size = 12))
  )

p_age_panel_17O


# ---- Save outputs ----


dir.create(here("figures", "age_domain"), recursive = TRUE, showWarnings = FALSE)

ggsave(
  here("figures", "age_domain", "BHB_age_panel_T47_d13C_d18Ow_D17Orsw.png"),
  p_age_panel_BHB_only,
  width = 12,
  height = 6,
  dpi = 600
)

ggsave(
  here("figures", "age_domain", "BHB_age_panel_CO2_SST_T47.png"),
  p_age_panel_temp_comparison,
  width = 12,
  height = 6,
  dpi = 600
)
