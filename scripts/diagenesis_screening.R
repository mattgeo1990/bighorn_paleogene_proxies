spar_altered_horizon_summary <- read_csv(
  here(
    "data",
    "processed",
    "spar_altered_horizon_summary.csv"
  )
)

BHB_multiproxy_final <- read_csv(
  here("data",
       "processed",
       "BHB_multiproxy_final.csv")
  )

# ---- Reconstruct d18Owater for spar / altered carbonate samples ----

names(spar_altered_horizon_summary)


spar_d18Ow_recon <- spar_altered_horizon_summary %>%
  filter(
    !is.na(mean_T47_C),
    !is.na(mean_d18Ocarb_vsmow)
  ) %>%
  mutate(
    alpha_calcite_water =
      exp(
        (
          18.03 * (1000 / (mean_T47_C + 273.15))
          - 32.42
        ) / 1000
      ),
    
    d18Ow_vsmow =
      ((mean_d18Ocarb_vsmow + 1000) / alpha_calcite_water) - 1000
  )


# ---- Quick look: micrite vs spar/altered T47-d18Ow space -----

micrite_T_d18Ow <- BHB_multiproxy_final %>%
  filter(
    T_recon_source == "IPL measured T47",
    !is.na(T_recon_C),
    !is.na(d18Ow_mean_vsmow)
  )

spar_T_d18Ow <- spar_d18Ow_recon %>%
  filter(
    !is.na(mean_T47_C),
    !is.na(d18Ow_vsmow)
  )

# ---- Combine micrite and spar/altered data for plotting ----
micrite_plot <- micrite_T_d18Ow %>%
  transmute(
    MLA_horizon_id,
    T47_C = IPLD47_mean_T47_C,
    T47_se_C = IPLD47_se_T47_C,
    d18Ocarb_vsmow,
    d18Ocarb_se_vsmow = d18Ocarb_se_used_vsmow,
    d18Ow_vsmow = d18Ow_mean_vsmow,
    carbonate_type = "micrite"
  )

spar_altered_plot <- spar_T_d18Ow %>%
  mutate(
    classification_text = str_to_lower(
      paste(sample_groups, sample_types)
    )
  ) %>%
  transmute(
    MLA_horizon_id,
    T47_C = mean_T47_C,
    T47_se_C = se_T47_C,
    d18Ocarb_vsmow = mean_d18Ocarb_vsmow,
    d18Ocarb_se_vsmow = se_d18Ocarb_vsmow,
    d18Ow_vsmow,
    carbonate_type = case_when(
      str_detect(classification_text, "fracture|spar") ~ "fracture spar",
      TRUE ~ "altered"
    )
  )

T47_isotope_plot_data <- bind_rows(
  micrite_plot,
  spar_altered_plot
) %>%
  mutate(
    carbonate_type = factor(
      carbonate_type,
      levels = c("micrite", "fracture spar", "altered")
    )
  )

# ---- Carbonate d18O versus T47 ----
# Okabe-Ito colorblind-friendly palette
carbonate_colors <- c(
  "micrite" = "#222222",
  "fracture spar" = "#0072B2",
  "altered" = "#D55E00"
)

# Shapes remain distinguishable when printed in grayscale
carbonate_shapes <- c(
  "micrite" = 16,         # circle
  "fracture spar" = 17,   # triangle
  "altered" = 15          # square
)

p_d18Ocarb_T47 <- T47_isotope_plot_data %>%
  filter(!is.na(T47_C), !is.na(d18Ocarb_vsmow)) %>%
  ggplot(
    aes(
      x = T47_C,
      y = d18Ocarb_vsmow,
      color = carbonate_type,
      shape = carbonate_type
    )
  ) +
  geom_errorbarh(
    aes(
      xmin = T47_C - T47_se_C,
      xmax = T47_C + T47_se_C
    ),
    height = 0,
    linewidth = 0.45,
    alpha = 0.55,
    na.rm = TRUE,
    show.legend = FALSE
  ) +
  geom_errorbar(
    aes(
      ymin = d18Ocarb_vsmow - d18Ocarb_se_vsmow,
      ymax = d18Ocarb_vsmow + d18Ocarb_se_vsmow
    ),
    width = 0,
    linewidth = 0.45,
    alpha = 0.55,
    na.rm = TRUE,
    show.legend = FALSE
  ) +
  geom_point(
    size = 3.2,
    alpha = 0.95,
    stroke = 0.7
  ) +
  scale_color_manual(values = carbonate_colors) +
  scale_shape_manual(values = carbonate_shapes) +
  scale_x_continuous(
    expand = expansion(mult = c(0.04, 0.06))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.06, 0.06))
  ) +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = expression(delta^18 * O[carbonate] ~ "(‰ VSMOW)"),
    color = NULL,
    shape = NULL
  ) +
  guides(
    color = guide_legend(nrow = 1),
    shape = guide_legend(nrow = 1)
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 10),
    legend.key.width = unit(1.2, "lines"),
    plot.margin = margin(6, 10, 6, 6)
  )

p_d18Ocarb_T47

# ---- Simple d18Owater evolution / diagenesis model ----

ggplot(
  BHB_d18Ow_recon %>%
    filter(
      !is.na(T_recon_C),
      !is.na(d18Ocarb_vsmow)
    ),
  aes(
    x = T_recon_C,
    y = d18Ocarb_vsmow
  )
) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()

library(tidyverse)

# ---- Forward burial-fluid scenarios in T47-d18Ocarb space ----
# Starting condition for burial trajectories
initial_T_C <- 30
initial_d18Ow_vsmow <- -4.5

calc_d18Ocarb_from_water_KO97 <- function(T_C, d18Ow_vsmow) {
  alpha <- exp((18.03 * (1000 / (T_C + 273.15)) - 32.42) / 1000)
  ((d18Ow_vsmow + 1000) * alpha) - 1000
}

burial_model <- expand_grid(
  T_C = seq(initial_T_C, 120, by = 1),
  scenario = c(
    "Open meteoric fluid",
    "Moderate rock buffering",
    "Strong rock buffering",
    "Evolved basin fluid"
  )
) %>%
  mutate(
    d18Ow_vsmow = case_when(
      scenario == "Open meteoric fluid" ~
        initial_d18Ow_vsmow,
      
      scenario == "Moderate rock buffering" ~
        initial_d18Ow_vsmow + (T_C - initial_T_C) * (6 / (120 - initial_T_C)),
      
      scenario == "Strong rock buffering" ~
        initial_d18Ow_vsmow + (T_C - initial_T_C) * (10 / (120 - initial_T_C)),
      
      scenario == "Evolved basin fluid" ~
        initial_d18Ow_vsmow + (T_C - initial_T_C) * (14 / (120 - initial_T_C))
    ),
    
    d18Ocarb_vsmow = calc_d18Ocarb_from_water_KO97(
      T_C,
      d18Ow_vsmow
    )
  )


obs <- T47_isotope_plot_data %>%
  filter(
    !is.na(T47_C),
    !is.na(d18Ocarb_vsmow)
  )

ggplot() +
  geom_line(
    data = burial_model,
    aes(
      x = T_C,
      y = d18Ocarb_vsmow,
      color = scenario
    ),
    linewidth = 1
  ) +
  geom_point(
    data = obs,
    aes(
      x = T47_C,
      y = d18Ocarb_vsmow,
      shape = carbonate_type
    ),
    size = 3,
    alpha = 0.8
  ) +
  geom_errorbar(
    data = obs,
    aes(
      x = T47_C,
      ymin = d18Ocarb_vsmow - d18Ocarb_se_vsmow,
      ymax = d18Ocarb_vsmow + d18Ocarb_se_vsmow
    ),
    width = 0,
    alpha = 0.35
  ) +
  geom_errorbarh(
    data = obs,
    aes(
      xmin = T47_C - T47_se_C,
      xmax = T47_C + T47_se_C,
      y = d18Ocarb_vsmow
    ),
    height = 0,
    alpha = 0.35
  ) +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = expression(delta^18 * O[carb] ~ "(‰ VSMOW)"),
    color = expression(delta^18 * O[water] ~ "scenario"),
    shape = "Carbonate type"
  ) +
  theme_classic(base_size = 14)
# Burial-fluid evolution scenarios used to illustrate how carbonate δ18O
# may respond to recrystallization under different diagenetic conditions.
#
# Open meteoric fluid:
#   δ18Ow remains constant (-8‰) during burial, representing an open system
#   continuously flushed by meteoric water.
#
# Moderate rock buffering:
#   δ18Ow becomes progressively enriched through water-rock interaction,
#   producing modest isotopic evolution of pore fluids.
#
# Strong rock buffering:
#   More extensive water-rock exchange drives larger increases in δ18Ow
#   as burial temperature increases.
#
# Evolved basin fluid:
#   Represents highly evolved burial fluids with strongly enriched δ18Ow,
#   approximating near-closed-system conditions and prolonged fluid-rock interaction.
#
# These trajectories are illustrative end-member models intended to visualize
# possible diagenetic pathways in T–δ18O space and are not calibrated burial
# history simulations.

# Starting condition for burial trajectories:
# The initial fluid composition is set to δ18Ow = -11‰ VSMOW, based on the
# mean LPEE meteoric-water estimate from Campbell et al. (2024).
# The initial crystallization temperature is set to 29°C, based on the mean
# Snell et al. (2013) carbonate Δ47 temperature after applying a -5°C
# radiation correction to approximate soil-temperature conditions.
#
# Together, these values define a plausible primary pedogenic carbonate
# starting point for exploring how carbonate δ18O would evolve during
# burial under different fluid-buffering scenarios.



# ---- Monte Carlo burial-fluid trajectories in T47-d18Ocarb space ----

set.seed(123)

n_paths <- 500

calc_d18Ocarb_from_water_KO97 <- function(T_C, d18Ow_vsmow) {
  alpha <- exp((18.03 * (1000 / (T_C + 273.15)) - 32.42) / 1000)
  ((d18Ow_vsmow + 1000) * alpha) - 1000
}

trajectory_params <- tibble(
  path_id = 1:n_paths,
  
  # plausible primary / early soil carbonate starting conditions
  initial_T_C = runif(n_paths, 15, 40),
  initial_d18Ow_vsmow = runif(n_paths, -8, -2),
  
  # maximum burial / recrystallization temperature
  max_T_C = runif(n_paths, 60, 120),
  
  # total enrichment in fluid d18O caused by water-rock interaction
  # 0 = open meteoric; larger values = stronger rock buffering / evolved fluids
  total_d18Ow_shift = runif(n_paths, 0, 14)
)

burial_mc <- trajectory_params %>%
  rowwise() %>%
  mutate(
    data = list(
      tibble(
        T_C = seq(initial_T_C, max_T_C, length.out = 100)
      ) %>%
        mutate(
          progress = (T_C - min(T_C)) / (max(T_C) - min(T_C)),
          
          d18Ow_vsmow =
            initial_d18Ow_vsmow +
            total_d18Ow_shift * progress,
          
          d18Ocarb_vsmow =
            calc_d18Ocarb_from_water_KO97(
              T_C = T_C,
              d18Ow_vsmow = d18Ow_vsmow
            )
        )
    )
  ) %>%
  ungroup() %>%
  select(path_id, initial_T_C, initial_d18Ow_vsmow, max_T_C, total_d18Ow_shift, data) %>%
  unnest(data)

obs <- T47_isotope_plot_data %>%
  filter(
    !is.na(T47_C),
    !is.na(d18Ocarb_vsmow)
  )

ggplot() +
  geom_line(
    data = burial_mc,
    aes(
      x = T_C,
      y = d18Ocarb_vsmow,
      group = path_id,
      color = total_d18Ow_shift
    ),
    alpha = 0.18,
    linewidth = 0.45
  ) +
  geom_errorbar(
    data = obs,
    aes(
      x = T47_C,
      ymin = d18Ocarb_vsmow - d18Ocarb_se_vsmow,
      ymax = d18Ocarb_vsmow + d18Ocarb_se_vsmow
    ),
    width = 0,
    alpha = 0.25
  ) +
  geom_errorbarh(
    data = obs,
    aes(
      xmin = T47_C - T47_se_C,
      xmax = T47_C + T47_se_C,
      y = d18Ocarb_vsmow
    ),
    height = 0,
    alpha = 0.25
  ) +
  geom_point(
    data = obs,
    aes(
      x = T47_C,
      y = d18Ocarb_vsmow,
      shape = carbonate_type
    ),
    size = 3,
    alpha = 0.85
  ) +
  scale_color_viridis_c(
    name = expression("Total " * delta^18 * O[water] * " shift (‰)")
  ) +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = expression(delta^18 * O[carb] ~ "(‰ VSMOW)"),
    shape = "Carbonate type"
  ) +
  theme_classic(base_size = 14)
# ---- Interactive forward model plot: burial evolution of d18Ocarb ----
library(plotly)

p_burial_model <- ggplot() +
  
  geom_line(
    data = burial_model,
    aes(
      x = T_C,
      y = d18Ocarb_vsmow,
      color = scenario,
      group = scenario
    ),
    linewidth = 1
  ) +
  
  geom_point(
    data = obs,
    aes(
      x = T47_C,
      y = d18Ocarb_vsmow,
      shape = carbonate_type,
      text = paste0(
        "Horizon: ", MLA_horizon_id,
        "<br>Type: ", carbonate_type,
        "<br>T47: ", round(T47_C, 1), " °C",
        "<br>d18Ocarb: ", round(d18Ocarb_vsmow, 2), "‰",
        "<br>d18Ow: ", round(d18Ow_vsmow, 2), "‰"
      )
    ),
    size = 3,
    alpha = 0.8
  ) +
  
  labs(
    x = "Δ47 temperature (°C)",
    y = "δ18Ocarb (‰ VSMOW)",
    color = "Scenario",
    shape = "Carbonate type"
  ) +
  
  theme_classic(base_size = 14)

ggplotly(
  p_burial_model,
  tooltip = "text"
)


# ---- Interactive: micrite vs spar/altered T47-d18Ow space ----

micrite_T_d18Ow <- BHB_d18Ow_recon %>%
  filter(
    T_recon_source == "IPL measured T47",
    !is.na(T_recon_C),
    !is.na(d18Ow_mean_vsmow)
  ) %>%
  mutate(
    carbonate_group = "Micrite",
    hover_text = paste0(
      "Horizon: ", MLA_horizon_id,
      "<br>Strat height: ", strat_height_m, " m",
      "<br>T47: ", round(T_recon_C, 1), " °C",
      "<br>d18Ow: ", round(d18Ow_mean_vsmow, 2), "‰ VSMOW",
      "<br>d18Ocarb source: ", d18Ocarb_source
    )
  )

spar_T_d18Ow <- spar_d18Ow_recon %>%
  filter(
    !is.na(mean_T47_C),
    !is.na(d18Ow_vsmow)
  ) %>%
  mutate(
    carbonate_group = sample_groups,
    hover_text = paste0(
      "Horizon: ", MLA_horizon_id,
      "<br>Strat height: ", strat_height_m, " m",
      "<br>Type: ", sample_groups,
      "<br>T47: ", round(mean_T47_C, 1), " °C",
      "<br>d18Ow: ", round(d18Ow_vsmow, 2), "‰ VSMOW",
      "<br>d18Ocarb: ", round(mean_d18Ocarb_vsmow, 2), "‰ VSMOW"
    )
  )

p_micrite_spar_d18Ow <- ggplot() +
  geom_point(
    data = micrite_T_d18Ow,
    aes(
      x = T_recon_C,
      y = d18Ow_mean_vsmow,
      text = hover_text
    ),
    size = 2.5,
    alpha = 0.75
  ) +
  geom_point(
    data = spar_T_d18Ow,
    aes(
      x = mean_T47_C,
      y = d18Ow_vsmow,
      shape = carbonate_group,
      text = hover_text
    ),
    size = 3.5,
    alpha = 0.95
  ) +
  labs(
    x = "T47 temperature (°C)",
    y = "d18Owater (‰ VSMOW)",
    shape = "Spar / altered type",
    title = "Micrite vs. spar/altered carbonate T47-d18Owater space"
  ) +
  theme_classic()

ggplotly(
  p_micrite_spar_d18Ow,
  tooltip = "text"
)
# model d18Ocarb from prior d18Owater and Temp --------
# Kim & O'Neil (1997) style fractionation
calc_d18Ocarb_from_water_KO97 <- function(T_C, d18Ow_vsmow) {
  alpha <- exp((18.03 * (1000 / (T_C + 273.15)) - 32.42) / 1000)
  ((d18Ow_vsmow + 1000) * alpha) - 1000
}

# Temperature range
temps <- tibble(
  T_C = seq(0, 120, by = 1)
) %>%
  mutate(
    d18Ow_vsmow =
      ((22.5 + 1000) /
         exp((18.03 * (1000 / (T_C + 273.15)) - 32.42) / 1000)
      ) - 1000
  )

ggplot(
  temps,
  aes(
    x = T_C,
    y = d18Ow_vsmow
  )
) +
  geom_line(linewidth = 1) +
  geom_hline(
    yintercept = -5,
    linetype = 2
  ) +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = expression(delta^18 * O[water] ~ "(‰ VSMOW)"),
    title = expression(
      delta^18 * O[carb] == 22.5 * "\u2030 VSMOW"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 50),
    breaks = seq(0, 50, by = 5)
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_line(color = "grey92")
  )

# Spar mixing sensitivity model ---------------------

library(tidyverse)

# Fixed micrite temperature
T_micrite <- 30

# Plausible spar temperatures
spar_temps <- c(45, 55, 85)

mixing_df <- expand_grid(
  T_spar = spar_temps,
  percent_spar = seq(0, 100, by = 1)
) %>%
  mutate(
    f_spar = percent_spar / 100,
    bulk_T = f_spar * T_spar +
      (1 - f_spar) * T_micrite
  )
ggplot(
  mixing_df,
  aes(
    x = percent_spar,
    y = bulk_T,
    color = factor(T_spar)
  )
) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 30,
    ymax = 60,
    alpha = 0.15
  ) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 10)
  ) +
  scale_y_continuous(
    breaks = seq(20, 100, by = 10)
  ) +
  labs(
    x = "Spar in analyzed powder (%)",
    y = expression("Bulk " * Delta[47] * " temperature (" * degree * "C)"),
    color = "Spar T (°C)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

# Cull data suspected of alteration -----

nuclear <- c(
  #tier 1 (heavy alteration):
  "PK95-SC-4",
  "PK95-SC-279",
  "PK95-SC-242",
  "PK95-SC-27",
  # tier 2 (moderate alteration): 
  "PK95-SC-176",
  "PK95-SC-246",
  "PK95-SC-118up",
  # tier 3 (possibly altered):
  "PK95-SC-187",
  "PK95-SC-160",
  "PK95-SC-6"
)

moderate <- c(
  #tier 1 (heavy alteration):
  "PK95-SC-4",
  "PK95-SC-279",
  "PK95-SC-242",
  "PK95-SC-27",
  # tier 2 (moderate alteration): 
  "PK95-SC-176",
  "PK95-SC-246",
  "PK95-SC-118up"
)

conservative <- c(
  #tier 1 (heavy alteration):
  "PK95-SC-4",
  "PK95-SC-279",
  "PK95-SC-242",
  "PK95-SC-27"
)


# Choose cull scheme
altered_cull_data <- nuclear 

BHB_multiproxy_with_temperature %>%
  filter(MLA_horizon_id %in% altered_cull_data) %>%
  select(
    MLA_horizon_id,
    strat_height_m,
    T_measured_C,
    T_model_C,
    T_model_lower95_C,
    T_model_upper95_C
  ) %>%
  arrange(strat_height_m)



