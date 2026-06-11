# 04_soilwater_d18O_reconstructions.R
# Purpose: Reconstruct soil-water d18O from carbonate d18O and temperature.

# ---- Load packages ----
library(tidyverse)
library(here)

# ---- Settings ----
n_mc <- 10000
set.seed(123)

# ---- Load data with modeled temperatures ----
BHB <- read_csv(
  here("data", "processed", "BHB_multiproxy_with_temperature.csv")
)

# ---- Add placeholder columns if future U-M data are not present yet ----
if (!"IPL_NuDog_d18Ocarb_VSMOW" %in% names(BHB)) {
  BHB <- BHB %>%
    mutate(IPL_NuDog_d18Ocarb_VSMOW = NA_real_)
}

# ---- Choose carbonate d18O source ----
# Priority:
#   Snell samples: use Snell d18Ocarb
#   PB samples: use Bowen d18Ocarb
#   PK95-SC samples: use Koch; if missing, use IPL_NuDog; if missing, use IPL17O

BHB_water_inputs <- BHB %>%
  mutate(
    is_snell = !is.na(Snell_sample_id) & !is.na(Snell_mean_d18Ocarb_vsmow),
    
    d18Ocarb_vsmow = case_when(
      is_snell ~ Snell_mean_d18Ocarb_vsmow,
      str_detect(MLA_horizon_id, "^PB") ~ Bowen_mean_d18Ocarb_vsmow,
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(Koch_mean_d18Ocarb_vsmow) ~ Koch_mean_d18Ocarb_vsmow,
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(IPL_NuDog_d18Ocarb_VSMOW) ~ IPL_NuDog_d18Ocarb_VSMOW,
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(IPL17O_mean_dp18Ocarb) ~ IPL17O_mean_dp18Ocarb,
      TRUE ~ NA_real_
    ),
    
    d18Ocarb_se_vsmow = case_when(
      is_snell ~ NA_real_,
      str_detect(MLA_horizon_id, "^PB") ~ Bowen_se_d18Ocarb_vsmow,
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(Koch_mean_d18Ocarb_vsmow) ~ Koch_se_d18Ocarb_vsmow,
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(IPL_NuDog_d18Ocarb_VSMOW) ~ NA_real_,
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(IPL17O_mean_dp18Ocarb) ~ IPL17O_se_dp18Ocarb,
      TRUE ~ NA_real_
    ),
    
    d18Ocarb_source = case_when(
      is_snell ~ "Snell",
      str_detect(MLA_horizon_id, "^PB") & !is.na(Bowen_mean_d18Ocarb_vsmow) ~ "Bowen",
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(Koch_mean_d18Ocarb_vsmow) ~ "Koch",
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(IPL_NuDog_d18Ocarb_VSMOW) ~ "IPL_NuDog",
      str_detect(MLA_horizon_id, "^PK95-SC") & !is.na(IPL17O_mean_dp18Ocarb) ~ "IPL17O",
      TRUE ~ NA_character_
    ),
    
    d18Ocarb_missing_flag = case_when(
      str_detect(MLA_horizon_id, "^PB") & is.na(Bowen_mean_d18Ocarb_vsmow) ~ "PB missing Bowen d18Ocarb",
      str_detect(MLA_horizon_id, "^PK95-SC") &
        is.na(Koch_mean_d18Ocarb_vsmow) &
        is.na(IPL_NuDog_d18Ocarb_VSMOW) &
        is.na(IPL17O_mean_dp18Ocarb) ~ "PK95-SC missing Koch, IPL_NuDog, and IPL17O d18Ocarb",
      is.na(d18Ocarb_vsmow) ~ "No selected d18Ocarb",
      TRUE ~ NA_character_
    )
  )

# ---- Choose temperature source ----
# Prioritize IPL measured T47. Otherwise use modeled/interpolated temperature.

BHB_water_inputs <- BHB_water_inputs %>%
  mutate(
    T_recon_C = case_when(
      !is.na(IPLD47_mean_T47_C) ~ IPLD47_mean_T47_C,
      TRUE ~ T_model_C
    ),
    
    T_recon_se_C = case_when(
      !is.na(IPLD47_mean_T47_C) & !is.na(IPLD47_se_T47_C) ~ IPLD47_se_T47_C,
      !is.na(IPLD47_mean_T47_C) & is.na(IPLD47_se_T47_C) ~ T_model_se_C,
      TRUE ~ T_model_se_C
    ),
    
    T_recon_source = case_when(
      !is.na(IPLD47_mean_T47_C) ~ "IPL measured T47",
      !is.na(T_model_C) ~ "modeled/interpolated T",
      TRUE ~ NA_character_
    )
  )

# ---- Fill missing carbonate SE values ----
default_d18Ocarb_se <- 0.15

BHB_water_inputs <- BHB_water_inputs %>%
  mutate(
    d18Ocarb_se_used_vsmow = case_when(
      !is.na(d18Ocarb_se_vsmow) & d18Ocarb_se_vsmow > 0 ~ d18Ocarb_se_vsmow,
      !is.na(d18Ocarb_vsmow) ~ default_d18Ocarb_se,
      TRUE ~ NA_real_
    )
  )

# ---- Calcite-water fractionation: Kim & O'Neil style ----
calcite_water_alpha18 <- function(T_C) {
  T_K <- T_C + 273.15
  exp((18.03 * (1000 / T_K) - 32.42) / 1000)
}

# ---- Convert carbonate d18O and temperature to water d18O ----
reconstruct_d18Ow <- function(d18Ocarb_vsmow, T_C) {
  alpha18 <- calcite_water_alpha18(T_C)
  ((1000 + d18Ocarb_vsmow) / alpha18) - 1000
}

# ---- Monte Carlo reconstruction for one row ----
mc_d18Ow_one <- function(d18Ocarb_mean, d18Ocarb_se, T_mean, T_se, n = 10000) {
  d18Ocarb_sim <- rnorm(n, mean = d18Ocarb_mean, sd = d18Ocarb_se)
  T_sim <- rnorm(n, mean = T_mean, sd = T_se)
  
  d18Ow_sim <- reconstruct_d18Ow(
    d18Ocarb_vsmow = d18Ocarb_sim,
    T_C = T_sim
  )
  
  tibble(
    d18Ow_mean_vsmow = mean(d18Ow_sim, na.rm = TRUE),
    d18Ow_median_vsmow = median(d18Ow_sim, na.rm = TRUE),
    d18Ow_sd_vsmow = sd(d18Ow_sim, na.rm = TRUE),
    d18Ow_lower95_vsmow = quantile(d18Ow_sim, 0.025, na.rm = TRUE),
    d18Ow_upper95_vsmow = quantile(d18Ow_sim, 0.975, na.rm = TRUE)
  )
}

# ---- Run reconstruction ----
BHB_d18Ow_recon <- BHB_water_inputs %>%
  filter(
    !is.na(d18Ocarb_vsmow),
    !is.na(d18Ocarb_se_used_vsmow),
    !is.na(T_recon_C),
    !is.na(T_recon_se_C)
  ) %>%
  mutate(
    mc = pmap(
      list(
        d18Ocarb_vsmow,
        d18Ocarb_se_used_vsmow,
        T_recon_C,
        T_recon_se_C
      ),
      ~ mc_d18Ow_one(
        d18Ocarb_mean = ..1,
        d18Ocarb_se = ..2,
        T_mean = ..3,
        T_se = ..4,
        n = n_mc
      )
    )
  ) %>%
  unnest(mc)

# ---- QC outputs ----
missing_d18Ocarb <- BHB_water_inputs %>%
  filter(!is.na(d18Ocarb_missing_flag)) %>%
  select(
    MLA_horizon_id,
    strat_height_m,
    d18Ocarb_missing_flag
  )

BHB_water_inputs %>%
  count(d18Ocarb_source)

missing_d18Ocarb

# Quick look plots ------

ggplot(
  BHB_d18Ow_recon,
  aes(x = d18Ow_mean_vsmow, y = strat_height_m)
) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 1500,
    ymax = 1540,
    fill = "grey70",
    alpha = 0.25
  ) +
  geom_errorbarh(
    aes(
      xmin = d18Ow_lower95_vsmow,
      xmax = d18Ow_upper95_vsmow
    ),
    height = 0,
    alpha = 0.5
  ) +
  geom_point(
    aes(shape = d18Ocarb_source, color = T_recon_source),
    size = 2.5,
    alpha = 0.85
  ) +
  annotate(
    "text",
    x = Inf,
    y = 1520,
    label = "PETM",
    hjust = 1.1,
    fontface = "bold"
  ) +
  scale_y_continuous(
    breaks = seq(500, 2300, by = 100)
  ) +
  scale_x_continuous(
    breaks = seq(-12, 2, by = 2)
  ) +
  labs(
    x = expression(delta^18 * O[sw] ~ "(‰ VSMOW)"),
    y = "Stratigraphic height (m)",
    shape = expression(delta^18 * O[carb] ~ "source"),
    color = "Temperature source"
  ) +
  theme_classic()


# ---- Save outputs ----
write_csv(
  BHB_water_inputs,
  here("data", "processed", "BHB_soilwater_inputs.csv")
)

write_csv(
  BHB_d18Ow_recon,
  here("data", "processed", "BHB_d18Ow_reconstruction.csv")
)

write_csv(
  missing_d18Ocarb,
  here("data", "processed", "BHB_missing_selected_d18Ocarb.csv")
)