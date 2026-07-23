# Legacy script: superseded by 05_reconstruct_CFB_soilwater.R
# Purpose: Probabilistically reconstruct soil-water D17O from carbonate
#          triple oxygen isotope data and reconstructed temperature.

# ---- Load packages ----
library(tidyverse)
library(here)

# ---- Settings ----
n_mc <- 10000
set.seed(123)

theta_carb <- 0.5250
lambda_ref <- 0.528

# ---- Load d18O soil-water reconstruction output ----
BHB_d18Ow_recon <- read_csv(
  here("data", "processed", "BHB_d18Ow_reconstruction.csv")
)

# ---- Keep horizons with IPL D17O carbonate data ----
BHB_D17O_inputs <- BHB_d18Ow_recon %>%
  filter(
    !is.na(IPL17O_mean_Dp17Ocarb),
    !is.na(IPL17O_mean_d17Ocarb),
    !is.na(IPL17O_mean_d18Ocarb),
    !is.na(T_recon_C),
    !is.na(T_recon_se_C)
  ) %>%
  mutate(
    D17Ocarb_se_used_permeg = pmax(
      case_when(
        !is.na(IPL17O_se_Dp17Ocarb_adj) & IPL17O_se_Dp17Ocarb_adj > 0 ~ IPL17O_se_Dp17Ocarb_adj,
        !is.na(IPL17O_se_Dp17Ocarb) & IPL17O_se_Dp17Ocarb > 0 ~ IPL17O_se_Dp17Ocarb,
        TRUE ~ 12
      ),
      12,
      na.rm = TRUE
    ),
    
    d17Ocarb_se_used = case_when(
      !is.na(IPL17O_se_dp17Ocarb) & IPL17O_se_dp17Ocarb > 0 ~ IPL17O_se_dp17Ocarb,
      TRUE ~ 0.15
    ),
    
    d18Ocarb_se_used = case_when(
      !is.na(IPL17O_se_dp18Ocarb) & IPL17O_se_dp18Ocarb > 0 ~ IPL17O_se_dp18Ocarb,
      TRUE ~ 0.15
    )
  )

# ---- Calcite-water fractionation ----
calcite_water_alpha18 <- function(T_C) {
  T_K <- T_C + 273.15
  exp((18.03 * (1000 / T_K) - 32.42) / 1000)
}

# ---- Reconstruct D17Owater from carbonate d17O, d18O, and T ----
compute_D17Owater <- function(d17O_carb, d18O_carb, T_C,
                              theta = theta_carb,
                              lambda = lambda_ref) {
  
  alpha18 <- calcite_water_alpha18(T_C)
  alpha17 <- exp(theta * log(alpha18))
  
  d17O_water <- ((1000 + d17O_carb) / alpha17) - 1000
  d18O_water <- ((1000 + d18O_carb) / alpha18) - 1000
  
  dp17O_water <- 1000 * log1p(d17O_water / 1000)
  dp18O_water <- 1000 * log1p(d18O_water / 1000)
  
  D17O_water_permeg <- (dp17O_water - lambda * dp18O_water) * 1000
  
  return(D17O_water_permeg)
}

# ---- Monte Carlo reconstruction for one horizon ----
mc_D17Owater_one <- function(d17O_mean, d17O_se,
                             d18O_mean, d18O_se,
                             T_mean, T_se,
                             n = 10000) {
  
  d17O_sim <- rnorm(n, mean = d17O_mean, sd = d17O_se)
  d18O_sim <- rnorm(n, mean = d18O_mean, sd = d18O_se)
  T_sim    <- rnorm(n, mean = T_mean, sd = T_se)
  
  D17Owater_sim <- compute_D17Owater(
    d17O_carb = d17O_sim,
    d18O_carb = d18O_sim,
    T_C = T_sim
  )
  
  tibble(
    D17Orsw_mean_permeg = mean(D17Owater_sim, na.rm = TRUE),
    D17Orsw_median_permeg = median(D17Owater_sim, na.rm = TRUE),
    D17Orsw_sd_permeg = sd(D17Owater_sim, na.rm = TRUE),
    D17Orsw_lower95_permeg = quantile(D17Owater_sim, 0.025, na.rm = TRUE),
    D17Orsw_upper95_permeg = quantile(D17Owater_sim, 0.975, na.rm = TRUE)
  )
}

# ---- Run Monte Carlo reconstruction ----
BHB_D17Orsw_recon <- BHB_D17O_inputs %>%
  mutate(
    mc = pmap(
      list(
        IPL17O_mean_d17Ocarb,
        d17Ocarb_se_used,
        IPL17O_mean_d18Ocarb,
        d18Ocarb_se_used,
        T_recon_C,
        T_recon_se_C
      ),
      ~ mc_D17Owater_one(
        d17O_mean = ..1,
        d17O_se = ..2,
        d18O_mean = ..3,
        d18O_se = ..4,
        T_mean = ..5,
        T_se = ..6,
        n = n_mc
      )
    )
  ) %>%
  unnest(mc)

# ---- Quick QC ----
BHB_D17Orsw_recon %>%
  summarise(
    n_horizons = n(),
    min_D17Orsw = min(D17Orsw_mean_permeg, na.rm = TRUE),
    max_D17Orsw = max(D17Orsw_mean_permeg, na.rm = TRUE),
    mean_D17Orsw = mean(D17Orsw_mean_permeg, na.rm = TRUE)
  )

# ---- Stratigraphic plot ----

ggplot(
  BHB_D17Orsw_recon,
  aes(
    x = D17Orsw_mean_permeg,
    y = strat_height_m
  )
) +
  geom_point(size = 3) +
  theme_classic()

ggplot(
  BHB_D17Orsw_recon,
  aes(x = D17Orsw_mean_permeg, y = strat_height_m)
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
      xmin = D17Orsw_mean_permeg - 1.96 * IPL17O_se_Dp17Ocarb_adj,
      xmax = D17Orsw_mean_permeg + 1.96 * IPL17O_se_Dp17Ocarb_adj
    ),
    height = 0,
    alpha = 0.5
  ) +
  geom_point(
    aes(color = T_recon_source),
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
    breaks = seq(500, 1800, by = 100)
  ) +
  labs(
    x = expression(Delta*"'"^17 * O[rsw] ~ "(per meg)"),
    y = "Stratigraphic height (m)",
    color = "Temperature source"
  ) +
  theme_classic()

ggplot(
  BHB_D17Orsw_recon,
  aes(x = D17Orsw_mean_permeg, y = strat_height_m)
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
      xmin = D17Orsw_lower95_permeg,
      xmax = D17Orsw_upper95_permeg
    ),
    height = 0,
    alpha = 0.5
  ) +
  geom_point(
    aes(color = T_recon_source),
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
    breaks = seq(500, 1700, by = 100)
  ) +
  labs(
    x = expression(Delta*"'"^17 * O[rsw] ~ "(per meg)"),
    y = "Stratigraphic height (m)",
    color = "Temperature source"
  ) +
  theme_classic()

# ---- Join D17O soil-water reconstructions back to full multiproxy table ----

BHB_multiproxy_final <- BHB_d18Ow_recon %>%
  left_join(
    BHB_D17Orsw_recon %>%
      select(
        MLA_horizon_id,
        D17Orsw_mean_permeg,
        D17Orsw_median_permeg,
        D17Orsw_sd_permeg,
        D17Orsw_lower95_permeg,
        D17Orsw_upper95_permeg
      ),
    by = "MLA_horizon_id"
  )


# ---- Save outputs ----
write_csv(
  BHB_multiproxy_final,
  here("data", "processed", "BHB_multiproxy_final.csv")
)

write_csv(
  BHB_D17O_inputs,
  here("data", "processed", "BHB_D17O_soilwater_inputs.csv")
)

write_csv(
  BHB_D17Orsw_recon,
  here("data", "processed", "BHB_D17Orsw_reconstruction.csv")
)
