library(tidyverse)
library(here)

# ---- Settings ----
BHB_paleolat <- 48
Kelson_paleolat <- 33   # approximate Big Bend / Tornillo paleolat; adjust if needed

# ---- Load LTG files ----
ltg <- bind_rows(
  read_csv(here("data", "raw", "PhanDA_LTG_Selandian&Danian.csv")) %>%
    mutate(stage_group = "Selandian-Danian"),
  read_csv(here("data", "raw", "PhanDA_LTG_Thanetian.csv")) %>%
    mutate(stage_group = "Thanetian"),
  read_csv(here("data", "raw", "PhanDA_LTG_Ypresian.csv")) %>%
    mutate(stage_group = "Ypresian")
)

# ---- Load Kelson processed data ----
Kelson_Tornillo_D47 <- read_csv(
  here("data", "processed", "Kelson_Tornillo_D47_strat_age_filtered.csv")
)

# ---- Prep Kelson primary/micrite temps ----
Kelson_temp_ltg <- Kelson_Tornillo_D47 %>%
  filter(
    petro_class == "Micrite",
    !is.na(Age_Ma),
    !is.na(T47_C)
  ) %>%
  mutate(
    paleolatitude = Kelson_paleolat,
    stage_group = case_when(
      Age_Ma > 56.00 & Age_Ma <= 59.24 ~ "Thanetian",
      Age_Ma >= 47.83 & Age_Ma <= 56.00 ~ "Ypresian",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(stage_group %in% c("Thanetian", "Ypresian"))

# ---- Prep BHB temperature data ----

BHB_multiproxy_final <- read.csv(here::here("data", "processed",
  "BHB_multiproxy_final.csv"))


BHB_temp_ltg <- BHB_multiproxy_final %>%
  filter(!is.na(Age_Ma)) %>%
  mutate(
    paleolatitude = BHB_paleolat,
    stage_group = case_when(
      Age_Ma > 56.00 & Age_Ma <= 59.24 ~ "Thanetian",
      Age_Ma >= 47.83 & Age_Ma <= 56.00 ~ "Ypresian",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(stage_group))

# ---- Plot Thanetian and Ypresian LTGs with BHB + Kelson temps ----
for (this_stage in c("Thanetian", "Ypresian")) {
  
  ltg_i <- ltg %>%
    filter(stage_group == this_stage)
  
  bhb_i <- BHB_temp_ltg %>%
    filter(
      stage_group == this_stage,
      !is.na(IPLD47_mean_T47_C)
    )
  
  kelson_i <- Kelson_temp_ltg %>%
    filter(stage_group == this_stage)
  
  y_lims <- range(
    c(
      ltg_i$LTG_05,
      ltg_i$LTG_95,
      bhb_i$IPLD47_mean_T47_C,
      bhb_i$IPLD47_mean_T47_C - bhb_i$IPLD47_se_T47_C,
      bhb_i$IPLD47_mean_T47_C + bhb_i$IPLD47_se_T47_C,
      kelson_i$T47_C,
      kelson_i$T47_C - kelson_i$T47_se_C,
      kelson_i$T47_C + kelson_i$T47_se_C
    ),
    na.rm = TRUE
  )
  
  p <- ggplot() +
    geom_ribbon(
      data = ltg_i,
      aes(x = Latitude, ymin = LTG_05, ymax = LTG_95),
      alpha = 0.15
    ) +
    geom_ribbon(
      data = ltg_i,
      aes(x = Latitude, ymin = LTG_16, ymax = LTG_84),
      alpha = 0.30
    ) +
    geom_line(
      data = ltg_i,
      aes(x = Latitude, y = LTG_50),
      linewidth = 1
    ) +
    
    # BHB temps
    geom_errorbar(
      data = bhb_i,
      aes(
        x = paleolatitude,
        ymin = IPLD47_mean_T47_C - IPLD47_se_T47_C,
        ymax = IPLD47_mean_T47_C + IPLD47_se_T47_C
      ),
      width = 1,
      alpha = 0.75,
      na.rm = TRUE
    ) +
    geom_point(
      data = bhb_i,
      aes(
        x = paleolatitude,
        y = IPLD47_mean_T47_C,
        shape = "BHB IPL micrite"
      ),
      size = 2,
      alpha = 0.85,
      na.rm = TRUE
    ) +
    
    # Kelson primary/micrite temps
    geom_errorbar(
      data = kelson_i,
      aes(
        x = paleolatitude,
        ymin = T47_C - T47_se_C,
        ymax = T47_C + T47_se_C
      ),
      width = 1,
      alpha = 0.75,
      na.rm = TRUE
    ) +
    geom_point(
      data = kelson_i,
      aes(
        x = paleolatitude,
        y = T47_C,
        shape = "Kelson Tornillo micrite"
      ),
      size = 2,
      alpha = 0.85,
      na.rm = TRUE
    ) +
    
    coord_cartesian(ylim = y_lims) +
    scale_shape_manual(
      values = c(
        "BHB IPL micrite" = 16,
        "Kelson Tornillo micrite" = 17
      )
    ) +
    labs(
      title = paste("Latitudinal temperature gradient:", this_stage),
      subtitle = paste0(
        "BHB n = ", nrow(bhb_i),
        "; Kelson micrite n = ", nrow(kelson_i)
      ),
      x = "Paleolatitude (°N)",
      y = "Temperature (°C)",
      shape = NULL,
      caption = "Line = LTG median; dark ribbon = 16–84%; light ribbon = 5–95%"
    ) +
    theme_classic()
  
  print(p)
}

# --- Plot soil water D17O vs d18O ---------------

# ---- Reconstructed soil-water isotope space by PETM interval ----

library(tidyverse)
library(here)

# PETM stratigraphic boundaries, inclusive
petm_base_m <- 1506
petm_top_m  <- 1543

# ---- Kelson et al. observed modern soil waters ----

kelson_soilwater <- read_csv(
  here(
    "data",
    "excel files",
    "jrkelson-CZ17O_soilwater-efc3bd3",
    "sw.csv"
  ),
  show_col_types = FALSE
) %>%
  # sw.csv contains a duplicated header as its first data row
  filter(Identifier_1 != "Identifier_1") %>%
  mutate(
    across(
      c(dp18O, dp18O_se, D17O_pmg, D17O_err),
      as.numeric
    )
  ) %>%
  filter(
    !is.na(dp18O),
    !is.na(D17O_pmg),
    
    # Retain the four modern soil-water field sites
    siteID.1 %in% c("MOJ", "JOR", "REY", "ESGR"),
    
    # Exclude the unrelated JOR Red Lake playa samples
    siteID.1 != "JOR" | siteID.2 == "CSAND"
  )

# ---- Reconstructed Bighorn Basin soil waters ----


BHB_water_isotope_space <- BHB_multiproxy_final %>%
  filter(
    !is.na(strat_height_m),
    !is.na(d18Ow_mean_vsmow),
    !is.na(D17Orsw_mean_permeg)
  ) %>%
  mutate(
    # Convert conventional δ18Owater to logarithmic δ′18Owater
    dp18Ow_mean = 1000 * log1p(d18Ow_mean_vsmow / 1000),
    
    # Convert the δ18Owater 95% bounds to δ′18Owater
    dp18Ow_lower95 = 1000 * log1p(d18Ow_lower95_vsmow / 1000),
    dp18Ow_upper95 = 1000 * log1p(d18Ow_upper95_vsmow / 1000),
    
    strat_interval = case_when(
      strat_height_m < petm_base_m ~ "Before PETM",
      between(strat_height_m, petm_base_m, petm_top_m) ~ "During PETM",
      strat_height_m > petm_top_m ~ "After PETM"
    ),
    
    strat_interval = factor(
      strat_interval,
      levels = c("Before PETM", "During PETM", "After PETM")
    )
  )

# ---- Interval mean Δ′17Owater and confidence bands ----

D17O_interval_summary <- BHB_water_isotope_space %>%
  group_by(strat_interval) %>%
  summarise(
    n = sum(!is.na(D17Orsw_mean_permeg)),
    mean_D17O = mean(D17Orsw_mean_permeg, na.rm = TRUE),
    sd_D17O = sd(D17Orsw_mean_permeg, na.rm = TRUE),
    se_D17O = sd_D17O / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    # Confidence intervals around the interval mean
    lower50 = mean_D17O - qt(0.75,  df = n - 1) * se_D17O,
    upper50 = mean_D17O + qt(0.75,  df = n - 1) * se_D17O,
    
    lower80 = mean_D17O - qt(0.90,  df = n - 1) * se_D17O,
    upper80 = mean_D17O + qt(0.90,  df = n - 1) * se_D17O,
    
    lower95 = mean_D17O - qt(0.975, df = n - 1) * se_D17O,
    upper95 = mean_D17O + qt(0.975, df = n - 1) * se_D17O
  )

D17O_interval_summary

# ---- Multipanel figure ----

# ---- Multipanel figure ----

reconstruction_color <- "#2166AC"
  reconstruction_dark  <- "#053061"
    reconstruction_fill  <- "#4393C3"
      
    p_D17O_d18Owater_PETM <- ggplot() +
      
      # Outermost interval: 95% CI
      geom_rect(
        data = D17O_interval_summary,
        aes(
          xmin = -Inf,
          xmax = Inf,
          ymin = lower95,
          ymax = upper95
        ),
        inherit.aes = FALSE,
        fill = reconstruction_color,
        alpha = 0.07,
        color = NA
      ) +
      
      # Middle interval: 80% CI
      geom_rect(
        data = D17O_interval_summary,
        aes(
          xmin = -Inf,
          xmax = Inf,
          ymin = lower80,
          ymax = upper80
        ),
        inherit.aes = FALSE,
        fill = reconstruction_color,
        alpha = 0.11,
        color = NA
      ) +
      
      # Innermost interval: 50% CI
      geom_rect(
        data = D17O_interval_summary,
        aes(
          xmin = -Inf,
          xmax = Inf,
          ymin = lower50,
          ymax = upper50
        ),
        inherit.aes = FALSE,
        fill = reconstruction_color,
        alpha = 0.16,
        color = NA
      ) +
      
      # Mean Δ′17Owater for each stratigraphic interval
      geom_hline(
        data = D17O_interval_summary,
        aes(yintercept = mean_D17O),
        inherit.aes = FALSE,
        color = reconstruction_dark,
        linewidth = 0.9
      ) +
      
      # Kelson modern soil-water x uncertainty
      geom_errorbarh(
        data = kelson_soilwater,
        aes(
          y = D17O_pmg,
          xmin = dp18O - 1.96 * dp18O_se,
          xmax = dp18O + 1.96 * dp18O_se
        ),
        inherit.aes = FALSE,
        height = 0,
        linewidth = 0.25,
        color = "grey72",
        alpha = 0.45
      ) +
      
      # Kelson modern soil-water y uncertainty
      geom_errorbar(
        data = kelson_soilwater,
        aes(
          x = dp18O,
          ymin = D17O_pmg - D17O_err,
          ymax = D17O_pmg + D17O_err
        ),
        inherit.aes = FALSE,
        width = 0,
        linewidth = 0.25,
        color = "grey72",
        alpha = 0.45
      ) +
      
      # Kelson modern soil-water observations
      geom_point(
        data = kelson_soilwater,
        aes(x = dp18O, y = D17O_pmg),
        inherit.aes = FALSE,
        shape = 16,
        size = 1.5,
        color = "grey60",
        alpha = 0.45
      ) +
      
      # Reconstructed δ′18Owater 95% uncertainty
      geom_errorbarh(
        data = BHB_water_isotope_space,
        aes(
          y = D17Orsw_mean_permeg,
          xmin = dp18Ow_lower95,
          xmax = dp18Ow_upper95
        ),
        inherit.aes = FALSE,
        height = 0,
        linewidth = 0.45,
        color = reconstruction_color,
        alpha = 0.75
      ) +
      
      # Reconstructed Δ′17Owater uncertainty
      geom_errorbar(
        data = BHB_water_isotope_space,
        aes(
          x = dp18Ow_mean,
          ymin = D17Orsw_mean_permeg - IPL17O_ci95_Dp17Ocarb_adj,
          ymax = D17Orsw_mean_permeg + IPL17O_ci95_Dp17Ocarb_adj
        ),
        inherit.aes = FALSE,
        width = 0,
        linewidth = 0.45,
        color = reconstruction_color,
        alpha = 0.75
      ) +
      
      # Reconstructed soil-water estimates
      geom_point(
        data = BHB_water_isotope_space,
        aes(x = dp18Ow_mean, y = D17Orsw_mean_permeg),
        inherit.aes = FALSE,
        shape = 21,
        size = 3,
        stroke = 0.7,
        color = reconstruction_dark,
        fill = reconstruction_fill
      ) +
      
      facet_wrap(
        ~ strat_interval,
        nrow = 1,
        drop = FALSE
      ) +
      
      labs(
        x = expression(delta*"'"^18*"O"[water]*" (‰ VSMOW)"),
        y = expression(Delta*"'"^17*"O"[water]*" (per meg)"),
        caption = paste0(
          "PETM interval: ", petm_base_m, "\u2013", petm_top_m,
          " m. Grey points are observed modern soil waters from Kelson et al. ",
          "Blue horizontal lines show interval means; shading shows ",
          "50%, 80%, and 95% confidence intervals."
        )
      ) +
      
      theme_classic(base_size = 11) +
      theme(
        strip.background = element_rect(
          fill = "grey92",
          color = "grey35",
          linewidth = 0.4
        ),
        strip.text = element_text(face = "bold"),
        panel.spacing = grid::unit(1, "lines"),
        axis.title = element_text(face = "bold"),
        plot.caption = element_text(
          size = 8.5,
          color = "grey35",
          hjust = 0
        )
      )
    
    p_D17O_d18Owater_PETM
    
    ggsave(
      here(
        "figures",
        "BHB_D17Owater_d18Owater_PETM_intervals.png"
      ),
      plot = p_D17O_d18Owater_PETM,
      width = 12,
      height = 4.5,
      dpi = 600
    )

    p_D17O_d18Owater_PETM
    


# ---- Save figure ----

    ggsave(
      here(
        "figures",
        "BHB_D17Owater_d18Owater_PETM_intervals.png"
      ),
      plot = p_D17O_d18Owater_PETM,
      width = 12,
      height = 4.5,
      dpi = 600
    )
    
    
    
    
# ---- Temperature versus reconstructed Δ′17O soil water ----
    # ---- Temperature versus reconstructed Δ′17O soil water ----
    # Vertical error bars use the carbonate Δ′17O 95% CI.
    
    BHB_T_D17Orsw <- BHB_multiproxy_final %>%
      filter(
        !is.na(T_recon_C),
        !is.na(T_recon_se_C),
        !is.na(D17Orsw_mean_permeg),
        !is.na(IPL17O_ci95_Dp17Ocarb_adj)
      )
    
    p_T_D17Orsw <- ggplot(
      BHB_T_D17Orsw,
      aes(
        x = T_recon_C,
        y = D17Orsw_mean_permeg
      )
    ) +
      
      # Temperature 95% uncertainty
      geom_errorbarh(
        aes(
          xmin = T_recon_C - 1.96 * T_recon_se_C,
          xmax = T_recon_C + 1.96 * T_recon_se_C
        ),
        height = 0,
        linewidth = 0.45,
        color = "#2166AC",
        alpha = 0.65
      ) +
      
      # Carbonate Δ′17O 95% CI centered on reconstructed Δ′17Owater
      geom_errorbar(
        aes(
          ymin = D17Orsw_mean_permeg -
            IPL17O_ci95_Dp17Ocarb_adj,
          ymax = D17Orsw_mean_permeg +
            IPL17O_ci95_Dp17Ocarb_adj
        ),
        width = 0,
        linewidth = 0.45,
        color = "#2166AC",
        alpha = 0.65
      ) +
      
      # Linear trend with 95% confidence ribbon
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        level = 0.95,
        color = "#053061",
        fill = "#4393C3",
        linewidth = 0.9,
        alpha = 0.20
      ) +
      
      geom_point(
        shape = 21,
        size = 3.2,
        stroke = 0.8,
        color = "#053061",
        fill = "#4393C3"
      ) +
      
      geom_point(
        shape = 21,
        size = 3.2,
        stroke = 0.8,
        color = "#053061",
        fill = "#4393C3"
      ) +
      
      labs(
        x = expression(
          "Reconstruction temperature (" * degree * C * ")"
        ),
        y = expression(
          Delta * "'"^17 * "O"[soil~water] * " (per meg)"
        )
      ) +
      
      theme_classic(base_size = 11) +
      
      theme(
        axis.title = element_text(face = "bold")
      )
    
    p_T_D17Orsw
    
    ggsave(
      here(
        "figures",
        "BHB_temperature_vs_D17Orsw.png"
      ),
      plot = p_T_D17Orsw,
      width = 7,
      height = 5.5,
      dpi = 600
    )
    
    # ---- Temperature versus carbonate Δ′17O ----
    
    BHB_T_D17Ocarb <- BHB_multiproxy_final %>%
      filter(
        !is.na(T_recon_C),
        !is.na(T_recon_se_C),
        !is.na(IPL17O_mean_Dp17Ocarb),
        !is.na(IPL17O_ci95_Dp17Ocarb_adj)
      )
    
    p_T_D17Ocarb <- ggplot(
      BHB_T_D17Ocarb,
      aes(
        x = T_recon_C,
        y = IPL17O_mean_Dp17Ocarb
      )
    ) +
      
      # Temperature 95% uncertainty
      geom_errorbarh(
        aes(
          xmin = T_recon_C - 1.96 * T_recon_se_C,
          xmax = T_recon_C + 1.96 * T_recon_se_C
        ),
        height = 0,
        linewidth = 0.45,
        color = "#2166AC",
        alpha = 0.65
      ) +
      
      # Carbonate Δ′17O 95% uncertainty
      geom_errorbar(
        aes(
          ymin = IPL17O_mean_Dp17Ocarb -
            IPL17O_ci95_Dp17Ocarb_adj,
          ymax = IPL17O_mean_Dp17Ocarb +
            IPL17O_ci95_Dp17Ocarb_adj
        ),
        width = 0,
        linewidth = 0.45,
        color = "#2166AC",
        alpha = 0.65
      ) +
      
      # Linear trend and 95% confidence ribbon
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        level = 0.95,
        color = "#053061",
        fill = "#4393C3",
        linewidth = 0.9,
        alpha = 0.20
      ) +
      
      geom_point(
        shape = 21,
        size = 3.2,
        stroke = 0.8,
        color = "#053061",
        fill = "#4393C3"
      ) +
      
      labs(
        x = expression(
          "Reconstruction temperature (" * degree * C * ")"
        ),
        y = expression(
          Delta * "'"^17 * "O"[carbonate] * " (per meg)"
        )
      ) +
      
      theme_classic(base_size = 11) +
      
      theme(
        axis.title = element_text(face = "bold")
      )
    
    p_T_D17Ocarb
    
    # Regression statistics
    T_D17Ocarb_model <- lm(
      IPL17O_mean_Dp17Ocarb ~ T_recon_C,
      data = BHB_T_D17Ocarb
    )
    
    summary(T_D17Ocarb_model)
    
    # Save figure
    ggsave(
      here(
        "figures",
        "BHB_temperature_vs_D17Ocarb.png"
      ),
      plot = p_T_D17Ocarb,
      width = 7,
      height = 5.5,
      dpi = 600
    )
    
    