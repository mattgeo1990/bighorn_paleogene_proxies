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
