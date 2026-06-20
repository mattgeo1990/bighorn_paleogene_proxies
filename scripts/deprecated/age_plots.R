### Load Libraries ----------
library(dplyr)
library(readr)
library(here)
library(isogeochem)
library(janitor)
library(ggplot2)
library(patchwork)
library(stringr)
library(plotly)
library(tidyr)
library(ggrepel)



BP_May2026_summary  <- read.csv(here("data", "raw", "BHB Paleogene Summary May 2026.csv"))

str(BP_May2026_summary)  

BP_May2026_summary <- subset(
  BP_May2026_summary,
  !grepl("SPAR", Sample.ID)
)

# Summarize by strat
BP_May2026_strat_summary <- BP_May2026_summary %>%
  group_by(Strat) %>%
  summarise(
    n = sum(!is.na(T.D47..Petersen)),
    mean_T47 = mean(T.D47..Petersen, na.rm = TRUE),
    sd_T47 = sd(T.D47..Petersen, na.rm = TRUE),
    se_T47 = sd_T47 / sqrt(n),
    .groups = "drop"
  ) %>%
  arrange(Strat)



# Strat Plot
ggplot(BP_May2026_strat_summary,
       aes(x = mean_T47, y = Strat)) +
  
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 1500, ymax = 1545,
           fill = "red", alpha = 0.15) +
  
  annotate("text",
           x = Inf, y = 1522.5,
           label = "PETM",
           hjust = 1.1,
           color = "red",
           fontface = "bold") +
  
  geom_errorbarh(
    aes(xmin = mean_T47 - sd_T47,
        xmax = mean_T47 + sd_T47),
    height = 0
  ) +
  
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = seq(15, 60, by = 5)
  ) +
  scale_y_continuous(
    breaks = seq(500, 1800, by = 100)
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# apply age model

CFB_age_model <- read.csv(here("data", "raw", "approx_age_mdl_PETM_PCB.csv"))

str(CFB_age_model)

BP_May2026_strat_summary <- BP_May2026_strat_summary %>%
  mutate(
    Age_Ma = approx(
      x = CFB_age_model$est_Depth_m_PCB_outcrop,
      y = CFB_age_model$Age_Ma,
      xout = Strat,
      rule = 1
    )$y
  )

ggplot(BP_May2026_strat_summary,
       aes(x = mean_T47, y = Age_Ma)) +
  
  geom_errorbarh(aes(xmin = mean_T47 - sd_T47,
                     xmax = mean_T47 + sd_T47),
                 height = 0) +
  
  geom_point(size = 3) +
  
  scale_y_reverse() +
  
  labs(
    x = expression(T[Delta47]~"Petersen ("*degree*C*")"),
    y = "Age (Ma)",
    title = expression(T[Delta47]~"by Age")
  ) +
  
  theme_bw(base_size = 14)

