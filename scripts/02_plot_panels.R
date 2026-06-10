# 02_stratplot_panels.R
# Create stratplot panels of proxy records

library(tidyverse)
library(here)
library(patchwork)
library(zoo)
library(grid)  # Needed for unit() in geom_segment()

# Load processed data
df <- read_csv(here("data", "processed", "PETM_combined.csv"))

#### AGE MODEL ------------

# Load modified Westerhold et al. (2017) age model
prelim_age_model <- read_csv(here("data", "raw", "approx_age_mdl_PETM_PCB.csv"))

# Filter to valid rows and fit linear model
model_data <- prelim_age_model %>%
  filter(!is.na(est_Depth_m_PCB_outcrop), !is.na(Age_Ma))

age_model <- lm(Age_Ma ~ est_Depth_m_PCB_outcrop, data = model_data)

# Create prediction over observed range (plot only)
depth_seq <- seq(min(model_data$est_Depth_m_PCB_outcrop),
                 max(model_data$est_Depth_m_PCB_outcrop),
                 length.out = 200)

pred_df <- tibble(
  est_Depth_m_PCB_outcrop = depth_seq,
  Age_Ma = predict(age_model, newdata = tibble(est_Depth_m_PCB_outcrop = depth_seq))
)

# Optional: plot age model (optional preview)
ggplot(model_data, aes(x = est_Depth_m_PCB_outcrop, y = Age_Ma)) +
  geom_point(color = "darkblue", size = 2) +
  geom_line(data = pred_df, aes(x = est_Depth_m_PCB_outcrop, y = Age_Ma), color = "red", linewidth = 1) +
  scale_x_reverse() +
  theme_minimal()

# Extrapolate full range (1300–1800 m)
depth_seq_full <- seq(1300, 1800, length.out = 501)

age_model_full <- tibble(
  est_Depth_m_PCB_outcrop = depth_seq_full,
  Age_Ma = predict(age_model, newdata = tibble(est_Depth_m_PCB_outcrop = depth_seq_full))
)

# Interpolate onto df$Strat_m_Bowen
df <- df %>%
  mutate(Age_Ma = approx(
    x = age_model_full$est_Depth_m_PCB_outcrop,
    y = age_model_full$Age_Ma,
    xout = strat_height_m,
    rule = 2
  )$y)

#### END AGE MODEL ----------

### Clumped age plot from other script --------

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



### Resuming current script -------


# Placeholder mammal turnover events
turnover_events <- tibble(
  label = c("Tiff-Wa", "Wa0", "Wa1", "Wa2"),
  Age_Ma = c(56.33, 56.0, 55.8, 55.6),
  x_pos = c(3.2, 3.2, 3.2, 3.2)
)

# Theme for panels
theme_petm <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.title.x = element_text(margin = margin(t = 5)),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

pt_size <- 2.5




# Δ17O panel
# Sort and compute rolling mean
df_sorted <- df %>%
  arrange(Age_Ma) %>%
  mutate(
    D17O_roll = zoo::rollmean(mean_D17Ocarb, k = 3, fill = NA, align = "center")
  )

# Plot Δ17O vs Age (rolling mean)

# Clean, sort, and filter NAs
df_sorted <- df %>%
  filter(!is.na(mean_D17Ocarb), !is.na(Age_Ma)) %>%
  arrange(Age_Ma) %>%
  mutate(
    D17O_roll = rollmean(mean_D17Ocarb, k = 5, fill = NA, align = "center")  # use k = 5 for more smoothing
  )

# Plot
df_sorted <- df %>%
  arrange(Age_Ma) %>%
  mutate(D17O_roll = zoo::rollmean(mean_D17Ocarb, k = 3, fill = NA, align = "center"))

# Use new data frame for plotting, with x mapped to D17O_roll
# Sort by Age_Ma and compute rolling mean for Δ17Ocarb
df_sorted <- df %>%
  mutate(D17O_roll = zoo::rollmean(mean_D17Ocarb, k = 3, fill = NA, align = "center")) %>%
  arrange(Age_Ma)  # this is what ensures clean left-to-right (x) lines

# Δ17O panel
p1 <- ggplot(df_sorted, aes(y = Age_Ma)) +
  geom_point(aes(x = mean_D17Ocarb), size = pt_size, color = "darkblue") +
  geom_line(aes(x = D17O_roll), linewidth = 0.5, color = "darkblue") +
  scale_y_reverse() +
  labs(
    x = expression(Delta*minute^17*O[carb]~"(per meg)"),
    y = "Age (Ma)", title = expression(Delta^17*O[carb])
  ) +
  geom_segment(data = turnover_events,
               aes(x = x_pos, xend = x_pos - 0.2,
                   y = Age_Ma, yend = Age_Ma),
               arrow = arrow(length = unit(0.1, "cm")),
               color = "black", linewidth = 0.3) +
  geom_text(data = turnover_events,
            aes(x = x_pos + 0.05, y = Age_Ma, label = label),
            hjust = 0, size = 3.3) +
  theme_petm
p1


# T47 panel
p2 <- ggplot(df, aes(x = T47_C, y = Age_Ma)) +
  geom_point(size = pt_size, color = "firebrick") +
  geom_line(linewidth = 0.5, color = "firebrick") +
  scale_y_reverse() +
  labs(
    x = expression(paste(Delta[47]," T (",degree,"C)")),
    y = NULL, title = expression(Delta[47]~"Temperature")
  ) +
  geom_segment(data = turnover_events,
               aes(x = 37.5, xend = 36.5,
                   y = Age_Ma, yend = Age_Ma),
               arrow = arrow(length = unit(0.1, "cm")),
               color = "black", linewidth = 0.3) +
  geom_text(data = turnover_events,
            aes(x = 37.6, y = Age_Ma, label = label),
            hjust = 0, size = 3.3) +
  theme_petm
p2

# δ13C panel
p3 <- ggplot(df, aes(x = d13C_carb, y = Age_Ma)) +
  geom_point(size = pt_size, color = "forestgreen") +
  geom_line(linewidth = 0.5, color = "forestgreen") +
  scale_y_reverse() +
  labs(
    x = expression(delta^13*C[carb]~"\u2030 VPDB"),
    y = NULL, title = expression(delta^13*C[carb])
  ) +
  geom_segment(data = turnover_events,
               aes(x = -8, xend = -8.4,
                   y = Age_Ma, yend = Age_Ma),
               arrow = arrow(length = unit(0.1, "cm")),
               color = "black", linewidth = 0.3) +
  geom_text(data = turnover_events,
            aes(x = -7.85, y = Age_Ma, label = label),
            hjust = 0, size = 3.3) +
  theme_petm
p3

# Combine panels
final_fig <- p1 + p2 + p3 + plot_layout(ncol = 3)

# Ensure output directory exists
if (!dir.exists(here("figures"))) dir.create(here("figures"))

# Export figure
ggsave(here("figures", "PETM_panel_v1.png"),
       final_fig, width = 11, height = 6, dpi = 600)
final_fig