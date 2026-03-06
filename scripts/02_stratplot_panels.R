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