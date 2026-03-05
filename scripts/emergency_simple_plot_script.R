# ---- Packages ----
library(tidyverse)   # ggplot2, dplyr, tibble, readr, etc.
library(here)        # project-root file paths
install.packages("janitor")
library(janitor)     # clean_names()
install.packages("patchwork")
library(zoo)         # rollmean()
library(patchwork)   # panel assembly

### LOAD DATA -------

# Load CU Δ47 data and clean column names
CU_data <- read_csv(here("data", "raw", "PETM_clumped.csv")) %>%
  janitor::clean_names() %>%
  rename(Strat_m_Bowen = strat_height_m) %>%
  mutate(d18Oc_SMOW = to_VSMOW(d18o_vpbd, eq = "IUPAC")) %>%
  round_depth()

# Bowen dataset
Bowen_d13C <- read.csv(
  here("data", "raw", "Bowen2001_IsotopeData.csv")
)


# Load Koch δ13C data and clean
koch <- read_csv(here("data", "raw", "Koch_SC_nodules_isotopes.csv")) %>%
  clean_names() %>%
  rename(Strat_m_Bowen = strat_m) %>%
  round_depth()

# MArine reference d13C
benthic_d13C <- read.csv(here("data", "raw", "ODP1262_BarnettEtAl2019.csv"))
benthic_d13C <- benthic_d13C %>%
  mutate(Age_Ma = Age..ka.BP. / 1000)

# Sample metadata (depths)
samples <- read.csv(here("data", "raw", "Polecat_samples.csv"))


# Corrected Δ′17O (IPL data, replicate-level)
IPL_17O <- read.csv(here("data", "processed", "IPL17O_PETM_summary.csv"))

combined <- read.csv(here("data", "processed", "combined_petm_data.csv"))

plot(combined$d18Oc_SMOW, combined$D17Orfw_permeg)

### SUMMARIZE -------



# READY TO PLOT OR APPLY AGE MODEL

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
# ggplot(model_data, aes(x = est_Depth_m_PCB_outcrop, y = Age_Ma)) +
#   geom_point(color = "darkblue", size = 2) +
#   geom_line(data = pred_df, aes(x = est_Depth_m_PCB_outcrop, y = Age_Ma), color = "red", linewidth = 1) +
#   scale_x_reverse() +
#   theme_minimal()

# Extrapolate full range (1300–1800 m)
depth_seq_full <- seq(1300, 1800, length.out = 501)

age_model_full <- tibble(
  est_Depth_m_PCB_outcrop = depth_seq_full,
  Age_Ma = predict(age_model, newdata = tibble(est_Depth_m_PCB_outcrop = depth_seq_full))
)

# Interpolate onto IPL_17O$Strat_m_Bowen
IPL_17O <- IPL_17O %>%
  mutate(Age_Ma = approx(
    x = age_model_full$est_Depth_m_PCB_outcrop,
    y = age_model_full$Age_Ma,
    xout = Strat_m_Bowen,
    rule = 2
  )$y)

#### END AGE MODEL ----------

### D17O Age PCB  ------------
IPL_17O <- IPL_17O %>%
  arrange(Age_Ma) %>%
  mutate(rollmean = zoo::rollmean(mean_D17Ocarb, k = 3,
                                   fill = NA, align = "center"))


p <- ggplot(IPL_17O, aes(x = mean_D17Ocarb, y = Age_Ma)) +
  # Error bars (no tips)
  geom_errorbarh(
    aes(xmin = mean_D17Ocarb - se_D17Ocarb,
        xmax = mean_D17Ocarb + se_D17Ocarb),
    height = 0,  # removes the horizontal caps
    na.rm = TRUE,
    color = "gray40"
  ) +
  # Points: singletons (red)
  geom_point(
    aes(x = mean_D17Ocarb, y = Age_Ma),
    color = "black",
    size = 2.5,
    alpha = 1
  ) +
  # Line
  #geom_path(
  #  aes(x = mean_D17Ocarb, y = Age_Ma),
  #  color = "blue",
  #  linewidth = 0.5,
  # na.rm = TRUE
  #) +
  # Axes and theme
  scale_y_reverse(
    limits = c(56.5, 55.2),
    breaks = seq(56.5, 55.2, by = -0.1)
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

p + points(Koch~Strat_m~Koch$T47)

# Save to file
ggsave(filename = here::here("figures", "D17O_age_plot.png"),
       plot = p,
       width = 4, height = 4, dpi = 600)

### D17O Depth PCB  ------------

D17O_PCB_depth <- ggplot(IPL_17O, aes(x = mean_D17Ocarb, y = Strat_m_Bowen)) +
  # Error bars
  geom_errorbarh(
    aes(xmin = mean_D17Ocarb - se_D17Ocarb,
        xmax = mean_D17Ocarb + se_D17Ocarb),
    height = 0,
    na.rm = TRUE,
    color = "gray40"
  ) +
  # Points: n > 1 (black)
  geom_point(
    aes(x = mean_D17Ocarb, y = Strat_m_Bowen),
    color = "black",
    size = 1.5,
    alpha = 0.6
  ) +
  # Rolling mean ribbon
  #geom_ribbon(
  #  data = rolled,
  #  aes(x = roll_mean, y = Strat_m_Bowen,
  #      xmin = D17O_lower, xmax = D17O_upper),
  #  fill = "steelblue", alpha = 0.3,
  #  inherit.aes = FALSE
  #) +
  # Rolling mean line
  #geom_path(
  #  data = rolled %>% arrange(Strat_m_Bowen),
  #  aes(x = roll_mean, y = Strat_m_Bowen),
  #  color = "steelblue",
  # size = 1,
  #  inherit.aes = FALSE
  #) +
  scale_x_continuous(
    limits = c(-150, -70),
    breaks = seq(-150, -70, by = 10)
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

D17O_PCB_depth

# Save to file
ggsave(filename = here::here("figures", "D17O_PETM_only.png"),
       plot = D17O_PETM_only,
       width = 4, height = 4, dpi = 600)

CU_data

#### Benthic d13C (ODP 1262 South Atlantic) -------

marine <- ggplot(benthic_d13C, aes(x = d13C_VPDB, y = Age_Ma)) +
  geom_point(
    data = benthic_d13C,
    aes(x = d13C_VPDB, y = Age_Ma),
    color = "black",
    size = 1.5,
    alpha = 0.6
  ) +
  # Axes and theme
  scale_y_reverse(
    limits = c(56.5, 55.2),
    breaks = seq(56.5, 55.2, by = -0.1)
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

marine


# Save to file
ggsave(filename = here::here("figures", "benthic_d13C_ODP1262.png"),
       plot = marine,
       width = 4, height = 4, dpi = 600)


### D47 Depth CU Boulder ---------
# Fresh plot of t47 with horizontal error bars
D47_CU_depth <- ggplot(CU_data, aes(x = t47_c, y = Strat_m_Bowen)) +
  geom_point(color = "black", size = 2) +
  geom_errorbarh(aes(xmin = t47_c - t47_2se_c, xmax = t47_c + t47_2se_c),
                 height = 0, color = "gray40") +
  labs(
    x = expression(paste("Clumped-isotope temperature (", degree*C, ")")),
    y = "Stratigraphic height (m)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

# Print it
D47_CU_depth




### Benthic temp (ODP 1262 South Atlantic) --------------

bwt <- ggplot(benthic_d13C, aes(x = BWT_degC, y = Age_Ma)) +
  geom_point(
    data = benthic_d13C,
    aes(x = BWT_degC, y = Age_Ma),
    color = "black",
    size = 1.5,
    alpha = 0.6
  ) +
  # Axes and theme
  scale_y_reverse(
    limits = c(56.5, 55.2),
    breaks = seq(56.5, 55.2, by = -0.1)
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

bwt


# Save to file
ggsave(filename = here::here("figures", "benthic_temp_ODP1262.png"),
       plot = bwt,
       width = 4, height = 4, dpi = 600)


#### PCB Koch d13C ----------------------

# Option 1: rename for clarity
age_model_clean <- age_model_full %>%
  filter(!is.na(est_Depth_m_PCB_outcrop) & !is.na(Age_Ma)) %>%
  rename(Strat_m_Bowen = est_Depth_m_PCB_outcrop)

koch <- koch %>%
  mutate(Age_Ma = approx(
    x = age_model_clean$Strat_m_Bowen,
    y = age_model_clean$Age_Ma,
    xout = Strat_m_Bowen,
    rule = 2
  )$y)


kochcarb <- ggplot(koch, aes(x = d13c_vpdb, y = Age_Ma)) +
  geom_point(
    data = koch,
    aes(x = d13c_vpdb, y = Age_Ma),
    color = "black",
    size = 1.5,
    alpha = 0.6
  ) +
  # Axes and theme
  scale_y_reverse(
    limits = c(56.5, 55.2),
    breaks = seq(56.5, 55.2, by = -0.1)
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

kochcarb


# Save to file
ggsave(filename = here::here("figures", "kochcarb.png"),
       plot = kochcarb,
       width = 4, height = 4, dpi = 600)
# ---- Age plots: Δ′17O (IPL) and T47 (Koch) ----

age_limits <- c(56.5, 55.2)
age_breaks <- seq(56.5, 55.2, by = -0.1)

p_D17O_age <- ggplot(IPL_17O, aes(x = mean_D17Ocarb, y = Age_Ma)) +
  geom_errorbarh(
    aes(xmin = mean_D17Ocarb - se_D17Ocarb,
        xmax = mean_D17Ocarb + se_D17Ocarb),
    height = 0, na.rm = TRUE, color = "gray40"
  ) +
  geom_point(color = "black", size = 2.0, alpha = 0.9) +
  scale_y_reverse(limits = age_limits, breaks = age_breaks) +
  labs(
    x = expression(Delta*minute^17*O[carb]~("per meg")),
    y = "Age (Ma)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

p_T47_age <- ggplot(subset(koch, !is.na(t47) & !is.na(Age_Ma)),
                    aes(x = t47, y = Age_Ma)) +
  geom_point(color = "black", size = 2.0, alpha = 0.8) +
  scale_y_reverse(limits = age_limits, breaks = age_breaks) +
  labs(
    x = expression(T[47]~("(°C)")),
    y = NULL
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.y = element_blank()
  )

panel_D17O_T47_age <- p_D17O_age | p_T47_age
panel_D17O_T47_age

ggsave(
  filename = here::here("figures", "panel_D17O_T47_age.png"),
  plot = panel_D17O_T47_age,
  width = 7, height = 4, dpi = 600
)
# ---- Add Koch δ18O(VSMOW) and δ13C(VPDB) age panels ----
# Assumes koch has columns:
#   d18Ocarb_SMOW (or d18Ocarb_VSMOW) and d13Ccarb_VPDB (or d13C_VPDB) and Age_Ma

# If your Koch columns are actually d18ocarb_vsmow and d13c_vpdb (from clean_names), you can alias them:
# koch <- koch %>%
#   mutate(
#     d18Ocarb_SMOW   = d18ocarb_vsmow,
#     d13Ccarb_VPDB   = d13c_vpdb
#   )

p_d18O_age <- ggplot(subset(koch, !is.na(d18ocarb_vsmow) & !is.na(Age_Ma)),
                     aes(x = d18ocarb_vsmow, y = Age_Ma)) +
  geom_point(color = "black", size = 2.0, alpha = 0.8) +
  scale_y_reverse(limits = age_limits, breaks = age_breaks) +
  labs(
    x = expression(delta^18*O[carb]~("\u2030 VSMOW")),
    y = NULL
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.y = element_blank()
  )

p_d13C_age <- ggplot(subset(koch, !is.na(d13c_vpdb) & !is.na(Age_Ma)),
                     aes(x = d13c_vpdb, y = Age_Ma)) +
  geom_point(color = "black", size = 2.0, alpha = 0.8) +
  scale_y_reverse(limits = age_limits, breaks = age_breaks) +
  labs(
    x = expression(delta^13*C[carb]~("\u2030 VPDB")),
    y = NULL
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.y = element_blank()
  )

panel_D17O_T47_d18O_d13C_age <- p_D17O_age | p_T47_age | p_d18O_age | p_d13C_age
panel_D17O_T47_d18O_d13C_age

ggsave(
  filename = here::here("figures", "panel_D17O_T47_d18O_d13C_age.png"),
  plot = panel_D17O_T47_d18O_d13C_age,
  width = 12, height = 4, dpi = 600
)

plot(koch$t47~koch$d18ocarb_vsmow)
#### PCB Koch d18O ----------------------

kochcarb18O <- ggplot(koch, aes(x = d18o_vpdb, y = Age_Ma)) +
  geom_point(
    data = koch,
    aes(x = d18o_vpdb, y = Age_Ma),
    color = "black",
    size = 1.5,
    alpha = 0.6
  ) +
  # Axes and theme
  scale_y_reverse(
    limits = c(56.5, 55.2),
    breaks = seq(56.5, 55.2, by = -0.1)
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

kochcarb18O


# Save to file
ggsave(filename = here::here("figures", "kochcarb18O.png"),
       plot = kochcarb,
       width = 4, height = 4, dpi = 600)


#### 17O by depth!!!! ------
PCB_17O_depth <- ggplot(IPL_17O, aes(x = mean_D17Ocarb, y = Strat_m_Bowen)) +
  # Error bars (no tips)
  geom_errorbarh(
    aes(xmin = mean_D17Ocarb - sd_D17Ocarb,
        xmax = mean_D17Ocarb + sd_D17Ocarb),
    height = 0,
    na.rm = TRUE,
    color = "gray40"
  ) +
  # Points: n > 1 (black)
  geom_point(
    data = subset(IPL_17O, n > 1),
    color = "black",
    size = 1.5,
    alpha = 0.6
  ) +
  # Points: singletons (red)
  geom_point(
    data = subset(IPL_17O, n == 1),
    color = "red",
    size = 2.5,
    alpha = 1
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

PCB_17O_depth 

# Save to file
ggsave(filename = here::here("figures", "D17O_depth_plot.png"),
       plot = PCB_17O_depth,
       width = 4, height = 4, dpi = 600)

#### 17O PETM only ------ INCOMPLETE CODE

PCB_17O_PETM_depth <- ggplot(IPL_17O, aes(x = mean_D17Ocarb, y = Strat_m_Bowen)) +
  # Error bars (no tips)
  geom_errorbarh(
    aes(xmin = mean_D17Ocarb - sd_D17Ocarb,
        xmax = mean_D17Ocarb + sd_D17Ocarb),
    height = 0,
    na.rm = TRUE,
    color = "gray40"
  ) +
  # Points: n > 1 (black)
  geom_point(
    data = subset(IPL_17O, n > 1),
    color = "black",
    size = 1.5,
    alpha = 0.6
  ) +
  # Points: singletons (red)
  geom_point(
    data = subset(IPL_17O, n == 1),
    color = "red",
    size = 2.5,
    alpha = 1
  ) +
  # Axis and theme
  scale_y_continuous(
    limits = c(1400, 1650),
    breaks = seq(1400, 1650, by = 25)
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

PCB_17O_PETM_depth

# Save to file
ggsave(filename = here::here("figures", "D17O_depth_plot.png"),
       plot = PCB_17O_depth,
       width = 4, height = 4, dpi = 600)

#### Bowen d13C depth --------
bowen_d13C_depth <- ggplot(Bowen_d13C, aes(x = d13C, y = Level)) +
  # Points:(black)
  geom_point(
    data = subset(Bowen_d13C),
    color = "black",
    size = 1.5,
    alpha = 0.6
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

bowen_d13C_depth

# Save to file
ggsave(filename = here::here("figures", "bowen_d13C_depth.png"),
       plot = bowen_d13C_depth,
       width = 4, height = 4, dpi = 600)


#### Koch d13C depth! -----

koch_d13C_depth <- ggplot(koch, aes(x = d13c_vpdb, y = Strat_m_Bowen)) +
  # Points:(black)
  geom_point(
    data = subset(koch),
    color = "black",
    size = 1.5,
    alpha = 0.6
  ) +
  # Axis and theme
  scale_y_continuous(
    limits = c(1200, 1800),
    breaks = seq(1200, 1800, by = 50)
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

koch_d13C_depth

# Save to file
ggsave(filename = here::here("figures", "koch_d13C_depth.png"),
       plot = koch_d13C_depth,
       width = 4, height = 4, dpi = 600)


#### Koch d13C depth! -----

koch_d18O_depth <- ggplot(koch, aes(x = d18o_vpdb, y = Strat_m_Bowen)) +
  # Points:(black)
  geom_point(
    data = subset(koch),
    color = "black",
    size = 1.5,
    alpha = 0.6
  ) +
  # Axis and theme
  scale_y_continuous(
    limits = c(1200, 1800),
    breaks = seq(1200, 1800, by = 50)
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

koch_d18O_depth

# Save to file
ggsave(filename = here::here("figures", "koch_d18O_depth.png"),
       plot = koch_d18O_depth,
       width = 4, height = 4, dpi = 600)


#### Biozone plot -------
library(ggplot2)
library(dplyr)
library(tibble)

# Input depth-based biozone data
biozones <- tribble(
  ~Zone, ~Top_m, ~Bottom_m,
  "Wa-5", 2200, 2200,   # Placeholder — no top defined
  "Wa-4", 2020, 2200,
  "Wa-3b", 1780, 2020,
  "Wa-3a", 1750, 1780,
  "Wa-2", 1645, 1750,
  "Wa-1", 1543, 1645,
  "Wa-0", 1506, 1543,
  "Cf-3", 1335, 1506,
  "Cf-2", 1180, 1335,
  "Cf-1", 885, 1180,
  "Ti-6", 820, 885,
  "Ti-5b", 655, 820,
  "Ti-5a", 530, 655,
  "Ti-4", 415, 530,
  "Ti-3", 215, 415,
  "Ti-2", 155, 215
)

# Filter to only the depth range you want to plot
biozones_trimmed <- biozones %>% 
  filter(Bottom_m >= 1200)

# Plot
biozone_plot <- ggplot(biozones_trimmed) +
  geom_rect(aes(xmin = 0, xmax = 1,
                ymin = Top_m, ymax = Bottom_m,
                fill = Zone),
            color = "black") +
  geom_text(aes(x = 0.5, y = (Top_m + Bottom_m) / 2, label = Zone),
            size = 3, color = "white") +
  scale_y_continuous(
    limits = c(1200, 1800),
    breaks = seq(1200, 1800, by = 50)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

biozone_plot

ggsave(filename = here::here("figures", "biozone_plot.png"),
       plot = biozone_plot,
       width = 4, height = 4, dpi = 600)

### d13C vs D17O ----------

koch <- koch %>% mutate(strat_m_rounded = round(Strat_m_Bowen, 1))
IPL_17O <- IPL_17O %>% mutate(strat_m_rounded = round(Strat_m_Bowen, 1))
koch_17O <- left_join(koch, IPL_17O, by = "strat_m_rounded")


# Join datasets by Strat_m_Bowen
koch_17O <- left_join(koch, IPL_17O, by = "Strat_m_Bowen")

ggplot(koch_17O, aes(x = d13c_vpdb, y = mean_D17Ocarb)) +
  geom_point(size = 2, alpha = 0.8, color = "blue") +
  labs(
    x = expression(delta^13*C[carb]~("\u2030 VPDB")),
    y = expression(Delta*minute^17*O[carb]~("per meg")),
    title = "Triple Oxygen vs δ13C (Koch samples with Δ′17O)"
  ) +
  theme_minimal()

ggplot(koch_17O, aes(x = d18o_vpdb, y = mean_D17Ocarb)) +
  geom_point(size = 2, alpha = 0.8, color = "blue") +
  labs(
    x = expression(delta^18*O[carb]~("\u2030 VPDB")),
    y = expression(Delta*minute^17*O[carb]~("per meg")),
    title = "Triple Oxygen vs δ18O (Koch samples with Δ′17O)"
  ) +
  theme_minimal()


# ==== AGE-BASED PLOTS ====
# p                  – Δ17O vs Age
# D17O_PETM_only     – Δ17O vs Age (PETM only w/ rolling mean)
# marine             – Benthic δ13C vs Age (ODP 1262)
# bwt                – Benthic BWT vs Age
# kochcarb           – Koch δ13C vs Age
# kochcarb18O        – Koch δ18O vs Age

# ==== DEPTH-BASED PLOTS ====
# PCB_17O_depth         – Δ17O vs Strat Depth
# PCB_17O_PETM_depth    – Δ17O vs Depth (PETM interval only)
# bowen_d13C_depth      – Bowen δ13C vs Depth
# koch_d13C_depth       – Koch δ13C vs Depth
# koch_d18O_depth       – Koch δ18O vs Depth
# D47_CU_depth          - CU D47 data

# ==== CROSSPLOTS ====
# (not saved as objects, but could assign these:)
# koch_d13C_vs_D17O     – δ13C vs Δ17O
# koch_d18O_vs_D17O     – δ18O vs Δ17O

### PANEL ASSEMBLY ----

# Load patchwork
library(patchwork)

# Set shared axis theme
shared_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

### PANEL ASSEMBLY ----

library(patchwork)

# Apply shared theme
shared_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

# === AGE-BASED PANEL (Horizontal) ===
panel_age <- (p + shared_theme) |
  (kochcarb + shared_theme) |
  (marine + shared_theme) +
  plot_layout(guides = "collect") &
  scale_y_reverse(limits = c(56.5, 55.2), breaks = seq(56.5, 55.2, -0.1))

ggsave(
  here("figures", "panel_age_combined_horizontal.png"),
  panel_age,
  width = 9, height = 4, dpi = 600
)


# === DEPTH-BASED PANEL (Horizontal) ===
# Define shared reversed y-axis scale (top = lower meter value)
shared_y <- scale_y_continuous(
  limits = c(1380, 1650),
  breaks = seq(1380, 1650, by = 100)
)

# Build panel with shared theme and shared y-scale
panel_depth <- (PCB_17O_depth + shared_theme + shared_y) |
  (D47_CU_depth + shared_theme + shared_y) |
  (bowen_d13C_depth + shared_theme + shared_y) +
  plot_layout(guides = "collect")

# print it 
panel_depth

ggsave(
  here("figures", "panel_depth_combined_horizontal.png"),
  panel_depth,
  width = 9, height = 4, dpi = 600
)