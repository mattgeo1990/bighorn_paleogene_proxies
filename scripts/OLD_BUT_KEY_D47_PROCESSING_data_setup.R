# 01_data_setup.R
# Create combined dataset for multipanel plotting

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

### Define Helper Functions ----------
round_depth <- function(df) {
  df %>% mutate(strat_height_m = round(strat_height_m, 1))
}

### Read Data ----------
samples     <- read.csv(here("data", "raw", "SandCoulee_Polecat_nodules.csv")) %>%
  rename(sample.ID = Sample_ID) 
R31_raw     <- read.csv(here("data", "raw", "R31_corData_linear.csv"))
R32_raw     <- read.csv(here("data", "raw", "R32_corData_linear.csv"))
R33_raw     <- read.csv(here("data", "raw", "R33_corrected.csv"))
R34_raw     <- read.csv(here("data", "raw", "R34_corrected.csv"))
R34_rfw     <- read.csv(here("data", "raw", "R34_reconstructed_waters_Excel.csv")) %>%
  rename(sample.ID = Sample_ID)

CU_data <- read_csv(here("data", "raw", "PETM_clumped.csv")) %>%
  #clean_names() %>%
  rename(strat_height_m = Strat_m_Bowen) %>%
  mutate(d18Oc_SMOW = to_VSMOW(d18Ocarb_VPBD, eq = "IUPAC")) %>%
  round_depth()

koch <- read_csv(here("data", "raw", "Koch_SC_nodules_isotopes.csv")) %>%
  #clean_names() %>%
  rename(strat_height_m = Strat_m) %>%
  round_depth()

Bowen2001 <- read.csv(
  here("data", "raw", "Bowen2001_IsotopeData.csv")
) %>%
  rename(strat_height_m = Level)






###Wrangle IPL Corrected IPL Clumped Data###


# Import corrected IPL clumped isotope data
IPL_clumped <- read_csv(
  here("data", "raw", "Matthew's BHB Carbs_5-8-26.csv")
) # these data include manual corrections by Ben Passey. These data are reliable

# Remove apostrophes introduced during import
names(IPL_clumped) <- str_remove_all(names(IPL_clumped), "'")

# Inspect imported data
table(IPL_clumped$SampleID)
str(IPL_clumped)



### Plot all replicates and identify 2025-2026 clumped session 1 vs clumped session 2
# Ben noted that session 1 seems to have produced higher temp values

# read in Ben's May 13 2026 summary (corrected data with strat levels and session #)

BP_May13_summary <- read.csv(here::here("data", "excel files", "BHB Paleogene Summary May 2026.csv"))

str(BP_May13_summary)

# Plot "T.D47..Petersen" vs "Strat". Color points based on "Session".
ggplot(BP_May13_summary, aes(x = T.D47..Petersen, y = Strat, color = Session)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Strat",
    color = "Session"
  ) +
  theme_classic()

# Average duplicate analyses within Sample × Session,
# keep only samples measured in both sessions

paired <- BP_May13_summary %>%
  filter(
    !is.na(Session),
    Session %in% c("Session 1", "Session 2"),
    !is.na(T.D47..Petersen),
    !is.na(Strat)
  ) %>%
  group_by(Strat, Session) %>%
  summarise(
    T_mean = mean(T.D47..Petersen, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Session,
    values_from = T_mean
  ) %>%
  filter(
    !is.na(`Session 1`),
    !is.na(`Session 2`)
  ) %>%
  mutate(
    dT = `Session 2` - `Session 1`
  )

print(paired)

ggplot(paired, aes(x = `Session 1`, y = `Session 2`)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_point(size = 3) +
  coord_equal() +
  labs(
    x = expression("Session 1 T" [Delta47] * " (" * degree * "C)"),
    y = expression("Session 2 T" [Delta47] * " (" * degree * "C)")
  ) +
  theme_classic()

paired %>%
  summarise(
    n_pairs = n(),
    mean_dT = mean(dT),
    median_dT = median(dT),
    sd_dT = sd(dT)
  )

ggplot(paired, aes(x = Strat, y = dT)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(size = 3) +
  coord_flip() +
  labs(
    x = "Strat",
    y = expression(Delta*"T = Session 2 - Session 1 ("*degree*"C)")
  ) +
  theme_classic()

# SC-27 is the only one with significant offest between session 1 mean and session 2 mean
# let's look at the data
BP_May13_summary$T.D47..Petersen[which(BP_May13_summary$Strat %in% "1525")]

# which strat levels only have data from one of the two sessions?
BP_May13_summary %>%
  filter(
    !is.na(Session),
    Session %in% c("Session 1", "Session 2")
  ) %>%
  group_by(Strat) %>%
  summarise(
    sessions = paste(sort(unique(Session)), collapse = ", "),
    n_sessions = n_distinct(Session),
    .groups = "drop"
  ) %>%
  filter(n_sessions == 1) %>%
  arrange(Strat)





### Prepare to summarize the data by strat level 

# Import sample metadata with stratigraphic positions
mod_sample_list <- read_csv(
  here("data", "raw", "SandCoulee_Polecat_nodules_IDflipped.csv")
)

# Summarize replicate-level IPL temperatures by sample
IPL_clumped_summary <- IPL_clumped %>%
  group_by(SampleID) %>%
  summarize(
    
    mean_T_D47_iCDES = mean(`T(D47)iCDES`, na.rm = TRUE),
    se_T_D47_iCDES = sd(`T(D47)iCDES`, na.rm = TRUE) /
      sqrt(sum(!is.na(`T(D47)iCDES`))),
    
    mean_TD47iCDES_plus = mean(`TD47iCDES+`, na.rm = TRUE),
    se_TD47iCDES_plus = sd(`TD47iCDES+`, na.rm = TRUE) /
      sqrt(sum(!is.na(`TD47iCDES+`))),
    
    mean_T_D47_CDES = mean(`T(D47)CDES`, na.rm = TRUE),
    se_T_D47_CDES = sd(`T(D47)CDES`, na.rm = TRUE) /
      sqrt(sum(!is.na(`T(D47)CDES`))),
    
    n = n(),
    .groups = "drop"
  )

# Add stratigraphic heights
IPL_clumped_summary <- IPL_clumped_summary %>%
  left_join(
    mod_sample_list %>%
      select(
        SampleID = Sample_ID_flip,
        Strat_m_Bowen
      ),
    by = "SampleID"
  )


# Some are missing, so here's a stop-gap measure
IPL_clumped_summary$Strat_m_Bowen[which(IPL_clumped_summary$SampleID %in% "PK95-SC-4")] <- 1570
IPL_clumped_summary$Strat_m_Bowen[which(IPL_clumped_summary$SampleID %in% "PK95-242")] <- 590
IPL_clumped_summary$Strat_m_Bowen[which(IPL_clumped_summary$SampleID %in% "PK95-SC-118up")] <- 1325
IPL_clumped_summary$Strat_m_Bowen[which(IPL_clumped_summary$SampleID %in% "PK95-SC-27")] <- 1525.0
IPL_clumped_summary$Strat_m_Bowen[which(IPL_clumped_summary$SampleID %in% "PK95-SC-27")] <- 1525.0


# EXPORT csv 
write_csv(IPL_clumped_summary, here::here("data", "processed", "IPL_Pg_BHB_D47_summary.csv"))



# IPL data: iCDESplus only, with 1 SE errors
ipl_plot_df <- IPL_clumped_summary %>%
  transmute(
    SampleID,
    Strat_m_Bowen,
    Temperature = mean_TD47iCDES_plus,
    Error = se_TD47iCDES_plus,
    Dataset = "IPL"
  ) %>%
  filter(!is.na(Strat_m_Bowen), !is.na(Temperature))

# CU clumped isotope temperatures for comparison
# T47_2SE_C is treated here as 2 SE
cu_plot_df <- CU_data %>%
  transmute(
    SampleID = Sample_ID,
    Strat_m_Bowen = strat_height_m,
    Temperature = T47_C,
    Error = T47_2SE_C,
    Dataset = "CU"
  ) %>%
  filter(!is.na(Strat_m_Bowen), !is.na(Temperature))

ggplot() +
  
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 1500,
    ymax = 1550,
    fill = "red",
    alpha = 0.12
  ) +
  
  geom_errorbarh(
    data = cu_plot_df,
    aes(
      xmin = Temperature - Error,
      xmax = Temperature + Error,
      y = Strat_m_Bowen,
      color = Dataset
    ),
    height = 0,
    linewidth = 0.35,
    alpha = 0.55,
    na.rm = TRUE
  ) +
  
  geom_point(
    data = cu_plot_df,
    aes(
      x = Temperature,
      y = Strat_m_Bowen,
      color = Dataset
    ),
    size = 1.8,
    shape = 1,
    stroke = 0.8,
    alpha = 0.9
  ) +
  
  geom_errorbarh(
    data = ipl_plot_df,
    aes(
      xmin = Temperature - Error,
      xmax = Temperature + Error,
      y = Strat_m_Bowen,
      color = Dataset
    ),
    height = 0,
    linewidth = 0.45,
    alpha = 0.7,
    na.rm = TRUE
  ) +
  
  geom_point(
    data = ipl_plot_df,
    aes(
      x = Temperature,
      y = Strat_m_Bowen,
      color = Dataset
    ),
    size = 2.6,
    alpha = 0.9
  ) +
  
  scale_x_continuous(
    limits = c(0, 90),
    breaks = seq(0, 90, 10),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  
  labs(
    x = expression(Delta[47]~"temperature (°C)"),
    y = "Stratigraphic height (m)",
    color = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.45),
    axis.ticks.length = unit(2.5, "mm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 11),
    legend.position = c(0.85, 0.15),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 9),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

ggsave(
  here::here("figures", "clumped_temperature_plot.png"),
  width = 7,
  height = 8
)






# Combine PB-only IPL and CU temperatures for sample-level comparison
jitter_df <- bind_rows(
  
  ipl_plot_df %>%
    filter(str_detect(SampleID, "^PB")) %>%
    select(
      SampleID,
      Calibration,
      Temperature,
      Error,
      Dataset
    ),
  
  cu_plot_df %>%
    filter(str_detect(SampleID, "^PB")) %>%
    select(
      SampleID,
      Calibration,
      Temperature,
      Error,
      Dataset
    )
)

# Use the same jitter position for points and error bars
jitter_pos <- position_jitter(
  width = 0.12,
  height = 0,
  seed = 1
)

# Jitter plot comparing CU and IPL temperatures for PB samples
ggplot(
  jitter_df,
  aes(
    x = SampleID,
    y = Temperature,
    color = Calibration
  )
) +
  
  geom_errorbar(
    aes(
      ymin = Temperature - Error,
      ymax = Temperature + Error
    ),
    position = jitter_pos,
    width = 0,
    linewidth = 0.35,
    alpha = 0.6,
    na.rm = TRUE
  ) +
  
  geom_point(
    aes(shape = Dataset),
    position = jitter_pos,
    size = 1.5,
    alpha = 0.75
  ) +
  
  scale_shape_manual(
    values = c(
      "IPL" = 16,
      "CU" = 1
    )
  ) +
  
  labs(
    x = "Sample",
    y = "Temperature (°C)",
    color = "Calibration",
    shape = "Dataset"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

## clumped spar mixing model -----
library(tibble)

# Endmember temperatures
T_spar <- 80
T_micrite <- 30

# Total analyzed powder mass (mg)
total_mass_mg <- 10

# Target bulk temperature
target_T <- 40

# Simple mass-based mixing model
mixing_df <- tibble(
  spar_mg = seq(0, total_mass_mg, by = 0.1),
  micrite_mg = total_mass_mg - spar_mg,
  f_spar = spar_mg / total_mass_mg,
  percent_spar = f_spar * 100,
  bulk_T = f_spar * T_spar + (1 - f_spar) * T_micrite
)

# Calculate spar required for target bulk temperature
f_spar_needed <- (target_T - T_micrite) / (T_spar - T_micrite)

spar_mg_needed <- f_spar_needed * total_mass_mg
micrite_mg_needed <- total_mass_mg - spar_mg_needed

# Print results
cat("Target bulk temperature:", target_T, "°C\n")
cat("Required spar mass:", round(spar_mg_needed, 2), "mg\n")
cat("Required micrite mass:", round(micrite_mg_needed, 2), "mg\n")
cat("Percent spar:", round(f_spar_needed * 100, 1), "%\n")

# View full mixing table
mixing_df


library(tibble)
library(ggplot2)

# Endmember temperatures
T_spar <- 80
T_micrite <- 30

# Total analyzed powder mass (mg)
total_mass_mg <- 10

# Mixing model
mixing_df <- tibble(
  spar_mg = seq(0, total_mass_mg, by = 0.1),
  f_spar = spar_mg / total_mass_mg,
  percent_spar = f_spar * 100,
  bulk_T = f_spar * T_spar + (1 - f_spar) * T_micrite
)

# Crossplot
ggplot(mixing_df, aes(x = percent_spar, y = bulk_T)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Spar (%)",
    y = "Bulk Temperature (°C)"
  ) +
  theme_classic(base_size = 12)


## its lready a mess lets throw in a random age plot here -----

matthews_age_model <- read.csv(here("data", "raw", "approx_age_mdl_PETM_PCB.csv"))

str(matthews_age_model)

# Interpolate ages for IPL_clumped_summary depths using the age model
# Treat est_Depth_m_PCB_outcrop as equivalent to Strat_m_Bowen

IPL_clumped_summary <- IPL_clumped_summary %>%
  mutate(
    Age_Ma = approx(
      x = matthews_age_model$est_Depth_m_PCB_outcrop,
      y = matthews_age_model$Age_Ma,
      xout = Strat_m_Bowen,
      rule = 1
    )$y
  )

# Inspect results
IPL_clumped_summary %>%
  select(SampleID, Strat_m_Bowen, Age_Ma)

### Wrangle IPL Corrected Δ17O Data ----------
R31_PETM <- R31_raw %>%
  filter(Type.2 == "Bighorn PETM") %>%
  mutate(reactor = "R31")

R32_PETM <- R32_raw %>%
  filter(Type.2 == "PETM Bighorn") %>%
  mutate(reactor = "R32")

R33_PETM <- R33_raw %>% mutate(reactor = "R33") %>%
  mutate(Date.Time = as.character.Date(Date.Time))

R34_PETM <- R34_raw %>% mutate(reactor = "R34")

IPL_PETM_corrected <- bind_rows(R31_PETM, R32_PETM, R33_PETM, R34_PETM)

# Compute summary stats per sample.ID
summary_stats <- IPL_PETM_corrected %>%
  group_by(sample.ID) %>%
  summarise(
    mean_D17Ocarb = mean(D17O.SMOWSLAP.per.meg.carbNorm, na.rm = TRUE),
    sd_D17Ocarb   = sd(D17O.SMOWSLAP.per.meg.carbNorm, na.rm = TRUE),
    se_D17Ocarb   = sd(D17O.SMOWSLAP.per.meg.carbNorm, na.rm = TRUE) / sqrt(n()),
    n             = n(),
    .groups       = "drop"
  )

# Add strat and rfw metadata
##R34_rfw doesn't have any matching sample IDS, 
##causing a bug where all the data goes away
IPL_PETM <- summary_stats %>%
  left_join(IPL_PETM, samples %>% select(sample.ID, Strat_m_Bowen), 
            by = "sample.ID") %>%
  left_join(IPL_PETM, R34_rfw, by = "sample.ID")

##Will's version
IPL_PETM <- summary_stats
left_join(IPL_PETM, samples, by = "sample.ID")

IPL_PETM2 = merge(IPL_PETM, samples, by = "sample.ID", all.x = TRUE) %>%
  rename(strat_height_m = Strat_m_Bowen)

# Save summary
write.csv(IPL_PETM2, here("data", "processed", "IPL17O_PETM_summary.csv"))

### Summarize Koch δ13C by strat level ----------
koch_summary <- koch %>%
  group_by(strat_height_m) %>%
  summarise(
    d13C_carb = mean(d13C_VPDB, na.rm = TRUE),
    n_d13C   = n(),
    .groups  = "drop"
  )

### Summarize Bowen δ13C by strat level ----------
Bowen2001 %>%
  #clean_names() %>%
  mutate(strat_height_m = round(strat_height_m, 1)) %>%  
  group_by(strat_height_m) %>%
  summarise(
    d13C_Bowen = mean(d13C_VPDB, na.rm = TRUE),
    n_Bowen    = n(),
    .groups    = "drop")

names(Bowen2001)



### Build Combined Dataset ----------
combined <- IPL_PETM2 %>%
  left_join(CU_data %>% select(strat_height_m, T47_C, d18Oc_SMOW), by = "strat_height_m") %>%
  left_join(koch_summary, by = "strat_height_m") %>%
  left_join(Bowen2001, by = "strat_height_m")

# Optional: write combined to file
write_csv(combined, here("data", "processed", "PETM_combined.csv"))


# Optional: check coverage of each variable
#combined %>%
#  summarise(
#    n_rows = n(),
 #   n_D17O = sum(!is.na(mean_D17Ocarb)),
#    n_T47 = sum(!is.na(T47_C)),
#    n_d13C = sum(!is.na(d13C_carb))
#  ) %>%
#  print()

# Export cleaned dataset
write_csv(combined, here("data", "processed", "combined_petm_data.csv"))



#### Emergency panel code --------

library(ggplot2)
library(patchwork)
library(dplyr)

# Common y-axis limits/breaks
y_lims <- c(500, 1800)
y_breaks <- seq(500, 1800, by = 100)

petm_rect <- annotate(
  "rect",
  xmin = -Inf, xmax = Inf,
  ymin = 1500, ymax = 1545,
  fill = "red", alpha = 0.15
)

panel_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank()
  )

# ---- D17O water / reconstructed fluid water ----
# Change x variable here if your column name is slightly different
p_D17Orfw <- ggplot(df, aes(x = D17Orfw, y = strat_height_m)) +
  petm_rect +
  geom_point(size = 2.5, na.rm = TRUE) +
  geom_path(na.rm = TRUE) +
  scale_y_continuous(limits = y_lims, breaks = y_breaks) +
  labs(
    x = expression(Delta^17*O[rfw]~"(per meg)"),
    y = "Stratigraphic height (m)",
    title = expression(Delta^17*O[rfw])
  ) +
  panel_theme +
  theme(axis.title.y = element_text())

# ---- D47 temperature ----
p_D47 <- ggplot(BP_May2026_strat_summary,
                aes(x = mean_T47, y = Strat)) +
  petm_rect +
  geom_errorbarh(
    aes(xmin = mean_T47 - sd_T47,
        xmax = mean_T47 + sd_T47),
    height = 0,
    na.rm = TRUE
  ) +
  geom_point(size = 2.5, na.rm = TRUE) +
  scale_x_continuous(breaks = seq(15, 60, by = 5)) +
  scale_y_continuous(limits = y_lims, breaks = y_breaks) +
  labs(
    x = expression(T[Delta47]~"("*degree*C*")"),
    title = expression(Delta[47])
  ) +
  panel_theme

# ---- d18O carb ----
p_d18O <- ggplot(df, aes(x = d18Oc_SMOW, y = strat_height_m)) +
  petm_rect +
  geom_point(size = 2.5, na.rm = TRUE) +
  geom_path(na.rm = TRUE) +
  scale_y_continuous(limits = y_lims, breaks = y_breaks) +
  labs(
    x = expression(delta^18*O[carb]~"\u2030 VSMOW"),
    title = expression(delta^18*O[carb])
  ) +
  panel_theme

# ---- d13C carb ----
p_d13C <- ggplot(df, aes(x = d13C_carb, y = strat_height_m)) +
  petm_rect +
  geom_point(size = 2.5, na.rm = TRUE) +
  geom_path(na.rm = TRUE) +
  scale_y_continuous(limits = y_lims, breaks = y_breaks) +
  labs(
    x = expression(delta^13*C[carb]~"\u2030 VPDB"),
    title = expression(delta^13*C[carb])
  ) +
  panel_theme

# ---- combine ----
final_panel <- p_D17Orfw + p_D47 + p_d18O + p_d13C +
  plot_layout(ncol = 4)

final_panel

ggsave(
  here::here("figures", "PETM_strat_proxy_panel_quick.png"),
  final_panel,
  width = 13,
  height = 7,
  dpi = 600
)
