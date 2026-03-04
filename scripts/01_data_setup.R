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

### Define Helper Functions ----------
round_depth <- function(df) {
  df %>% mutate(strat_height_m = round(strat_height_m, 1))
}

### Read Data ----------
samples     <- read.csv(here("data", "raw", "SandCoulee_Polecat_nodules.csv"))
R31_raw     <- read.csv(here("data", "raw", "R31_corData_linear.csv"))
R32_raw     <- read.csv(here("data", "raw", "R32_corData_linear.csv"))
R33_raw     <- read.csv(here("data", "raw", "R33_corrected.csv"))
R34_raw     <- read.csv(here("data", "raw", "R34_corrected.csv"))
R34_rfw     <- read.csv(here("data", "raw", "R34_reconstructed_waters_Excel.csv"))

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
)

### Wrangle IPL Corrected Δ17O Data ----------
R31_PETM <- R31_raw %>%
  filter(Type.2 == "Bighorn PETM") %>%
  mutate(reactor = "R31")

R32_PETM <- R32_raw %>%
  filter(Type.2 == "PETM Bighorn") %>%
  mutate(reactor = "R32")

R33_PETM <- R33_raw %>% mutate(reactor = "R33")
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
IPL_PETM <- summary_stats %>%
  left_join(samples %>% select(Sample_ID, Strat_m_Bowen), by = c("sample.ID" = "Sample_ID")) %>%
  left_join(R34_rfw, by = "sample.ID")

# Save summary
write_csv(IPL_PETM, here("data", "processed", "IPL17O_PETM_summary.csv"))

### Summarize Koch δ13C by strat level ----------
koch_summary <- koch %>%
  group_by(Strat_m_Bowen) %>%
  summarise(
    d13C_carb = mean(d13c_vpdb, na.rm = TRUE),
    n_d13C   = n(),
    .groups  = "drop"
  )

### Summarize Bowen δ13C by strat level ----------
Bowen2001 %>%
  clean_names() %>%
  mutate(level = round(level, 1)) %>%  
  group_by(level) %>%
  summarise(
    d13C_Bowen = mean(d13c_vpdb, na.rm = TRUE),
    n_Bowen    = n(),
    .groups    = "drop"
  ) %>%
  rename(Strat_m_Bowen = level)

names(Bowen2001)



### Build Combined Dataset ----------
combined <- IPL_PETM %>%
  left_join(CU_data %>% select(Strat_m_Bowen, T47_C = t47_c, d18Oc_SMOW), by = "Strat_m_Bowen") %>%
  left_join(koch_summary, by = "Strat_m_Bowen") %>%
  left_join(Bowen2001, by = "Strat_m_Bowen")

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
