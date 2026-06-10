# 01_data_setup.R
# Create processed datasets for PETM/BHB plotting and analysis

### Load libraries ----
library(dplyr)
library(readr)
library(here)
library(isogeochem)
library(stringr)


### Helper functions ----
round_depth <- function(df) {
  df %>%
    mutate(strat_height_m = round(strat_height_m, 1))
}


### Read raw data ----
samples <- read_csv(here("data", "raw", "SandCoulee_Polecat_nodules.csv")) %>%
  rename(sample.ID = Sample_ID)

mod_sample_list <- read_csv(
  here("data", "raw", "SandCoulee_Polecat_nodules_IDflipped.csv")
)

R31_raw <- read_csv(here("data", "raw", "R31_corData_linear.csv"))
R32_raw <- read_csv(here("data", "raw", "R32_corData_linear.csv"))
R33_raw <- read_csv(here("data", "raw", "R33_corrected.csv"))
R34_raw <- read_csv(here("data", "raw", "R34_corrected.csv"))

R34_rfw <- read_csv(
  here("data", "raw", "R34_reconstructed_waters_Excel.csv")
) %>%
  rename(sample.ID = Sample_ID)

CU_data <- read_csv(here("data", "raw", "PETM_clumped.csv")) %>%
  rename(strat_height_m = Strat_m_Bowen) %>%
  mutate(d18Oc_SMOW = to_VSMOW(d18Ocarb_VPBD, eq = "IUPAC")) %>%
  round_depth()

koch <- read_csv(here("data", "raw", "Koch_SC_nodules_isotopes.csv")) %>%
  rename(strat_height_m = Strat_m) %>%
  round_depth()

Bowen2001 <- read_csv(here("data", "raw", "Bowen2001_IsotopeData.csv")) %>%
  rename(strat_height_m = Level) %>%
  mutate(strat_height_m = round(strat_height_m, 1))

matthews_age_model <- read_csv(
  here("data", "raw", "approx_age_mdl_PETM_PCB.csv")
)


### Wrangle IPL clumped isotope data ----
IPL_clumped <- read_csv(
  here("data", "raw", "Matthew's BHB Carbs_5-8-26.csv")
)

names(IPL_clumped) <- str_remove_all(names(IPL_clumped), "'")

IPL_clumped_summary <- IPL_clumped %>%
  group_by(SampleID) %>%
  summarize(
    mean_T_D47_iCDES = mean(`T(D47)iCDES`, na.rm = TRUE),
    se_T_D47_iCDES =
      sd(`T(D47)iCDES`, na.rm = TRUE) /
      sqrt(sum(!is.na(`T(D47)iCDES`))),
    
    mean_TD47iCDES_plus = mean(`TD47iCDES+`, na.rm = TRUE),
    se_TD47iCDES_plus =
      sd(`TD47iCDES+`, na.rm = TRUE) /
      sqrt(sum(!is.na(`TD47iCDES+`))),
    
    mean_T_D47_CDES = mean(`T(D47)CDES`, na.rm = TRUE),
    se_T_D47_CDES =
      sd(`T(D47)CDES`, na.rm = TRUE) /
      sqrt(sum(!is.na(`T(D47)CDES`))),
    
    n = n(),
    .groups = "drop"
  ) %>%
  left_join(
    mod_sample_list %>%
      select(
        SampleID = Sample_ID_flip,
        Strat_m_Bowen
      ),
    by = "SampleID"
  )


### Manual stratigraphic fixes ----
IPL_clumped_summary <- IPL_clumped_summary %>%
  mutate(
    Strat_m_Bowen = case_when(
      SampleID == "PK95-SC-4"     ~ 1570,
      SampleID == "PK95-242"      ~ 590,
      SampleID == "PK95-SC-118up" ~ 1325,
      SampleID == "PK95-SC-27"    ~ 1525,
      TRUE ~ Strat_m_Bowen
    )
  )


### Add age model to IPL clumped summary ----
IPL_clumped_summary <- IPL_clumped_summary %>%
  mutate(
    Age_Ma = approx(
      x = matthews_age_model$est_Depth_m_PCB_outcrop,
      y = matthews_age_model$Age_Ma,
      xout = Strat_m_Bowen,
      rule = 1
    )$y
  )


### Export IPL clumped summary ----
write_csv(
  IPL_clumped_summary,
  here("data", "processed", "IPL_Pg_BHB_D47_summary.csv")
)


### Wrangle IPL corrected Δ17O data ----
R31_PETM <- R31_raw %>%
  filter(Type.2 == "Bighorn PETM") %>%
  mutate(reactor = "R31")

R32_PETM <- R32_raw %>%
  filter(Type.2 == "PETM Bighorn") %>%
  mutate(reactor = "R32")

R33_PETM <- R33_raw %>%
  mutate(
    reactor = "R33",
  )

R34_PETM <- R34_raw %>%
  mutate(reactor = "R34")

IPL_PETM_corrected <- bind_rows(
  R31_PETM,
  R32_PETM,
  R33_PETM,
  R34_PETM
)

IPL17O_summary <- IPL_PETM_corrected %>%
  group_by(sample.ID) %>%
  summarize(
    mean_D17Ocarb = mean(D17O.SMOWSLAP.per.meg.carbNorm, na.rm = TRUE),
    sd_D17Ocarb   = sd(D17O.SMOWSLAP.per.meg.carbNorm, na.rm = TRUE),
    se_D17Ocarb   =
      sd(D17O.SMOWSLAP.per.meg.carbNorm, na.rm = TRUE) /
      sqrt(sum(!is.na(D17O.SMOWSLAP.per.meg.carbNorm))),
    n = sum(!is.na(D17O.SMOWSLAP.per.meg.carbNorm)),
    .groups = "drop"
  ) %>%
  left_join(samples, by = "sample.ID") %>%
  rename(strat_height_m = Strat_m_Bowen)


### Export IPL Δ17O summary ----
write_csv(
  IPL17O_summary,
  here("data", "processed", "IPL17O_PETM_summary.csv")
)


### Summarize Koch δ13C by stratigraphic level ----
koch_summary <- koch %>%
  group_by(strat_height_m) %>%
  summarize(
    d13C_carb = mean(d13C_VPDB, na.rm = TRUE),
    n_d13C = n(),
    .groups = "drop"
  )


### Summarize Bowen 2001 δ13C by stratigraphic level ----
Bowen_summary <- Bowen2001 %>%
  group_by(strat_height_m) %>%
  summarize(
    d13C_Bowen = mean(d13C_VPDB, na.rm = TRUE),
    n_Bowen = n(),
    .groups = "drop"
  )


### Build combined plotting dataset ----
combined <- IPL17O_summary %>%
  left_join(
    CU_data %>%
      select(strat_height_m, T47_C, T47_2SE_C, d18Oc_SMOW),
    by = "strat_height_m"
  ) %>%
  left_join(koch_summary, by = "strat_height_m") %>%
  left_join(Bowen_summary, by = "strat_height_m")


### Export combined dataset ----
write_csv(
  combined,
  here("data", "processed", "PETM_combined.csv")
)

write_csv(
  combined,
  here("data", "processed", "combined_petm_data.csv")
)


### Optional coverage check ----
combined %>%
  summarize(
    n_rows = n(),
    n_D17O = sum(!is.na(mean_D17Ocarb)),
    n_T47 = sum(!is.na(T47_C)),
    n_d13C_Koch = sum(!is.na(d13C_carb)),
    n_d13C_Bowen = sum(!is.na(d13C_Bowen))
  ) %>%
  print()