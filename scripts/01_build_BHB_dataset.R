
# 01_build_BHB_dataset.R
# Purpose: Load, clean, summarize, and combine BHB proxy datasets
# Outputs: Clean proxy-specific summaries + BHB multiproxy master by strat


# 1. Setup ----------------------------------------------------
# Load required packages

library(tidyverse)
library(here)
library(isogeochem)
library(plotly)

# source(here::here("scripts", "00_setup.R"))

# Helper function to round stratigraphic depths to the nearest 0.1 m
round_depth <- function(df) {
  df %>% mutate(strat_height_m = round(strat_height_m, 1))
}

# Helper function to check for multiple sample IDs at the same stratigraphic level
#
# Purpose:
# Identify stratigraphic levels that contain more than one unique
# horizon/sample identifier. This is useful before merging datasets
# by stratigraphic height because multiple IDs at the same level may
# require averaging, selection of a representative sample, or special
# handling during joins.

check_levels <- function(df, strat_col, id_col) {
  
  df %>%
    group_by({{ strat_col }}) %>%
    summarise(
      n_ids = n_distinct({{ id_col }}),
      ids = paste(sort(unique({{ id_col }})), collapse = ", "),
      .groups = "drop"
    ) %>%
    filter(n_ids > 1) %>%
    arrange(desc(n_ids))
}


# 2. Load raw data --------------------------------------------

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data model notes
#
# MLA_sample_id:
#   Project-standard sample identifier.
#   For IPL datasets, this refers to a powder analyzed in the IPL.
#   For external datasets, this is my standardized name for the
#   individual sample/entity reported by the original study.
#
# MLA_horizon_id:
#   Project-standard stratigraphic horizon identifier.
#   Multiple samples may belong to one horizon, and multiple horizons
#   may occur within one paleosol or locality.
#
# Important:
#   MLA_sample_id and MLA_horizon_id are not always one-to-one.
#   Source-specific spreadsheets/text should be consulted when distinctions
#   among locality, paleosol, horizon, nodule, and powder matter.

# Stratigraphic framework:
#   Unless otherwise noted, strat_height_m is meters above the K–Pg boundary
#   in the northern Bighorn Basin composite section.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load IPL sample list (PB and SC samples)

IPL_sample_list <- read.csv(
  here("data", "raw", "SandCoulee_Polecat_nodules.csv")
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IPL triple oxygen isotope dataset
#
# Compiled project spreadsheet integrating IPL Δ′17O data from multiple
# reactors. Reactor 31/32 data were reduced with the IPL R pipeline;
# Reactor 33/34 data were manually corrected/standardized in Excel
# by Matthew Allen using the same reporting framework.

IPL_D17O_data <- read.csv(
  here("data", "raw", "all_data_pre-June2026_PaleogeneBHB_IPL17O_standardized_columns.csv")
)

names(IPL_D17O_data)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IPL clumped isotope (Δ47) dataset
#
# Cleaned, corrected, and simplified IPL Δ47 spreadsheet prepared by
# Ben Passey for Bighorn Basin Paleogene carbonate samples.

IPL_D47_data <- read.csv(
  here("data", "raw", "IPL_D47_BHB_Pg_Summary_June2026.csv")
)

names(IPL_D47_data)
table(IPL_D47_data$MLA_horizon_id)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CU Boulder Δ47 dataset
#
# Havranek/CU Boulder data are already summarized at the sample level.
# Retain raw-style CU_data for inspection, then create CU_summary with
# source-prefixed column names for downstream merging.

CU_data <- read_csv(
  here("data", "raw", "PETM_clumped.csv")
) %>%
  mutate(
    d18Oc_SMOW = to_VSMOW(d18Ocarb_VPBD, eq = "IUPAC")
  ) %>%
  round_depth()

CU_summary <- CU_data %>%
  transmute(
    CU_sample_id = CU_Sample_ID,
    MLA_sample_id,
    MLA_horizon_id,
    strat_height_m,
    
    CU_mean_d13Ccarb_vpdb   = d13C_VPBD,
    CU_mean_d18Ocarb_vpdb   = d18Ocarb_VPBD,
    CU_se_d18Ocarb_vpdb     = d18O_SE_VPBD,
    CU_mean_d18Ocarb_vsmow  = d18Oc_SMOW,
    
    CU_mean_D47_ICDES       = D47_ICDES,
    CU_se_D47_ICDES         = D47_SE_ICDES,
    CU_mean_T47_C           = T47_C,
    CU_2se_T47_C            = T47_2SE_C,
    
    CU_mean_d18Ow_vsmow     = d18Ow_VSMOW,
    CU_se_d18Ow_vsmow       = d18Ow_SE_VSMOW
  )


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Snell et al. (2013) data are already summarized at the sample level.
# Create a source-prefixed summary object for downstream merging while
# preserving provenance.


Snell2013 <- read_csv(
  here("data", "raw", "SnellEtAl2013_summary.csv")
) %>%
  round_depth()

Snell2013_summary <- Snell2013 %>%
  transmute(
    Snell_sample_id = Snell_Sample_ID,
    MLA_sample_id,
    MLA_horizon_id,
    strat_height_m,
    Snell_Age_Ma = Age_Ma,
    Snell_sample_type = Sample_Type,
    
    Snell_mean_d18Ow_vsmow   = Average_d18Ow_permil_SMOW,
    Snell_mean_d18Ocarb_vsmow = Average_d18Oc_permil_SMOW,
    Snell_mean_d13Ccarb_vpdb  = Average_d13Cc_permil_PDB,
    
    Snell_mean_D47 = Average_D47_permil,
    Snell_se_D47   = D47_1se,
    
    Snell_mean_T47_C = Average_Temp_C,
    Snell_se_T47_C   = Temp_1se
  )


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Paul Koch (1992, 1995) carbonate isotope dataset
#
# Replicate-level paleosol carbonate δ13C and δ18O data.
# Will be summarized by MLA_horizon_id before merging.

koch <- read_csv(
  here("data", "raw", "Koch_SC_nodules_isotopes.csv")
) %>%
  mutate(
    Koch_d13Ccarb_vpdb  = d13C_VPDB,
    Koch_d18Ocarb_vpdb  = d18Ocarb_VPDB,
    Koch_d18Ocarb_vsmow = d18Ocarb_VSMOW
  ) %>%
  select(
    Koch_Sample,
    Locality,
    MLA_sample_id,
    MLA_horizon_id,
    strat_height_m,
    Koch_d13Ccarb_vpdb,
    Koch_d18Ocarb_vpdb,
    Koch_d18Ocarb_vsmow
  ) %>%
  round_depth()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bowen et al. (2001) carbonate isotope dataset
#
# Replicate-level Polecat Bench paleosol carbonate δ13C and δ18O data.
# Will be screened for altered values and summarized by MLA_horizon_id
# before merging.

bowen <- read_csv(
  here("data", "raw", "Bowen2001_IsotopeData.csv")
) %>%
  mutate(
    Bowen_d13Ccarb_vpdb     = as.numeric(d13C_VPDB),
    Bowen_d18Ocarb_vpdb     = as.numeric(d18Ocarb_VPDB),
    Bowen_d18Ocarb_vsmow    = as.numeric(d18Ocarb_VSMOW)
  ) %>%
  select(
    Bowen_Soil_ID,
    Bowen_Sample_ID,
    MLA_sample_id,
    MLA_horizon_id,
    strat_height_m,
    Bowen_d13Ccarb_vpdb,
    Bowen_d18Ocarb_vpdb,
    Bowen_d18Ocarb_vsmow
  ) %>%
  round_depth()


# 3. Clean IPL Δ′17O data -------------------------------------

# check for mismatch issues or outliers
# Plot a histogram of mismatch values to visualize the overall distribution
hist(IPL_D17O_data$X33_mismatch,
     main = "Histogram of X33 Mismatch",
     xlab = "X33 Mismatch",
     col = "skyblue",
     border = "white")

# Scatterplot to examine relationship between mismatch and Δ′17O values

p <- ggplot(
  IPL_D17O_data,
  aes(
    x = Dp17Ocarb_permeg_final_correction,
    y = X33_mismatch,
    text = paste(
      "IPL_num:", IPL_num,
      "<br>Sample:", MLA_sample_id,
      "<br>Horizon:", MLA_horizon_id,
      "<br>Δ′17O:", round(Dp17Ocarb_permeg_final_correction, 1),
      "<br>Mismatch:", round(X33_mismatch, 3)
    )
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.8,
    color = "gray40"
  ) +
  labs(
    title = "Mismatch vs Δ′17Ocarb",
    x = "Δ′17Ocarb (per meg)",
    y = "X33 Mismatch"
  ) +
  theme_classic()

ggplotly(
  p,
  tooltip = "text"
)
# Manually exclude analyses with anomalously high X33 mismatch.
# These exclusions are based on analytical QC, not on Δ′17Ocarb value.
high_mismatch <- c("5699", "5841")

# Inspect excluded high-mismatch analyses.
IPL_D17O_data %>%
  filter(IPL_num %in% high_mismatch) %>%
  select(
    IPL_num,
    MLA_sample_id,
    MLA_horizon_id,
    Dp17Ocarb_permeg_final_correction,
    X33_mismatch
  ) %>%
  arrange(desc(abs(X33_mismatch)))

# Remove high-mismatch analyses before summary statistics.
IPL_D17O_data_clean <- IPL_D17O_data %>%
  filter(!IPL_num %in% high_mismatch)

# Examine relationship between analytical mismatch and Δ′17Ocarb.
plot(
  IPL_D17O_data_clean$X33_mismatch ~
    IPL_D17O_data_clean$Dp17Ocarb_permeg_final_correction,
  main = "Mismatch vs Δ′17Ocarb",
  xlab = expression(Delta * minute^17 * O[carb] ~ "(per meg)"),
  ylab = "X33 Mismatch",
  pch = 19,
  col = "gray40"
)

# Identify statistical Δ′17Ocarb outliers for inspection only.
# These are not removed automatically because statistical outlier status
# alone is not sufficient evidence of analytical failure.
outliers <- boxplot.stats(
  IPL_D17O_data_clean$Dp17Ocarb_permeg_final_correction
)$out

# Inspect analyses flagged by the 1.5×IQR rule.
IPL_D17O_data_clean %>%
  filter(Dp17Ocarb_permeg_final_correction %in% outliers) %>%
  select(
    IPL_num,
    MLA_sample_id,
    MLA_horizon_id,
    Dp17Ocarb_permeg_final_correction,
    X33_mismatch
  ) %>%
  arrange(Dp17Ocarb_permeg_final_correction)

# Retain SC-49 and SC-50 for now.
# SC-49 has two low-mismatch analyses after removal of the high-mismatch run,
# and SC-50 is a single extreme value with no independent QC reason for exclusion.
# Additional SC-50 replicates are pending.

# Inspect PB-00-02-09L, which shows poor replicate agreement.
IPL_D17O_data_clean %>%
  filter(MLA_horizon_id == "PB-00-02-09L") %>%
  select(
    IPL_num,
    MLA_sample_id,
    MLA_horizon_id,
    Dp17Ocarb_permeg_final_correction,
    X33_mismatch
  )

# Temporarily omit the high Δ′17Ocarb replicate from PB-00-02-09L.
# The paired replicate is ~-120 per meg, so this sample requires additional
# replication before the less negative value is treated as representative.
omit_low_confidence <- c("5780")

IPL_D17O_data_final <- IPL_D17O_data_clean %>%
  filter(!IPL_num %in% omit_low_confidence)

# Exclusion justification:
# Analyses were excluded only when independent QC information supported removal.
# High-mismatch analyses were removed because anomalous X33 mismatch indicates
# reduced confidence in the analytical correction. Statistical Δ′17Ocarb outliers
# identified by the 1.5×IQR rule were inspected but not automatically removed.
# SC-49 and SC-50 were retained because their unusual Δ′17Ocarb values are not
# accompanied by clear analytical failure. One PB-00-02-09L replicate was
# temporarily omitted because the two-replicate spread is large and additional
# replication is needed to determine which value best represents the sample.
# All excluded analyses remain archived in the raw dataset.
# Remove singletons ?
# IPL17O_summary <- IPL17O_summary %>% filter(n > 1)


# 4. Clean IPL Δ47 data ---------------------------------------

# Count analyses per Sample.ID to identify special qualifiers in sample names (e.g.,"SPAR")
table(IPL_D47_data$MLA_sample_id)

# Extract analyses of sparitic calcite 
# These samples contain "SPAR" in the Sample.ID and are
# treated separately from the primary paleosol dataset.
IPL_D47_SPAR_data <- IPL_D47_data %>%
  filter(grepl("SPAR", MLA_sample_id))

# Extract analyses of micritic or microsparitic calcite
# Any sample containing "SPAR" is excluded.
IPL_D47_primary_data <- IPL_D47_data %>%
  filter(!grepl("SPAR", MLA_sample_id))

# Plot replicate-level Δ47-derived temperatures versus stratigraphic position
#
# Purpose:
# Visualize all individual Δ47 temperature replicates from primary (non-SPAR)
# pedogenic carbonate samples to assess temperature variability with stratigraphy
# and identify potential outliers or stratigraphic trends.
ggplot(IPL_D47_primary_data,
       aes(x = T47_preferred,
           y = strat_height_m)) +
  geom_point(
    size = 2,
    alpha = 0.7,
    color = "gray40"
  ) +
  labs(
    x = expression("Temperature from " * Delta[47] * " (" * degree * "C)"),
    y = "Stratigraphic Level (m)",
    title = expression("Replicate-Level " * T[47] * " Values")
  ) +
  theme_classic()


# Plot Δ47 temperatures by sample ID
#
# Purpose:
# Evaluate within-sample reproducibility and identify samples with unusually
# high or low replicate temperatures. Useful for spotting potential outliers,
# mixed populations, or problematic analyses before calculating sample means.
ggplot(IPL_D47_data,
       aes(
         x = reorder(
           MLA_sample_id,
           T47_preferred,
           FUN = median,
           na.rm = TRUE
         ),
         y = T47_preferred,
         color = Session
       )) +
  
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8
  ) +
  
  coord_flip() +
  
  labs(
    x = "Sample ID",
    y = expression(T[Delta47] ~ "(" * degree * "C)"),
    color = "Session"
  ) +
  
  theme_classic()


# Compare Δ47-derived temperatures between analytical sessions by strat level
#
# Purpose:
# Identify stratigraphic levels analyzed in both Session 1 and Session 2,
# calculate the mean temperature for each strat level within each session,
# and quantify the temperature difference between sessions.
#
# Positive dT values indicate Session 1 produced hotter temperatures than
# Session 2 at the same stratigraphic level.

paired <- IPL_D47_data %>%
  
  # Keep only analyses from Session 1 and Session 2 with valid temperatures
  # and valid stratigraphic positions
  filter(
    !is.na(Session),
    Session %in% c("Session 1", "Session 2"),
    !is.na(T47_preferred),
    !is.na(strat_height_m)
  ) %>%
  
  # Group replicate analyses by stratigraphic level and session
  group_by(strat_height_m, Session) %>%
  
  # Calculate mean temperature for each strat level in each session
  summarise(
    T_mean = mean(T47_preferred, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  
  # Reshape from long to wide format so each strat level occupies one row
  # with separate columns for Session 1 and Session 2 mean temperatures
  pivot_wider(
    names_from = Session,
    values_from = T_mean
  ) %>%
  
  # Retain only strat levels measured in both sessions
  filter(
    !is.na(`Session 1`),
    !is.na(`Session 2`)
  ) %>%
  
  # Calculate temperature offset between sessions
  mutate(
    dT = `Session 1` - `Session 2`
  )


# Inspect paired strat-level session comparison
#
# Interpretation:
#   dT > 0  -> Session 1 hotter
#   dT < 0  -> Session 2 hotter
print(paired)

summary(paired$dT)

paired %>%
  summarise(
    n_pairs = n(),
    mean_dT = mean(dT, na.rm = TRUE),
    median_dT = median(dT, na.rm = TRUE),
    sd_dT = sd(dT, na.rm = TRUE)
  )


# Paired t-test
#
# Tests whether the mean strat-level temperature difference between Session 1
# and Session 2 is significantly different from zero.
t.test(
  paired$`Session 1`,
  paired$`Session 2`,
  paired = TRUE
)


# 1:1 plot of Session 1 versus Session 2 temperatures
#
# Points above the dashed line indicate Session 1 is hotter.
# Points below the dashed line indicate Session 2 is hotter.
ggplot(paired,
       aes(x = `Session 2`,
           y = `Session 1`)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "red"
  ) +
  geom_point(size = 3) +
  coord_equal() +
  labs(
    x = expression("Session 2 " * T[Delta47] * " (" * degree * "C)"),
    y = expression("Session 1 " * T[Delta47] * " (" * degree * "C)")
  ) +
  theme_classic()


# Stratigraphic plot of session offsets
#
# Useful for identifying whether the session offset is systematic or driven by
# one/few stratigraphic levels.
ggplot(paired,
       aes(x = dT,
           y = strat_height_m)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "red"
  ) +
  geom_point(size = 3) +
  labs(
    x = expression(Delta * "T = Session 1 - Session 2 (" * degree * "C)"),
    y = "Stratigraphic Level (m)"
  ) +
  theme_classic()

# OPTIONAL: Restrict dataset to Session 2 only
#
# Session comparison indicates a small tendency toward warmer temperatures in
# Session 1, although the offset is not statistically significant. The full
# dataset is retained for all primary analyses.
#
# Uncomment the code below to perform sensitivity analyses using only
# Session 2 measurements.

# UN-COMMENT IF WANT TO REMOVE SESSION 1 DATA
# D47_primary_data <- D47_primary_data %>%
#   filter(Session == "Session 2")



# Inspect replicate-level Δ47 temperatures by sample
#
# Purpose: identify samples with unusually high within-sample
# scatter prior to summarization.

ggplot(IPL_D47_primary_data,
       aes(x = reorder(MLA_sample_id,
                       T47_preferred,
                       median),
           y = T47_preferred)) +
  geom_point() +
  coord_flip()


# Calculate within-sample temperature variability
sample_sd <- IPL_D47_primary_data %>%
  group_by(MLA_sample_id) %>%
  summarise(
    n = n(),
    mean_T = mean(T47_preferred),
    sd_T = sd(T47_preferred),
    range_T = max(T47_preferred) -
      min(T47_preferred)
  )

# Summarize distribution of within-sample variability
summary(sample_sd$sd_T)
summary(sample_sd$range_T)

# Identify samples with the largest temperature ranges
sample_sd %>%
  arrange(desc(range_T))


# Plot within-sample standard deviation for each sample
ggplot(sample_sd,
       aes(x = reorder(MLA_sample_id, sd_T),
           y = sd_T)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Sample ID",
    y = "Within-sample SD (°C)"
  ) +
  theme_classic()

# Plot individual Δ47 temperature replicates by sample, colored by session
#
# Purpose:
# Identify whether high within-sample scatter is driven by:
#   (1) a single outlying replicate,
#   (2) a systematic Session 1 vs Session 2 offset, or
#   (3) broadly scattered measurements across all replicates.

ggplot(IPL_D47_primary_data,
       aes(
         x = reorder(
           MLA_sample_id,
           T47_preferred,
           FUN = median,
           na.rm = TRUE
         ),
         y = T47_preferred,
         color = Session
       )) +
  
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8
  ) +
  
  coord_flip() +
  
  labs(
    x = "Sample ID",
    y = expression(T[Delta47] ~ "(" * degree * "C)"),
    color = "Session"
  ) +
  theme_classic()


# 5. Clean Snell et al. (2013) Δ47 data --------------------------------
str(Snell2013_summary)

# Separate Snell et al. (2013) samples by sample type
#
# Paleosol micrite = primary paleoclimate target
# Altered carbonate, fracture spar, and bivalves are retained but separated
# because they are not equivalent to primary paleosol micrite.

Snell2013_summary <- Snell2013_summary %>%
  mutate(
    Snell_sample_group = case_when(
      Snell_sample_type == "Paleosol Micrite" ~ "Primary paleosol micrite",
      Snell_sample_type == "Altered Paleosol carbonate" ~ "Altered paleosol carbonate",
      Snell_sample_type == "Fracture Spar" ~ "Fracture spar",
      Snell_sample_type == "Bivalve Fossil" ~ "Bivalve fossil",
      TRUE ~ "Other"
    )
  )

# Check grouping
table(Snell2013_summary$Snell_sample_group)

# Plot Snell et al. (2013) Δ47 temperatures through time
#
# Sample types are represented by point shape so primary paleosol micrites can
# be visually distinguished from altered carbonates, fracture spars, and
# bivalve fossils.

ggplot(Snell2013_summary,
       aes(x = Snell_mean_T47_C,
           y = Snell_Age_Ma,
           shape = Snell_sample_type)) +
  
  geom_point(size = 3) +
  
  scale_y_reverse() +
  
  labs(
    x = expression(T[47] ~ "(" * degree * "C)"),
    y = "Age (Ma)",
    shape = "Sample Type",
    title = expression("Snell et al. (2013) " * Delta[47] * " Temperatures")
  ) +
  
  guides(
    shape = guide_legend(nrow = 2)
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )


# Plot Snell et al. (2013) Δ47 temperature versus reconstructed water δ18O
#
# Purpose:
# Evaluate whether higher Δ47 temperatures are associated with anomalous
# reconstructed δ18Owater values, which could indicate alteration or
# recrystallization.

ggplot(Snell2013_summary,
       aes(x = Snell_mean_T47_C,
           y = Snell_mean_d18Ow_vsmow,
           shape = Snell_sample_type)) +
  
  geom_errorbarh(
    aes(
      xmin = Snell_mean_T47_C - Snell_se_T47_C,
      xmax = Snell_mean_T47_C + Snell_se_T47_C
    ),
    height = 0,
    alpha = 0.5
  ) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  labs(
    x = expression(T[47] ~ "(" * degree * "C)"),
    y = expression(delta^{18} * O[water] ~ "(‰ VSMOW)"),
    shape = "Sample Type"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "bottom"
  )

# Separate Snell et al. (2013) data by sample type
#
# Paleosol micrites are the primary paleoclimate archive and will eventually
# be merged with IPL and CU Δ47 datasets for comparison.
#
# Altered paleosol carbonates, fracture spars, and bivalves are retained as
# separate datasets for evaluating diagenetic and non-pedogenic signatures.

Snell2013_micrite <- Snell2013_summary %>%
  filter(Snell_sample_type == "Paleosol Micrite")

Snell2013_altered <- Snell2013_summary %>%
  filter(Snell_sample_type == "Altered Paleosol carbonate")

Snell2013_spar <- Snell2013_summary %>%
  filter(Snell_sample_type == "Fracture Spar")

Snell2013_bivalve <- Snell2013_summary %>%
  filter(Snell_sample_type == "Bivalve Fossil")

# Inspect sample counts
nrow(Snell2013_micrite)
nrow(Snell2013_altered)
nrow(Snell2013_spar)
nrow(Snell2013_bivalve)

# Verify separation
table(Snell2013_summary$Snell_sample_type)


# 6. Clean CU Boulder Δ47 data --------------------------------
str(CU_summary)


# Plot CU Δ47 temperatures versus stratigraphic position
#
# Purpose:
# Visualize CU temperature trends through the stratigraphic section and identify
# unusually warm or cool samples that may warrant further investigation.

ggplot(CU_summary,
       aes(x = CU_mean_T47_C,
           y = strat_height_m)) +
  
  geom_errorbarh(
    aes(
      xmin = CU_mean_T47_C - CU_2se_T47_C,
      xmax = CU_mean_T47_C + CU_2se_T47_C
    ),
    height = 0,
    alpha = 0.6
  ) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  labs(
    x = expression(T[47] ~ "(" * degree * "C)"),
    y = "Stratigraphic Height (m)"
  ) +
  
  theme_classic()


# Plot CU Δ47 temperature versus reconstructed soil-water δ18O
#
# Purpose:
# Evaluate whether high temperatures are associated with anomalous δ18Ow
# values that could indicate recrystallization, diagenesis, or alteration.

ggplot(CU_summary,
       aes(x = CU_mean_T47_C,
           y = CU_mean_d18Ow_vsmow)) +
  
  # Horizontal temperature uncertainty
  geom_errorbarh(
    aes(
      xmin = CU_mean_T47_C - CU_2se_T47_C,
      xmax = CU_mean_T47_C + CU_2se_T47_C
    ),
    height = 0,
    alpha = 0.5
  ) +
  
  # Vertical water-isotope uncertainty
  geom_errorbar(
    aes(
      ymin = CU_mean_d18Ow_vsmow - CU_se_d18Ow_vsmow,
      ymax = CU_mean_d18Ow_vsmow + CU_se_d18Ow_vsmow
    ),
    width = 0,
    alpha = 0.5
  ) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  labs(
    x = expression(T[47] ~ "(" * degree * "C)"),
    y = expression(delta^{18} * O[water] ~ "(‰ VSMOW)")
  ) +
  
  theme_classic()

# 7. Clean δ13C and δ18O carbonate data ------------------------


# Koch 1992/1995 data ~~~~~~~~~~~~~~

str(koch)

# Classify Koch samples as SPAR or primary micrite
koch <- koch %>%
  mutate(
    Koch_sample_type = ifelse(
      grepl("SPAR", MLA_sample_id),
      "SPAR",
      "MICRITE"
    )
  )

# Carbon vs. oxygen isotope crossplot
ggplot(koch,
       aes(x = Koch_d18Ocarb_vpdb,
           y = Koch_d13Ccarb_vpdb,
           shape = Koch_sample_type)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_reverse() +
  labs(
    x = expression(delta^18 * O[carb] ~ "(‰ VPDB)"),
    y = expression(delta^13 * C[carb] ~ "(‰ VPDB)"),
    shape = "Sample Type"
  ) +
  theme_classic()

# Separate SPAR and primary micrite samples
koch_spar <- koch %>%
  filter(Koch_sample_type == "SPAR")

koch_primary <- koch %>%
  filter(Koch_sample_type == "MICRITE")

# Quick stratigraphic checks for primary Koch δ18O values
ggplot(koch_primary,
       aes(x = Koch_d18Ocarb_vpdb,
           y = strat_height_m)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    x = expression(delta^18 * O[carb] ~ "(‰ VPDB)"),
    y = "Stratigraphic Height (m)"
  ) +
  theme_classic()

ggplot(koch_primary,
       aes(x = Koch_d18Ocarb_vsmow,
           y = strat_height_m)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    x = expression(delta^18 * O[carb] ~ "(‰ VSMOW)"),
    y = "Stratigraphic Height (m)"
  ) +
  theme_classic()

# Replicate δ18O values by horizon
ggplot(
  koch_primary %>% filter(!is.na(Koch_d18Ocarb_vsmow)),
  aes(
    x = reorder(
      MLA_horizon_id,
      Koch_d18Ocarb_vsmow,
      FUN = median,
      na.rm = TRUE
    ),
    y = Koch_d18Ocarb_vsmow
  )
) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.75) +
  coord_flip() +
  labs(
    x = "Horizon ID",
    y = expression(delta^18 * O[carb] ~ "(‰ VSMOW)")
  ) +
  theme_classic()

# Check sample counts
nrow(koch)
nrow(koch_spar)
nrow(koch_primary)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bowen et al. (2001) PB data

str(bowen)

# Quick sanity-check plots for Bowen carbonate isotope data

ggplot(bowen,
       aes(x = Bowen_d13Ccarb_vpdb,
           y = strat_height_m)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    x = expression(delta^13 * C[carb] ~ "(‰ VPDB)"),
    y = "Stratigraphic Height (m)"
  ) +
  theme_classic()

ggplot(bowen,
       aes(x = Bowen_d18Ocarb_vsmow,
           y = strat_height_m)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    x = expression(delta^18 * O[carb] ~ "(‰ VSMOW)"),
    y = "Stratigraphic Height (m)"
  ) +
  theme_classic()

ggplot(bowen,
       aes(x = Bowen_d18Ocarb_vpdb,
           y = strat_height_m)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    x = expression(delta^18 * O[carb] ~ "(‰ VPDB)"),
    y = "Stratigraphic Height (m)"
  ) +
  theme_classic()

ggplot(bowen,
       aes(x = Bowen_d18Ocarb_vsmow)) +
  geom_histogram(bins = 20) +
  labs(
    x = expression(delta^18 * O[carb] ~ "(‰ VSMOW)"),
    y = "Count"
  ) +
  theme_classic()

ggplot(bowen,
       aes(x = Bowen_d18Ocarb_vpdb)) +
  geom_histogram(bins = 20) +
  labs(
    x = expression(delta^18 * O[carb] ~ "(‰ VPDB)"),
    y = "Count"
  ) +
  theme_classic()


# Remove likely altered samples:
# 1) Exclude samples with δ18Ocarb < -11‰ VPDB
# 2) Exclude samples identified as fracture spar

bowen_primary <- bowen %>%
  filter(
    Bowen_d18Ocarb_vpdb > -11,
    !grepl("SPAR", MLA_sample_id)
  )

bowen_spar <- bowen %>%
  filter(
    Bowen_d18Ocarb_vpdb <= -11 |
      grepl("SPAR", MLA_sample_id)
  )


# Histogram of retained carbonate δ18O values
ggplot(bowen_primary,
       aes(x = Bowen_d18Ocarb_vpdb)) +
  geom_histogram(bins = 20) +
  labs(
    x = expression(delta^18 * O[carb] ~ "(‰ VPDB)"),
    y = "Count",
    title = expression("Retained Bowen " * delta^18 * O[carb] ~ "Values")
  ) +
  theme_classic()


# Stratigraphic distribution of retained carbonate δ18O values
ggplot(bowen_primary,
       aes(x = Bowen_d18Ocarb_vpdb,
           y = strat_height_m)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    x = expression(delta^18 * O[carb] ~ "(‰ VPDB)"),
    y = "Stratigraphic Height (m)",
    title = expression("Retained Bowen " * delta^18 * O[carb] ~ "Values")
  ) +
  theme_classic()

# Check sample counts
nrow(bowen)
nrow(bowen_primary)
nrow(bowen_spar)




# 8. Summarize each dataset by meter level ---------------

## 8.1 Summarize IPL Δ′17O data ---------------------------------
#
# Calculate horizon-level mean carbonate triple oxygen isotope values
# and associated uncertainties. Analyses are grouped by MLA_horizon_id,
# producing one row per stratigraphic horizon.
#
# Error treatment:
# Observed replicate standard deviations are calculated for each horizon.
# Because replicate scatter may underestimate true analytical uncertainty
# when sample size is small, a minimum Δ′17O standard deviation of
# 12 per meg is imposed. This value approximates long-term analytical
# reproducibility and is used whenever:
#   (1) only a single analysis exists for a horizon, or
#   (2) observed replicate SD is less than 12 per meg.
#
# Adjusted standard errors and 95% confidence intervals are calculated
# from this minimum uncertainty threshold.

generic_sd <- 12

IPL17O_summary <- IPL_D17O_data_final %>%
  group_by(MLA_horizon_id) %>%
  summarise(
    
    # δ′17O carbonate
    IPL17O_mean_dp17Ocarb = mean(dp17Ocarb_SMOWSLAP, na.rm = TRUE),
    IPL17O_sd_dp17Ocarb   = sd(dp17Ocarb_SMOWSLAP, na.rm = TRUE),
    IPL17O_se_dp17Ocarb   = IPL17O_sd_dp17Ocarb /
      sqrt(sum(!is.na(dp17Ocarb_SMOWSLAP))),
    IPL17O_n_dp17Ocarb    = sum(!is.na(dp17Ocarb_SMOWSLAP)),
    
    # δ′18O carbonate
    IPL17O_mean_dp18Ocarb = mean(dp18Ocarb_SMOWSLAP, na.rm = TRUE),
    IPL17O_sd_dp18Ocarb   = sd(dp18Ocarb_SMOWSLAP, na.rm = TRUE),
    IPL17O_se_dp18Ocarb   = IPL17O_sd_dp18Ocarb /
      sqrt(sum(!is.na(dp18Ocarb_SMOWSLAP))),
    IPL17O_n_dp18Ocarb    = sum(!is.na(dp18Ocarb_SMOWSLAP)),
    
    # Δ′17O carbonate
    IPL17O_mean_Dp17Ocarb = mean(
      Dp17Ocarb_permeg_final_correction,
      na.rm = TRUE
    ),
    IPL17O_sd_Dp17Ocarb = sd(
      Dp17Ocarb_permeg_final_correction,
      na.rm = TRUE
    ),
    IPL17O_se_Dp17Ocarb = IPL17O_sd_Dp17Ocarb /
      sqrt(sum(!is.na(Dp17Ocarb_permeg_final_correction))),
    IPL17O_n_Dp17Ocarb = sum(
      !is.na(Dp17Ocarb_permeg_final_correction)
    ),
    
    IPL17O_n_analyses = n(),
    
    .groups = "drop"
  ) %>%
  
  # Attach stratigraphic height after summarizing so each horizon receives
  # one stratigraphic position from the project sample list.
  left_join(
    IPL_sample_list %>%
      select(MLA_horizon_id, strat_height_m) %>%
      distinct(),
    by = "MLA_horizon_id"
  ) %>%
  
  mutate(
    
    # Impose minimum Δ′17O uncertainty of 12 per meg
    IPL17O_sd_Dp17Ocarb_adj = ifelse(
      IPL17O_n_Dp17Ocarb < 2 |
        IPL17O_sd_Dp17Ocarb < generic_sd,
      generic_sd,
      IPL17O_sd_Dp17Ocarb
    ),
    
    IPL17O_se_Dp17Ocarb_adj =
      IPL17O_sd_Dp17Ocarb_adj /
      sqrt(IPL17O_n_Dp17Ocarb),
    
    IPL17O_ci95_Dp17Ocarb_adj =
      IPL17O_se_Dp17Ocarb_adj * 1.96,
    
    # Convert mean δ′ values back to conventional δ notation
    IPL17O_mean_d17Ocarb =
      1000 * (exp(IPL17O_mean_dp17Ocarb / 1000) - 1),
    
    IPL17O_mean_d18Ocarb =
      1000 * (exp(IPL17O_mean_dp18Ocarb / 1000) - 1)
  )


## 8.2 Summarize IPL Δ47 data ----------------------------

# Summarize IPL D47 data from micrite/microspar
#
# Produces a sample-level Δ47 temperature summary table for downstream
# analyses. Primary samples exclude analyses identified as SPAR.
mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

min_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

max_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

summarize_IPLD47 <- function(data, sample_type) {
  data %>%
    group_by(MLA_sample_id, MLA_horizon_id) %>%
    summarise(
      strat_height_m = mean_or_na(strat_height_m),
      
      IPLD47_n_T47 = sum(!is.na(T47_preferred)),
      IPLD47_mean_T47_C = mean_or_na(T47_preferred),
      IPLD47_sd_T47_C = sd(T47_preferred, na.rm = TRUE),
      IPLD47_se_T47_C = IPLD47_sd_T47_C / sqrt(IPLD47_n_T47),
      
      IPLD47_min_T47_C = min_or_na(T47_preferred),
      IPLD47_max_T47_C = max_or_na(T47_preferred),
      IPLD47_range_T47_C =
        IPLD47_max_T47_C - IPLD47_min_T47_C,
      
      # Preserve corrected mineral isotope data
      IPL_NuDog_d13Ccarb_VPDB = mean_or_na(
        IPL_NuDog_d13Ccarb_VPDB
      ),
      
      IPL_NuDog_d18Ocarb_VPDB = mean_or_na(
        IPL_NuDog_d18Ocarb_VPDB
      ),
      
      IPL_NuDog_d18Ocarb_VSMOW = mean_or_na(
        IPL_NuDog_d18Ocarb_VSMOW
      ),
      
      IPLD47_sessions = paste(
        sort(unique(na.omit(Session))),
        collapse = ", "
      ),
      
      IPLD47_IPLnums = paste(
        sort(unique(na.omit(IPLnum))),
        collapse = ", "
      ),
      
      .groups = "drop"
    ) %>%
    mutate(
      IPLD47_se_T47_C = if_else(
        IPLD47_n_T47 > 0,
        IPLD47_se_T47_C,
        NA_real_
      ),
      IPLD47_sample_type = sample_type
    ) %>%
    arrange(strat_height_m)
}

# Primary micrite/microspar summary
IPLD47_primary_summary <- summarize_IPLD47(
  IPL_D47_primary_data,
  sample_type = "primary"
)

# Spar summary
IPLD47_spar_summary <- summarize_IPLD47(
  IPL_D47_SPAR_data,
  sample_type = "spar"
)

IPLD47_primary_summary
IPLD47_spar_summary

# Export summary tables
write_csv(
  IPLD47_primary_summary,
  here("data", "processed", "IPLD47_primary_summary.csv")
)

write_csv(
  IPLD47_spar_summary,
  here("data", "processed", "IPLD47_spar_summary.csv")
)

## 8.3 Summarize Koch data ----------------------

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Summarize Koch primary carbonate data

koch_summary <- koch_primary %>%
  group_by(MLA_horizon_id, strat_height_m) %>%
  summarise(
    Koch_mean_d13Ccarb_vpdb = mean(Koch_d13Ccarb_vpdb, na.rm = TRUE),
    Koch_se_d13Ccarb_vpdb   = sd(Koch_d13Ccarb_vpdb, na.rm = TRUE) /
      sqrt(sum(!is.na(Koch_d13Ccarb_vpdb))),
    Koch_n_d13Ccarb         = sum(!is.na(Koch_d13Ccarb_vpdb)),
    
    Koch_mean_d18Ocarb_vpdb = mean(Koch_d18Ocarb_vpdb, na.rm = TRUE),
    Koch_se_d18Ocarb_vpdb   = sd(Koch_d18Ocarb_vpdb, na.rm = TRUE) /
      sqrt(sum(!is.na(Koch_d18Ocarb_vpdb))),
    Koch_n_d18Ocarb_vpdb    = sum(!is.na(Koch_d18Ocarb_vpdb)),
    
    Koch_mean_d18Ocarb_vsmow = mean(Koch_d18Ocarb_vsmow, na.rm = TRUE),
    Koch_se_d18Ocarb_vsmow   = sd(Koch_d18Ocarb_vsmow, na.rm = TRUE) /
      sqrt(sum(!is.na(Koch_d18Ocarb_vsmow))),
    Koch_n_d18Ocarb_vsmow    = sum(!is.na(Koch_d18Ocarb_vsmow)),
    
    .groups = "drop"
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Summarize Koch fracture spar data

koch_spar_summary <- koch_spar %>%
  group_by(MLA_horizon_id, strat_height_m) %>%
  summarise(
    KochSPAR_mean_d13Ccarb_vpdb = mean(Koch_d13Ccarb_vpdb, na.rm = TRUE),
    KochSPAR_se_d13Ccarb_vpdb   = sd(Koch_d13Ccarb_vpdb, na.rm = TRUE) /
      sqrt(sum(!is.na(Koch_d13Ccarb_vpdb))),
    KochSPAR_n_d13Ccarb         = sum(!is.na(Koch_d13Ccarb_vpdb)),
    
    KochSPAR_mean_d18Ocarb_vpdb = mean(Koch_d18Ocarb_vpdb, na.rm = TRUE),
    KochSPAR_se_d18Ocarb_vpdb   = sd(Koch_d18Ocarb_vpdb, na.rm = TRUE) /
      sqrt(sum(!is.na(Koch_d18Ocarb_vpdb))),
    KochSPAR_n_d18Ocarb_vpdb    = sum(!is.na(Koch_d18Ocarb_vpdb)),
    
    KochSPAR_mean_d18Ocarb_vsmow = mean(Koch_d18Ocarb_vsmow, na.rm = TRUE),
    KochSPAR_se_d18Ocarb_vsmow   = sd(Koch_d18Ocarb_vsmow, na.rm = TRUE) /
      sqrt(sum(!is.na(Koch_d18Ocarb_vsmow))),
    KochSPAR_n_d18Ocarb_vsmow    = sum(!is.na(Koch_d18Ocarb_vsmow)),
    
    .groups = "drop"
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Quick checks

ggplot(koch_summary,
       aes(x = Koch_mean_d18Ocarb_vsmow,
           y = strat_height_m)) +
  geom_errorbarh(
    aes(
      xmin = Koch_mean_d18Ocarb_vsmow - Koch_se_d18Ocarb_vsmow,
      xmax = Koch_mean_d18Ocarb_vsmow + Koch_se_d18Ocarb_vsmow
    ),
    height = 0
  ) +
  geom_point(size = 2) +
  labs(
    x = expression(delta^18 * O[carb] ~ "(‰ VSMOW)"),
    y = "Stratigraphic Height (m)"
  ) +
  theme_classic()

ggplot(koch_summary,
       aes(x = Koch_mean_d13Ccarb_vpdb,
           y = strat_height_m)) +
  geom_errorbarh(
    aes(
      xmin = Koch_mean_d13Ccarb_vpdb - Koch_se_d13Ccarb_vpdb,
      xmax = Koch_mean_d13Ccarb_vpdb + Koch_se_d13Ccarb_vpdb
    ),
    height = 0
  ) +
  geom_point(size = 2) +
  labs(
    x = expression(delta^13 * C[carb] ~ "(‰ VPDB)"),
    y = "Stratigraphic Height (m)"
  ) +
  theme_classic()



## 8.4 Summarize Bowen data ----------------------

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Summarize Bowen primary carbonate data

bowen_summary <- bowen_primary %>%
  group_by(MLA_horizon_id, strat_height_m) %>%
  summarise(
    Bowen_mean_d13Ccarb_vpdb = mean(Bowen_d13Ccarb_vpdb, na.rm = TRUE),
    Bowen_se_d13Ccarb_vpdb   = sd(Bowen_d13Ccarb_vpdb, na.rm = TRUE) /
      sqrt(sum(!is.na(Bowen_d13Ccarb_vpdb))),
    Bowen_n_d13Ccarb         = sum(!is.na(Bowen_d13Ccarb_vpdb)),
    
    Bowen_mean_d18Ocarb_vpdb = mean(Bowen_d18Ocarb_vpdb, na.rm = TRUE),
    Bowen_se_d18Ocarb_vpdb   = sd(Bowen_d18Ocarb_vpdb, na.rm = TRUE) /
      sqrt(sum(!is.na(Bowen_d18Ocarb_vpdb))),
    Bowen_n_d18Ocarb_vpdb    = sum(!is.na(Bowen_d18Ocarb_vpdb)),
    
    Bowen_mean_d18Ocarb_vsmow = mean(Bowen_d18Ocarb_vsmow, na.rm = TRUE),
    Bowen_se_d18Ocarb_vsmow   = sd(Bowen_d18Ocarb_vsmow, na.rm = TRUE) /
      sqrt(sum(!is.na(Bowen_d18Ocarb_vsmow))),
    Bowen_n_d18Ocarb_vsmow    = sum(!is.na(Bowen_d18Ocarb_vsmow)),
    
    .groups = "drop"
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Summarize Bowen excluded / altered carbonate data

#bowen_excluded_summary <- bowen_excluded %>%
#  group_by(MLA_horizon_id, strat_height_m) %>%
#  summarise(
#    BowenExcluded_mean_d13Ccarb_vpdb = mean(Bowen_d13Ccarb_vpdb, na.rm = TRUE),
#    BowenExcluded_se_d13Ccarb_vpdb   = sd(Bowen_d13Ccarb_vpdb, na.rm = TRUE) /
#      sqrt(sum(!is.na(Bowen_d13Ccarb_vpdb))),
#    BowenExcluded_n_d13Ccarb         = sum(!is.na(Bowen_d13Ccarb_vpdb)),
    
#    BowenExcluded_mean_d18Ocarb_vpdb = mean(Bowen_d18Ocarb_vpdb, na.rm = TRUE),
#    BowenExcluded_se_d18Ocarb_vpdb   = sd(Bowen_d18Ocarb_vpdb, na.rm = TRUE) /
#      sqrt(sum(!is.na(Bowen_d18Ocarb_vpdb))),
#    BowenExcluded_n_d18Ocarb_vpdb    = sum(!is.na(Bowen_d18Ocarb_vpdb)),
    
#    BowenExcluded_mean_d18Ocarb_vsmow = mean(Bowen_d18Ocarb_vsmow, na.rm = TRUE),
#    BowenExcluded_se_d18Ocarb_vsmow   = sd(Bowen_d18Ocarb_vsmow, na.rm = TRUE) /
#     sqrt(sum(!is.na(Bowen_d18Ocarb_vsmow))),
#    BowenExcluded_n_d18Ocarb_vsmow    = sum(!is.na(Bowen_d18Ocarb_vsmow)),
    
#   .groups = "drop"
# )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Quick checks

ggplot(bowen_summary,
       aes(x = Bowen_mean_d18Ocarb_vsmow,
           y = strat_height_m)) +
  geom_errorbarh(
    aes(
      xmin = Bowen_mean_d18Ocarb_vsmow - Bowen_se_d18Ocarb_vsmow,
      xmax = Bowen_mean_d18Ocarb_vsmow + Bowen_se_d18Ocarb_vsmow
    ),
    height = 0
  ) +
  geom_point(size = 2) +
  labs(
    x = expression(delta^18 * O[carb] ~ "(‰ VSMOW)"),
    y = "Stratigraphic Height (m)"
  ) +
  theme_classic()

ggplot(bowen_summary,
       aes(x = Bowen_mean_d13Ccarb_vpdb,
           y = strat_height_m)) +
  geom_errorbarh(
    aes(
      xmin = Bowen_mean_d13Ccarb_vpdb - Bowen_se_d13Ccarb_vpdb,
      xmax = Bowen_mean_d13Ccarb_vpdb + Bowen_se_d13Ccarb_vpdb
    ),
    height = 0
  ) +
  geom_point(size = 2) +
  labs(
    x = expression(delta^13 * C[carb] ~ "(‰ VPDB)"),
    y = "Stratigraphic Height (m)"
  ) +
  theme_classic()

ggplot(bowen_summary,
       aes(x = Bowen_mean_d18Ocarb_vpdb,
           y = Bowen_mean_d13Ccarb_vpdb)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_x_reverse() +
  labs(
    x = expression(delta^18 * O[carb] ~ "(‰ VPDB)"),
    y = expression(delta^13 * C[carb] ~ "(‰ VPDB)")
  ) +
  theme_classic()

## 8.5 Summarize spar/microspar data --------

# ---- Combine IPL spar + Snell altered/spar data ----

spar_altered_combined <- bind_rows(
  
  IPL_D47_SPAR_data %>%
    transmute(
      source = "IPL D47",
      sample_group = "SPAR",
      sample_type = "SPAR",
      MLA_sample_id,
      MLA_horizon_id,
      strat_height_m,
      d13Ccarb_vpdb = IPL_NuDog_d13Ccarb_VPDB,
      d18Ocarb_vpdb = IPL_NuDog_d18Ocarb_VPDB,
      d18Ocarb_vsmow = IPL_NuDog_d18Ocarb_VSMOW,
      D47 = D47.CDES,
      D47_carb_corr = D47.CDES.Carb.Corr,
      T47_C = T.D47..Petersen,
      T47_preferred,
      T47_se_C = NA_real_
    ),
  
  Snell2013_altered %>%
    transmute(
      source = "Snell et al. 2013",
      sample_group = Snell_sample_group,
      sample_type = Snell_sample_type,
      MLA_sample_id,
      MLA_horizon_id,
      strat_height_m,
      d13Ccarb_vpdb = Snell_mean_d13Ccarb_vpdb,
      d18Ocarb_vpdb = NA_real_,
      d18Ocarb_vsmow = Snell_mean_d18Ocarb_vsmow,
      D47 = Snell_mean_D47,
      D47_carb_corr = NA_real_,
      T47_C = Snell_mean_T47_C,
      T47_preferred = Snell_mean_T47_C,
      T47_se_C = Snell_se_T47_C
    ),
  
  Snell2013_spar %>%
    transmute(
      source = "Snell et al. 2013",
      sample_group = Snell_sample_group,
      sample_type = Snell_sample_type,
      MLA_sample_id,
      MLA_horizon_id,
      strat_height_m,
      d13Ccarb_vpdb = Snell_mean_d13Ccarb_vpdb,
      d18Ocarb_vpdb = NA_real_,
      d18Ocarb_vsmow = Snell_mean_d18Ocarb_vsmow,
      D47 = Snell_mean_D47,
      D47_carb_corr = NA_real_,
      T47_C = Snell_mean_T47_C,
      T47_preferred = Snell_mean_T47_C,
      T47_se_C = Snell_se_T47_C
    )
) %>%
  arrange(strat_height_m, MLA_horizon_id, source, MLA_sample_id)

# 9. Combine summaries by horizon / strat ----------------------


# Check whether each summary has one row per MLA_horizon_id
check_horizon_dupes <- function(df) {
  df %>%
    count(MLA_horizon_id) %>%
    filter(n > 1)
}

check_horizon_dupes(IPL17O_summary)
check_horizon_dupes(IPLD47_primary_summary)
check_horizon_dupes(CU_summary)
check_horizon_dupes(Snell2013_micrite)
check_horizon_dupes(koch_summary)
check_horizon_dupes(bowen_summary)

summary_list <- list(
  IPL17O = IPL17O_summary,
  IPLD47 = IPLD47_primary_summary %>%
    select(-MLA_sample_id),
  CU = CU_summary %>%
    select(-MLA_sample_id),
  Snell = Snell2013_micrite %>%
    select(-MLA_sample_id),
  Koch = koch_summary,
  Bowen = bowen_summary
)

lapply(summary_list, function(x) nrow(x))


BHB_multiproxy_summary <- summary_list %>%
  purrr::reduce(
    full_join,
    by = c("MLA_horizon_id", "strat_height_m")
  ) %>%
  arrange(strat_height_m)

check_horizon_dupes(BHB_multiproxy_summary)

BHB_multiproxy_summary %>%
 filter(MLA_horizon_id == "PK95-SC-295") %>%
 select(MLA_horizon_id, strat_height_m)

# Summarize spar / altered data by horizon 

spar_altered_horizon_summary <- spar_altered_combined %>%
  group_by(MLA_horizon_id) %>%
  summarise(
    strat_height_m = mean(strat_height_m, na.rm = TRUE),
    
    sources = paste(sort(unique(source)), collapse = "; "),
    sample_groups = paste(sort(unique(sample_group)), collapse = "; "),
    sample_types = paste(sort(unique(sample_type)), collapse = "; "),
    
    n_analyses = n(),
    n_T47 = sum(!is.na(T47_C)),
    n_d18Ocarb = sum(!is.na(d18Ocarb_vsmow)),
    
    mean_T47_C = mean(T47_C, na.rm = TRUE),
    sd_T47_C = sd(T47_C, na.rm = TRUE),
    se_T47_C = sd_T47_C / sqrt(n_T47),
    min_T47_C = min(T47_C, na.rm = TRUE),
    max_T47_C = max(T47_C, na.rm = TRUE),
    
    mean_d18Ocarb_vsmow = mean(d18Ocarb_vsmow, na.rm = TRUE),
    sd_d18Ocarb_vsmow = sd(d18Ocarb_vsmow, na.rm = TRUE),
    se_d18Ocarb_vsmow = sd_d18Ocarb_vsmow / sqrt(n_d18Ocarb),
    min_d18Ocarb_vsmow = min(d18Ocarb_vsmow, na.rm = TRUE),
    max_d18Ocarb_vsmow = max(d18Ocarb_vsmow, na.rm = TRUE),
    
    mean_d13Ccarb_vpdb = mean(d13Ccarb_vpdb, na.rm = TRUE),
    sd_d13Ccarb_vpdb = sd(d13Ccarb_vpdb, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    across(
      where(is.numeric),
      ~ if_else(is.nan(.x) | is.infinite(.x), NA_real_, .x)
    )
  ) %>%
  arrange(strat_height_m)

# 10. Age Model --------------

# Assign ages to all horizons using the composite Bighorn Basin age model.
#
# matthews_age_model contains age-depth tie points derived from a combination
# of:
#   (1) the astronomically calibrated Polecat Bench core chronology of
#       Westerhold et al. (2018), projected onto Polecat Bench outcrop
#       stratigraphic heights using the onset of the PETM carbon isotope
#       excursion (CIE) as a shared datum, and
#   (2) the paleomagnetic age model of Secord et al. (2006), which provides
#       age control outside the interval directly constrained by Westerhold
#       et al. (2018).
#
# Ages for individual horizons are obtained by linear interpolation between
# adjacent age-depth tie points. Horizons falling outside the modeled range
# are assigned NA (rule = 1).

age_priors <- read.csv(
  here("data", "raw", "PCB-CFB_age_priors.csv")
)

str(age_priors)

# ---- Continuous composite age model ------

age_priors <- age_priors %>%
  mutate(
    est_Depth_m_PCB_outcrop = as.numeric(est_Depth_m_PCB_outcrop),
    Age_Ma_best_estimate = as.numeric(Age_Ma_best_estimate)
  ) %>%
  filter(
    !is.na(est_Depth_m_PCB_outcrop),
    !is.na(Age_Ma_best_estimate)
  ) %>%
  arrange(est_Depth_m_PCB_outcrop)

BHB_multiproxy_summary <- BHB_multiproxy_summary %>%
  mutate(
    Age_Ma = approx(
      x = age_priors$est_Depth_m_PCB_outcrop,
      y = age_priors$Age_Ma_best_estimate,
      xout = strat_height_m,
      rule = 1
    )$y
  )

# Quick-look age model diagnostics

# Age-depth model
ggplot() +
  geom_line(
    data = tibble(
      strat_height_m = seq(
        min(age_priors$est_Depth_m_PCB_outcrop),
        max(age_priors$est_Depth_m_PCB_outcrop),
        by = 1
      )
    ) %>%
      mutate(
        Age_Ma = approx(
          x = age_priors$est_Depth_m_PCB_outcrop,
          y = age_priors$Age_Ma_best_estimate,
          xout = strat_height_m,
          rule = 1
        )$y
      ),
    aes(x = Age_Ma, y = strat_height_m),
    linewidth = 1
  ) +
  geom_point(
    data = age_priors,
    aes(x = Age_Ma_best_estimate, y = est_Depth_m_PCB_outcrop),
    size = 2
  ) +
  scale_x_reverse() +
  labs(
    title = "Composite age-depth model",
    subtitle = "Piecewise linear interpolation through stratigraphic age priors",
    x = "Age (Ma)",
    y = "Stratigraphic height above K-Pg (m)"
  ) +
  theme_classic()


# ---- Apply composite age model to BHB multiproxy summary 

age_priors_clean <- age_priors %>%
  mutate(
    est_Depth_m_PCB_outcrop = as.numeric(est_Depth_m_PCB_outcrop),
    Age_Ma_best_estimate = as.numeric(Age_Ma_best_estimate)
  ) %>%
  filter(
    !is.na(est_Depth_m_PCB_outcrop),
    !is.na(Age_Ma_best_estimate)
  ) %>%
  arrange(est_Depth_m_PCB_outcrop)

BHB_multiproxy_summary <- BHB_multiproxy_summary %>%
  mutate(
    Age_Ma = approx(
      x = age_priors_clean$est_Depth_m_PCB_outcrop,
      y = age_priors_clean$Age_Ma_best_estimate,
      xout = strat_height_m,
      rule = 1
    )$y
  )

# Quick check
BHB_multiproxy_summary %>%
  select(MLA_horizon_id, strat_height_m, Age_Ma) %>%
  arrange(strat_height_m)

ggplot(
  BHB_multiproxy_summary %>%
    filter(
      !is.na(Age_Ma),
      !is.na(Koch_mean_d13Ccarb_vpdb)
    ),
  aes(
    x = Age_Ma,
    y = Koch_mean_d13Ccarb_vpdb
  )
) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(
    method = "loess",
    span = 0.2,
    se = FALSE,
    linewidth = 1.2,
    color = "black"
  ) +
  scale_x_reverse() +
  labs(
    x = "Age (Ma)",
    y = "Koch d13Ccarb (per mil VPDB)",
    title = "Koch paleosol carbonate d13C in age space"
  ) +
  theme_classic()

# ---- Add geologic stage from age model ----
BHB_multiproxy_summary <- BHB_multiproxy_summary %>%
  mutate(
    stage = case_when(
      Age_Ma >= 66.0 & Age_Ma < 61.6 ~ "Danian",
      Age_Ma >= 61.6 & Age_Ma < 59.2 ~ "Selandian",
      Age_Ma >= 59.2 & Age_Ma < 56.0 ~ "Thanetian",
      Age_Ma >= 56.0 & Age_Ma < 47.8 ~ "Ypresian",
      TRUE ~ NA_character_
    ),
    stage = factor(
      stage,
      levels = c("Danian", "Selandian", "Thanetian", "Ypresian")
    )
  )
# 11. Export summarized primary datasets  ---------------------

write_csv(
  IPL17O_summary,
  here("data", "processed", "IPL17O_summary.csv")
)

write_csv(
  IPLD47_primary_summary,
  here("data", "processed", "IPLD47_primary_summary.csv")
)

write_csv(
  CU_summary,
  here("data", "processed", "CU_summary.csv")
)

write_csv(
  Snell2013_micrite,
  here("data", "processed", "Snell2013_micrite_summary.csv")
)

write_csv(
  koch_summary,
  here("data", "processed", "koch_summary.csv")
)

write_csv(
  bowen_summary,
  here("data", "processed", "bowen_summary.csv")
)

write_csv(
  BHB_multiproxy_summary,
  here("data", "processed", "BHB_multiproxy_summary.csv")
)

# Export spar / altered carbonate datasets

write_csv(
  spar_altered_combined,
  here(
    "data",
    "processed",
    "spar_altered_combined.csv"
  )
)

write_csv(
  spar_altered_horizon_summary,
  here(
    "data",
    "processed",
    "spar_altered_horizon_summary.csv"
  )
)

# 12. Quick checks ---------------------------------------------

glimpse(BHB_multiproxy_summary)

BHB_multiproxy_summary %>%
  summarise(
    n_horizons = n(),
    min_strat = min(strat_height_m, na.rm = TRUE),
    max_strat = max(strat_height_m, na.rm = TRUE)
  )

# quick check on d13C 
paired_d13c <- BHB_multiproxy_summary %>%
  filter(
    !is.na(Koch_mean_d13Ccarb_vpdb),
    !is.na(IPL_NuDog_d13Ccarb_VPDB)
  )

ggplot(
  paired_d13c,
  aes(
    x = Koch_mean_d13Ccarb_vpdb,
    y = IPL_NuDog_d13Ccarb_VPDB
  )
) +
  geom_point(size = 3, alpha = 0.8) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    color = "gray50"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black"
  ) +
  labs(
    x = "Koch d13Ccarb (‰ VPDB)",
    y = "IPL NuDog d13Ccarb (‰ VPDB)",
    title = "Koch vs. IPL NuDog carbonate d13C"
  ) +
  theme_classic()

# -------------------------------------------------------------------
# Interactive d13Ccarb stratigraphic plot by source
# -------------------------------------------------------------------

library(tidyverse)
library(plotly)

d13C_by_source <- BHB_multiproxy_summary %>%
  
  select(
    MLA_horizon_id,
    strat_height_m,
    
    Koch_mean_d13Ccarb_vpdb,
    Koch_se_d13Ccarb_vpdb,
    
    Bowen_mean_d13Ccarb_vpdb,
    Bowen_se_d13Ccarb_vpdb,
    
    CU_mean_d13Ccarb_vpdb,
    Snell_mean_d13Ccarb_vpdb,
    IPL_NuDog_d13Ccarb_VPDB
  ) %>%
  
  pivot_longer(
    cols = c(
      Koch_mean_d13Ccarb_vpdb,
      Bowen_mean_d13Ccarb_vpdb,
      CU_mean_d13Ccarb_vpdb,
      Snell_mean_d13Ccarb_vpdb,
      IPL_NuDog_d13Ccarb_VPDB
    ),
    names_to = "source_raw",
    values_to = "d13Ccarb_vpdb"
  ) %>%
  
  mutate(
    source = recode(
      source_raw,
      Koch_mean_d13Ccarb_vpdb  = "Koch",
      Bowen_mean_d13Ccarb_vpdb = "Bowen et al. (2001)",
      CU_mean_d13Ccarb_vpdb    = "CU Boulder",
      Snell_mean_d13Ccarb_vpdb = "Snell et al. (2013)",
      IPL_NuDog_d13Ccarb_VPDB  = "IPL NuDog"
    ),
    
    # Only Koch and Bowen currently have d13C SE columns
    d13C_se_vpdb = case_when(
      source_raw == "Koch_mean_d13Ccarb_vpdb" ~
        Koch_se_d13Ccarb_vpdb,
      
      source_raw == "Bowen_mean_d13Ccarb_vpdb" ~
        Bowen_se_d13Ccarb_vpdb,
      
      TRUE ~ NA_real_
    ),
    
    d13C_lower = d13Ccarb_vpdb - d13C_se_vpdb,
    d13C_upper = d13Ccarb_vpdb + d13C_se_vpdb,
    
    source = factor(
      source,
      levels = c(
        "Koch",
        "Bowen et al. (2001)",
        "CU Boulder",
        "Snell et al. (2013)",
        "IPL NuDog"
      )
    ),
    
    hover_text = paste0(
      "Source: ", source,
      "<br>Horizon: ", MLA_horizon_id,
      "<br>Stratigraphic height: ",
      round(strat_height_m, 1), " m",
      "<br>δ13Ccarb: ",
      round(d13Ccarb_vpdb, 2), "‰ VPDB",
      ifelse(
        is.na(d13C_se_vpdb),
        "",
        paste0(
          "<br>SE: ±",
          round(d13C_se_vpdb, 2),
          "‰"
        )
      )
    )
  ) %>%
  
  filter(
    !is.na(d13Ccarb_vpdb),
    !is.na(strat_height_m)
  )

# Optional source-coverage summary
d13C_source_summary <- d13C_by_source %>%
  group_by(source) %>%
  summarise(
    n = n(),
    n_horizons = n_distinct(MLA_horizon_id),
    minimum_height_m = min(strat_height_m),
    maximum_height_m = max(strat_height_m),
    mean_d13Ccarb_vpdb = mean(d13Ccarb_vpdb),
    sd_d13Ccarb_vpdb = sd(d13Ccarb_vpdb),
    .groups = "drop"
  )

print(d13C_source_summary)

# -------------------------------------------------------------------
# Static ggplot used as the basis for the interactive figure
# -------------------------------------------------------------------

p_d13C_by_source <- ggplot(
  d13C_by_source,
  aes(
    x = d13Ccarb_vpdb,
    y = strat_height_m,
    colour = source,
    shape = source,
    text = hover_text
  )
) +
  
  # Horizontal SE bars where uncertainty is available
  geom_errorbarh(
    aes(
      xmin = d13C_lower,
      xmax = d13C_upper
    ),
    height = 0,
    linewidth = 0.35,
    alpha = 0.45,
    na.rm = TRUE
  ) +
  
  geom_point(
    size = 2.2,
    alpha = 0.80
  ) +
  
  # PETM interval; remove this layer if not wanted
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 1500,
    ymax = 1540,
    fill = "grey70",
    alpha = 0.20
  ) +
  
  scale_x_continuous(
    limits = c(-17, -4),
    breaks = seq(-16, -4, by = 2),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  
  scale_y_continuous(
    limits = c(0, 2300),
    breaks = seq(0, 2300, by = 200),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    x = "δ13Ccarb (‰ VPDB)",
    y = "Stratigraphic height (m)",
    colour = "Source",
    shape = "Source"
  ) +
  
  theme_classic(
    base_family = "Arial",
    base_size = 11
  ) +
  
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.4
    )
  )

# -------------------------------------------------------------------
# Convert to interactive Plotly object
# -------------------------------------------------------------------

p_d13C_by_source_interactive <- ggplotly(
  p_d13C_by_source,
  tooltip = "text",
  width = 850,
  height = 900
) %>%
  
  layout(
    hovermode = "closest",
    legend = list(
      title = list(text = "<b>Source</b>"),
      itemclick = "toggle",
      itemdoubleclick = "toggleothers"
    ),
    margin = list(
      l = 80,
      r = 180,
      b = 70,
      t = 25
    )
  ) %>%
  
  config(
    displaylogo = FALSE,
    responsive = TRUE,
    modeBarButtonsToRemove = c(
      "lasso2d",
      "select2d"
    )
  )

p_d13C_by_source_interactive

# 13. Data Availability ------

# ---- UM sample coverage jitter plots ----

UM_sample_coverage <- IPL_sample_list %>%
  filter(
    Status == "UM",
    !is.na(strat_height_m)
  ) %>%
  mutate(sample_set = "UM samples")

p_UM_coverage_full <- ggplot(
  UM_sample_coverage,
  aes(x = sample_set, y = strat_height_m)
) +
  geom_jitter(
    width = 0.0,
    height = 0,
    size = 2.5,
    alpha = 0.75
  ) +
  scale_y_continuous(
    breaks = seq(500, 2300, by = 100)
  ) +
  labs(
    title = "UM Sample Coverage",
    x = NULL,
    y = "Stratigraphic height (m)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

p_UM_coverage_zoom <- p_UM_coverage_full +
  coord_cartesian(ylim = c(1200, 1800)) +
  scale_y_continuous(
    breaks = seq(1200, 1800, by = 50)
  ) +
  labs(
    title = "UM Sample Coverage: P–E Boundary Interval"
  )

p_UM_coverage_full
p_UM_coverage_zoom





# IPL / CU 

# ---- 1. Horizon-level availability summary 

data_availability <- BHB_multiproxy_summary %>%
  mutate(
    has_IPL17O = !is.na(IPL17O_mean_Dp17Ocarb),
    has_IPLD47 = !is.na(IPLD47_mean_T47_C),
    has_CUD47  = !is.na(CU_mean_T47_C),
    has_any_D47 = has_IPLD47 | has_CUD47,
    
    data_combo = case_when(
      has_IPL17O & has_IPLD47 & has_CUD47 ~ "IPL 17O + IPL D47 + CU D47",
      has_IPL17O & has_IPLD47             ~ "IPL 17O + IPL D47",
      has_IPL17O & has_CUD47              ~ "IPL 17O + CU D47",
      has_IPLD47 & has_CUD47              ~ "IPL D47 + CU D47",
      has_IPL17O                          ~ "IPL 17O only",
      has_IPLD47                          ~ "IPL D47 only",
      has_CUD47                           ~ "CU D47 only",
      TRUE                                ~ "No IPL/CU isotope data"
    )
  ) %>%
  select(
    MLA_horizon_id,
    strat_height_m,
    data_combo,
    
    has_IPL17O,
    IPL17O_n_analyses,
    IPL17O_mean_Dp17Ocarb,
    IPL17O_sd_Dp17Ocarb,
    IPL17O_se_Dp17Ocarb,
    IPL17O_mean_d18Ocarb,
    
    has_IPLD47,
    IPLD47_n_T47,
    IPLD47_mean_T47_C,
    IPLD47_sd_T47_C,
    IPLD47_se_T47_C,
    IPLD47_range_T47_C,
    IPLD47_sessions,
    
    has_CUD47,
    CU_sample_id,
    CU_mean_T47_C,
    CU_2se_T47_C,
    
    has_any_D47
  ) %>%
  arrange(strat_height_m)

# View table
print(data_availability, n = 200)

# Save table
write_csv(
  data_availability,
  here::here("data", "processed", "IPL_CU_data_availability_by_horizon.csv")
)


## 17O but no D47, and D47 but no 17O ----

IPL17O_no_D47 <- data_availability %>%
  filter(has_IPL17O, !has_any_D47) %>%
  arrange(strat_height_m)

D47_no_IPL17O <- data_availability %>%
  filter(has_any_D47, !has_IPL17O) %>%
  arrange(strat_height_m)

IPL17O_and_D47 <- data_availability %>%
  filter(has_IPL17O, has_any_D47) %>%
  arrange(strat_height_m)

print(IPL17O_no_D47, n = 100)
print(D47_no_IPL17O, n = 100)
print(IPL17O_and_D47, n = 100)

write_csv(IPL17O_no_D47, here::here("data", "processed", "IPL17O_no_D47.csv"))
write_csv(D47_no_IPL17O, here::here("data", "processed", "D47_no_IPL17O.csv"))
write_csv(IPL17O_and_D47, here::here("data", "processed", "IPL17O_and_D47.csv"))

# ---- UM samples still missing IPL D47 ----

UM_missing_IPLD47 <- IPL_sample_list %>%
  filter(
    Status == "UM",
    !is.na(MLA_horizon_id),
    !is.na(strat_height_m)
  ) %>%
  distinct(
    MLA_horizon_id,
    strat_height_m,
    Status,
    D47_IPL,
    X17O_IPL,
    Reps_needed,
    mg_needed,
    MLA_notes
  ) %>%
  left_join(
    IPLD47_primary_summary %>%
      distinct(
        MLA_horizon_id,
        has_IPLD47 = IPLD47_mean_T47_C,
        IPLD47_n_T47,
        IPLD47_mean_T47_C,
        IPLD47_se_T47_C
      ),
    by = "MLA_horizon_id"
  ) %>%
  mutate(
    has_IPLD47 = !is.na(has_IPLD47)
  ) %>%
  filter(!has_IPLD47) %>%
  arrange(strat_height_m)

print.data.frame(UM_missing_IPLD47)

UM_missing_IPLD47 %>%
  filter(!has_IPLD47) %>%
  arrange(strat_height_m) %>%
  select(MLA_horizon_id, strat_height_m)

# quick count
UM_missing_IPLD47 %>%
  summarise(
    n_missing_horizons = n(),
    min_strat = min(strat_height_m, na.rm = TRUE),
    max_strat = max(strat_height_m, na.rm = TRUE)
  )
## Quick counts ----

availability_counts <- data_availability %>%
  count(data_combo, sort = TRUE)

print(availability_counts)

write_csv(
  availability_counts,
  here::here("data", "processed", "IPL_CU_availability_counts.csv")
)

library(plotly)

# ---- Interactive jitter plot: D47 availability ----
library(plotly)
library(ggplot2)

p_D47 <- data_availability %>%
  mutate(
    D47_status = if_else(has_any_D47, "Has D47", "No D47"),
    hover_text = paste0(
      "Horizon: ", MLA_horizon_id,
      "<br>Strat height: ", round(strat_height_m,1),
      "<br>D47 status: ", D47_status,
      "<br>IPL T47: ", round(IPLD47_mean_T47_C,1),
      "<br>CU T47: ", round(CU_mean_T47_C,1)
    )
  ) %>%
  ggplot(
    aes(
      x = 1,
      y = strat_height_m,
      color = D47_status,
      text = hover_text
    )
  ) +
  geom_jitter(width = 0.15, height = 0) +
  scale_color_manual(
    values = c(
      "Has D47" = "black",
      "No D47" = "grey80"
    )
  ) +
  labs(
    title = "D47 coverage",
    x = "",
    y = "Stratigraphic height (m)"
  ) +
  theme_classic()+ scale_y_continuous(
    breaks = seq(500, 2300, by = 100)
  )


ggplotly(p_D47, tooltip = "text")


# ---- Interactive jitter plot: D17O availability ----
p_D17O_availability <- data_availability %>%
  mutate(
    D17O_status = if_else(has_IPL17O, "Has D17O", "No D17O"),
    hover_text = paste0(
      "Horizon: ", MLA_horizon_id,
      "<br>Strat height: ", strat_height_m, " m",
      "<br>D17O status: ", D17O_status,
      "<br>n analyses: ", IPL17O_n_analyses,
      "<br>Mean Dp17Ocarb: ", round(IPL17O_mean_Dp17Ocarb, 1),
      "<br>SD Dp17Ocarb: ", round(IPL17O_sd_Dp17Ocarb, 1),
      "<br>Mean d18Ocarb: ", round(IPL17O_mean_d18Ocarb, 2),
      "<br>Data combo: ", data_combo
    )
  ) %>%
  ggplot(aes(
    x = D17O_status,
    y = strat_height_m,
    color = D17O_status,
    text = hover_text
  )) +
  geom_jitter(
    width = 0.18,
    height = 0,
    size = 2.4,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = c(
      "Has D17O" = "black",
      "No D17O" = "grey75"
    )
  ) +
  scale_y_continuous(
    breaks = seq(500, 2300, by = 100)
  ) +
  labs(
    title = "D17O availability by horizon",
    x = NULL,
    y = "Stratigraphic height (m)"
  ) +
  theme_classic() +
  theme(legend.position = "none")

plotly::ggplotly(p_D17O_availability, tooltip = "text")