
# 01_build_BHB_dataset.R
# Purpose: Load, clean, summarize, and combine BHB proxy datasets
# Outputs: Clean proxy-specific summaries + BHB multiproxy master by strat


# 1. Setup ----------------------------------------------------

# source(here::here("scripts", "00_setup.R"))

# Helper function to round stratigraphic depths to the nearest 0.1 m
round_depth <- function(df) {
  df %>% mutate(strat_height_m = round(strat_height_m, 1))
}

# 2. Load raw data --------------------------------------------

# IPL 17O data
IPL_D17O_data  <- read.csv(here("data", "raw", "all_data_pre-Oct2025_PaleogeneBHB_IPL17O_standardized_columns.csv"))

names(IPL_D17O_data)


# IPL D47 data

# Load the raw IPL clumped isotope dataset
IPL_D47_data <- read.csv(
  here("data", "raw", "BHB Paleogene Summary May 2026.csv")
)

# Display column names to inspect dataset structure
names(IPL_D47_data)

# Some replicate analyses were assigned sample IDs ending in
# "-01", "-02", or "-03" (e.g., PK95-SC-179-01, PK95-SC-179-02).
# Remove these suffixes so all replicates share a common Sample.ID
# and can be grouped together during summarization.
#
# Only the specific suffixes -01, -02, and -03 are removed.
# Other suffixes (e.g., -04, -05, etc.) are retained.
IPL_D47_data$Sample.ID <- sub(
  "-(01|02|03)$",
  "",
  IPL_D47_data$Sample.ID
)

# Count the number of analyses associated with each Sample.ID
# after replicate suffixes have been standardized.
table(IPL_D47_data$Sample.ID)


# CU Bouder data (Havranek, 2023)

CU_data <- read_csv(here("data", "raw", "PETM_clumped.csv")) %>%
  #clean_names() %>%
  rename(strat_height_m = Strat_m_Bowen) %>%
  mutate(d18Oc_SMOW = to_VSMOW(d18Ocarb_VPBD, eq = "IUPAC")) %>%
  round_depth()

# UC Santa Cruz data (Snell et al., 2013)

Snell2013 <- read_csv(here("data", "raw", "SnellEtAl2013_summary.csv")) 


# Paul Koch 1992/1995 data
koch <- read_csv(here("data", "raw", "Koch_SC_nodules_isotopes.csv")) %>%
  #clean_names() %>%
  rename(strat_height_m = Strat_m) %>%
  round_depth()

# Bowen et al. (2001) PB data
# Bowen carbonate data
Bowen <- read.csv(here("data", "raw", "Bowen2001_IsotopeData.csv")) %>%
  rename(strat_height = Level) %>%
  mutate(
    d13C_vpdb_bowen       = as.numeric(d13C_VPDB),
    d18Ocarb_vpdb_bowen    = as.numeric(d18Ocarb_VPDB),
    d18Ocarb_vsmow_bowen   = as.numeric(d18Ocarb_VSMOW)
  )


# 3. Clean IPL Δ′17O data -------------------------------------


# check for mismatch issues or outliers
# Plot a histogram of mismatch values to visualize the overall distribution
hist(IPL_D17O_data$X33_mismatch,
     main = "Histogram of X33 Mismatch",
     xlab = "X33 Mismatch",
     col = "skyblue",
     border = "white")

# Scatterplot to examine relationship between mismatch and Δ′17O values
plot(IPL_D17O_data$X33_mismatch ~ IPL_D17O_data$Dp17Ocarb_permeg_final_correction,
     main = "Mismatch vs Δ′17Ocarb",
     xlab = expression(Delta * minute^17 * O[carb] ~ "(per meg)"),
     ylab = "X33 Mismatch",
     pch = 19, col = "gray40")

# Manually flag samples with high mismatch values (by IPL_num)
high_mismatch <- c("5699", "5841", "5830")  


# Filter out high mismatch samples by IPL_num
IPL_D17O_data_clean <- IPL_D17O_data[!IPL_D17O_data$IPL_num %in% high_mismatch, ]

# Scatterplot to examine relationship between mismatch and Δ′17O values
plot(IPL_D17O_data_clean$X33_mismatch ~ IPL_D17O_data_clean$Dp17Ocarb_permeg_final_correction,
     main = "Mismatch vs Δ′17Ocarb",
     xlab = expression(Delta * minute^17 * O[carb] ~ "(per meg)"),
     ylab = "X33 Mismatch",
     pch = 19, col = "gray40")


# Identify Δ'17O outliers automatically using boxplot rule
outliers <- boxplot.stats(IPL_D17O_data_clean$Dp17Ocarb_permeg_final_correction)$out

# Remove Δ'17O outliers
IPL_D17O_data_final <- IPL_D17O_data_clean[!IPL_D17O_data_clean$Dp17Ocarb_permeg_final_correction %in% outliers, ]

# Outlier Identification and Removal Justification
#
# Outliers in Δ′17Ocarb values were identified statistically using the 
# standard 1.5×IQR rule implemented in boxplot.stats(). This is a widely 
# accepted exploratory approach across the geosciences and other fields 
# for detecting anomalous data points that fall beyond 1.5 times the 
# interquartile range (IQR) from the 25th or 75th percentile.
#
# The flagged values were inspected alongside analytical metadata 
# (reactor mismatch, replicate consistency, and data quality flags).
# Outliers were removed only when independent evidence (e.g., high 
# X33_mismatch or analytical notes) supported their exclusion as 
# likely analytical artifacts rather than true geochemical variability.
#
# This approach ensures that the remaining dataset represents 
# reproducible Δ′17Ocarb values within analytical uncertainty while 
# preserving natural variability. All removed values are archived and 
# documented for reproducibility.


# Remove singletons ?
# IPL17O_summary <- IPL17O_summary %>% filter(n > 1)



# 4. Clean IPL Δ47 data ---------------------------------------

# Count analyses per Sample.ID to identify special qualifiers in sample names (e.g.,"SPAR")
table(IPL_D47_data$Sample.ID)

# Extract analyses of sparitic calcite 
# These samples contain "SPAR" in the Sample.ID and are
# treated separately from the primary paleosol dataset.
SPAR_D47_data <- IPL_D47_data %>%
  filter(grepl("SPAR", Sample.ID))

# Retain only analyses of micritic or microsparitic calcite
# Any sample containing "SPAR" is excluded.
D47_primary_data <- IPL_D47_data %>%
  filter(!grepl("SPAR", Sample.ID))

# Plot replicate-level Δ47-derived temperatures versus stratigraphic position
#
# Purpose:
# Visualize all individual Δ47 temperature replicates from primary (non-SPAR)
# pedogenic carbonate samples to assess temperature variability with stratigraphy
# and identify potential outliers or stratigraphic trends.
ggplot(D47_primary_data,
       aes(x = T.D47..Petersen,
           y = Strat)) +
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
           Sample.ID,
           T.D47..Petersen,
           FUN = median,
           na.rm = TRUE
         ),
         y = T.D47..Petersen,
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
    !is.na(T.D47..Petersen),
    !is.na(Strat)
  ) %>%
  
  # Group replicate analyses by stratigraphic level and session
  group_by(Strat, Session) %>%
  
  # Calculate mean temperature for each strat level in each session
  summarise(
    T_mean = mean(T.D47..Petersen, na.rm = TRUE),
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
           y = Strat)) +
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

ggplot(D47_primary_data,
       aes(x = reorder(Sample.ID,
                       T.D47..Petersen,
                       median),
           y = T.D47..Petersen)) +
  geom_point() +
  coord_flip()

sample_sd <- D47_primary_data %>%
  group_by(Sample.ID) %>%
  summarise(
    n = n(),
    mean_T = mean(T.D47..Petersen),
    sd_T = sd(T.D47..Petersen),
    range_T = max(T.D47..Petersen) -
      min(T.D47..Petersen)
  )

summary(sample_sd$sd_T)
summary(sample_sd$range_T)

sample_sd %>%
  arrange(desc(range_T))


ggplot(sample_sd,
       aes(x = reorder(Sample.ID, sd_T),
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

ggplot(D47_primary_data,
       aes(
         x = reorder(
           Sample.ID,
           T.D47..Petersen,
           FUN = median,
           na.rm = TRUE
         ),
         y = T.D47..Petersen,
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

ggplot(D47_primary_data,
       aes(x = Sample.ID,
           y = T.D47..Petersen,
           color = Session)) +
  geom_jitter(width = 0.15) +
  coord_flip() +
  theme_classic()


# 5. Clean Snell et al. (2013) Δ47 data --------------------------------

str(Snell2013)

# Separate Snell et al. (2013) samples by sample type
#
# Paleosol micrite = primary paleoclimate target
# Altered carbonate, fracture spar, and bivalves are retained but separated
# because they are not equivalent to primary paleosol micrite.

Snell2013 <- Snell2013 %>%
  mutate(
    Sample_Group = case_when(
      Sample_Type == "Paleosol Micrite" ~ "Primary paleosol micrite",
      Sample_Type == "Altered Paleosol carbonate" ~ "Altered paleosol carbonate",
      Sample_Type == "Fracture Spar" ~ "Fracture spar",
      Sample_Type == "Bivalve Fossil" ~ "Bivalve fossil",
      TRUE ~ "Other"
    )
  )

# Check grouping
table(Snell2013$Sample_Group)

# Plot Snell et al. (2013) Δ47 temperatures through time
#
# Sample types are represented by point shape so primary paleosol micrites can
# be visually distinguished from altered carbonates, fracture spars, and
# bivalve fossils.
ggplot(Snell2013,
       aes(x = Average_Temp_C,
           y = Age_Ma,
           shape = Sample_Type)) +
  
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

ggplot(Snell2013,
       aes(x = Average_Temp_C,
           y = Average_d18Ow_permil_SMOW,
           shape = Sample_Type)) +
  
  geom_errorbarh(
    aes(
      xmin = Average_Temp_C - Temp_1se,
      xmax = Average_Temp_C + Temp_1se
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
# be merged with the IPL Δ47 dataset for comparison.
#
# Altered paleosol carbonates, fracture spars, and bivalves are retained as
# separate datasets for evaluating diagenetic and non-pedogenic signatures.

Snell2013_micrite <- Snell2013 %>%
  filter(Sample_Type == "Paleosol Micrite")

Snell2013_altered <- Snell2013 %>%
  filter(Sample_Type == "Altered Paleosol carbonate")

Snell2013_spar <- Snell2013 %>%
  filter(Sample_Type == "Fracture Spar")

Snell2013_bivalve <- Snell2013 %>%
  filter(Sample_Type == "Bivalve Fossil")

# Inspect sample counts
nrow(Snell2013_micrite)
nrow(Snell2013_altered)
nrow(Snell2013_spar)
nrow(Snell2013_bivalve)

# Verify separation
table(Snell2013$Sample_Type)


# 6. Clean CU Boulder Δ47 data --------------------------------

str(CU_data)


# Plot Δ47 temperatures versus stratigraphic position
#
# Purpose:
# Visualize temperature trends through the stratigraphic section and identify
# unusually warm or cool samples that may warrant further investigation.

ggplot(CU_data,
       aes(x = T47_C,
           y = strat_height_m)) +
  
  geom_errorbarh(
    aes(
      xmin = T47_C - T47_2SE_C,
      xmax = T47_C + T47_2SE_C
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

# Plot Δ47 temperature versus reconstructed soil-water δ18O
#
# Purpose:
# Evaluate whether high temperatures are associated with anomalous δ18Ow
# values that could indicate recrystallization, diagenesis, or alteration.
ggplot(CU_data,
       aes(x = T47_C,
           y = d18Ow_VSMOW)) +
  
  # Horizontal temperature uncertainty
  geom_errorbarh(
    aes(
      xmin = T47_C - T47_2SE_C,
      xmax = T47_C + T47_2SE_C
    ),
    height = 0,
    alpha = 0.5
  ) +
  
  # Vertical water-isotope uncertainty
  geom_errorbar(
    aes(
      ymin = d18Ow_VSMOW - d18Ow_SE_VSMOW,
      ymax = d18Ow_VSMOW + d18Ow_SE_VSMOW
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

# Koch 1992/1995 data

str(koch)

plot(koch$d18Ocarb_VPDB, koch$d13C_VPDB)

plot(koch$d18Ocarb_VPDB, koch$strat_height_m)

# Plot replicate δ18O values by horizon
#
# Purpose:
# Visualize within-horizon variability in Koch et al. paleosol carbonate δ18O.
# Each point is an individual analysis; points are grouped by horizon_id.

ggplot(
  koch %>% filter(!is.na(d18Ocarb_VSMOW)),
  aes(
    x = reorder(
      horizon_id,
      d18Ocarb_VSMOW,
      FUN = median,
      na.rm = TRUE
    ),
    y = d18Ocarb_VSMOW
  )
) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.75
  ) +
  coord_flip() +
  labs(
    x = "Horizon ID",
    y = expression(delta^{18} * O[carb] ~ "(‰ VSMOW)")
  ) +
  theme_classic()

# Bowen PB data

# Quick sanity-check plots for Bowen carbonate isotope data

# δ13Ccarb vs. stratigraphic height
ggplot(Bowen, aes(x = d13C_vpdb_bowen, y = strat_height)) +
  geom_point(size = 2) +
  labs(x = expression(delta^13*C[carb]~"(‰ VPDB)"), y = "Strat. height (m)") +
  theme_bw()

# δ18Ocarb VSMOW vs. stratigraphic height
ggplot(Bowen, aes(x = d18Ocarb_vsmow_bowen, y = strat_height)) +
  geom_point(size = 2) +
  labs(x = expression(delta^18*O[carb]~"(‰ VSMOW)"), y = "Strat. height (m)") +
  theme_bw()

# δ18Ocarb VPDB vs. stratigraphic height
ggplot(Bowen, aes(x = d18Ocarb_vpdb_bowen, y = strat_height)) +
  geom_point(size = 2) +
  labs(x = expression(delta^18*O[carb]~"(‰ VPDB)"), y = "Strat. height (m)") +
  theme_bw()

# Histogram of δ18Ocarb VSMOW values
ggplot(Bowen, aes(x = d18Ocarb_vsmow_bowen)) +
  geom_histogram(bins = 20) +
  labs(x = expression(delta^18*O[carb]~"(‰ VSMOW)"), y = "Count") +
  theme_bw()

# Histogram of δ18Ocarb VPDB values
ggplot(Bowen, aes(x = d18Ocarb_vpdb_bowen)) +
  geom_histogram(bins = 20) +
  labs(x = expression(delta^18*O[carb]~"(‰ VPDB)"), y = "Count") +
  theme_bw()

# Remove likely diagenetically altered samples.
# Bowen et al. noted that altered carbonates tend to have
# anomalously low δ18O values; here we conservatively exclude
# samples with δ18Ocarb < -11‰ VPDB.

Bowen <- Bowen %>%
  filter(d18Ocarb_vpdb_bowen > -11)


# Histogram of retained carbonate δ18O values
ggplot(Bowen, aes(x = d18Ocarb_vpdb_bowen)) +
  geom_histogram(bins = 20) +
  labs(
    x = expression(delta^18*O[carb]~"(‰ VPDB)"),
    y = "Count",
    title = expression(delta^18*O[carb]~"Distribution")
  ) +
  theme_bw()


# Stratigraphic distribution of retained carbonate δ18O values
ggplot(Bowen,
       aes(x = d18Ocarb_vpdb_bowen,
           y = strat_height)) +
  geom_point(size = 2) +
  labs(
    x = expression(delta^18*O[carb]~"(‰ VPDB)"),
    y = "Stratigraphic Height (m)",
    title = expression(delta^18*O[carb]~"vs. Stratigraphic Height")
  ) +
  theme_bw()

# 8. Summarize each dataset by meter level ---------------

# IPL 17O data
# Group by horizon_id and summarize relevant columns
generic_sd <- 12

IPL17O_summary <- IPL_D17O_data_final %>%
  group_by(horizon_id) %>%
  summarise(
    mean_dp17Ocarb = mean(dp17Ocarb_SMOWSLAP, na.rm = TRUE),
    sd_dp17Ocarb   = sd(dp17Ocarb_SMOWSLAP, na.rm = TRUE),
    se_dp17Ocarb   = sd_dp17Ocarb / sqrt(n()),
    
    mean_dp18Ocarb = mean(dp18Ocarb_SMOWSLAP, na.rm = TRUE),
    sd_dp18Ocarb   = sd(dp18Ocarb_SMOWSLAP, na.rm = TRUE),
    se_dp18Ocarb   = sd_dp18Ocarb / sqrt(n()),
    
    mean_Dp17Ocarb = mean(Dp17Ocarb_permeg_final_correction, na.rm = TRUE),
    sd_Dp17Ocarb   = sd(Dp17Ocarb_permeg_final_correction, na.rm = TRUE),
    se_Dp17Ocarb   = sd_Dp17Ocarb / sqrt(n()),
    
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    # Use generic SD of 12 if n < 2 OR observed SD < 12
    sd_Dp17Ocarb_adj = ifelse(n < 2 | sd_Dp17Ocarb < generic_sd, generic_sd, sd_Dp17Ocarb),
    se_Dp17Ocarb_adj = sd_Dp17Ocarb_adj / sqrt(n),
    ci95_Dp17Ocarb_adj = se_Dp17Ocarb_adj * 1.96
  )


# Convert delta-prime values to delta values and add to IPL17O_summary
IPL17O_summary <- IPL17O_summary %>%
  mutate(
    mean_d17Ocarb = 1000 * (exp(mean_dp17Ocarb / 1000) - 1),
    mean_d18Ocarb = 1000 * (exp(mean_dp18Ocarb / 1000) - 1)
  )



# Summarize IPL primary (non-SPAR) Δ47 data by sample
#
# Produces a sample-level temperature summary table for downstream analyses.
D47_primary_summary <- D47_primary_data %>%
  group_by(Sample.ID) %>%
  summarise(
    Strat = mean(Strat, na.rm = TRUE),
    
    n_T47 = sum(!is.na(T.D47..Petersen)),
    mean_T47 = mean(T.D47..Petersen, na.rm = TRUE),
    sd_T47 = sd(T.D47..Petersen, na.rm = TRUE),
    se_T47 = sd_T47 / sqrt(n_T47),
    
    min_T47 = min(T.D47..Petersen, na.rm = TRUE),
    max_T47 = max(T.D47..Petersen, na.rm = TRUE),
    range_T47 = max_T47 - min_T47,
    
    sessions = paste(
      sort(unique(na.omit(Session))),
      collapse = ", "
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    sample_type = "primary"
  ) %>%
  arrange(Strat)


# Summarize IPL SPAR Δ47 data by sample
#
# Kept separate from the primary micrite/microspar dataset.
D47_spar_summary <- SPAR_D47_data %>%
  group_by(Sample.ID) %>%
  summarise(
    Strat = mean(Strat, na.rm = TRUE),
    
    n_T47 = sum(!is.na(T.D47..Petersen)),
    mean_T47 = mean(T.D47..Petersen, na.rm = TRUE),
    sd_T47 = sd(T.D47..Petersen, na.rm = TRUE),
    se_T47 = sd_T47 / sqrt(n_T47),
    
    min_T47 = min(T.D47..Petersen, na.rm = TRUE),
    max_T47 = max(T.D47..Petersen, na.rm = TRUE),
    range_T47 = max_T47 - min_T47,
    
    sessions = paste(
      sort(unique(na.omit(Session))),
      collapse = ", "
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    sample_type = "spar"
  ) %>%
  arrange(Strat)

# Inspect summaries
D47_primary_summary

D47_spar_summary

# Export summary tables
write_csv(
  D47_primary_summary,
  here("data", "processed", "IPL_D47_primary_summary.csv")
)

write_csv(
  D47_spar_summary,
  here("data", "processed", "IPL_D47_spar_summary.csv")
)







# summarize the koch data

# Summarize Koch data by horizon_id
koch_summary <- koch %>%
  group_by(horizon_id, strat_height_m) %>%  # Standardized ID for joins
  summarise(
    # δ13Ccarb summary
    d13Ccarb_mean_koch = mean(d13C_VPDB, na.rm = TRUE),
    d13Ccarb_se_koch   = sd(d13C_VPDB, na.rm = TRUE) / sqrt(sum(!is.na(d13C_VPDB))),
    n_d13Ccarb_koch    = sum(!is.na(d13C_VPDB)),
    
    # δ18Ocarb (VPDB) summary
    d18Ocarb_vpdb_mean_koch = mean(d18Ocarb_VPDB, na.rm = TRUE),
    d18Ocarb_vpdb_se_koch   = sd(d18Ocarb_VPDB, na.rm = TRUE) / sqrt(sum(!is.na(d18Ocarb_VPDB))),
    n_d18Ocarb_vpdb_koch    = sum(!is.na(d18Ocarb_VPDB)),
    
    # δ18Ocarb (VSMOW) summary
    d18Ocarb_vsmow_mean_koch = mean(d18Ocarb_VSMOW, na.rm = TRUE),
    d18Ocarb_vsmow_se_koch   = sd(d18Ocarb_VSMOW, na.rm = TRUE) / sqrt(sum(!is.na(d18Ocarb_VSMOW))),
    n_d18Ocarb_vsmow_koch    = sum(!is.na(d18Ocarb_VSMOW)),
    
    .groups = "drop"
  )

# Quick sanity-check plot of mean carbonate δ18O values vs. stratigraphic height
ggplot(koch_summary,
       aes(x = d18Ocarb_vsmow_mean_koch,
           y = strat_height_m)) +
  geom_point(size = 2) +
  geom_errorbarh(
    aes(
      xmin = d18Ocarb_vsmow_mean_koch - d18Ocarb_vsmow_se_koch,
      xmax = d18Ocarb_vsmow_mean_koch + d18Ocarb_vsmow_se_koch
    ),
    height = 0
  ) +
  labs(
    x = expression(delta^18 * O[carb] ~ "(‰ VSMOW)"),
    y = "Stratigraphic Height (m)",
    title = expression(delta^18 * O[carb] ~ "vs. Stratigraphic Height")
  ) +
  theme_bw()

# quick plot for sanity check, error is d13Ccarb_se_koch 
plot(koch_summary$d13Ccarb_mean_koch, koch_summary$strat_height)

# Quick sanity-check plot of mean carbonate δ13C values vs. stratigraphic height
ggplot(koch_summary,
       aes(x = d13Ccarb_mean_koch,
           y = strat_height_m)) +
  geom_errorbarh(
    aes(
      xmin = d13Ccarb_mean_koch - d13Ccarb_se_koch,
      xmax = d13Ccarb_mean_koch + d13Ccarb_se_koch
    ),
    height = 0
  ) +
  geom_point(size = 3) +
  labs(
    x = expression(delta^13 * C[carb] ~ "(‰ VPDB)"),
    y = "Stratigraphic Height (m)",
    title = expression("Mean " * delta^13 * C[carb] ~ "Values")
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )




# Bowen 
# Summarize Bowen data by strat height, output column = Sample_ID
bowen_summary <- Bowen %>%
  group_by(horizon_id, strat_height) %>%
  summarise(
    # δ13Ccarb summary
    d13Ccarb_mean_bowen = mean(d13C_vpdb_bowen, na.rm = TRUE),
    d13Ccarb_se_bowen   = sd(d13C_vpdb_bowen, na.rm = TRUE) / sqrt(sum(!is.na(d13C_vpdb_bowen))),
    n_d13Ccarb_bowen    = sum(!is.na(d13C_vpdb_bowen)),
    
    # δ18Ocarb (VPDB) summary
    d18Ocarb_vpdb_mean_bowen = mean(d18Ocarb_vpdb_bowen, na.rm = TRUE),
    d18Ocarb_vpdb_se_bowen  = sd(d18Ocarb_vpdb_bowen, na.rm = TRUE) / sqrt(sum(!is.na(d18Ocarb_vpdb_bowen))),
    n_d18Ocarb_vpdb_bowen    = sum(!is.na(d18Ocarb_vpdb_bowen)),
    
    # δ18Ocarb (VSMOW) summary
    d18Ocarb_vsmow_mean_bowen = mean(d18Ocarb_vsmow_bowen, na.rm = TRUE),
    d18Ocarb_vsmow_se_bowen   = sd(d18Ocarb_vsmow_bowen, na.rm = TRUE) / sqrt(sum(!is.na(d18Ocarb_vsmow_bowen))),
    n_d18Ocarb_vsmow_bowen    = sum(!is.na(d18Ocarb_vsmow_bowen)),
    
    .groups = "drop"
  )

# 9. Combine summaries by horizon / strat ----------------------




# 10. Save proxy-specific processed outputs ---------------------

write_csv(IPL_D17O_summary, here("data", "processed", "IPL_D17O_summary.csv"))
write_csv(IPL_D47_summary,  here("data", "processed", "IPL_D47_summary.csv"))
write_csv(CU_D47_summary,   here("data", "processed", "CU_D47_summary.csv"))
write_csv(carb_summary,     here("data", "processed", "Koch_Bowen_carb_summary.csv"))


# 11. Save combined BHB dataset ---------------------------------

write_csv(
  BHB_multiproxy_by_strat,
  here("data", "processed", "BHB_multiproxy_by_strat.csv")
)


# 12. Quick checks ---------------------------------------------

glimpse(BHB_multiproxy_by_strat)

BHB_multiproxy_by_strat %>%
  summarise(
    n_horizons = n(),
    min_strat = min(Strat_m, na.rm = TRUE),
    max_strat = max(Strat_m, na.rm = TRUE)
  )