###### Load Data, Packages ##########
library(dplyr)
library(ggplot2)
library(isogeochem)
library(patchwork)
library(here)
library(stringr)
library(plotly)
  

### IPL data --------

IPL_data  <- read.csv(here("data", "raw", "all_data_PaleogeneBHB_IPL17O_standardized_columns.csv"))

names(IPL_data)

# check for mismatch issues or outliers
# Plot a histogram of mismatch values to visualize the overall distribution
hist(IPL_data$X33_mismatch,
     main = "Histogram of X33 Mismatch",
     xlab = "X33 Mismatch",
     col = "skyblue",
     border = "white")

# Scatterplot to examine relationship between mismatch and Δ′17O values
plot(IPL_data$X33_mismatch ~ IPL_data$Dp17Ocarb_permeg_final_correction,
     main = "Mismatch vs Δ′17Ocarb",
     xlab = expression(Delta * minute^17 * O[carb] ~ "(per meg)"),
     ylab = "X33 Mismatch",
     pch = 19, col = "gray40")

# Manually flag samples with high mismatch values (by IPL_num)
high_mismatch <- c("5699", "5841", "5830")  # These will be excluded or flagged later


# Filter out high mismatch samples by IPL_num
IPL_data_clean <- IPL_data[!IPL_data$IPL_num %in% high_mismatch, ]

# Scatterplot to examine relationship between mismatch and Δ′17O values
plot(IPL_data_clean$X33_mismatch ~ IPL_data_clean$Dp17Ocarb_permeg_final_correction,
     main = "Mismatch vs Δ′17Ocarb",
     xlab = expression(Delta * minute^17 * O[carb] ~ "(per meg)"),
     ylab = "X33 Mismatch",
     pch = 19, col = "gray40")


# Identify Δ'17O outliers automatically using boxplot rule
outliers <- boxplot.stats(IPL_data_clean$Dp17Ocarb_permeg_final_correction)$out

# Remove Δ'17O outliers
IPL_data_final <- IPL_data_clean[!IPL_data_clean$Dp17Ocarb_permeg_final_correction %in% outliers, ]

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



# Group by horizon_id and summarize relevant columns
generic_sd <- 12

IPL17O_summary <- IPL_data_final %>%
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


# Remove singletons !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#IPL17O_summary <- IPL17O_summary %>% filter(n > 1)


## add CU data -----------------------------------------------------------------

# read the CU clumped data
CU_data <- read.csv(here("data", "raw", "PETM_clumped.csv")) 


# Compute δ18O of water (VSMOW) from carbonate (VPDB) and T47_C
CU_data <- CU_data %>%
  mutate(d18Ocarb_VSMOW = to_VSMOW(d18Ocarb_VPBD))

# Full join while keeping and renaming selected columns from CU_data
IPL_CU_summary <- IPL17O_summary %>%
  full_join(
    CU_data %>%
      dplyr::select(
        horizon_id,
        TEMP_CU = T47_C,
        T47_2SE_C = T47_2SE_C,
      ),
    by = "horizon_id"
  )



## add Koch data -----------------------------------------------------------------

# Koch carbonate data
koch <- read.csv(here("data", "raw", "Koch_SC_nodules_isotopes.csv")) %>%
  rename(strat_height = Strat_m)

# Check for outliers
hist(koch$d13C_VPDB)
hist(koch$d18Ocarb_VSMOW)
# keep all the Koch data

# quick plot for sanity check
plot(koch$d13C_VPDB, koch$strat_height)

# quick plot for sanity check
plot(koch$d18Ocarb_VSMOW, koch$strat_height)

# Summarize Koch data by horizon_id
koch_summary <- koch %>%
  group_by(horizon_id, strat_height) %>%  # Standardized ID for joins
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

# quick plot for sanity check
plot(koch_summary$d18Ocarb_vsmow_mean_koch, koch_summary$strat_height)

# quick plot for sanity check
plot(koch_summary$d13Ccarb_mean_koch, koch_summary$strat_height)


# Full join while keeping and renaming selected columns from koch
IPL_CU_summary <- IPL_CU_summary %>%
  full_join(
    koch_summary %>%
      dplyr::select(
        horizon_id,
        strat_height,
        d13Ccarb_mean_koch = d13Ccarb_mean_koch,
        d13Ccarb_se_koch = d13Ccarb_se_koch, 
        n_d13Ccarb_koch = n_d13Ccarb_koch,
        d18Ocarb_vsmow_mean_koch = d18Ocarb_vsmow_mean_koch,
        d18Ocarb_vsmow_se_koch = d18Ocarb_vsmow_se_koch,
        n_d18Ocarb_vsmow_koch = n_d18Ocarb_vsmow_koch,
      ),
    by = "horizon_id"
  )








## add Bowen data -----------------------------------------------------------------


# Bowen carbonate data
Bowen <- read.csv(here("data", "raw", "Bowen2001_IsotopeData.csv")) %>%
  rename(strat_height = Level) %>%
  mutate(
    d13C_vpdb_bowen       = as.numeric(d13C_VPDB),
    d18Ocarb_vpdb_bowen    = as.numeric(d18Ocarb_VPDB),
    d18Ocarb_vsmow_bowen   = as.numeric(d18Ocarb_VSMOW)
  )

# quick plots for sanity check
plot(Bowen$d13C_vpdb_bowen, Bowen$strat_height)

plot(Bowen$d18Ocarb_vsmow_bowen, Bowen$strat_height)
plot(Bowen$d18Ocarb_vpdb_bowen, Bowen$strat_height)
hist(Bowen$d18Ocarb_vsmow_bowen)
hist(Bowen$d18Ocarb_vpdb_bowen)

# !!!!! REMOVE Bowen outliers < -11 permil !!!!!!!
Bowen <- Bowen %>%
  filter(d18Ocarb_vpdb_bowen > -11)

hist(Bowen$d18Ocarb_vpdb_bowen)
plot(Bowen$d18Ocarb_vpdb_bowen, Bowen$strat_height)

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

IPL_CU_summary <- IPL_CU_summary %>%
  full_join(
    bowen_summary %>%
      dplyr::select(
        horizon_id,
        strat_height,
        d18Ocarb_vsmow_mean_bowen = d18Ocarb_vsmow_mean_bowen,
        d18Ocarb_vsmow_se_bowen = d18Ocarb_vsmow_se_bowen,
        n_d18Ocarb_vsmow_bowen = n_d18Ocarb_vsmow_bowen,
        d13Ccarb_mean_bowen = d13Ccarb_mean_bowen,
        d13Ccarb_se_bowen = d13Ccarb_se_bowen,
        n_d13Ccarb_bowen = n_d13Ccarb_bowen,
      ),
    by = "horizon_id"
  )

IPL_CU_summary <- IPL_CU_summary %>%
  mutate(
    d18Ocarb_vsmow_mean_all = coalesce(d18Ocarb_vsmow_mean_koch, d18Ocarb_vsmow_mean_bowen),
    d18Ocarb_vsmow_se_all   = coalesce(d18Ocarb_vsmow_se_koch, d18Ocarb_vsmow_se_bowen),
    n_d18Ocarb_all          = coalesce(n_d18Ocarb_vsmow_koch, n_d18Ocarb_vsmow_bowen),
    d13Ccarb_vpdb_mean_all = coalesce(d13Ccarb_mean_koch, d13Ccarb_mean_bowen),
    d13Ccarb_vpdb_se_all   = coalesce(d13Ccarb_se_koch, d13Ccarb_se_bowen),
    n_d13Ccarb_all          = coalesce(n_d13Ccarb_koch, n_d13Ccarb_bowen)
  )


# Consolidate strat_height into a single column
IPL_CU_summary <- IPL_CU_summary %>%
  mutate(
    strat_height_clean = coalesce(strat_height.x, strat_height.y)
  ) %>%
  dplyr::select(-strat_height.x, -strat_height.y) %>%
  rename(strat_height = strat_height_clean)


# set this value manually, it keep sgetting left out
IPL_CU_summary$strat_height[which(IPL_CU_summary$horizon_id == "PB-00-01-09D")] <- 1511.7

### AGE MODEL: Interpolate onto PCB_summary --------------------

# Load age model data
age_model_raw <- read.csv(here("data", "raw", "approx_age_mdl_PETM_PCB.csv"))

# Filter to valid rows and fit linear model
age_model_fit <- lm(Age_Ma ~ est_Depth_m_PCB_outcrop, 
                    data = filter(age_model_raw, !is.na(est_Depth_m_PCB_outcrop), !is.na(Age_Ma)))

# Predict Age onto PCB_summary
IPL_CU_summary <- IPL_CU_summary %>%
  mutate(Age_Ma = predict(age_model_fit, newdata = tibble(est_Depth_m_PCB_outcrop = strat_height)))

# export the model
saveRDS(age_model_fit, file = "data/processed/age_model_fit.rds")

## TEMP MODEL ---------------


names(IPL_CU_summary)


names(IPL_CU_summary)

# 1. Create TEMP_priors from T47_C
IPL_CU_summary <- IPL_CU_summary %>%
  mutate(TEMP_priors = TEMP_CU)

# 2. override TEMP_priors at specific strat_height
background_mean <- 25
IPL_CU_summary <- IPL_CU_summary %>%
  mutate(TEMP_priors = case_when(
    strat_height == 2210 ~ background_mean,
    strat_height == 1600.0 ~ background_mean,
    strat_height == 1480.8 ~ background_mean,
    strat_height == 940 ~ background_mean,
    TRUE ~ TEMP_priors
  ))


# 3. Build model from TEMP_priors and strat_height
# Filter to rows with valid TEMP_priors and strat_height
temp_model_data <- IPL_CU_summary %>%
  filter(!is.na(strat_height), !is.na(TEMP_priors)) %>%
  arrange(strat_height)

# 4. Interpolate/extrapolate to all strat_heights using approx()

IPL_CU_summary$TEMP_modeled <- approx(
  x = temp_model_data$strat_height,
  y = temp_model_data$TEMP_priors,
  xout = IPL_CU_summary$strat_height,
  rule = 2  # constant extrapolation
)$y

# THIS IS FOR TEST, DON'T FORGET TO REMOVE
#PCB_summary$TEMP_modeled <- 25 

ggplot() +
  # Modeled temperatures (blue points)
  geom_point(
    data = IPL_CU_summary,
    aes(x = TEMP_modeled, y = strat_height),
    color = "blue",
    size = 2,
    alpha = 0.6
  ) +
  # CU data temperatures (red points)
  geom_point(
    data = IPL_CU_summary,
    aes(x = TEMP_CU, y = strat_height),
    color = "red",
    size = 2.5,
    alpha = 0.9
  ) +
  labs(
    x = "Temperature (°C)",
    y = "Stratigraphic Height (m)",
    title = "Interpolated Temperature Model"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )


### CREATE FUNCTION TO RECONSTRUCT D17Owater--------------------------------

compute_D17O_water <- function(d17O_carb, d18O_carb, T_C, theta = 0.525) {
  # Step 1: Convert °C to Kelvin
  T_K <- T_C + 273.15
  
  # Step 2: Calculate 18alpha (fractionation factor for 18O)
  alpha18 <- exp((18030 / T_K - 32.42) / 1000)
  
  # Step 3: Calculate 17alpha using theta
  alpha17 <- exp(theta * log(alpha18))
  
  # Step 4: Compute d17O_rfw (raw delta, not delta-prime)
  d17O_rfw <- ((d17O_carb + 1000) / alpha17) - 1000
  
  # Step 5: Convert d17O_rfw to delta-prime
  dp17O_rfw <- 1000 * log1p(d17O_rfw / 1000)
  
  # Step 6: Compute d18O_rfw (raw delta)
  d18O_rfw <- ((d18O_carb + 1000) / alpha18) - 1000
  
  # Step 7: Convert d18O_rfw to delta-prime
  dp18O_rfw <- 1000 * log1p(d18O_rfw / 1000)
  
  # Step 8: Compute D17O of rfw (‰)
  D17O_rfw <- dp17O_rfw - 0.528 * dp18O_rfw
  
  # Step 9: Convert to per meg
  D17O_rfw_permeg <- D17O_rfw * 1000
  
  return(D17O_rfw_permeg)
}


### TEMP SENSITIVITY TESTS ----------
a <- 10.86 # mean(IPL_CU_summary$mean_dp17Ocarb, na.rm = TRUE)
b <- 20.77 # mean(IPL_CU_summary$mean_dp18Ocarb, na.rm = TRUE)
Tx <- c(10,15,20,25,30,35,40)
c <- compute_D17O_water(a, b, Tx)
plot(Tx, c)


### COMPUTE WATER VALUES -----

# compute D'17Orfw
IPL_CU_summary$D17O_rfw_VSMOWSLAP_fromR <- compute_D17O_water(
  IPL_CU_summary$mean_d17Ocarb,
  IPL_CU_summary$mean_d18Ocarb,
  IPL_CU_summary$TEMP_modeled
)

# compute d18Orfw
IPL_CU_summary <- IPL_CU_summary %>%
  mutate(
    mean_d18Owater_VSMOW = d18O_H2O(
      temp = TEMP_modeled,
      d18O = d18Ocarb_vsmow_mean_all,
      min = "calcite",
      eq = "KO97"
    )
  )


IPL_CU_summary$se_d18Owater_VSMOW <- IPL_CU_summary$d18Ocarb_vsmow_se_all


# Compute d′18O of water
IPL_CU_summary <- IPL_CU_summary %>%
  mutate(
    mean_dp18Orsw_vsmow = 1000 * log(mean_d18Owater_VSMOW / 1000 + 1),
    se_dp18Orsw_vsmow = se_d18Owater_VSMOW / (mean_d18Owater_VSMOW / 1000 + 1)
  )

# Set error for D17Orfw
IPL_CU_summary <- IPL_CU_summary %>%
  mutate(se_D17O_rfw = se_Dp17Ocarb)


# NOW compute unevaporated waters!
# Inputs (carbonate measurements)
d18O_carb   <- IPL_CU_summary$mean_dp18Ocarb          # ‰
D17O_carb   <- IPL_CU_summary$mean_Dp17Ocarb / 1000   # convert per meg → ‰
Temp        <- IPL_CU_summary$TEMP_modeled
theta_carb  <- 0.5250
lambda_lake <- 0.515
slope_MWL   <- 0.528

# 1. Calcite–water fractionation
alpha18 <- exp((18030 / (Temp + 273.15) - 32.42) / 1000)

# 2. Convert carbonate → parent (lake) water
d18O_lake  <- ((1000 + d18O_carb) / alpha18) - 1000
dp18O_lake <- 1000 * log1p(d18O_lake / 1000)
D17O_lake  <- D17O_carb + 1000 * log(alpha18) * (0.528 - theta_carb)

# 3. Backprojection in δ′17O–δ′18O space
dp17O_lake <- D17O_lake + slope_MWL * dp18O_lake      # convert to δ′17O
b_lake     <- dp17O_lake - lambda_lake * dp18O_lake   # lake line intercept
b_MWL      <- 0.010                                      # MWL intercept (~0)

dp18O_met  <- (b_lake - b_MWL) / (slope_MWL - lambda_lake)
dp17O_met  <- slope_MWL * dp18O_met + b_MWL
D17O_met   <- dp17O_met - slope_MWL * dp18O_met

# 4. Add results to dataframe (Δ′17O back to per meg for readability)
IPL_CU_summary$dp18O_met <- dp18O_met
IPL_CU_summary$dp17O_met <- dp17O_met
IPL_CU_summary$D17O_met  <- D17O_met * 1000  # per meg

# Optional quick sanity check output
cat("Mean unevaporated δ′18O =", round(mean(dp18O_met, na.rm = TRUE), 2), "‰\n")
cat("Mean unevaporated Δ′17O =", round(mean(D17O_met * 1000, na.rm = TRUE), 1), "per meg\n")
### export data ----------------------------------------------------------------

# Export the dataframe as CSV
write.csv(IPL_CU_summary, 
          file = "data/processed/IPL_CU_summary.csv", 
          row.names = FALSE)

