# 06_process_reference_datasets.R
# Purpose:
#   Import, clean, summarize, visualize, and export published datasets used
#   to compare with the northern Bighorn Basin multiproxy record.
#
# Reference datasets:
#   1. Harper et al. (2024): marine CO2 and SST
#   2. Kelson et al.: Tornillo Basin carbonate clumped isotopes
#   3. Fricke et al. (1998): mammal and gar apatite isotope records
#   4. Wing et al. (2000): leaf-margin MAT estimates
#   5. GTS2020: updated BHB age-model tie points
#
# Main outputs:
#   - Harper2024_CO2_SST_processed.csv
#   - Kelson_Tornillo_D47_processed.csv
#   - Kelson_Tornillo_D47_strat_age_filtered.csv
#   - FrickeEtAl1998_specimen_means.csv
#   - FrickeEtAl1998_zone_means.csv
#   - FrickeEtAl1998_temperature_change.csv
#   - WingEtAl2000_LMA_MAT_processed.csv
#   - BHB_GTS2020_tiepoints.csv
#   - BHB_GTS2020_age_model.csv

# ---- Load packages ----
library(tidyverse)
library(here)
library(zoo)
library(readr)

# ---- Create output directories ----
processed_dir <- here("data", "processed")
figure_dir <- here("figures", "reference_datasets")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# - 1. Harper et al. (2024): CO2 and SST reference records -----


Harper2024_CO2_raw <- read_csv(
  here("data", "raw", "HarperEtAl2024_co2_out.csv"),
  show_col_types = FALSE
)

Harper2024_SST_raw <- read_csv(
  here("data", "raw", "HarperEtAl2024_sst_out.csv"),
  show_col_types = FALSE
)

# Join the records by age, calculate three-point anchors, and then fit
# lightly smoothed splines for plotting broad trends.
Harper2024_CO2_SST <- Harper2024_CO2_raw %>%
  transmute(
    Age_Ma = age / 1000,
    Harper2024_mean_CO2_ppm = mean,
    Harper2024_CO2_lower95_ppm = `2.50%`,
    Harper2024_CO2_upper95_ppm = `97.50%`
  ) %>%
  left_join(
    Harper2024_SST_raw %>%
      transmute(
        Age_Ma = age / 1000,
        Harper2024_mean_SST_C = mean,
        Harper2024_SST_lower95_C = `2.50%`,
        Harper2024_SST_upper95_C = `97.50%`
      ),
    by = "Age_Ma"
  ) %>%
  arrange(Age_Ma) %>%
  mutate(
    Harper2024_SST_anchor_C = rollapply(
      Harper2024_mean_SST_C,
      width = 3,
      FUN = mean,
      fill = NA_real_,
      align = "center",
      partial = TRUE,
      na.rm = TRUE
    ),
    Harper2024_CO2_anchor_ppm = rollapply(
      Harper2024_mean_CO2_ppm,
      width = 3,
      FUN = mean,
      fill = NA_real_,
      align = "center",
      partial = TRUE,
      na.rm = TRUE
    )
  )

Harper_SST_spline <- smooth.spline(
  x = Harper2024_CO2_SST$Age_Ma,
  y = Harper2024_CO2_SST$Harper2024_SST_anchor_C,
  spar = 0.35
)

Harper_CO2_spline <- smooth.spline(
  x = Harper2024_CO2_SST$Age_Ma,
  y = Harper2024_CO2_SST$Harper2024_CO2_anchor_ppm,
  spar = 0.35
)

Harper2024_CO2_SST <- Harper2024_CO2_SST %>%
  mutate(
    Harper2024_SST_smooth_C = as.numeric(
      predict(Harper_SST_spline, x = Age_Ma)$y
    ),
    Harper2024_CO2_smooth_ppm = as.numeric(
      predict(Harper_CO2_spline, x = Age_Ma)$y
    )
  ) %>%
  arrange(desc(Age_Ma))


# --- 2. Kelson Tornillo / Big Bend carbonate D47 dataset ----


# Samples excluded only from stratigraphic and age-domain plots.
# They remain in the complete processed dataset.
Kelson_exclude_strat_age <- c(
  "BB-TF3-14-003",
  "BB-TF2-14-036",
  "PS3-Bk",
  "BB-TF2-14-002",
  "BB12-077",
  "BB-TF3-14-012nod2",
  "BB-TF2-14-030"
)

Kelson_Tornillo_raw <- read_csv(
  here("data", "raw", "kelson_tornillo_D47.csv"),
  show_col_types = FALSE
)

Kelson_Tornillo_D47 <- Kelson_Tornillo_raw %>%
  transmute(
    sample_id = Sample_ID,
    Age_Ma = `Age (Ma)`,
    epoch = Epoch,
    section = Section,
    petro_type = Type,
    
    d13Ccarb_vpdb = d13Ccarb,
    d13Ccarb_sd = `d13C SD`,
    d13Ccarb_se = d13C_SE,
    
    d18Ocarb_vpdb = d18Ocarb,
    d18Ocarb_vsmow = 1.03091 * d18Ocarb + 30.91,
    d18Ocarb_sd = `d18O SD`,
    d18Ocarb_se = d18O_SE,
    
    D47 = `∆47 a`,
    D47_sd = `∆47_SD`,
    D47_se = `∆47 SE`,
    n_replicates = `n (# of replicates)`,
    
    T47_C = `T(∆47) b`,
    T47_se_C = `T(∆47) SE`,
    T47_95CI_C = `T(∆47) 95 % CI`,
    
    d18Ow_vsmow = d18Owater,
    d18Ow_error = `d18Ow error`
  ) %>%
  mutate(
    petro_class = case_when(
      str_detect(str_to_lower(petro_type), "micrite") ~ "Micrite",
      str_detect(str_to_lower(petro_type), "spar") ~ "Spar",
      str_detect(str_to_lower(petro_type), "radial") ~ "Radial",
      TRUE ~ petro_type
    )
  ) %>%
  filter(!is.na(Age_Ma)) %>%
  arrange(desc(Age_Ma))

Kelson_Tornillo_D47_strat_age <- Kelson_Tornillo_D47 %>%
  filter(!sample_id %in% Kelson_exclude_strat_age)

Kelson_Tornillo_micrite_strat_age <- Kelson_Tornillo_D47_strat_age %>%
  filter(
    petro_class == "Micrite",
    !is.na(T47_C),
    !is.na(Age_Ma)
  )

# ---3. Fricke et al. (1998): biogenic apatite isotope records ------- 

# Geographic scope: Clarks Fork Basin and McCullough Peaks,
# northern Bighorn Basin.
#
# Important:
#   - published_age_Ma uses the original Fricke et al. age model.
#   - within_specimen_height_mm records position along a tooth or scale;
#     it is NOT stratigraphic height.
#   - Table 1 does not identify a section for each specimen.

Fricke1998_raw <- read_csv(
  here("data", "raw", "FrickeEtAl1998_Table1_isotope_data.csv"),
  show_col_types = FALSE
)

Fricke1998_temperature_change <- read_csv(
  here("data", "raw", "FrickeEtAl1998_Table2_temperature_change.csv"),
  show_col_types = FALSE
)

Fricke1998_samples <- Fricke1998_raw %>%
  mutate(
    taxon = factor(taxon, levels = c("Coryphodon", "Gar")),
    lmz = factor(
      lmz,
      levels = c(
        "Tiffanian", "CF-3", "CF-2", "Wa-0", "Wa-1",
        "Wa-3", "Wa-4", "Wa-5", "Wa-6"
      ),
      ordered = TRUE
    ),
    dataset = "Fricke et al. (1998)",
    proxy_type = if_else(
      !is.na(d18Ophosphate_VSMOW),
      "Biogenic phosphate d18O",
      NA_character_
    ),
    geographic_scope = "Northern Bighorn Basin",
    stratigraphic_section = NA_character_,
    section_note = paste(
      "Samples are from either the Clarks Fork Basin or McCullough Peaks,",
      "but Table 1 does not provide specimen-level section assignments."
    )
  )

# First calculate one mean per independent specimen. This prevents a densely
# sampled Coryphodon tooth from receiving more weight than another tooth.
Fricke1998_specimen_means <- Fricke1998_samples %>%
  group_by(taxon, specimen_id, lmz, published_age_Ma) %>%
  summarise(
    n_d18Op = sum(!is.na(d18Ophosphate_VSMOW)),
    mean_d18Op_VSMOW = if_else(
      n_d18Op > 0,
      mean(d18Ophosphate_VSMOW, na.rm = TRUE),
      NA_real_
    ),
    sd_d18Op_VSMOW = if_else(
      n_d18Op > 1,
      sd(d18Ophosphate_VSMOW, na.rm = TRUE),
      NA_real_
    ),
    
    n_d18Oc = sum(!is.na(d18Ocarb_VSMOW)),
    mean_d18Oc_VSMOW = if_else(
      n_d18Oc > 0,
      mean(d18Ocarb_VSMOW, na.rm = TRUE),
      NA_real_
    ),
    sd_d18Oc_VSMOW = if_else(
      n_d18Oc > 1,
      sd(d18Ocarb_VSMOW, na.rm = TRUE),
      NA_real_
    ),
    
    n_d13C = sum(!is.na(d13Ccarb_VPDB)),
    mean_d13Ccarb_VPDB = if_else(
      n_d13C > 0,
      mean(d13Ccarb_VPDB, na.rm = TRUE),
      NA_real_
    ),
    sd_d13Ccarb_VPDB = if_else(
      n_d13C > 1,
      sd(d13Ccarb_VPDB, na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  )

# Then summarize independent specimens within each taxon and land-mammal zone.
Fricke1998_zone_means <- Fricke1998_specimen_means %>%
  group_by(taxon, lmz) %>%
  summarise(
    published_age_Ma = mean(published_age_Ma, na.rm = TRUE),
    n_specimens_d18Op = sum(!is.na(mean_d18Op_VSMOW)),
    zone_mean_d18Op_VSMOW = if_else(
      n_specimens_d18Op > 0,
      mean(mean_d18Op_VSMOW, na.rm = TRUE),
      NA_real_
    ),
    zone_sd_d18Op_VSMOW = if_else(
      n_specimens_d18Op > 1,
      sd(mean_d18Op_VSMOW, na.rm = TRUE),
      NA_real_
    ),
    zone_se_d18Op_VSMOW = if_else(
      n_specimens_d18Op > 1,
      zone_sd_d18Op_VSMOW / sqrt(n_specimens_d18Op),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    d18Op_95ci_lower = if_else(
      n_specimens_d18Op > 1,
      zone_mean_d18Op_VSMOW -
        qt(0.975, df = n_specimens_d18Op - 1) * zone_se_d18Op_VSMOW,
      NA_real_
    ),
    d18Op_95ci_upper = if_else(
      n_specimens_d18Op > 1,
      zone_mean_d18Op_VSMOW +
        qt(0.975, df = n_specimens_d18Op - 1) * zone_se_d18Op_VSMOW,
      NA_real_
    )
  ) %>%
  filter(is.finite(zone_mean_d18Op_VSMOW))

# Convenient taxon-specific plotting tables.
Fricke1998_Coryphodon <- Fricke1998_zone_means %>%
  filter(taxon == "Coryphodon")

Fricke1998_gar <- Fricke1998_zone_means %>%
  filter(taxon == "Gar")

# --- 4. Wing et al. (2000): leaf-margin MAT estimates ----
# The published meter levels are section-specific and should not be treated as
# equivalent to the Polecat Bench composite stratigraphic-height scale.

Wing2000_MAT <- read_csv(
  here("data", "raw", "Wing_et_al_2000_LMA_MAT.csv"),
  show_col_types = FALSE
) %>%
  arrange(desc(published_age_model_2_Ma)) %>%
  mutate(
    dataset = "Wing et al. (2000)",
    proxy_type = "Leaf-margin MAT",
    geographic_scope = "Bighorn Basin"
  )

# --- 4. Fricke and Wing, 2004: d18Op/LMA temp and d18Owater estimates North America -------


# Load Fricke & Wing (2004) dataset
fw <- read_csv(
  here("data", "raw", "Fricke&Wing2004_MAT_d18Owater.csv")
)

# Plot MAT vs paleolatitude
ggplot(
  fw,
  aes(
    x = Paleolatitude,
    y = MAAT,
    shape = Interval == "PETM"
  )
) +
  geom_point(size = 3.5, stroke = 1) +
  scale_shape_manual(
    values = c(16, 17),
    labels = c("Background", "PETM"),
    name = NULL
  ) +
  labs(
    x = "Paleolatitude (°N)",
    y = "Mean annual temperature (°C)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top"
  )

# 5. Kelson et al. (2026) D17O of modern soil waters -----------

modern_soilwater <- read_csv(
  here(
    "data",
    "excel files",
    "jrkelson-CZ17O_soilwater-efc3bd3",
    "sw.csv"
  ),
  show_col_types = FALSE
) %>%
  filter(
    !is.na(dp18O),
    !is.na(D17O_pmg),
    
    # Follow the filtering used in the source manuscript:
    siteID.1 %in% c("MOJ", "JOR", "REY", "ESGR"),
    siteID.1 != "JOR" | siteID.2 == "CSAND"
  )


# --- Quick diagnostic plots ----

p_Wing2000_MAT <- ggplot(
  Wing2000_MAT,
  aes(x = published_age_model_2_Ma, y = MAT_C)
) +
  geom_errorbar(
    aes(
      ymin = MAT_C - MAT_error_C,
      ymax = MAT_C + MAT_error_C
    ),
    width = 0,
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  scale_x_reverse() +
  labs(
    x = "Published age (Ma)",
    y = "Mean annual temperature (°C)",
    title = "Wing et al. (2000) leaf-margin MAT"
  ) +
  theme_classic(base_size = 11)

p_Wing2000_MAT 

p_Fricke1998_d18Op <- ggplot(
  Fricke1998_zone_means,
  aes(
    x = published_age_Ma,
    y = zone_mean_d18Op_VSMOW,
    shape = taxon,
    group = taxon
  )
) +
  geom_errorbar(
    aes(
      ymin = d18Op_95ci_lower,
      ymax = d18Op_95ci_upper
    ),
    width = 0,
    linewidth = 0.35,
    alpha = 0.6,
    na.rm = TRUE
  ) +
  geom_line(linewidth = 0.5, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  scale_x_reverse() +
  labs(
    x = "Published age (Ma)",
    y = expression(delta^18 * O[phosphate] ~ "(‰ VSMOW)"),
    shape = "Taxon",
    title = "Fricke et al. (1998) biogenic phosphate"
  ) +
  theme_classic(base_size = 11)

p_Fricke1998_d18Op 

p_Kelson_T47_d18Ocarb <- Kelson_Tornillo_D47 %>%
  filter(!is.na(T47_C), !is.na(d18Ocarb_vsmow)) %>%
  ggplot(aes(x = T47_C, y = d18Ocarb_vsmow, shape = petro_class)) +
  geom_errorbarh(
    aes(
      xmin = T47_C - T47_se_C,
      xmax = T47_C + T47_se_C
    ),
    height = 0,
    linewidth = 0.35,
    alpha = 0.35,
    na.rm = TRUE
  ) +
  geom_errorbar(
    aes(
      ymin = d18Ocarb_vsmow - d18Ocarb_se,
      ymax = d18Ocarb_vsmow + d18Ocarb_se
    ),
    width = 0,
    linewidth = 0.35,
    alpha = 0.35,
    na.rm = TRUE
  ) +
  geom_point(size = 2.2, alpha = 0.8) +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = expression(delta^18 * O[carb] ~ "(‰ VSMOW)"),
    shape = "Petrography"
  ) +
  theme_classic(base_size = 11)

p_Kelson_T47_d18Ocarb

p_Kelson_micrite_T47_age <- ggplot(
  Kelson_Tornillo_micrite_strat_age,
  aes(x = T47_C, y = Age_Ma)
) +
  geom_errorbarh(
    aes(
      xmin = T47_C - T47_se_C,
      xmax = T47_C + T47_se_C
    ),
    height = 0,
    linewidth = 0.35,
    alpha = 0.45,
    na.rm = TRUE
  ) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_y_reverse() +
  scale_x_continuous(
    breaks = seq(10, 70, by = 10),
    limits = c(10, 70),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = "Age (Ma)",
    title = "Kelson Tornillo micrite D47 temperatures"
  ) +
  theme_classic(base_size = 11)

# Display quick checks when the script is run interactively.
if (interactive()) {
  print(p_Wing2000_MAT)
  print(p_Fricke1998_d18Op)
  print(p_Kelson_T47_d18Ocarb)
  print(p_Kelson_micrite_T47_age)
  
  print(Kelson_Tornillo_D47 %>% count(epoch, petro_class))
  print(Kelson_Tornillo_D47_strat_age %>% count(epoch, petro_class))
}

# -- Export processed datasets for downstream use ----
# All downstream scripts should load these processed files rather than
# re-reading or re-cleaning the original source data.

write_csv(
  Harper2024_CO2_SST,
  file.path(processed_dir, "Harper2024_CO2_SST_processed.csv")
)

write_csv(
  Kelson_Tornillo_D47,
  file.path(processed_dir, "Kelson_Tornillo_D47_processed.csv")
)

write_csv(
  Kelson_Tornillo_D47_strat_age,
  file.path(processed_dir, "Kelson_Tornillo_D47_strat_age_filtered.csv")
)

write_csv(
  Kelson_Tornillo_micrite_strat_age,
  file.path(processed_dir, "Kelson_Tornillo_micrite_strat_age_filtered.csv")
)

write_csv(
  Fricke1998_samples,
  file.path(processed_dir, "FrickeEtAl1998_samples_processed.csv")
)

write_csv(
  Fricke1998_specimen_means,
  file.path(processed_dir, "FrickeEtAl1998_specimen_means.csv")
)

write_csv(
  Fricke1998_zone_means,
  file.path(processed_dir, "FrickeEtAl1998_zone_means.csv")
)

write_csv(
  Fricke1998_temperature_change,
  file.path(processed_dir, "FrickeEtAl1998_temperature_change.csv")
)

write_csv(
  Wing2000_MAT,
  file.path(processed_dir, "WingEtAl2000_LMA_MAT_processed.csv")
)

# Optional confirmation when run interactively.
if (interactive()) {
  message("Processed reference datasets exported to: ", processed_dir)
}

