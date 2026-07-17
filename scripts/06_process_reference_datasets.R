# 06_process_reference_datasets.R
# Purpose:
#   Process external reference datasets for comparison with BHB records.
#
# Outputs:
#   - Harper2024_CO2_SST_processed.csv
#   - Kelson_Tornillo_D47_processed.csv
#   - Kelson_Tornillo_D47_strat_age_filtered.csv
#   - BHB_GTS2020_tiepoints.csv
#   - BHB_GTS2020_age_model.csv

# ---- Load packages ----
library(tidyverse)
library(here)
library(zoo)

# --- Output directories ----
dir.create(here("data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("figures", "reference_datasets"), recursive = TRUE, showWarnings = FALSE)

# ---- 1. Harper et al. (2024) CO2 and SST reference records -----------


Harper2024_CO2 <- read_csv(
  here("data", "raw", "HarperEtAl2024_co2_out.csv")
)

Harper2024_SST <- read_csv(
  here("data", "raw", "HarperEtAl2024_sst_out.csv")
)

Harper2024_CO2_SST <- Harper2024_CO2 %>%
  transmute(
    Age_Ma = age / 1000,
    Harper2024_mean_CO2_ppm = mean,
    Harper2024_CO2_lower95_ppm = `2.50%`,
    Harper2024_CO2_upper95_ppm = `97.50%`
  ) %>%
  left_join(
    Harper2024_SST %>%
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
      fill = NA,
      align = "center",
      partial = TRUE
    ),
    Harper2024_CO2_anchor_ppm = rollapply(
      Harper2024_mean_CO2_ppm,
      width = 3,
      FUN = mean,
      fill = NA,
      align = "center",
      partial = TRUE
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

# --- 2. Kelson Tornillo / Big Bend D47 reference dataset ----

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
  here("data", "raw", "kelson_tornillo_D47.csv")
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

# Filtered version for strat/age plots only
Kelson_Tornillo_D47_strat_age <- Kelson_Tornillo_D47 %>%
  filter(!sample_id %in% Kelson_exclude_strat_age)


# --- 3. Fricke et al. (1998) biogenic apatite isotope records ---------
# Northern Bighorn Basin: Clarks Fork Basin and McCullough Peaks
#
# Important:
# - published_age_Ma uses the original Fricke et al. age model.
# - within_specimen_height_mm is position along a tooth or scale,
#   NOT stratigraphic height.
# - Individual specimens cannot be assigned confidently to either
#   Clarks Fork Basin or McCullough Peaks from Table 1 alone.


library(tidyverse)
library(here)

# Load extracted sample-level isotope data ----

Fricke1998_raw <- read_csv(
  here(
    "data",
    "raw",
    "FrickeEtAl1998_Table1_isotope_data.csv"
  ),
  show_col_types = FALSE
)

# Load published temperature-change estimates ----

Fricke1998_temperature_change <- read_csv(
  here(
    "data",
    "raw",
    "FrickeEtAl1998_Table2_temperature_change.csv"
  ),
  show_col_types = FALSE
)

# Clean sample-level dataset ----

Fricke1998_samples <- Fricke1998_raw %>%
  mutate(
    taxon = factor(
      taxon,
      levels = c("Coryphodon", "Gar")
    ),
    
    lmz = factor(
      lmz,
      levels = c(
        "Tiffanian",
        "CF-3",
        "CF-2",
        "Wa-0",
        "Wa-1",
        "Wa-3",
        "Wa-4",
        "Wa-5",
        "Wa-6"
      ),
      ordered = TRUE
    ),
    
    dataset = "Fricke et al. (1998)",
    
    proxy_type = case_when(
      !is.na(d18Ophosphate_VSMOW) ~ "Biogenic phosphate d18O",
      TRUE ~ NA_character_
    ),
    
    geographic_scope = "Northern Bighorn Basin",
    
    stratigraphic_section = NA_character_,
    
    section_note = paste(
      "Samples are from either the Clarks Fork Basin",
      "or McCullough Peaks, but Table 1 does not provide",
      "a specimen-level section assignment."
    )
  )


# Specimen-level summaries
#
# Coryphodon teeth contain multiple serial samples. Calculate
# one mean per tooth before calculating land-mammal-zone means,
# so densely sampled teeth do not receive disproportionate weight.
#
# Gar scales are already individual sample observations.


Fricke1998_specimen_means <- Fricke1998_samples %>%
  group_by(
    taxon,
    specimen_id,
    lmz,
    published_age_Ma
  ) %>%
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

# Land-mammal-zone summaries
#
# These summaries treat individual teeth or scales as the
# independent observational units.

Fricke1998_zone_means <- Fricke1998_specimen_means %>%
  group_by(
    taxon,
    lmz
  ) %>%
  summarise(
    published_age_Ma = mean(
      published_age_Ma,
      na.rm = TRUE
    ),
    
    n_specimens_d18Op = sum(
      !is.na(mean_d18Op_VSMOW)
    ),
    
    mean_d18Op_VSMOW = mean(
      mean_d18Op_VSMOW,
      na.rm = TRUE
    ),
    
    sd_d18Op_VSMOW = if_else(
      n_specimens_d18Op > 1,
      sd(mean_d18Op_VSMOW, na.rm = TRUE),
      NA_real_
    ),
    
    se_d18Op_VSMOW = if_else(
      n_specimens_d18Op > 1,
      sd_d18Op_VSMOW / sqrt(n_specimens_d18Op),
      NA_real_
    ),
    
    d18Op_95ci_lower = if_else(
      n_specimens_d18Op > 1,
      mean_d18Op_VSMOW -
        qt(0.975, df = n_specimens_d18Op - 1) *
        se_d18Op_VSMOW,
      NA_real_
    ),
    
    d18Op_95ci_upper = if_else(
      n_specimens_d18Op > 1,
      mean_d18Op_VSMOW +
        qt(0.975, df = n_specimens_d18Op - 1) *
        se_d18Op_VSMOW,
      NA_real_
    ),
    
    .groups = "drop"
  ) %>%
  filter(
    is.finite(mean_d18Op_VSMOW)
  )

# ---- Optional: separate taxon-specific plotting tables ----

Fricke1998_Coryphodon <- Fricke1998_zone_means %>%
  filter(taxon == "Coryphodon")

Fricke1998_gar <- Fricke1998_zone_means %>%
  filter(taxon == "Gar")

# ---- Export processed reference tables ----

write_csv(
  Fricke1998_specimen_means,
  here(
    "data",
    "processed",
    "FrickeEtAl1998_specimen_means.csv"
  )
)

write_csv(
  Fricke1998_zone_means,
  here(
    "data",
    "processed",
    "FrickeEtAl1998_zone_means.csv"
  )
)

write_csv(
  Fricke1998_temperature_change,
  here(
    "data",
    "processed",
    "FrickeEtAl1998_temperature_change.csv"
  )
)


library(tidyverse)
library(here)

# ---- Load Wing et al. (2000) MAT ----
Wing_MAT <- read_csv(
  here("data", "raw", "Wing_et_al_2000_LMA_MAT.csv")
)


# ---- Quick checks ----

ggplot(
  Wing_MAT,
  aes(
    x = published_age_model_2_Ma,
    y = MAT_C
  )
) +
  geom_errorbar(
    aes(
      ymin = MAT_C - MAT_error_C,
      ymax = MAT_C + MAT_error_C
    ),
    width = 0
  ) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.5) +
  scale_x_reverse() +
  labs(
    x = "Published age (Ma)",
    y = "Mean annual temperature (°C)"
  ) +
  theme_classic()

ggplot(
  Fricke1998_zone_means,
  aes(
    x = published_age_Ma,
    y = mean_d18Op_VSMOW,
    color = taxon
  )
) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.6) +
  scale_x_reverse() +
  labs(
    x = "Published age (Ma)",
    y = expression(delta^18 * O[phosphate] ~ "(\u2030 VSMOW)")
  ) +
  theme_classic()


Kelson_Tornillo_D47 %>%
  count(epoch, petro_class)

Kelson_Tornillo_D47_strat_age %>%
  count(epoch, petro_class)

summary(Kelson_Tornillo_D47$T47_C)
summary(Kelson_Tornillo_D47$d18Ocarb_vpdb)
summary(Kelson_Tornillo_D47$d18Ocarb_vsmow)

# ---- Plot: Kelson D47 temperature vs d18Ocarb ----
# Full dataset retained here.

p_Kelson_T47_d18Ocarb <- Kelson_Tornillo_D47 %>%
  filter(!is.na(T47_C), !is.na(d18Ocarb_vsmow)) %>%
  ggplot(aes(
    x = T47_C,
    y = d18Ocarb_vsmow,
    shape = petro_class
  )) +
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
  geom_point(size = 2.2, alpha = 0.80) +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = expression(delta^18 * O[carb] ~ "(‰ VSMOW)"),
    shape = "Petrography"
  ) +
  theme_classic(base_size = 11)

p_Kelson_T47_d18Ocarb

# ---- Plot: Kelson micrite-only D47 temperature vs age ----
# Excluded samples removed here.

Kelson_Tornillo_micrite_strat_age <- Kelson_Tornillo_D47_strat_age %>%
  filter(
    petro_class == "Micrite",
    !is.na(T47_C),
    !is.na(Age_Ma)
  )

p_Kelson_micrite_T47_age <- ggplot(Kelson_Tornillo_micrite_strat_age) +
  geom_errorbarh(
    aes(
      xmin = T47_C - T47_se_C,
      xmax = T47_C + T47_se_C,
      y = Age_Ma
    ),
    height = 0,
    linewidth = 0.35,
    alpha = 0.45,
    na.rm = TRUE
  ) +
  geom_point(
    aes(x = T47_C, y = Age_Ma),
    size = 1.8,
    alpha = 0.75
  ) +
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

p_Kelson_micrite_T47_age

# ============================================================
# 3. GTS2020 BHB age-model update
# ============================================================

GTS2020_tiepoints <- tribble(
  ~tie_point, ~strat_height_m, ~Age_Ma, ~source,
  "C26n",     NA_real_,        59.237, "GTS2020",
  "C25r",     NA_real_,        58.959, "GTS2020",
  "C25n",     NA_real_,        57.656, "GTS2020",
  "C24r",     NA_real_,        57.101, "GTS2020",
  "C24n.3n",  NA_real_,        53.900, "GTS2020"
)

BHB_multiproxy_summary <- read_csv(
  here("data", "processed", "BHB_multiproxy_summary.csv")
)

BHB_GTS2020_tiepoints_clean <- GTS2020_tiepoints %>%
  filter(!is.na(strat_height_m), !is.na(Age_Ma)) %>%
  arrange(strat_height_m)

BHB_GTS2020_age_model <- BHB_multiproxy_summary %>%
  distinct(MLA_horizon_id, strat_height_m) %>%
  filter(!is.na(strat_height_m)) %>%
  arrange(strat_height_m) %>%
  mutate(
    Age_Ma_GTS2020 = approx(
      x = BHB_GTS2020_tiepoints_clean$strat_height_m,
      y = BHB_GTS2020_tiepoints_clean$Age_Ma,
      xout = strat_height_m,
      rule = 1
    )$y
  )

p_BHB_GTS2020_age_model <- ggplot() +
  geom_path(
    data = BHB_GTS2020_age_model,
    aes(x = Age_Ma_GTS2020, y = strat_height_m),
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  geom_point(
    data = BHB_GTS2020_tiepoints_clean,
    aes(x = Age_Ma, y = strat_height_m),
    size = 2
  ) +
  geom_text(
    data = BHB_GTS2020_tiepoints_clean,
    aes(x = Age_Ma, y = strat_height_m, label = tie_point),
    hjust = -0.05,
    size = 3
  ) +
  scale_x_reverse() +
  labs(
    x = "Age (Ma; GTS2020)",
    y = "Stratigraphic height (m)",
    title = "BHB GTS2020 age-depth model"
  ) +
  theme_classic(base_size = 11)

p_BHB_GTS2020_age_model

# ============================================================
# 4. Save processed outputs
# ============================================================

write_csv(
  Harper2024_CO2_SST,
  here("data", "processed", "Harper2024_CO2_SST_processed.csv")
)

write_csv(
  Kelson_Tornillo_D47,
  here("data", "processed", "Kelson_Tornillo_D47_processed.csv")
)

write_csv(
  Kelson_Tornillo_D47_strat_age,
  here("data", "processed", "Kelson_Tornillo_D47_strat_age_filtered.csv")
)

write_csv(
  GTS2020_tiepoints,
  here("data", "processed", "BHB_GTS2020_tiepoints.csv")
)

write_csv(
  BHB_GTS2020_age_model,
  here("data", "processed", "BHB_GTS2020_age_model.csv")
)

ggsave(
  here("figures", "reference_datasets", "Kelson_T47_vs_d18Ocarb.png"),
  p_Kelson_T47_d18Ocarb,
  width = 5,
  height = 4,
  dpi = 600
)

ggsave(
  here("figures", "reference_datasets", "Kelson_micrite_T47_vs_age.png"),
  p_Kelson_micrite_T47_age,
  width = 4,
  height = 6,
  dpi = 600
)

ggsave(
  here("figures", "reference_datasets", "BHB_GTS2020_age_model.png"),
  p_BHB_GTS2020_age_model,
  width = 5,
  height = 6,
  dpi = 600
)