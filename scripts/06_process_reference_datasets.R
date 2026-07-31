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
#   5. Barnet et al. (2019): ODP Site 1262 benthic stable isotopes
#   6. GTS2020: updated BHB age-model tie points
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

# Prevent unsaved exploratory ggplot expressions from creating Rplots.pdf
# during non-interactive pipeline runs. Explicit ggsave() calls are unaffected.
opened_null_graphics_device <- FALSE
if (!interactive()) {
  grDevices::pdf(file = NULL)
  opened_null_graphics_device <- TRUE
}

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

# - 1b. Barnet et al. (2019): South Atlantic benthic isotopes -----
#
# This is the complete, orbitally tuned ODP Site 1262 record archived by
# Barnet et al. in PANGAEA (dataset DOI: 10.1594/PANGAEA.884585; article DOI:
# 10.1029/2019PA003556). The original PANGAEA text file is retained verbatim
# under data/raw so that the repository preserves its citation, metadata,
# license, source-study labels, and the exact values used here.
#
# The PANGAEA file contains a 29-line metadata block followed by a tabular
# header. Ages are supplied in ka BP and are converted to Ma. `d18O_corrected`
# is the published value corrected to Cibicidoides; `BWT_C` is the temperature
# estimate reported by Barnet et al., rather than a temperature recalculated
# in this project.

Barnet2019_raw <- read_tsv(
  here("data", "raw", "BarnetEtAl2019_PANGAEA_884585.tab"),
  skip = 29,
  name_repair = "minimal",
  show_col_types = FALSE
)

if (ncol(Barnet2019_raw) != 11) {
  stop(
    "Barnet et al. (2019) source does not have the expected 11 columns; ",
    "check the archived PANGAEA file before proceeding."
  )
}

# Assign portable ASCII names by verified PANGAEA column order. This avoids
# locale-dependent corruption of the source's delta, per-mil, and degree
# symbols when the pipeline is run under a C locale.
names(Barnet2019_raw) <- c(
  "event", "sample_label", "depth_top_m", "d13C_benthic_vpdb",
  "d18O_benthic_vpdb", "depth_mbsf", "depth_composite_mcd", "age_ka_bp",
  "source_record", "d18O_corrected_vpdb", "bottom_water_temperature_C"
)

Barnet2019_ODP1262 <- Barnet2019_raw %>%
  transmute(
    site = "ODP Site 1262",
    event,
    sample_label,
    depth_mbsf,
    depth_composite_mcd,
    Age_Ma = age_ka_bp / 1000,
    d13C_benthic_vpdb,
    d18O_benthic_vpdb,
    d18O_corrected_vpdb,
    bottom_water_temperature_C,
    source_record,
    dataset_doi = "10.1594/PANGAEA.884585",
    article_doi = "10.1029/2019PA003556"
  ) %>%
  filter(
    !is.na(Age_Ma),
    between(Age_Ma, 52, 68)
  ) %>%
  arrange(Age_Ma)

if (
  nrow(Barnet2019_ODP1262) == 0 ||
    !all(c(52.3, 67.1) >= range(Barnet2019_ODP1262$Age_Ma) - 0.1) ||
    !all(c(52.3, 67.1) <= range(Barnet2019_ODP1262$Age_Ma) + 0.1)
) {
  stop("Barnet et al. (2019) import failed its expected age-range check.")
}

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


# Samples excluded from the published paleoclimate reconstruction and/or the
# project's Paleocene-Eocene age-domain comparison remain in the complete
# inventory. The reasons are assigned explicitly below. Kelson et al. (2018)
# identify two intrusion-reset nodules, PS3 mixed material, and the radial
# black-paleosol nodule as problematic. The two Cretaceous samples are outside
# this project's target interval. BB-TF3-14-012nod2 was already excluded in the
# legacy project code; it is a paired nodule at the same age/horizon, but the
# original project rationale was not documented, so it is not called altered.

Kelson_Tornillo_raw <- read_csv(
  here("data", "raw", "kelson_tornillo_D47.csv"),
  show_col_types = FALSE
) %>%
  # Normalize Unicode source headers immediately so the pipeline is stable
  # across operating-system locale settings.
  rename(
    D47_source = 12,
    D47_sd_source = 13,
    D47_se_source = 14,
    T47_source_C = 16,
    T47_se_source_C = 17,
    T47_95CI_source_C = 18
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
    
    D47 = D47_source,
    D47_sd = D47_sd_source,
    D47_se = D47_se_source,
    n_replicates = `n (# of replicates)`,
    
    T47_C = T47_source_C,
    T47_se_C = T47_se_source_C,
    T47_95CI_C = T47_95CI_source_C,
    
    d18Ow_vsmow = d18Owater,
    d18Ow_error = `d18Ow error`
  ) %>%
  mutate(
    petro_class = case_when(
      str_detect(str_to_lower(petro_type), "micrite") ~ "Micrite",
      str_detect(str_to_lower(petro_type), "spar") ~ "Spar",
      str_detect(str_to_lower(petro_type), "radial") ~ "Radial",
      TRUE ~ petro_type
    ),
    kelson_exclusion_reason = case_when(
      sample_id %in% c("BB-TF3-14-003", "BB-TF2-14-036") ~
        paste(
          "Thermally reset near local igneous intrusion;",
          "excluded from paleoclimate reconstruction by Kelson et al. (2018)"
        ),
      sample_id == "PS3-Bk" ~
        paste(
          "Bulk nodule with probable micrite-spar mixing and T47 > 55 C;",
          "excluded by Kelson et al. (2018)"
        ),
      sample_id == "BB-TF2-14-030" ~
        paste(
          "Radial carbonate from black paleosol;",
          "excluded from paleoclimate reconstruction by Kelson et al. (2018)"
        ),
      sample_id %in% c("BB-TF2-14-002", "BB12-077") ~
        "Cretaceous; outside project Paleocene-Eocene age-domain target",
      sample_id == "BB-TF3-14-012nod2" ~
        paste(
          "Legacy project exclusion of paired nodule at the same 53.9 Ma",
          "horizon; original rationale not documented"
        ),
      petro_class != "Micrite" ~
        "Non-micrite material; excluded from primary temperature model",
      TRUE ~ NA_character_
    ),
    used_in_tornillo_temperature_model =
      petro_class == "Micrite" & is.na(kelson_exclusion_reason),
    kelson_screening_basis = if_else(
      used_in_tornillo_temperature_model,
      paste(
        "Primary micrite retained after Kelson et al. (2018) and legacy",
        "project sample exclusions"
      ),
      kelson_exclusion_reason
    )
  ) %>%
  filter(!is.na(Age_Ma)) %>%
  arrange(desc(Age_Ma))

Kelson_Tornillo_D47_strat_age <- Kelson_Tornillo_D47 %>%
  filter(used_in_tornillo_temperature_model)

Kelson_Tornillo_micrite_strat_age <- Kelson_Tornillo_D47_strat_age %>%
  filter(
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

# --- 5. Fricke and Wing (2004): North American temperature estimates ----


# Load Fricke & Wing (2004) dataset
fw <- read_csv(
  here("data", "raw", "Fricke&Wing2004_MAT_d18Owater.csv"),
  show_col_types = FALSE
)

# The source table provides two Bighorn Basin estimates tied to mammalian
# biozones rather than sample-specific numerical ages. Preserve that coarse
# temporal resolution instead of treating the values as precisely dated.
# Wa-0 is bracketed by the project PETM interval; Wa-6 is bracketed using the
# taxon-specific ages carried by the Fricke et al. (1998) reference table.
FrickeWing2004_BHB_MAAT <- fw %>%
  filter(str_detect(Locality, fixed("Bighorn Basin"))) %>%
  transmute(
    dataset = "Fricke and Wing (2004)",
    proxy_type = "Isotope/LMA MAAT estimate",
    geographic_scope = "Bighorn Basin",
    interval = Interval,
    biozone = NALMA,
    paleolatitude_deg_n = Paleolatitude,
    temperature_C = parse_number(as.character(MAAT)),
    Age_Ma = case_when(
      biozone == "Wa-0" ~ 55.90,
      biozone == "Wa-6" ~ mean(c(52.83, 53.01)),
      TRUE ~ NA_real_
    ),
    age_younger_ma = case_when(
      biozone == "Wa-0" ~ 55.75,
      biozone == "Wa-6" ~ 52.83,
      TRUE ~ NA_real_
    ),
    age_older_ma = case_when(
      biozone == "Wa-0" ~ 55.93,
      biozone == "Wa-6" ~ 53.01,
      TRUE ~ NA_real_
    ),
    age_basis = case_when(
      biozone == "Wa-0" ~
        "Project Wa-0/PETM datum; interval shown as 55.93-55.75 Ma",
      biozone == "Wa-6" ~
        "Range of Fricke et al. (1998) taxon-specific Wa-6 ages",
      TRUE ~ "No numerical age assignment"
    )
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

#-- 6.) ZB20a Astronomical Forcing at 47 Degrees North ---------------------
#
# The compact raw input is a 1-kyr extract of the pre-computed ZB20a(1,1)
# astronomical solution distributed with Kocken and Zeebe (2026):
# https://doi.org/10.1029/2025PA005287
# https://github.com/japhir/paleoinsolation
#
# Daily-mean top-of-atmosphere insolation is calculated at 47 degrees north
# for Northern Hemisphere summer solstice (true solar longitude = 90 degrees)
# using S0 = 1360.7 W m-2. The equations follow the authors' R implementation.
# Ages are the solution's astronomical ages relative to J2000 and are plotted
# on the same Ma axis as GTS2020-positioned proxy ages. No phase shift or tuning
# to the BHB data is applied. ZB20a eccentricity is geologically supported over
# this interval, but exact precession/obliquity phase becomes less secure before
# approximately 58 Ma; the processed table retains that interpretation flag.

ZB20a_orbital <- read_csv(
  here("data", "raw", "ZB20a_1-1_Thanetian_Ypresian_1kyr.csv"),
  show_col_types = FALSE
) %>%
  filter(
    if_all(
      c(eccentricity, obliquity_rad, longitude_perihelion_rad),
      ~ !is.na(.x)
    )
  )

daily_insolation <- function(
    eccentricity,
    obliquity_rad,
    longitude_perihelion_rad,
    latitude_rad,
    solar_longitude_rad,
    solar_constant_w_m2 = 1360.7) {
  true_anomaly <- solar_longitude_rad - longitude_perihelion_rad
  earth_sun_distance <-
    (1 - eccentricity^2) /
    (1 + eccentricity * cos(true_anomaly))
  sin_declination <- sin(obliquity_rad) * sin(solar_longitude_rad)
  cos_declination <- sqrt(1 - sin_declination^2)
  sin_lat_sin_dec <- sin(latitude_rad) * sin_declination
  cos_lat_cos_dec <- cos(latitude_rad) * cos_declination
  cos_hour_angle <- pmin(
    pmax(-1, -sin_lat_sin_dec / cos_lat_cos_dec),
    1
  )
  hour_angle <- acos(cos_hour_angle)
  sin_hour_angle <- sqrt(1 - cos_hour_angle^2)

  solar_constant_w_m2 / (pi * earth_sun_distance^2) *
    (hour_angle * sin_lat_sin_dec +
       cos_lat_cos_dec * sin_hour_angle)
}

BHB_insolation_47N <- ZB20a_orbital %>%
  transmute(
    Age_Ma = age_ma,
    astronomical_time_kyr_j2000 = time_kyr_j2000,
    eccentricity,
    obliquity_deg = obliquity_rad * 180 / pi,
    climatic_precession,
    summer_solstice_insolation_w_m2 = daily_insolation(
      eccentricity = eccentricity,
      obliquity_rad = obliquity_rad,
      longitude_perihelion_rad = longitude_perihelion_rad,
      latitude_rad = 47 * pi / 180,
      solar_longitude_rad = pi / 2,
      solar_constant_w_m2 = 1360.7
    ),
    paleolatitude_deg_n = 47,
    solar_longitude_deg = 90,
    solar_constant_w_m2 = 1360.7,
    astronomical_solution = "ZB20a(1,1)",
    chronology_alignment =
      "Native astronomical age plotted on GTS2020 Ma axis; no phase shift",
    phase_interpretation = if_else(
      Age_Ma > 58,
      paste(
        "Precession/obliquity phase less secure;",
        "interpret exact insolation peaks cautiously"
      ),
      "Geologically supported ZB20a interval"
    )
  ) %>%
  arrange(desc(Age_Ma))

#-- 7.) Kelson et al. (2026) D17O of Modern Soil Waters --------------------

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

#-- 8.) Assemble Non-CFB Soil-Carbonate Reference Records ------------------

# Script 01 is the authoritative cleaning location for Snell and Koch. This
# script consumes those processed summaries rather than repeating their source
# cleaning. Only non-CFB observations enter the regional reference table.
Snell_soilcarb_processed <- read_csv(
  here("data", "processed", "SnellEtAl2013_soilcarb_summary.csv"),
  show_col_types = FALSE
)

Koch_soilcarb_processed <- read_csv(
  here("data", "processed", "Koch_soilcarb_summary.csv"),
  show_col_types = FALSE
)

BHB_regional_soilcarb_reference_summary <- bind_rows(
  Snell_soilcarb_processed %>%
    filter(section_id != "CFB") %>%
    transmute(
      dataset = "Snell et al. (2013)",
      section_id,
      MLA_horizon_id,
      strat_height_m,
      published_age_ma = snell2013_age_ma,
      d13Ccarb_vpdb = Snell_mean_d13Ccarb_vpdb,
      d18Ocarb_vsmow = Snell_mean_d18Ocarb_vsmow,
      T47_C = Snell_mean_T47_C,
      T47_se_C = Snell_se_T47_C,
      published_d18Owater_vsmow = Snell_mean_d18Ow_vsmow
    ),
  Koch_soilcarb_processed %>%
    filter(section_id != "CFB") %>%
    transmute(
      dataset = "Koch et al.",
      section_id,
      MLA_horizon_id,
      strat_height_m,
      published_age_ma = NA_real_,
      d13Ccarb_vpdb = Koch_mean_d13Ccarb_vpdb,
      d18Ocarb_vsmow = Koch_mean_d18Ocarb_vsmow,
      T47_C = NA_real_,
      T47_se_C = NA_real_,
      published_d18Owater_vsmow = NA_real_
    )
) %>%
  arrange(section_id, strat_height_m)

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

write_csv(
  FrickeWing2004_BHB_MAAT,
  file.path(processed_dir, "FrickeWing2004_BHB_MAAT_processed.csv")
)

write_csv(
  BHB_insolation_47N,
  file.path(processed_dir, "BHB_ZB20a_summer_insolation_47N.csv")
)

write_csv(
  BHB_regional_soilcarb_reference_summary,
  file.path(processed_dir, "BHB_regional_soilcarb_reference_summary.csv")
)

write_csv(
  Barnet2019_ODP1262,
  file.path(processed_dir, "BarnetEtAl2019_ODP1262_benthic_isotopes.csv")
)

# Optional confirmation when run interactively.
if (interactive()) {
  message("Processed reference datasets exported to: ", processed_dir)
}

if (opened_null_graphics_device) {
  grDevices::dev.off()
}
