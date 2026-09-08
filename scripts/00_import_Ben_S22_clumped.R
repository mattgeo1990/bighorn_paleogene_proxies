# Import Ben Passey's final MATLAB corrections for Nu Dog Session 22.
#
# Temperature work uses carbonate-residual-corrected D47 on I-CDES (column AN).
# Dual-clumped work uses pure-CDES D47 and D48 (columns Z and AC). The workbook
# Petersen temperatures are retained for audit only; preferred temperatures are
# recalculated with Anderson et al. (2021), with no additional acid correction.

library(tidyverse)
library(readxl)
library(here)
source(here("scripts", "helpers", "D47_temperature_calibrations.R"))

source_file <- here(
  "data", "excel files", "Nu_Dog_Clump_Session22_BHB.xlsx"
)
if (!file.exists(source_file)) stop("Missing Ben Session 22 workbook: ", source_file)

strip_matlab_quotes <- function(x) {
  x %>% str_remove("^'") %>% str_remove("'$") %>% str_replace_all("''", "'")
}

read_results_sheet <- function(sheet_name, session_id) {
  x <- read_excel(source_file, sheet = sheet_name)
  expected <- c(
    "'Type1'", "'Status'", "'SampleID'", "'AnalyID'", "'d13CPDB'",
    "'d18OPDBmin'", "'D47_90 CDES'", "'D48_90 CDES'",
    "'D47_90_iCDES carbcorr'", "'TD47iCDES carb corr'"
  )
  if (!all(expected %in% names(x))) {
    stop(sheet_name, " does not contain the expected final-result columns.")
  }

  tibble(
    correction_session = session_id,
    workbook_sheet = sheet_name,
    workbook_row = seq_len(nrow(x)) + 1L,
    matlab_order = x[[1]],
    Type_1 = strip_matlab_quotes(x[["'Type1'"]]),
    Type_2 = strip_matlab_quotes(x[["'Type2'"]]),
    User = strip_matlab_quotes(x[["'User'"]]),
    ben_status = str_to_lower(strip_matlab_quotes(x[["'Status'"]])),
    Sample_Name = strip_matlab_quotes(x[["'SampleID'"]]),
    analysis_id = strip_matlab_quotes(x[["'AnalyID'"]]),
    analysis_datetime_excel = x[["'Date_Time'"]],
    acid_temperature_C = x[["'AcidTemp'"]],
    d13Ccarb_VPDB = x[["'d13CPDB'"]],
    d18Ocarb_VPDB = x[["'d18OPDBmin'"]],
    D47_CDES = x[["'D47_90 CDES'"]],
    D48_CDES = x[["'D48_90 CDES'"]],
    D47_iCDES_carb_corr = x[["'D47_90_iCDES carbcorr'"]],
    T47_Petersen_workbook_C = x[["'TD47iCDES carb corr'"]]
  ) %>%
    mutate(
      IPLnum = as.integer(str_match(analysis_id, "IPL-CI-([0-9]+)")[, 2]),
      analysis_datetime = as.POSIXct(
        (analysis_datetime_excel - 25569) * 86400,
        origin = "1970-01-01", tz = "UTC"
      ),
      d18Ocarb_VSMOW = 1.03091 * d18Ocarb_VPDB + 30.91,
      T47_Anderson2021_C = anderson_2021_T47_C(D47_iCDES_carb_corr),
      correction_source = paste0(
        "Ben Passey Session ", session_id,
        " MATLAB; D47 I-CDES carbonate-residual correction"
      )
    )
}

ben_S22_all_results <- bind_rows(
  read_results_sheet("S22A Results", "22A"),
  read_results_sheet("S22B Results", "22B")
)

legacy_crosswalk <- read_csv(
  here("data", "raw", "IPL_D47_BHB_Pg_Summary_June2026.csv"),
  show_col_types = FALSE
) %>%
  select(IPLnum, MLA_sample_id, MLA_horizon_id, strat_height_m) %>%
  distinct(IPLnum, .keep_all = TRUE)

horizon_registry <- read_csv(
  here("data", "raw", "SandCoulee_Polecat_nodules.csv"),
  show_col_types = FALSE
)
names(horizon_registry)[1] <- "MLA_horizon_id"
horizon_registry <- horizon_registry %>%
  select(MLA_horizon_id, registry_strat_height_m = strat_height_m) %>%
  distinct()

normalize_bhb_sample <- function(x) {
  x %>%
    str_remove("^IPL-CI-[0-9]+\\s+") %>%
    str_replace(regex("^PK95-242", ignore_case = TRUE), "PK95-SC-242") %>%
    str_remove(regex("-0[1-9]$", ignore_case = TRUE))
}

ben_S22_BHB_samples <- ben_S22_all_results %>%
  filter(
    Type_1 == "Sample",
    str_detect(Sample_Name, regex("^(IPL-CI-[0-9]+ )?(PK95|PB-)",
                                  ignore_case = TRUE))
  ) %>%
  left_join(legacy_crosswalk, by = "IPLnum") %>%
  mutate(
    normalized_sample_name = normalize_bhb_sample(Sample_Name),
    MLA_sample_id = coalesce(MLA_sample_id, normalized_sample_name),
    MLA_horizon_id = coalesce(
      MLA_horizon_id,
      normalized_sample_name %>%
        str_remove(regex("-(SPAR|clky|drk|lte)$", ignore_case = TRUE)) %>%
        str_replace(regex("^PK95-SC-80-n-[12]$", ignore_case = TRUE),
                    "PK95-SC-80")
    ),
    material_type = if_else(
      str_detect(MLA_sample_id, regex("-SPAR$", ignore_case = TRUE)),
      "secondary_fill_spar", "host_matrix"
    ),
    passes_analytical_screen = ben_status == "include" &
      is.finite(D47_iCDES_carb_corr) & is.finite(T47_Anderson2021_C),
    passes_fabric_screen = material_type == "host_matrix",
    passes_primary_temperature_screen =
      passes_analytical_screen & passes_fabric_screen,
    primary_screen_reason = case_when(
      ben_status != "include" ~ paste("Ben status:", ben_status),
      !is.finite(D47_iCDES_carb_corr) ~ "Missing final D47",
      !is.finite(T47_Anderson2021_C) ~ "Invalid Anderson temperature",
      material_type == "secondary_fill_spar" ~
        "Secondary fill spar; retained for paragenetic comparison",
      TRUE ~ "Accepted host-matrix analysis"
    )
  ) %>%
  left_join(horizon_registry, by = "MLA_horizon_id") %>%
  mutate(strat_height_m = coalesce(strat_height_m, registry_strat_height_m)) %>%
  select(-registry_strat_height_m)

if (any(is.na(ben_S22_BHB_samples$IPLnum))) {
  stop("Could not recover IPL number for one or more Bighorn analyses.")
}

write_csv(
  ben_S22_all_results,
  here("data", "processed", "Ben_S22_all_results_audit.csv")
)
write_csv(
  ben_S22_BHB_samples,
  here("data", "processed", "Ben_S22_BHB_run_level.csv")
)

message(
  "Imported ", nrow(ben_S22_BHB_samples), " Bighorn sample analyses: ",
  sum(ben_S22_BHB_samples$passes_primary_temperature_screen),
  " pass the primary host-matrix temperature screen."
)
