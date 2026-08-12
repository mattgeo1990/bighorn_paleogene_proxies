# 00_correct_IPL_triple_oxygen_reactors.R
#
# Purpose
# -------
# Recalculate IPL triple-oxygen results from the six reactor workbooks,
# convert O2 measurements to carbonate-mineral values following the manual
# BP-validated workflow, compile project samples in the historical standardized
# column format, and reconcile the new calculation against legacy IPL/manual
# results.
#
# Source workbooks are read only. Outputs are written to data/processed unless
# IPL17O_OUTPUT_DIR is set. The source directory can be overridden with
# IPL17O_REACTOR_DIR.

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(readr)
  library(readxl)
  library(stringr)
  library(tidyr)
})

# ---- Configuration -------------------------------------------------------

reactor_dir <- Sys.getenv(
  "IPL17O_REACTOR_DIR",
  unset = paste0(
    "/Users/allen/Library/CloudStorage/OneDrive-Personal/MLA Work/",
    "NSF_EARPF/Data/D17O and D47 results/all_IPL_data_7-19-26"
  )
)

output_dir <- Sys.getenv(
  "IPL17O_OUTPUT_DIR",
  unset = here("data", "processed")
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

reactor_files <- tibble::tribble(
  ~reactor_id, ~file_name,
  31L, "Cap17O Compiled REACTOR THIRTY ONE.xlsx",
  32L, "Cap17O Compiled REACTOR THIRTY TWO.xlsx",
  33L, "Cap17O Compiled REACTOR THIRTY THREE.xlsx",
  34L, "Cap17O Compiled REACTOR THIRTY FOUR.xlsx",
  35L, "Cap17O Compiled REACTOR THIRTY FIVE.xlsx",
  36L, "Cap17O Compiled REACTOR THIRTY SIX_7-19-26.xlsx"
) %>%
  mutate(
    path = file.path(reactor_dir, file_name),
    raw_path = if_else(
      reactor_id == 33L,
      paste0(
        "/Users/allen/Library/CloudStorage/OneDrive-Personal/MLA Work/",
        "NSF_EARPF/Data/Data Corrections/R33_DataCorrection.xlsx"
      ),
      path
    ),
    raw_skip = if_else(reactor_id == 33L, 2L, 0L)
  )

missing_workbooks <- unique(c(
  reactor_files$path[!file.exists(reactor_files$path)],
  reactor_files$raw_path[!file.exists(reactor_files$raw_path)]
))
if (length(missing_workbooks) > 0) {
  stop(
    "Missing reactor workbook(s):\n",
    paste0("  - ", missing_workbooks, collapse = "\n")
  )
}

# Historical project tables provide sample/horizon identifiers and independent
# comparison values. They are never used to calculate the new BP correction.
legacy_path <- here(
  "data", "raw", "all_data_PgBHB_IPL17O_standardized_columns.csv"
)
comparison_path <- here(
  "data", "raw", "all_data_PgBHB_IPL17O_CODEX_AUTOMATED_CORRECTIONS.csv"
)

legacy <- read_csv(legacy_path, show_col_types = FALSE) %>%
  rename_with(~ str_remove(.x, "^\\ufeff")) %>%
  mutate(IPL_num = as.integer(IPL_num))

comparison <- read_csv(comparison_path, show_col_types = FALSE) %>%
  mutate(IPL_num = as.integer(IPL_num))

# ---- Constants from the BP-validated manual workflow --------------------

lambda_D17O <- 0.528

# Huth et al. (2022) conversion recorded in R34_DataCorrection_bp_validated.
alpha18_mineral_O2 <- 0.9918723
theta_mineral_O2 <- 0.5224019071720026
alpha17_mineral_O2 <- alpha18_mineral_O2 ^ theta_mineral_O2

standard_assignments <- tibble::tribble(
  ~standard_key, ~known_dp18O_gas, ~accepted_Dp17O_mineral_permeg,
  "102GCAZ01", 23.91, -67,
  "IAEAC1",     36.30, -100
)

normalize_standard_name <- function(x) {
  str_to_upper(x) %>% str_replace_all("[^A-Z0-9]", "")
}

num <- function(x) suppressWarnings(as.numeric(x))

cell_num <- function(tbl, row, col) {
  num(tbl[[col]][row])
}

# ---- Read and independently reproduce the workbook O2 reduction ---------

read_reactor <- function(reactor_id, path, raw_path, raw_skip) {
  raw <- suppressWarnings(read_excel(
    raw_path,
    sheet = "All Data",
    skip = raw_skip,
    .name_repair = "unique_quiet"
  ))
  names(raw) <- names(raw) %>%
    str_replace(fixed("d'17O Final, as O2"), "d'17O Final") %>%
    str_replace(fixed("d'18O Final, as O2"), "d'18O Final") %>%
    str_replace(fixed("D17O per meg, as O2"), "D17O per meg")
  smow <- suppressWarnings(read_excel(
    path,
    sheet = "SMOW",
    col_names = FALSE,
    .name_repair = "minimal"
  ))

  # Excel coordinates used by all six compiled reactor workbooks:
  # Z4/AA4: measured SMOW mass-33/34 baseline
  # AN6/AN12: zero-intercept SLAP transfer slopes
  # AN14: lambda used by the workbook (normally 0.528)
  smow_d33 <- cell_num(smow, 4, 26)
  smow_d34 <- cell_num(smow, 4, 27)
  slope17 <- cell_num(smow, 6, 40)
  slope18 <- cell_num(smow, 12, 40)
  workbook_lambda <- cell_num(smow, 14, 40)

  if (any(!is.finite(c(smow_d33, smow_d34, slope17, slope18)))) {
    stop("Missing SMOW-SLAP calibration constants in Reactor ", reactor_id)
  }
  if (!is.finite(workbook_lambda)) workbook_lambda <- lambda_D17O

  out <- raw %>%
    transmute(
      reactor_id = reactor_id,
      IPL_num = as.integer(num(.data[["IPL num"]])),
      Name_full = as.character(.data[["NAME"]]),
      comments = as.character(.data[["Comments"]]),
      type_1 = as.character(.data[["Type 1"]]),
      type_2 = as.character(.data[["Type 2"]]),
      d17O = num(.data[["d17O"]]),
      d.17O = num(.data[["d'17O"]]),
      d17O.err = num(.data[["d17O err"]]),
      d18O = num(.data[["d18O"]]),
      d.18O = num(.data[["d'18O"]]),
      d18O.err = num(.data[["d18O err"]]),
      CAP.17O = num(.data[["CAP 17O"]]),
      CAP17O.err = num(.data[["CAP17O err"]]),
      d33 = num(.data[["d33"]]),
      d33.err = num(.data[["d33 err"]]),
      d34 = num(.data[["d34"]]),
      d34.err = num(.data[["d34 err"]]),
      d35 = num(.data[["d35"]]),
      d35.err = num(.data[["d35 err"]]),
      d36 = num(.data[["d36"]]),
      d36.err = num(.data[["d36 err"]]),
      `33_mismatch` = num(.data[["33 mismatch R2"]]),
      `34_mismatch` = num(.data[["34 mismatch R2"]]),
      flag_major = num(.data[["flag.major"]]),
      flag_analysis = num(.data[["flag.analysis"]]),
      workbook_dp17O_O2 = num(.data[["d'17O Final"]]),
      workbook_dp18O_O2 = num(.data[["d'18O Final"]]),
      workbook_Dp17O_O2_permeg = num(.data[["D17O per meg"]])
    ) %>%
    filter(!is.na(IPL_num), is.finite(d33), is.finite(d34)) %>%
    mutate(
      d17O_SMOW = 1000 * (((1 + d33 / 1000) /
        (1 + smow_d33 / 1000)) - 1),
      d18O_SMOW = 1000 * (((1 + d34 / 1000) /
        (1 + smow_d34 / 1000)) - 1),
      d17O_SMOWSLAP = d17O_SMOW * slope17,
      d18O_SMOWSLAP = d18O_SMOW * slope18,
      dp17O_O2_SMOWSLAP = 1000 * log1p(d17O_SMOWSLAP / 1000),
      dp18O_O2_SMOWSLAP = 1000 * log1p(d18O_SMOWSLAP / 1000),
      Dp17O_O2_SMOWSLAP2 = 1000 * (
        dp17O_O2_SMOWSLAP - workbook_lambda * dp18O_O2_SMOWSLAP
      ),
      standard_key = normalize_standard_name(type_2),
      calibration_eligible =
        standard_key %in% standard_assignments$standard_key &
        coalesce(flag_major, 0) == 0 &
        coalesce(flag_analysis, 0) == 0 &
        is.finite(dp18O_O2_SMOWSLAP)
    )

  calibration <- out %>%
    filter(calibration_eligible) %>%
    inner_join(standard_assignments, by = "standard_key")

  n_standard_types <- n_distinct(calibration$standard_key)
  if (nrow(calibration) < 2) {
    stop("Reactor ", reactor_id, " lacks enough accepted carbonate standards.")
  }

  # Ben Passey's validated procedure regresses known gas delta-prime-18O on
  # measured gas delta-prime-18O. Equal standard-type weighting prevents a
  # reactor with many IAEA-C1 analyses from overwhelming 102-GC-AZ01.
  calibration_balanced <- calibration %>%
    group_by(standard_key) %>%
    mutate(regression_weight = 1 / n()) %>%
    ungroup()

  if (n_standard_types >= 2) {
    d18_fit <- lm(
      known_dp18O_gas ~ dp18O_O2_SMOWSLAP,
      data = calibration_balanced,
      weights = regression_weight
    )
    d18_intercept <- unname(coef(d18_fit)[1])
    d18_slope <- unname(coef(d18_fit)[2])
    d18_method <- "two-standard linear transfer"
  } else {
    # Reactor 31 contains IAEA-C1 but not 102-GC-AZ01. A slope cannot be
    # estimated from one material, so retain unit slope and correct its mean
    # to the accepted gas value. This limitation is explicit in the output.
    d18_slope <- 1
    d18_intercept <- mean(
      calibration_balanced$known_dp18O_gas -
        calibration_balanced$dp18O_O2_SMOWSLAP
    )
    d18_method <- "one-standard intercept-only transfer"
  }

  convert_to_mineral <- function(df) {
    df %>%
      mutate(
        dp18O_O2_corrected =
          d18_intercept + d18_slope * dp18O_O2_SMOWSLAP,
        d18O_correction = dp18O_O2_corrected - dp18O_O2_SMOWSLAP,
        # Preserve the measured gas Delta-prime-17O while correcting delta-18O.
        dp17O_O2_corrected =
          dp17O_O2_SMOWSLAP + lambda_D17O * d18O_correction,
        d17O_O2_corrected = 1000 * expm1(dp17O_O2_corrected / 1000),
        d18O_O2_corrected = 1000 * expm1(dp18O_O2_corrected / 1000),
        d17Ocarb = (1000 + d17O_O2_corrected) * alpha17_mineral_O2 - 1000,
        d18Ocarb = (1000 + d18O_O2_corrected) * alpha18_mineral_O2 - 1000,
        dp17Ocarb_SMOWSLAP = 1000 * log1p(d17Ocarb / 1000),
        dp18Ocarb_SMOWSLAP = 1000 * log1p(d18Ocarb / 1000),
        Dp17Ocarb_pre_standard = 1000 * (
          dp17Ocarb_SMOWSLAP - lambda_D17O * dp18Ocarb_SMOWSLAP
        )
      )
  }

  calibration_mineral <- convert_to_mineral(calibration_balanced) %>%
    mutate(
      standard_residual =
        accepted_Dp17O_mineral_permeg - Dp17Ocarb_pre_standard
    )

  # Give the two standard materials equal influence, matching their equal role
  # in the d18O regression and avoiding dependence on unequal replicate counts.
  standard_offset_table <- calibration_mineral %>%
    group_by(standard_key) %>%
    summarise(
      n = n(),
      mean_measured = mean(Dp17Ocarb_pre_standard),
      accepted = first(accepted_Dp17O_mineral_permeg),
      offset = accepted - mean_measured,
      sd = sd(Dp17Ocarb_pre_standard),
      .groups = "drop"
    )

  standard_offset <- mean(standard_offset_table$offset)

  corrected <- convert_to_mineral(out) %>%
    mutate(
      Dp17Ocarb_permeg_final_correction =
        Dp17Ocarb_pre_standard + standard_offset,
      calibration_d18_slope = d18_slope,
      calibration_d18_intercept = d18_intercept,
      calibration_standard_offset = standard_offset,
      calibration_n = nrow(calibration),
      calibration_method = paste(
        "BP validated; equal standard-type weighting;", d18_method
      )
    )

  calibration_summary <- tibble(
    reactor_id = reactor_id,
    smow_d33 = smow_d33,
    smow_d34 = smow_d34,
    slap_slope17 = slope17,
    slap_slope18 = slope18,
    lambda = workbook_lambda,
    d18_slope = d18_slope,
    d18_intercept = d18_intercept,
    d18_method = d18_method,
    standard_offset = standard_offset,
    n_calibration = nrow(calibration),
    n_102 = sum(calibration$standard_key == "102GCAZ01"),
    n_IAEA_C1 = sum(calibration$standard_key == "IAEAC1")
  )

  list(
    corrected = corrected,
    calibration = calibration_summary,
    standard_offsets = standard_offset_table %>% mutate(reactor_id = reactor_id)
  )
}

reactor_results <- Map(
  read_reactor,
  reactor_files$reactor_id,
  reactor_files$path,
  reactor_files$raw_path,
  reactor_files$raw_skip
)

all_corrected <- bind_rows(lapply(reactor_results, `[[`, "corrected"))
calibration_summary <- bind_rows(lapply(reactor_results, `[[`, "calibration"))
standard_offset_summary <- bind_rows(
  lapply(reactor_results, `[[`, "standard_offsets")
)

# ---- Compile the project analyses in the historical standardized format --

id_lookup <- legacy %>%
  select(IPL_num, MLA_sample_id, MLA_horizon_id) %>%
  distinct(IPL_num, .keep_all = TRUE)

manual_lookup <- comparison %>%
  select(
    IPL_num,
    Dp17Ocarb_permeg_MANUAL_correction,
    previous_automated = Dp17Ocarb_permeg_final_correction,
    previous_accepted = Dp17Ocarb_permeg_ACCEPTED
  )

compiled <- all_corrected %>%
  inner_join(id_lookup, by = "IPL_num") %>%
  left_join(manual_lookup, by = "IPL_num") %>%
  mutate(
    manual_automated_offset =
      Dp17Ocarb_permeg_MANUAL_correction -
        Dp17Ocarb_permeg_final_correction,
    Dp17Ocarb_permeg_ACCEPTED = coalesce(
      Dp17Ocarb_permeg_MANUAL_correction,
      Dp17Ocarb_permeg_final_correction
    )
  ) %>%
  arrange(reactor_id, IPL_num)

standardized_columns <- c(
  "MLA_sample_id", "MLA_horizon_id", "reactor_id", "IPL_num",
  "Name_full", "comments", "d17O", "d.17O", "d17O.err", "d18O",
  "d.18O", "d18O.err", "CAP.17O", "CAP17O.err", "d33", "d33.err",
  "d34", "d34.err", "d35", "d35.err", "d36", "d36.err",
  "33_mismatch", "34_mismatch", "dp17O_O2_SMOWSLAP",
  "dp18O_O2_SMOWSLAP", "Dp17O_O2_SMOWSLAP2",
  "dp17Ocarb_SMOWSLAP", "dp18Ocarb_SMOWSLAP",
  "Dp17Ocarb_permeg_final_correction", "final_correction?"
)

compiled_standardized <- compiled %>%
  mutate(`final_correction?` = is.finite(Dp17Ocarb_permeg_final_correction)) %>%
  select(all_of(standardized_columns))

# ---- Reconciliation and validation --------------------------------------

reference_origin <- tibble(
  reactor_id = 31:36,
  reference_origin = c(
    "IPL R pipeline", "IPL R pipeline", "manual Excel", "BP manual Excel",
    "compiled workbook/no independent manual", "compiled workbook/no manual"
  )
)

reconciliation <- compiled %>%
  left_join(reference_origin, by = "reactor_id") %>%
  transmute(
    reactor_id,
    IPL_num,
    MLA_sample_id,
    MLA_horizon_id,
    reference_origin,
    new_automated = Dp17Ocarb_permeg_final_correction,
    legacy_standardized = previous_accepted,
    legacy_manual = Dp17Ocarb_permeg_MANUAL_correction,
    previous_automated,
    new_minus_legacy = new_automated - legacy_standardized,
    new_minus_previous_automated = new_automated - previous_automated,
    abs_new_minus_legacy = abs(new_minus_legacy),
    within_5_permeg = abs_new_minus_legacy <= 5,
    workbook_O2_difference =
      Dp17O_O2_SMOWSLAP2 - workbook_Dp17O_O2_permeg,
    calibration_method,
    calibration_d18_slope,
    calibration_d18_intercept,
    calibration_standard_offset
  )

reconciliation_summary <- reconciliation %>%
  group_by(reactor_id, reference_origin) %>%
  summarise(
    n = n(),
    mean_difference = mean(new_minus_legacy, na.rm = TRUE),
    median_difference = median(new_minus_legacy, na.rm = TRUE),
    rmse = sqrt(mean(new_minus_legacy ^ 2, na.rm = TRUE)),
    max_abs_difference = max(abs_new_minus_legacy, na.rm = TRUE),
    n_within_5 = sum(within_5_permeg, na.rm = TRUE),
    cached_O2_mean_difference = mean(workbook_O2_difference, na.rm = TRUE),
    cached_O2_max_abs_difference = max(abs(workbook_O2_difference), na.rm = TRUE),
    .groups = "drop"
  )

o2_check <- reconciliation %>%
  filter(is.finite(workbook_O2_difference)) %>%
  summarise(max_abs = max(abs(workbook_O2_difference))) %>%
  pull(max_abs)

if (!is.finite(o2_check) || o2_check > 0.01) {
  message(
    "Fresh O2 recalculation differs from cached workbook formula results by up to ",
    signif(o2_check, 5), " per meg. Several source workbooks contain cached ",
    "external-link formulas; the fresh raw-data calculation is retained."
  )
}

if (any(!reconciliation$within_5_permeg, na.rm = TRUE)) {
  warning(
    "Some BP-protocol results differ from legacy values by >5 per meg. ",
    "These are retained and explicitly flagged in the reconciliation table."
  )
}

# ---- Outputs -------------------------------------------------------------

write_csv(
  compiled_standardized,
  file.path(output_dir, "IPL17O_all_reactors_automated_standardized.csv"),
  na = ""
)
write_csv(
  reconciliation,
  file.path(output_dir, "IPL17O_automated_reconciliation.csv"),
  na = ""
)
write_csv(
  reconciliation_summary,
  file.path(output_dir, "IPL17O_automated_reconciliation_summary.csv"),
  na = ""
)
write_csv(
  calibration_summary,
  file.path(output_dir, "IPL17O_reactor_calibration_summary.csv"),
  na = ""
)
write_csv(
  standard_offset_summary,
  file.path(output_dir, "IPL17O_standard_offset_summary.csv"),
  na = ""
)

message("Wrote IPL triple-oxygen correction products to: ", output_dir)
print(reconciliation_summary, n = Inf)
