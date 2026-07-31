as_number <- function(x) suppressWarnings(as.numeric(as.character(x)))

petersen_2019_ben <- function(D47, A = 38300, B = 0.258,
                              acid_fractionation = 0.088) {
  denominator <- D47 - B + acid_fractionation
  out <- rep(NA_real_, length(denominator))
  valid <- is.finite(denominator) & denominator > 0
  out[valid] <- sqrt(A / denominator[valid]) - 273.15
  out
}

classify_ben_reference <- function(dat) {
  dat %>%
    mutate(
      sample_upper = str_to_upper(str_squish(replace_na(Sample_Name, ""))),
      type2_upper = str_to_upper(str_squish(replace_na(Type.2, ""))),
      standard_key = str_remove(type2_upper, "\\s+HEATED$"),
      ben_reference_class = case_when(
        Type.1 == "EquilibratedGas" ~ "equilibrated gas",
        Type.1 == "HeatedGas" |
          str_detect(sample_upper, "ONLINE.*HEATED GAS") ~ "heated gas",
        str_detect(sample_upper, "HEATED") |
          str_detect(type2_upper, "HEATED") ~ "heated carbonate",
        Type.1 == "Standard" ~ "carbonate standard",
        TRUE ~ "unknown"
      )
    )
}

flag_extreme_standard_outliers <- function(dat, threshold = 3.5) {
  dat %>%
    group_by(standard_key) %>%
    mutate(
      standard_group_n = sum(is.finite(standard_residual)),
      standard_residual_median = median(standard_residual, na.rm = TRUE),
      standard_residual_mad = mad(
        standard_residual, center = standard_residual_median,
        constant = 1.4826, na.rm = TRUE
      ),
      standard_outlier_score = case_when(
        !is.finite(standard_residual) ~ NA_real_,
        is.finite(standard_residual_mad) & standard_residual_mad > 0 ~
          abs(standard_residual - standard_residual_median) /
            standard_residual_mad,
        TRUE ~ 0
      ),
      standard_extreme_outlier = standard_group_n >= 4 &
        standard_outlier_score > threshold
    ) %>%
    ungroup()
}

fit_ben_correction <- function(dat, outlier_threshold = 4,
                               fixed_parameters = NULL) {
  accepted_carbonates <- tribble(
    ~standard_key, ~D47_accepted,
    "102-GC-AZ01", 0.6120,
    "ETH-1", 0.2052,
    "ETH-2", 0.2085,
    "ETH-3", 0.6132,
    "ETH-4", 0.4505,
    "IAEA-C1", 0.3018
  )

  working <- dat %>%
    mutate(
      included = str_to_lower(str_squish(replace_na(ignoreAnalysis, ""))) ==
        "include",
      Counter = as_number(Counter),
      IPLnum = as_number(IPLnum),
      d47 = as_number(d47),
      D47_raw = as_number(D472)
    ) %>%
    classify_ben_reference()

  slope_rows <- working %>%
    filter(
      included,
      ben_reference_class %in% c("equilibrated gas", "heated carbonate"),
      is.finite(d47), is.finite(D47_raw)
    )

  slope_models <- slope_rows %>%
    group_by(ben_reference_class) %>%
    group_nest() %>%
    mutate(
      n = map_int(data, nrow),
      model = map(data, ~ lm(D47_raw ~ d47, data = .x)),
      slope = map_dbl(model, ~ unname(coef(.x)[["d47"]])),
      intercept = map_dbl(model, ~ unname(coef(.x)[["(Intercept)"]]))
    ) %>%
    select(-data, -model)

  required_slopes <- c("equilibrated gas", "heated carbonate")
  if (!all(required_slopes %in% slope_models$ben_reference_class)) {
    stop("Both equilibrated-gas and heated-carbonate slopes are required.")
  }
  estimated_linearity_slope <- mean(
    slope_models$slope[match(required_slopes,
                             slope_models$ben_reference_class)]
  )
  linearity_slope <- if (is.null(fixed_parameters)) {
    estimated_linearity_slope
  } else {
    fixed_parameters$linearity_slope
  }

  projected <- working %>%
    mutate(D47_d47_0 = D47_raw - linearity_slope * d47)

  transfer_models <- working %>%
    filter(
      included,
      ben_reference_class %in%
        c("equilibrated gas", "heated gas", "heated carbonate"),
      is.finite(d47), is.finite(D47_raw)
    ) %>%
    group_by(ben_reference_class) %>%
    group_nest() %>%
    mutate(
      n = map_int(data, nrow),
      model = map(data, ~ lm(D47_raw ~ d47, data = .x)),
      D47_measured = map_dbl(
        model, ~ unname(coef(.x)[["(Intercept)"]])
      ),
      D47_accepted = if_else(
        ben_reference_class == "equilibrated gas", 0.91985, 0.02659
      )
    ) %>%
    select(-data, -model)

  required_transfer <- c("equilibrated gas", "heated gas", "heated carbonate")
  if (!all(required_transfer %in% transfer_models$ben_reference_class)) {
    stop("Ben transfer requires equilibrated gas, heated gas, and heated carbonate.")
  }
  transfer_model <- lm(D47_accepted ~ D47_measured, data = transfer_models)
  estimated_transfer_slope <- unname(
    coef(transfer_model)[["D47_measured"]]
  )
  estimated_transfer_intercept <- unname(
    coef(transfer_model)[["(Intercept)"]]
  )
  transfer_slope <- if (is.null(fixed_parameters)) {
    estimated_transfer_slope
  } else {
    fixed_parameters$transfer_slope
  }
  transfer_intercept <- if (is.null(fixed_parameters)) {
    estimated_transfer_intercept
  } else {
    fixed_parameters$transfer_intercept
  }

  transferred <- projected %>%
    mutate(D47_transferred = transfer_intercept + transfer_slope * D47_d47_0) %>%
    left_join(accepted_carbonates, by = "standard_key") %>%
    mutate(
      standard_residual = if_else(
        included & ben_reference_class == "carbonate standard" &
          is.finite(D47_accepted),
        D47_transferred - D47_accepted,
        NA_real_
      )
    )

  standard_rows <- transferred %>%
    filter(is.finite(standard_residual)) %>%
    flag_extreme_standard_outliers(outlier_threshold)

  retained_standard_means <- standard_rows %>%
    filter(!standard_extreme_outlier) %>%
    group_by(standard_key) %>%
    summarise(
      n_retained = n(),
      mean_standard_residual = mean(standard_residual),
      .groups = "drop"
    )

  if (!nrow(retained_standard_means)) {
    stop("No carbonate standards remain after outlier filtering.")
  }
  estimated_carbonate_offset <- mean(
    retained_standard_means$mean_standard_residual
  )
  carbonate_offset <- if (is.null(fixed_parameters)) {
    estimated_carbonate_offset
  } else {
    fixed_parameters$carbonate_offset
  }

  outlier_flags <- standard_rows %>%
    select(
      Counter, IPLnum, Sample_Name, Sample, standard_key,
      standard_residual, standard_residual_median, standard_residual_mad,
      standard_outlier_score, standard_extreme_outlier
    )

  results <- transferred %>%
    left_join(
      outlier_flags %>%
        select(Counter, standard_extreme_outlier, standard_outlier_score),
      by = "Counter"
    ) %>%
    mutate(
      standard_extreme_outlier = replace_na(standard_extreme_outlier, FALSE),
      correction_qc_include = included & !standard_extreme_outlier,
      D47_corrected = D47_transferred - carbonate_offset,
      T47_C = petersen_2019_ben(D47_corrected),
      correction_method = if_else(
        is.null(fixed_parameters),
        "Ben manual-workbook procedure (parameters refit)",
        "Ben manual-workbook procedure (stored parameters)"
      ),
      linearity_slope = linearity_slope,
      transfer_slope = transfer_slope,
      transfer_intercept = transfer_intercept,
      carbonate_offset = carbonate_offset,
      temperature_A = 38300,
      temperature_B = 0.258,
      acid_fractionation = 0.088
    )

  parameters <- tibble(
    parameter = c(
      "linearity_slope", "transfer_slope", "transfer_intercept",
      "carbonate_offset", "standard_outlier_MAD_threshold",
      "estimated_linearity_slope", "estimated_transfer_slope",
      "estimated_transfer_intercept", "estimated_carbonate_offset"
    ),
    value = c(
      linearity_slope, transfer_slope, transfer_intercept,
      carbonate_offset, outlier_threshold, estimated_linearity_slope,
      estimated_transfer_slope, estimated_transfer_intercept,
      estimated_carbonate_offset
    )
  )

  list(
    results = results,
    parameters = parameters,
    slope_models = slope_models,
    transfer_means = transfer_models,
    standard_rows = standard_rows,
    retained_standard_means = retained_standard_means
  )
}

read_ben_parameters <- function(path) {
  sheet <- readxl::read_excel(
    path, sheet = "Manual Correction", col_names = FALSE,
    .name_repair = "minimal"
  )
  list(
    linearity_slope = as_number(sheet[[34]][154]),
    transfer_slope = as_number(sheet[[31]][199]),
    transfer_intercept = as_number(sheet[[31]][200]),
    carbonate_offset = as_number(sheet[[42]][142])
  )
}

read_ben_reported <- function(path) {
  sheet <- readxl::read_excel(
    path, sheet = "Manual Correction", col_names = FALSE,
    .name_repair = "minimal"
  )
  tibble(
    IPLnum = as_number(sheet[[10]]),
    Sample_Name = str_squish(as.character(sheet[[6]])),
    D47_Ben = as_number(sheet[[44]]),
    T47_Ben_C = as_number(sheet[[45]])
  ) %>%
    filter(is.finite(IPLnum), is.finite(D47_Ben)) %>%
    group_by(IPLnum, Sample_Name) %>%
    summarise(
      D47_Ben = first(D47_Ben),
      T47_Ben_C = first(T47_Ben_C),
      .groups = "drop"
    )
}

validate_ben_overlap <- function(full_dat, ben_workbook,
                                 outlier_threshold = 4) {
  ben_sheet <- readxl::read_excel(
    ben_workbook, sheet = "Manual Correction", col_names = FALSE,
    .name_repair = "minimal"
  )
  ben_samples <- str_squish(as.character(ben_sheet[[7]]))
  ben_samples <- unique(ben_samples[!is.na(ben_samples) & ben_samples != ""])

  validation_fit <- fit_ben_correction(
    full_dat %>% filter(str_squish(Sample) %in% ben_samples),
    outlier_threshold = outlier_threshold,
    fixed_parameters = read_ben_parameters(ben_workbook)
  )
  ben_reported <- read_ben_reported(ben_workbook)

  comparison <- validation_fit$results %>%
    transmute(
      Counter, IPLnum, Sample_Name,
      D47_pipeline = D47_corrected,
      T47_pipeline_C = T47_C
    ) %>%
    inner_join(
      ben_reported,
      by = c("IPLnum", "Sample_Name")
    ) %>%
    mutate(
      D47_pipeline_minus_Ben = D47_pipeline - D47_Ben,
      T47_pipeline_minus_Ben_C = T47_pipeline_C - T47_Ben_C
    )

  metrics <- comparison %>%
    summarise(
      n_overlap = n(),
      mean_D47_difference = mean(D47_pipeline_minus_Ben, na.rm = TRUE),
      RMSE_D47 = sqrt(mean(D47_pipeline_minus_Ben^2, na.rm = TRUE)),
      max_abs_D47_difference = max(
        abs(D47_pipeline_minus_Ben), na.rm = TRUE
      ),
      n_temperature = sum(is.finite(T47_pipeline_minus_Ben_C)),
      mean_temperature_difference_C = mean(
        T47_pipeline_minus_Ben_C, na.rm = TRUE
      ),
      RMSE_temperature_C = sqrt(mean(
        T47_pipeline_minus_Ben_C^2, na.rm = TRUE
      )),
      max_abs_temperature_difference_C = max(
        abs(T47_pipeline_minus_Ben_C), na.rm = TRUE
      )
    )

  list(
    fit = validation_fit,
    comparison = comparison,
    metrics = metrics
  )
}

run_ben_pipeline <- function(input_path, ben_workbook, output_dir,
                             minimum_IPL = 4800,
                             outlier_threshold = 4) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  full_dat <- readr::read_csv(input_path, show_col_types = FALSE,
                              na = c("", "NA"))
  production_input <- full_dat %>%
    filter(as_number(IPLnum) >= minimum_IPL)

  production <- fit_ben_correction(
    production_input,
    outlier_threshold = outlier_threshold,
    fixed_parameters = read_ben_parameters(ben_workbook)
  )
  validation <- validate_ben_overlap(
    full_dat, ben_workbook, outlier_threshold = outlier_threshold
  )

  readr::write_csv(
    production$results,
    file.path(output_dir, "corrected_run_level_data.csv"), na = ""
  )
  readr::write_csv(
    production$parameters,
    file.path(output_dir, "ben_correction_parameters.csv"), na = ""
  )
  readr::write_csv(
    production$standard_rows,
    file.path(output_dir, "standard_outlier_audit.csv"), na = ""
  )
  readr::write_csv(
    validation$comparison,
    file.path(output_dir, "ben_overlap_comparison.csv"), na = ""
  )
  readr::write_csv(
    validation$metrics,
    file.path(output_dir, "ben_overlap_validation_metrics.csv"), na = ""
  )

  list(production = production, validation = validation)
}
