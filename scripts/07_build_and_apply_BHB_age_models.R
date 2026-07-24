# 07_build_and_apply_BHB_age_models.R
# Purpose:
#   Construct reproducible age-depth models for individual Bighorn Basin
#   sections and apply those models to processed proxy datasets that contain
#   both section_id and native strat_height_m.
#
# Data model:
#   - BHB_section_chronostratigraphy.xlsx stores where events occur.
#   - absolute_age_priors.csv stores numerical ages shared among sections.
#   - Dated ashes retain their numerical ages in the chronostratigraphy file.
#   - Proxy measurements never influence the age-depth models.

#-- 1.) Setup ---------------------------------------------------------------

library(tidyverse)
library(here)
library(readxl)
source(here("scripts", "helpers", "save_figure_variants.R"))

processed_dir <- here("data", "processed")
figure_dir <- here("figures", "age_models")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

#-- 2.) Load Chronostratigraphic Events and Absolute Ages ------------------

chronostrat_events <- read_excel(
  here(
    "data",
    "excel files",
    "BHB_section_chronostratigraphy.xlsx"
  ),
  sheet = "chronostrat_events"
) %>%
  mutate(
    strat_position_m = as.numeric(strat_position_m),
    strat_position_min_m = as.numeric(strat_position_min_m),
    strat_position_max_m = as.numeric(strat_position_max_m),
    mean_age_ma = as.numeric(mean_age_ma),
    age_error_ma = as.numeric(age_error_ma),
    candidate_correlation_tie = as.logical(candidate_correlation_tie)
  )

section_registry <- read_excel(
  here(
    "data",
    "excel files",
    "BHB_section_chronostratigraphy.xlsx"
  ),
  sheet = "section_registry"
) %>%
  mutate(
    section_base_m = as.numeric(section_base_m),
    section_top_m = as.numeric(section_top_m)
  )

absolute_age_priors <- read_csv(
  here("data", "raw", "absolute_age_priors.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    mean_age_ma = as.numeric(mean_age_ma),
    age_error_ma = as.numeric(age_error_ma)
  )

#-- 3.) Standardize Event Keys ---------------------------------------------

# Magnetochron names are already standardized in both files. Wa0 labels carry
# additional local descriptions, so they are mapped to one explicit project
# key. All other correlations remain undated unless they carry a local age.
chronostrat_events <- chronostrat_events %>%
  mutate(
    absolute_age_prior_id = case_when(
      event_name == "Base C24n.1r" ~ "BASE_C24N_1R",
      event_name == "Base C24n.2n" ~ "BASE_C24N_2N",
      event_name == "Base C24n.2r" ~ "BASE_C24N_2R",
      event_name == "Base C24n.3n" ~ "BASE_C24N_3N",
      event_name == "Base C24r" ~ "BASE_C24R",
      event_name == "Base C25n" ~ "BASE_C25N",
      event_name == "Base C25r" ~ "BASE_C25R",
      event_name == "Base C26n" ~ "BASE_C26N",
      event_name == "Base C26r" ~ "BASE_C26R",
      str_detect(event_name, regex("Base Wa0", ignore_case = TRUE)) ~
        "BASE_WA0_PETM_ONSET",
      TRUE ~ NA_character_
    )
  )

#-- 4.) Assemble Dated Tie Points by Section -------------------------------

section_age_tiepoints <- chronostrat_events %>%
  left_join(
    absolute_age_priors %>%
      transmute(
        absolute_age_prior_id,
        shared_mean_age_ma = mean_age_ma,
        shared_age_error_ma = age_error_ma,
        shared_age_error_metric = age_error_metric,
        shared_timescale = timescale,
        shared_calibration_basis = calibration_basis
      ),
    by = "absolute_age_prior_id"
  ) %>%
  mutate(
    prior_age_source = case_when(
      !is.na(mean_age_ma) ~ "section_specific_dated_event",
      !is.na(shared_mean_age_ma) ~ "shared_absolute_age_prior",
      TRUE ~ NA_character_
    ),
    prior_age_ma = coalesce(mean_age_ma, shared_mean_age_ma),
    prior_age_error_ma = coalesce(age_error_ma, shared_age_error_ma),
    prior_age_error_metric = coalesce(
      age_error_metric,
      shared_age_error_metric
    )
  ) %>%
  filter(
    candidate_correlation_tie,
    !is.na(strat_position_m),
    !is.na(prior_age_ma)
  ) %>%
  arrange(section_id, strat_position_m) %>%
  select(
    section_id,
    event_id,
    event_name,
    event_type,
    strat_position_m,
    strat_position_min_m,
    strat_position_max_m,
    position_scale,
    prior_age_ma,
    prior_age_error_ma,
    prior_age_error_metric,
    prior_age_source,
    absolute_age_prior_id,
    position_status,
    confidence,
    source_citation,
    source_locator,
    notes
  ) %>%
  mutate(
    # The published Belt Ash mean overlaps the base-C26n age only when its
    # analytical uncertainty is considered; using both central values would
    # create a local age reversal. Preserve the ash in the audit table, but do
    # not force its central value through the deterministic interpolation.
    use_in_deterministic_model = event_id != "CFB_BELT_ASH",
    model_exclusion_reason = case_when(
      event_id == "CFB_BELT_ASH" ~ paste(
        "Central age conflicts with overlying GTS2020 base C26n;",
        "retained as an uncertainty-bearing dated observation"
      ),
      TRUE ~ NA_character_
    )
  )

section_model_tiepoints <- section_age_tiepoints %>%
  filter(use_in_deterministic_model)

#-- 5.) Validate Section Model Eligibility ---------------------------------

section_model_status <- section_registry %>%
  select(section_id, section_name, section_base_m, section_top_m) %>%
  left_join(
    section_model_tiepoints %>%
      group_by(section_id) %>%
      summarise(
        n_age_tiepoints = n(),
        n_unique_strat_positions = n_distinct(strat_position_m),
        age_order_is_monotonic =
          all(diff(prior_age_ma[order(strat_position_m)]) < 0),
        lowest_dated_position_m = min(strat_position_m),
        highest_dated_position_m = max(strat_position_m),
        .groups = "drop"
      ),
    by = "section_id"
  ) %>%
  mutate(
    n_age_tiepoints = replace_na(n_age_tiepoints, 0L),
    n_unique_strat_positions = replace_na(n_unique_strat_positions, 0L),
    model_eligible =
      n_unique_strat_positions >= 2 &
      replace_na(age_order_is_monotonic, FALSE),
    model_status = case_when(
      n_unique_strat_positions < 2 ~ "insufficient_absolute_control",
      !replace_na(age_order_is_monotonic, FALSE) ~
        "non_monotonic_tiepoints",
      TRUE ~ "eligible"
    )
  )

eligible_section_ids <- section_model_status %>%
  filter(model_eligible) %>%
  pull(section_id)

#-- 6.) Define the Reusable Section-Age Predictor --------------------------

predict_section_age <- function(strat_height_m, tiepoints) {
  tiepoints <- tiepoints %>%
    arrange(strat_position_m) %>%
    distinct(strat_position_m, .keep_all = TRUE)

  x <- tiepoints$strat_position_m
  y <- tiepoints$prior_age_ma
  n_ties <- length(x)

  if (n_ties < 2) {
    return(tibble(
      Age_Ma = rep(NA_real_, length(strat_height_m)),
      age_model_position = "insufficient_absolute_control",
      distance_to_nearest_prior_m = NA_real_,
      age_control_distance_index = NA_real_
    ))
  }

  Age_Ma <- approx(x, y, xout = strat_height_m, rule = 1)$y

  below <- !is.na(strat_height_m) & strat_height_m < x[1]
  above <- !is.na(strat_height_m) & strat_height_m > x[n_ties]

  lower_slope <- (y[2] - y[1]) / (x[2] - x[1])
  upper_slope <-
    (y[n_ties] - y[n_ties - 1]) /
    (x[n_ties] - x[n_ties - 1])

  Age_Ma[below] <- y[1] + lower_slope * (strat_height_m[below] - x[1])
  Age_Ma[above] <-
    y[n_ties] + upper_slope * (strat_height_m[above] - x[n_ties])

  distance_to_nearest_prior_m <- map_dbl(
    strat_height_m,
    ~ if (is.na(.x)) NA_real_ else min(abs(.x - x))
  )

  age_model_position <- case_when(
    is.na(strat_height_m) ~ NA_character_,
    strat_height_m < x[1] ~ "extrapolated_below",
    strat_height_m > x[n_ties] ~ "extrapolated_above",
    strat_height_m %in% x ~ "tie_point",
    TRUE ~ "interpolated"
  )

  # This is a relative control-distance index, not a formal age uncertainty.
  # It equals zero at tie points, reaches one at the midpoint between bounding
  # controls, and increases with distance during extrapolation.
  age_control_distance_index <- map_dbl(strat_height_m, function(z) {
    if (is.na(z)) return(NA_real_)
    if (z < x[1]) return((x[1] - z) / (x[2] - x[1]))
    if (z > x[n_ties]) {
      return((z - x[n_ties]) / (x[n_ties] - x[n_ties - 1]))
    }
    if (z %in% x) return(0)

    lower_i <- findInterval(z, x, all.inside = TRUE)
    w <- (z - x[lower_i]) / (x[lower_i + 1] - x[lower_i])
    4 * w * (1 - w)
  })

  tibble(
    Age_Ma,
    age_model_position,
    distance_to_nearest_prior_m,
    age_control_distance_index
  )
}

#-- 7.) Build Section Age-Model Grids --------------------------------------

section_age_model_grid <- section_model_status %>%
  filter(model_eligible) %>%
  transmute(
    section_id,
    strat_height_m = map2(
      section_base_m,
      section_top_m,
      ~ seq(.x, .y, by = 1)
    )
  ) %>%
  unnest(strat_height_m) %>%
  group_by(section_id) %>%
  group_modify(~ bind_cols(
    .x,
    predict_section_age(
      .x$strat_height_m,
      section_model_tiepoints %>% filter(section_id == .y$section_id)
    )
  )) %>%
  ungroup() %>%
  mutate(age_model_id = "BHB_section_linear_GTS2020_v1")

#-- 8.) Apply Models to Standardized BHB Proxy Tables ----------------------

apply_section_age_models <- function(data) {
  required_columns <- c("section_id", "strat_height_m")
  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns) > 0) {
    stop(
      "Cannot apply section ages; missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  data %>%
    group_by(section_id) %>%
    group_modify(function(.x, .y) {
      this_section <- .y$section_id

      if (!this_section %in% eligible_section_ids) {
        return(.x %>% mutate(
          Age_Ma = NA_real_,
          age_model_position = "insufficient_absolute_control",
          distance_to_nearest_prior_m = NA_real_,
          age_control_distance_index = NA_real_,
          age_model_id = NA_character_
        ))
      }

      predictions <- predict_section_age(
        .x$strat_height_m,
        section_model_tiepoints %>% filter(section_id == this_section)
      )

      bind_cols(.x, predictions) %>%
        mutate(age_model_id = "BHB_section_linear_GTS2020_v1")
    }) %>%
    ungroup() %>%
    relocate(section_id)
}

age_application_manifest <- tribble(
  ~input_file, ~output_file,
  "CFB_soilcarb_isotope_summary.csv",
  "CFB_soilcarb_isotope_summary_age_calibrated.csv",
  "CFB_soilcarb_with_temperature.csv",
  "CFB_soilcarb_with_temperature_age_calibrated.csv",
  "CFB_soilwater_reconstruction_summary.csv",
  "CFB_soilwater_reconstruction_summary_age_calibrated.csv",
  "BHB_regional_soilcarb_reference_summary.csv",
  "BHB_regional_soilcarb_reference_summary_age_calibrated.csv",
  "SnellEtAl2013_soilcarb_summary.csv",
  "SnellEtAl2013_soilcarb_summary_age_calibrated.csv",
  "Koch_soilcarb_summary.csv",
  "Koch_soilcarb_summary_age_calibrated.csv"
)

age_application_log <- pmap_dfr(
  age_application_manifest,
  function(input_file, output_file) {
    input_path <- file.path(processed_dir, input_file)

    if (!file.exists(input_path)) {
      return(tibble(
        input_file,
        output_file,
        n_rows = NA_integer_,
        n_ages_assigned = NA_integer_,
        status = "input_missing"
      ))
    }

    age_calibrated <- read_csv(input_path, show_col_types = FALSE) %>%
      select(-any_of(c(
        "Age_Ma",
        "age_model_position",
        "distance_to_nearest_prior_m",
        "age_control_distance_index",
        "age_model_id"
      ))) %>%
      apply_section_age_models()

    write_csv(age_calibrated, file.path(processed_dir, output_file))

    tibble(
      input_file,
      output_file,
      n_rows = nrow(age_calibrated),
      n_ages_assigned = sum(!is.na(age_calibrated$Age_Ma)),
      status = "written"
    )
  }
)

#-- 9.) Export Models, Validation Tables, and Diagnostics ------------------

write_csv(
  section_age_tiepoints,
  file.path(processed_dir, "BHB_section_age_tiepoints.csv")
)

write_csv(
  section_model_tiepoints,
  file.path(processed_dir, "BHB_section_age_model_tiepoints.csv")
)

write_csv(
  section_model_status,
  file.path(processed_dir, "BHB_section_age_model_status.csv")
)

write_csv(
  section_age_model_grid,
  file.path(processed_dir, "BHB_section_age_model_grid.csv")
)

write_csv(
  age_application_log,
  file.path(processed_dir, "BHB_age_application_log.csv")
)

age_model_diagnostic <- ggplot(
  section_age_model_grid,
  aes(x = Age_Ma, y = strat_height_m)
) +
  geom_path(linewidth = 0.8) +
  geom_point(
    data = section_model_tiepoints %>%
      filter(section_id %in% eligible_section_ids),
    aes(x = prior_age_ma, y = strat_position_m),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    size = 2
  ) +
  facet_wrap(~ section_id, scales = "free_y") +
  scale_x_reverse() +
  labs(
    x = "Age (Ma)",
    y = "Native stratigraphic position (m)",
    title = "Section-specific Bighorn Basin age models",
    subtitle = "Points are dated controls; lines use piecewise-linear interpolation/extrapolation"
  ) +
  theme_classic()

save_figure_variants(
  age_model_diagnostic, figure_dir, "BHB_section_age_models",
  10, 7, presentation_width = 12
)

print(section_model_status %>% select(section_id, n_age_tiepoints, model_status))
print(age_application_log)
