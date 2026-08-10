# 04_model_CFB_temperatures.R
# Purpose: Fit the CFB stratigraphic temperature model under every cumulative
#          alteration-screening scenario, compare scenario sensitivity, and
#          apply one explicitly selected model to all CFB carbonate horizons.

#-- 1.) Setup ---------------------------------------------------------------
library(tidyverse)
library(here)
library(zoo)
source(here("scripts", "helpers", "save_figure_variants.R"))
source(
  here(
    "scripts", "helpers",
    "BHB_d18O_alteration_probability.R"
  )
)

dir.create(here("figures", "temperature_models"), recursive = TRUE,
           showWarnings = FALSE)

# This is the only setting that determines which scenario supplies the
# production T_model_* columns used by the soil-water reconstruction.
# The talk model uses an explicit, observation-level geochemical screen while
# all observations remain available to the diagenesis and screening figures.
primary_screening_scenario <- "talk_geochemical_screen"

talk_d18Ocarb_min_vsmow <- 20
talk_temperature_max_C <- 50

screening_scenarios <- tribble(
  ~screening_scenario_id, ~screening_scenario, ~screen_column,
  "none", "All data", "none",
  "talk_geochemical_screen",
    "Talk screen: d18Ocarb >= 20 per mil and T <= 50 C",
    "talk_geochemical_screen",
  "exclude_high_likelihood", "Exclude high likelihood",
    "exclude_high_likelihood",
  "exclude_moderate_or_higher", "Exclude moderate or higher",
    "exclude_moderate_or_higher",
  "exclude_any_alteration_indication", "Exclude any indication",
    "exclude_any_alteration_indication"
)

if (!primary_screening_scenario %in% screening_scenarios$screening_scenario_id) {
  stop("Unknown primary_screening_scenario: ", primary_screening_scenario)
}

#-- 2.) Load CFB Data and Upstream Screening Flags -------------------------
CFB_soilcarb_isotope_summary <- read_csv(
  here("data", "processed", "CFB_soilcarb_isotope_summary.csv"),
  show_col_types = FALSE
)

if (any(CFB_soilcarb_isotope_summary$section_id != "CFB")) {
  stop("Temperature modeling input must contain only section_id == 'CFB'.")
}

CFB_temperature_screening_flags <- read_csv(
  here("data", "processed", "CFB_temperature_screening_flags.csv"),
  show_col_types = FALSE
)
BHB_d18O_probability_parameters <- read_csv(
  here(
    "data", "processed",
    "BHB_d18O_alteration_probability_parameters.csv"
  ),
  show_col_types = FALSE
)
BHB_d18Ocarb_reference_mean_vsmow <-
  BHB_d18O_probability_parameters$reference_mean_d18Ocarb_vsmow[[1]]

screen_columns <- screening_scenarios$screen_column %>%
  keep(~ str_starts(.x, "exclude_"))
probability_columns <- c(
  "p_altered_preservation",
  "p_altered_preservation_lower_sensitivity",
  "p_altered_preservation_upper_sensitivity",
  "alteration_evidence_class",
  "probability_model_version"
)
required_flag_columns <- c(
  "MLA_horizon_id", "alteration_likelihood", screen_columns,
  probability_columns
)
missing_flag_columns <- setdiff(
  required_flag_columns, names(CFB_temperature_screening_flags)
)

if (length(missing_flag_columns) > 0) {
  stop(
    "CFB_temperature_screening_flags.csv is missing: ",
    paste(missing_flag_columns, collapse = ", "),
    "\nRun 03_screen_CFB_clumped_diagenesis.R first."
  )
}

if (anyDuplicated(CFB_temperature_screening_flags$MLA_horizon_id)) {
  stop("Screening flags contain duplicate MLA_horizon_id values.")
}

CFB_soilcarb_isotope_summary <- CFB_soilcarb_isotope_summary %>%
  select(-any_of(c(required_flag_columns[-1], "exclude_from_temp_model"))) %>%
  left_join(
    CFB_temperature_screening_flags %>%
      select(all_of(required_flag_columns)),
    by = "MLA_horizon_id"
  ) %>%
  mutate(
    across(all_of(screen_columns), ~ replace_na(.x, FALSE)),
    alteration_likelihood = replace_na(
      alteration_likelihood, "no_indication"
    )
  )

#-- 3.) Define Shared Temperature-Model Functions --------------------------
make_temperature_observations <- function(data, screen_column) {
  if (str_starts(screen_column, "exclude_")) {
    data <- data %>% filter(!.data[[screen_column]])
  }

  observations <- data %>%
    select(
      section_id, MLA_horizon_id, strat_height_m,
      all_of(probability_columns),
      IPLD47_mean_T47_C, IPLD47_se_T47_C,
      IPL_NuDog_d18Ocarb_VSMOW,
      CU_mean_T47_C, CU_2se_T47_C,
      CU_mean_d18Ocarb_vsmow,
      Snell_mean_T47_C, Snell_se_T47_C,
      Snell_mean_d18Ocarb_vsmow
    ) %>%
    pivot_longer(
      cols = c(IPLD47_mean_T47_C, CU_mean_T47_C, Snell_mean_T47_C),
      names_to = "source", values_to = "T_C"
    ) %>%
    mutate(
      source = recode(
        source,
        IPLD47_mean_T47_C = "U-M",
        CU_mean_T47_C = "CU",
        Snell_mean_T47_C = "Caltech"
      ),
      T_se_C = case_when(
        source == "U-M" ~ IPLD47_se_T47_C,
        source == "CU" ~ CU_2se_T47_C / 2,
        source == "Caltech" ~ Snell_se_T47_C
      ),
      d18Ocarb_vsmow = case_when(
        source == "U-M" ~ IPL_NuDog_d18Ocarb_VSMOW,
        source == "CU" ~ CU_mean_d18Ocarb_vsmow,
        source == "Caltech" ~ Snell_mean_d18Ocarb_vsmow
      ),
      p_altered_preservation = calc_d18O_alteration_probability(
        d18Ocarb_vsmow,
        BHB_d18Ocarb_reference_mean_vsmow
      ),
      p_altered_preservation_lower_sensitivity =
        p_altered_preservation,
      p_altered_preservation_upper_sensitivity =
        p_altered_preservation,
      alteration_evidence_class = case_when(
        p_altered_preservation < 0.20 ~ "low",
        p_altered_preservation < 0.50 ~ "limited",
        p_altered_preservation < 0.80 ~ "substantial",
        TRUE ~ "strong"
      ),
      probability_model_version = "BHB_d18O_trajectory_index_v2"
    ) %>%
    select(
      section_id, MLA_horizon_id, strat_height_m,
      all_of(probability_columns), source, T_C, T_se_C, d18Ocarb_vsmow
    ) %>%
    filter(!is.na(T_C), !is.na(strat_height_m)) %>%
    mutate(
      passes_d18Ocarb_screen =
        is.finite(d18Ocarb_vsmow) &
        d18Ocarb_vsmow >= talk_d18Ocarb_min_vsmow,
      passes_temperature_screen =
        is.finite(T_C) & T_C <= talk_temperature_max_C,
      passes_talk_temperature_model_screen =
        passes_d18Ocarb_screen & passes_temperature_screen,
      talk_screen_exclusion_reason = case_when(
        !is.finite(d18Ocarb_vsmow) ~ "Missing d18Ocarb VSMOW",
        d18Ocarb_vsmow < talk_d18Ocarb_min_vsmow &
          T_C > talk_temperature_max_C ~
            "d18Ocarb < 20 per mil VSMOW and T > 50 C",
        d18Ocarb_vsmow < talk_d18Ocarb_min_vsmow ~
          "d18Ocarb < 20 per mil VSMOW",
        T_C > talk_temperature_max_C ~ "T > 50 C",
        TRUE ~ NA_character_
      )
    )

  if (screen_column == "talk_geochemical_screen") {
    observations <- observations %>%
      filter(passes_talk_temperature_model_screen)
  }

  valid_se <- observations$T_se_C[
    !is.na(observations$T_se_C) & observations$T_se_C > 0
  ]
  fallback_se <- median(valid_se, na.rm = TRUE)
  if (!is.finite(fallback_se)) stop("No valid temperature uncertainties found.")

  observations %>%
    mutate(
      T_se_imputed = is.na(T_se_C) | T_se_C <= 0,
      T_se_C = if_else(T_se_imputed, fallback_se, T_se_C),
      weight = 1 / T_se_C^2
    )
}

fit_temperature_scenario <- function(scenario_id, scenario_label,
                                     screen_column, prediction_grid) {
  observations <- make_temperature_observations(
    CFB_soilcarb_isotope_summary, screen_column
  )

  # Combine co-located temperature estimates using an inverse-variance-
  # weighted mean. This is appropriate when U-M, CU, and Caltech observations
  # from the same horizon are treated as independent estimates of one shared
  # horizon temperature and their reported analytical SEs adequately describe
  # their relative precision. Weighting by 1 / SE^2 allows more precise
  # measurements to contribute more strongly without discarding less precise
  # measurements. Under this fixed-effect assumption, sqrt(1 / sum(weight))
  # is the propagated analytical SE of the combined horizon mean.
  #
  # Important limitation: this SE does not account for systematic differences
  # among laboratories, calibrations, analyzed carbonate material, or excess
  # disagreement among source estimates. It should therefore be interpreted
  # as analytical precision conditional on a common-temperature assumption,
  # not as the complete uncertainty in the original soil temperature.
  horizons <- observations %>%
    group_by(section_id, MLA_horizon_id, strat_height_m) %>%
    summarise(
      T_measured_C = weighted.mean(T_C, weight, na.rm = TRUE),
      T_measured_se_C = sqrt(1 / sum(weight, na.rm = TRUE)),
      n_T_obs = n(),
      T_sources = paste(sort(unique(source)), collapse = ", "),
      p_altered_preservation =
        mean(p_altered_preservation, na.rm = TRUE),
      p_altered_preservation_lower_sensitivity =
        mean(p_altered_preservation_lower_sensitivity, na.rm = TRUE),
      p_altered_preservation_upper_sensitivity =
        mean(p_altered_preservation_upper_sensitivity, na.rm = TRUE),
      alteration_evidence_class = case_when(
        p_altered_preservation < 0.20 ~ "low",
        p_altered_preservation < 0.50 ~ "limited",
        p_altered_preservation < 0.80 ~ "substantial",
        TRUE ~ "strong"
      ),
      probability_model_version = first(probability_model_version),
      .groups = "drop"
    ) %>%
    arrange(strat_height_m) %>%
    mutate(
      T_model_anchor_C = rollapply(
        T_measured_C, width = 3, FUN = mean, fill = NA,
        align = "center", partial = TRUE
      )
    )

  spline_fit <- smooth.spline(
    x = horizons$strat_height_m,
    y = horizons$T_model_anchor_C,
    spar = 0.35
  )
  prediction <- predict(spline_fit, x = prediction_grid$strat_height_m)
  fitted_at_data <- predict(spline_fit, x = horizons$strat_height_m)$y
  residuals <- horizons$T_measured_C - fitted_at_data
  model_min_strat_height_m <- min(horizons$strat_height_m, na.rm = TRUE)
  model_max_strat_height_m <- max(horizons$strat_height_m, na.rm = TRUE)

  # Approximate local prediction uncertainty combines analytical uncertainty
  # from the three nearest measured horizons with their local model residuals.
  local_prediction_se <- function(x, k = 3) {
    nearest <- horizons %>%
      mutate(
        distance_m = abs(strat_height_m - x),
        model_residual_C = residuals
      ) %>%
      arrange(distance_m) %>%
      slice_head(n = k)
    sqrt(mean(
      nearest$T_measured_se_C^2 + nearest$model_residual_C^2,
      na.rm = TRUE
    ))
  }

  model <- prediction_grid %>%
    mutate(
      temperature_model_position = case_when(
        strat_height_m < model_min_strat_height_m ~ "below_observed_range",
        strat_height_m > model_max_strat_height_m ~ "above_observed_range",
        TRUE ~ "interpolated"
      ),
      T_model_C = if_else(
        temperature_model_position == "interpolated",
        as.numeric(prediction$y),
        NA_real_
      ),
      T_model_se_C = map2_dbl(
        strat_height_m,
        temperature_model_position,
        ~ if (.y == "interpolated") local_prediction_se(.x) else NA_real_
      ),
      T_model_lower95_C = T_model_C - 1.96 * T_model_se_C,
      T_model_upper95_C = T_model_C + 1.96 * T_model_se_C,
      model_min_strat_height_m = model_min_strat_height_m,
      model_max_strat_height_m = model_max_strat_height_m,
      screening_scenario_id = scenario_id,
      screening_scenario = scenario_label,
      n_excluded_horizons = case_when(
        screen_column == "none" ~ 0L,
        screen_column == "talk_geochemical_screen" ~
          as.integer(
            n_distinct(
              make_temperature_observations(
                CFB_soilcarb_isotope_summary, "none"
              )$MLA_horizon_id
            ) - n_distinct(observations$MLA_horizon_id)
          ),
        TRUE ~ as.integer(
          sum(CFB_soilcarb_isotope_summary[[screen_column]], na.rm = TRUE)
        )
      ),
      uncertainty_method = paste(
        "Approximate 95% interval: analytical SE plus spline residual",
        "variance among three nearest measured horizons"
      ),
      extrapolation_policy =
        "No extrapolation; predictions outside retained data range are NA"
    ) %>%
    left_join(
      horizons,
      by = c("section_id", "MLA_horizon_id", "strat_height_m")
    ) %>%
    mutate(
      has_measured_T47 = !is.na(T_measured_C),
      n_T_obs = replace_na(n_T_obs, 0L)
    )

  list(observations = observations, horizons = horizons, model = model)
}

#-- 4.) Fit Every Screening Scenario ---------------------------------------
prediction_grid <- CFB_soilcarb_isotope_summary %>%
  distinct(section_id, MLA_horizon_id, strat_height_m) %>%
  filter(!is.na(strat_height_m)) %>%
  arrange(strat_height_m)

scenario_fits <- pmap(
  screening_scenarios,
  ~ fit_temperature_scenario(..1, ..2, ..3, prediction_grid)
)
names(scenario_fits) <- screening_scenarios$screening_scenario_id

CFB_temperature_scenario_models <- map_dfr(scenario_fits, "model") %>%
  mutate(
    screening_scenario = factor(
      screening_scenario,
      levels = screening_scenarios$screening_scenario
    )
  )

CFB_temperature_scenario_summary <- CFB_temperature_scenario_models %>%
  group_by(screening_scenario_id, screening_scenario) %>%
  summarise(
    n_excluded_horizons = first(n_excluded_horizons),
    n_temperature_horizons = sum(has_measured_T47),
    model_min_strat_height_m = first(model_min_strat_height_m),
    model_max_strat_height_m = first(model_max_strat_height_m),
    n_interpolated_horizons = sum(
      temperature_model_position == "interpolated"
    ),
    n_outside_model_range = sum(
      temperature_model_position != "interpolated"
    ),
    min_model_T_C = min(T_model_C, na.rm = TRUE),
    max_model_T_C = max(T_model_C, na.rm = TRUE),
    .groups = "drop"
  )

print(CFB_temperature_scenario_summary)

#-- 5.) Select and Apply the Production Scenario ---------------------------
primary_fit <- scenario_fits[[primary_screening_scenario]]
temp_obs <- scenario_fits[["none"]]$observations %>%
  mutate(
    used_in_primary_temperature_model =
      passes_talk_temperature_model_screen
  )
temp_horizon <- primary_fit$horizons
CFB_temperature_model <- primary_fit$model

selected_screen_column <- screening_scenarios %>%
  filter(screening_scenario_id == primary_screening_scenario) %>%
  pull(screen_column)

CFB_soilcarb_with_temperature <- CFB_soilcarb_isotope_summary %>%
  left_join(
    temp_obs %>%
      group_by(section_id, MLA_horizon_id, strat_height_m) %>%
      summarise(
        n_temperature_observations = n(),
        n_temperature_observations_passing_screen =
          sum(used_in_primary_temperature_model),
        exclude_from_temp_model =
          !any(used_in_primary_temperature_model),
        .groups = "drop"
      ),
    by = c("section_id", "MLA_horizon_id", "strat_height_m")
  ) %>%
  mutate(
    exclude_from_temp_model = replace_na(exclude_from_temp_model, FALSE),
    primary_screening_scenario = primary_screening_scenario
  ) %>%
  left_join(
    CFB_temperature_model,
    by = c("section_id", "MLA_horizon_id", "strat_height_m")
  )

#-- 6.) Plot Scenario Sensitivity with 95% Error Ribbons -------------------
scenario_colors <- c(
  "All data" = "#000000",
  "Exclude high likelihood" = "#0072B2",
  "Exclude moderate or higher" = "#E69F00",
  "Exclude any indication" = "#D55E00"
)

p_temperature_scenarios <- ggplot(
  CFB_temperature_scenario_models %>%
    filter(temperature_model_position == "interpolated"),
  aes(y = strat_height_m, group = screening_scenario)
) +
  annotate(
    "rect", xmin = -Inf, xmax = Inf, ymin = 1500, ymax = 1540,
    fill = "grey70", alpha = 0.25
  ) +
  geom_ribbon(
    aes(
      xmin = T_model_lower95_C, xmax = T_model_upper95_C,
      fill = screening_scenario
    ),
    alpha = 0.10, color = NA
  ) +
  geom_path(aes(x = T_model_C, color = screening_scenario), linewidth = 1) +
  scale_color_manual(values = scenario_colors, drop = FALSE) +
  scale_fill_manual(values = scenario_colors, drop = FALSE) +
  scale_x_continuous(
    breaks = seq(10, 75, by = 5), limits = c(10, 75),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = expression(Delta[47] * "-derived temperature (" * degree * "C)"),
    y = "CFB stratigraphic height (m)",
    color = "Alteration screen", fill = "Alteration screen",
    title = "Sensitivity of the CFB temperature model to diagenetic screening",
    subtitle = "Ribbons are approximate 95% intervals for each scenario"
  ) +
  theme_classic(base_size = 18) +
  theme(legend.position = "top", legend.box = "vertical")

save_figure_variants(
  p_temperature_scenarios, here("figures", "temperature_models"),
  "CFB_temperature_screening_scenarios", 8, 8,
  presentation_width = 6
)

# A faceted version avoids overlap and makes the uncertainty of each scenario
# independently legible while retaining identical x and y scales.
p_temperature_scenarios_faceted <- p_temperature_scenarios +
  facet_wrap(~ screening_scenario, ncol = 2) +
  guides(color = "none", fill = "none") +
  theme(legend.position = "none")

save_figure_variants(
  p_temperature_scenarios_faceted, here("figures", "temperature_models"),
  "CFB_temperature_screening_scenarios_faceted", 10, 9,
  presentation_width = 6
)

#-- 7.) Export Temperature Products ----------------------------------------
write_csv(
  temp_obs,
  here("data", "processed", "CFB_temperature_observations.csv")
)
write_csv(
  temp_horizon,
  here("data", "processed", "CFB_temperature_horizon_summary.csv")
)
write_csv(
  CFB_temperature_model,
  here("data", "processed", "CFB_temperature_model.csv")
)
write_csv(
  CFB_temperature_scenario_models,
  here("data", "processed", "CFB_temperature_scenario_models.csv")
)
write_csv(
  CFB_temperature_scenario_summary,
  here("data", "processed", "CFB_temperature_scenario_summary.csv")
)
write_csv(
  CFB_soilcarb_with_temperature,
  here("data", "processed", "CFB_soilcarb_with_temperature.csv")
)
