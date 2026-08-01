# 03_screen_CFB_clumped_diagenesis.R
# Purpose: Document horizon-level alteration assessments, generate cumulative
#          temperature-screening scenarios, and evaluate paired D47-D48 data.
#
# This script does not delete observations or choose the primary temperature
# model. It creates reproducible flags that the temperature-modeling script
# applies under several explicitly named sensitivity scenarios.

#-- 1.) Setup ---------------------------------------------------------------
library(tidyverse)
library(here)
source(here("scripts", "helpers", "save_figure_variants.R"))
source(
  here(
    "scripts", "helpers",
    "BHB_d18O_alteration_probability.R"
  )
)

dir.create(here("figures", "diagenetic_screening"), recursive = TRUE,
           showWarnings = FALSE)

#-- 2.) Load the Authoritative CFB Horizon Roster --------------------------
CFB_soilcarb_isotope_summary <- read_csv(
  here("data", "processed", "CFB_soilcarb_isotope_summary.csv"),
  show_col_types = FALSE
)

CFB_horizons <- CFB_soilcarb_isotope_summary %>%
  filter(section_id == "CFB") %>%
  distinct(section_id, MLA_horizon_id, strat_height_m)

if (anyDuplicated(CFB_horizons$MLA_horizon_id)) {
  stop("CFB horizon identifiers are not unique in the screening input.")
}

#-- 3.) Record the Qualitative Alteration Assessment -----------------------
# These classifications preserve the assessment developed in the original
# diagenesis_screening.R analysis. The evidence and classification can be
# revised here without changing any raw or integrated isotope data.
alteration_assessment <- tribble(
  ~MLA_horizon_id, ~alteration_likelihood,
  "PK95-SC-4",     "high",
  "PK95-SC-279",   "high",
  "PK95-SC-242",   "high",
  "PK95-SC-27",    "high",
  "PK95-SC-176",   "moderate",
  "PK95-SC-246",   "moderate",
  "PK95-SC-118up", "moderate",
  "PK95-SC-187",   "possible",
  "PK95-SC-160",   "possible",
  "PK95-SC-6",     "possible"
) %>%
  mutate(
    screening_basis = paste(
      "Qualitative alteration assessment retained from the original",
      "diagenesis-screening analysis; see diagnostic products."
    )
  )

unmatched_assessments <- anti_join(
  alteration_assessment, CFB_horizons, by = "MLA_horizon_id"
)

if (nrow(unmatched_assessments) > 0) {
  warning(
    "Alteration assessments did not match the CFB horizon roster: ",
    paste(unmatched_assessments$MLA_horizon_id, collapse = ", ")
  )
}

# Retain one row for every CFB horizon. Unflagged horizons are explicitly
# classified as having no current indication of alteration.
CFB_temperature_screening_flags <- CFB_horizons %>%
  left_join(alteration_assessment, by = "MLA_horizon_id") %>%
  mutate(
    alteration_likelihood = replace_na(
      alteration_likelihood, "no_indication"
    ),
    screening_basis = replace_na(
      screening_basis,
      "No alteration indication assigned in the current screening assessment."
    ),
    exclude_high_likelihood = alteration_likelihood == "high",
    exclude_moderate_or_higher =
      alteration_likelihood %in% c("high", "moderate"),
    exclude_any_alteration_indication =
      alteration_likelihood %in% c("high", "moderate", "possible")
  ) %>%
  arrange(strat_height_m)

write_csv(
  CFB_temperature_screening_flags,
  here("data", "processed", "CFB_temperature_screening_flags.csv")
)

#-- 4.) Summarize and Plot the Screening Scenarios ------------------------
screening_scenario_summary <- tribble(
  ~screening_scenario, ~screen_column,
  "All data", "none",
  "Exclude high likelihood", "exclude_high_likelihood",
  "Exclude moderate or higher", "exclude_moderate_or_higher",
  "Exclude any indication", "exclude_any_alteration_indication"
) %>%
  mutate(
    n_total_horizons = nrow(CFB_temperature_screening_flags),
    n_excluded_horizons = map_int(
      screen_column,
      ~ if (.x == "none") 0L else
        sum(CFB_temperature_screening_flags[[.x]], na.rm = TRUE)
    ),
    n_retained_horizons = n_total_horizons - n_excluded_horizons
  )

write_csv(
  screening_scenario_summary,
  here("data", "processed", "CFB_temperature_screening_scenario_summary.csv")
)

p_screening_flags <- CFB_temperature_screening_flags %>%
  filter(alteration_likelihood != "no_indication") %>%
  mutate(
    alteration_likelihood = factor(
      alteration_likelihood,
      levels = c("possible", "moderate", "high")
    )
  ) %>%
  ggplot(aes(x = alteration_likelihood, y = strat_height_m,
             color = alteration_likelihood)) +
  geom_point(size = 3) +
  geom_text(aes(label = MLA_horizon_id), hjust = -0.1, size = 3,
            show.legend = FALSE) +
  scale_color_manual(values = c(
    possible = "#E69F00", moderate = "#D55E00", high = "#A50026"
  )) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.65))) +
  labs(
    x = "Assigned alteration likelihood",
    y = "CFB stratigraphic height (m)",
    color = "Alteration likelihood",
    title = "Horizons flagged by the CFB clumped-isotope screen"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

save_figure_variants(
  p_screening_flags, here("figures", "diagenetic_screening"),
  "CFB_temperature_screening_flags", 7, 6, presentation_width = 6
)

#-- 5.) Plot T47-d18Ocarb Space with Parent-Water Contours -----------------
# Following the T47-d18O framework used by Henkes et al. (2018), pair each
# clumped-isotope temperature with carbonate d18O and overlay contours of
# constant calculated parent-water d18O. Movement across contours represents
# a change in inferred water composition, whereas movement along a contour can
# be produced by temperature-dependent calcite-water fractionation at constant
# water d18O. These curves are equilibrium reference contours, not fitted
# temporal trajectories or a kinetic burial-history model.
#
# References:
#   Henkes et al. (2018), Earth Planet. Sci. Lett. 490, 40-50,
#   doi:10.1016/j.epsl.2018.02.001.
#   Calcite-water fractionation follows Kim & O'Neil (1997), the same
#   relationship used in 05_reconstruct_CFB_soilwater.R.

calc_d18Ocarb_from_water_KO97 <- function(T_C, d18Owater_vsmow) {
  alpha_calcite_water <- exp(
    (18.03 * (1000 / (T_C + 273.15)) - 32.42) / 1000
  )
  ((d18Owater_vsmow + 1000) * alpha_calcite_water) - 1000
}

CFB_d18O_T47_observations <- bind_rows(
  CFB_soilcarb_isotope_summary %>%
    transmute(
      section_id, MLA_horizon_id, strat_height_m,
      source = "U-M micrite", carbonate_type = "Pedogenic micrite",
      T47_C = IPLD47_mean_T47_C,
      T47_se_C = IPLD47_se_T47_C,
      d18Ocarb_vsmow = IPL_NuDog_d18Ocarb_VSMOW,
      d18Ocarb_se_vsmow = NA_real_
    ),
  CFB_soilcarb_isotope_summary %>%
    transmute(
      section_id, MLA_horizon_id, strat_height_m,
      source = "CU micrite", carbonate_type = "Pedogenic micrite",
      T47_C = CU_mean_T47_C,
      T47_se_C = CU_2se_T47_C / 2,
      d18Ocarb_vsmow = CU_mean_d18Ocarb_vsmow,
      # VSMOW = 1.03091 * VPDB + 30.91, so its SE scales by 1.03091.
      d18Ocarb_se_vsmow = 1.03091 * CU_se_d18Ocarb_vpdb
    ),
  CFB_soilcarb_isotope_summary %>%
    transmute(
      section_id, MLA_horizon_id, strat_height_m,
      source = "Caltech micrite", carbonate_type = "Pedogenic micrite",
      T47_C = Snell_mean_T47_C,
      T47_se_C = Snell_se_T47_C,
      d18Ocarb_vsmow = Snell_mean_d18Ocarb_vsmow,
      d18Ocarb_se_vsmow = NA_real_
    ),
  read_csv(
    here("data", "processed", "spar_altered_horizon_summary.csv"),
    show_col_types = FALSE
  ) %>%
    transmute(
      section_id = replace_na(section_id, "CFB"),
      MLA_horizon_id, strat_height_m,
      source = sources,
      carbonate_type = if_else(
        str_detect(str_to_lower(paste(sample_groups, sample_types)), "spar"),
        "Spar", "Altered carbonate"
      ),
      T47_C = mean_T47_C,
      T47_se_C = se_T47_C,
      d18Ocarb_vsmow = mean_d18Ocarb_vsmow,
      d18Ocarb_se_vsmow = se_d18Ocarb_vsmow
    )
) %>%
  filter(
    section_id == "CFB",
    !is.na(T47_C),
    !is.na(d18Ocarb_vsmow)
  ) %>%
  distinct(source, MLA_horizon_id, T47_C, d18Ocarb_vsmow, .keep_all = TRUE) %>%
  left_join(
    CFB_temperature_screening_flags %>%
      select(MLA_horizon_id, alteration_likelihood),
    by = "MLA_horizon_id"
  ) %>%
  mutate(
    alteration_likelihood = replace_na(
      alteration_likelihood, "not_assessed"
    )
  )

d18Owater_contour_values <- seq(-15, 10, by = 5)

CFB_d18Owater_equilibrium_contours <- expand_grid(
  T47_C = seq(10, 120, by = 0.5),
  d18Owater_vsmow = d18Owater_contour_values
) %>%
  mutate(
    d18Ocarb_vsmow = calc_d18Ocarb_from_water_KO97(
      T47_C, d18Owater_vsmow
    ),
    contour_label = paste0(
      "starting d18Ow = ", d18Owater_vsmow, " per mil"
    )
  )

contour_labels <- CFB_d18Owater_equilibrium_contours %>%
  filter(T47_C == 110)

carbonate_shapes <- c(
  "Pedogenic micrite" = 21,
  "Spar" = 24,
  "Altered carbonate" = 22
)

p_CFB_d18O_T47_contours <- ggplot() +
  geom_line(
    data = CFB_d18Owater_equilibrium_contours,
    aes(T47_C, d18Ocarb_vsmow, group = d18Owater_vsmow),
    color = "grey55", linewidth = 0.55
  ) +
  geom_text(
    data = contour_labels,
    aes(T47_C, d18Ocarb_vsmow, label = contour_label),
    hjust = -0.05, size = 2.8, color = "grey35"
  ) +
  geom_errorbarh(
    data = CFB_d18O_T47_observations,
    aes(
      xmin = T47_C - T47_se_C,
      xmax = T47_C + T47_se_C,
      y = d18Ocarb_vsmow
    ),
    height = 0, linewidth = 0.35, alpha = 0.45, na.rm = TRUE
  ) +
  geom_errorbar(
    data = CFB_d18O_T47_observations,
    aes(
      x = T47_C,
      ymin = d18Ocarb_vsmow - d18Ocarb_se_vsmow,
      ymax = d18Ocarb_vsmow + d18Ocarb_se_vsmow
    ),
    width = 0, linewidth = 0.35, alpha = 0.45, na.rm = TRUE
  ) +
  geom_point(
    data = CFB_d18O_T47_observations,
    aes(
      T47_C, d18Ocarb_vsmow,
      shape = carbonate_type,
      fill = alteration_likelihood
    ),
    size = 2.8, color = "black", stroke = 0.6, alpha = 0.9
  ) +
  scale_shape_manual(values = carbonate_shapes, drop = FALSE) +
  scale_fill_manual(
    values = c(
      no_indication = "white", possible = "#F0E442",
      moderate = "#E69F00", high = "#D55E00",
      not_assessed = "grey70"
    ),
    labels = c(
      no_indication = "No indication", possible = "Possible",
      moderate = "Moderate", high = "High",
      not_assessed = "Not assessed"
    ),
    na.value = "grey70"
  ) +
  scale_x_continuous(
    limits = c(10, 130), breaks = seq(20, 120, by = 20),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = expression(delta^18 * O[carbonate] ~ "(per mil VSMOW)"),
    shape = "Carbonate material",
    fill = "Alteration likelihood",
    title = expression("CFB carbonate " * delta^18 * O * " versus " * T[47]),
    subtitle =
      "Grey curves show constant calculated parent-water d18O (per mil VSMOW)",
    caption = paste(
      "T47-d18O framework: Henkes et al. (2018).",
      "Calcite-water fractionation: Kim and O'Neil (1997)."
    )
  ) +
  theme_classic(base_size = 12) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 3)),
    shape = guide_legend(override.aes = list(fill = "white", size = 3))
  ) +
  theme(
    legend.position = "top",
    legend.box = "vertical",
    plot.margin = margin(6, 65, 6, 6)
  )

write_csv(
  CFB_d18O_T47_observations,
  here("data", "processed", "CFB_d18O_T47_screening_observations.csv")
)
write_csv(
  CFB_d18Owater_equilibrium_contours,
  here("data", "processed", "CFB_d18Owater_equilibrium_contours.csv")
)
save_figure_variants(
  p_CFB_d18O_T47_contours, here("figures", "diagenetic_screening"),
  "CFB_d18O_T47_contours", 9, 7.5, presentation_width = 6
)

#-- 6.) Restore the Monte Carlo Burial-Fluid Trajectory Plot ---------------
# This sensitivity experiment is retained from the original exploratory
# diagenesis_screening.R script. It draws plausible primary conditions,
# maximum burial/recrystallization temperatures, and progressive enrichment
# of burial-fluid d18O. Each path then calculates the equilibrium carbonate
# d18O expected as temperature and fluid composition evolve together.
#
# These paths are illustrative end members, not posterior simulations, fitted
# histories, or probabilities of alteration. Their purpose is to show which
# combinations of heating and water-rock evolution could enter the observed
# T47-d18Ocarb field. A fixed seed and an exported parameter table make the
# experiment exactly reproducible.
set.seed(123)

n_burial_paths <- 500

CFB_burial_trajectory_parameters <- tibble(
  path_id = seq_len(n_burial_paths),
  initial_T_C = runif(n_burial_paths, 15, 40),
  initial_d18Owater_vsmow = runif(n_burial_paths, -8, -2),
  max_T_C = runif(n_burial_paths, 60, 120),
  total_d18Owater_shift = runif(n_burial_paths, 0, 14)
)

CFB_burial_mc_trajectories <- CFB_burial_trajectory_parameters %>%
  mutate(
    trajectory = pmap(
      list(initial_T_C, initial_d18Owater_vsmow, max_T_C,
           total_d18Owater_shift),
      function(initial_T_C, initial_d18Owater_vsmow, max_T_C,
               total_d18Owater_shift) {
        tibble(T47_C = seq(initial_T_C, max_T_C, length.out = 100)) %>%
          mutate(
            reaction_progress = (T47_C - initial_T_C) /
              (max_T_C - initial_T_C),
            d18Owater_vsmow = initial_d18Owater_vsmow +
              total_d18Owater_shift * reaction_progress,
            d18Ocarb_vsmow = calc_d18Ocarb_from_water_KO97(
              T47_C, d18Owater_vsmow
            )
          )
      }
    )
  ) %>%
  unnest(trajectory)

p_CFB_d18O_T47_monte_carlo <- ggplot() +
  geom_line(
    data = CFB_burial_mc_trajectories,
    aes(
      x = T47_C,
      y = d18Ocarb_vsmow,
      group = path_id,
      color = total_d18Owater_shift
    ),
    alpha = 0.18, linewidth = 0.45
  ) +
  geom_errorbarh(
    data = CFB_d18O_T47_observations,
    aes(
      xmin = T47_C - T47_se_C,
      xmax = T47_C + T47_se_C,
      y = d18Ocarb_vsmow
    ),
    height = 0, linewidth = 0.35, alpha = 0.35, na.rm = TRUE
  ) +
  geom_errorbar(
    data = CFB_d18O_T47_observations,
    aes(
      x = T47_C,
      ymin = d18Ocarb_vsmow - d18Ocarb_se_vsmow,
      ymax = d18Ocarb_vsmow + d18Ocarb_se_vsmow
    ),
    width = 0, linewidth = 0.35, alpha = 0.35, na.rm = TRUE
  ) +
  geom_point(
    data = CFB_d18O_T47_observations,
    aes(T47_C, d18Ocarb_vsmow, shape = carbonate_type),
    fill = "white", color = "black", size = 2.8, stroke = 0.7
  ) +
  scale_shape_manual(values = carbonate_shapes, drop = FALSE) +
  scale_color_viridis_c(
    option = "C",
    name = expression("Total " * delta^18 * O[water] * " shift (per mil)")
  ) +
  scale_x_continuous(
    limits = c(10, 125), breaks = seq(20, 120, by = 20),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = expression(delta^18 * O[carbonate] ~ "(per mil VSMOW)"),
    shape = "Carbonate material",
    title = expression("Monte Carlo burial-fluid paths in " * T[47] *
                         "-" * delta^18 * O * " space"),
    subtitle = paste(
      "500 illustrative paths span heating and progressive burial-fluid",
      "18O enrichment; paths are not fitted histories"
    ),
    caption =
      "Calcite-water fractionation: Kim and O'Neil (1997); random seed = 123."
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.box = "vertical"
  )

write_csv(
  CFB_burial_trajectory_parameters,
  here("data", "processed", "CFB_burial_trajectory_parameters.csv")
)

save_figure_variants(
  p_CFB_d18O_T47_monte_carlo, here("figures", "diagenetic_screening"),
  "CFB_d18O_T47_Monte_Carlo_burial_trajectories", 9, 7.5,
  presentation_width = 6
)

#-- 7.) Summarize Monte Carlo Paths as Nested Simulation Envelopes --------
# Interpolate every simulated path onto a common 1-degree temperature grid,
# then summarize carbonate d18O across all paths that span each temperature.
# The resulting bands are central quantile envelopes of the specified Monte
# Carlo parameter space. They are not frequentist confidence intervals and do
# not measure the probability that the true burial path lies within a band.
CFB_burial_mc_common_grid <- CFB_burial_trajectory_parameters %>%
  mutate(
    trajectory = pmap(
      list(initial_T_C, initial_d18Owater_vsmow, max_T_C,
           total_d18Owater_shift),
      function(initial_T_C, initial_d18Owater_vsmow, max_T_C,
               total_d18Owater_shift) {
        T_grid_C <- seq(ceiling(initial_T_C), floor(max_T_C), by = 1)
        tibble(T47_C = T_grid_C) %>%
          mutate(
            reaction_progress = (T47_C - initial_T_C) /
              (max_T_C - initial_T_C),
            d18Owater_vsmow = initial_d18Owater_vsmow +
              total_d18Owater_shift * reaction_progress,
            d18Ocarb_vsmow = calc_d18Ocarb_from_water_KO97(
              T47_C, d18Owater_vsmow
            )
          )
      }
    )
  ) %>%
  select(path_id, trajectory) %>%
  unnest(trajectory)

CFB_burial_mc_envelope_summary <- CFB_burial_mc_common_grid %>%
  group_by(T47_C) %>%
  summarise(
    n_paths = n_distinct(path_id),
    mean_d18Ocarb_vsmow = mean(d18Ocarb_vsmow),
    median_d18Ocarb_vsmow = median(d18Ocarb_vsmow),
    lower95_d18Ocarb_vsmow = quantile(d18Ocarb_vsmow, 0.025),
    upper95_d18Ocarb_vsmow = quantile(d18Ocarb_vsmow, 0.975),
    lower80_d18Ocarb_vsmow = quantile(d18Ocarb_vsmow, 0.10),
    upper80_d18Ocarb_vsmow = quantile(d18Ocarb_vsmow, 0.90),
    lower50_d18Ocarb_vsmow = quantile(d18Ocarb_vsmow, 0.25),
    upper50_d18Ocarb_vsmow = quantile(d18Ocarb_vsmow, 0.75),
    .groups = "drop"
  ) %>%
  # Avoid presenting unstable tails supported by only a few simulated paths.
  filter(n_paths >= 100)

p_CFB_d18O_T47_monte_carlo_envelopes <- ggplot(
  CFB_burial_mc_envelope_summary,
  aes(x = T47_C)
) +
  geom_ribbon(
    aes(
      ymin = lower95_d18Ocarb_vsmow,
      ymax = upper95_d18Ocarb_vsmow,
      fill = "95%"
    ),
    alpha = 0.28
  ) +
  geom_ribbon(
    aes(
      ymin = lower80_d18Ocarb_vsmow,
      ymax = upper80_d18Ocarb_vsmow,
      fill = "80%"
    ),
    alpha = 0.42
  ) +
  geom_ribbon(
    aes(
      ymin = lower50_d18Ocarb_vsmow,
      ymax = upper50_d18Ocarb_vsmow,
      fill = "50%"
    ),
    alpha = 0.60
  ) +
  geom_line(
    aes(y = mean_d18Ocarb_vsmow),
    color = "#3F007D", linewidth = 1.2
  ) +
  geom_errorbarh(
    data = CFB_d18O_T47_observations,
    aes(
      xmin = T47_C - T47_se_C,
      xmax = T47_C + T47_se_C,
      y = d18Ocarb_vsmow
    ),
    inherit.aes = FALSE,
    height = 0, linewidth = 0.35, alpha = 0.35, na.rm = TRUE
  ) +
  geom_errorbar(
    data = CFB_d18O_T47_observations,
    aes(
      x = T47_C,
      ymin = d18Ocarb_vsmow - d18Ocarb_se_vsmow,
      ymax = d18Ocarb_vsmow + d18Ocarb_se_vsmow
    ),
    inherit.aes = FALSE,
    width = 0, linewidth = 0.35, alpha = 0.35, na.rm = TRUE
  ) +
  geom_point(
    data = CFB_d18O_T47_observations,
    aes(T47_C, d18Ocarb_vsmow, shape = carbonate_type),
    inherit.aes = FALSE,
    fill = "white", color = "black", size = 2.8, stroke = 0.7
  ) +
  scale_fill_manual(
    values = c("95%" = "#CBC9E2", "80%" = "#9E9AC8", "50%" = "#6A51A3"),
    breaks = c("50%", "80%", "95%"),
    name = "Central simulation interval"
  ) +
  scale_shape_manual(values = carbonate_shapes, drop = FALSE) +
  scale_x_continuous(
    limits = c(10, 125), breaks = seq(20, 120, by = 20),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = expression(delta^18 * O[carbonate] ~ "(per mil VSMOW)"),
    shape = "Carbonate material",
    title = expression("Mean Monte Carlo burial path and " * delta^18 *
                         O * " simulation envelopes"),
    subtitle = paste0(
      "Line is the ensemble mean; nested bands contain the central 50%, 80%,\n",
      "and 95% of simulated paths at each temperature"
    ),
    caption = paste(
      "Intervals summarize the specified uniform parameter ranges and are",
      "not fitted confidence intervals; at least 100 paths support each temperature."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.box = "vertical"
  )

write_csv(
  CFB_burial_mc_envelope_summary,
  here("data", "processed", "CFB_burial_trajectory_envelope_summary.csv")
)

save_figure_variants(
  p_CFB_d18O_T47_monte_carlo_envelopes,
  here("figures", "diagenetic_screening"),
  "CFB_d18O_T47_Monte_Carlo_mean_and_envelopes", 9, 7.5,
  presentation_width = 6
)

#-- 8.) Evaluate Paired D47-D48 Measurements -------------------------------
dual_clumped_raw <- read_csv(
  here(
    "data", "processed",
    "Nu_Dog_Clump_Session22_Oct 2025-July2026_Matlab-2_MLAcleaned.csv"
  ),
  show_col_types = FALSE
)

IPL_D47_lookup <- read_csv(
  here("data", "raw", "IPL_D47_BHB_Pg_Summary_June2026.csv"),
  show_col_types = FALSE
)

CFB_dual_clumped_analyses <- dual_clumped_raw %>%
  filter(
    Type.1 == "Sample", ignoreAnalysis == "include",
    !is.na(D472), !is.na(D484), !is.na(D47.err3), !is.na(D48.err5),
    D47.err3 > 0, D48.err5 > 0
  ) %>%
  transmute(
    IPLnum, analysis_name = Sample_Name, DateTime,
    D47 = D472, D47_se = D47.err3,
    D48 = D484, D48_se = D48.err5
  ) %>%
  inner_join(
    IPL_D47_lookup %>%
      select(IPLnum, MLA_sample_id, MLA_horizon_id, strat_height_m,
             T47_preferred),
    by = "IPLnum"
  ) %>%
  inner_join(
    CFB_temperature_screening_flags %>%
      select(MLA_horizon_id, alteration_likelihood),
    by = "MLA_horizon_id"
  )

CFB_D47_D48_screening_summary <- CFB_dual_clumped_analyses %>%
  mutate(D47_weight = 1 / D47_se^2, D48_weight = 1 / D48_se^2) %>%
  group_by(MLA_horizon_id, strat_height_m, alteration_likelihood) %>%
  summarise(
    n_analyses = n(),
    D47_mean = weighted.mean(D47, D47_weight, na.rm = TRUE),
    D47_se = sqrt(1 / sum(D47_weight, na.rm = TRUE)),
    D48_mean = weighted.mean(D48, D48_weight, na.rm = TRUE),
    D48_se = sqrt(1 / sum(D48_weight, na.rm = TRUE)),
    .groups = "drop"
  )

write_csv(
  CFB_D47_D48_screening_summary,
  here("data", "processed", "CFB_D47_D48_screening_summary.csv")
)

D47_D48_scenario_data <- map2_dfr(
  screening_scenario_summary$screening_scenario,
  screening_scenario_summary$screen_column,
  function(scenario_label, flag_column) {
    scenario_data <- CFB_D47_D48_screening_summary
    if (flag_column != "none") {
      excluded <- CFB_temperature_screening_flags %>%
        filter(.data[[flag_column]]) %>%
        pull(MLA_horizon_id)
      scenario_data <- scenario_data %>%
        filter(!MLA_horizon_id %in% excluded)
    }
    scenario_data %>% mutate(screening_scenario = scenario_label)
  }
) %>%
  mutate(
    screening_scenario = factor(
      screening_scenario,
      levels = screening_scenario_summary$screening_scenario
    )
  )

scenario_colors <- c(
  "All data" = "#000000",
  "Exclude high likelihood" = "#0072B2",
  "Exclude moderate or higher" = "#E69F00",
  "Exclude any indication" = "#D55E00"
)

p_D47_D48 <- ggplot() +
  geom_point(
    data = CFB_dual_clumped_analyses,
    aes(D47, D48), color = "grey65", size = 1.1, alpha = 0.25
  ) +
  geom_smooth(
    data = D47_D48_scenario_data,
    aes(D47_mean, D48_mean, color = screening_scenario,
        fill = screening_scenario, group = screening_scenario),
    method = "lm", formula = y ~ x, se = TRUE, level = 0.95,
    linewidth = 1, alpha = 0.10
  ) +
  geom_errorbarh(
    data = CFB_D47_D48_screening_summary,
    aes(xmin = D47_mean - D47_se, xmax = D47_mean + D47_se, y = D48_mean),
    height = 0, color = "grey45", linewidth = 0.4
  ) +
  geom_errorbar(
    data = CFB_D47_D48_screening_summary,
    aes(x = D47_mean, ymin = D48_mean - D48_se,
        ymax = D48_mean + D48_se),
    width = 0, color = "grey45", linewidth = 0.4
  ) +
  geom_point(
    data = CFB_D47_D48_screening_summary,
    aes(D47_mean, D48_mean), shape = 21, fill = "white", size = 2.6
  ) +
  scale_color_manual(values = scenario_colors, drop = FALSE) +
  scale_fill_manual(values = scenario_colors, drop = FALSE) +
  labs(
    x = expression(Delta[47]), y = expression(Delta[48]),
    color = "Screening scenario", fill = "Screening scenario",
    title = expression("CFB dual clumped-isotope screening: " * Delta[47] *
                         " versus " * Delta[48]),
    subtitle = "Scenario fits include 95% confidence intervals"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

save_figure_variants(
  p_D47_D48, here("figures", "diagenetic_screening"),
  "CFB_D47_D48_screening_scenarios", 8, 7, presentation_width = 6
)

#-- 9.) Estimate a d18O-Only Alteration-Trajectory Index -------------------
# The continuous value below is a transparent screening index, not a trained
# classifier or calibrated posterior probability. It is based only on the
# position of each point along the low-d18Ocarb direction relative to the
# pooled BHB pedogenic-micrite distribution. Qualitative petrographic classes,
# D47-D48 residuals, and temperature plausibility are not inputs.

calc_d18Owater_from_carb_KO97 <- function(T_C, d18Ocarb_vsmow) {
  alpha_calcite_water <- exp(
    (18.03 * (1000 / (T_C + 273.15)) - 32.42) / 1000
  )
  ((d18Ocarb_vsmow + 1000) / alpha_calcite_water) - 1000
}

set.seed(20260723)

# Propagate T47 and carbonate-d18O analytical uncertainty into inferred
# parent-water d18O. Missing SE values receive explicit conservative defaults.
CFB_isotopic_thermal_consistency <- CFB_d18O_T47_observations %>%
  filter(carbonate_type == "Pedogenic micrite") %>%
  mutate(
    T47_se_for_mc = if_else(
      is.finite(T47_se_C) & T47_se_C > 0, T47_se_C, 3
    ),
    d18Ocarb_se_for_mc = if_else(
      is.finite(d18Ocarb_se_vsmow) & d18Ocarb_se_vsmow > 0,
      d18Ocarb_se_vsmow, 0.15
    ),
    mc = pmap(
      list(T47_C, T47_se_for_mc, d18Ocarb_vsmow, d18Ocarb_se_for_mc),
      function(T47_C, T47_se_C, d18Ocarb_vsmow, d18Ocarb_se_vsmow) {
        T_draw <- pmax(rnorm(5000, T47_C, T47_se_C), 0.1)
        d18Ocarb_draw <- rnorm(
          5000, d18Ocarb_vsmow, d18Ocarb_se_vsmow
        )
        calc_d18Owater_from_carb_KO97(T_draw, d18Ocarb_draw)
      }
    ),
    inferred_d18Owater_median = map_dbl(mc, median, na.rm = TRUE),
    inferred_d18Owater_sd = map_dbl(mc, sd, na.rm = TRUE)
  ) %>%
  select(
    MLA_horizon_id, source,
    inferred_d18Owater_median, inferred_d18Owater_sd
  )

# Combine laboratories at horizon level. Between-laboratory spread is retained
# in the total uncertainty instead of treating several laboratories as wholly
# independent measurements of preservation.
CFB_isotopic_thermal_horizon <- CFB_isotopic_thermal_consistency %>%
  group_by(MLA_horizon_id) %>%
  summarise(
    n_isotopic_thermal_sources = n(),
    inferred_d18Owater_between_source_sd = if_else(
      n() > 1, sd(inferred_d18Owater_median, na.rm = TRUE), 0
    ),
    inferred_d18Owater_measurement_sd = sqrt(mean(
      inferred_d18Owater_sd^2, na.rm = TRUE
    )),
    inferred_d18Owater_median = median(
      inferred_d18Owater_median, na.rm = TRUE
    ),
    inferred_d18Owater_total_sd = sqrt(
      inferred_d18Owater_measurement_sd^2 +
        inferred_d18Owater_between_source_sd^2
    ),
    .groups = "drop"
  ) %>%
  mutate(
    # The broad soil-water prior is N(-5, 3.5 per mil). Measurement uncertainty
    # is added in quadrature before expressing inconsistency as a z score.
    isotopic_thermal_z = (
      inferred_d18Owater_median - (-5)
    ) / sqrt(3.5^2 + inferred_d18Owater_total_sd^2),
    isotopic_thermal_z_lower_sensitivity = (
      inferred_d18Owater_median - (-5)
    ) / sqrt(4.5^2 + inferred_d18Owater_total_sd^2),
    isotopic_thermal_z_upper_sensitivity = (
      inferred_d18Owater_median - (-5)
    ) / sqrt(2.5^2 + inferred_d18Owater_total_sd^2)
  )

# Fit an iterative errors-in-variables approximation to the horizon means.
# Each iteration includes propagated x uncertainty through the current slope.
# The standardized residual includes D47 and D48 analytical uncertainty plus
# robust residual scatter. It is a consistency diagnostic, not proof that an
# observation following the line is primary.
fit_D47_D48_consistency <- function(data, n_iter = 8) {
  initial_fit <- lm(D48_mean ~ D47_mean, data = data)
  slope <- unname(coef(initial_fit)[["D47_mean"]])

  for (iteration in seq_len(n_iter)) {
    variance_i <- data$D48_se^2 + slope^2 * data$D47_se^2
    fit <- lm(
      D48_mean ~ D47_mean,
      data = data,
      weights = 1 / pmax(variance_i, .Machine$double.eps)
    )
    slope <- unname(coef(fit)[["D47_mean"]])
  }

  residual_raw <- residuals(fit)
  robust_scatter <- 1.4826 * median(
    abs(residual_raw - median(residual_raw)), na.rm = TRUE
  )
  if (!is.finite(robust_scatter) || robust_scatter == 0) {
    robust_scatter <- sqrt(mean(residual_raw^2, na.rm = TRUE))
  }

  data %>%
    mutate(
      D48_fitted = predict(fit, newdata = data),
      D47_D48_residual = D48_mean - D48_fitted,
      D47_D48_residual_sd = sqrt(
        D48_se^2 + slope^2 * D47_se^2 + robust_scatter^2
      ),
      D47_D48_std_residual =
        D47_D48_residual / D47_D48_residual_sd
    )
}

CFB_D47_D48_consistency <- fit_D47_D48_consistency(
  CFB_D47_D48_screening_summary
) %>%
  select(
    MLA_horizon_id, D47_mean, D47_se, D48_mean, D48_se,
    D48_fitted, D47_D48_residual, D47_D48_residual_sd,
    D47_D48_std_residual
  )

# Define the reference distribution from all available BHB pedogenic-micrite
# T47 observations, not only the new IPL measurements. CFB laboratory
# observations remain separate points; the MCP values extend the geographic
# coverage without duplicating CFB horizons.
BHB_regional_T47_reference <- read_csv(
  here(
    "data", "processed",
    "BHB_regional_soilcarb_reference_summary.csv"
  ),
  show_col_types = FALSE
) %>%
  filter(
    section_id == "MCP",
    is.finite(T47_C),
    is.finite(d18Ocarb_vsmow)
  ) %>%
  transmute(
    section_id,
    MLA_horizon_id,
    source = dataset,
    carbonate_type = "Pedogenic micrite",
    T47_C,
    T47_se_C,
    d18Ocarb_vsmow,
    d18Ocarb_se_vsmow = NA_real_
  )

BHB_d18O_probability_reference <- bind_rows(
  CFB_d18O_T47_observations %>%
    filter(
      carbonate_type == "Pedogenic micrite",
      is.finite(T47_C),
      is.finite(d18Ocarb_vsmow)
    ) %>%
    select(
      section_id, MLA_horizon_id, source, carbonate_type,
      T47_C, T47_se_C, d18Ocarb_vsmow, d18Ocarb_se_vsmow
    ),
  BHB_regional_T47_reference
)

BHB_d18Ocarb_reference_mean_vsmow <- mean(
  BHB_d18O_probability_reference$d18Ocarb_vsmow,
  na.rm = TRUE
)

BHB_d18O_probability_parameters <- tibble(
  probability_model_version = "BHB_d18O_trajectory_index_v2",
  n_reference_observations = nrow(BHB_d18O_probability_reference),
  reference_mean_d18Ocarb_vsmow = BHB_d18Ocarb_reference_mean_vsmow,
  altered_anchor_d18Ocarb_vsmow = 20,
  probability_at_reference_mean = 0.05,
  probability_at_altered_anchor = 0.95,
  reference_population = paste(
    "All BHB pedogenic-micrite observations with paired T47 and",
    "d18Ocarb: IPL, CU, Caltech/CFB, and Snell/MCP"
  )
)

CFB_d18O_horizon_probability <- CFB_d18O_T47_observations %>%
  filter(
    carbonate_type == "Pedogenic micrite",
    is.finite(d18Ocarb_vsmow)
  ) %>%
  group_by(MLA_horizon_id) %>%
  summarise(
    alteration_reference_d18Ocarb_vsmow =
      mean(d18Ocarb_vsmow, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    p_altered_preservation = calc_d18O_alteration_probability(
      alteration_reference_d18Ocarb_vsmow,
      BHB_d18Ocarb_reference_mean_vsmow
    )
  )

CFB_temperature_screening_flags <- CFB_temperature_screening_flags %>%
  left_join(CFB_isotopic_thermal_horizon, by = "MLA_horizon_id") %>%
  left_join(CFB_D47_D48_consistency, by = "MLA_horizon_id") %>%
  left_join(CFB_d18O_horizon_probability, by = "MLA_horizon_id") %>%
  mutate(
    petrographic_prior_basis = NA_character_,
    p_altered_preservation_lower_sensitivity = p_altered_preservation,
    p_altered_preservation_upper_sensitivity = p_altered_preservation,
    alteration_evidence_class = case_when(
      p_altered_preservation < 0.20 ~ "low",
      p_altered_preservation < 0.50 ~ "limited",
      p_altered_preservation < 0.80 ~ "substantial",
      TRUE ~ "strong"
    ),
    p_climate_inconsistent = NA_real_,
    probability_model_version = "BHB_d18O_trajectory_index_v2",
    screening_basis = paste(
      screening_basis,
      "Alteration probability is a d18Ocarb-only trajectory index:",
      "5% at the pooled BHB pedogenic-micrite mean and 95% at",
      "20 per mil VSMOW. Qualitative petrographic classes, D47-D48,",
      "and temperature plausibility are excluded from this probability."
    )
  ) %>%
  arrange(strat_height_m)

# Assign every plotted point directly from its own d18Ocarb value. Carbonate
# material labels remain available as plot shapes, but they do not contribute
# to the probability estimate.
CFB_d18O_T47_observations <- CFB_d18O_T47_observations %>%
  select(-any_of(c(
    "p_altered_preservation",
    "p_altered_preservation_lower_sensitivity",
    "p_altered_preservation_upper_sensitivity",
    "alteration_evidence_class",
    "probability_model_version"
  ))) %>%
  mutate(
    probability_basis = paste(
      "Point-level d18Ocarb trajectory index relative to pooled BHB",
      "pedogenic-micrite mean"
    ),
    p_altered_preservation = calc_d18O_alteration_probability(
      d18Ocarb_vsmow,
      BHB_d18Ocarb_reference_mean_vsmow
    ),
    p_altered_preservation_lower_sensitivity = p_altered_preservation,
    p_altered_preservation_upper_sensitivity = p_altered_preservation,
    probability_model_version = "BHB_d18O_trajectory_index_v2",
    alteration_evidence_class = case_when(
      p_altered_preservation < 0.20 ~ "low",
      p_altered_preservation < 0.50 ~ "limited",
      p_altered_preservation < 0.80 ~ "substantial",
      TRUE ~ "strong"
    )
  )

write_csv(
  BHB_d18O_probability_parameters,
  here("data", "processed", "CFB_alteration_probability_parameters.csv")
)
write_csv(
  BHB_d18O_probability_parameters,
  here("data", "processed", "BHB_d18O_alteration_probability_parameters.csv")
)
write_csv(
  CFB_temperature_screening_flags,
  here("data", "processed", "CFB_temperature_screening_flags.csv")
)
write_csv(
  CFB_d18O_T47_observations,
  here("data", "processed", "CFB_d18O_T47_screening_observations.csv")
)

alteration_probability_colors <- c(
  "#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"
)

p_CFB_alteration_probability_T_d18O <- ggplot() +
  geom_hline(
    yintercept = BHB_d18Ocarb_reference_mean_vsmow,
    color = "#2166AC", linewidth = 0.7, linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 20,
    color = "#B2182B", linewidth = 0.7, linetype = "dashed"
  ) +
  geom_line(
    data = CFB_d18Owater_equilibrium_contours,
    aes(T47_C, d18Ocarb_vsmow, group = d18Owater_vsmow),
    color = "grey72", linewidth = 0.5
  ) +
  geom_label(
    data = contour_labels,
    aes(
      T47_C, d18Ocarb_vsmow,
      label = paste0("d18Ow = ", d18Owater_vsmow, " per mil")
    ),
    hjust = 1.04, size = 2.25, color = "grey35",
    fill = scales::alpha("white", 0.78),
    label.size = 0, label.padding = unit(0.06, "lines")
  ) +
  geom_errorbarh(
    data = CFB_d18O_T47_observations,
    aes(
      xmin = T47_C - T47_se_C, xmax = T47_C + T47_se_C,
      y = d18Ocarb_vsmow
    ),
    height = 0, color = "grey45", linewidth = 0.3, alpha = 0.45,
    na.rm = TRUE
  ) +
  geom_point(
    data = CFB_d18O_T47_observations,
    aes(
      T47_C, d18Ocarb_vsmow, shape = carbonate_type,
      fill = p_altered_preservation
    ),
    color = "black", size = 3, stroke = 0.65
  ) +
  scale_shape_manual(values = carbonate_shapes, drop = FALSE) +
  scale_fill_gradientn(
    colors = alteration_probability_colors,
    limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::label_percent(accuracy = 1),
    name = "Modeled probability\nof alteration"
  ) +
  scale_x_continuous(
    limits = c(10, 130), breaks = seq(20, 120, by = 20),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = expression(delta^18 * O[carbonate] ~ "(per mil VSMOW)"),
    shape = "Carbonate material",
    title = expression(
      "CFB alteration probability in " * T[47] * "-" *
        delta^18 * O * " space"
    ),
    subtitle = paste0(
      "Point fill is based only on d18Ocarb: 5% at the pooled BHB mean (",
      round(BHB_d18Ocarb_reference_mean_vsmow, 2),
      " per mil) and 95% at 20 per mil"
    ),
    caption = paste(
      "Dashed blue line = pooled BHB micrite mean; dashed red line =",
      "95% alteration-probability anchor. This is a trajectory index,",
      "not a calibrated posterior probability."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", legend.box = "vertical")

save_figure_variants(
  p_CFB_alteration_probability_T_d18O,
  here("figures", "diagenetic_screening"),
  "CFB_alteration_probability_T47_d18O",
  9, 7.5, presentation_width = 6
)

BHB_d18O_T47_probability_observations <- bind_rows(
  CFB_d18O_T47_observations %>%
    mutate(plot_dataset = if_else(
      str_detect(source, regex("U-M|IPL", ignore_case = TRUE)),
      "This study", "Published CFB"
    )),
  BHB_regional_T47_reference %>%
    mutate(
      p_altered_preservation = calc_d18O_alteration_probability(
        d18Ocarb_vsmow,
        BHB_d18Ocarb_reference_mean_vsmow
      ),
      plot_dataset = "Published MCP"
    )
)

p_BHB_alteration_probability_T_d18O <- ggplot() +
  annotate(
    "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 20,
    fill = "#B2182B", alpha = 0.045
  ) +
  geom_hline(
    yintercept = BHB_d18Ocarb_reference_mean_vsmow,
    color = "#2166AC", linewidth = 0.8, linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 20,
    color = "#B2182B", linewidth = 0.8, linetype = "dashed"
  ) +
  geom_line(
    data = CFB_d18Owater_equilibrium_contours,
    aes(T47_C, d18Ocarb_vsmow, group = d18Owater_vsmow),
    color = "grey78", linewidth = 0.45
  ) +
  geom_label(
    data = contour_labels,
    aes(
      T47_C, d18Ocarb_vsmow,
      label = paste0("d18Ow = ", d18Owater_vsmow, " per mil")
    ),
    hjust = 1.04, size = 2.25, color = "grey35",
    fill = scales::alpha("white", 0.78),
    label.size = 0, label.padding = unit(0.06, "lines")
  ) +
  geom_errorbarh(
    data = BHB_d18O_T47_probability_observations,
    aes(
      xmin = T47_C - T47_se_C, xmax = T47_C + T47_se_C,
      y = d18Ocarb_vsmow
    ),
    height = 0, color = "grey45", linewidth = 0.3, alpha = 0.45,
    na.rm = TRUE
  ) +
  geom_point(
    data = BHB_d18O_T47_probability_observations,
    aes(
      T47_C, d18Ocarb_vsmow,
      shape = carbonate_type,
      fill = p_altered_preservation
    ),
    color = "black", size = 3, stroke = 0.65
  ) +
  scale_shape_manual(values = carbonate_shapes, drop = FALSE) +
  scale_fill_gradientn(
    colors = alteration_probability_colors,
    limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::label_percent(accuracy = 1),
    name = "d18O trajectory\nP(altered)"
  ) +
  scale_x_continuous(
    limits = c(10, 130), breaks = seq(20, 120, by = 20),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    breaks = seq(10, 30, by = 5),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(ylim = c(10, 30), clip = "on") +
  labs(
    x = expression(Delta[47] * " temperature (" * degree * "C)"),
    y = expression(delta^18 * O[carbonate] ~ "(per mil VSMOW)"),
    shape = "Carbonate material",
    title = expression(
      "BHB " * T[47] * "-" * delta^18 * O[carbonate] *
        " alteration-trajectory index"
    ),
    subtitle = paste0(
      "5% at pooled BHB micrite mean (",
      round(BHB_d18Ocarb_reference_mean_vsmow, 2),
      " per mil); 95% at 20 per mil"
    ),
    caption = paste0(
      "Index uses d18Ocarb only; petrography, D47-D48 residuals, and\n",
      "temperature plausibility are excluded."
    )
  ) +
  theme_classic(base_size = 12) +
  guides(
    shape = guide_legend(
      order = 1, nrow = 1,
      override.aes = list(fill = "white", size = 3)
    ),
    fill = guide_colorbar(
      order = 2,
      direction = "horizontal",
      barwidth = grid::unit(7.5, "cm"),
      barheight = grid::unit(0.35, "cm"),
      title.position = "top"
    )
  ) +
  theme(
    plot.title = element_text(size = 11),
    plot.subtitle = element_text(size = 8.5),
    plot.caption = element_text(size = 7.5, lineheight = 1.05),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.justification = "center",
    legend.margin = margin(1, 1, 1, 1),
    plot.margin = margin(5, 8, 5, 5)
  )

save_figure_variants(
  p_BHB_alteration_probability_T_d18O,
  here("figures", "diagenetic_screening"),
  "BHB_alteration_probability_T47_d18O",
  manuscript_width = 6,
  manuscript_height = 6,
  presentation_width = 6,
  presentation_height = 6
)
