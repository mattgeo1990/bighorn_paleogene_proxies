# 02_analyze_CFB_carbonate_agreement.R
# Purpose: Evaluate agreement among carbonate-isotope datasets immediately
#          after construction of the integrated CFB soil-carbonate table.

# ---- Load packages ----
library(tidyverse)
library(here)
source(here("scripts", "helpers", "save_figure_variants.R"))

# ---- Load processed dataset ----
CFB <- read_csv(
  here("data", "processed", "CFB_soilcarb_isotope_summary.csv"),
  show_col_types = FALSE
)

# ============================================================
# PART 1: Compare d18Ocarb estimates among datasets
# ============================================================

# ---- Build long-format carbonate dataset ----

d18Ocarb_long <- CFB %>%
  select(
    MLA_horizon_id,
    strat_height_m,
    
    Koch_mean_d18Ocarb_vsmow,
    Bowen_mean_d18Ocarb_vsmow,
    CU_mean_d18Ocarb_vsmow,
    Snell_mean_d18Ocarb_vsmow
  ) %>%
  pivot_longer(
    cols = c(
      Koch_mean_d18Ocarb_vsmow,
      Bowen_mean_d18Ocarb_vsmow,
      CU_mean_d18Ocarb_vsmow,
      Snell_mean_d18Ocarb_vsmow
    ),
    names_to = "source",
    values_to = "d18Ocarb_vsmow"
  ) %>%
  mutate(
    source = case_when(
      source == "Koch_mean_d18Ocarb_vsmow"  ~ "Koch",
      source == "Bowen_mean_d18Ocarb_vsmow" ~ "Bowen",
      source == "CU_mean_d18Ocarb_vsmow"    ~ "CU",
      source == "Snell_mean_d18Ocarb_vsmow" ~ "Snell"
    )
  ) %>%
  filter(!is.na(d18Ocarb_vsmow))

# ---- How many horizons have multiple datasets? ----

cat("\n=============================\n")
cat("MULTI-DATASET HORIZONS\n")
cat("=============================\n")

d18Ocarb_long %>%
  group_by(MLA_horizon_id) %>%
  summarise(
    n_datasets = n(),
    .groups = "drop"
  ) %>%
  count(n_datasets)

# ---- Variability among datasets for same horizon ----

d18Ocarb_agreement <- d18Ocarb_long %>%
  group_by(
    MLA_horizon_id,
    strat_height_m
  ) %>%
  summarise(
    n_datasets = n(),
    
    mean_d18Ocarb = mean(
      d18Ocarb_vsmow,
      na.rm = TRUE
    ),
    
    sd_d18Ocarb = sd(
      d18Ocarb_vsmow,
      na.rm = TRUE
    ),
    
    range_d18Ocarb =
      max(d18Ocarb_vsmow, na.rm = TRUE) -
      min(d18Ocarb_vsmow, na.rm = TRUE),
    
    min_d18Ocarb = min(
      d18Ocarb_vsmow,
      na.rm = TRUE
    ),
    
    max_d18Ocarb = max(
      d18Ocarb_vsmow,
      na.rm = TRUE
    ),
    
    sources = paste(
      sort(unique(source)),
      collapse = ", "
    ),
    
    .groups = "drop"
  ) %>%
  filter(n_datasets > 1) %>%
  arrange(desc(range_d18Ocarb))

# ---- Summary statistics of inter-dataset variability ----

cat("\n=============================\n")
cat("INTER-DATASET AGREEMENT\n")
cat("=============================\n")
d18Ocarb_agreement <- d18Ocarb_long %>%
  group_by(
    MLA_horizon_id,
    strat_height_m
  ) %>%
  summarise(
    n_datasets = n(),
    
    mean_d18Ocarb = mean(d18Ocarb_vsmow),
    
    sd_d18Ocarb = sd(d18Ocarb_vsmow),
    
    range_d18Ocarb =
      max(d18Ocarb_vsmow) -
      min(d18Ocarb_vsmow),
    
    sources = paste(
      sort(unique(source)),
      collapse = ", "
    ),
    
    .groups = "drop"
  ) %>%
  filter(n_datasets > 1) %>%
  arrange(desc(range_d18Ocarb))

# ---- Show worst disagreements ----

cat("\n=============================\n")
cat("TOP 20 LARGEST DISAGREEMENTS\n")
cat("=============================\n")

d18Ocarb_agreement %>%
  select(
    MLA_horizon_id,
    strat_height_m,
    n_datasets,
    range_d18Ocarb,
    sd_d18Ocarb,
    sources
  ) %>%
  slice_head(n = 20)


# test for systematic offset -----

bowen_cu_pairs <- d18Ocarb_long %>%
  filter(source %in% c("Bowen", "CU")) %>%
  select(
    MLA_horizon_id,
    source,
    d18Ocarb_vsmow
  ) %>%
  pivot_wider(
    names_from = source,
    values_from = d18Ocarb_vsmow
  ) %>%
  filter(
    !is.na(Bowen),
    !is.na(CU)
  ) %>%
  mutate(
    difference = Bowen - CU
  )

summary(bowen_cu_pairs$difference)

mean(bowen_cu_pairs$difference)

sd(bowen_cu_pairs$difference)

t.test(bowen_cu_pairs$difference)

bowen_cu_vpdb <- CFB %>%
  select(
    MLA_horizon_id,
    strat_height_m,
    Bowen_mean_d18Ocarb_vpdb,
    CU_mean_d18Ocarb_vpdb
  ) %>%
  filter(
    !is.na(Bowen_mean_d18Ocarb_vpdb),
    !is.na(CU_mean_d18Ocarb_vpdb)
  ) %>%
  mutate(
    difference = Bowen_mean_d18Ocarb_vpdb -
      CU_mean_d18Ocarb_vpdb,
    abs_difference = abs(difference)
  )

summary(bowen_cu_vpdb$difference)

mean(bowen_cu_vpdb$difference)

sd(bowen_cu_vpdb$difference)

t.test(bowen_cu_vpdb$difference)
# ---- Save agreement table ----

write_csv(
  d18Ocarb_agreement,
  here(
    "data",
    "processed",
    "d18Ocarb_dataset_agreement.csv"
  )
)

# ============================================================
# PART 2: Compare Delta47 temperatures among laboratories
# ============================================================
#
# Laboratory attribution:
#   U-M     = University of Michigan; measurements generated in this study.
#   CU      = University of Colorado; previously published measurements.
#   Caltech = Snell et al. (2013). Their supplement states that all sample
#             gases were cleaned and measured on the MAT 253 at Caltech,
#             although some carbonate digestion occurred at UCSC.
#
# Comparisons are made only when the exact same MLA_horizon_id has results
# from more than one laboratory. Nearby stratigraphic levels are not treated
# as replicates. Temperature is used rather than raw Delta47 because the
# source datasets use different historical reference frames/calibrations.

agreement_figure_dir <- here(
  "figures", "carbonate_agreement", "interlaboratory_D47"
)
dir.create(agreement_figure_dir, recursive = TRUE, showWarnings = FALSE)

D47_temperature_long <- CFB %>%
  select(
    section_id, MLA_horizon_id, strat_height_m,
    IPLD47_mean_T47_C, IPLD47_se_T47_C,
    CU_mean_T47_C, CU_2se_T47_C,
    Snell_mean_T47_C, Snell_se_T47_C
  ) %>%
  pivot_longer(
    cols = c(IPLD47_mean_T47_C, CU_mean_T47_C, Snell_mean_T47_C),
    names_to = "laboratory",
    values_to = "T47_C"
  ) %>%
  mutate(
    laboratory = recode(
      laboratory,
      IPLD47_mean_T47_C = "U-M",
      CU_mean_T47_C = "CU",
      Snell_mean_T47_C = "Caltech"
    ),
    T47_se_C = case_when(
      laboratory == "U-M" ~ IPLD47_se_T47_C,
      laboratory == "CU" ~ CU_2se_T47_C / 2,
      laboratory == "Caltech" ~ Snell_se_T47_C
    ),
    study_status = if_else(
      laboratory == "U-M", "This study", "Published data"
    )
  ) %>%
  select(
    section_id, MLA_horizon_id, strat_height_m,
    laboratory, study_status, T47_C, T47_se_C
  ) %>%
  filter(is.finite(T47_C))

D47_multilab_horizons <- D47_temperature_long %>%
  group_by(section_id, MLA_horizon_id, strat_height_m) %>%
  filter(n_distinct(laboratory) > 1) %>%
  ungroup()

# Create every available laboratory pair within a shared horizon. The
# standardized difference divides the pairwise temperature difference by the
# combined independent analytical SE. Values with absolute z <= 1.96 agree at
# the approximate 95% level, conditional on the reported analytical errors.
D47_pairwise_agreement <- D47_multilab_horizons %>%
  inner_join(
    D47_multilab_horizons,
    by = c("section_id", "MLA_horizon_id", "strat_height_m"),
    suffix = c("_1", "_2"),
    relationship = "many-to-many"
  ) %>%
  filter(laboratory_1 < laboratory_2) %>%
  mutate(
    laboratory_pair = paste(laboratory_1, laboratory_2, sep = " - "),
    difference_C = T47_C_1 - T47_C_2,
    absolute_difference_C = abs(difference_C),
    combined_se_C = sqrt(T47_se_C_1^2 + T47_se_C_2^2),
    standardized_difference_z = difference_C / combined_se_C,
    agrees_within_95pct_analytical =
      abs(standardized_difference_z) <= 1.96,
    pair_mean_C = (T47_C_1 + T47_C_2) / 2
  ) %>%
  arrange(desc(absolute_difference_C))

D47_pairwise_summary <- D47_pairwise_agreement %>%
  group_by(laboratory_pair) %>%
  summarise(
    n_shared_horizons = n(),
    mean_difference_C = mean(difference_C),
    median_difference_C = median(difference_C),
    sd_difference_C = sd(difference_C),
    mean_absolute_difference_C = mean(absolute_difference_C),
    rmse_C = sqrt(mean(difference_C^2)),
    proportion_agreeing_within_95pct_analytical =
      mean(agrees_within_95pct_analytical),
    .groups = "drop"
  )

cat("\n=============================\n")
cat("INTERLABORATORY DELTA47-TEMPERATURE AGREEMENT\n")
cat("=============================\n")
print(D47_pairwise_agreement)
print(D47_pairwise_summary)

write_csv(
  D47_temperature_long,
  here("data", "processed", "CFB_D47_temperature_by_laboratory.csv")
)
write_csv(
  D47_pairwise_agreement,
  here("data", "processed", "CFB_D47_interlaboratory_pairwise_agreement.csv")
)
write_csv(
  D47_pairwise_summary,
  here("data", "processed", "CFB_D47_interlaboratory_agreement_summary.csv")
)

if (nrow(D47_pairwise_agreement) > 0) {
  p_D47_interlab_one_to_one <- ggplot(
    D47_pairwise_agreement,
    aes(T47_C_1, T47_C_2)
  ) +
    geom_abline(slope = 1, intercept = 0, color = "grey55", linetype = 2) +
    geom_errorbar(
      aes(ymin = T47_C_2 - T47_se_C_2, ymax = T47_C_2 + T47_se_C_2),
      width = 0, color = "grey35", alpha = 0.65
    ) +
    geom_errorbarh(
      aes(xmin = T47_C_1 - T47_se_C_1, xmax = T47_C_1 + T47_se_C_1),
      height = 0, color = "grey35", alpha = 0.65
    ) +
    geom_point(aes(fill = laboratory_pair), shape = 21, size = 3) +
    facet_wrap(~ laboratory_pair) +
    coord_equal() +
    labs(
      x = expression("Laboratory 1 " * Delta[47] * " temperature (" *
                       degree * "C)"),
      y = expression("Laboratory 2 " * Delta[47] * " temperature (" *
                       degree * "C)"),
      fill = "Laboratory pair",
      title = "Independent clumped-isotope temperatures from shared horizons",
      subtitle = "Error bars are reported 1 SE; dashed line is 1:1"
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")

  p_D47_interlab_difference <- ggplot(
    D47_pairwise_agreement,
    aes(pair_mean_C, difference_C)
  ) +
    geom_hline(yintercept = 0, color = "grey55", linetype = 2) +
    geom_errorbar(
      aes(
        ymin = difference_C - combined_se_C,
        ymax = difference_C + combined_se_C
      ),
      width = 0, color = "grey35", alpha = 0.65
    ) +
    geom_point(aes(fill = laboratory_pair), shape = 21, size = 3) +
    facet_wrap(~ laboratory_pair, scales = "free") +
    labs(
      x = expression("Pair-mean " * Delta[47] * " temperature (" *
                       degree * "C)"),
      y = "Laboratory 1 - laboratory 2 (degrees C)",
      fill = "Laboratory pair",
      title = "Interlaboratory temperature differences",
      subtitle = "Error bars show the combined independent analytical SE"
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")

  save_figure_variants(
    p_D47_interlab_one_to_one, agreement_figure_dir,
    "CFB_D47_interlaboratory_one_to_one", 7, 6,
    presentation_width = 6
  )
  save_figure_variants(
    p_D47_interlab_difference, agreement_figure_dir,
    "CFB_D47_interlaboratory_difference_plot", 7, 6,
    presentation_width = 6
  )
}
