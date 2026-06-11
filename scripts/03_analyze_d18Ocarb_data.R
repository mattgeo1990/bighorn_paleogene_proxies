# 03__analyze_d18Ocarb_data.R
# Purpose: Evaluate carbonate oxygen isotope variability among datasets

# ---- Load packages ----
library(tidyverse)
library(here)

# ---- Load processed dataset ----
BHB <- read_csv(
  here("data", "processed", "BHB_multiproxy_with_temperature.csv")
)

# ============================================================
# PART 1: Compare d18Ocarb estimates among datasets
# ============================================================

# ---- Build long-format carbonate dataset ----

d18Ocarb_long <- BHB %>%
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

bowen_cu_vpdb <- BHB %>%
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