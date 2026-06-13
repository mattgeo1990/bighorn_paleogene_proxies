library(tidyverse)
library(here)

# ---- Load data ----
BP_corrections <- read_csv(
  here("data", "junkyard data", "IPL_D47_Bens_corrections.csv")
)

names(BP_corrections)

WN_corrections <- read_csv(
  here("data", "junkyard data", "IPL_D47_Wills_corrections.csv")
)

names(WN_corrections)
BP_summary <- read_csv(
  here("data", "junkyard data", "BHB Paleogene Summary May 2026-1.csv")
)

names(BP_summary)
# ---- Clean Will sample IDs ----
WN_corrections_clean <- WN_corrections %>%
  mutate(
    MLA_horizon_id = gsub("'", "", `'SampleID'`),
    MLA_horizon_id = sub("-(01|02|03)$", "", MLA_horizon_id)
  )

# ---- Check IPLnum overlap ----
bp_only_IPLnums <- setdiff(BP_corrections$IPLnum, WN_corrections_clean$IPLnum)
wn_only_IPLnums <- setdiff(WN_corrections_clean$IPLnum, BP_corrections$IPLnum)

bp_only_IPLnums
wn_only_IPLnums

# ---- Compare shared corrections ----
comparison <- inner_join(
  BP_corrections %>%
    select(IPLnum, BP_T47 = T47_preferred),
  WN_corrections_clean %>%
    select(IPLnum, WN_T47 = T47_preferred),
  by = "IPLnum"
) %>%
  mutate(
    diff_C = WN_T47 - BP_T47
  )

ggplot(comparison, aes(x = BP_T47, y = WN_T47)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    color = "red"
  ) +
  geom_point(size = 2) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "BP T47 preferred (°C)",
    y = "WN T47 preferred (°C)"
  )

comparison %>%
  summarise(
    n = n(),
    mean_diff_C = mean(diff_C, na.rm = TRUE),
    sd_diff_C = sd(diff_C, na.rm = TRUE),
    rmse_C = sqrt(mean(diff_C^2, na.rm = TRUE)),
    max_abs_diff_C = max(abs(diff_C), na.rm = TRUE)
  )

# ---- Pull WN-only corrections and add strat height ----
WN_only <- WN_corrections_clean %>%
  filter(IPLnum %in% wn_only_IPLnums) %>%
  select(
    IPLnum,
    MLA_horizon_id,
    IPLD47_mean_T47_C = T47_preferred
  ) %>%
  left_join(
    IPL_sample_list %>%
      select(MLA_horizon_id, strat_height_m),
    by = "MLA_horizon_id"
  ) %>%
  arrange(IPLnum)

WN_only

# ---- Check unmatched MLA horizon IDs ----
WN_only %>%
  filter(is.na(strat_height_m))

# ---- Save output ----
write_csv(
  WN_only,
  here("data", "processed", "WN_only_corrections.csv")
)