# 06_process_reference_datasets.R
# Purpose: Process external reference datasets for comparison with BHB records.


# ---- Load packages ----
library(tidyverse)
library(here)

# ---- Load Harper et al. (2024) reference datasets ----


# Harper et al. (2024) used new and published planktic foraminiferal
# geochemical data from Pacific ODP Sites 1209 and 1210, together with a
# multiproxy Bayesian hierarchical model, to reconstruct atmospheric CO2
# and SST from ~59 to 53 Ma. Their record spans long-term late Paleocene–
# early Eocene warming, the PETM, and ETM-2.
#
# Source:
# Harper, D.T., Hönisch, B., Bowen, G.J., Zeebe, R.E., Haynes, L.L.,
# Penman, D.E., and Zachos, J.C. (2024). Long- and short-term coupling of
# sea surface temperature and atmospheric CO2 during the late Paleocene and
# early Eocene. PNAS.
#
# Raw input files:
#   HarperEtAl2024_co2_out.csv
#   HarperEtAl2024_sst_out.csv
#
# Notes:
#   - Raw ages are reported in kyr before present / kyr-style units.
#   - Ages are converted to Ma by dividing by 1000.
#   - Only the posterior mean and 95% credible interval bounds are retained.
#   - CO2 and SST records share the same age grid and are merged by Age_Ma.


Harper2024_CO2 <- read.csv(
  here::here("data", "raw", "HarperEtAl2024_co2_out.csv")
)

Harper2024_SST <- read.csv(
  here::here("data", "raw", "HarperEtAl2024_sst_out.csv")
)

# Merge CO2 and SST records 
# Keep one age column, convert kyr to Ma, and retain only mean and 95% CI bounds.

Harper2024_CO2_SST <- Harper2024_CO2 %>%
  transmute(
    Age_Ma = age / 1000,
    Harper2024_mean_CO2_ppm = mean,
    Harper2024_CO2_lower95_ppm = X2.50.,
    Harper2024_CO2_upper95_ppm = X97.50.
  ) %>%
  left_join(
    Harper2024_SST %>%
      transmute(
        Age_Ma = age / 1000,
        Harper2024_mean_SST_C = mean,
        Harper2024_SST_lower95_C = X2.50.,
        Harper2024_SST_upper95_C = X97.50.
      ),
    by = "Age_Ma"
  ) %>%
  arrange(desc(Age_Ma))

# Quick checks 

str(Harper2024_CO2_SST)

summary(Harper2024_CO2_SST$Age_Ma)
summary(Harper2024_CO2_SST$Harper2024_mean_CO2_ppm)
summary(Harper2024_CO2_SST$Harper2024_mean_SST_C)



# CO2 versus age
p_Harper_CO2_age <- ggplot(
  Harper2024_CO2_SST,
  aes(x = Harper2024_mean_CO2_ppm, y = Age_Ma)
) +
  geom_errorbarh(
    aes(
      xmin = Harper2024_CO2_lower95_ppm,
      xmax = Harper2024_CO2_upper95_ppm
    ),
    height = 0,
    alpha = 0.35
  ) +
  geom_path(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_y_reverse() +
  labs(
    title = expression("Harper et al. (2024) atmospheric CO"[2]),
    x = expression("Atmospheric CO"[2] ~ "(ppm)"),
    y = "Age (Ma)"
  ) +
  theme_classic()

print(p_Harper_CO2_age)


# SST versus age
p_Harper_SST_age <- ggplot(
  Harper2024_CO2_SST,
  aes(x = Harper2024_mean_SST_C, y = Age_Ma)
) +
  geom_errorbarh(
    aes(
      xmin = Harper2024_SST_lower95_C,
      xmax = Harper2024_SST_upper95_C
    ),
    height = 0,
    alpha = 0.35
  ) +
  geom_path(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_y_reverse() +
  labs(
    title = "Harper et al. (2024) Pacific sea-surface temperature",
    x = expression("Sea-surface temperature (" * degree * "C)"),
    y = "Age (Ma)"
  ) +
  theme_classic()

print(p_Harper_SST_age)



# Save processed reference dataset

write_csv(
  Harper2024_CO2_SST,
  here("data", "processed", "Harper2024_CO2_SST_processed.csv")
)