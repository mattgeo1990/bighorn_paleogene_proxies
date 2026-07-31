# 14_build_Judd2024_LPEE_GMST_1Ma.R
# Purpose: Regrid the published Judd et al. (2024) PhanDA GMST
#          stage/event posterior percentiles to 1-Ma bins across 59-52 Ma.
#
# Important interpretation:
#   PhanDA reports three independent time slices across this interval:
#   Thanetian (59.24-56 Ma), PETM (56-55.7 Ma), and Ypresian
#   (55.7-48.07 Ma). This script does not create new 1-Ma information.
#   It represents those published estimates on a regular 1-Ma grid using
#   temporal-overlap weights. Corresponding posterior quantiles are combined
#   across source slices (a comonotonic dependence assumption) when a 1-Ma bin
#   overlaps more than one source interval.

#-- 1.) Setup ---------------------------------------------------------------
library(tidyverse)
library(here)

source(here("scripts", "helpers", "save_figure_variants.R"))

source_url <- paste0(
  "https://raw.githubusercontent.com/EJJudd/PhanDA/main/",
  "5_Outputs/PhanDA_GMSTandCO2_percentiles.csv"
)

raw_path <- here(
  "data", "raw", "JuddEtAl2024_PhanDA_GMST_CO2_percentiles.csv"
)
summary_path <- here(
  "data", "processed", "JuddEtAl2024_PhanDA_GMST_1Ma_bins.csv"
)
weights_path <- here(
  "data", "processed", "JuddEtAl2024_PhanDA_GMST_1Ma_overlap_weights.csv"
)
figure_dir <- here("figures", "global_temperature", "Judd2024_PhanDA")

dir.create(dirname(raw_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)

# Downloading is performed by R and the published file is retained verbatim.
if (!file.exists(raw_path)) {
  download.file(source_url, raw_path, mode = "wb", quiet = FALSE)
}

#-- 2.) Read the published PhanDA percentiles -------------------------------
phanda <- read_csv(raw_path, show_col_types = FALSE) %>%
  mutate(
    source_interval = factor(
      Stage,
      levels = c("Thanetian", "PETM", "Ypresian")
    )
  ) %>%
  filter(
    !is.na(source_interval),
    LowerAge > 52,
    UpperAge < 59
  ) %>%
  transmute(
    source_interval,
    source_younger_ma = UpperAge,
    source_older_ma = LowerAge,
    GMST_05, GMST_16, GMST_50, GMST_84, GMST_95
  )

if (nrow(phanda) != 3) {
  stop(
    "Expected exactly three PhanDA source intervals across 59-52 Ma; found ",
    nrow(phanda), "."
  )
}

#-- 3.) Define 1-Ma bins and temporal-overlap weights -----------------------
gmst_bins <- tibble(
  bin_older_ma = seq(59, 53, by = -1),
  bin_younger_ma = bin_older_ma - 1,
  age_ma = (bin_older_ma + bin_younger_ma) / 2
)

overlap_weights <- crossing(
  gmst_bins,
  phanda
) %>%
  mutate(
    overlap_ma = pmax(
      0,
      pmin(bin_older_ma, source_older_ma) -
        pmax(bin_younger_ma, source_younger_ma)
    )
  ) %>%
  filter(overlap_ma > 0) %>%
  group_by(bin_older_ma, bin_younger_ma, age_ma) %>%
  mutate(
    weight = overlap_ma / sum(overlap_ma),
    source_composition = paste(
      paste0(source_interval, " (", scales::percent(weight), ")"),
      collapse = "; "
    )
  ) %>%
  ungroup()

coverage_check <- overlap_weights %>%
  group_by(bin_older_ma, bin_younger_ma) %>%
  summarise(total_overlap_ma = sum(overlap_ma), .groups = "drop")

if (any(abs(coverage_check$total_overlap_ma - 1) > 1e-10)) {
  stop("At least one requested 1-Ma bin lacks complete PhanDA coverage.")
}

#-- 4.) Create overlap-weighted 1-Ma posterior summaries -------------------
# Combining corresponding quantiles assumes that posterior ranks co-vary
# perfectly among source intervals. This preserves each published PhanDA
# quantile exactly in bins that draw on a single source interval.
gmst_1ma_summary <- overlap_weights %>%
  group_by(bin_older_ma, bin_younger_ma, age_ma, source_composition) %>%
  summarise(
    GMST_05 = sum(weight * GMST_05),
    GMST_16 = sum(weight * GMST_16),
    GMST_50 = sum(weight * GMST_50),
    GMST_84 = sum(weight * GMST_84),
    GMST_95 = sum(weight * GMST_95),
    .groups = "drop"
  ) %>%
  mutate(
    bin_older_ma,
    bin_younger_ma,
    age_ma,
    source_composition,
    method = paste(
      "Overlap-weighted regridding of published PhanDA stage/event",
      "posterior percentiles under comonotonic dependence;",
      "not an independent 1-Ma data assimilation"
    )
  )

write_csv(gmst_1ma_summary, summary_path)

overlap_weights %>%
  select(
    bin_older_ma, bin_younger_ma, age_ma,
    source_interval, source_younger_ma, source_older_ma,
    overlap_ma, weight
  ) %>%
  write_csv(weights_path)

#-- 5.) Plot ---------------------------------------------------------------
petm_limits <- c(56, 55.7)

p_gmst_1ma <- ggplot(gmst_1ma_summary, aes(age_ma, GMST_50)) +
  annotate(
    "rect",
    xmin = petm_limits[2], xmax = petm_limits[1],
    ymin = -Inf, ymax = Inf,
    fill = "#D73027", alpha = 0.08
  ) +
  geom_rect(
    aes(
      xmin = bin_younger_ma, xmax = bin_older_ma,
      ymin = GMST_05, ymax = GMST_95
    ),
    fill = "#5B7FA3", alpha = 0.18
  ) +
  geom_rect(
    aes(
      xmin = bin_younger_ma, xmax = bin_older_ma,
      ymin = GMST_16, ymax = GMST_84
    ),
    fill = "#356A9A", alpha = 0.30
  ) +
  geom_step(
    color = "#173F5F", linewidth = 0.85, direction = "mid"
  ) +
  geom_point(
    shape = 21, size = 2.5, stroke = 0.6,
    color = "#173F5F", fill = "white"
  ) +
  scale_x_reverse(
    limits = c(59, 52),
    breaks = 59:52,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Age (Ma)",
    y = expression("GMST (" * degree * "C)"),
    title = "PhanDA GMST across the LPEE",
    subtitle = paste(
      "1-Ma overlap-weighted grid of Judd et al. (2024) stage/event",
      "posteriors; points are bin centers"
    ),
    caption = paste(
      "Dark/light ribbons: 16-84% and 5-95%. Pink: PhanDA PETM interval.",
      "The 1-Ma grid does not add temporal resolution beyond",
      "Thanetian, PETM, and Ypresian source estimates."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9.5),
    plot.caption = element_text(size = 8, hjust = 0),
    plot.margin = margin(5, 7, 5, 7)
  )

save_figure_variants(
  plot = p_gmst_1ma,
  base_dir = figure_dir,
  stem = "Judd2024_PhanDA_GMST_LPEE_1Ma",
  manuscript_width = 7.2,
  manuscript_height = 4.8,
  presentation_width = 10,
  presentation_height = 5.5
)

message("Wrote: ", summary_path)
message("Wrote: ", weights_path)
message("Wrote figure variants under: ", figure_dir)
