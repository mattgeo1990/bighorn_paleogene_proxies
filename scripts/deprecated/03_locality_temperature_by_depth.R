# 03_locality_temperature_by_depth.R
# Summarize clumped_data temperatures by locality and plot each calibration by depth.

library(dplyr)
library(ggplot2)
library(here)
library(readr)
library(readxl)
library(stringr)

standardize_clumped_locality <- function(sample_name) {
  if (is.na(sample_name) || sample_name == "") {
    return(NA_character_)
  }

  if (str_detect(sample_name, "SC-118up")) {
    return("SC-118 upper-PK95")
  }

  if (str_detect(sample_name, "(^|-)242(-|$)") && !str_detect(sample_name, "SC-")) {
    return("SC-242-PK95")
  }

  locality_number <- str_extract(sample_name, "(?<=SC-)[0-9]+")

  if (is.na(locality_number)) {
    return(NA_character_)
  }

  paste0("SC-", locality_number, "-PK95")
}

clean_sample_id <- function(sample_id) {
  str_remove_all(sample_id, "^'+|'+$")
}

build_depth_lookup <- function() {
  metadata_depths <- read_csv(
    here("data", "raw", "SandCoulee_Polecat_nodules.csv"),
    show_col_types = FALSE
  ) %>%
    transmute(
      locality_id = Sample_ID,
      depth_meta_m = Strat_m_Bowen
    ) %>%
    filter(str_detect(locality_id, "^SC-"), !is.na(depth_meta_m)) %>%
    group_by(locality_id) %>%
    summarise(depth_meta_m = first(depth_meta_m), .groups = "drop")

  koch_depths <- read_csv(
    here("data", "raw", "Koch_SC_nodules_isotopes.csv"),
    show_col_types = FALSE
  ) %>%
    transmute(
      locality_id = case_when(
        horizon_id == "SC-118-UP-PK95" ~ "SC-118 upper-PK95",
        TRUE ~ horizon_id
      ),
      depth_koch_m = Strat_m
    ) %>%
    filter(!is.na(locality_id), !is.na(depth_koch_m)) %>%
    group_by(locality_id) %>%
    summarise(depth_koch_m = first(depth_koch_m), .groups = "drop")

  full_join(metadata_depths, koch_depths, by = "locality_id") %>%
    mutate(
      # SC-80 differs between the metadata table (1492 m) and Koch (1495 m).
      # The Koch table ties the SC-80 horizon directly to 1495 m, so use that here.
      depth_locality_m = case_when(
        locality_id == "SC-80-PK95" ~ depth_koch_m,
        !is.na(depth_meta_m) ~ depth_meta_m,
        TRUE ~ depth_koch_m
      ),
      depth_source = case_when(
        locality_id == "SC-80-PK95" ~ "Koch_SC_nodules_isotopes.csv (SC-80 override)",
        !is.na(depth_meta_m) ~ "SandCoulee_Polecat_nodules.csv",
        !is.na(depth_koch_m) ~ "Koch_SC_nodules_isotopes.csv",
        TRUE ~ NA_character_
      )
    ) %>%
    select(locality_id, depth_locality_m, depth_source)
}

build_locality_plot <- function(summary_data, model_name, x_limits, output_path) {
  model_slug <- tolower(model_name)

  plot_data <- summary_data %>%
    filter(temperature_model == model_name) %>%
    mutate(locality_label = str_remove(locality_id, "-PK95$"))

  plot_obj <- ggplot(plot_data, aes(x = mean_temp_c, y = depth_locality_m)) +
    geom_errorbar(
      aes(
        xmin = mean_temp_c - two_se_temp_c,
        xmax = mean_temp_c + two_se_temp_c
      ),
      orientation = "y",
      width = 0,
      color = "gray45",
      linewidth = 0.5,
      na.rm = TRUE
    ) +
    geom_point(color = "firebrick", size = 2.8) +
    geom_text(
      aes(label = locality_label),
      nudge_x = diff(x_limits) * 0.03,
      hjust = 0,
      size = 3,
      check_overlap = TRUE
    ) +
    scale_x_continuous(
      limits = x_limits,
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = paste("Average temperature by locality:", model_name),
      subtitle = "Error bars show +/- 2 SE across locality measurements",
      x = expression(paste("Mean clumped-isotope temperature (", degree, "C)")),
      y = "Locality stratigraphic height (m)"
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray90"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.margin = margin(5.5, 72, 5.5, 5.5)
    )

  ggsave(
    filename = output_path,
    plot = plot_obj,
    width = 8.5,
    height = 6,
    dpi = 600
  )

  message("Wrote ", model_slug, " plot to: ", output_path)
}

clumped_data <- read_csv(
  here("clumped_data.csv"),
  show_col_types = FALSE
)

clumped_data_filtered <- clumped_data %>%
  filter(!str_detect(`Sample Name`, regex("spar", ignore_case = TRUE)))

clumped_long <- bind_rows(
  clumped_data_filtered %>%
    transmute(
      sample_name = `Sample Name`,
      temperature_model = "Peterson",
      temperature_c = `T(47) Peterson`
    ),
  clumped_data_filtered %>%
    transmute(
      sample_name = `Sample Name`,
      temperature_model = "Anderson",
      temperature_c = `T(47) Anderson`
    )
) %>%
  mutate(
    locality_id = vapply(sample_name, standardize_clumped_locality, character(1))
  ) %>%
  filter(!is.na(temperature_c))

missing_localities <- clumped_long %>%
  filter(is.na(locality_id)) %>%
  distinct(sample_name) %>%
  pull(sample_name)

if (length(missing_localities) > 0) {
  stop(
    "Unable to assign locality IDs for: ",
    paste(missing_localities, collapse = ", ")
  )
}

depth_lookup <- build_depth_lookup()

locality_temperature_summary <- clumped_long %>%
  group_by(temperature_model, locality_id) %>%
  summarise(
    n_measurements = n(),
    mean_temp_c = mean(temperature_c),
    sd_temp_c = if_else(n() > 1, sd(temperature_c), NA_real_),
    min_temp_c = min(temperature_c),
    max_temp_c = max(temperature_c),
    .groups = "drop"
  ) %>%
  left_join(depth_lookup, by = "locality_id") %>%
  arrange(temperature_model, depth_locality_m)

replicate_mask <- locality_temperature_summary$n_measurements > 1

locality_temperature_summary <- locality_temperature_summary %>%
  mutate(
    se_temp_c = NA_real_,
    two_se_temp_c = NA_real_,
    ci95_half_width_c = NA_real_
  )

locality_temperature_summary$se_temp_c[replicate_mask] <-
  locality_temperature_summary$sd_temp_c[replicate_mask] /
  sqrt(locality_temperature_summary$n_measurements[replicate_mask])

locality_temperature_summary$two_se_temp_c[replicate_mask] <-
  2 * locality_temperature_summary$se_temp_c[replicate_mask]

locality_temperature_summary$ci95_half_width_c[replicate_mask] <-
  qt(
    0.975,
    df = locality_temperature_summary$n_measurements[replicate_mask] - 1
  ) * locality_temperature_summary$se_temp_c[replicate_mask]

missing_depths <- locality_temperature_summary %>%
  filter(is.na(depth_locality_m)) %>%
  distinct(locality_id) %>%
  pull(locality_id)

if (length(missing_depths) > 0) {
  stop(
    "Missing locality depths for: ",
    paste(missing_depths, collapse = ", ")
  )
}

summary_output <- here("data", "processed", "locality_temperature_summary_clumped_data.csv")
sample_error_output <- here("data", "processed", "sample_temperature_error_table_clumped_data.csv")
peterson_plot_output <- here("figures", "locality_temperature_by_depth_peterson.png")
anderson_plot_output <- here("figures", "locality_temperature_by_depth_anderson.png")

write_csv(locality_temperature_summary, summary_output)

sample_temperature_error_table <- clumped_data_filtered %>%
  transmute(
    sample_name = `Sample Name`,
    locality_id = vapply(`Sample Name`, standardize_clumped_locality, character(1)),
    peterson_temp_c = `T(47) Peterson`,
    anderson_temp_c = `T(47) Anderson`
  ) %>%
  group_by(sample_name, locality_id) %>%
  summarise(
    n_rows = n(),
    peterson_n = sum(!is.na(peterson_temp_c)),
    peterson_mean_c = mean(peterson_temp_c, na.rm = TRUE),
    peterson_sd_c = if_else(peterson_n > 1, sd(peterson_temp_c, na.rm = TRUE), NA_real_),
    anderson_n = sum(!is.na(anderson_temp_c)),
    anderson_mean_c = mean(anderson_temp_c, na.rm = TRUE),
    anderson_sd_c = if_else(anderson_n > 1, sd(anderson_temp_c, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  arrange(locality_id, sample_name)

peterson_mask <- sample_temperature_error_table$peterson_n > 1
anderson_mask <- sample_temperature_error_table$anderson_n > 1

sample_temperature_error_table <- sample_temperature_error_table %>%
  mutate(
    peterson_se_c = NA_real_,
    peterson_two_se_c = NA_real_,
    peterson_ci95_half_width_c = NA_real_,
    anderson_se_c = NA_real_,
    anderson_two_se_c = NA_real_,
    anderson_ci95_half_width_c = NA_real_
  )

sample_temperature_error_table$peterson_se_c[peterson_mask] <-
  sample_temperature_error_table$peterson_sd_c[peterson_mask] /
  sqrt(sample_temperature_error_table$peterson_n[peterson_mask])

sample_temperature_error_table$peterson_two_se_c[peterson_mask] <-
  2 * sample_temperature_error_table$peterson_se_c[peterson_mask]

sample_temperature_error_table$peterson_ci95_half_width_c[peterson_mask] <-
  qt(0.975, df = sample_temperature_error_table$peterson_n[peterson_mask] - 1) *
  sample_temperature_error_table$peterson_se_c[peterson_mask]

sample_temperature_error_table$anderson_se_c[anderson_mask] <-
  sample_temperature_error_table$anderson_sd_c[anderson_mask] /
  sqrt(sample_temperature_error_table$anderson_n[anderson_mask])

sample_temperature_error_table$anderson_two_se_c[anderson_mask] <-
  2 * sample_temperature_error_table$anderson_se_c[anderson_mask]

sample_temperature_error_table$anderson_ci95_half_width_c[anderson_mask] <-
  qt(0.975, df = sample_temperature_error_table$anderson_n[anderson_mask] - 1) *
  sample_temperature_error_table$anderson_se_c[anderson_mask]

write_csv(sample_temperature_error_table, sample_error_output)

matthew_bhb <- read_excel(
  here("data", "processed", "Matthew's BHB Carbs.xlsx"),
  .name_repair = "minimal"
) %>%
  transmute(
    sample_name = clean_sample_id(`'SampleID'`),
    sample_type = `'Type1'`,
    status = `'Status'`,
    t_d47_cdes_c = `'T(D47)CDES'`,
    td47_icdes_plus_c = `'TD47iCDES+'`,
    t_d47_icdes_c = `'T(D47)iCDES'`
  ) %>%
  filter(
    sample_type == "'Sample'",
    status == "'include'",
    !str_detect(sample_name, regex("spar", ignore_case = TRUE))
  )

matthew_bhb_long <- bind_rows(
  matthew_bhb %>%
    transmute(
      sample_name = sample_name,
      temperature_model = "T(D47)CDES",
      temperature_c = t_d47_cdes_c
    ),
  matthew_bhb %>%
    transmute(
      sample_name = sample_name,
      temperature_model = "TD47iCDES+",
      temperature_c = td47_icdes_plus_c
    ),
  matthew_bhb %>%
    transmute(
      sample_name = sample_name,
      temperature_model = "T(D47)iCDES",
      temperature_c = t_d47_icdes_c
    )
) %>%
  mutate(
    locality_id = vapply(sample_name, standardize_clumped_locality, character(1))
  ) %>%
  filter(!is.na(temperature_c))

unmatched_matthew_samples <- matthew_bhb_long %>%
  filter(is.na(locality_id)) %>%
  distinct(sample_name) %>%
  pull(sample_name)

if (length(unmatched_matthew_samples) > 0) {
  message(
    "Dropped workbook sample IDs without locality/depth matches: ",
    paste(unmatched_matthew_samples, collapse = ", ")
  )
}

matthew_bhb_long <- matthew_bhb_long %>%
  filter(!is.na(locality_id))

matthew_bhb_summary <- matthew_bhb_long %>%
  group_by(temperature_model, locality_id) %>%
  summarise(
    n_measurements = n(),
    mean_temp_c = mean(temperature_c),
    sd_temp_c = if_else(n() > 1, sd(temperature_c), NA_real_),
    min_temp_c = min(temperature_c),
    max_temp_c = max(temperature_c),
    .groups = "drop"
  ) %>%
  left_join(depth_lookup, by = "locality_id") %>%
  arrange(depth_locality_m, temperature_model)

matthew_replicate_mask <- matthew_bhb_summary$n_measurements > 1

matthew_bhb_summary <- matthew_bhb_summary %>%
  mutate(
    se_temp_c = NA_real_,
    two_se_temp_c = NA_real_,
    ci95_half_width_c = NA_real_
  )

matthew_bhb_summary$se_temp_c[matthew_replicate_mask] <-
  matthew_bhb_summary$sd_temp_c[matthew_replicate_mask] /
  sqrt(matthew_bhb_summary$n_measurements[matthew_replicate_mask])

matthew_bhb_summary$two_se_temp_c[matthew_replicate_mask] <-
  2 * matthew_bhb_summary$se_temp_c[matthew_replicate_mask]

matthew_bhb_summary$ci95_half_width_c[matthew_replicate_mask] <-
  qt(
    0.975,
    df = matthew_bhb_summary$n_measurements[matthew_replicate_mask] - 1
  ) * matthew_bhb_summary$se_temp_c[matthew_replicate_mask]

missing_matthew_depths <- matthew_bhb_summary %>%
  filter(is.na(depth_locality_m)) %>%
  distinct(locality_id) %>%
  pull(locality_id)

if (length(missing_matthew_depths) > 0) {
  stop(
    "Missing workbook locality depths for: ",
    paste(missing_matthew_depths, collapse = ", ")
  )
}

matthew_summary_output <- here(
  "data",
  "processed",
  "matthews_bhb_carbs_ai_aj_ak_locality_summary.csv"
)
matthew_plot_output <- here(
  "figures",
  "matthews_bhb_carbs_ai_aj_ak_by_depth.png"
)

write_csv(matthew_bhb_summary, matthew_summary_output)

matthew_x_lower <- min(
  ifelse(
    is.na(matthew_bhb_summary$two_se_temp_c),
    matthew_bhb_summary$mean_temp_c,
    matthew_bhb_summary$mean_temp_c - matthew_bhb_summary$two_se_temp_c
  )
)
matthew_x_upper <- max(
  ifelse(
    is.na(matthew_bhb_summary$two_se_temp_c),
    matthew_bhb_summary$mean_temp_c,
    matthew_bhb_summary$mean_temp_c + matthew_bhb_summary$two_se_temp_c
  )
)
matthew_x_padding <- max(2, (matthew_x_upper - matthew_x_lower) * 0.08)
matthew_x_limits <- c(
  matthew_x_lower - matthew_x_padding,
  matthew_x_upper + (2 * matthew_x_padding)
)

model_offsets <- c(
  "T(D47)CDES" = -8,
  "TD47iCDES+" = 0,
  "T(D47)iCDES" = 8
)

matthew_plot_data <- matthew_bhb_summary %>%
  mutate(
    temperature_model = factor(
      temperature_model,
      levels = c("T(D47)CDES", "TD47iCDES+", "T(D47)iCDES")
    ),
    depth_plot_m = depth_locality_m + unname(model_offsets[as.character(temperature_model)])
  )

matthew_label_data <- matthew_plot_data %>%
  group_by(locality_id, depth_locality_m) %>%
  summarise(
    label_x = max(mean_temp_c + if_else(is.na(two_se_temp_c), 0, two_se_temp_c)),
    .groups = "drop"
  ) %>%
  mutate(
    label_x = label_x + diff(matthew_x_limits) * 0.04,
    locality_label = str_remove(locality_id, "-PK95$")
  )

matthew_plot <- ggplot(
  matthew_plot_data,
  aes(x = mean_temp_c, y = depth_plot_m, color = temperature_model)
) +
  geom_errorbar(
    aes(
      xmin = mean_temp_c - two_se_temp_c,
      xmax = mean_temp_c + two_se_temp_c
    ),
    orientation = "y",
    width = 0,
    linewidth = 0.5,
    na.rm = TRUE
  ) +
  geom_point(size = 2.6) +
  geom_text(
    data = matthew_label_data,
    aes(x = label_x, y = depth_locality_m, label = locality_label),
    inherit.aes = FALSE,
    color = "black",
    hjust = 0,
    size = 3,
    check_overlap = TRUE
  ) +
  scale_x_continuous(
    limits = matthew_x_limits,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_color_manual(
    values = c(
      "T(D47)CDES" = "#1b9e77",
      "TD47iCDES+" = "#d95f02",
      "T(D47)iCDES" = "#7570b3"
    ),
    name = "Workbook column"
  ) +
  labs(
    title = "Matthew's BHB Carbs AI:AK temperatures by locality",
    subtitle = "Workbook columns AI, AJ, and AK summarized by locality; error bars show +/- 2 SE",
    x = expression(paste("Mean clumped-isotope temperature (", degree, "C)")),
    y = "Locality stratigraphic height (m)"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.margin = margin(5.5, 90, 5.5, 5.5)
  )

ggsave(
  filename = matthew_plot_output,
  plot = matthew_plot,
  width = 9,
  height = 6.5,
  dpi = 600
)

x_lower <- min(
  ifelse(
    is.na(locality_temperature_summary$two_se_temp_c),
    locality_temperature_summary$mean_temp_c,
    locality_temperature_summary$mean_temp_c - locality_temperature_summary$two_se_temp_c
  )
)
x_upper <- max(
  ifelse(
    is.na(locality_temperature_summary$two_se_temp_c),
    locality_temperature_summary$mean_temp_c,
    locality_temperature_summary$mean_temp_c + locality_temperature_summary$two_se_temp_c
  )
)
x_padding <- max(2, (x_upper - x_lower) * 0.08)
x_limits <- c(x_lower - x_padding, x_upper + (2 * x_padding))

build_locality_plot(
  locality_temperature_summary,
  model_name = "Peterson",
  x_limits = x_limits,
  output_path = peterson_plot_output
)

build_locality_plot(
  locality_temperature_summary,
  model_name = "Anderson",
  x_limits = x_limits,
  output_path = anderson_plot_output
)

print(locality_temperature_summary, n = Inf)
message("Wrote summary to: ", summary_output)
message("Wrote sample error table to: ", sample_error_output)
message("Wrote Matthew workbook summary to: ", matthew_summary_output)
message("Wrote Matthew workbook plot to: ", matthew_plot_output)
