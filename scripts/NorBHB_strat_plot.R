
# NorBHB_strat_plot.R
# Creates one reusable ggplot object named `p_strat`.
# Source this script, then combine `p_strat` with other ggplot objects
# using patchwork, cowplot, or another layout package.

library(ggplot2)
library(dplyr)
library(readr)
library(here)
library(stringr)
library(tibble)

# -------------------------------------------------------------------
# 1. FILE PATHS
# -------------------------------------------------------------------

units_file <- here(
  "data",
  "raw",
  "NorBHB_strat_units.csv"
)

annotations_file <- here(
  "data",
  "raw",
  "NorBHB_strat_annotations.csv"
)

columns_file <- here(
  "data",
  "raw",
  "NorBHB_strat_columns.csv"
)

required_files <- c(
  units_file,
  annotations_file,
  columns_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required CSV file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

# -------------------------------------------------------------------
# 2. LOAD DATA
# -------------------------------------------------------------------

strat_units <- read_csv(
  units_file,
  show_col_types = FALSE,
  na = c("", "NA")
)

strat_annotations <- read_csv(
  annotations_file,
  show_col_types = FALSE,
  na = c("", "NA")
)

strat_columns <- read_csv(
  columns_file,
  show_col_types = FALSE,
  na = c("", "NA")
)

# -------------------------------------------------------------------
# DISPLAY OPTIONS
# -------------------------------------------------------------------

show_formation  <- TRUE
show_member     <- FALSE
show_chron      <- TRUE
show_lma        <- FALSE
show_biozone    <- TRUE

show_samples    <- TRUE
show_boundaries <- TRUE
show_headers    <- TRUE
show_grainsize  <- FALSE

# -------------------------------------------------------------------
# 3. VALIDATE INPUT TABLES
# -------------------------------------------------------------------

required_unit_cols <- c(
  "unit_id",
  "base_m",
  "top_m",
  "lithology",
  "grain_class",
  "unit_label",
  "notes"
)

required_annotation_cols <- c(
  "depth_m",
  "type",
  "label",
  "notes"
)

required_column_cols <- c(
  "base_m",
  "top_m",
  "label",
  "column"
)

missing_unit_cols <- setdiff(
  required_unit_cols,
  names(strat_units)
)

missing_annotation_cols <- setdiff(
  required_annotation_cols,
  names(strat_annotations)
)

missing_column_cols <- setdiff(
  required_column_cols,
  names(strat_columns)
)

if (length(missing_unit_cols) > 0) {
  stop(
    "NorBHB_strat_units.csv is missing column(s): ",
    paste(missing_unit_cols, collapse = ", ")
  )
}

if (length(missing_annotation_cols) > 0) {
  stop(
    "NorBHB_strat_annotations.csv is missing column(s): ",
    paste(missing_annotation_cols, collapse = ", ")
  )
}

if (length(missing_column_cols) > 0) {
  stop(
    "NorBHB_strat_columns.csv is missing column(s): ",
    paste(missing_column_cols, collapse = ", ")
  )
}

if (any(
  strat_units$top_m <= strat_units$base_m,
  na.rm = TRUE
)) {
  stop(
    "Every lithologic unit must have top_m greater than base_m."
  )
}

if (any(
  strat_columns$top_m <= strat_columns$base_m,
  na.rm = TRUE
)) {
  stop(
    "Every stratigraphic-column interval must have top_m greater than base_m."
  )
}

# -------------------------------------------------------------------
# 4. STANDARDIZE AND SORT DATA
# -------------------------------------------------------------------

strat_units <- strat_units %>%
  mutate(
    lithology = str_to_lower(str_trim(lithology)),
    grain_class = str_to_lower(str_trim(grain_class))
  ) %>%
  arrange(base_m)

strat_annotations <- strat_annotations %>%
  mutate(
    type = str_to_lower(str_trim(type)),
    label = str_trim(label)
  ) %>%
  arrange(depth_m)

strat_columns <- strat_columns %>%
  mutate(
    column = str_trim(column),
    label = str_trim(label)
  ) %>%
  arrange(column, base_m)

# -------------------------------------------------------------------
# 5. PLOT LIMITS
# -------------------------------------------------------------------

strat_y_limits <- c(
  min(
    strat_units$base_m,
    strat_columns$base_m,
    strat_annotations$depth_m,
    na.rm = TRUE
  ),
  max(
    strat_units$top_m,
    strat_columns$top_m,
    strat_annotations$depth_m,
    na.rm = TRUE
  )
)

strat_y_breaks <- seq(
  floor(strat_y_limits[1] / 50) * 50,
  ceiling(strat_y_limits[2] / 50) * 50,
  by = 50
)

# -------------------------------------------------------------------
# 6. LITHOLOGY AND GRAIN-SIZE CONTROLS
# -------------------------------------------------------------------

grain_widths <- c(
  mud = 0.18,
  vf = 0.32,
  fine = 0.46,
  medium = 0.62,
  coarse = 0.78,
  coal = 0.24,
  ash = 0.28,
  covered = 0.42
)

lithology_fills <- c(
  mudstone = "grey82",
  sandstone = "grey97",
  coal = "grey15",
  ash = "grey60",
  covered = "white"
)

strat_units_plot <- strat_units %>%
  mutate(
    xmin = 0,
    xmax = unname(grain_widths[grain_class]),
    xmax = coalesce(xmax, 0.18)
  )

# -------------------------------------------------------------------
# 7. PREPARE ANNOTATIONS
# -------------------------------------------------------------------

strat_annotations_plot <- strat_annotations %>%
  mutate(
    plot_label = case_when(
      is.na(notes) | notes == "" ~ label,
      TRUE ~ paste0(label, " — ", notes)
    )
  )

strat_samples_plot <- strat_annotations_plot %>%
  filter(
    type %in% c(
      "sample",
      "fossil",
      "flora",
      "quarry"
    )
  )

strat_boundaries_plot <- strat_annotations_plot %>%
  filter(
    type %in% c(
      "ash",
      "biozone",
      "boundary",
      "unconformity",
      "marker bed",
      "marker_bed"
    )
  )

# -------------------------------------------------------------------
# 8. PREPARE LEFT-HAND STRATIGRAPHIC COLUMNS
# -------------------------------------------------------------------

# Edit this order to control the left-to-right arrangement.
column_order <- c(
  "Formation",
  "Member",
  "Chron",
  "LMA",
  "Biozone"
)

column_order <- c()

if (show_formation) column_order <- c(column_order, "Formation")
if (show_member)    column_order <- c(column_order, "Member")
if (show_chron)     column_order <- c(column_order, "Chron")
if (show_lma)       column_order <- c(column_order, "LMA")
if (show_biozone)   column_order <- c(column_order, "Biozone")

# Any additional column names not listed above are appended automatically.
extra_columns <- setdiff(
  unique(strat_columns$column),
  column_order
)

column_order <- c(
  column_order,
  extra_columns
)

column_width <- 0.30
column_gap <- 0.02
lithology_gap <- 0.10

column_positions <- tibble(
  column = column_order,
  column_number = seq_along(column_order)
) %>%
  mutate(
    xmax = -lithology_gap -
      (length(column_order) - column_number) *
      (column_width + column_gap),
    xmin = xmax - column_width,
    xmid = (xmin + xmax) / 2
  )

strat_columns_plot <- strat_columns %>%
  mutate(
    column = factor(
      column,
      levels = column_order
    )
  ) %>%
  left_join(
    column_positions,
    by = "column"
  ) %>%
  mutate(
    ymid = (base_m + top_m) / 2
  )

if (any(is.na(strat_columns_plot$xmin))) {
  stop(
    "At least one value in strat_columns$column could not be assigned ",
    "a plotting position."
  )
}

strat_x_min <- min(
  strat_columns_plot$xmin,
  na.rm = TRUE
)

strat_x_max <- 2.15

# -------------------------------------------------------------------
# 9. GRAIN-SIZE HEADER
# -------------------------------------------------------------------

grain_headers <- tibble(
  x = c(
    0.03,
    grain_widths["vf"],
    grain_widths["fine"],
    grain_widths["medium"],
    grain_widths["coarse"]
  ),
  label = c(
    "Mud",
    "vf",
    "f",
    "m",
    "c"
  )
)

header_bottom <- strat_y_limits[2] + 4
header_top <- strat_y_limits[2] + 18
header_mid <- (header_bottom + header_top) / 2

# -------------------------------------------------------------------
# 10. REUSABLE GGPLOT OBJECT
# -------------------------------------------------------------------

p_strat <- ggplot() +

  # Left-hand formation/member/chron/LMA/biozone columns
  geom_rect(
    data = strat_columns_plot,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = base_m,
      ymax = top_m
    ),
    fill = "white",
    colour = "black",
    linewidth = 0.3
  ) +

  geom_text(
    data = strat_columns_plot,
    aes(
      x = xmid,
      y = ymid,
      label = label
    ),
    angle = 90,
    size = 2.45,
    fontface = "bold"
  ) +

  # Lithologic bodies
  geom_rect(
    data = strat_units_plot,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = base_m,
      ymax = top_m,
      fill = lithology
    ),
    colour = "black",
    linewidth = 0.25
  ) +

  # Covered intervals
  geom_rect(
    data = filter(
      strat_units_plot,
      lithology == "covered"
    ),
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = base_m,
      ymax = top_m
    ),
    fill = NA,
    colour = "grey35",
    linewidth = 0.35,
    linetype = "dashed"
  ) +

  # Samples, fossils, flora, and quarries
  geom_segment(
    data = strat_samples_plot,
    aes(
      x = 0.84,
      xend = 1.00,
      y = depth_m,
      yend = depth_m
    ),
    linewidth = 0.35
  ) +

  geom_point(
    data = strat_samples_plot,
    aes(
      x = 1.00,
      y = depth_m,
      shape = type
    ),
    size = 1.25
  ) +

  geom_text(
    data = strat_samples_plot,
    aes(
      x = 1.06,
      y = depth_m,
      label = plot_label
    ),
    hjust = 0,
    size = 2.35,
    check_overlap = TRUE
  ) +

  # Marker beds and boundaries
  geom_hline(
    data = strat_boundaries_plot,
    aes(
      yintercept = depth_m,
      linetype = type
    ),
    linewidth = 0.35
  ) +

  geom_text(
    data = strat_boundaries_plot,
    aes(
      x = 1.06,
      y = depth_m,
      label = label
    ),
    hjust = 0,
    vjust = -0.35,
    size = 2.25,
    check_overlap = TRUE
  ) +

  # Left-column headings
  geom_text(
    data = column_positions,
    aes(
      x = xmid,
      y = strat_y_limits[2],
      label = column
    ),
    angle = 90,
    hjust = -0.08,
    vjust = 0.5,
    size = 2.7,
    fontface = "bold",
    inherit.aes = FALSE
  ) +

  # Meters heading immediately left of lithology
  annotate(
    "text",
    x = -0.04,
    y = strat_y_limits[2],
    label = "Meters",
    angle = 90,
    hjust = -0.08,
    vjust = 0.5,
    size = 2.7,
    fontface = "bold"
  ) +

  # Grain-size header boxes
  geom_rect(
    data = grain_headers,
    aes(
      xmin = x - 0.075,
      xmax = x + 0.075,
      ymin = header_bottom,
      ymax = header_top
    ),
    fill = "white",
    colour = "black",
    linewidth = 0.25,
    inherit.aes = FALSE
  ) +

  geom_text(
    data = grain_headers,
    aes(
      x = x,
      y = header_mid,
      label = label
    ),
    size = 2.1,
    inherit.aes = FALSE
  ) +

  scale_fill_manual(
    values = lithology_fills,
    drop = FALSE
  ) +

  scale_shape_manual(
    values = c(
      sample = 16,
      fossil = 17,
      flora = 15,
      quarry = 18
    )
  ) +

  scale_linetype_manual(
    values = c(
      ash = "solid",
      biozone = "dashed",
      boundary = "dashed",
      unconformity = "longdash",
      `marker bed` = "solid",
      marker_bed = "solid"
    )
  ) +

  scale_y_continuous(
    limits = strat_y_limits,
    breaks = strat_y_breaks,
    expand = expansion(mult = c(0, 0))
  ) +

  scale_x_continuous(
    limits = c(
      strat_x_min - 0.08,
      strat_x_max
    ),
    expand = expansion(mult = c(0, 0))
  ) +

  labs(
    x = NULL,
    y = NULL
  ) +

  coord_cartesian(
    clip = "off"
  ) +

  theme_classic(
    base_size = 9
  ) +

  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),

    axis.title.y = element_blank(),

    legend.position = "none",

    plot.margin = margin(
      t = 50,
      r = 110,
      b = 5.5,
      l = 5.5
    )
  )

# Display the stratigraphic panel
p_strat

# -------------------------------------------------------------------
# DOWNSTREAM EXAMPLE
#
# library(patchwork)
#
# final_plot <-
#   p_strat +
#   p_d13C +
#   p_temperature +
#   plot_layout(
#     widths = c(1.9, 1, 1)
#   )
#
# Adjacent panels should use:
#
# scale_y_continuous(
#   limits = strat_y_limits,
#   breaks = strat_y_breaks,
#   expand = expansion(mult = c(0, 0))
# )
# -------------------------------------------------------------------
