
# polecat_strat_plot.R
# Creates one reusable ggplot object named `p_strat`.
# Source this script, then add `p_strat` to a patchwork/cowplot multipanel figure.

library(ggplot2)
library(dplyr)
library(tibble)
library(readr)
library(here)
# -------------------------------------------------------------------
# 1. LOAD STRATIGRAPHIC DATA FROM CSV FILES
# -------------------------------------------------------------------

library(readr)
library(dplyr)
library(here)

# Define file paths
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

# Check that all files exist before loading
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

# Load data
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
# 2. VALIDATE COLUMN STRUCTURE
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

# -------------------------------------------------------------------
# 3. VALIDATE STRATIGRAPHIC VALUES
# -------------------------------------------------------------------

if (any(
  strat_units$top_m <= strat_units$base_m,
  na.rm = TRUE
)) {
  stop(
    "Every stratigraphic unit must have top_m greater than base_m."
  )
}

if (any(
  strat_columns$top_m <= strat_columns$base_m,
  na.rm = TRUE
)) {
  stop(
    "Every stratigraphic column interval must have top_m greater than base_m."
  )
}

# -------------------------------------------------------------------
# 4. SORT DATA
# -------------------------------------------------------------------

strat_units <- strat_units %>%
  arrange(base_m)

strat_annotations <- strat_annotations %>%
  arrange(depth_m)

strat_columns <- strat_columns %>%
  arrange(column, base_m)
# -------------------------------------------------------------------
# 2. PREPARE DATA AND PLOT CONTROLS
# -------------------------------------------------------------------

# Plot limits based on all three stratigraphic tables
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

# Width assigned to each grain-size class
grain_widths <- c(
  mud     = 0.30,
  fine    = 0.58,
  medium  = 0.75,
  coarse  = 0.92,
  coal    = 0.38,
  ash     = 0.44,
  covered = 0.52
)

# Fill assigned to each lithology
lithology_fills <- c(
  mudstone  = "grey82",
  sandstone = "grey97",
  coal      = "grey15",
  ash       = "grey60",
  covered   = "white"
)

# Prepare lithologic intervals
strat_units_plot <- strat_units %>%
  mutate(
    lithology = stringr::str_to_lower(lithology),
    grain_class = stringr::str_to_lower(grain_class),
    xmin = 0,
    xmax = unname(grain_widths[grain_class]),
    xmax = dplyr::coalesce(xmax, 0.30)
  )

# Standardize annotation types
strat_annotations_plot <- strat_annotations %>%
  mutate(
    type = stringr::str_to_lower(type),
    plot_label = case_when(
      is.na(notes) | notes == "" ~ label,
      TRUE ~ paste0(label, " — ", notes)
    )
  )

# Point annotations
strat_samples_plot <- strat_annotations_plot %>%
  filter(type %in% c(
    "sample",
    "fossil",
    "flora",
    "quarry"
  ))

# Horizontal boundaries
strat_boundaries_plot <- strat_annotations_plot %>%
  filter(type %in% c(
    "ash",
    "biozone",
    "boundary",
    "unconformity",
    "marker bed",
    "marker_bed"
  ))

# -------------------------------------------------------------------
# PREPARE FORMATION / MEMBER / BIOZONE COLUMNS
# -------------------------------------------------------------------

# The strat_columns CSV should contain:
# xmin, xmax, base_m, top_m, label, column

required_column_plot_cols <- c(
  "xmin",
  "xmax",
  "base_m",
  "top_m",
  "label",
  "column"
)

missing_column_plot_cols <- setdiff(
  required_column_plot_cols,
  names(strat_columns)
)

if (length(missing_column_plot_cols) > 0) {
  stop(
    "NorBHB_strat_columns.csv is missing column(s): ",
    paste(missing_column_plot_cols, collapse = ", ")
  )
}

strat_columns_plot <- strat_columns %>%
  mutate(
    xmid = (xmin + xmax) / 2,
    ymid = (base_m + top_m) / 2
  ) %>%
  arrange(column, base_m)

# One heading position for each side column
column_positions <- strat_columns_plot %>%
  group_by(column) %>%
  summarise(
    xmin = min(xmin, na.rm = TRUE),
    xmax = max(xmax, na.rm = TRUE),
    xmid = (xmin + xmax) / 2,
    .groups = "drop"
  )

# Expand the plot far enough to show all side columns
strat_x_max <- max(
  c(2.15, strat_columns_plot$xmax),
  na.rm = TRUE
)

# -------------------------------------------------------------------
# 3. REUSABLE GGPLOT OBJECT
# -------------------------------------------------------------------

p_strat <- ggplot() +
  
  # ---------------------------------------------------------------
# Lithologic bodies
# ---------------------------------------------------------------

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
  
  # Dashed outline for covered intervals
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
  
  # ---------------------------------------------------------------
# Samples, fossils, and flora
# ---------------------------------------------------------------

geom_segment(
  data = strat_samples_plot,
  aes(
    x = 0.88,
    xend = 1.02,
    y = depth_m,
    yend = depth_m
  ),
  linewidth = 0.35
) +
  
  geom_point(
    data = strat_samples_plot,
    aes(
      x = 1.02,
      y = depth_m,
      shape = type
    ),
    size = 1.25
  ) +
  
  geom_text(
    data = strat_samples_plot,
    aes(
      x = 1.08,
      y = depth_m,
      label = plot_label
    ),
    hjust = 0,
    size = 2.35,
    check_overlap = TRUE
  ) +
  
  # ---------------------------------------------------------------
# Marker beds and boundaries
# ---------------------------------------------------------------

geom_hline(
  data = strat_boundaries_plot,
  aes(
    yintercept = depth_m,
    linetype = type
  ),
  linewidth = 0.35
) +
  
  # Optional labels beside important boundaries
  geom_text(
    data = strat_boundaries_plot,
    aes(
      x = 1.08,
      y = depth_m,
      label = label
    ),
    hjust = 0,
    vjust = -0.35,
    size = 2.25,
    check_overlap = TRUE
  ) +
  
  # ---------------------------------------------------------------
# Formation, member, zone, or other interval columns
# ---------------------------------------------------------------

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
      y = (base_m + top_m) / 2,
      label = label
    ),
    angle = 90,
    size = 2.4
  ) +
  
  # Column headings
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
    size = 2.5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  
  # ---------------------------------------------------------------
# Scales
# ---------------------------------------------------------------

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
    limits = c(-0.08, strat_x_max + 0.08),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    x = NULL,
    y = "Stratigraphic height (m)"
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
    
    legend.position = "none",
    
    plot.margin = margin(
      t = 35,
      r = 10,
      b = 5.5,
      l = 5.5
    )
  )

# Display the stratigraphic panel
p_strat
