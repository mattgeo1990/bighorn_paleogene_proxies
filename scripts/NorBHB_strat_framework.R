
# NorBHB_strat_framework.R
#
# Creates a compact stratigraphic framework plot named `p_strat_framework`
# using:
#   - NorBHB_strat_columns.csv
#   - NorBHB_strat_annotations.csv
#
# Columns displayed:
#   Height | Formation | Chron | Mammal biozone | Sample strip
#
# Intended final width: ~3.5 inches.
# Arial is used throughout.
#
# The sample strip has no jitter. All samples occur at the same x position.
# Vertical uncertainty bars are drawn when uncertainty information is present.

library(ggplot2)
library(dplyr)
library(readr)
library(here)
library(stringr)
library(tibble)

# -------------------------------------------------------------------
# 1. FILE PATHS
# -------------------------------------------------------------------

columns_file <- here(
  "data",
  "raw",
  "NorBHB_strat_columns.csv"
)

annotations_file <- here(
  "data",
  "raw",
  "NorBHB_strat_annotations.csv"
)

required_files <- c(
  columns_file,
  annotations_file
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

strat_columns <- read_csv(
  columns_file,
  show_col_types = FALSE,
  na = c("", "NA")
)

strat_annotations <- read_csv(
  annotations_file,
  show_col_types = FALSE,
  na = c("", "NA")
)

# -------------------------------------------------------------------
# 3. VALIDATE INPUT TABLES
# -------------------------------------------------------------------

required_column_cols <- c(
  "base_m",
  "top_m",
  "label",
  "column"
)

required_annotation_cols <- c(
  "depth_m",
  "type",
  "label"
)

missing_column_cols <- setdiff(
  required_column_cols,
  names(strat_columns)
)

missing_annotation_cols <- setdiff(
  required_annotation_cols,
  names(strat_annotations)
)

if (length(missing_column_cols) > 0) {
  stop(
    "NorBHB_strat_columns.csv is missing column(s): ",
    paste(missing_column_cols, collapse = ", ")
  )
}

if (length(missing_annotation_cols) > 0) {
  stop(
    "NorBHB_strat_annotations.csv is missing column(s): ",
    paste(missing_annotation_cols, collapse = ", ")
  )
}

if (any(
  strat_columns$top_m <= strat_columns$base_m,
  na.rm = TRUE
)) {
  stop(
    "Every interval in NorBHB_strat_columns.csv must have ",
    "top_m greater than base_m."
  )
}

# -------------------------------------------------------------------
# 4. DISPLAY CONTROLS
# -------------------------------------------------------------------

font_family <- "Arial"

# Final export dimensions
figure_width_in  <- 1.5
figure_height_in <- 3.5

# Height scale
height_break_m <- 200

# Show annotation labels beside sample points?
show_annotation_labels <- FALSE

# Which annotation types should appear in the sample strip?
sample_types <- c(
  "sample",
  "carbonate",
  "pedogenic carbonate",
  "enamel",
  "mammal enamel",
  "gar",
  "gar scale",
  "gar scales",
  "fossil",
  "flora"
)

# Stratigraphic columns to retain, in left-to-right order
column_order <- c(
    "Epoch",
"Formation",
  "Chron",
  "Biozone"
)

# Header labels shown above the columns
column_headers <- c(
  Formation = "Formation",
  Epoch     = "Epoch",
  Chron     = "Chron",
  Biozone   = "Mammal\nbiozone"
)

# Column widths in arbitrary plotting units
column_widths <- c(
  Formation = 0.4,
  Epoch     = 0.4,
  Chron     = 0.4,
  Biozone   = 0.5
)

column_gap <- 0.00
sample_strip_gap <- 0.12
sample_strip_width <- 0.34

# -------------------------------------------------------------------
# 5. STANDARDIZE STRATIGRAPHIC COLUMNS
# -------------------------------------------------------------------

strat_columns <- strat_columns %>%
  mutate(
    column = str_trim(column),
    label = str_trim(label)
  ) %>%
  filter(column %in% column_order) %>%
  arrange(
    factor(column, levels = column_order),
    base_m
  )

missing_requested_columns <- setdiff(
  column_order,
  unique(strat_columns$column)
)

if (length(missing_requested_columns) > 0) {
  stop(
    "The following requested stratigraphic column(s) are absent from ",
    "NorBHB_strat_columns.csv: ",
    paste(missing_requested_columns, collapse = ", ")
  )
}

# Calculate horizontal positions for each framework column
column_positions <- tibble(
  column = column_order,
  width = unname(column_widths[column_order])
) %>%
  mutate(
    xmin = lag(cumsum(width + column_gap), default = 0),
    xmax = xmin + width,
    xmid = (xmin + xmax) / 2
  )

strat_columns_plot <- strat_columns %>%
  left_join(
    column_positions,
    by = "column"
  ) %>%
  mutate(
    ymid = (base_m + top_m) / 2,
    label_angle = if_else(
      column %in% c("Formation", "Chron"),
      90,
      0
    )
  )

# -------------------------------------------------------------------
# 6. PREPARE ANNOTATIONS AND UNCERTAINTIES
# -------------------------------------------------------------------

strat_annotations <- strat_annotations %>%
  mutate(
    type = str_to_lower(str_trim(type)),
    label = str_trim(label)
  ) %>%
  filter(type %in% sample_types) %>%
  arrange(depth_m)

# Support several common uncertainty-column conventions.
#
# Preferred:
#   lower_m and upper_m
#
# Also supported:
#   depth_min_m and depth_max_m
#   uncertainty_m (symmetric ± uncertainty)
#   strat_uncertainty_m (symmetric ± uncertainty)

if (all(c("lower_m", "upper_m") %in% names(strat_annotations))) {

  strat_annotations <- strat_annotations %>%
    mutate(
      depth_lower_m = lower_m,
      depth_upper_m = upper_m
    )

} else if (all(
  c("depth_min_m", "depth_max_m") %in% names(strat_annotations)
)) {

  strat_annotations <- strat_annotations %>%
    mutate(
      depth_lower_m = depth_min_m,
      depth_upper_m = depth_max_m
    )

} else if ("uncertainty_m" %in% names(strat_annotations)) {

  strat_annotations <- strat_annotations %>%
    mutate(
      depth_lower_m = depth_m - uncertainty_m,
      depth_upper_m = depth_m + uncertainty_m
    )

} else if ("strat_uncertainty_m" %in% names(strat_annotations)) {

  strat_annotations <- strat_annotations %>%
    mutate(
      depth_lower_m = depth_m - strat_uncertainty_m,
      depth_upper_m = depth_m + strat_uncertainty_m
    )

} else {

  message(
    "No recognized uncertainty columns were found. ",
    "Points will be plotted without vertical uncertainty bars."
  )

  strat_annotations <- strat_annotations %>%
    mutate(
      depth_lower_m = depth_m,
      depth_upper_m = depth_m
    )
}

# -------------------------------------------------------------------
# 7. PLOT LIMITS AND SAMPLE-STRIP POSITIONS
# -------------------------------------------------------------------

strat_y_limits <- c(
  min(
    strat_columns_plot$base_m,
    strat_annotations$depth_lower_m,
    na.rm = TRUE
  ),
  max(
    strat_columns_plot$top_m,
    strat_annotations$depth_upper_m,
    na.rm = TRUE
  )
)

strat_y_breaks <- seq(
  floor(strat_y_limits[1] / height_break_m) * height_break_m,
  ceiling(strat_y_limits[2] / height_break_m) * height_break_m,
  by = height_break_m
)

framework_right <- max(column_positions$xmax)

sample_strip_xmin <- framework_right + sample_strip_gap
sample_strip_xmax <- sample_strip_xmin + sample_strip_width
sample_x <- (sample_strip_xmin + sample_strip_xmax) / 2

strat_annotations_plot <- strat_annotations %>%
  mutate(
    sample_x = sample_x
  )

# -------------------------------------------------------------------
# 8. REUSABLE GGPLOT OBJECT
# -------------------------------------------------------------------
# Text-size controls, in points
text_sizes <- list(
  formation_label = 6.8,
  chron_label = 6.5,
  biozone_label = 6.3,
  axis_tick = 6.5
)

# Pull the Fort Union–Willwood boundary from the Formation column
formation_boundary_m <- strat_columns_plot %>%
  filter(
    column == "Formation",
    str_to_lower(label) == "willwood"
  ) %>%
  pull(base_m)

if (length(formation_boundary_m) != 1) {
  stop(
    "Could not uniquely identify the base of the Willwood Formation."
  )
}

formation_xmin <- column_positions %>%
  filter(column == "Formation") %>%
  pull(xmin)

formation_xmax <- column_positions %>%
  filter(column == "Formation") %>%
  pull(xmax)

p_strat_framework <- ggplot() +
  
  # Formation, Epoch, Chron, and Biozone boxes
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
    linewidth = 0.28
  ) +
  
  # Formation labels
  geom_text(
    data = strat_columns_plot %>%
      filter(
        column == "Formation",
        !is.na(label),
        label != "",
        label != "(blank)"
      ),
    aes(
      x = xmid,
      y = ymid,
      label = label
    ),
    family = font_family,
    size = text_sizes$formation_label / ggplot2::.pt,
    angle = 90,
    lineheight = 0.92
  ) +
  
  # Epoch labels
  geom_text(
    data = strat_columns_plot %>%
      filter(
        column == "Epoch",
        !is.na(label),
        label != "",
        label != "(blank)"
      ),
    aes(
      x = xmid,
      y = ymid,
      label = label
    ),
    family = font_family,
    size = text_sizes$formation_label / ggplot2::.pt,
    angle = 90,
    lineheight = 0.92
  ) +
  
  # Chron labels
  geom_text(
    data = strat_columns_plot %>%
      filter(
        column == "Chron",
        !is.na(label),
        label != "",
        label != "(blank)"
      ),
    aes(
      x = xmid,
      y = ymid,
      label = label
    ),
    family = font_family,
    size = text_sizes$chron_label / ggplot2::.pt,
    angle = 90,
    lineheight = 0.92
  ) +
  
  # Mammal-biozone labels
  geom_text(
    data = strat_columns_plot %>%
      filter(
        column == "Biozone",
        !is.na(label),
        label != "",
        label != "(blank)"
      ),
    aes(
      x = xmid,
      y = ymid,
      label = label
    ),
    family = font_family,
    size = text_sizes$biozone_label / ggplot2::.pt,
    angle = 0,
    lineheight = 0.92
  ) +
  
  
  scale_y_continuous(
    limits = strat_y_limits,
    breaks = strat_y_breaks,
    expand = expansion(mult = c(0, 0)),
    position = "left"
  ) +
  
  scale_x_continuous(
    limits = c(
      -0.15,
      max(column_positions$xmax)
    ),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  coord_cartesian(
    ylim = strat_y_limits,
    clip = "on",
    expand = FALSE
  ) +
  
  theme_classic(
    base_family = font_family,
    base_size = text_sizes$axis_tick
  ) +
  
  theme(
    axis.title = element_blank(),
    
    axis.text.y = element_text(
      family = font_family,
      size = text_sizes$axis_tick,
      colour = "black",
      margin = margin(r = 2)
    ),
    
    axis.ticks.y = element_line(
      linewidth = 0.25
    ),
    
    axis.ticks.length.y = grid::unit(
      1.2,
      "mm"
    ),
    
    axis.line.y = element_line(
      linewidth = 0.28
    ),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    
    legend.position = "none",
    
    plot.margin = margin(
      t = 0,
      r = 3,
      b = 3,
      l = 2,
      unit = "pt"
    )
  )

p_strat_framework

# -------------------------------------------------------------------
# 9. OPTIONAL EXPORT
# -------------------------------------------------------------------


ggsave(
  filename = here(
    "figures",
    "NorBHB_strat_framework.png"
  ),
  plot = p_strat_framework,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = here(
    "figures",
    "NorBHB_strat_framework.svg"
  ),
  plot = p_strat_framework,
  width = 1.5,
  height = 3,
  units = "in",
  device = svglite::svglite
)

# -------------------------------------------------------------------
# DOWNSTREAM EXAMPLE
#
# library(patchwork)
#
# final_figure <-
#   p_map +
#   p_strat_framework +
#   p_geochem +
#   plot_layout(
#     widths = c(1.2, 0.8, 1.4)
#   )
# -------------------------------------------------------------------
