# CFB_chronostrat_panels.R
# Reusable chronostratigraphic side panels for CFB figures.
#
# The authoritative interval definitions remain in
# data/raw/NorBHB_strat_columns.csv. The age-domain panel transforms those
# native meter boundaries through the same CFB age model used for samples.

library(tidyverse)
library(here)

read_CFB_chronostrat_columns <- function() {
  read_csv(
    here("data", "raw", "NorBHB_strat_columns.csv"),
    show_col_types = FALSE,
    na = c("", "NA")
  ) %>%
    mutate(
      column = str_trim(column),
      label = str_trim(label)
    ) %>%
    filter(column %in% c("Epoch", "Formation", "Chron", "Biozone"))
}

chronostrat_fill <- function(column, label) {
  case_when(
    column == "Chron" & str_detect(label, "n$") ~ "#202020",
    column == "Chron" ~ "white",
    column == "Epoch" & label == "Paleocene" ~ "#F6E8A6",
    column == "Epoch" & label == "Eocene" ~ "#F7B267",
    column == "Formation" ~ "#D9D9D9",
    column == "Biozone" ~ "#F7F7F7",
    TRUE ~ "white"
  )
}

chronostrat_text_color <- function(column, label) {
  if_else(column == "Chron" & str_detect(label, "n$"), "white", "black")
}

prepare_CFB_chronostrat_strat <- function(strat_limits) {
  read_CFB_chronostrat_columns() %>%
    mutate(
      ymin = pmax(base_m, min(strat_limits)),
      ymax = pmin(top_m, max(strat_limits))
    ) %>%
    filter(ymax > ymin) %>%
    mutate(
      ymid = (ymin + ymax) / 2,
      fill_value = chronostrat_fill(column, label),
      text_color = chronostrat_text_color(column, label),
      column = factor(
        column,
        levels = c("Epoch", "Formation", "Chron", "Biozone")
      ),
      x = as.numeric(column)
    )
}

prepare_CFB_chronostrat_age <- function(age_lookup, age_limits) {
  lookup <- age_lookup %>%
    filter(is.finite(strat_height_m), is.finite(Age_Ma)) %>%
    distinct(strat_height_m, Age_Ma) %>%
    arrange(strat_height_m)

  if (nrow(lookup) < 2) {
    stop("At least two CFB meter-age pairs are required.")
  }

  meter_to_age <- function(x) {
    approx(
      x = lookup$strat_height_m,
      y = lookup$Age_Ma,
      xout = x,
      rule = 2,
      ties = mean
    )$y
  }

  read_CFB_chronostrat_columns() %>%
    mutate(
      age_at_base_ma = meter_to_age(base_m),
      age_at_top_ma = meter_to_age(top_m),
      ymin = pmax(
        pmin(age_at_base_ma, age_at_top_ma),
        min(age_limits)
      ),
      ymax = pmin(
        pmax(age_at_base_ma, age_at_top_ma),
        max(age_limits)
      )
    ) %>%
    filter(ymax > ymin) %>%
    mutate(
      ymid = (ymin + ymax) / 2,
      fill_value = chronostrat_fill(column, label),
      text_color = chronostrat_text_color(column, label),
      column = factor(
        column,
        levels = c("Epoch", "Formation", "Chron", "Biozone")
      ),
      x = as.numeric(column)
    )
}

build_chronostrat_panel <- function(intervals, y_scale, title = NULL) {
  column_levels <- levels(intervals$column)
  column_labels <- recode(
    column_levels,
    "Biozone" = "Mammal\nbiozone",
    .default = column_levels
  )

  ggplot(intervals) +
    geom_rect(
      aes(
        xmin = x - 0.5, xmax = x + 0.5,
        ymin = ymin, ymax = ymax,
        fill = fill_value
      ),
      color = "black", linewidth = 0.28
    ) +
    geom_text(
      aes(
        x = x, y = ymid, label = label, color = text_color,
        angle = if_else(column == "Formation", 90, 0)
      ),
      size = 2.25, lineheight = 0.88
    ) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_x_continuous(
      breaks = seq_along(column_levels),
      labels = column_labels,
      position = "top",
      expand = c(0, 0)
    ) +
    y_scale +
    labs(x = NULL, y = NULL, title = title) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 9) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      axis.text.x.top = element_text(size = 7, face = "bold"),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 9, face = "bold"),
      plot.margin = margin(5, 1, 5, 2)
    )
}

build_CFB_chronostrat_strat_panel <- function(
    strat_limits,
    breaks = seq(min(strat_limits), max(strat_limits), by = 200)) {
  intervals <- prepare_CFB_chronostrat_strat(strat_limits)
  build_chronostrat_panel(
    intervals,
    scale_y_continuous(
      limits = strat_limits,
      breaks = breaks,
      expand = expansion(mult = c(0.01, 0.02))
    )
  )
}

build_CFB_chronostrat_age_panel <- function(
    age_lookup,
    age_limits,
    breaks = seq(
      ceiling(min(age_limits) * 2) / 2,
      floor(max(age_limits) * 2) / 2,
      by = 0.5
    )) {
  intervals <- prepare_CFB_chronostrat_age(age_lookup, age_limits)
  build_chronostrat_panel(
    intervals,
    scale_y_reverse(
      limits = age_limits,
      breaks = breaks,
      expand = expansion(mult = c(0.01, 0.02))
    )
  )
}
