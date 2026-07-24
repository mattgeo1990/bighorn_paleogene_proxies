# save_figure_variants.R
# Shared export utilities for production figures.
#
# Every production figure is written in two forms:
#   manuscript/   publication layout at an explicitly supplied size
#   presentation/ slide/web layout, 6 x 6 inches by default or 12 x 6 when wide
#
# Each form is exported as PNG, PDF, and editable SVG. Existing legacy exports
# may remain during the transition, but downstream manuscripts and talks should
# use the files in these explicit variant folders.

prepare_figure_variant_dirs <- function(base_dir) {
  dirs <- list(
    manuscript = file.path(base_dir, "manuscript"),
    presentation = file.path(base_dir, "presentation")
  )
  purrr::walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  dirs
}

add_presentation_theme <- function(plot) {
  slide_theme <- ggplot2::theme(
    axis.title = ggplot2::element_text(size = 13),
    axis.text = ggplot2::element_text(size = 11, color = "black"),
    legend.title = ggplot2::element_text(size = 11),
    legend.text = ggplot2::element_text(size = 10),
    plot.title = ggplot2::element_text(size = 14, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 11),
    plot.caption = ggplot2::element_text(size = 8, hjust = 0),
    plot.margin = ggplot2::margin(5, 7, 5, 7)
  )

  if (inherits(plot, "patchwork")) {
    plot & slide_theme
  } else {
    plot + slide_theme
  }
}

export_plot_triplet <- function(
    plot, directory, stem, width, height, dpi = 500
) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    file.path(directory, paste0(stem, ".png")),
    plot, width = width, height = height, dpi = dpi, bg = "white"
  )
  ggplot2::ggsave(
    file.path(directory, paste0(stem, ".pdf")),
    plot, width = width, height = height,
    device = grDevices::pdf, useDingbats = FALSE, bg = "white"
  )
  ggplot2::ggsave(
    file.path(directory, paste0(stem, ".svg")),
    plot, width = width, height = height,
    device = svglite::svglite, bg = "white"
  )
}

save_figure_variants <- function(
    plot,
    base_dir,
    stem,
    manuscript_width,
    manuscript_height,
    presentation_plot = NULL,
    presentation_width = NULL,
    presentation_height = 6,
    manuscript_dpi = 500,
    presentation_dpi = 400
) {
  dirs <- prepare_figure_variant_dirs(base_dir)

  if (is.null(presentation_width)) {
    aspect <- manuscript_width / manuscript_height
    presentation_width <- if (aspect > 1.35) 12 else 6
  }

  if (is.null(presentation_plot)) {
    presentation_plot <- add_presentation_theme(plot)
  }

  export_plot_triplet(
    plot, dirs$manuscript, stem,
    manuscript_width, manuscript_height, manuscript_dpi
  )
  export_plot_triplet(
    presentation_plot, dirs$presentation, stem,
    presentation_width, presentation_height, presentation_dpi
  )

  invisible(list(
    manuscript_dir = dirs$manuscript,
    presentation_dir = dirs$presentation
  ))
}
