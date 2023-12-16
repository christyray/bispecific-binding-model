
# LOAD PACKAGES -----------------------------------------------------------

library(tidyverse)
library(cli)
library(scales)
library(sciscales)
library(ggh4x)
library(crthemes)
library(paletteer)
library(patchwork)

# FUNCTIONS ---------------------------------------------------------------

# Plot Function -----------------------------------------------------------

# Generate the main text and supplemental figures for simulations with varying
# antibody and receptor concentrations
vary_concentrations <- function(data_file = NULL, fig_file = NULL,
                          fig_path = NULL, poster = FALSE,
                          suppress_warn = TRUE, dpi = 600) {

  # Set display name for current figure
  name <- "Vary Concentrations"

  # Exit with message if required input not provided
  if (any(c(is.null(data_file), is.null(fig_file), is.null(fig_path)))) {
    cli_alert_info("Skipping {.field {name}}")
    return(NULL)
  }

  # Helper function to create full figure file path from the specific name of
  # the desired plot
  # Since `fig_path` is not a function input, the function will determine its
  # values from the variables in the global environment (i.e., the workspace)
  figpath <- function(fig_file, suffix = NULL) {
    if (!is.null(suffix)) {
      fullfile <- paste(fig_file, suffix, sep = "-")
    } else {
      fullfile <- fig_file
    }
    here::here(paste0(fig_path, fullfile, ".png"))
  }

  # Suppress warnings when saving plots if requested
  if (suppress_warn) {
    savefig <- function(...) {
      suppressWarnings(crsave(dpi = dpi, ...))
    }
  } else {
    savefig <- function(...) {
      crsave(dpi = dpi, ...)
    }
  }

  # Print update message
  cli_progress_step(name)

  # DATA IMPORT -------------------------------------------------------------

  # Varying antibody and receptor concentrations for different antibodies
  data <- import_data(data_file, type = "sim")

  # Convert the antibody and receptor concentrations to log concentrations, and
  # filter the data to only the necessary values for plotting
  data <- data$occupied |>
    log_conc() |>
    filter(
      Antibody %in% c("BS1", "Toci_H2"),
      Occupied %in% c("Antibody", "Receptor"),
      Species %in% c("Free", "Binary", "Ternary", "Total")
    )

  # PLOTTING ----------------------------------------------------------------

  # Heat map with varying total antibody and receptor concentrations
  vary_heatmap <- data |>
    filter(
      Antibody == "BS1",
      Occupied == "Receptor",
      Species %in% c("Binary", "Ternary", "Total"),
      Ab.Total >= -2
    ) |>
    plot_heatmap(
      x = Ab.Total,
      y = Recep.Total,
      xlab = "[BS1] (nM)",
      ylab = bquote("[R]"[total] ~ "(#/Cell)")
    )

  # Supplemental heat map with varying total antibody and receptor
  # concentrations for the combination of tocilizumab and 10H2
  vary_heatmap_tocih2 <- data |>
    filter(
      Antibody == "Toci_H2",
      Occupied == "Receptor",
      Species %in% c("Binary", "Ternary", "Total"),
      Ab.Total >= -2
    ) |>
    plot_heatmap(
      x = Ab.Total,
      y = Recep.Total,
      xlab = bquote("[Ab]"[total] ~ "(nM)"),
      ylab = bquote("[R]"[total] ~ "(#/Cell)"),
      palette = "pals::ocean.algae"
    )

  # Line plots with varying antibody concentrations and lines for receptor
  # concentration
  bound_free_ab <- data |>
    filter(
      Antibody == "BS1",
      Ab.Total >= -2,
      Recep.Total %~% c(2, 3, 4, 5, 6, 7),
      Species %in% c("Free", "Binary", "Ternary", "Total")
    ) |>
    mutate(Recep.Total = factor(Recep.Total)) |>
    plot_recep_ab(
      x = Ab.Total,
      color = Recep.Total,
      xlab = "[BS1] (nM)",
      colorlab = bquote("[R]"[total] ~ "(#/Cell)"),
      palette = "Sunset"
    )

  # Line plots with varying receptor concentrations and lines for antibody
  # concentration
  bound_free_recep <- data |>
    filter(
      Antibody == "BS1",
      Ab.Total %~% c(-2, -1, 0, 1, 2, 3),
      Species %in% c("Free", "Binary", "Ternary", "Total")
    ) |>
    mutate(Ab.Total = factor(Ab.Total)) |>
    plot_recep_ab(
      x = Recep.Total,
      color = Ab.Total,
      xlab = bquote("[R]"[total] ~ "(#/Cell)"),
      colorlab = "[BS1] (nM)",
      palette = "Sunset"
    )

  # COMBINE AND SAVE --------------------------------------------------------

  # Save the poster version if requested; otherwise, save the manuscript version
  if (poster) {

    # Save a condensed version of the figure for poster presentations
    layout <- c(
      area(1, 1, 3, 20),
      area(4, 2, 6, 19)
    )

    vary_conc_poster <-
      (bound_free_ab) /
      (vary_heatmap) +
      plot_layout(design = layout) &
      tag_plots()

    savefig(
      vary_conc_poster,
      figpath(fig_file),
      ratio = 1.85,
      width = 15,
      scaling = 1.2
    )
  } else {

    # Main text figure
    layout <- c(
      area(1, 1, 3, 20),
      area(4, 1, 6, 20),
      area(7, 2, 10, 19)
    )

    vary_conc <-
      (bound_free_ab) /
      (bound_free_recep) /
      (vary_heatmap) +
      plot_layout(design = layout) &
      tag_plots(tags = list(c("A", "C", "B", "D", "E")))

    savefig(
      vary_conc,
      figpath(fig_file),
      ratio = 1.5,
      width = 12.25
    )

    # Supplemental figure
    savefig(
      vary_heatmap_tocih2,
      figpath(fig_file, "mAbs"),
      ratio = 3,
      width = 9
    )
  }

  # Return a list with the figures if the output is assigned
  figures <- list(vary_conc, vary_heatmap_tocih2)
  invisible(figures)
}
