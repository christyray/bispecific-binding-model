
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

# Generate the main text and supplemental figures for the simulations restricted
# to monovalent binding
monovalent_only <- function(data_file = NULL, fig_file = NULL,
                          fig_path = NULL, suppress_warn = TRUE, dpi = 600) {

  # Set display name for current figure
  name <- "Monovalent Only"

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

  # Varying antibody and receptor concentrations for monovalent binding only
  data <- import_data(data_file, type = "sim")

  # Combine the occupancy results for all of the antibodies into one data frame,
  # convert the antibody and receptor concentrations to log concentrations, and
  # filter the data to only the necessary values for plotting
  # `modify_depth()` selects the "occupied" tables from inside each of the 1st
  # level lists
  data <- data$occupied |>
    log_conc() |>
    filter(
      Antibody %in% c("BS1", "Toci_H2"),
      Occupied %in% c("Antibody", "Receptor"),
      Species %in% c("Free", "Binary", "Ternary", "Total")
    )

  # PLOTTING ----------------------------------------------------------------

  # Wrapper function to generate the monovalent figure panels depending on which
  # antibodies are being plotted
  plot_monovalent <- function(data, antibody) {

    # Filter the data and set panel labels and colors depending on antibody
    if (antibody == "BS1") {
      data <- data |> filter(Antibody == "BS1")
      ab_lab <- "[BS1] (nM)"
      palettes <- list(
        heatmap = "pals::ocean.deep",
        lines = "BluYl",
        compare = c(
          palette_gen(3, palette = "Sunset", range = c(0.15, 0.9)),
          palette_gen(1, "BluYl", range = c(0.6, 0.6))
        )
      )

    } else if (antibody == "Toci_H2") {
      data <- data |> filter(Antibody == "Toci_H2")
      ab_lab <- bquote("[Ab]"[total] ~ "(nM)")
      palettes <- list(
        heatmap = "pals::ocean.algae",
        lines = "BurgYl",
        compare = c(
          palette_gen(3, palette = "BluYl", range = c(0.15, 0.75)),
          palette_gen(1, "BurgYl", range = c(0.6, 0.6))
        )
      )
    }

    # Varying Concentrations --------------------------------------------------

    # Heat map with varying total antibody and receptor concentrations
    vary_heatmap <- data |>
      filter(
        ParamID == 2,
        Occupied == "Receptor",
        Species == "Total",
        Ab.Total >= -2
      ) |>
      plot_heatmap(
        x = Ab.Total,
        y = Recep.Total,
        xlab = ab_lab,
        ylab = bquote("[R]"[total] ~ "(#/Cell)"),
        palette = palettes$heatmap,
        nbreaks = 4
      ) +
      facet_grid2()

    # Line plots with varying antibody concentrations and lines for receptor
    # concentration
    bound_free_ab <- data |>
      filter(
        ParamID == 2,
        Ab.Total >= -2,
        Recep.Total %~% c(2, 3, 4, 5, 6, 7),
        Species %in% c("Free", "Total")
      ) |>
      mutate(Recep.Total = factor(Recep.Total)) |>
      plot_recep_ab(
        x = Ab.Total,
        color = Recep.Total,
        xlab = ab_lab,
        colorlab = bquote("[R]"[total] ~ "(#/Cell)"),
        palette = palettes$lines
      )

    # Line plots with varying receptor concentrations and lines for antibody
    # concentration
    bound_free_recep <- data |>
      filter(
        ParamID == 2,
        Ab.Total %~% c(-2, -1, 0, 1, 2, 3),
        Recep.Total >= 3,
        Species %in% c("Free", "Total")
      ) |>
      mutate(Ab.Total = factor(Ab.Total)) |>
      plot_recep_ab(
        x = Recep.Total,
        color = Ab.Total,
        xlab = bquote("[R]"[total] ~ "(#/Cell)"),
        colorlab = ab_lab,
        palette = palettes$lines
      )

    # Comparison of Monovalent and Bivalent Binding ---------------------------

    # Specify the factor levels and labels for the receptor concentrations
    recep_levels <- c(3, 4, 5, 6, 7)
    recep_labels <- str_glue('"[R]"[total] == 10^{val}', val = recep_levels)

    # Compare the monovalent and bivalent complex formation over varying
    # antibody concentrations at different receptor levels
    mono_bi <- data |>
      filter(
        Ab.Total >= -2,
        Recep.Total %~% recep_levels,
        Occupied == "Receptor",
        ((ParamID == 1 & Species %in% c("Binary", "Ternary", "Total")) |
           (ParamID == 2 & Species == "Total"))
      ) |>
      mutate(
        Recep.Total = factor_number(
          Recep.Total,
          levels = recep_levels,
          labels = recep_labels
        ),
        Type = case_when(
          ParamID == 1 & Species == "Binary" ~ "Bivalent, Binary",
          ParamID == 1 & Species == "Ternary" ~ "Bivalent, Ternary",
          ParamID == 1 & Species == "Total" ~ "Bivalent, Total",
          ParamID == 2 ~ "Monovalent"
        )
      ) |>
      proper_names() |>
      ggplot(aes(x = Ab.Total, y = Frac, color = Type)) +
      geom_line() +
      scale_x_continuous(
        name = ab_lab,
        labels = label_math(expr = 10^.x)
      ) +
      scale_color_manual(name = "", values = palettes$compare) +
      labs(y = "Fractional Occupancy") +
      facet_grid2(
        . ~ Recep.Total,
        strip = strip_vanilla(clip = "off"),
        labeller = label_parsed
      ) +
      theme_cr()

    # Return a list with all of the created plot panels
    panels <- list(
      heatmap = vary_heatmap,
      ab = bound_free_ab,
      recep = bound_free_recep,
      compare = mono_bi
    )
    return(panels)
  }

  # Create the figure panels for BS1 and for the combination of mAbs
  bs1_panels <- plot_monovalent(data, antibody = "BS1")
  tocih2_panels <- plot_monovalent(data, antibody = "Toci_H2")

  # COMBINE AND SAVE --------------------------------------------------------

  # Main text figure
  mono_only <-
    ((bs1_panels$recep) + bs1_panels$heatmap +
       plot_layout(widths = c(3, 6, 4), guides = "keep")) /
    bs1_panels$compare &
    tag_plots()

  savefig(
    mono_only,
    figpath(fig_file),
    ratio = 2,
    width = 12.5
  )

  # Supplemental figures
  savefig(
    bs1_panels$ab,
    figpath(fig_file, "Ab"),
    ratio = 3,
    width = 9
  )

  mono_only_tocih2 <-
    ((tocih2_panels$recep) + tocih2_panels$heatmap +
       plot_layout(widths = c(3, 6, 4), guides = "keep")) /
    tocih2_panels$compare &
    tag_plots()

  savefig(
    mono_only_tocih2,
    figpath(fig_file, "Toci-H2"),
    ratio = 2,
    width = 12.5
  )

  # Return a list with the figures if the output is assigned
  figures <- list(mono_only, bs1_panels$ab, mono_only_tocih2)
  invisible(figures)
}
