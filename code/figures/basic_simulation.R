
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

# Generate the main text and supplemental figures for the simulations over time
basic_simulation <- function(data_file = NULL, fig_file = NULL,
                          fig_path = NULL, suppress_warn = TRUE, dpi = 600) {

  # Set display name for current figure
  name <- "Basic Simulation"

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

  # Binding over time for all antibodies
  binding_time <- import_data(data_file, type = "sim")

  # Convert the antibody and receptor concentrations to log concentrations
  binding_time$occupied <- log_conc(binding_time$occupied)

  # PLOTTING ----------------------------------------------------------------

  # Plot Options ------------------------------------------------------------

  # Define vector of colors to choose from
  bind_colors <- c("#03618C", "#983434", "#37BDFB", "#D58585",
                   "#906EAF", "#E6B450")

  # Antibody-Receptor Complexes over Time -----------------------------------

  # Wrapper function to plot the bound concentrations versus time for different
  # antibodies and concentrations
  plot_binding <- function(ab, conc, cells = c("6R", "8R", "6R8R"),
                           cell_label = "6R", colors = bind_colors,
                           ab_label = "BS1") {

    # Create the data frame for the washout time label
    washout_label <- data.frame(Time = 1.8, Bound = 1.3e6, Cell = cell_label) |>
      proper_names()

    # Set the receptors to be filtered out based on which antibody is being
    # plotted
    if (ab == "Toci") {
      rm_recep <- "IL8R"
    } else if (ab == "H2") {
      rm_recep <- "IL6R"
    } else {
      rm_recep <- ""
    }

    # Create plot
    binding_time$occupied |>
      filter(
        Antibody == ab,
        Ab.Total %~% conc,
        Cell %in% cells,
        Time <= 4,
        Occupied == "Receptor",
        !Species %in% c("Free", "Binary", "Ternary", "Total", rm_recep)
      ) |>
      rename_species() |>
      proper_names() |>
      ggplot(aes(x = Time, y = Bound, color = Species)) +
      geom_line(alpha = 0.8) +
      geom_vline(xintercept = 2, linetype = "dashed", color = "gray20") +
      geom_text(
        data = washout_label,
        label = "Washout",
        color = "gray20",
        hjust = "right",
        fontface = "italic",
        size = 4
      ) +
      scale_y_continuous(
        name = "Concentration (#/Cell)",
        labels = label_sci(digits = 3),
        breaks = breaks_extended(6),
        expand = expand1()
      ) +
      scale_color_manual(values = colors) +
      xlab("Time (hr)") +
      facet_grid2(
        Antibody ~ Cell,
        strip = strip_vanilla(clip = "off"),
        labeller = labeller(Antibody = ab_label, Cell = label_parsed)
      ) +
      theme_cr()
  }

  # Generate Figures --------------------------------------------------------

  # BS1 at original concentration
  bound_bs1 <-
    plot_binding(ab = "BS1", conc = 2, colors = bind_colors[1:5]) +
    facet_grid2(
      . ~ Cell,
      strip = strip_vanilla(clip = "off"),
      labeller = label_parsed
    )

  # Tocilizumab
  bound_toci <-
    plot_binding(
      ab = "Toci",
      conc = 2,
      cells = c("6R", "6R8R"),
      colors = bind_colors[c(1, 3, 6)],
      ab_label = c("Tocilizumab" = "Tocilizumab, [Ab] = 100 nM")
    )

  # 10H2
  bound_h2 <-
    plot_binding(
      ab = "H2",
      conc = 2,
      cells = c("8R", "6R8R"),
      cell_label = "8R",
      colors = bind_colors[c(2, 4, 6)],
      ab_label = c("10H2" = "10H2, [Ab] = 100 nM")
    )

  # BS1 at lower concentration
  bound_bs1_low <-
    plot_binding(
      ab = "BS1",
      conc = 1,
      colors = bind_colors[1:5],
      ab_label = c("BS1" = "BS1, [Ab] = 10 nM")
    )

  # COMBINE AND SAVE --------------------------------------------------------

  # Main text figure
  savefig(
    bound_bs1,
    figpath(fig_file),
    ratio = 2.25,
    width = 8
  )

  # Supplemental figure
  bound_supp <-
    bound_bs1_low / bound_toci / bound_h2 &
    tag_plots()

  savefig(
    bound_supp,
    figpath(fig_file, "Antibodies"),
    ratio = 0.75,
    width = 9
  )

  # Return a list with the figures if the output is assigned
  figures <- list(bound_bs1, bound_supp)
  invisible(figures)
}
