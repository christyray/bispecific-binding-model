
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

# Generate the main text and supplemental figures for the binding curves
binding_curve <- function(data_file = NULL, fig_file = NULL,
                          fig_path = NULL, poster = FALSE,
                          suppress_warn = TRUE, dpi = 600) {

  # Set display name for current figure
  name <- "Binding Curve"

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

  # Simulation output for binding curve including high concentrations
  sims <- import_data(data_file, type = "sim")
  norm <- sims$norm

  # Experimental data
  exp_file <- "Flow-MFI_"
  exps <- import_files(exp_file, type = "exp")
  exps <- clean_data(exps, type = "exp")

  # DATA ORGANIZATION -------------------------------------------------------

  # Log-scale Ab concentrations
  norm <- log_conc(norm)

  # Determine antibodies and cell lines present in the simulation output and the
  # normalization basis used
  ab <- levels(fct_drop(norm$Antibody))
  cells <- levels(fct_drop(norm$Cell))
  basis <- levels(fct_drop(norm$Basis))

  # Filter the experimental data to the same antibodies, cell lines, and
  # normalization basis used for the model output; convert the concentrations to
  # log scale; and calculate the mean and standard deviation for each set of
  # replicates
  exps <- convert_exp(exps, ab = ab, cell = cells, basis = basis)

  # Determine the maximum concentration used in the experimental data, to mark
  # any model simulations above the experimental concentration range
  max_ab <- max(exps$Conc)

  # PLOTTING ----------------------------------------------------------------

  # Binding Curve for Higher Concentrations with Experimental Data ----------

  # Filter the model output to just the normalized bound antibody, and label any
  # initial antibody concentrations that are above the experimental range
  norm_ab <- norm |>
    filter(Interest == "receptors") |>
    mutate(HighConc = ifelse(Ab.Total > max_ab, TRUE, FALSE)) |>
    proper_names()

  binding <- ggplot() +
    # Solid lines for model simulations within the experimental concentration
    # range, dashed lines for model simulations outside the experimental
    # concentration range, points for the individual data points, and error bars
    # for the standard error of the replicates
    geom_line(
      data = filter(norm_ab, HighConc == FALSE),
      aes(x = Ab.Total, y = Value, color = Antibody, group = Antibody,
          linetype = HighConc),
      alpha = 0.6
    ) +
    geom_line(
      data = filter(norm_ab, HighConc == TRUE),
      aes(x = Ab.Total, y = Value, color = Antibody, group = Antibody,
          linetype = HighConc),
      alpha = 0.6
    ) +
    geom_errorbar(
      data = exps,
      aes(x = Conc, y = Mean, ymin = Mean - SE, ymax = Mean + SE),
      linewidth = 0.55,
      width = 0.5,
      color = "gray40"
    ) +
    geom_point(
      data = exps,
      aes(x = Conc, y = Mean, color = Antibody, shape = Antibody),
      alpha = 1
    ) +
    coord_cartesian(xlim = c(-2, 6)) +
    scale_x_continuous(
      name = bquote("[Ab]"[total] ~ "(nM)"),
      breaks = c(-2, 0, 2, 4, 6),
      labels = label_math(expr = 10^.x)
    ) +
    scale_y_continuous(name = "Normalized Ab Binding") +
    # Setting the same name for the scales will combine them into one legend
    scale_color_cr(name = "Antibody", palette = "antibodies") +
    scale_shape_manual(name = "Antibody", values = c(15, 16, 17, NA)) +
    scale_linetype_discrete(
      name = "Concentrations",
      labels = c("Experimental", "Simulated")
    ) +
    facet_grid2(
      . ~ Cell,
      labeller = label_parsed,
      strip = strip_vanilla(clip = "off")
    ) +
    theme_cr(plot_scale = 1.2) +
    # Since color and shape are combined into a single legend, they both have to
    # be assigned to the same place in the order
    guides(
      color = guide_legend(order = 1),
      shape = guide_legend(order = 1),
      linetype = guide_legend(order = 2)
    )

  # Binding Curve Divided between Complexes ---------------------------------

  binding_complexes <- norm |>
    filter(
      Antibody %in% c("BS1", "Toci_H2"),
      Time %~% 2.25
    ) |>
    mutate(Interest = fct_recode(Interest, complex = "receptors")) |>
    proper_names() |>
    mutate(Interest = fct_recode(Interest, "Total Bound" = "Total Complexes")) |>
    ggplot(aes(x = Ab.Total, y = Value, color = Interest)) +
    geom_line() +
    scale_x_continuous(
      name = bquote("[Ab]"[total] ~ "(nM)"),
      breaks = seq(-2, 6, by = 2),
      labels = label_math(expr = 10^.x)
    ) +
    scale_color_cr(name = "", palette = "contrast3") +
    labs(y = "Normalized Binding") +
    facet_grid2(
      Antibody ~ Cell,
      strip = strip_vanilla(clip = "off"),
      labeller = labeller(Cell = label_parsed)
    ) +
    theme_cr()

  # COMBINE AND SAVE --------------------------------------------------------

  # Save the poster version if requested; otherwise, save the manuscript version
  if (poster) {

    # Save a version of the figure with slightly less space for posters
    binding_poster <- binding + theme(legend.spacing = unit(12, "pt"))

    savefig(
      binding_poster,
      figpath(fig_file),
      ratio = 2.75,
      width = 15,
      scaling = 1.45
    )
  } else {
    savefig(
      binding,
      figpath(fig_file),
      ratio = 2.75,
      width = 11
    )
  }

  # Supplemental figure
  savefig(
    binding_complexes,
    figpath(fig_file, "Complexes"),
    ratio = 1.75,
    width = 9
  )

  # Return a list with the figures if the output is assigned
  figures <- list(binding, binding_complexes)
  invisible(figures)
}
