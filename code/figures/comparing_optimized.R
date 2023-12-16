
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

# Compare the binding curves for all the optimized parameter sets
comparing_optimized <- function(data_file = NULL, fig_file = NULL,
                                fig_path = NULL, suppress_warn = TRUE,
                                dpi = 600) {

  # Set display name for current figure
  name <- "Supplemental All Parameter Sets"

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

  # Simulations from all optimized parameter sets
  data <- import_data(data_file$norms, type = "sim")

  # Experimental data for comparison
  exps <- import_data(data_file$experimental, type = "exp")

  # Merge Data --------------------------------------------------------------

  # Only need to keep the antibodies, receptors, parameters, and calculations
  # from a single data set because they are the same for all
  # `map_at()` only maps over the selected list elements
  df <- c("antibodies", "receptors", "parameters", "calculations")
  output <- map_at(
    data, df,
    \(x) x |> filter(FileID == "bs1-data") |> select(-FileID)
  )

  # Process the merged output and normalized tables
  # `separate_wider_delim()` splits the normalization scheme names into separate
  # columns, and the scheme types are converted to factors by `short_names()`
  output$norm <- data$norm |> select(-FileID)
  output$output <- data$output |>
    separate_wider_delim(
      cols = FileID,
      delim = "-",
      names = c("Basis", "NormCalc")
    ) |>
    short_names()

  # Select the table with the normalized output
  norms <- output$norm

  # Filter Data -------------------------------------------------------------

  # Determine antibodies and cell lines present in the simulation output
  ab <- levels(fct_drop(norms$Antibody))
  cells <- levels(fct_drop(norms$Cell))

  # Filter experimental data to just the antibodies and cell lines present
  exps <- filter(exps, Antibody %in% ab, Cell %in% cells)

  # Set the number of replicates for the standard error calculation
  nrep <- max(as.numeric(levels(exps$Rep)))

  # Select the correct normalization basis, convert the concentrations to log
  # scale, and calculate mean and standard error for each group of replicates
  exps <- exps |>
    mutate(Conc = log10(.data$Conc)) |>
    group_by(.data$Conc, .data$Antibody, .data$Cell, .data$Basis) |>
    summarize(
      Mean = mean(.data$Value),
      SE = sd(.data$Value)/sqrt(nrep),
      .groups = "drop"
    ) |>
    proper_names()

  # Convert the concentrations to log scale and set the names for the plot
  norms <- norms |>
    filter(Interest == "receptors") |>
    log_conc() |>
    proper_names()

  # PLOTTING ----------------------------------------------------------------

  # Normalized binding versus antibody concentration, facetted by antibody and
  # normalization basis, comparison of experimental data and simulation output
  compared <-
    ggplot(data = norms) +
    geom_line(
      aes(
        x = Ab.Total,
        y = Value,
        group = interaction(ParamID, NormCalc),
        color = Antibody,
        linetype = "Model Output"
      ),
      alpha = 0.05,
      linewidth = 0.5,
    ) +
    geom_point(
      data = exps,
      aes(x = Conc, y = Mean, fill = "Experimental Data"),
      color = "grey20",
      size = 0.75,
    ) +
    geom_errorbar(
      data = exps,
      aes(x = Conc, y = Mean, ymin = Mean - SE, ymax = Mean + SE),
      color = "grey20",
      width = 0.4
    ) +
    scale_x_continuous(
      name = bquote("[Ab]"[total] ~ "(nM)"),
      labels = label_math(expr = 10^.x)
    ) +
    scale_y_continuous(name = "Normalized Ab Binding") +
    scale_color_cr(name = "", palette = "antibodies") +
    coord_cartesian(xlim = c(-2.5, 3.5), ylim = c(0, 2)) +
    facet_nested(
      Basis + Antibody ~ Cell,
      labeller = labeller(.rows = label_value, .cols = label_parsed),
      strip = strip_nested(clip = "off")
    ) +
    theme_cr(font_scale = 0.9) +
    theme(legend.position = "bottom") +
    guides(
      color = "none",
      linetype = guide_legend(
        title = NULL,
        order = 2,
        override.aes = list(alpha = 1)
      ),
      fill = guide_legend(title = NULL, order = 1)
    )

  # COMBINE AND SAVE --------------------------------------------------------

  # Save the figure for the manuscript
  savefig(
    compared,
    figpath(fig_file),
    ratio = 0.75
  )

  # Return a list with the figures if the output is assigned
  figures <- list(compared)
  invisible(figures)
}
