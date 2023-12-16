
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

# Generate the figure for the local and global sensitivity analyses
univariate_sensitivity <- function(data_file = NULL, fig_file = NULL,
                                   fig_path = NULL, poster = FALSE,
                                   suppress_warn = TRUE, dpi = 600) {

  # Set display name for current figure
  name <- "Univariate Sensitivity"

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

  # PLOTTING ----------------------------------------------------------------

  # Local Sensitivity -------------------------------------------------------

  # Set labels for varied parameters in the sensitivity analysis
  local_labels <- list(
    '"[BS1]"' = "Ab",
    '"[IL-6R]"' = "IL6R",
    '"[IL-8R]"' = "IL8R",
    'k["on,6R"]' = "kon6R",
    'k["on,8R"]' = "kon8R",
    'k["on,6R*"]' = "kon6Rprime",
    'k["on,8R*"]' = "kon8Rprime",
    'k["off,6R"]' = "koff6R",
    'k["off,8R"]' = "koff8R"
  )

  # Filter the occupancy data to just the final concentrations for the ternary
  # and total complexes, and relabel the factors to match the output table
  local_occupied <- data$occupied |>
    filter(
      FileID == "local",
      Species %in% c("Ternary", "Total"),
      Occupied == "Receptor",
      Calc == "conct",
      Time == 2
    ) |>
    mutate(
      Interest = case_when(
        Species == "Ternary" ~ "ternary",
        Species == "Total" ~ "complexes"
      ),
      Add = FALSE
    ) |>
    select(!Occupied)

  # Filter the output table to just the area under the curve and peak
  # concentrations for the ternary and total complexes, and add columns matching
  # the occupancy table before appending the two tables together
  local_output <- data$output |>
    filter(
      FileID == "local",
      Calc %in% c("auc", "peak"),
      (Interest == "present" & Species == "BS1_6R_8R") |
        (Interest == "complexes")
    ) |>
    mutate(Bound = Value, Frac = 0, .keep = "unused") |>
    mutate(Interest = case_when(
      Interest == "present" ~ "ternary",
      Interest == "complexes" ~ "complexes"
    )) |>
    bind_rows(local_occupied)

  # Determine the values for each output with the baseline parameters and
  # concentrations, and select only the columns with the relevant output
  baseline <- local_output |>
    filter(ID == 1) |>
    select(!c("FileID", "ID", "ConcID", "ParamID", "OutputID", "Antibody",
              "Cell", "Ab.BS1", "Recep.Total", "Time"))

  # Calculate the sensitivity of the model to each parameter value by
  # calculating the percentage change in each output relative to the baseline
  # value and then dividing by the percent change in the parameter
  compare <- local_output |>
    select(!c("Ab.BS1", "Recep.Total")) |>
    # Join the baseline values onto the output table
    left_join(
      baseline,
      by = c("Interest", "Calc", "Add", "Species"),
      suffix = c("_sim", "_base")
    ) |>
    # Pivot the table into long form with separate rows for the values of each
    # parameter for each simulation
    pivot_longer(
      cols = ends_with(c("_base", "_sim")),
      names_to = c("Variable", ".value"),
      names_sep = "_"
    ) |>
    # Calculate the sensitivity as the percent change in the output divided by
    # the percent change in the parameter value (set as 10 percent in the
    # simulation)
    # Only the parameters that were actually varied will have values in this
    # calculation
    mutate(delta = (sim - base) / base / 0.1, .keep = "unused") |>
    # Convert the table back to wide form with separate columns for each varied
    # parameter and concentration
    pivot_wider(
      names_from = Variable,
      values_from = delta
    ) |>
    # Add up the relative change in each input and filter to only the rows where
    # just one parameter is varied (i.e., when the relative change is 1 for that
    # parameter and 0 for all the others)
    # Necessary because the MATLAB code is not set up to only vary one parameter
    # at a time, so there are extra simulations where multiple parameters were
    # changed
    mutate(Varied = add_across("(Ab|Recep|Param)\\.")) |>
    filter(Varied %~% 1) |>
    select(!Varied) |>
    # Convert the table back to long form with separate rows for each parameter
    # and concentration for each simulation, and select only the rows with the
    # specific parameter being varied in the simulation
    pivot_longer(
      cols = matches("(Ab|Recep|Param)\\."),
      names_to = "Parameter",
      values_to = "Varied"
    ) |>
    filter(Varied %~% 1)

  # Generate the heat map with the local sensitivity of variables of interest
  local_plot <- compare |>
    # Relabel the output and parameter columns with the desired labels for the
    # figure, and reverse the output label order so it looks correct from top to
    # bottom in the final figure
    mutate(
      Output = factor(
        case_when(
          Interest == "ternary" & Calc == "auc" ~ "Ternary AUC",
          Interest == "ternary" & Calc == "conct" ~ "Ternary\nFinal Conc",
          Interest == "ternary" & Calc == "peak" ~ "Ternary\nPeak Conc",
          Interest == "complexes" & Calc == "auc" ~ "Total AUC",
          Interest == "complexes" & Calc == "conct" ~ "Total\nFinal Conc",
          Interest == "complexes" & Calc == "peak" ~ "Total\nPeak Conc"
        ),
        levels = c("Ternary AUC", "Ternary\nFinal Conc", "Ternary\nPeak Conc",
                   "Total AUC", "Total\nFinal Conc", "Total\nPeak Conc")
      ),
      Output = fct_rev(Output),
      Parameter = factor(
        Parameter,
        levels = c("Ab.Total", "Recep.IL6R", "Recep.IL8R",
                   "Param.kon6R", "Param.kon8R", "Param.kon6Rprime",
                   "Param.kon8Rprime", "Param.koff6R", "Param.koff8R"),
        labels = c("Ab", "IL6R", "IL8R", "kon6R", "kon8R",
                   "kon6Rprime", "kon8Rprime", "koff6R", "koff8R")
      ),
      Parameter = fct_recode(Parameter, !!!local_labels)
    ) |>
    # Create the heat map from the sensitivity data
    ggplot(aes(x = Parameter, y = Output, fill = Bound)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_x_discrete(name = "", labels = label_parse(), expand = expand0()) +
    scale_y_discrete(name = NULL, expand = expand0()) +
    scale_fill_paletteer_c(
      name = "Sensitivity",
      palette = "pals::ocean.balance",
      guide = heatmap_legend(plot_scale = 0.8),
      limits = c(-round(max(compare$Bound), 2), round(max(compare$Bound), 2))
    ) +
    theme_cr()

  # Global Sensitivity ------------------------------------------------------

  # Select the baseline values for each parameter for comparison
  global_params <- data$occupied |>
    filter(FileID == "global", ParamID == 1) |>
    distinct(pick(starts_with("Param.")))

  # Set the labels for the species facets
  global_labels <- list(
    "Ternary Complexes" = "Ternary",
    "Total Bound" = "Total"
  )

  # Calculate the change in each parameter value and filter to the specific
  # simulations for plotting
  global_occupied <- data$occupied |>
    filter(FileID == "global") |>
    # Calculate the fold change in each parameter from the baseline
    mutate(across(starts_with("Param."), ~ log10(.x/global_params$.x))) |>
    # Convert the concentrations to log scale
    mutate(Ab.Total = log10(Ab.Total)) |>
    # Filter to just the occupied receptor for a limited range of antibody
    # concentrations and for just one receptor concentration at a 1:1 ratio,
    # plotting just the ternary and total complexes at the final time point
    filter(
      Occupied == "Receptor",
      Ab.Total %~% c(-2, -1, 0, 1, 2, 3),
      (Recep.IL6R / Recep.IL8R) %~% 1,
      Recep.Total %~% 1e5,
      Species %in% c("Ternary", "Total"),
      Calc == "end"
    ) |>
    # Convert the table to long form with separate rows for each parameter
    pivot_longer(
      cols = starts_with("Param."),
      names_to = "Param",
      values_to = "Value",
      names_prefix = "Param."
    ) |>
    # Filter to only the rows with the baseline parameters and the individual
    # parameters being varied in each simulation
    filter(ID == 1 | Value != 0)

  # Create the line plot from the global sensitivity data
  global_colors <- palette_gen(6, palette = "Sunset")
  global_plot <- global_occupied |>
    mutate(Ab.Total = factor(Ab.Total)) |>
    mutate(Species = fct_recode(Species, !!!global_labels)) |>
    proper_names() |>
    ggplot(aes(x = Value, y = Frac, color = Ab.Total)) +
    geom_line(alpha = 0.8) +
    scale_x_continuous(
      name = bquote(Delta ~ log[10]("Parameter")),
    ) +
    scale_color_manual(
      name = "[BS1] (nM)",
      values = global_colors,
      labels = label_math(expr = 10^.x)
    ) +
    labs(y = "Fractional Occupancy") +
    facet_grid2(
      Species ~ Param,
      strip = strip_vanilla(clip = "off"),
      labeller = labeller(Param = label_parsed)
    ) +
    theme_cr()

  # COMBINE AND SAVE --------------------------------------------------------

  # Set the plot layout
  uni_sens <- local_plot +
    global_plot +
    plot_layout(heights = c(8, 12)) +
    plot_annotation(tag_levels = "A")

  uni_sens[[2]] <- uni_sens[[2]] +
    theme(axis.title.y = element_text(vjust = -12))

  # Save the poster version if requested; otherwise, save the manuscript version
  if (poster) {

    # Save a condensed version of the figure for poster presentations
    savefig(
      global_plot,
      figpath(fig_file),
      ratio = 3.5,
      width = 15,
      scaling = 1.6
    )

  } else {

    # Save the figure for the manuscript
    savefig(
      uni_sens,
      figpath(fig_file),
      ratio = 1.25,
      width = 10
    )
  }

  # Return a list with the figures if the output is assigned
  figures <- list(uni_sens)
  invisible(figures)
}
