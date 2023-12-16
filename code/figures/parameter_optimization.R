
# LOAD PACKAGES -----------------------------------------------------------

library(tidyverse)
library(cli)
library(scales)
library(sciscales)
library(ggh4x)
library(ggrepel)
library(crthemes)
library(paletteer)
library(patchwork)

# FUNCTIONS ---------------------------------------------------------------

# Helper Functions --------------------------------------------------------

# Helper function to match values between data frames
# Inputs are data frame to search, tibble containing the columns and values to
# filter on (for `rowwise()`, can be obtained with `cur_group()`), column to
# search in filtered data frame, value to search for in filtered data frame, and
# column to return value from
# Similar concept to VLOOKUP() from Excel
match_val <- function(df, match, value_col, value, return_col) {

  # Filter data frame to only include the matching values from the input table
  df <- inner_join(df, match, by = names(match))

  # Select the column containing data values and determine the index where the
  # data value is closest to the input value
  values <- df[[value_col]]
  idx <- which.min(abs(values - value))

  # Select the value in the return column from the row where the data value was
  # closest to the search value
  out <- df[[return_col]][[idx]]
  return(out)
}

# Helper function to calculate the density value for a particular parameter
# value depending on how the data is grouped/faceted
# Inputs are the optimal parameter values and the groups that should be used for
# determining the density, needs to match the facets used on the plot
density_y <- function(df, groups, ...) {

  # Group the input data frame using the provided groups
  group_df <- df |>
    group_by(across(all_of(groups))) |>
    select(all_of(groups), "Optimal")

  # Calculate the density distribution using the same method as is used by
  # ggplot - necessary to determine the y value for a particular parameter
  # Based on `https://stackoverflow.com/q/69699009`

  # Necessary to pre-specify the anonymous function for map to be able to pass
  # the additional arguments in
  # Based on `https://stackoverflow.com/a/64867742`
  dens_opt <- function(x, ...) density(x$Optimal, ...)

  # Nest the grouped data frame to collect all of the parameter values for one
  # group in separate data frames, map the density function across each
  # separate data frame (calculates a separate density distribution for each
  # group, just like the faceted plots), pull the X and Y values from the
  # calculated density distribution into tibbles, remake a table with the new
  # density values, and scale the density values to a range of 0 to 1 (same as
  # the `after_stat(scaled)` from the ggplot)
  dens <- group_df |>
    nest(data = c("Optimal")) |>
    mutate(Density = map(.data$data, dens_opt, ...)) |>
    mutate(Density = map(.data$Density, ~ as_tibble(.x[c("x", "y")]))) |>
    select(-"data") |>
    unnest("Density") |>
    mutate(y = .data$y / max(.data$y))

  # Find the parameter values in each group that give the minimal cost
  minval <- df |>
    group_by(across(all_of(groups))) |>
    filter(Cost == min(.data$Cost)) |>
    select(all_of(groups), "Optimal")

  # Compare the parameter values with minimal cost from each group to the
  # density distributions to determine the Y value for each optimal parameter
  matched <- minval |>
    rowwise() |>
    mutate(
      y = match_val(dens, cur_group(), "x", Optimal, "y"),
      Label = as.character(to_scientific(10^Optimal))
    ) |>
    ungroup()
  return(matched)
}

# Helper function to join the data for just the "BS1, Data" normalization onto
# the data for all of the normalizations; used for the plots where "BS1, Data"
# is compared to the results for all normalizations
# If no argument is given for `all_df`, the function will use the same data
# frame for both
join_bs1 <- function(bs1_df, all_df = NULL) {

  # Use the BS1 data frame if no data frame was given for all normalizations
  all_df <- all_df %||% bs1_df

  bs1_df |>
    filter(.data$Basis == "BS1", .data$NormCalc == "Data") |>
    bind_rows(all_df, .id = "Facet") |>
    mutate(Facet = factor(
      .data$Facet,
      levels = c("2", "1"),
      labels = c("all", "bs1")
    ))
}

# Plot Function -----------------------------------------------------------

# Generate the main text and supplemental figures for the parameter optimization
parameter_optimization <- function(data_file = NULL, fig_file = NULL,
                                   fig_path = NULL, suppress_warn = TRUE,
                                   dpi = 600) {

  # Set display name for current figure
  name <- "Parameter Optimization"

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

  # Optimized parameters from all of the different optimization options
  data <- import_files(data_file, type = "opt")
  data <- clean_data(data, type = "opt")

  # DATA ORGANIZATION -------------------------------------------------------

  # Remove optimizations that are not useful results
  # Flag == 0 removes optimizations that stopped because of the number of
  # function iterations (meaning that the optimization never converged and that
  # the tolerance and step size were still high)
  # The near functions check if the any of the optimized parameters are the same
  # as the initial guesses (within a small relative tolerance) to remove
  # optimizations that never moved

  # Filter on different variables depending on type of optimization done
  # For the newer data, using `match_near()` from `utilities.R` to select values
  # that are the same with relative tolerance, also rejecting an optimization if
  # *any* of the values did not change compared to the original
  cut_rows <- data |>
    dplyr::filter(
      Flag == 0 |
        Cost > 1000 |
        match_near(Param.kon6R, Initial.kon6R, 1e-6) |
        match_near(Param.kon8R, Initial.kon8R, 1e-6) |
        match_near(Param.kon6Rprime, Initial.kon6Rprime, 1e-6) |
        match_near(Param.koff6R, Initial.koff6R, 1e-6) |
        match_near(Param.koff8R, Initial.koff8R, 1e-6)
    )
  data <- dplyr::setdiff(data, cut_rows)

  # Log-transform the initial guesses and optimal values and pivot the data into
  # long form with separate columns for the initial and optimal values
  data <- data |>
    mutate(across(starts_with(c("Initial", "Param")), log10)) |>
    pivot_longer(
      starts_with(c("Initial", "Param")),
      names_to = "Param",
      values_to = "Value"
    ) |>
    separate(Param, c("Type", "Rate")) |>
    pivot_wider(
      names_from = "Type",
      values_from = "Value"
    ) |>
    rename(Optimal = Param, Param = Rate, NormCalc = Calc) |>
    short_names()

  # PLOTTING ----------------------------------------------------------------

  # Plot Options ------------------------------------------------------------

  # Select colors to use for figures
  cost_colors <- paletteer_d("basetheme::clean")[c(1,2,4)]
  norm_colors <- paletteer_d("RColorBrewer::Set1", 4)
  compare_colors <- c(paletteer_d("basetheme::clean")[2], norm_colors[1])

  # Set the range to expand the axes of the grouped cost plots by
  expand_range <- c(0.2, 0.1)

  # Set the facet labels for panels faceted by parameter and normalization type
  facet_labels <- labeller(
    Facet = c("all" = "All Norm", "bs1" = "BS1 Norm"),
    Param = label_parsed
  )

  # Set the ranges and breaks for the X-axes for each parameter

  # Initial = the range of initial guesses used in the optimization
  facet_initial <- list(
    scale_x_continuous(breaks = seq(-8, -2, by = 3), limits = c(-9, -2)),
    scale_x_continuous(breaks = seq(-8, -2, by = 3), limits = c(-9, -2)),
    scale_x_continuous(breaks = seq(-12, -5, by = 3), limits = c(-13, -5)),
    scale_x_continuous(breaks = seq(-6, -2, by = 2), limits = c(-7, -2)),
    scale_x_continuous(breaks = seq(-6, -2, by = 2), limits = c(-7, -2))
  )

  # Optimal = the lower and upper bounds used for the optimization
  facet_optimal <- list(
    scale_x_continuous(breaks = seq(-10, 0, by = 5), limits = c(-11, 0)),
    scale_x_continuous(breaks = seq(-10, 0, by = 5), limits = c(-11, 0)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-8, 0, by =  4), limits = c(-9, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-9, 0))
  )

  # Narrow = +/- one order magnitude from the optimal value
  facet_narrow <- list(
    scale_x_continuous(breaks = seq(-6, -5, by = 1), limits = c(-6.23, -4.23)),
    scale_x_continuous(breaks = seq(-6, -5, by = 1), limits = c(-6.04, -4.04)),
    scale_x_continuous(breaks = seq(-8, -6, by = 1), limits = c(-8.09, -6.09)),
    scale_x_continuous(breaks = seq(-8, -6, by = 1), limits = c(-7.91, -5.91)),
    scale_x_continuous(breaks = seq(-5, -4, by = 1), limits = c(-5.25, -3.25)),
    scale_x_continuous(breaks = seq(-5, -4, by = 1), limits = c(-5.20, -3.20))
  )

  # Initial Guesses with Optimal Values -------------------------------------

  # Rename the data and add column to group the points into three groups by cost
  initial <- data |>
    drop_na() |>
    proper_names() |>
    mutate(
      Cost3 = ifelse(Cost < 2.5, 1, ifelse(Cost < 5, 2, 3)),
      Cost3 = factor(
        Cost3,
        levels = c(1, 2, 3),
        labels = c("< 2.5", "2.5 - 5", "> 5")
      ),
      Facet = factor("all", levels = c("all", "bs1"))
    )

  # Plot the relationship of the optimal values with the initial guesses used to
  # determine how dependent the optimization is on the specific initial guesses

  # Compare all normalization schemes with the "BS1, Data" normalization
  initial_bs1 <- join_bs1(initial) |>
    ggplot(aes(x = Initial, y = Optimal, color = Cost3, alpha = Facet)) +
    geom_point(shape = 16, size = 1.2) +
    facet_grid2(
      Facet ~ Param,
      scales = "free_x",
      labeller = facet_labels,
      strip = strip_vanilla(clip = "off")
    ) +
    facetted_pos_scales(x = facet_initial) +
    scale_y_continuous(breaks = breaks_extended(n = 5)) +
    scale_alpha_manual(values = c(0.1, 0.25)) +
    scale_color_manual(name = "Cost", values = cost_colors) +
    labs(
      x = bquote(log[10] ~ "Initial Guess"),
      y = bquote(log[10] ~ "Optimal Value")
    ) +
    theme_cr(font_scale = 0.95) +
    guides(color = guide_legend(override.aes = list(alpha = 1)), alpha = "none")

  # Facet by normalization scheme used
  initial_norms <- initial |>
    ggplot(aes(x = Initial, y = Optimal, color = Cost3)) +
    geom_point(shape = 16, alpha = 0.2, size = 1.2) +
    facet_nested(
      Basis + NormCalc ~ Param,
      scales = "free_x",
      labeller = label_parsed,
      strip = strip_nested(clip = "off")
    ) +
    facetted_pos_scales(x = facet_initial) +
    scale_y_continuous(breaks = breaks_extended(n = 6)) +
    scale_color_manual(name = "Cost", values = cost_colors) +
    labs(
      x = bquote(log[10] ~ "Initial Guess"),
      y = bquote(log[10] ~ "Optimal Value")
    ) +
    theme_cr() +
    guides(color = guide_legend(override.aes = list(alpha = 1)))

  # Density of Optimal Values -----------------------------------------------

  # Set the smoothing factor for the density distributions
  adjust <- 1

  # Calculate scaled density Y values for the minimal cost parameters
  dens <- data |>
    density_y("Param", adjust = adjust) |>
    proper_names()
  dens_norms <- data |>
    density_y(c("Basis", "NormCalc", "Param"), adjust = adjust) |>
    proper_names()

  # Scaled density of the optimal parameter values - gives a visualization of
  # how tight or widely spread the optimized values for each individual
  # parameter are

  # Compare all normalization schemes with the "BS1, Data" normalization
  param_dens_base <- join_bs1(data) |>
    proper_names() |>
    ggplot(aes(x = Optimal, y = after_stat(scaled), color = Facet)) +
    geom_density(linewidth = 1, adjust = adjust) +
    geom_point(
      data = join_bs1(dens_norms, dens),
      mapping = aes(x = Optimal, y = y),
      color = "gray20"
    ) +
    scale_color_manual(values = compare_colors) +
    facet_grid2(
      Facet ~ Param,
      scales = "free_x",
      labeller = facet_labels,
      strip = strip_vanilla(clip = "off")
    ) +
    labs(x = bquote(log[10] ~ "Parameter Value")) +
    scale_y_continuous(
      name = "Scaled Density",
      breaks = c(0, 0.5, 1),
      expand = expansion(mult = c(0, 0.1))
    ) +
    theme_cr() +
    guides(color = "none")

  # Apply the x-axis limits and breaks for the main text figure
  param_dens_bs1 <-
    param_dens_base +
    facetted_pos_scales(x = facet_optimal)

  # Supplemental panel with narrower X-axis range (+/- one order of magnitude
  # around the optimal value)
  param_dens_narrow <-
    param_dens_base +
    facetted_pos_scales(x = facet_narrow)

  # Facet by normalization scheme used

  # In `geom_label_repel()`, `color` is the text and border color, `fill` is the
  # label fill (the extra digits at the end of the hex code make it slightly
  # transparent), `nudge_x` moves each label along the X axis (same units as the
  # data), `ylim` sets borders on the Y axis that the label won't cross, `size`
  # sets the text size, `label.padding` sets the space between the label text
  # and the surrounding box, `label.size` sets the label border width, and
  # `min.segment.length` sets the shortest distance to draw a line between the
  # label and the point; `force_pull` sets the attraction between the data
  # labels and the points
  param_dens_norms <- data |>
    proper_names() |>
    ggplot(aes(
      x = Optimal, y = after_stat(scaled),
      color = interaction(Basis, NormCalc, lex.order = TRUE)
    )) +
    geom_density(linewidth = 1, adjust = adjust) +
    geom_point(data = dens_norms, aes(x = Optimal, y = y), color = "gray20") +
    geom_label_repel(
      data = dens_norms,
      aes(x = Optimal, y = y, label = Label),
      color = "gray20",
      fill = "#FFFFFF8C",
      family = "Roboto Regular",
      size = 3.5,
      parse = TRUE,
      point.size = 1,
      label.padding = 0.15,
      label.size = 0.3,
      min.segment.length = 1,
      nudge_y = -0.3
    ) +
    scale_color_manual(values = norm_colors) +
    facet_nested(
      Basis + NormCalc ~ Param,
      scales = "free_x",
      labeller = label_parsed,
      strip = strip_nested(clip = "off")
    ) +
    labs(x = bquote(log[10] ~ "Parameter Value")) +
    facetted_pos_scales(x = facet_optimal) +
    scale_y_continuous(name = "Scaled Density", breaks = c(0, 0.5, 1)) +
    theme_cr() +
    guides(color = "none")

  # Cost for Parameter Groups -----------------------------------------------

  # Group parameters by same cost and similar parameter values
  # Select just the columns with parameter values, round the Cost to three
  # decimal points to allow similar points to be grouped together, and convert
  # the table to long form
  group_cost <- data |>
    select(ID, Basis, NormCalc, Cost, Param, Optimal) |>
    mutate(Cost = signif(Cost, 3), Optimal = round(Optimal, 1)) |>
    group_by(Basis, NormCalc, Param, Optimal, Cost) |>
    summarize(Group = cur_group_id(), Count = n(), .groups = "drop_last") |>
    ungroup() |>
    mutate(Facet = factor("all", levels = c("all", "bs1")))

  # Plot the optimized parameter values where the size of the points varies
  # depending on the number of parameter sets that share the same value

  # Compare all normalization schemes with the "BS1, Data" normalization
  group_cost_base <- join_bs1(group_cost) |>
    proper_names() |>
    ggplot(aes(
      x = Optimal, y = log10(Cost), size = Count, alpha = Count, color = Facet
    )) +
    geom_point(shape = 16) +
    scale_y_continuous(
      name = bquote(log[10] ~ "Cost"),
      expand = expansion(expand_range)
    ) +
    scale_size_continuous(range = c(1, 6)) +
    scale_alpha_continuous(range = c(0.2, 0.5), guide = "none") +
    scale_color_manual(values = compare_colors) +
    facet_grid2(
      Facet ~ Param,
      scales = "free_x",
      labeller = facet_labels,
      strip = strip_vanilla(clip = "off")
    ) +
    labs(x = bquote(log[10] ~ "Parameter Value")) +
    theme_cr() +
    guides(color = "none")

  # Apply the x-axis limits and breaks for the main text figure
  group_cost_bs1 <-
    group_cost_base +
    facetted_pos_scales(x = facet_optimal)

  # Supplemental panel with narrower X-axis range (+/- one order of magnitude
  # around the optimal value)
  group_cost_narrow <-
    group_cost_base +
    facetted_pos_scales(x = facet_narrow)

  # Facet by normalization scheme used
  group_cost_norms <- group_cost |>
    proper_names() |>
    ggplot(aes(
      x = Optimal, y = log10(Cost), size = Count, alpha = Count,
      color = interaction(Basis, NormCalc, lex.order = TRUE)
    )) +
    geom_point(shape = 16) +
    scale_y_continuous(
      name = bquote(log[10] ~ "Cost"),
      expand = expansion(expand_range)
    ) +
    scale_size_continuous(range = c(1, 6)) +
    scale_alpha_continuous(range = c(0.3, 0.6), guide = "none") +
    scale_color_manual(values = norm_colors) +
    facet_nested(
      Basis + NormCalc ~ Param,
      scales = "free_x",
      labeller = label_parsed,
      strip = strip_nested(clip = "off")
    ) +
    labs(x = bquote(log[10] ~ "Parameter Value")) +
    facetted_pos_scales(x = facet_optimal) +
    theme_cr() +
    guides(color = "none")

  # Cost for Different Normalizations ---------------------------------------

  # Cumulative distribution of the cost of the optimal values, creates curves of
  # fraction of optimal sets versus cost showing what fraction of optimal sets
  # from each normalization are under the specific cost
  cost_curve <- data |>
    group_by(Basis, NormCalc) |>
    mutate(CostBelow = cume_dist(Cost)) |>
    distinct(Basis, NormCalc, Cost, CostBelow) |>
    proper_names() |>
    ggplot(aes(
      x = log10(Cost),
      y = CostBelow,
      color = interaction(Basis, NormCalc, sep = ", ", lex.order = TRUE)
    )) +
    geom_line() +
    scale_y_continuous(
      name = "Fraction of Optimal Sets",
      limits = c(-0.01, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      expand = expansion(c(0, 0.05))
    ) +
    scale_color_manual(name = "Normalization", values = norm_colors) +
    labs(x = bquote(log[10] ~ "Cost")) +
    theme_cr(font_scale = 0.9)

  # COMBINE AND SAVE --------------------------------------------------------

  # Main text figure
  param_optim <-
    initial_bs1 /
    param_dens_bs1 /
    group_cost_bs1 &
    tag_plots()

  savefig(
    param_optim,
    figpath(fig_file),
    ratio = 0.8,
    width = 9
  )

  # Supplemental figures
  savefig(
    initial_norms,
    figpath(fig_file, "Initial"),
    ratio = 1.5,
    width = 8
  )

  param_optim_narrow <-
    param_dens_narrow /
    group_cost_narrow &
    tag_plots()

  savefig(
    param_optim_narrow,
    figpath(fig_file, "Narrow"),
    ratio = 1.2,
    width = 9
  )

  param_optim_norms <-
    param_dens_norms /
    group_cost_norms &
    tag_plots()

  savefig(
    param_optim_norms,
    figpath(fig_file, "Facets"),
    ratio = 1,
    width = 9
  )

  savefig(
    cost_curve,
    figpath(fig_file, "Cost"),
    ratio = 1.5,
    width = 6
  )

  # Return a list with the figures if the output is assigned
  figures <- list(
    param_optim,
    initial_norms,
    param_optim_narrow,
    param_optim_norms,
    cost_curve
  )
  invisible(figures)
}
