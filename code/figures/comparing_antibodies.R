
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

# Relative Binding Plot Function ------------------------------------------

# Helper function to generate the heatmaps of the relative binding between BS1
# and the mAbs for specific antibody concentrations and plot limits and breaks
relative <- function(data, conc, limits, breaks, labels,
                     palette = "pals::ocean.curl") {

  relative_plot <-
    # Filter the data to the specific antibody concentrations, label the
    # antibody concentrations for the different panels, and label the complex
    # types
    data |>
    filter(Ab.Total %~% conc) |>
    mutate(
      Ab.Total = factor_number(
        Ab.Total,
        conc,
        str_glue('"[Ab]"[total] == {10^val} ~ "nM"', val = conc)
      ),
      Occupied = fct_recode(
        Occupied,
        "Bound IL-6R" = "IL6R",
        "Bound IL-8R" = "IL8R",
        "Total Bound Receptor" = "Receptor"
      )
    ) |>
    # Plot the relative binding heatmap for varying IL6R and IL8R concentrations
    ggplot(aes(x = Recep.IL6R, y = Recep.IL8R, fill = Ratio)) +
    geom_raster() +
    scale_x_continuous(
      name = "[IL-6R] (#/Cell)",
      labels = label_math(expr = 10^.x),
      expand = expand0()
    ) +
    scale_y_continuous(
      name = "[IL-8R] (#/Cell)",
      labels = label_math(expr = 10^.x),
      expand = expand0()
    ) +
    scale_fill_paletteer_c(
      name = "Relative Binding",
      palette = palette,
      limits = limits,
      breaks = breaks,
      labels = labels,
      guide = heatmap_legend(plot_scale = 0.7)
    ) +
    facet_grid2(
      Ab.Total ~ Occupied,
      strip = strip_vanilla(clip = "off"),
      labeller = labeller(Ab.Total = label_parsed)
    ) +
    theme_cr() +
    theme(
      panel.spacing = unit(20, "pt"),
      legend.title = element_text(margin = margin(b = 6))
    )

  # If only one antibody concentration is given, facet by the complex type only
  if (length(conc) == 1) {
    relative_plot <- relative_plot +
      facet_grid2(. ~ Occupied, strip = strip_vanilla(clip = "off"))
  }

  return(relative_plot)
}

# Main Plot Function ------------------------------------------------------

# Generate the main text and supplemental figures for the simulations comparing
# the monospecific and bispecific antibodies
comparing_antibodies <- function(data_file = NULL, fig_file = NULL,
                                 fig_path = NULL, poster = FALSE,
                                 suppress_warn = TRUE, dpi = 600) {

  # Set display name for current figure
  name <- "Comparing Antibodies"

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
  # `modify_depth()` selects the "occupied" tables from inside each of the 1st
  # level lists
  data <- data$occupied |>
    log_conc() |>
    filter(
      Antibody %in% c("BS1", "Toci_H2"),
      Occupied %in% c("Receptor", "IL6R", "IL8R"),
      Species %in% c("Free", "Binary", "Ternary", "Total")
    )

  # PLOTTING ----------------------------------------------------------------

  # Specify narrower range for receptor concentrations for figures in text
  recep_range <- c(3, 6)

  # Specify the factor levels and labels for the supplement figures with
  # multiple antibody concentrations
  supp_levels <- c(0, 1, 2)
  supp_labels <- str_glue('"[Ab]"[total] == {10^val} ~ "nM"', val = supp_levels)

  # Varying IL6R and IL8R Concentrations ------------------------------------

  # Main text heatmap with varying IL6R and IL8R concentrations
  vary_recep_compare <- data |>
    filter(
      Ab.Total %~% 1,
      (Recep.IL6R >= recep_range[1] & Recep.IL6R <= recep_range[2]),
      (Recep.IL8R >= recep_range[1] & Recep.IL8R <= recep_range[2]),
      Occupied == "Receptor",
      Species %in% c("Binary", "Ternary", "Total")
    ) |>
    plot_heatmap(
      x = Recep.IL6R,
      y = Recep.IL8R,
      xlab = "[IL-6R] (#/Cell)",
      ylab = "[IL-8R] (#/Cell)",
      nbreaks = 4
    ) +
    facet_grid2(
      Antibody ~ Species,
      strip = strip_vanilla(clip = "off")
    ) +
    guides(fill = heatmap_legend(plot_scale = 1))

  # Supplemental heatmaps with varying IL6R and IL8R concentrations at
  # multiple antibody concentration levels
  vary_recep_compare_multi <- data |>
    filter(
      Ab.Total %~% supp_levels,
      Occupied == "Receptor",
      Species %in% c("Binary", "Ternary", "Total")
    ) |>
    mutate(Ab.Total = factor_number(
      Ab.Total,
      levels = supp_levels,
      labels = supp_labels
    )) |>
    plot_heatmap(
      x = Recep.IL6R,
      y = Recep.IL8R,
      xlab = "[IL-6R] (#/Cell)",
      ylab = "[IL-8R] (#/Cell)"
    ) +
    facet_nested(
      Ab.Total + Antibody ~ Species,
      strip = strip_nested(clip = "off"),
      labeller = labeller(
        Ab.Total = label_parsed,
        Antibody = c(
          "BS1" = "BS1",
          "Tocilizumab + 10H2" = "mAbs"
        )
      )
    ) +
    guides(fill = heatmap_legend(plot_scale = 1)) +
    coord_fixed()

  # Bound Receptor with Varying Excess Receptor Concentration ---------------

  # Set the discrete receptor concentrations to use for the plot lines
  il6r_fixed <- 3
  il8r_vary <- c(2, 3, 4, 5, 6, 7)

  # Plot the amount of each receptor bound for each antibody with varying IL8R
  # concentration and constant IL6R concentration
  vary_il8r <- data |>
    filter(
      Recep.IL6R %~% il6r_fixed,
      Recep.IL8R %~% il8r_vary,
      Occupied %in% c("IL6R", "IL8R"),
      Species == "Total"
    ) |>
    mutate(
      Recep.Ratio = factor_number(
        Recep.IL8R - Recep.IL6R,
        levels = il8r_vary - il6r_fixed
      ),
      Occupied = fct_recode(
        Occupied,
        "Bound IL-6R" = "IL6R",
        "Bound IL-8R" = "IL8R"
      )
    ) |>
    proper_names() |>
    ggplot(aes(x = Ab.Total, y = Frac, color = Recep.Ratio)) +
    geom_line() +
    scale_x_continuous(
      name = bquote("[Ab]"[total] ~ "(nM)"),
      labels = label_math(expr = 10^.x)
    ) +
    scale_color_manual(
      name = "[IL-8R]:[IL-6R]",
      values = palette_gen(length(il8r_vary), palette = "Emrld"),
      labels = \(x) 10^as.numeric(x)
    ) +
    facet_grid2(Antibody ~ Occupied, strip = strip_vanilla(clip = "off")) +
    labs(y = "Fraction of Receptor Bound") +
    theme_cr() +
    theme(panel.spacing = unit(20, "pt"))

  # Set the discrete receptor concentrations to use for the plot lines
  il6r_vary <- c(2, 3, 4, 5, 6, 7)
  il8r_fixed <- 3

  # Plot the amount of each receptor bound for each antibody with varying IL6R
  # concentration and constant IL8R concentration for the supplement
  vary_il6r <- data |>
    filter(
      Recep.IL6R %~% il6r_vary,
      Recep.IL8R %~% il8r_fixed,
      Occupied %in% c("IL6R", "IL8R"),
      Species == "Total"
    ) |>
    mutate(
      Recep.Ratio = factor_number(
        Recep.IL6R - Recep.IL8R,
        levels = il6r_vary - il8r_fixed
      ),
      Occupied = fct_recode(
        Occupied,
        "Bound IL-6R" = "IL6R",
        "Bound IL-8R" = "IL8R"
      )
    ) |>
    proper_names() |>
    ggplot(aes(x = Ab.Total, y = Frac, color = Recep.Ratio)) +
    geom_line() +
    scale_x_continuous(
      name = bquote("[Ab]"[total] ~ "(nM)"),
      labels = label_math(expr = 10^.x)
    ) +
    scale_color_manual(
      name = "[IL-6R]:[IL-8R]",
      values = palette_gen(length(il6r_vary), palette = "Teal"),
      labels = \(x) 10^as.numeric(x)
    ) +
    facet_grid2(Antibody ~ Occupied, strip = strip_vanilla(clip = "off")) +
    labs(y = "Fraction of Receptor Bound") +
    theme_cr() +
    theme(panel.spacing = unit(20, "pt"))

  # Relative Binding between Antibodies -------------------------------------

  # Filter the data to select the bound IL6R, IL8R, and total receptor, pivot so
  # each antibody is in a separate column, and calculate the ratio of occupied
  # receptor with BS1 to occupied receptor with the mAbs
  data_relative <- data |>
    filter(Species == "Total") |>
    select(Antibody, Ab.Total, Recep.IL6R, Recep.IL8R, Occupied, Frac) |>
    distinct(Antibody, Ab.Total, Recep.IL6R, Recep.IL8R, Occupied,
             .keep_all = TRUE) |>
    pivot_wider(
      names_from = Antibody,
      values_from = Frac
    ) |>
    mutate(Ratio = log10(BS1 / Toci_H2))

  # Determine the limits to use based on the maximum relative binding for each
  # antibody concentration
  limit <- data_relative |>
    filter(
      (Recep.IL6R >= recep_range[1] & Recep.IL6R <= recep_range[2]),
      (Recep.IL8R >= recep_range[1] & Recep.IL8R <= recep_range[2])
    ) |>
    summarize(Max = max(abs(Ratio)), .by = Ab.Total)
  limit1 <- limit |> filter(Ab.Total %~% 1) |> pull(Max)
  # limit3 <- limit |> filter(Ab.Total %~% supp_levels) |> pull(Max) |> max()

  # Generate heatmaps for the main text comparing the relative binding between
  # BS1 and the mAbs for varying receptor expression
  relative_ab <- data_relative |>
    filter(
      (Recep.IL6R >= recep_range[1] & Recep.IL6R <= recep_range[2]),
      (Recep.IL8R >= recep_range[1] & Recep.IL8R <= recep_range[2])
    ) |>
    relative(
      conc = 1,
      breaks = log10(c(1/2, 1, 2)),
      limits = c(-log10(2.001), log10(2.001)),
      labels = c("mAbs 2-fold greater", "Equal", "BsAb 2-fold greater")
    )

  # Generate heatmaps for the supplemental section comparing the relative
  # binding for multiple antibody concentrations
  relative_ab_multi <- relative(
    data_relative,
    conc = supp_levels,
    breaks = log10(c(1/10, 1/3, 1, 3, 10)),
    limits = c(-1.01, 1),
    labels = c(
      "mAbs 10-fold greater",
      "mAbs 3-fold greater",
      "Equal",
      "BsAb 3-fold greater",
      "BsAb 10-fold greater"
    )
  )

  # COMBINE AND SAVE --------------------------------------------------------

  # Main text figure
  layout <- c(
    area(1, 1, 5, 7),
    area(1, 8, 5, 10),
    area(6, 2, 8, 9)
  )

  compare_ab <-
    vary_recep_compare +
    vary_il8r +
    relative_ab +
    plot_layout(design = layout) &
    tag_plots()

  # Save the poster version if requested; otherwise, save the manuscript version
  if (poster) {

    # Save a version of the figure with larger text for poster presentations
    compare_ab_poster <- compare_ab &
      apply_scaling(font_scale = 1.05, set_margin = FALSE) &
      theme(panel.spacing = unit(20, "pt"))

    savefig(
      compare_ab_poster,
      figpath(fig_file),
      ratio = 1.75,
      width = 15
    )

  } else {

    # Save main text figure
    savefig(
      compare_ab,
      figpath(fig_file),
      ratio = 1.75,
      width = 15
    )

    # Save supplemental figures
    savefig(
      vary_recep_compare_multi,
      figpath(fig_file, "Concentrations"),
      ratio = 0.7,
      width = 8.5
    )

    savefig(
      vary_il6r,
      figpath(fig_file, "IL6R"),
      ratio = 1.5,
      width = 7
    )

    savefig(
      relative_ab_multi,
      figpath(fig_file, "Relative"),
      ratio = 1.25,
      width = 9
    )
  }

  # Return a list with the figures if the output is assigned
  figures <- list(compare_ab, vary_recep_compare_multi, relative_ab_multi)
  invisible(figures)
}
