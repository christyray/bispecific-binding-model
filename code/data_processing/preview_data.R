here::i_am("code/data_processing/preview_data.R")

# SOURCE FUNCTIONS --------------------------------------------------------

# Import helper functions to import and organize the simulation data
source(here::here("code/data_processing/initialization.R"))

# PLOTTING ----------------------------------------------------------------

# Function to quickly generate a concentration versus time figure for varying
# antibody and receptor concentration; takes input of occupancy data output from
# the `occupancy()` function and an assortment of plot options

# `data` = simulation data passed through the `occupancy()` calculation function
# `fraction` = logical value; dictates if fractional occupancy (default) or
# total bound concentration should be plotted
# `ab` = which antibody to plot; used to filter the `Antibody` column
# `cell` = which cell line to plot; used to filter the `Cell` column
# `conc` = range of antibody concentrations to include; defaults to all
# `recep` = range of receptor concentrations to include; defaults to all
# `endtime` = maximum time value for X axis in hours
# `ratio` = ratio of IL6R to IL8R concentrations to plot; defaults to 1
# `individual` = logical value; dictates if the species should be grouped into
# "Free", "Binary", "Ternary", "Bound" (default) or if individual species should
# be plotted
# `labels` = list of receptor concentrations and their labels to be passed to
# `factor_number()`; should have two elements: "levels" for the concentration
# values and "labels" for the labels for those concentrations; generated
# automatically if it is not provided
# `free_scales` = option for the `scales` argument of `facet_grid`; options are
# "fixed" (default), "free_x", "free_y", or "free"
# `independent` = option for the `independent` argument of `facet_grid`; options
# are "none" (default), "x", "y", or "all"
# `output` = option for the function to output a `ggplot` figure ("figure",
# default) or the underlying data used for the figure ("data")

# Example:
# prefix <- "2023-04-18_base"
# data <- import_files(prefix, type = "sim")
# data <- clean_data(data, type = "sim")
# vary_total(data$occupied)

vary_total <- function(data, fraction = TRUE, ab = "BS1", cell = "6R8R",
                       conc = c(-Inf, Inf), recep = c(-Inf, Inf),
                       endtime = 6, ratio = 1, individual = FALSE,
                       labels = NULL, free_scales = "fixed",
                       independent = "none", output = "figure") {


  # Options -----------------------------------------------------------------
  # Setup function variables based on the input argument

  # Plot the fractional occupancy or the total bound concentration, based on the
  # `fraction` argument
  if (fraction) {
    value <- "Frac"
    yaxis <- "Fractional Occupancy"
  } else {
    value <- "Bound"
    yaxis <- "Concentration (#/Cell)"
  }

  # Specify which species to plot based on the `individual` argument, either the
  # summed free molecules and complexes or the individual species
  all_species <- levels(data$Species)
  summed_species <- c("Free", "Binary", "Ternary", "Total")
  if (individual) {
    keep_species <- all_species[!all_species %in% summed_species]
    species_levels <- keep_species
  } else {
    keep_species <- summed_species
    species_levels <- c("Free Ab", "Free R", keep_species)
  }

  # If a list of factor levels and labels for the receptor concentrations was
  # not provided, generate the levels and labels from the values in the
  # Recep.Total column; note that this can be slightly slow for very large data
  # frames
  if (is.null(labels)) {
    concentrations <- unique(data$Recep.Total)
    labels <- list(levels = concentrations, labels = log10(concentrations))
  }

  # Match the "Antibody" column to the possible free antibodies present in the
  # system; necessary because the antibody column does not dictate monovalent vs
  # bivalent, but the species do, also needed for the combination of multiple
  # antibodies
  antibodies <- list(
    "Toci" = c("TociM", "TociB"),
    "H2" = c("H2M", "H2B"),
    "BS1" = "BS1",
    "Toci_H2" = c("TociM", "TociB", "H2M", "H2B")
  )

  # Data Organization -------------------------------------------------------
  # Filter the data based on the function options and calculate columns

  data <- data |>
    # Filter the data to only the desired antibodies, cell line, species (summed
    # or individual), concentrations, and time range
    dplyr::filter(
      (Species %in% c(antibodies[[ab]], "Free") & Occupied == "Antibody" |
         Occupied == "Receptor"),
      Species %in% keep_species,
      Antibody == ab,
      Cell ==  cell,
      Ab.Total >= conc[1],
      Ab.Total <= conc[2],
      Recep.Total >= recep[1],
      Recep.Total <= recep[2],
      Time <= endtime
    ) |>
    # Convert the receptor concentrations to factors to make it easier to use
    # discrete colors for each concentration, calculate the receptor ratio,
    # select the correct value column for the data being plotted (fractional
    # occupancy or total concentration), and label the free antibody and
    # receptor separately
    dplyr::mutate(
      Recep.Total = factor_number(Recep.Total, labels$levels, labels$labels),
      Ratio = Recep.IL6R / Recep.IL8R,
      Value = .data[[value]],
      Species = dplyr::case_when(
        .data$Species == "Free" & .data$Occupied == "Antibody" ~ "Free Ab",
        .data$Species == "Free" & .data$Occupied == "Receptor" ~ "Free R",
        .default = .data$Species
      )
    ) |>
    # Re-factor the species column in the correct order for the figure
    dplyr::mutate(Species = factor(Species, levels = species_levels)) |>
    # Filter the receptor ratio to only the input ratio or any rows where the
    # ratio could not be calculated because only one receptor was present
    dplyr::filter(Ratio %~% ratio | Ratio %in% c(0, Inf))

  # Create Plot -------------------------------------------------------------
  # Create a ggplot of the data if requested

  if (output == "data") {
    # If only the underlying data was requested, output the data and exit the
    # function without plotting
    return(data)
  } else if (output == "figure") {

    # Use the `palette_gen()` function from `utilities.R` to generate a color
    # palette with different colors for each receptor concentration; then,
    # select only the colors needed for the concentrations being plotted
    # This makes the colors consistent for the same receptor concentration
    # regardless of which concentrations were selected for the plot
    colors <- palette_gen(nlevels(data$Recep.Total))
    colors <- colors[unique(data$Recep.Total)]

    # X-axis = time (hr); y-axis = concentration (fractional or total); color =
    # receptor concentration; facet rows = species (summed or individual); facet
    # columns = antibody concentration
    data |>
      proper_names() |>
      ggplot2::ggplot(ggplot2::aes(x = Time, y = Value, color = Recep.Total)) +
      ggplot2::geom_line() +
      ggplot2::scale_x_continuous(name = "Time (hr)", limits = c(0, endtime)) +
      ggplot2::scale_y_continuous(
        name = yaxis,
        labels = sciscales::label_sci()
      ) +
      ggplot2::scale_color_manual(
        name = "[R] (#/Cell)",
        values = colors,
        labels = scales::label_math(expr = 10^.x)
      ) +
      ggh4x::facet_grid2(
        Species ~ Ab.Total,
        labeller = ggplot2::labeller(
          Ab.Total = function(conc) paste(conc, "nM")
        ),
        strip = ggh4x::strip_vanilla(clip = "off"),
        scales = free_scales,
        independent = independent
      ) +
      crthemes::theme_cr()
  }
}
