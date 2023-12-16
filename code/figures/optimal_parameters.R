here::i_am("code/figures/optimal_parameters.R")

# SOURCE FUNCTIONS --------------------------------------------------------

# Import helper functions to import and organize the simulation data
source(here::here("code/data_processing/initialization.R"))

# LOAD PACKAGES -----------------------------------------------------------

library(tidyverse)
library(NbClust)
library(ggforce)
library(scales)
library(sciscales)
library(ggh4x)
library(ggrepel)
library(crthemes)
library(ggthemes)
library(paletteer)
library(patchwork)
library(gt)
library(gtExtras)

# USER OPTIONS ------------------------------------------------------------

# Option to save figures when script is run
saveQ <- TRUE

# Set file name for naming output figures
fn <- "optimal-parameters"
figpath <- paste0("output/figures/", fn)

# Wrapper function to only save figures when requested and to suppress warnings
crsave_nowarn <- function(...) {
  if (saveQ) {
    suppressWarnings(crsave(...))
  }
}

# DATA IMPORT -------------------------------------------------------------

# Optimized parameters from all of the different optimization options
prefix <- "optimization"
output <- import_files(prefix, type = "opt")
output <- clean_opt(output)

# DATA ORGANIZATION -------------------------------------------------------

# Remove optimizations that are not useful results
# Flag == 0 removes optimizations that stopped because of the number of function
# iterations (meaning that the optimization never converged and that the
# tolerance and step size were still high)
# The near functions check if the optimized parameters are the same as the
# initial guesses (within a small tolerance) to remove optimizations that never
# moved

# Filter on different variables depending on which type of optimization was done
# For the newer data, using `match_near()` from `utilities.R` to select values
# that are the same with relative tolerance, also rejecting an optimization if
# *any* of the values did not change compared to the original
if ("Param.koff6R" %in% names(output)) {
  cut_rows <- output |>
    filter(
      Flag == 0 |
        Cost > 1000 |
        match_near(Param.kon6R, Initial.kon6R, 1e-6) |
        match_near(Param.kon8R, Initial.kon8R, 1e-6) |
        match_near(Param.kon6Rprime, Initial.kon6Rprime, 1e-6) |
        match_near(Param.koff6R, Initial.koff6R, 1e-6) |
        match_near(Param.koff8R, Initial.koff8R, 1e-6)
    )
} else {
  cut_rows <- output |>
    filter(Flag == 0 | (
      near(Param.kon6R, Initial.kon6R) &
        near(Param.kon8R, Initial.kon8R) &
        near(Param.kon6Rprime, Initial.kon6Rprime)
    ))
}

output <- setdiff(output, cut_rows)

# Write the filtered parameter values out to a single CSV file; replacing
# periods in column names to be compatible with MATLAB
if (saveQ) {
  out_fn <- "optimal_parameters.csv"
  output_save <- output |>
    rename_with(~ str_replace(.x, "\\.", "_"))
  write_csv(output_save, here("output", out_fn))
}

# Rename the columns with the parameter values for simplicity and capitalize
# the factor names
output <- output |>
  mutate(
    Basis = fct_recode(Basis, BS1 = "bs1", Ab = "ab"),
    Calc = fct_recode(Calc, Data = "data", Max = "max")) |>
  rename(NormCalc = Calc)

# Convert the second binding steps to nM-1 s-1, log-transform the initial
# guesses and optimal values and pivot the data into long form with separate
# columns for the initial and optimal values
alpha <- 8.3029e-07      # nM/(# receptors/cell)
output <- output |>
  mutate(across(contains("prime"), \(x) x / alpha, .names = "{.col}.nM")) |>
  mutate(across(starts_with(c("Initial", "Param")), log10)) |>
  pivot_longer(
    starts_with(c("Initial", "Param")),
    names_to = "Param",
    values_to = "Value"
  ) |>
  separate_wider_delim(
    Param,
    delim = ".",
    names = c("Type", "Rate", "Units"),
    too_few = "align_start"
  ) |>
  pivot_wider(
    names_from = "Type",
    values_from = "Value"
  ) |>
  rename(Optimal = Param, Param = Rate) |>
  short_names()

# Determine the standard deviation of the log-transformed optimal parameters
output_sd <- output |>
  summarize(LogSD = sd(Optimal), SD = sd(10 ^ Optimal), .by = c(Param, Units))

# Remove the converted rate constants from the data frame before plotting
output <- output |>
  filter(is.na(Units)) |>
  select(-Units)

# PLOTTING ----------------------------------------------------------------

# Define Plot Colors ------------------------------------------------------

# Select colors to use for figures
colors3 <- paletteer_d("basetheme::clean")[c(1,2,4)]
colors4 <- paletteer_d("RColorBrewer::Set1", 4)
colors2 <- c(paletteer_d("basetheme::clean")[2], colors4[1])

# Set the range to expand the cost points by
expand_range <- c(0.2, 0.1)

# Initial Guesses with Optimal Values -------------------------------------

# Rename the data and add a column to group the points into three groups by cost
initial <- output |>
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

# Determine the number and percentage of optimal sets in each cost group
initial_count <- initial |>
  distinct(ID, Opt, Basis, NormCalc, Cost3) |>
  summarize(Count = n(), .by = Cost3) |>
  mutate(Percent = Count / sum(Count))

# Plot the relationship of the optimal values with the initial guesses used to
# determine how dependent the optimization is on the specific initial guesses
initial_optimal <- initial |>
  ggplot(aes(x = Initial, y = Optimal, color = Cost3)) +
  geom_point(shape = 16, alpha = 0.1, size = 1.2) +
  facet_grid2(
    . ~ Param,
    scales = "free_x",
    labeller = label_parsed,
    strip = strip_vanilla(clip = "off")
  ) +
  scale_y_continuous(breaks = breaks_extended(n = 6)) +
  scale_color_manual(name = "Cost", values = colors3) +
  labs(
    x = bquote(log[10] ~ "Initial Guess"),
    y = bquote(log[10] ~ "Optimal Value")
  ) +
  theme_cr() +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

# Same plot as previous but with the points colored and the plot faceted by the
# normalization options
initial_faceted <- initial |>
  ggplot(aes(x = Initial, y = Optimal, color = Cost3)) +
  geom_point(shape = 16, alpha = 0.2, size = 1.2) +
  facet_nested(
    Basis + NormCalc ~ Param,
    scales = "free_x",
    labeller = label_parsed,
    strip = strip_nested(clip = "off")
  ) +
  scale_y_continuous(breaks = breaks_extended(n = 6)) +
  scale_color_manual(name = "Cost", values = colors3) +
  labs(
    x = bquote(log[10] ~ "Initial Guess"),
    y = bquote(log[10] ~ "Optimal Value")
  ) +
  theme_cr() +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

# Combine the plots of all normalization methods with just the plots from the
# BS1 normalization into a single figure

# Create plot of initial guesses versus cost for normalization to BS1-data
initial_bs1 <- initial |>
  filter(Basis == "BS1", NormCalc == "Data") |>
  mutate(Facet = factor("bs1", levels = c("all", "bs1"))) |>
  bind_rows(initial) |>
  ggplot(aes(x = Initial, y = Optimal, color = Cost3, alpha = Facet)) +
  geom_point(shape = 16, size = 1.2) +
  facet_grid2(
    Facet ~ Param,
    scales = "free_x",
    labeller = labeller(
      Facet = c("all" = "All Norm", "bs1" = "BS1 Norm"),
      Param = label_parsed
    ),
    strip = strip_vanilla(clip = "off")
  ) +
  facetted_pos_scales(x = list(
    scale_x_continuous(breaks = seq(-8, -2, by = 3), limits = c(-9, -2)),
    scale_x_continuous(breaks = seq(-8, -2, by = 3), limits = c(-9, -2)),
    scale_x_continuous(breaks = seq(-12, -5, by = 3), limits = c(-12, -5)),
    scale_x_continuous(breaks = seq(-6, -2, by = 2), limits = c(-7, -2)),
    scale_x_continuous(breaks = seq(-6, -2, by = 2), limits = c(-7, -2))
  )) +
  scale_y_continuous(breaks = breaks_extended(n = 5)) +
  scale_alpha_manual(values = c(0.1, 0.25)) +
  scale_color_manual(name = "Cost", values = colors3) +
  labs(
    x = bquote(log[10] ~ "Initial Guess"),
    y = bquote(log[10] ~ "Optimal Value")
  ) +
  theme_cr(font_scale = 0.95) +
  guides(color = guide_legend(override.aes = list(alpha = 1)), alpha = "none")

# Save plots
crsave_nowarn(
  initial_optimal,
  path = here(paste0(figpath, "_initial-optimal.png")),
  ratio = 3.5,
  width = 8
)
crsave_nowarn(
  initial_faceted,
  path = here(paste0(figpath, "_initial-faceted.png")),
  ratio = 1.5,
  width = 8
)
crsave_nowarn(
  initial_bs1,
  path = here(paste0(figpath, "_initial-bs1.png")),
  ratio = 2.5,
  width = 8
)

# Cost for Different Normalizations ---------------------------------------

# Lollipop plot of the cost of each set of normalized values for the different
# normalization types
cost_lollipop <- output |>
  filter(Param == "kon6R") |>
  group_by(Basis, NormCalc) |>
  mutate(Rank = rank(Cost, ties.method = "first")) |>
  proper_names() |>
  ggplot(aes(
    x = Rank,
    y = Cost,
    color = interaction(Basis, NormCalc, lex.order = TRUE)
  )) +
  geom_point(size = 1) +
  geom_segment(
    aes(x = .data$Rank, xend = .data$Rank, y = 0, yend = .data$Cost),
    linewidth = 0.4
  ) +
  scale_color_manual(name = "", guide = "none", values = colors4) +
  scale_y_continuous(expand = expand1()) +
  facet_nested(Basis + NormCalc ~ ., strip = strip_nested(clip = "off")) +
  theme_cr()

# Cumulative distribution of the cost of the optimal values, creates curves of
# fraction of optimal sets versus cost showing what fraction of optimal sets
# from each normalization are under the specific cost
cost_curve <- output |>
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
  scale_color_manual(name = "Normalization", values = colors4) +
  labs(x = bquote(log[10] ~ "Cost")) +
  theme_cr(font_scale = 0.9)

# Save plots
crsave_nowarn(
  cost_lollipop,
  path = here(paste0(figpath, "_cost-lollipop.png")),
  ratio = 1,
  width = 6
)
crsave_nowarn(
  cost_curve,
  path = here(paste0(figpath, "_cost-curve.png")),
  ratio = 1.5,
  width = 5.5
)

# Cost for Parameter Groups -----------------------------------------------

# Group parameters by same cost and similar parameter values
# Select just the columns with parameter values, round the Cost to three decimal
# points to allow similar points to be grouped together, and convert the table
# to long form
group_cost <- output |>
  select(ID, Basis, NormCalc, Cost, Param, Optimal) |>
  mutate(Cost = signif(Cost, 3), Optimal = round(Optimal, 1)) |>
  group_by(Basis, NormCalc, Param, Optimal, Cost) |>
  summarize(Group = cur_group_id(), Count = n(), .groups = "drop_last") |>
  ungroup() |>
  mutate(Facet = factor("all", levels = c("all", "bs1")))

# Plot the parameter groups where the size of the points varies depending on
# group size
group_cost_points <- group_cost |>
  proper_names() |>
  ggplot(aes(x = Optimal, y = log10(Cost), size = Count, alpha = Count)) +
  geom_point(shape = 16, color = "gray20") +
  scale_alpha_continuous(range = c(0.25, 0.5), guide = "none") +
  scale_size_continuous(range = c(0.75, 6)) +
  scale_y_continuous(
    name = bquote(log[10] ~ "Cost"),
    expand = expansion(expand_range)
  ) +
  facet_nested(
    . ~ Param,
    scales = "free_x",
    labeller = label_parsed,
    strip = strip_vanilla(clip = "off")
  ) +
  labs(x = bquote(log[10] ~ "Parameter Value")) +
  facetted_pos_scales(x = list(
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-8, 0, by =  4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0))
  )) +
  theme_cr() +
  guides(color = "none")

# Plot the parameter groups where the size of the points varies depending on
# group size and facet/color based on normalization options
group_cost_facets <- group_cost |>
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
  scale_color_manual(values = colors4) +
  facet_nested(
    Basis + NormCalc ~ Param,
    scales = "free_x",
    labeller = label_parsed,
    strip = strip_nested(clip = "off")
  ) +
  labs(x = bquote(log[10] ~ "Parameter Value")) +
  facetted_pos_scales(x = list(
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-8, 0, by =  4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0))
  )) +
  theme_cr() +
  guides(color = "none")

# Combine the plots of all normalization methods with just the plots from the
# BS1 normalization into a single figure
group_cost_bs1_facets <- group_cost |>
  filter(Basis == "BS1", NormCalc == "Data") |>
  mutate(Facet = factor("bs1", levels = c("all", "bs1"))) |>
  bind_rows(group_cost) |>
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
  scale_color_manual(values = colors2) +
  facet_grid2(
    Facet ~ Param,
    scales = "free_x",
    labeller = labeller(
      Facet = c("all" = "All Norm", "bs1" = "BS1 Norm"),
      Param = label_parsed
    ),
    strip = strip_vanilla(clip = "off")
  ) +
  labs(x = bquote(log[10] ~ "Parameter Value")) +
  facetted_pos_scales(x = list(
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-8, 0, by =  4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0))
  )) +
  theme_cr() +
  guides(color = "none")

# Save plots
crsave_nowarn(
  group_cost_points,
  path = here(paste0(figpath, "_group-cost-points.png")),
  ratio = 4.25,
  width = 9
)
crsave_nowarn(
  group_cost_facets,
  path = here(paste0(figpath, "_group-cost-facets.png")),
  ratio = 2,
  width = 9
)

# Density of Optimal Values -----------------------------------------------

# Helper function to match values between data frames
# Inputs are data frame to search, tibble containing the columns and values to
# filter on (for `rowwise()`, can be obtained with `cur_group()`), column to
# search in filtered data frame, value to search for in filtered data frame, and
# column to return value from
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
  # Based on `https://stackoverflow.com/questions/69699009/how-to-extract-the-
  # density-value-from-ggplot-in-r`

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
  minval <- output |>
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

# Set the smoothing factor for the density distributions
adjust <- 1

# Calculate scaled density Y values for the minimal cost parameters
dens <- output |>
  density_y("Param", adjust = adjust) |>
  proper_names()
dens_norm <- output |>
  density_y(c("Basis", "NormCalc", "Param"), adjust = adjust) |>
  proper_names()

# Scaled density of the optimal parameter values - gives a visualization of how
# tight or widely spread the optimized values for each individual parameter are

# In `geom_label_repel()`, `color` is the text and border color, `fill` is the
# label fill (the extra digits at the end of the hex code make it slightly
# transparent), `nudge_x` moves each label along the X axis (same units as the
# data), `ylim` sets borders on the Y axis that the label won't cross, `size`
# sets the text size, `label.padding` sets the space between the label text and
# the surrounding box, `label.size` sets the label border width, and
# `min.segment.length` sets the shortest distance to draw a line between the
# label and the point
param_dens <- output |>
  proper_names() |>
  ggplot(aes(x = Optimal, y = after_stat(scaled))) +
  geom_density(
    adjust = adjust,
    linewidth = 1,
    color = pal_select("contrast3", colors = "green")) +
  geom_point(data = dens, aes(x = Optimal, y = y), color = "gray20") +
  geom_label_repel(
    data = dens,
    aes(x = Optimal, y = y, label = Label),
    color = "gray20",
    fill = "#FFFFFF99",
    family = "Roboto Regular",
    nudge_x = c(6, 6, -6, -6, -6, -6),
    ylim = c(0,0.95),
    size = 3.5,
    parse = TRUE,
    label.padding = 0.15,
    label.size = 0.3,
    min.segment.length = 1
  ) +
  scale_y_continuous(
    name = "Scaled Density",
    breaks = c(0, 0.5, 1),
    expand = expansion(mult = c(0, 0.1))
  ) +
  facet_grid2(
    . ~ Param,
    scales = "free_x",
    labeller = label_parsed,
    strip = strip_vanilla(clip = "off")
  ) +
  labs(x = bquote(log[10] ~ "Parameter Value")) +
  facetted_pos_scales(x = list(
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-8, 0, by =  4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0))
  )) +
  theme_cr()

# Same plot as previous but also colored and faceted by normalization options
# `force_pull` sets the attraction between the data labels and the points
param_dens_norm <- output |>
  proper_names() |>
  ggplot(aes(
    x = Optimal, y = after_stat(scaled),
    color = interaction(Basis, NormCalc, lex.order = TRUE)
  )) +
  geom_density(linewidth = 1, adjust = adjust) +
  geom_point(data = dens_norm, aes(x = Optimal, y = y), color = "gray20") +
  geom_label_repel(
    data = dens_norm,
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
  scale_color_manual(values = colors4) +
  facet_nested(
    Basis + NormCalc ~ Param,
    scales = "free_x",
    labeller = label_parsed,
    strip = strip_nested(clip = "off")
  ) +
  labs(x = bquote(log[10] ~ "Parameter Value")) +
  facetted_pos_scales(x = list(
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-8, 0, by =  4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0))
  )) +
  scale_y_continuous(name = "Scaled Density", breaks = c(0, 0.5, 1)) +
  theme_cr() +
  guides(color = "none")

# Same data as the previous plots but with the plot of all of the normalization
# methods combined into the same figure as just the BS1 normalization

# Set up the combined data frames for the plotting
output_facet <- output |>
  filter(Basis == "BS1", NormCalc == "Data") |>
  bind_rows(output, .id = "Facet") |>
  mutate(
    Facet = factor(Facet, levels = c("2", "1"), labels = c("all", "bs1"))
  )

dens_facet <- dens_norm |>
  filter(Basis == "BS1", NormCalc == "Data") |>
  bind_rows(dens, .id = "Facet") |>
  mutate(
    Facet = factor(Facet, levels = c("2", "1"), labels = c("all", "bs1"))
  )

param_dens_facet <- output_facet |>
  proper_names() |>
  ggplot(aes(x = Optimal, y = after_stat(scaled), color = Facet)) +
  geom_density(linewidth = 1, adjust = adjust) +
  geom_point(data = dens_facet, aes(x = Optimal, y = y), color = "gray20") +
  scale_color_manual(values = colors2) +
  facet_grid2(
    Facet ~ Param,
    scales = "free_x",
    labeller = labeller(
      Facet = c("all" = "All Norm", "bs1" = "BS1 Norm"),
      Param = label_parsed
    ),
    strip = strip_vanilla(clip = "off")
  ) +
  labs(x = bquote(log[10] ~ "Parameter Value")) +
  facetted_pos_scales(x = list(
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-8, 0, by =  4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0))
  )) +
  scale_y_continuous(
    name = "Scaled Density",
    breaks = c(0, 0.5, 1),
    expand = expansion(mult = c(0, 0.1))
  ) +
  theme_cr() +
  guides(color = "none")

# Save plots
crsave_nowarn(
  param_dens,
  path = here(paste0(figpath, "_param-dens.png")),
  ratio = 4.5,
  width = 9
)
crsave_nowarn(
  param_dens_norm,
  path = here(paste0(figpath, "_param-dens-norm.png")),
  ratio = 2,
  width = 9
)
crsave_nowarn(
  param_dens_facet,
  path = here(paste0(figpath, "_param-dens-bs1.png")),
  ratio = 3,
  width = 9
)

# Density and Cost Together -----------------------------------------------

# Combine the cost points plot and the scaled density plot for the normalization
# to the BS1 data into a single figure

# Make the grouped cost versus parameter value figure for the BS1 normalized
# data separately
group_cost_bs1 <- group_cost |>
  filter(Basis == "BS1", NormCalc == "Data") |>
  proper_names() |>
  ggplot(aes(x = Optimal, y = log10(Cost), size = Count, alpha = Count)) +
  geom_point(shape = 16, color = colors4[1]) +
  scale_alpha_continuous(range = c(0.2, 0.4), guide = "none") +
  scale_size_continuous(range = c(1, 6)) +
  scale_y_continuous(
    name = bquote(log[10] ~ "Cost"),
    expand = expansion(expand_range)
  ) +
  facet_nested(
    . ~ Param,
    scales = "free_x",
    labeller = label_parsed,
    strip = strip_vanilla(clip = "off")
  ) +
  facetted_pos_scales(x = list(
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0))
  )) +
  xlab(bquote(log[10] ~ "Parameter Value")) +
  theme_cr() +
  guides(color = "none")

# Filter the density data frame to just the BS1-data normalization
dens_norm_bs1 <- dens_norm |> filter(Basis == "BS1", NormCalc == "Data")

# Make the scaled parameter value density figure for the BS1 normalized data
# separately; first without the labeled parameter values
param_dens_bs1_nolabel <- output |>
  filter(Basis == "BS1", NormCalc == "Data") |>
  proper_names() |>
  ggplot(aes(x = Optimal, y = after_stat(scaled))) +
  geom_density(adjust = adjust, linewidth = 1, color = colors4[1]) +
  geom_point(data = dens_norm_bs1, aes(x = Optimal, y = y), color = "gray20") +
  scale_y_continuous(
    name = "Scaled Density",
    breaks = c(0, 0.5, 1),
    expand = expansion(mult = c(0, 0.1))
  ) +
  facet_grid2(
    . ~ Param,
    scales = "free_x",
    labeller = label_parsed,
    strip = strip_vanilla(clip = "off")
  ) +
  facetted_pos_scales(x = list(
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-15, -3, by = 6), limits = c(-15, -3)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0)),
    scale_x_continuous(breaks = seq(-8, 0, by = 4), limits = c(-8, 0))
  )) +
  xlab(bquote(log[10] ~ "Parameter Value")) +
  theme_cr()

# Add the labeled parameter values to the scaled density plots
param_dens_bs1 <- param_dens_bs1_nolabel +
  geom_label_repel(
    data = dens_norm_bs1,
    aes(x = Optimal, y = y, label = Label),
    color = "gray20",
    fill = "#FFFFFF99",
    family = "Roboto Regular",
    nudge_x = c(6, 6, -6, -6, -6, -6),
    ylim = c(0,0.95),
    size = 3.5,
    parse = TRUE,
    label.padding = 0.15,
    label.size = 0.3,
    min.segment.length = 1
  )

# Combine the figures together (with and without labels) and save
cost_dens_bs1 <- (group_cost_bs1 + theme_cut_bottom()) /
  (param_dens_bs1 + theme_cut_top() + plot_layout(tag_level = "new")) +
  plot_layout(guides = "collect", nrow = 2, heights = c(3, 2)) &
  theme(axis.title = element_text(size = 10))

cost_dens_bs1_nolabel <- (group_cost_bs1 + theme_cut_bottom()) /
  (param_dens_bs1_nolabel + theme_cut_top() + plot_layout(tag_level = "new")) +
  plot_layout(guides = "collect", nrow = 2, heights = c(3, 2))

crsave_nowarn(
  cost_dens_bs1,
  path = here(paste0(figpath, "_cost-dens-bs1.png")),
  ratio = 3,
  width = 9
)

crsave_nowarn(
  cost_dens_bs1_nolabel,
  path = here(paste0(figpath, "_cost-dens-bs1_nolabel.png")),
  ratio = 3,
  width = 9
)
