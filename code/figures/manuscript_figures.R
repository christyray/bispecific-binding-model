here::i_am("code/figures/manuscript_figures.R")

# SOURCE FUNCTIONS --------------------------------------------------------

# Import helper functions to import and organize the simulation data
source(here::here("code/data_processing/initialization.R"))

# Import plotting functions that are shared between multiple figures
source(here::here("code/figures/common_functions.R"))

# LOAD PACKAGES -----------------------------------------------------------

library(tidyverse)
library(cli)
library(scales)
library(crthemes)
library(sciscales)
library(ggh4x)
library(paletteer)
library(ggthemes)
library(patchwork)

# USER OPTIONS ------------------------------------------------------------

# saveqQ = option to save the generated figures
# posterQ = option to save the poster-sized versions of the figures
# setupQ = option to only run the package and function loading sections

saveQ <- TRUE
posterQ <- FALSE
setupQ <- FALSE

# Set file name for naming output figures
manupath <- "output/figures/manuscript/"
posterpath <- "output/figures/presentations/"

# Set resolution of saved images
dpi <- 150

if (posterQ) {
  fig_path <- posterpath
} else {
  fig_path <- manupath
}

# Setting which figure functions to run
figures <- c(
  "parameter_optimization" = TRUE,
  "binding_curve" = TRUE,
  "basic_simulation" = TRUE,
  "vary_concentrations" = TRUE,
  "monovalent_only" = TRUE,
  "comparing_antibodies" = TRUE,
  "univariate_sensitivity" = TRUE,
  "comparing_optimized" = TRUE
)

# If running in "setup mode", print message and stop script
if (setupQ) {
  # `suppress_warn` is a variable inside the plot functions that is not defined
  # elsewhere, so it is helpful to pre-define it before testing code inside the
  # functions
  suppress_warn <- TRUE

  rlang::abort("Figure generation setup complete!")
}

# Print Script Options ----------------------------------------------------

# Print warning that manuscript and poster versions cannot be generated together
if (saveQ && posterQ) {
  cli_abort(c(
    "Manuscript and poster figures must be created separately.",
    "x" = "{.code saveQ} and {.code posterQ} are both {.strong TRUE}."
  ))
}

# Set the display names of the figure functions
figure_names <- c(
  "Parameter Optimization",
  "Binding Curve",
  "Basic Simulation",
  "Vary Concentrations",
  "Monovalent Only",
  "Comparing Antibodies",
  "Univariate Sensitivity",
  "Supplemental All Parameter Sets"
)

# Print summary of options before script is run
cli_h1("Script Options")
cli_alert_info("Project path: {.file {here()}}")

if (saveQ) {
  cli_alert_success("Save manuscript figures: {.strong {saveQ}}")
  cli_alert_info("Manuscript figure path: {.file {manupath}}")
} else {
  cli_alert_warning("Save manuscript figures: {.strong {saveQ}}")
}
if (posterQ) {
  cli_alert_success("Save poster figures: {.strong {posterQ}}")
  cli_alert_info("Poster figure path: {.file {posterpath}}")
} else {
  cli_alert_warning("Save poster figures: {.strong {posterQ}}")
}

if (saveQ || posterQ) {
  cli_h1("Figure Selection")
  cli_ul(figure_names[figures])
}

cli_h1("Generating Figures")

# PARAMETER OPTIMIZATION --------------------------------------------------

# Run the script with the figure generation functions
source(here::here("code/figures/parameter_optimization.R"))

# Set the data file to use for the figures and the figure name
data_file <- "optimization"
fig_file <- "Parameter-Optimization"

# Generate and save the figures
if (figures["parameter_optimization"]) {
  parameter_optimization(data_file, fig_file, fig_path, dpi = dpi)
} else {
  parameter_optimization()
}

# BINDING CURVE -----------------------------------------------------------

# Run the script with the figure generation functions
source(here::here("code/figures/binding_curve.R"))

# Set the data file to use for the figures and the figure name
data_file <- "binding-curve"
fig_file <- "Binding-Curve"

# Generate and save the figures
if (figures["binding_curve"]) {
  binding_curve(data_file, fig_file, fig_path, poster = posterQ, dpi = dpi)
} else {
  binding_curve()
}

# BASIC SIMULATION --------------------------------------------------------

# Run the script with the figure generation functions
source(here::here("code/figures/basic_simulation.R"))

# Set the data file to use for the figures and the figure name
data_file <- "time"
fig_file <- "Basic-Simulation"

# Calculate the fractional occupancy if it doesn't exist
walk(data_file, \(x) save_occupancy(x, type = "sim"))

# Generate and save the figures
if (figures["basic_simulation"]) {
  basic_simulation(data_file, fig_file, fig_path, dpi = dpi)
} else {
  basic_simulation()
}

# VARY CONCENTRATIONS -----------------------------------------------------

# Run the script with the figure generation functions
source(here::here("code/figures/vary_concentrations.R"))

# Set the data file(s) to use for the figures and the figure name
data_file <- "concentration"
fig_file <- "Vary-Concentrations"

# Calculate the fractional occupancy if it doesn't exist
walk(data_file, \(x) save_occupancy(x, type = "sim"))

# Generate and save the figures
if (figures["vary_concentrations"]) {
  vary_concentrations(data_file, fig_file, fig_path,
                      poster = posterQ, dpi = dpi)
} else {
  vary_concentrations()
}

# MONOVALENT ONLY ---------------------------------------------------------

# Run the script with the figure generation functions
source(here::here("code/figures/monovalent_only.R"))

# Set the data file(s) to use for the figures and the figure name
data_file <- "monovalent"
fig_file <- "Monovalent-Only"

# Calculate the fractional occupancy if it doesn't exist
walk(data_file, \(x) save_occupancy(x, type = "sim"))

# Generate and save the figures
if (figures["monovalent_only"]) {
  monovalent_only(data_file, fig_file, fig_path, dpi = dpi)
} else {
  monovalent_only()
}

# COMPARING ANTIBODIES ----------------------------------------------------

# Run the script with the figure generation functions
source(here::here("code/figures/comparing_antibodies.R"))

# Set the data file(s) to use for the figures and the figure name
data_file <- c(
  "compare-recep",
  "compare-ab"
)
fig_file <- "Comparing-Antibodies"

# Calculate the fractional occupancy if it doesn't exist
walk(data_file, \(x) save_occupancy(x, type = "sim"))

# Generate and save the figures
if (figures["comparing_antibodies"]) {
  comparing_antibodies(data_file, fig_file, fig_path,
                       poster = posterQ, dpi = dpi)
} else {
  comparing_antibodies()
}

# UNIVARIATE SENSITIVITY --------------------------------------------------

# Run the script with the figure generation functions
source(here::here("code/figures/univariate_sensitivity.R"))

# Set the data file(s) to use for the figures and the figure name
data_file <- list(
  local = "local",
  global = "global"
)
fig_file <- "Univariate-Sensitivity"

# Calculate the fractional occupancy if it doesn't exist
walk(data_file, \(x) save_occupancy(x, type = "sim"))

# Generate and save the figures
if (figures["univariate_sensitivity"]) {
  univariate_sensitivity(data_file, fig_file, fig_path,
                         poster = posterQ, dpi = dpi)
} else {
  univariate_sensitivity()
}

# SUPPLEMENTAL OPTIMIZED FIT ----------------------------------------------

# Run the script with the figure generation functions
source(here::here("code/figures/comparing_optimized.R"))

# Set the data file(s) to use for the figures and the figure name
prefix <- "compare-opt-"
norms <- c("bs1-data", "bs1-max", "ab-data", "ab-max")

# Using `set_names()` allows the norm types to passed in alongside the data file
# paths
norms <- set_names(norms) |> map_chr(\(x) paste0(prefix, x))

data_file <- list(norms = norms, experimental = "Flow-MFI_")
fig_file <- "Binding-Curve-All-Sets"

# Generate and save the figures
if (figures["comparing_optimized"]) {
  comparing_optimized(data_file, fig_file, fig_path, dpi = dpi)
} else {
  comparing_optimized()
}

# Print final message
cli_progress_step("All figures complete!")
