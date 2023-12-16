
# SOURCE FUNCTIONS --------------------------------------------------------

# Import helper functions to relabel the data values
source(here::here("code/data_processing/relabel_data.R"))
# Specifically imports `short_names()` which is my standard function for
# renaming simulation variable labels (like the antibody names, parameter
# labels, etc.)

# Import the occupancy calculation for the `save_occupancy()` function
source(here::here("code/data_processing/calculations.R"))

# HELPER FUNCTIONS --------------------------------------------------------

# Wrapper for `dplyr::inner_join` with my common settings applied
# relationship = "many-to-many" makes the join expect multiple values in `x` to
# match multiple values in `y`, which is typically the case when I am joining
# two tables by the ID values and am using the join to expand out the data into
# a single table
# unmatched = "error" makes the join report an error if any columns are dropped
# instead of silently dropping them; works because I am using filter steps to
# remove unwanted data, and I expect my join steps to not drop any columns
join_many <- function(df1, df2, by, ...) {
  dplyr::inner_join(
    df1, df2, by,
    relationship = "many-to-many",
    unmatched = "error",
    ...
  )
}

# Determine all of the file names that match a given prefix; returns a list
# where `$files` is all of the file names and `$type` is the simulation type
find_files <- function(prefix, type = NULL){

  # If type is not provided, attempt to determine the type from the file name
  if (is.null(type)) {
    # Extract the portion of the file name directly after the first underscore
    type <- stringr::str_extract(prefix, "(?<=_)[^_]*(?=_?)")

    # Set the type based on the extracted portion of the file name
    type <- switch(
      type,
      "opt" = "opt",
      "sim" = "sim",
      "base" = "sim",
      "sim-opt" = "sim",
      "univariate" = "sim",
      "Binding" = "exp"
    )
  }

  # Set file path based on file type
  if (type == "sim" || type == "opt") {
    folder <- "output"
  } else if (type == "exp") {
    folder <- "data"
  }

  # Make regex pattern to match the file names and find matching files
  # Including `.rds` files for importing files saved from R (e.g., occupancy)
  pattern <- paste0(prefix, ".*.(csv|rds)$")
  files <- list.files(path = here::here(folder), pattern = pattern)

  # If there are no files that match the regex pattern, throw an error
  if (length(files) == 0) {
    stop("No files matched the input file prefix.")
  }

  # Create the output list
  list(files = files, type = type)
}

# Helper function to set all of the column types corresponding to the different
# column names; accepts a default column type for any column names that aren't
# present in the formats list, by default, uses "?" to guess the type for any
# columns that aren't in the list
col_formats <- function(col_names, default = "?") {

  # Inputs for the `col_types` argument of `read_csv()`
  formats <- list(
    ID = "f",
    ConcID = "f",
    ParamID = "f",
    OutputID = "f",
    Ab = "f",
    Antibody = "f",
    Recep = "f",
    Cell = "f",
    Species = "f",
    Param = "f",
    Period = "f",
    Interest = "f",
    Basis = "f",
    Calc = "f",
    Occupied = "f",
    Add = "l",
    Conc = "d",
    Start = "d",
    End = "d",
    Time = "d",
    Value = "d",
    Bound = "d",
    Frac = "d"
  )

  # Determine which column names aren't present in the column types list, and
  # set those column types to the default type
  missing <- col_names[!col_names %in% names(formats)]
  formats[missing] <- default

  # Return the column type for each column name in the list
  formats[col_names]
}

# Merge the lists of data frames from separate files that are generated from
# mapping `import_files()` and `clean_data()` over a list of input prefixes
# Outputs a single list with separate data frames for each clean data table,
# in the same format as the output from the import functions when a single
# prefix is used
merge_data <- function(df_list) {

  # If input is not a nested list, nest the data frames inside a list to match
  # correct syntax
  # Occurs for "exp" and "opt" type data because those functions yield a single
  # data frame of clean data
  # Necessary form: level 1 = separate file prefixes, level 2 = separate data
  # frames from each prefix
  if ("tbl" %in% class(df_list[[1]])) {
    df_list <- purrr::map(df_list, \(x) list(data = x))
  }

  # Store the original order of the data frames in the list
  orig <- names(df_list[[1]])

  # Merge the data frames in the list
  data <- df_list |>
    # Sort the list of data tables in order by name because `pmap` iterates
    # over list elements by position
    purrr::map(\(x) x[order(names(x))]) |>
    # Row-bind the data frames of each list by position (i.e., df_list[[1]][[1]]
    # will bind with df_list[[2]][[1]], df_list[[3]][[1]], etc.)
    purrr::pmap(\(...) dplyr::bind_rows(..., .id = "FileID")) |>
    # Convert the "FileID" column to a factor
    purrr::map(\(x) dplyr::mutate(x, FileID = factor(FileID)))

  # Restore the original data frame order
  data <- data[orig]

  # Unlist outputs that only have a single data frame to restore the original
  # output format on data frames that had to be nested in a list
  if (length(data) == 1) {
    data <- data[[1]]
  }

  return(data)
}

# DATA IMPORT -------------------------------------------------------------

# Utility function to import the simulation output files and standardize the
# table format
import_sim <- function(file) {

  # Setup -------------------------------------------------------------------
  # Setup variables needed for data organization

  # Columns in the data frames that are used for simulation identification
  # instead of holding data
  id_cols <- c(
    "ID", "ConcID", "ParamID", "OutputID",
    "Ab", "Recep",
    "Period", "Start", "End",
    "Species", "Cell", "Param", "Basis", "Calc", "Add",
    "Conc", "Time")

  # Import ------------------------------------------------------------------
  # Import the data file as a tibble

  # If the file to be imported is an `.rds` file (e.g., occupancy), import the
  # file and exit the function
  if (stringr::str_detect(file, "\\.rds$")) {
    sim <- readRDS(file = here::here("output", file))
    return(sim)
  }

  # Determine column names in the imported file
  col_names <- names(readr::spec_csv(here::here("output", file))$cols)

  # Read in the data using the correct string for `col_types`
  sim <- readr::read_csv(
    file = here::here("output", file),
    col_types = col_formats(col_names, default = "?")
  )

  # Organization ------------------------------------------------------------
  # Organize the data into a standard format for plotting

  # Remove duplicated time points
  if ("Time" %in% names(sim)) {
    sim <- sim |>
      dplyr::distinct(
        dplyr::across(tidyselect::any_of(id_cols)),
        .keep_all = TRUE
      ) |>
      dplyr::mutate(Time = .data$Time / 3600)
  }

  # Rename non-existent species
  # Can occur where "ternary" was one of the outputs but no ternary complexes
  # were possible (e.g., BS1 in 6R+ cells) or where antibody could not bind to a
  # receptor type (e.g., Toci in 8R+ cells)
  # To handle NA, first explicitly change NA to a factor level so fct_collapse
  # recognizes it, then change it to the desired factor
  # Also need to replace "*2+" because it is an artifact from making complexes
  # when binding is not possible (e.g., Toci-8R-8R); first have to add it to the
  # factor levels so fct_collapse doesn't throw an "Unknown levels" warning
  # It is also possible for "*2+" to appear at the end of a name (e.g., BS1 in
  # 6R+ cells where no ternary complexes formed), so the final step strips "+2*"
  # in those cases
  # Moving the "None" factor to the end so data that was saved prior to the
  # correction to the MATLAB species names will match sorting of the newer data
  if ("Species" %in% names(sim)) {
    sim <- sim |>
      dplyr::mutate(
        Species = forcats::fct_na_value_to_level(Species),
        Species = forcats::fct_expand(Species, "+2*"),
        Species = forcats::fct_collapse(Species, None = c(NA, "+2*")),
        Species = forcats::fct_relabel(
          Species,
          ~ stringr::str_replace(.x, "\\+2\\*$", "")
        ),
        Species = forcats::fct_relevel(Species, "None", after = Inf)
      )
  }

  return(sim)
}

# Utility function to import the experimental data files and standardize the
# table format
import_exp <- function(file) {

  # Setup -------------------------------------------------------------------
  # Setup variables needed for data organization

  # Antibodies and cell lines present in the experimental data
  ab <- c("Toci", "H2", "BS1", "BS2")
  cells <- c("6R", "8R", "6R8R", "Neg")

  # Import ------------------------------------------------------------------
  # Import the data file as a tibble

  # Import the data file as a tibble
  exp <- readr::read_csv(
    file = here::here("data", file),
    show_col_types = FALSE
  )

  # Organization ------------------------------------------------------------
  # Organize the data into a standard format for plotting

  # Select only columns that do not start with "h" to drop the data from the
  # humanized antibodies
  exp <- exp |>
    dplyr::select(!tidyselect::starts_with("hBS"))

  # Convert the experimental data to long form, and separate the antibodies,
  # cell lines, and replicates into separate columns
  exp <- exp |>
    tidyr::pivot_longer(
      cols = !"Conc",
      names_to = c("Species", "Cell", "Rep"),
      names_sep = "_",
      names_transform = list(
        Species = ~ readr::parse_factor(.x, levels = ab),
        Cell = ~ readr::parse_factor(.x, levels = cells),
        Rep = ~ readr::parse_factor(.x, levels = c("1", "2", "3"))
      ),
      values_to = "Value"
    )

  # Determine the normalization type from the file name
  norm_type <- stringr::str_extract(file, "Ab|BS1(?=-Norm)")

  # Set the norm type to BS1 if there was no information in the file name
  if (is.na(norm_type)) {
    norm_type <- "BS1"
  }

  # Add a column for the normalization type
  exp <- exp |>
    dplyr::mutate(Basis = factor(norm_type, levels = c("Ab", "BS1")))

  return(exp)
}

# Utility function to import the optimization results and standardize the table
# format
import_opt <- function(file) {

  # Setup -------------------------------------------------------------------
  # Setup variables needed for data organization

  # Inputs for the `col_types` argument of `read_csv()`
  formats <- c(
    initial = "fdddddddd",
    optimal = "fdddddddd",
    results = "fdfd"
  )

  # Optimization types performed
  opt_types <- c("on", "on-off", "on-off2", "koff3", "koff4", "koff5",
                 "koff6", "koff7")

  # Import ------------------------------------------------------------------
  # Import the data file as a tibble

  # Determine the variable being imported from the file name
  suffix <- gsub(".*_([^_].*).csv$", "\\1", file)
  suffix <- gsub("-(.*)$", "", suffix)

  # Read in the data using the correct string for `col_types`
  data <- readr::read_csv(
    file = here::here("output", file),
    col_types = formats[[suffix]]
  )

  # Organization ------------------------------------------------------------
  # Organize the data into a standard format for plotting

  # Replace the parameter regular expressions with simpler names
  # arrange() after the renaming sorts the parameters in the correct order
  # before converting them back to columns
  prefix <- c("initial" = "Initial.", "optimal" = "Param.")
  if (suffix == "initial" || suffix == "optimal") {
    data <- data |>
      tidyr::pivot_longer(
        cols = !"ID",
        names_to = "Param",
        values_to = "Value",
        names_transform = ~ readr::parse_factor(.x)
      ) |>
      short_names() |>
      dplyr::arrange(.data$ID, .data$Param) |>
      tidyr::pivot_wider(
        names_from = "Param",
        values_from = "Value",
        names_prefix = prefix[[suffix]]
      )
  }

  # Classify Data -----------------------------------------------------------
  # Set the optimization type, normalization basis, and normalization
  # calculation for each import file

  # Determine the optimization parameters from the file name - looking for the
  # substring before the final underscore in the file name
  opt_params <- stringr::str_extract(file, "[^_]*(?=_[^_]*$)")

  # Classify the optimization type
  if (stringr::str_detect(opt_params, "kon-koff2")) {
    opt_type <- "on-off2"
  } else if (stringr::str_detect(opt_params, "koff[0-9]$")) {
    opt_type <- stringr::str_extract(opt_params, "koff[0-9]$")
  } else if (stringr::str_detect(opt_params, "^kon-koff")) {
    opt_type <- "on-off"
  } else {
    opt_type <- "on-off2"
  }

  # Classify the normalization basis
  opt_basis <- stringr::str_extract(opt_params, "(ab|bs1)")

  # Classify the normalization calculation
  opt_calc <- ifelse(stringr::str_detect(opt_params, "max$"), "max", "data")

  # Add columns for the optimization parameters to the imported tables
  data <- data |>
    dplyr::mutate(
      Opt = factor(opt_type, levels = opt_types),
      Basis = factor(opt_basis, levels = c("ab", "bs1")),
      Calc = factor(opt_calc, levels = c("data", "max"))
    ) |>
    dplyr::select(
      "ID",
      "Opt",
      "Basis",
      "Calc",
      tidyselect::everything()
    )

  return(data)
}

# Wrapper function to import all files with a particular prefix
import_files <- function(prefix, type = NULL) {

  # Determine the list of files to import and the file type
  imports <- find_files(prefix = prefix, type = type)
  files <- imports$files
  type <- imports$type

  # Import the data files based on data type
  if (type == "sim") {

    # Import data tibbles into a list
    data <- lapply(files, import_sim)

    # Name each table based on the variable name in MATLAB
    suffix <- gsub(".*_([^_].*).(csv|rds)$", "\\1", files)
    names(data) <- suffix

  } else if (type == "exp") {

    # Import data tibbles into a list
    data <- lapply(files, import_exp)

    # Concatenate the tibbles into a single tibble
    data <- do.call(rbind, data)

  } else if (type == "opt") {

    # Determine the variables being imported from the file names
    suffix <- gsub(".*_([^_].*).csv$", "\\1", files)
    suffix <- gsub("-(.*)$", "", suffix)
    variables <- unique(suffix)

    # Do not import the saved parameters since they are saved from the R
    # analysis, not from the MATLAB optimizations
    variables <- variables[variables != "params"]

    # Import all of the files for each variable separately
    # This is functionally equivalent to a nested for loop. The inner `map`
    # applies the `import_opt()` function to all of the files of a particular
    # variable (e.g., "optimal"). The `reduce` step binds all of the tables for
    # that variable into a single table. The outer `map` applies the inner `map`
    # and `reduce` steps to every variable in the file set. The end result is a
    # list with a single table for each imported variable
    data <- purrr::map(
      variables,
      ~ purrr::reduce(
        purrr::map(files[suffix == .], import_opt),
        dplyr::bind_rows
      )
    )

    # Name each table based on the variable name in MATLAB
    names(data) <- variables
  }
  return(data)
}

# DATA ORGANIZATION -------------------------------------------------------

# Utility function to make tibbles with the variables used in each simulation
clean_sim <- function(sims) {

  # Setup -------------------------------------------------------------------
  # Setup variables needed for data organization

  # Define names of antibodies, receptors, and cell lines for matching
  ab_names <- c("Toci", "H2", "BS1", "Toci_H2")
  recep_names <- c("IL6R", "IL8R")

  # Initialize output list
  output <- list()

  # Initial Concentrations --------------------------------------------------
  # Determine the initial antibody concentration in each simulation
  # The normalization function has to round the total antibody concentrations to
  # compare values without floating-point error in MATLAB, so they also need to
  # be rounded here
  output$antibodies <-
    join_many(sims$id, sims$yin, by = "ConcID") |>
    dplyr::filter(
      .data$Species %in% ab_names,
      .data$Period == "initial"
    ) |>
    dplyr::select(tidyselect::all_of(
      c("ID", "ConcID", "Ab", "Species", "Conc")
    )) |>
    dplyr::distinct() |>
    tidyr::pivot_wider(
      names_from = "Species",
      values_from = "Conc",
      names_prefix = "Ab.",
      values_fill = 0
    ) |>
    dplyr::rename(Antibody = "Ab") |>
    short_names() |>
    dplyr::mutate(
      Ab.Total = rowSums(dplyr::across(tidyselect::ends_with(ab_names))),
      Ab.Total = signif(.data$Ab.Total, digits = 12)
    )

  # Receptors ---------------------------------------------------------------
  # Determine the amount of receptor present in different cell lines
  output$receptors <-
    join_many(sims$id, sims$yin, by = "ConcID") |>
    dplyr::filter(
      .data$Species %in% recep_names,
      .data$Period == "initial"
    ) |>
    dplyr::select(tidyselect::all_of(
      c("ID", "ConcID", "Recep", "Species", "Conc")
    )) |>
    dplyr::distinct() |>
    tidyr::pivot_wider(
      names_from = "Species",
      values_from = "Conc",
      names_prefix = "Recep.",
      values_fill = 0
    ) |>
    dplyr::rename(Cell = "Recep") |>
    short_names() |>
    dplyr::mutate(
      Recep.Total = .data$Recep.IL6R + .data$Recep.IL8R
    )

  # Parameters --------------------------------------------------------------
  # Make a tibble with the parameters used for each simulation
  output$parameters <-
    join_many(sims$id, sims$params, by = "ParamID") |>
    dplyr::select(tidyselect::all_of(c("ID", "ParamID", "Param", "Value"))) |>
    dplyr::distinct() |>
    short_names() |>
    tidyr::pivot_wider(
      names_from = "Param",
      values_from = "Value",
      names_prefix = "Param."
    )

  # Output Calculations -----------------------------------------------------
  # Determine the output calculations used for each output value
  output_cut <- sims$output |>
    dplyr::select(!tidyselect::all_of(c("Time")))
  out_cut <- sims$out |>
    dplyr::select(tidyselect::all_of(c("ID", "OutputID")))

  output$calculations <-
    dplyr::inner_join(out_cut, output_cut, by = "OutputID") |>
    dplyr::distinct()

  # Output Values -----------------------------------------------------------
  # Join all of the simulation tables onto the output table
  output$output <-
    dplyr::inner_join(sims$out, output$antibodies, by = "ID") |>
    dplyr::inner_join(output$receptors, by = c("ID", "ConcID")) |>
    dplyr::inner_join(output$parameters, by = "ID") |>
    dplyr::inner_join(output$calculations, by = c("ID", "OutputID", "Calc")) |>
    dplyr::select(
      "ID",
      "ConcID",
      "ParamID",
      "OutputID",
      "Antibody",
      "Cell",
      tidyselect::starts_with("Ab.", ignore.case = FALSE),
      tidyselect::starts_with("Recep.", ignore.case = FALSE),
      tidyselect::starts_with("Param.", ignore.case = FALSE),
      "Interest",
      "Calc",
      "Add",
      "Species",
      "Time",
      "Value"
    )

  # Join the relevant simulation tables onto the normalized output table if it
  # is given in the results
  if ("norm" %in% names(sims)) {

    # Add a column for OutputID if it is not already present in the data,
    # necessary for compatibility with data saved prior to 2023-03-07
    # Find which OutputID corresponds to the output used to normalize the output
    # (the total bound antibody concentration), and then perform a cross join to
    # add that OutputID as a column in the normalized output table
    if (!"OutputID" %in% names(sims$norm)) {
      sims$norm <- output$calculations |>
        dplyr::distinct(dplyr::pick(!"ID")) |>
        dplyr::filter(
          .data$Interest == "receptors",
          .data$Add == TRUE,
          .data$Calc %in% c("conc", "conct")
        ) |>
        dplyr::select("OutputID") |>
        dplyr::cross_join(sims$norm)
    }

    output$norm <- sims$norm |>
      dplyr::rename(
        "Ab.Total" = "Conc",
        "Antibody" = "Species",
        "NormCalc" = "Calc"
      ) |>
      short_names() |>
      dplyr::mutate(Ab.Total = signif(.data$Ab.Total, digits = 12)) |>
      dplyr::inner_join(
        output$antibodies,
        by = c("ID", "Antibody", "Ab.Total")
      ) |>
      dplyr::inner_join(output$receptors, by = c("ID", "ConcID", "Cell")) |>
      dplyr::inner_join(output$parameters, by = "ID") |>
      dplyr::left_join(output$calculations, by = c("ID", "OutputID")) |>
      dplyr::select(
        "ID",
        "ConcID",
        "ParamID",
        "OutputID",
        "Antibody",
        "Cell",
        tidyselect::starts_with("Ab.", ignore.case = FALSE),
        tidyselect::starts_with("Recep.", ignore.case = FALSE),
        tidyselect::starts_with("Param.", ignore.case = FALSE),
        "Interest",
        "Calc",
        "Add",
        "Basis",
        "NormCalc",
        "Time",
        "Value"
      ) |>
      dplyr::mutate(Value = replace(.data$Value, is.nan(.data$Value), 0))
  }

  # Join the relevant simulation tables onto the cost table if it is given
  if ("cost" %in% names(sims)) {
    output$cost <- sims$cost |>
      dplyr::inner_join(sims$id, by = "ParamID") |>
      dplyr::inner_join(output$parameters, by = c("ID", "ParamID")) |>
      dplyr::select(
        "ID",
        "ParamID",
        tidyselect::starts_with("Param.", ignore.case = FALSE),
        "Value"
      ) |>
      dplyr::distinct() |>
      dplyr::arrange(.data$ID)
  }

  # Include the occupancy data in the output if it was imported
  if ("occupied" %in% names(sims)) {
    output$occupied <- sims$occupied
  }

  return(output)
}

# Utility function make the experimental data into a tibble matching the
# simulation output
clean_exp <- function(exp) {

  # Experimental Data -------------------------------------------------------
  # Clean experimental data columns

  # Convert concentration to nM, convert normalized values to decimals, reorder
  # the table columns, and sort the table rows
  exp <- exp |>
    dplyr::mutate(Conc = 10^.data$Conc * 10^9, Value = .data$Value/100) |>
    dplyr::rename(Antibody = "Species") |>
    dplyr::select(c("Antibody", "Cell", "Conc", "Rep", "Basis", "Value")) |>
    short_names() |>
    dplyr::arrange(
      .data$Basis,
      .data$Antibody,
      .data$Cell,
      .data$Conc,
      .data$Rep
    )

  return(exp)
}

# Utility function to combine all of the optimization variables and results into
# a single table
clean_opt <- function(opt) {

  # Combine Tables ----------------------------------------------------------
  # Combine data tables into a single table
  output <-
    dplyr::inner_join(
      opt$initial,
      opt$optimal,
      by = c("ID", "Opt", "Basis", "Calc")
    ) |>
    dplyr::inner_join(opt$results, by = c("ID", "Opt", "Basis", "Calc"))

  return(output)
}

# Wrapper function to take in the imported data and the type of data and apply
# the correct cleaning function to the data
clean_data <- function(data, type = NULL) {

  # If type is not provided, attempt to determine the type from the data
  if (is.null(type)) {
    # Determine the names of the lists/columns in the data
    data_names <- names(data)

    # Set the type based on the lists/columns present
    if ("initial" %in% data_names) {
      type <- "opt"
    } else if ("id" %in% data_names) {
      type <- "sim"
    } else if ("Conc" %in% data_names) {
      type <- "exp"
    } else {
      stop("`type` was not given and could not be determined from the data.")
    }
  }

  # Clean the imported data based on which type it is
  if (type == "sim") {
    data <- clean_sim(data)
  } else if (type == "exp") {
    data <- clean_exp(data)
  } else if (type == "opt") {
    data <- clean_opt(data)
  }

  return(data)
}

# MERGE MULTIPLE FILES ----------------------------------------------------

# Imports data from each of the file prefixes, cleans the data, and merges the
# data from each separate prefix into a single list of data frames
import_data <- function(prefixes, type = NULL) {

  # Import data from each given file prefix
  # Uses `map()` to loop through multiple files if supplied as inputs
  data <- purrr::map(prefixes, \(x) import_files(x, type = type))

  # Iterate through each data set and clean data
  data <- purrr::map(data, \(x) clean_data(x, type = type))

  # Merge the separate data from each prefix into a single list
  data <- merge_data(data)
  return(data)
}

# SAVE PROCESSED DATA -----------------------------------------------------

# Save the calculated fractional occupancy data to avoid repeating slow
# processing steps
# `include` argument dictates which molecule types should be used as the bases
# for the occupancy calculation; "all" uses all individual and summed species,
# "totals" uses only the summed antibodies and receptors, and "species" uses
# only the individual species as the bases
save_occupancy <- function(prefix, type = NULL, include = "all") {

  # Determine the file name to save the data under
  # Determine the names of all of the files that were imported
  imports <- find_files(prefix = prefix, type = type)

  # Select just the portion of the file name that does not correspond to the
  # saved variable
  basename <- stringr::str_extract(
    imports$files,
    "(.*)_.*\\.(csv|rds)$",
    group = 1
  )

  # Throw an error if the basename is not unique
  if (length(unique(basename)) > 1) {
    stop("Files from more than one simulation were provided.")
  }

  # Construct the file name from the base name for the simulation
  fn <- paste0(unique(basename), "_occupied.rds")

  # Check if occupancy file exists and only generate new file if it doesn't
  file_match <- list.files(
    path = here::here("output"),
    pattern = paste0("^", fn, "$")
  )

  if (length(file_match) == 0) {

    # Import and clean the data, then calculate the fractional occupancy of the
    # processed data
    data <- import_data(prefixes = prefix, type = type)
    occupied <- occupancy(data$output, include = include)

    # Save the occupancy data in the constructed file name
    saveRDS(occupied, here::here("output", fn))

  } else {
    occupied <- list()
  }

  # Return occupancy data frame if requested
  invisible(occupied)
}
