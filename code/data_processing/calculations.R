
# CALCULATIONS ------------------------------------------------------------

# General -----------------------------------------------------------------

# Generate a list with all possible antibodies, receptors, and antibody-
# receptor combinations, sorted in the standard ordering; primarily used for
# generating the factor levels for output data frames
# Takes input of the species present in the data (e.g., from
# `levels(data$Species)`) and the maximum number of receptors to include in the
# complexes
create_molecules <- function(species, maxR = 2) {

  # Extract just the antibody name from the species and convert to factor to
  # preserve the ordering
  antibody <- unique(stringr::str_extract(species, "^[^_]*(?=_)"))
  antibody <- antibody[!is.na(antibody)]
  antibody <- factor(antibody, levels = antibody)

  # Extract the receptors present from the species names
  # The `_` in `do.call()` refers to the input from the pipe; using `c()` to
  # join the lists
  receptor <- species |>
    stringr::str_extract_all("(?<=_)[0-9]R") |>
    do.call(c, args = _) |>
    unique()

  # Create all possible permutations of bound complexes; outputs a list where
  # each element is a data frame of all receptor combinations for that number of
  # receptors (e.g., [[1]] = one receptor, [[2]] = two receptors)
  complexes <- lapply(
    1:maxR,
    function(x) as.data.frame(gtools::permutations(
      length(receptor),
      x,
      receptor,
      repeats.allowed = TRUE
    ))
  )

  # Combine the separate data frames from lapply into a single table and replace
  # the NA values with empty strings
  # Then, duplicate the complexes for each antibody in the set and create the
  # complex names by joining the antibody to the receptors
  complexes <- complexes |>
    do.call(dplyr::bind_rows, args = _) |>
    tidyr::unite(Complex, sep = "_", na.rm = TRUE) |>
    # `uncount()` duplicates each row N times, where N = `length(antibody)`
    tidyr::uncount(length(antibody)) |>
    # Add a column with the antibodies, using grouping to match each antibody
    # once with each possible complex type
    dplyr::mutate(Antibody = antibody, .by = "Complex", .before = 1) |>
    dplyr::mutate(Complex = paste0(.data$Antibody, "_", .data$Complex)) |>
    dplyr::arrange(.data$Antibody)

  # Make a list with antibodies, receptors, and bound complexes for the output
  molecules <- list(
    antibodies = levels(antibody),
    receptors = paste0("IL", receptor),
    complexes = unique(complexes$Complex)
  )
  return(molecules)
}

# Conversion from nM to # molecules/cell using the default volume and cell
# number from the flow cytometry experiments; volume = mL
nM2recep <- function(x, volume = 0.2, cells = 1e5) {
  # molecules / mol = 6.022e23
  # mol / nmol = 1/1e9
  # L / mL = 1/1e3
  # Volume = mL
  x * 6.022e23 / 1e12 * volume / cells
}

# Convert the antibody and receptor concentration columns to log concentrations
log_conc <- function(df) {
  df |> dplyr::mutate(dplyr::across(
    tidyselect::starts_with(c("Ab.", "Recep.")),
    log10
  ))
}

# Experimental Data -------------------------------------------------------

# Helper function to filter experimental data to correct normalization basis,
# convert concentrations to log scale, and calculate the mean and standard
# error for each group of replicates
convert_exp <- function(exp, ab = NULL, cell = NULL, basis = "BS1") {

  # If antibodies or cells are not given, use all values present in the data
  if (is.null(ab)) {
    ab <- levels(exp$Antibody)
  }

  if (is.null(cell)) {
    cell <- levels(exp$Cell)
  }

  # Fill in empty values for any antibodies, cell lines, concentrations,
  # replicates, or normalization bases that are missing from the data
  exp <- exp |>
    dplyr::mutate(Rep = forcats::fct_drop(.data$Rep)) |>
    tidyr::complete(
      .data$Antibody,
      .data$Cell,
      .data$Conc,
      .data$Rep,
      .data$Basis
    )

  # Filter the data to just the antibodies, cell lines, and normalization basis
  # used in the model simulations
  exp <- exp |>
    dplyr::filter(
      .data$Antibody %in% ab,
      .data$Cell %in% cell,
      .data$Basis == basis
    )

  # Convert the concentrations to log scale and calculate the mean and standard
  # error for each group of replicates
  exp <- exp |>
    dplyr::mutate(Conc = log10(.data$Conc)) |>
    dplyr::group_by(.data$Conc, .data$Antibody, .data$Cell) |>
    dplyr::summarize(
      N = max(as.numeric(.data$Rep)),
      Mean = mean(.data$Value),
      SE = sd(.data$Value)/sqrt(N),
      .groups = "drop"
    ) |>
    proper_names()

   return(exp)
}

# Fractional Occupancy ----------------------------------------------------

# Wrapper function finding all columns that match the regular expression pattern
# and then summing all values in each row (with NA = 0)
# Used to calculate the total concentration of each different molecule type
# (free, binary, ternary, total) based on matching the column names
add_across <- function(pattern) {
  rowSums(dplyr::pick(tidyselect::matches(pattern)), na.rm = TRUE)
}

# Calculate the free concentration, bound concentration, and fractional
# occupancy of a given antibody or receptor in the system across all complexes
# it can participate in, including binary, ternary, and total complexes
individual_occupancy <- function(output, basis = "recep") {

  # Setup -------------------------------------------------------------------
  # Set the molecules, search patterns, and columns for calculating occupancy

  # Determine which molecules are present in the simulation output
  species <- levels(output$Species)
  species <- species[species != "None"]
  freeAb <- species[!stringr::str_detect(species, "[0-9]R$")]
  freeR <- species[stringr::str_detect(species, "^IL[0-9]R$")]

  # Check if the given basis exists as one of the individual species in the data
  # or if it is one of the special cases; otherwise, throw an error
  other <- c("6R", "8R", "IL6R", "IL8R",
             "R", "Recep", "Receptor", "Ab", "Antibody")
  all_options <- c(freeAb, freeR, other)

  if (!any(stringr::str_equal(all_options, basis, ignore_case = TRUE))) {
    stop(stringr::str_wrap(
      "`basis` argument must be one of the species in the simulation or one of
      the special cases: `6R`, `8R`, `Recep`, `Ab`.",
      exdent = 4
    ))
  }

  # Set the search string and concentration column to use as the denominator
  # based on which antibody or receptor is being used, and specify which type of
  # molecule it is for calculating the free concentration
  # `(?i)` in the regex pattern makes the search case-insensitive
  if (stringr::str_detect(basis, "(?i)^(IL)?(6R|8R)$")) {
    basis <- stringr::str_replace(basis, "^IL", "")
    pattern <- paste0("(?i)", basis)
    column <- paste0("Recep.IL", stringr::str_to_upper(basis))
    molecule <- "recep"
  } else if (stringr::str_detect(basis, "(?i)^(r|recep|receptor)$")) {
    pattern <- "[0-9]R"
    column <- "Recep.Total"
    molecule <- "recep"
  } else if (stringr::str_detect(basis, "(?i)^(ab|antibody)$")) {
    pattern <- paste0("(", paste0(freeAb, collapse = "|"), ")")
    column <- "Ab.Total.Cell"
    molecule <- "ab"
  } else {
    pattern <- species[stringr::str_equal(species, basis, ignore_case = TRUE)]
    column <- paste0("Ab.", pattern, ".Cell")
    molecule <- "ab"
  }

  # Define the simulation outputs that can be used for calculating occupancy
  # Defined in order of importance so redundant outputs can be filtered out of
  # the data set prior to performing the calculations
  types <- c("present", "receptors", "all", "inputs", "binary", "ternary")
  calculations <- c("conc", "conct", "end", "peak")

  # Calculation -------------------------------------------------------------
  # Calculate total concentration and fractional occupancy

  # Determine the how much of the given antibody or receptor in the system is in
  # a particular form, including individual species and binary, ternary, and
  # total complexes, and calculate the fractional occupancy of each form
  output <- output |>
    # Filter the input data frame to only simulations where all molecules were
    # included in the output for only relevant outputs (any of the bound
    # concentrations), and filter to only species where the given molecule is
    # present
    dplyr::filter(
      .data$Interest %in% types,
      .data$Calc %in% calculations,
      .data$Add == FALSE,
      stringr::str_detect(.data$Species, pattern)
    ) |>
    # Convert the molecules of interest column into a factored variable with all
    # of the possible levels, with the factor levels sorted by order of
    # importance of that interest type (e.g., "present" is the most important
    # because it contains every molecule in the system)
    dplyr::mutate(Interest = factor(.data$Interest, levels = types)) |>
    # Group the data by Simulation ID, Calculation, Time, and Species to prepare
    # to filter out redundant data
    dplyr::group_by(.data$ID, .data$Calc, .data$Time, .data$Species)

  # Filter the data to only one row for each group to remove redundant data ID,
  # Calculation, Time, and Species together fully define a unique output, and
  # outputs with the same ID, Calc, Time, and Species but different molecules of
  # interest will have same Bound and Frac values and aren't needed

  # This step will keep rows from different molecules of interest if necessary
  # to fully define the system; e.g., if "binary" and "ternary" are both
  # included in the outputs, "Toci_6R" from "binary" and "Toci_6R_6R" from
  # ternary will both be kept in the data

  # Only filtering the data if it actually will remove rows - this step is very
  # slow for data frames with a large number of groups, so will avoid it when
  # possible; if the number of groups is equal to the number of rows, there is
  # no redundant data so nothing will be filtered
  if (dplyr::n_groups(output) != nrow(output)) {
    output <- output |>
      dplyr::arrange(.data$Interest, .by_group = TRUE) |>
      dplyr::slice_head(n = 1)
  }

  # Ungroup the data after the filter step and start the occupancy calculation
  output |>
    dplyr::ungroup() |>
    # Determine how many of the input antibody or receptor are present in each
    # species (e.g., TociB_6R_6R has two IL6R, and BS1_6R_8R has one IL6R), and
    # scale the complex concentrations by the number of that molecule in the
    # complex; also determining how many total receptors are present for
    # distinguishing binary and ternary complexes
    dplyr::mutate(
      nmatch = stringr::str_count(.data$Species, pattern),
      Bound = Value * nmatch,
      nrecep = stringr::str_count(.data$Species, "[0-9]R")
    ) |>
    dplyr::select(!c("OutputID", "Add", "Interest", "nmatch", "Value")) |>
    # Pivot the data so each species and complex is in a separate column to
    # prepare for calculating the total bound molecule, including the number of
    # receptors in the column name to separate free/binary/ternary; setting any
    # non-existent species to NA to filter out later
    tidyr::pivot_wider(
      names_from = c("Species", "nrecep"),
      names_sep = ".",
      values_from = "Bound",
      values_fill = NA
    ) |>
    # Convert the initial antibody and free antibody concentrations to #/cell
    # to be consistent with the bound antibody
    dplyr::mutate(
      dplyr::across(tidyselect::starts_with("Ab."), \(Conc) nM2recep(Conc),
                    .names = "{.col}.Cell"),
      dplyr::across(tidyselect::ends_with(".0"), \(Conc) nM2recep(Conc))
    ) |>
    # Add the concentrations of each complex type to determine the total free
    # and bound molecule in the different forms; using regular expressions to
    # select all of the columns that correspond to a specific complex type, then
    # using `rowSums()` to add all of the values in each row of those columns
    # (with all NA values set to 0)
    # For the free molecule, using the free antibody columns or the free
    # receptor columns depending on which basis was given, and setting the value
    # to NA if the output did not include the free molecules
    dplyr::mutate(
      Free = dplyr::case_when(
        molecule == "ab" && length(freeAb) > 0 ~ add_across("\\.0$"),
        molecule == "recep" && length(freeR) > 0 ~ add_across("IL[0-9]R\\.1$"),
        .default = NA
      ),
      Binary = add_across("_[0-9]R\\.1$"),
      Ternary = add_across("_[0-9]R\\.2$"),
      Total = add_across("_[0-9]R\\.[0-9]$"),
    ) |>
    # Pivot the data back to long form with separate rows for each complex type
    tidyr::pivot_longer(
      c(tidyselect::matches("\\.[0-9]$"),
        "Free", "Binary", "Ternary", "Total"),
      names_to = "Species",
      values_to = "Bound",
      names_transform = \(x) factor(stringr::str_replace(x, "\\.[0-9]$", ""))
    ) |>
    # Calculate the fractional occupancy of the input molecule type, using the
    # initial concentration of the molecule specified by "basis"
    dplyr::mutate(
      Frac = .data$Bound / .data[[column]]
    ) |>
    # Drop the antibody concentration columns in #/cell because they were only
    # needed for the fractional calculation
    dplyr::select(!tidyselect::ends_with(".Cell"))
}

# Calculate concentration and fractional occupancy for all antibody and receptor
# types and combine all of the occupancy tables for the different molecule types
# into a single table
# `complete` argument dictates whether or not missing species and complexes
# should be included in the output (e.g., Toci_8R_8R or ternary complexes when
# BS1 is used in 6R+ cells)
# `include` argument dictates which molecule types should be used as the bases
# for the occupancy calculation; "all" uses all individual and summed species,
# "totals" uses only the summed antibodies and receptors, and "species" uses
# only the individual species as the bases
occupancy <- function(output, complete = FALSE, include = "all") {

  # Setup -------------------------------------------------------------------
  # Set the molecules for the occupancy bases, and determine which calculations
  # can be performed based on the data available

  # Confirm that the output data table was passed in instead of the list of data
  # frames
  if (!all(c("ID", "Species", "Value") %in% names(output))) {
    stop(stringr::str_wrap(
      "The input to `occupancy()` must be the simulation output data frame,
      not the complete list of data frames.",
      exdent = 4
    ))
  }

  # Define the possible antibodies, receptors, and complexes in the data set
  # species in the data set
  species <- levels(output$Species)
  molecules <- create_molecules(species)

  # Determine which of the possible species are present in the data and set the
  # order for the factored Species column
  mol_vec <- unlist(molecules, use.names = FALSE)
  species <- mol_vec[mol_vec %in% species]
  antibodies <- molecules$antibodies[molecules$antibodies %in% species]
  receptors <- molecules$receptors
  complexes <- c(species, "Free", "Binary", "Ternary", "Total")

  # Define the species to use as bases for the fractional occupancy based on the
  # value of the `include` argument
  if (include == "all") {
    # Using all free antibodies and receptors present in the data and the total
    # antibody and receptor concentrations
    basis <- c(antibodies, "Antibody", receptors, "Receptor")
  } else if (include == "totals") {
    # Using only the total antibody and receptor concentrations
    basis <- c("Antibody", "Receptor")
  } else if (include == "species") {
    # Using all free antibodies and receptors present in the data
    basis <- c(antibodies, receptors)
  } else {
    # Throw an error if none of the specified options were given
    stop(stringr::str_wrap(
      "The value for `include` must be one of: \"all\", \"totals\",
      \"species\".",
      exdent = 4
    ))
  }

  # Determine which output calculations were performed by scanning the relevant
  # columns of the output data
  calculations <- output |>
    dplyr::select(c("Interest", "Add")) |>
    dplyr::distinct() |>
    dplyr::filter(.data$Add == FALSE)

  # Remove bases that cannot be calculated based on the output calculations that
  # were performed; e.g., the "Antibody" basis requires that the free antibody
  # concentrations were included in the output
  if (!any(c("present", "inputs", "all") %in% calculations[["Interest"]])) {
    # If none of the outputs included free antibody, cannot use it as a basis
    basis <- basis[!basis == "Antibody"]
  }

  # Specify the options for the calculation progress bar
  progress_bar <- list(name = "Occupancy", type = "tasks", show_after = 2)

  # Calculation -------------------------------------------------------------
  # Calculate total concentration and fractional occupancy

  # Calculate the occupancy of each species from the basis list separately, then
  # bind the occupancy tables together into a single large table and re-factor
  # the "Species" and "Occupied" columns into the right order
  occupied <-
    purrr::map(
      basis,
      \(x) individual_occupancy(output, basis = x),
      .progress = progress_bar
    ) |>
    rlang::set_names(basis) |>
    purrr::list_rbind(names_to = "Occupied") |>
    dplyr::mutate(
      Species = factor(.data$Species, levels = complexes),
      Occupied = factor(.data$Occupied, levels = basis)
    )

  # If the missing species and complexes should be included in the data set
  if (complete) {

    # Set which columns in the data set have data that needs to be copied into
    # the added rows
    # These columns will be used with the `nesting()` function so only
    # combinations of these columns that occur in the data will be used in the
    # complete data set; for the "Species" and "Occupied" columns, all unique
    # combinations will be expanded out, even if they are not currently present
    # in the data
    cols <- names(occupied)
    cols <- cols[!cols %in% c("Species", "Occupied", "Bound", "Frac")]

    # Replace missing values with 0 and fill in all missing species and
    # complexes in the output data
    occupied <- occupied |>
      # Replace all NA values in all columns with 0
      dplyr::mutate(dplyr::across(
        tidyselect::everything(),
        \(x) tidyr::replace_na(x, 0)
      )) |>
      # Complete all missing combinations of "Species" and "Occupied" and fill
      # in 0s where there is missing data
      tidyr::complete(
        tidyr::nesting(!!!rlang::syms(cols)),
        .data$Species,
        .data$Occupied,
        fill = list(Bound = 0, Frac = 0)
      )
  } else {
    # If the missing data is not needed in the output, drop all rows with NA
    # values
    # NA was set as the "Bound" and "Frac" values for species that were not
    # present in the data set
    occupied <- occupied |>
      tidyr::drop_na(c("Bound", "Frac"))
  }

  # Reorder the columns and sort the data into the desired output order
  occupied |>
    dplyr::relocate(c("Occupied", "Species"), .before = "Bound") |>
    dplyr::arrange(ID, Calc, Time, Occupied, Species)
}
