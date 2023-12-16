here::i_am("code/data_processing/compare_data.R")

# SOURCE FUNCTIONS --------------------------------------------------------

# Import helper functions to import and organize the simulation data
source(here::here("code/data_processing/initialization.R"))

# HELPER FUNCTIONS --------------------------------------------------------

# Determine the symmetric difference between the two data sets in the values in
# the filtering column; returns all of the values that are not present in both
# data sets so the data sets can be filtered to compare only the simulations
# that were done in both sets
symdiff_col <- function(A, B, col) {

  # If the column is present in both data frames, return the values that need
  # to be filtered out of both data sets
  if (col %in% names(A) && col %in% names(B)) {
    values <- dplyr::symdiff(A[[col]], B[[col]])
  } else {
    values <- NULL
  }
  return(values)
}

# Filter out specified values from a column of a data frame and return ID values
# to be removed from the other data frames from the simulation; used as part of
# the filter_inputs() function to remove extra inputs from all data frames in
# simulation
filter_id <- function(data, col, values) {

  # If the column is present in the data frame, filter to the given `values`
  # from the column and return the ID values that need to be removed; otherwise,
  # return an empty table
  if (col %in% names(data)) {
    id <- data |>
      dplyr::filter(.data[[col]] %in% values) |>
      dplyr::select("ID")
  } else {
    id <- tibble::tibble()
  }
  return(id)
}

# Filter out specified values from a column of a data frame, then remove the
# filtered out IDs from all data frames in the input list and drop all unused
# factor levels from all data frames in the input list; useful when one of the
# data sets had extra simulations that need to be removed entirely from the
# comparison
filter_inputs <- function(data, col, values) {

  # Determine which ID values correspond to the values to be removed
  # Search all tables in the list for the values to be removed, and return the
  # corresponding ID values, then combine the ID values into a single table and
  # keep only the unique values
  id <- lapply(data, filter_id, col = col, values = values)
  id <- do.call(rbind, id) |> dplyr::distinct()

  # Remove those ID values from all tables in the list
  data <- lapply(data, dplyr::filter, !.data$ID %in% id$ID)

  # Drop all unused factor levels from all tables in the list
  data <- lapply(data, function(df) {
    dplyr::mutate(
      df,
      dplyr::across(tidyselect::where(is.factor), ~ forcats::fct_drop(.x))
    )
  })
  return(data)
}

# Replace factor levels with a sequence of numbers, used along with `mutate` to
# re-number that ID values after filtering out inputs
recode_id_col <- function(col) {
  forcats::fct_relabel(col, ~ as.character(1:length(levels(col))))
}

# Re-number all ID columns in a data frame, used with `lapply` to re-number the
# ID columns in all tables in a set of simulation data
recode_id <- function(df) {
  dplyr::mutate(
    df,
    dplyr::across(tidyselect::ends_with("ID"), ~ recode_id_col(.x))
  )
}

# Round numeric data to remove small precision differences between values
round_column <- function(data, digits) {

  # Round all numeric columns to the given significant figures before comparing
  data <- data |>
    dplyr::mutate(dplyr::across(
      tidyselect::where(is.numeric),
      ~ signif(.x, digits = digits)
    ))
}

# Check if two data frames are equal and display a progress message after the
# check with information
all_equal_msg <- function(A, B, var, truncate = TRUE) {
  eq <- all.equal(A, B)

  # Displays a success message if both data frames are equal or information
  # about where they are unequal
  if (isTRUE(eq)) {
    print(glue::glue("`{var}` complete: Both data sets are equal."))
  } else {
    # Truncate the output of unequal rows to a reasonable length if requested
    if (isTRUE(truncate)) {
      eq <- substr(eq, start = 1, stop = 500)
      eq <- glue::glue('{eq}... (truncated at 500 characters)')
    }
    print(glue::glue('`{var}` complete: \n{eq}'))
  }

  # Display max error between calculated values if present
  if ("Value" %in% names(A)) {
    re <- max(abs((A$Value - B$Value) / A$Value), na.rm = TRUE)
    re <- format(re, digits = 4)
    print(glue::glue('Maximum relative error in `{var}` = {re}\n\n'))
  }
}

# COMPARE DATA ------------------------------------------------------------

# Import and clean data from different simulations, filter out extra data if
# necessary, and verify that the two simulations are equal
compare_data <- function(old, new, type = NULL, col = NULL, values = NULL,
                         fixID = TRUE, digits = NULL, truncate = TRUE) {

  # Import and Clean Data ---------------------------------------------------
  # Import the data from both simulations and organize the data into standard
  # tables for comparison
  old_data <- import_data(old, type = type)
  new_data <- import_data(new, type = type)

  # Filter Extra Data -------------------------------------------------------
  # Remove rows with extra data that is not present in both simulations

  # If columns to filter on were provided, filter out the IDs that correspond to
  # the given `values` in each column in both simulations; useful if additional
  # inputs or outputs were present in one simulation
  if (!is.null(col)) {
    for (i in seq_along(col)) {

      # If no values were provided for the column, filter all of the values in
      # that column that are not present in both data sets
      if (is.null(values[[i]])) {
        # Search all tables in the list for the filtering column and return the
        # values that are not present in that column in both data sets, then
        # combine the values into a single table and keep only unique entries
        values[[i]] <- lapply(
          names(old_data),
          function(x) { symdiff_col(old_data[[x]], new_data[[x]], col[[i]]) }
        )
        values[[i]] <- do.call(c, values[[i]]) |> unique()
      }

      # Filter out the extra values from both data sets
      old_data <- filter_inputs(old_data, col = col[[i]], values = values[[i]])
      new_data <- filter_inputs(new_data, col = col[[i]], values = values[[i]])
    }

    # Re-number the ID columns after filtering out IDs - optional and purely so
    # that the ID columns do not generate extra messages about not being equal
    # when the rest of the columns match
    if (isTRUE(fixID)) {
      old_data <- lapply(old_data, recode_id)
      new_data <- lapply(new_data, recode_id)
    }
  }

  # Round Numeric Data ------------------------------------------------------
  # Round numeric data to remove small precision differences between values

  # If a number of significant figures was provided, round the numeric data to
  # the provided number of significant figures before comparing
  if (!is.null(digits)) {
    old_data <- lapply(old_data, round_column, digits = digits)
    new_data <- lapply(new_data, round_column, digits = digits)
  }

  # Compare Data Tables -----------------------------------------------------
  # Verify equality of the simulations

  # Compare each data table for equality between the simulations and display a
  # progress message
  if (is.data.frame(old_data) && is.data.frame(new_data)) {
    all_equal_msg(old_data, new_data, type, truncate = truncate)
  } else {
    # If the data is comprised of multiple tables in a list, compare each table
    # separately
    # Only compare tables that are present in both lists
    tables_both <- intersect(names(old_data), names(new_data))
    for (var in tables_both) {
      all_equal_msg(old_data[[var]], new_data[[var]], var, truncate = truncate)
    }

    # Print a message for any tables that were not present in both lists
    tables_old <- setdiff(names(old_data), names(new_data))
    tables_new <- setdiff(names(new_data), names(old_data))

    if (length(tables_old) > 0) {
      print(glue::glue("Tables present in old data only: ",
                       "{paste(tables_old, collapse = \", \")}"))
    }
    if (length(tables_new) > 0 ) {
      print(glue::glue("Tables present in new data only: ",
                       "{paste(tables_new, collapse = \", \")}"))
    }
  }

  # Return the data from both simulations for manual checking
  data <- list(old = old_data, new = new_data)
  invisible(data)
}
