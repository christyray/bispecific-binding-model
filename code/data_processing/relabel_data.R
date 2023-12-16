
# NAME DEFINITIONS --------------------------------------------------------

# Define the short variable names to be available to the functions - used for
# collapsing the factor levels into standardized shorter names
short_name_definitions <- function() {

  # Initialize output object
  defs <- list()

  # Factor Levels -----------------------------------------------------------
  # Set up the variables for the factor levels and the collapsed factor names

  # Columns = Possible names of the columns containing factors to be collapsed
  # and renamed
  # Names = all possible names that could be used for the levels of a given
  # model variable
  # Types = list that assigns each of the possible names to a specific standard
  # name that is used across the R analysis and plotting code

  defs$ab <- list(
    columns = c("Antibody", "antibody", "Ab", "ab"),
    names = c(
      "TociM", "TociB", "Toci", "toci",
      "H2M", "H2B", "H2", "h2",
      "BS1", "bs1", "BS2M", "BS2B", "BS2", "bs2",
      "TociM_H2M", "TociB_H2M", "TociM_H2B", "TociB_H2B", "Toci_H2", "toci_h2"
    ),
    types = list(
      "Toci" = c("TociM", "TociB", "Toci", "toci"),
      "H2" = c("H2M", "H2B", "H2", "h2"),
      "BS1" = c("BS1", "bs1"),
      "Toci_H2" = c("TociM_H2M", "TociB_H2M", "TociM_H2B",
                    "TociB_H2B", "Toci_H2", "toci_h2"),
      "BS2" = c("BS2M", "BS2B", "BS2", "bs2")
    )
  )

  defs$recep <- list(
    columns = c("Receptor", "receptor", "Recep", "recep"),
    names = c("IL6R", "IL8R", "6R", "8R", "il6r", "il8r", "6r", "8r",
              "IL-6R", "IL-8R", "il-6r", "il-8r"),
    types = list(
      "IL6R" = c("IL6R", "6R", "il6r", "6r", "IL-6R", "il-6r"),
      "IL8R" = c("IL8R", "8R", "il8r", "8r", "IL-8R", "il-8r")
    )
  )

  defs$cell <- list(
    columns = c("Cell", "cell"),
    names = c(
      "6R", "8R",
      "6R-6R", "6R-8R", "8R-8R", "6R6R", "6R8R", "8R8R",
      "6R-6R-8R", "6R-8R-8R", "6R-6R-8R-8R", "6R6R8R", "6R8R8R", "6R6R8R8R",
      "6R_8R", "6R-6R_8R-8R", "6R-6R-8R_6R-8R-8R", "6R-6R_6R-6R", "8R-8R_8R-8R"
    ),
    types = list(
      "6R" = c("6R", "6R-6R", "6R-6R_6R-6R", "6R6R"),
      "8R" = c("8R", "8R-8R", "8R-8R_8R-8R", "8R8R"),
      "6R8R" = c(
        "6R-8R", "6R-6R-8R", "6R-8R-8R", "6R-6R-8R-8R",
        "6R_8R", "6R-6R_8R-8R", "6R-6R-8R_6R-8R-8R",
        "6R8R", "6R6R8R", "6R8R8R", "6R6R8R8R"
      )
    )
  )

  defs$k <- list(
    columns = c("Parameter", "parameter", "Param", "param"),
    names = c(
      "kon_Ab_6R$", "kon6R", "kon_Ab_8R$", "kon8R",
      "kon_Ab.*_[68]R_6R$", "kon6Rprime", "kon_Ab.*_[68]R_8R$", "kon8Rprime",
      "koff_Ab.*_6R$", "koff6R", "koff_Ab.*_8R$", "koff8R",
      "koff_Ab", "koff",
      "initial_6R", "initial_8R", "initial_6R_prime", "initial_8R_prime",
      "initial_off", "initial_off6R", "initial_off8R"
    ),
    types = list(
      "kon6R" = c("kon_Ab_6R$", "kon6R", "initial_6R"),
      "kon8R" = c("kon_Ab_8R$", "kon8R", "initial_8R"),
      "kon6Rprime" = c("kon_Ab.*_[68]R_6R$", "kon6Rprime", "initial_6R_prime"),
      "kon8Rprime" = c("kon_Ab.*_[68]R_8R$", "kon8Rprime", "initial_8R_prime"),
      "koff" = c("koff_Ab", "koff", "initial_off"),
      "koff6R" = c("koff_Ab.*_6R$", "koff6R", "initial_off6R"),
      "koff8R" = c("koff_Ab.*_8R$", "koff8R", "initial_off8R")
    )
  )

  defs$basis <- list(
    columns = c("Basis", "basis"),
    names = c("bs1", "BS1", "ab", "Ab"),
    types = list("BS1" = c("bs1", "BS1"), "Ab" = c("ab", "Ab"))
  )

  defs$calc <- list(
    columns = c("NormCalc", "normcalc"),
    names = c("data", "Data", "max", "Max"),
    types = list("Data" = c("data", "Data"), "Max" = c("max", "Max"))
  )

  defs$interest <- list(
    columns = c("Interest", "interest"),
    names = c("receptors", "Receptors", "binary", "Binary",
              "ternary", "Ternary", "complex", "Complex",
              "complexes", "Complexes", "total", "Total",
              "present", "Present"),
    types = list(
      "receptors" = c("receptors", "Receptors"),
      "binary" = c("binary", "Binary"),
      "ternary" = c("ternary", "Ternary"),
      "complex" = c("complex", "Complex",
                    "complexes", "Complexes", "total", "Total"),
      "present" = c("present", "Present")
    )
  )

  defs$calculations <- list(
    columns = c("Calc", "calc", "Calculation", "calculation"),
    names = c("peak", "Peak", "auc", "AUC", "conct", "ConcT", "conc", "Conc",
              "end", "End"),
    types = list(
      "peak" = c("peak", "Peak"),
      "auc" = c("auc", "AUC"),
      "conct" = c("conct", "ConcT"),
      "conc" = c("conc", "Conc"),
      "end" = c("end", "End")
    )
  )

  return(defs)
}

# Define the formatted variable names for plotting - used to replace the
# standardized short names after data analysis
proper_name_definitions <- function() {

  # Initialize output object
  defs <- list()

  # Factor Levels -----------------------------------------------------------
  # Set up the variables for the factor levels and the collapsed factor names

  # Columns = Possible names of the columns containing factors to be collapsed
  # and renamed
  # Names = the standardized names used for the levels of a given model variable
  # Types = list that assigns each of the standard variable names to its
  # plotmath formatted name for labeling in figures

  defs$ab <- list(
    columns = c("Antibody", "antibody", "Ab", "ab"),
    names = c("Toci", "H2", "BS1", "Toci_H2", "BS2"),
    types = list(
      "Tocilizumab" = "Toci",
      "10H2" = "H2",
      "BS1" = "BS1",
      "Tocilizumab + 10H2" = "Toci_H2",
      "BS2" = "BS2"
    )
  )

  defs$cell <- list(
    columns = c("Cell", "cell"),
    names = c("6R", "8R", "6R8R"),
    types = list(
      '"IL-6R"^{"+"} ~ Cells' = "6R",
      '"IL-8R"^{"+"} ~ Cells' = "8R",
      '"IL-6R"^{"+"} ~ "IL-8R"^{"+"} ~ Cells' = "6R8R"
    )
  )

  defs$k <- list(
    columns = c("Parameter", "parameter", "Param", "param"),
    names = c("kon6R", "kon8R", "kon6Rprime", "kon8Rprime",
              "koff", "koff6R", "koff8R"),
    types = list(
      'k["on,6R"]' = "kon6R",
      'k["on,8R"]' = "kon8R",
      'k["on,6R*"]' = "kon6Rprime",
      'k["on,8R*"]' = "kon8Rprime",
      'k["off"]' = "koff",
      'k["off,6R"]' = "koff6R",
      'k["off,8R"]' = "koff8R"
    )
  )

  defs$interest <- list(
    columns = c("Interest", "interest"),
    names = c("receptors", "binary", "ternary", "complex", "present"),
    types = list(
      "Bound Antibodies" = "receptors",
      "Binary Complexes" = "binary",
      "Ternary Complexes" = "ternary",
      "Total Complexes" = "complex",
      "All Species" = "present"
    )
  )

  defs$calculations <- list(
    columns = c("Calc", "calc", "Calculation", "calculation"),
    names = c("peak", "auc", "conct", "conc", "end"),
    types = list(
      "Peak Concentration" = "peak",
      "Area Under Concentration Curve" = "auc",
      "Concentration at Time Point" = "conct",
      "Concentration" = "conc",
      "Concentration at End Time" = "end"
    )
  )

  return(defs)
}

# HELPER FUNCTIONS --------------------------------------------------------

# Convert levels of a variable into the short names for that variable type
rename_factors <- function(data, var, name) {

  # Setup -------------------------------------------------------------------
  # Import the lists of factor levels and matching collapsed factor names
  # depending on which was requested
  if (name == "short") {
    defs <- short_name_definitions()
  } else if (name == "proper") {
    defs <- proper_name_definitions()
  }

  # Collapse Factor Levels --------------------------------------------------
  # Determine which columns in the data frame correspond to the model variable
  # being renamed, convert the column into a factor variable with all possible
  # variable labels used as the factor levels, collapse the factor levels into
  # the standardized variable levels

  # Columns = possible names of the columns containing factors to be collapsed
  # and renamed
  # Names = all possible names that could be used for the levels of a given
  # model variable
  # Types = list that assigns each of the possible names to a specific standard
  # name that is used across the R analysis and plotting code

  data <- dplyr::mutate(data, dplyr::across(
    tidyselect::any_of(defs[[var]]$columns),
    ~ forcats::fct_collapse(
      factor(.x, levels = defs[[var]]$names),
      !!!defs[[var]]$types
    )
  ))
  # !!! turns a list of inputs into distinct arguments since fct_recode()
  # requires separate named strings instead of a named vector; similar concept
  # to do.call()

  return(data)
}

# Find the proper name corresponding to the input short name(s) for filtering
# data frames that have already been converted to "proper names"
match_proper <- function(short_name, type = NULL) {

  # Setup -------------------------------------------------------------------
  # Import the lists of factor levels and matching collapsed factor names
  defs <- proper_name_definitions()

  # If type is not provided, determine it from the short name
  if (is.null(type)) {

    # Find the index of the variable list that contains the short name
    # `detect_index()` finds the first list index where the function is true
    idx <- purrr::detect_index(defs, ~ short_name %in% .x$names)

    # Get the first option of the possible column names for the variable that
    # contains the matching short name to use for the later search
    type <- defs[[idx]]$columns[1]
  }

  # Filter the proper names to just the names that correspond to the requested
  # variable type

  # `keep` will keep only the list elements where the given function is true
  # In this case, only keeps elements where the value given for `type` is
  # present in the list of column names for that variable

  # `flatten` removes extraneous list levels so only the single variable that
  # matched `type` will remain
  defs <- purrr::flatten(purrr::keep(defs, ~ type %in% .x$columns))

  # Match Names -------------------------------------------------------------
  # Return the "proper name" corresponding to the short name input(s)

  proper_name <- names(defs$types[defs$types %in% short_name])
  return(proper_name)
}

# RENAME FACTORS ----------------------------------------------------------

# Relabel the factors in order and with the short names for convenience
short_names <- function(data) {

  # Renaming ----------------------------------------------------------------
  # Rename the factor levels and collapse the factors into the correct levels

  # Determine all of the variables that have standard names available
  vars <- names(short_name_definitions())

  # For each available variable, rename the factors in the matching columns to
  # the standard names
  for (var in vars) {
    data <- rename_factors(data, var, "short")
  }

  return(data)
}

# Relabel the factors in order and with the correct full names for plotting
proper_names <- function(data) {

  # Renaming ----------------------------------------------------------------
  # Rename the factor levels into the proper names for plotting

  # Determine all of the variables that have proper names available
  vars <- names(proper_name_definitions())

  # First, pass the data through `short_names()` to standardize the names
  data <- short_names(data)

  # For each available variable, rename the factors in the matching columns to
  # the proper names
  for (var in vars) {
    data <- rename_factors(data, var, "proper")
  }

  return(data)
}

# Relabel the species column when necessary
# Not including in the other renaming steps because it is not frequently needed
# and it only works when the data does not include the summed complexes
rename_species <- function(data) {

  # Generate the names of all possible antibodies, receptors, and
  # antibody-receptor complexes in the system, then combine the molecules into a
  # single vector and append the summed species to the end
  orig <- create_molecules(levels(data$Species))
  orig <- c(orig$receptors, orig$antibodies, orig$complexes)
  orig <- c(orig, "Free", "Binary", "Ternary", "Total")

  # Create the formatted labels for each species
  # Removing the M or B from the antibody names, using the full names for Toci
  # and H2, adding dashes after "IL", and replacing underscores with dashes to
  # make the formatted names
  new <- orig
  new <- stringr::str_remove(new, "(M|B)$")
  new <- stringr::str_replace(new, "(M|B)_", "_")
  new <- stringr::str_replace(new, "Toci", "Tocilizumab")
  new <- stringr::str_replace(new, "H2", "10H2")
  new <- stringr::str_replace_all(new, "(IL)(6|8)", "\\1-\\2")
  new <- stringr::str_replace_all(new, "_", "-")

  # Create a named vector to use in `fct_collapse()` where the names are the new
  # factor levels and the values are the original factors
  collapse <- orig
  names(collapse) <- new

  # Re-level the Species column using the new factor levels and labels
  data <- dplyr::mutate(
    data,
    Species = forcats::fct_collapse(factor(Species, levels = orig), !!!collapse)
  )

  return(data)
}
