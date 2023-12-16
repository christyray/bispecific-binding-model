
# COMPARE VALUES PAIRWISE WITHIN TOLERANCE --------------------------------

# Compare two numbers or vectors (pairwise) with relative tolerance; similar to
# `near()` from `dplyr` except that it uses relative tolerance instead of
# absolute tolerance
# Returns logical vector the same length as x and y indicating where the values
# are equal (within the tolerance)
match_near <- function(x, y, tol = 1e-12) {
  # Using 1e-12 as the default because MATLAB's `ismembertol()` uses 1e-12 as
  # tolerance for double-precision numbers

  # Determine tolerance using the given relative tolerance value
  # `pmax()` determines the parallel maximum - compares matching elements from
  # x and y and returns the maximum at each position
  tol <- pmax(abs(x), abs(y)) * tol

  # Find matches that are within the calculated tolerance
  abs(x - y) <= tol
}

# FIND VALUES WITHIN TOLERANCE --------------------------------------------

# Find all values in x that are anywhere within in y within a relative tolerance
# Returns logical vector the same length as x indicating which values of x are
# present somewhere within y (within the tolerance)
find_near <- function(x, y, tol = 1e-12) {
  # Using 1e-12 as the default because MATLAB's `ismembertol()` uses 1e-12 as
  # tolerance for double-precision numbers

  # Use relative tolerance to determine tolerance for differences in input
  # vectors
  tol <- max(c(abs(x), abs(y))) * tol

  # Pre-allocate output logical vector
  out <- logical(length(x))

  # For each value in the first vector, determine if it is within the tolerance
  # of any values in the second vector
  for (i in 1:length(x)) {
    out[i] <- any(abs(x[i] - y) <= tol)
  }
  return(out)
}

# Infix operator shortcut for `find_near()`
# Note that user-defined infix operators have a higher order of precedence than
# most other operations, including arithmetic operations
`%~%` <- function(x, y, tol = 1e-12) {
  find_near(x, y, tol = tol)
}

# ROUND TO VALUES IN LIST -------------------------------------------------

# Round a set of numbers to the closest values from an arbitrary list; helpful
# to remove rounding errors when simulations were given a specific set of inputs
# Adapted from StackOverflow post: https://stackoverflow.com/questions/12861061/
round_list <- function(values, round_to) {

  # List of numbers to compare to needs to be sorted for `findInterval`
  round_to <- sort(round_to)

  # `findInterval` returns the interval in round_to that each value belongs to;
  # i.e., interval is the index of the closest number smaller than the value and
  # interval + 1 is the index of the closest number larger than the value
  interval <- findInterval(values, round_to)

  # Use the interval to return the closest numbers on either side of the value;
  # need to include -Inf and Inf for numbers that are outside of the range of
  # the rounding list
  low <- c(-Inf, round_to)[interval + 1]
  high <- c(round_to, Inf)[interval + 1]

  # Determine where the values are closer to the higher number
  # Absolute value is not necessary because low < values < high
  idx_high <- high - values < values - low

  # Replace the low values with the high values where the high values are better
  low[idx_high] <- high[idx_high]
  return(low)
}

# CONVERT FACTOR AND NUMERIC ----------------------------------------------

# Convert factored variable to numeric
# Wrapper function for the base R method of conversion to make it clearer in
# model functions
fct2num <- function(x) {
  as.numeric(levels(x))[x]
}

# Convert a numeric vector to a factored variable with a specific set of levels;
# useful for variables like the ratio of IL6R to IL8R where set of values are
# known
# Attempts to convert factor and character inputs to numeric and throws an error
# if conversion introduces NAs
# Uses the `round_list()` function from to first round each value to the closest
# number within the set of possible levels, then maps each value to the matching
# factor level and label
factor_number <- function(x, levels, labels = levels) {

  # Convert factor and character inputs to numeric
  tryCatch({
    if (class(x) == "factor") {
      x <- as.numeric(levels(x))[x]
    } else if (class(x) == "character") {
      x <- as.numeric(x)
    }
  }, warning = function(warn, name, call = rlang::caller_env(n = 4)) {

    cli::cli_abort(c(
      "{.var x} must be coercible to {.cls numeric}.",
      "x" = "Conversion of {.var x} introduced NAs."
    ), call = call)
  })

  factor(
    round_list(x, round_to = levels),
    levels = levels,
    labels = labels
  )
}

# REPLACE ALL NA ----------------------------------------------------------

# Replace all NA values in a data frame with a specific value (default 0)
replace_all_na <- function(df, replace = 0) {
  # Replace all NA values with 0s across all columns of the data frame
  df |>
    dplyr::mutate(dplyr::across(
      tidyselect::everything(),
      \(x) tidyr::replace_na(x, replace)
    ))
}

# CREATE COLOR PALETTE ----------------------------------------------------

# Generate colors for plots from a continuous color scale
# `colour_ramp` creates a function that maps the interval [0,1] to the color
# space of the input color palette, and that function can be used with a
# sequence of numbers to get a discrete color palette of a particular length

# Preview the generated palette with `scales::show_col()`

# `n` = the number of colors to generate
# `palette` = the name of the color palette to select colors from; must be from
# a package included in the `paletteer` package
# `package` = the package the named `palette` is from; must be a package
# included in the `paletteer` package
# `range` = the interval to map the colors to, within [0,1]; a wider range
# provides more variation in the colors and more extreme ends of the color
# spectrum
palette_gen <- function(n, palette = "PurpOr", package = "rcartocolor",
                        range = c(0.1, 1)) {

  # Concatenate the package and palette unless a package was listed with the
  # palette name
  if (!stringr::str_detect(palette, "::")) {
    palette <- paste0(package, "::", palette)
  }

  # Generate the specified number of colors from the palette
  scales::colour_ramp(paletteer::paletteer_d(palette))(
    seq(range[[1]], range[[2]], length.out = n)
  )
}

# TAG PLOTS ---------------------------------------------------------------

# Add annotation tags to `patchwork` plots and add extra space in front of
# Y-axis title

# `plot_scale` = the scale factor used in `theme_cr()`
# `tags` = the type of letters or numbers to use for the annotation tags
tag_plots <- function(plot_scale = 1, tags = "A") {
  font_size <- plot_scale * 12

  return(list(
    patchwork::plot_annotation(tag_levels = tags),
    ggplot2::theme(axis.title.y = ggplot2::element_text(
      margin = ggplot2::margin(l = font_size, r = font_size * 0.75)
    ))
  ))
}
