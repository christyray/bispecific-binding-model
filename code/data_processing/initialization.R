
# INITIALIZATION ----------------------------------------------------------
# Load the necessary packages and R helper functions together to make them
# available for the data analysis and plotting code

# Packages ----------------------------------------------------------------

# `here` package makes it easier to specify relative paths to code/data files
library(here)

# Functions ---------------------------------------------------------------

# `relabel_data` contains functions to map the antibodies, complexes, cells, and
# other variables in the data to a standardized set of names; `short_names()`
# converts all data to standard, easy-to-reference names; `proper_names()`
# converts the short names to prettier versions for the final figures
source(here("code/data_processing/relabel_data.R"))

# `import_data` contains functions to import and tidy model simulation data;
# `import_files()` imports all data associated with a simulation; the
# `clean_*()` functions tidy the imported data and convert it to long form
source(here("code/data_processing/import_data.R"))

# `calculations` contains functions to perform frequently used calculations on
# the data; `receptor_occupancy()` calculates the fractional occupancy of a
# specific receptor type; `occupancy()` calculates the fractional occupancy of
# each receptor type and combines all of the results into a single table
source(here("code/data_processing/calculations.R"))

# `utilities` contains other helper functions for miscellaneous data tasks;
# `match_near()` works similarly to my `isequaltol()` function in MATLAB and
# determines pairwise equality between two vectors within a relative tolerance;
# `find_near()` works similarly to `ismembertol()` in MATLAB and determines if
# values within one vector are present within another vector within a relative
# tolerance; `%~%` is a binary operator shortcut for `find_near()`
source(here("code/data_processing/utilities.R"))
