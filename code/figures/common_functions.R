
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

# Define the labels to use for the different complex types
species_labels <- list(
  "Binary Complexes" = "Binary",
  "Ternary Complexes" = "Ternary",
  "Total Bound" = "Total"
)

# Create Heatmap Plots ----------------------------------------------------

# Generic heat map function that uses consistent formatting and styling
plot_heatmap <- function(df, x, y, xlab, ylab, palette = "pals::ocean.matter",
                         nbreaks = 7) {
  df |>
    proper_names() |>
    mutate(Species = fct_recode(Species, !!!species_labels)) |>
    ggplot(aes(x = {{ x }}, y = {{ y }}, fill = Frac)) +
    geom_raster() +
    scale_x_continuous(
      name = xlab,
      expand = expand0(),
      labels = label_math(expr = 10^.x),
      breaks = breaks_extended(n = nbreaks)
    ) +
    scale_y_continuous(
      name = ylab,
      expand = expand0(),
      labels = label_math(expr = 10^.x)
    ) +
    scale_fill_paletteer_c(
      name = "Occupancy",
      palette = palette,
      guide = heatmap_legend(plot_scale = 0.7),
      limits = c(-0.005, 1.005)
    ) +
    facet_grid2(
      . ~ Species,
      strip = strip_vanilla(clip = "off")
    ) +
    theme_cr() +
    theme(
      panel.spacing = unit(20, "pt"),
      legend.title = element_text(margin = margin(b = 6))
    )
}

# Create Line Plots -------------------------------------------------------

# Helper function for just the plotting code for the line plots of the
# fractional occupancy of the antibodies and receptors
plot_lines <- function(df, x, color, xlab, colorlab, ylab,
                           palette = "Sunset") {
  df |>
    proper_names() |>
    ggplot(aes(x = {{ x }}, y = Frac, color = {{ color }})) +
    geom_line(alpha = 0.8) +
    scale_x_continuous(
      name = xlab,
      labels = label_math(expr = 10^.x)
    ) +
    scale_y_continuous(name = ylab, limits = c(0, 1)) +
    scale_color_manual(
      name = colorlab,
      values = palette_gen(6, palette = palette),
      labels = label_math(expr = 10^.x)
    ) +
    facet_grid2(
      . ~ Species,
      strip = strip_vanilla(clip = "off")
    ) +
    theme_cr()
}

# Combine Line Plots ------------------------------------------------------

# Create and combine line plots for the antibodies and the receptors
plot_recep_ab <- function(df, x, color, xlab, colorlab, palette = "Sunset") {

  ab_panel <- df |>
    filter(Occupied == "Antibody", Species == "Free") |>
    mutate(Species = fct_recode(Species, "Free Antibody" = "Free")) |>
    plot_lines(
      x = {{ x }},
      color = {{ color }},
      xlab = xlab,
      colorlab = colorlab,
      ylab = "Fraction of Total Antibody",
      palette = palette
    ) +
    guides(color = "none")

  recep_panel <- df |>
    filter(
      Occupied == "Receptor",
      Species %in% c("Free", "Binary", "Ternary", "Total")
    ) |>
    mutate(Species = fct_recode(Species, !!!species_labels)) |>
    mutate(Species = fct_recode(Species, "Free Receptor" = "Free")) |>
    plot_lines(
      x = {{ x }},
      color = {{ color }},
      xlab = xlab,
      colorlab = colorlab,
      ylab = "Fraction of Total Receptor",
      palette = palette
    )

  nplots <- nlevels(fct_drop(df$Species))

  ab_panel + recep_panel +
    plot_layout(widths = c(1, nplots))
}
