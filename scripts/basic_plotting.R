# Function to define the custom theme
eb <- element_blank()

theme_custom <- function() {
  theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.major = eb,
      panel.grid.minor = eb,
      strip.background = eb,
      strip.text.x = element_text(size = 16, face = "bold"),
      panel.border = element_rect(colour = "black"),
      axis.title = element_text(size = 15, face = "bold"),
      axis.title.x = element_text(
        size = 23,
        face = "bold",
        colour = "black"
      ),
      axis.title.y = element_text(
        size = 21,
        face = "bold",
        colour = "black"
      ),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 17, face = "bold"),
      legend.box = "horizontal",
      legend.direction = "horizontal",
      axis.text.y = element_text(size = 21, face = "bold"),
      axis.text.x = element_text(size = 21, face = "bold")
    )
}

### custom colors
get_viridis_colors <- function(
    df_to_plot,
    plot_by,
    pallette_type = "D",
    color_direction,
    begin_point) {
  # Define the number of samples
  number_of_samples <- df_to_plot %>%
    distinct({{ plot_by }}) %>%
    nrow()
  # Create a custom color palette based on an existing one (e.g., viridis)
  custom_palette <- viridisLite::viridis(
    number_of_samples,
    option = pallette_type,
    direction = color_direction,
    begin = begin_point
  )

  return(custom_palette)
}

## colours specific to nifh clusters
# nifh_cluster_colours <- c(
#   "1A" = "steelblue",
#   "1J/1K" = "chocolate4",
#   "1O/1P" = "magenta",
#   # "3" = "goldenrod2",
#   "3" = "red",
#   "1G" = "chocolate1",
#   "1B" = "green2",
#   "other" = "lightgrey",
#   "unknown" = "darkgrey"
# )

nifh_cluster_colours <- c(
  "1A" = "steelblue",
  "1J/1K" = "chocolate4",
  "1O/1P" = "magenta",
  # "3" = "goldenrod2",
  "3" = "red",
  "1G" = "chocolate1",
  "1B" = "green2",
  "other" = "lightgrey",
  "unknown" = "darkgrey",
  "4" = "black",
  "2" = "blue"
)

# Function to create a custom dot plot
create_custom_plot_point <- function(
    data1, data2, x1, y1, x2, y2,
    colour, group, size, label, legend_position) {
  ggplot() +
    geom_point(
      data = data1,
      aes(x = {{ x1 }}, y = {{ y1 }}, colour = {{ colour }}, group = {{ group }}, size = {{ size }}),
      na.rm = TRUE,
      shape = 1,
      stroke = 1
    ) +
    geom_text(
      data = data2,
      aes(x = {{ x2 }}, y = {{ y2 }}, label = {{ label }}),
      size = 10,
      fontface = "bold",
      vjust = -1
    ) +
    labs(
      y = expression(bold(paste(
        "PO"[4]^"3-" * " (µmol L"^"-1", ")"
      ))),
      x = expression(bold(paste("SST (˚C)"))),
      title = "",
      shape = "hemisphere",
      fill = "Taxa"
    ) +
    scale_size_continuous(
      breaks = c(0, 0.1, 0.33, 0.66, 1.00),
      range = c(.1, 10),
      limits = c(0, 1)
    ) +
    scale_color_fermenter(
      type = "div",
      palette = "Spectral",
      direction = -1,
      breaks = c(-11, -9, -7)
    ) +
    ylim(-0.05, 1.95) +
    xlim(-2, 30.4) +
    theme_custom() +
    theme(legend.position = legend_position)
}

# Function to create a simple custom bar plot
bar_plot <- function(
    df, x, y,
    fill = NULL,
    fill_pallete = NULL,
    width = NULL,
    fill_lab = NULL,
    colour_lab = NULL,
    title_lab = NULL,
    subtitle_lab = NULL,
    x_lab, y_lab = "Number of samples",
    legend_position = "bottom",
    legend_direction = "horizontal") {
  gg <- ggplot({{ df }}, aes(x = {{ x }}, y = {{ y }}, fill = {{ fill }})) +
    geom_bar(stat = "identity", colour = "#3a3838", width = width) +
    labs(
      title = title_lab,
      subtitle = subtitle_lab,
      fill = fill_lab,
      colour = colour_lab,
      x = x_lab,
      y = y_lab
    ) +
    # theme_bw() +
    # theme_minimal() +
    scale_fill_manual(values = fill_pallete) +
    # scale_fill_manual(values = viridis_color_pallete) +
    theme_custom() +
    theme(
      legend.position = legend_position,
      legend.direction = legend_direction
    )
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
  # guides(fill = "none")
  return(gg)
}

# Function to create a simple custom histogram plot
# histogram_plot <- function(
#     df, y, fill = NULL, binwith = 1,
#     x_lab = "Number of samples", y_lab) {
#   # FIXME: gg <- ggplot({{ df }}, aes(y = {{ y }}, fill = as.factor({{ fill }}))) +
#   gg <- ggplot({{ df }}, aes(y = {{ y }}, fill = {{ fill }})) +
#     geom_histogram(binwidth = binwith, colour = "black") +
#     labs(
#       x = x_lab,
#       y = y_lab
#     ) +
#     # theme_bw() +
#     # theme_minimal() +
#     theme_custom() +
#     # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
#     guides(fill = "none")
#   return(gg)
# }

histogram_plot <- function(
    df, y, fill = NULL, binwith = 1,
    x_lab = "Number of samples", y_lab) {
  # FIXME: gg <- ggplot({{ df }}, aes(y = {{ y }}, fill = as.factor({{ fill }}))) +
  gg <- ggplot({{ df }}, aes(y = {{ y }}, fill = {{ fill }})) +
    geom_histogram(binwidth = binwith, colour = "black") +
    labs(
      x = x_lab,
      y = y_lab
    ) +
    # theme_bw() +
    # theme_minimal() +
    theme_custom() #+
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
  # guides(fill = "none")
  return(gg)
}

histogram_plot_x_or_y <- function(
    df, aes_var,
    x = NULL, y = NULL, fill = NULL, binwidth = 1,
    x_lab = "Number of samples", y_lab = "Number of samples") {
  aes_mapping <- if ({{ aes_var }} == "x") {
    aes(x = {{ x }}, fill = {{ fill }})
  } else if ({{ aes_var }} == "y") {
    aes(y = {{ y }}, fill = {{ fill }})
  } else {
    stop("Invalid value for aes_var. Use 'x' or 'y'.")
  }

  gg <- ggplot({{ df }}, aes_mapping) +
    geom_histogram(binwidth = binwidth, colour = "black") +
    labs(
      x = x_lab,
      y = y_lab
    ) +
    theme_custom() #+
  # guides(fill = "none")

  return(gg)
}


### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
cat("Done loading script!!!")
cat("Woooooooohooooooo!!!")
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
