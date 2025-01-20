#' Module Docstring:
#'
#' This script contains a collection of functions designed to perform various data processing tasks using ggplot2 and tidyverse packages in R. Below are the main functionalities of the script:
#'
#' - Custom Theme Function:
#'     - `theme_custom`: Define a custom theme for ggplot2 plots.
#'
#' - Color Functions:
#'     - `get_viridis_colors`: Generate a custom color palette based on the viridis color scheme.
#'     - `nifh_cluster_colours`: Define custom colors for nifh clusters.
#'
#' - Plotting Functions:
#'     - `create_custom_plot_point`: Create a custom dot plot.
#'     - `bar_plot`: Create a simple custom bar plot.
#'     - `histogram_plot`: Create a simple custom histogram plot.
#'     - `histogram_plot_x_or_y`: Create a custom histogram plot with options to specify x or y axis.
#'
#' These functions can be used individually to customize plot aesthetics and create visually appealing visualizations.
#'
#' @import tidyverse
#' @importFrom viridisLite viridis
#'

#' Define custom colors for nifh clusters
#'
#' This function defines custom colors for nifh clusters.
#'
#' @return A named vector of custom colors.
#'
nifh_cluster_colours <- c(
  "1A" = "steelblue",
  "1J/1K" = "chocolate4",
  "1O/1P" = "magenta",
  "3" = "red",
  "1G" = "chocolate1",
  "1B" = "green2",
  "other" = "lightgrey",
  # "unknown" = "darkgrey",
  "4" = "black",
  "2" = "blue"
)


# Function to define the custom theme
eb <- element_blank()

#' Define a custom theme for ggplot2 plots
#'
#' This function defines a custom theme for ggplot2 plots, including modifications to plot titles, axes, grids, and legends.
#'
#' @return A ggplot2 theme object.
#'
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
      # legend.title = element_text(size = 20), # size of legend title
      # legend.key.size = unit(2, "lines"), # size of legend keys
      # legend.text = element_text(size = 17), # size of legend text
      # legend.spacing.x = unit(0.5, "lines"), # horizontal spacing between legend elements
      # legend.spacing.y = unit(0.5, "lines"), # vertical spacing between legend elements
      legend.box = "horizontal",
      legend.direction = "horizontal",
      axis.text.y = element_text(size = 21, face = "bold"),
      axis.text.x = element_text(size = 21, face = "bold")
    )
}

#' Generate custom color palette based on viridis color scheme
#'
#' This function generates a custom color palette based on the viridis color scheme.
#'
#' @param df_to_plot Dataframe to plot.
#' @param plot_by Variable to plot.
#' @param pallette_type Type of viridis color palette (e.g., "D").
#' @param color_direction Direction of color scheme.
#' @param begin_point Starting point of the color scheme.
#'
#' @return A vector of custom color palette.
#'
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

#' Generate a colorblind-safe palette for nifH clusters
#'
#' This function generates a colorblind-safe palette based on the RColorBrewer "Paired" palette
#' and assigns the colors to nifH cluster names.
#'
#' @param nifh_cluster_colours A named vector of nifH cluster names and their corresponding colors.
#'
#' @return A named vector of colorblind-safe colors assigned to nifH cluster names.
#'
#' @examples
#' nifh_cluster_colours <- c("Cluster1" = "#FF0000", "Cluster2" = "#00FF00", "Cluster3" = "#0000FF")
#' nifh_cluster_colours_colbldsafe <- generate_nifh_palette(nifh_cluster_colours)
#'
#' @import RColorBrewer
#'
#' @export
generate_nifh_palette <- function(nifh_cluster_colours) {
  # Generate palette
  pal <- c(RColorBrewer::brewer.pal(12, "Paired")[c(8, 2, 11, 9, 5, 4, 10, 3)], "#777777")

  # Assign colors to nifH cluster names
  names(pal) <- names(nifh_cluster_colours)

  return(pal)
}


#' Create a custom dot plot
#'
#' This function creates a custom dot plot using ggplot2.
#'
#' @param data1 Dataframe containing points for the plot.
#' @param data2 Dataframe containing labels for the points.
#' @param x1 X-axis variable for data1.
#' @param y1 Y-axis variable for data1.
#' @param x2 X-axis variable for data2.
#' @param y2 Y-axis variable for data2.
#' @param colour Variable for point color.
#' @param group Variable for grouping points.
#' @param size Variable for point size.
#' @param label Variable for point labels.
#' @param legend_position Position of the legend.
#'
#' @return A ggplot2 object representing the custom dot plot.
#'
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

#' Create a simple custom bar plot
#'
#' This function creates a simple custom bar plot using ggplot2.
#'
#' @param df Dataframe containing plot data.
#' @param x X-axis variable.
#' @param y Y-axis variable.
#' @param fill Fill variable for bars.
#' @param fill_pallete Color palette for fill variable.
#' @param width Width of bars.
#' @param fill_lab Label for fill legend.
#' @param colour_lab Label for color legend.
#' @param title_lab Title label for the plot.
#' @param subtitle_lab Subtitle label for the plot.
#' @param x_lab X-axis label.
#' @param y_lab Y-axis label.
#' @param legend_position Position of the legend.
#' @param legend_direction Direction of the legend.
#' @param print_out Boolean, determines if plot is printed. Default is FALSE
#'
#' @return A ggplot2 object representing the custom bar plot.
#'
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
    legend_direction = "horizontal",
    n_row = 1,
    x_axis_angle = FALSE,
    facet = FALSE,
    print_out = FALSE) {
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
    scale_fill_manual(values = fill_pallete) +
    theme_custom() +
    theme(
      legend.position = legend_position,
      legend.direction = legend_direction,
      legend.key.size = unit(2, "lines"), # size of legend keys
      legend.text = element_text(size = 16), # size of legend text
      legend.title = element_text(size = 25), # size of legend title
      legend.spacing.x = unit(0.5, "lines"), # horizontal spacing between legend elements
      legend.spacing.y = unit(0.5, "lines") # vertical spacing between legend elements
    )

  # Conditionally set x-axis angle
  if (x_axis_angle) {
    gg <- gg +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 17, face = "bold"),
        axis.title.x = element_blank()
      )
  }

  gg <- gg + guides(fill = guide_legend(nrow = n_row))

  if (facet) {
    gg <- gg + facet_wrap(~facet_by)
  }

  if (print_out) {
    print(gg) # print plot if prompted
  }

  return(gg)
}

#' Create a simple custom histogram plot
#'
#' This function creates a simple custom histogram plot using ggplot2.
#'
#' @param df Dataframe containing plot data.
#' @param y Y-axis variable.
#' @param fill Fill variable for bars.
#' @param binwidth Width of bins for histogram.
#' @param x_lab X-axis label.
#' @param y_lab Y-axis label.
#'
#' @return A ggplot2 object representing the custom histogram plot.
#'
histogram_plot <- function(
    df, y, fill = NULL, binwidth = 1,
    x_lab = "Number of samples", y_lab) {
  gg <- ggplot({{ df }}, aes(y = {{ y }}, fill = {{ fill }})) +
    geom_histogram(binwidth = binwidth, colour = "black") +
    labs(
      x = x_lab,
      y = y_lab
    ) +
    theme_custom()
  return(gg)
}

#' Create a custom histogram plot with options to specify x or y axis
#'
#' This function creates a custom histogram plot with options to specify x or y axis using ggplot2.
#'
#' @param df Dataframe containing plot data.
#' @param aes_var Variable to map to x or y axis.
#' @param x X-axis variable.
#' @param y Y-axis variable.
#' @param fill Fill variable for bars.
#' @param binwidth Width of bins for histogram.
#' @param x_lab X-axis label.
#' @param y_lab Y-axis label.
#'
#' @return A ggplot2 object representing the custom histogram plot.
#'
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
    theme_custom()

  return(gg)
}


## function for plotting
scatter_line_plot <- function(
    df,
    x,
    y,
    colour = NULL,
    colour_pallete = NULL,
    group = NULL,
    colour_lab = NULL,
    fill_lab = NULL,
    fill_pallete = NULL,
    fill = NULL,
    title_lab = NULL,
    subtitle_lab = NULL,
    shape_lab = NULL,
    x_lab,
    y_lab = "Number of samples",
    legend_position = "bottom",
    legend_direction = "horizontal") {
  plot <- ggplot(
    df,
    aes(
      y = {{ y }},
      x = {{ x }},
      colour = {{ colour }},
      group = {{ group }},
      # shape = hemi
    )
  ) +
    geom_point(
      aes(
        # shape=df_type
      ),
      na.rm = T, size = 2, stroke = 1.5, # fill = 'black'
      shape = 1 ## _# you have to comment this out if you want to do anything with the shapes, e.g., differentiate hemispheres
    ) +
    geom_smooth(
      aes(
        fill = {{ fill }}
      ),
      method = "loess", formula = y ~ x, se = T, alpha = 0.3, show.legend = F, linewidth = 2
    ) +
    scale_color_manual(values = colour_pallete) +
    scale_fill_manual(values = fill_pallete) +
    labs(
      title = title_lab,
      subtitle = subtitle_lab,
      # x = expression(bold(paste(x_lab))),
      x = eval(substitute(expression(bold(paste(x_lab))))),
      y = eval(substitute(expression(bold(paste(y_lab))))),
      colour = colour_lab,
      shape = shape_lab,
      fill = fill_lab
    ) +
    guides(
      colour = guide_legend(
        nrow = 1, byrow = TRUE,
        override.aes = list(size = 8, stroke = 5)
      ),
      shape = guide_legend(nrow = 1, byrow = TRUE),
      # shape = guide_legend(override.aes = list(size = 10))
    ) +
    ylim(0, 1.0) +
    theme_custom() +
    theme(
      legend.position = legend_position,
      legend.direction = legend_direction,
      legend.key.size = unit(3, "lines"), # size of legend keys
      # legend.key.height = unit(3, "cm"),
      # legend.key.width = unit(3, "cm"),
      legend.text = element_text(size = 20), # size of legend text
      legend.title = element_text(size = 26), # size of legend title
      legend.spacing.x = unit(0.5, "lines"), # horizontal spacing between legend elements
      legend.spacing.y = unit(0.5, "lines") # vertical spacing between legend elements
    )

  return(plot)
}


#' Save a ggplot object as an image file
#'
#' This function saves a specified ggplot object as an image file with customizable parameters.
#' If no plot is specified, it will save the last plot created.
#'
#' @param plot_obj The ggplot object to save. If not provided, the last plot will be used.
#' @param output_dir The directory where the file will be saved. Default is 'files_out_path'.
#' @param filename The name of the output file (without extension). Default is "Samples_per_studyID_ocean_fill_bar_facet_nucelic".
#' @param file_ext The file extension for the output. Default is 'plot_ext'.
#' @param height The height of the output image in inches. Default is 8.5.
#' @param width The width of the output image in inches. Default is 14.
#' @param dpi The resolution of the output image. Default is 300.
#' @param device The device to use for plotting. Default is NULL (automatically chosen based on file extension).
#'
#' @return NULL
#'
#' @examples
#' # Assuming you have a ggplot object named 'my_plot'
#' save_custom_plot(plot = my_plot)
#' save_custom_plot(plot = my_plot, height = 10, width = 16, dpi = 600)
#' 
#' # To save the last plot created
#' save_custom_plot()
#'
#' @export
save_custom_plot <- function(plot_obj = NULL,
                             output_dir,
                             filename,
                             file_ext = plot_ext,
                             height = 8.5,
                             width = 14,
                             dpi = 300,
                             device = NULL) {
  
  # Use the provided plot or the last plot if not specified
  if (is.null(plot_obj)) {
    plot_obj <- last_plot()
    cat("No plot specified. Using the last plot created.\n")
  }
  
  # Construct the full file path
  file_path <- file.path(output_dir, paste0(filename, file_ext))
  
  # Save the plot
  ggsave(plot = plot_obj,
         filename = file_path,
         height = height,
         width = width,
         units = "in",
         dpi = dpi,
         device = get(device))
  
  # Print a success message
  cat("Plot saved successfully: ", file_path, "\n")
}


### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
cat("Done loading script!!!\n")
cat("Woooooooohooooooo!!!\n\n")
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
