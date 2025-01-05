#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(argparser)
})

#' Source a file with error handling and path validation
#'
#' @param file_path Character string specifying the path to the file to be sourced
#' @return Invisible NULL. Prints status messages.
#' @examples
#' source_file("/path/to/your/file.R")
source_file <- function(file_path) {
  for (file in file_path) {
    # Check if file exists
    if (!file.exists(file)) {
      stop("File does not exist: ", file, "\n")
    }

    # Try to source the file
    tryCatch(
      {
        source(file)
        cat("Successfully sourced : ", file, "\n")
      },
      error = function(e) {
        stop("Error sourcing : ", file_path, ":", conditionMessage(e), "\n")
      },
      warning = function(w) {
        cat("Warning while sourcing : ", file_path, ":", conditionMessage(w), "\n")
      }
    )
  }

  cat("Finished sourcing files. \n")
}


# Source needed files
files_to_source <- c(
  "/Users/mo/Projects/nifH_amp_project/myWork/scripts/functions.R",
  "/Users/mo/Projects/nifH_amp_project/myWork/scripts/basic_plotting.R"
)

setup_parser <- function() {
  parser <- arg_parser("Process input files and generate world plot of each study for nifH database")
  parser <- add_argument(
    parser = parser, "--files",
    help = "CVS list of files to read in",
    default = "cmap_coloc"
  )
  parser <- add_argument(
    parser = parser,
    arg = "--input_path",
    help = "Input directory path",
    default = "../analysis/out_files"
  )
  parser <- add_argument(
    parser = parser,
    arg = "--output_path",
    help = "Output directory path",
    default = "../analysis/plots"
  )
  parser <- add_argument(parser,
    arg = "--plot_ext",
    help = "Extension to add to each plot, setting the type of figure produced",
    default = ".jpeg"
  )
  parser <- add_argument(parser,
    arg = "--plot_device",
    help = "Device to add to each plot, setting how the figure is produced.
    Mainly important if generating PDF figures",
    default = "jpeg"
  )

  return(parser)
}

# cmap_coloc <- read_csv("/Users/mo/Projects/nifH_amp_project/myWork/analysis/out_files/cmap_coloc.csv")

#' Parse command-line arguments
#'
#' Parses the command-line arguments using the provided parser.
#'
#' @param parser An argument parser object created by setup_parser()
#'
#' @return A list containing:
#'    \item{files_to_read}{A character vector of file names to read}
#'    \item{files_in_path}{The input directory path}
#'    \item{files_to_read}{The output directory path}
#' @importFrom argparser parse_args
#' @export
#'
#' @examples
#' parser <- setup_parser()
#' args <- parse_arg(parser)
parse_arg <- function(parser) {
  # Parse the arguments
  argv <- parse_args(parser)

  # Convert the comma-separated string to a vector
  files_to_read <- strsplit(argv$files, ",")[[1]]
  files_in_path <- argv$input_path
  files_out_path <- argv$output_path
  plot_ext <- argv$plot_ext
  plot_device <- argv$plot_device

  return(list(
    files_to_read = files_to_read,
    files_in_path = files_in_path,
    files_out_path = files_out_path,
    plot_ext = plot_ext,
    plot_device = plot_device
  ))
}


#' Define a custom theme for ggplot2 plots
#'
#' This function defines a custom theme for ggplot2 plots, including modifications to plot titles, axes, grids, and legends.
#'
#' @return A ggplot2 theme object.
#'
theme_custom_maps <- function() {
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 13, face = "bold"),
    panel.border = element_rect(colour = "black"),
    axis.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 25, face = "bold", colour = "black"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 13, face = "bold"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold")
  )
}


### - load function for plotting world maps
create_custom_worldmap_plot <- function(
    data,
    x_var,
    y_var,
    color_var,
    legend_nrow,
    world,
    fill_palette = NULL,
    legend_position = "bottom") {
  plot <- ggplot(data) +
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "black", fill = "black", size = 0.2
    ) +
    geom_point(
      aes(
        x = {{ x_var }}, y = {{ y_var }}, color = {{ color_var }}
      ),
      data = data # %>% mutate(studyID = str_replace_all(studyID, "_", " "))
      ,
      alpha = 0.15,
      size = 4, shape = 1, stroke = 3
    ) +
    # scale_color_manual(values = getPalette(palette_count)) +
    scale_fill_manual(values = fill_palette) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      y = "Latitude",
      x = "Longitude",
      colour = "Study ID"
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank()) +
    theme_custom_maps() +
    theme(legend.position = legend_position) +
    guides(colour = guide_legend(nrow = legend_nrow, byrow = TRUE, override.aes = list(alpha = 1))) +
    coord_quickmap()

  return(plot)
}


### write function where you do not need quotes to enter information
filter_df <- function(x, plt_flt) {
  cat("Input Variables:\n")
  cat("x:", deparse(substitute(x)), "\n")
  cat("plt_flt:", deparse(substitute(plt_flt)), "\n")

  filtered_data <- x %>%
    filter({{ plt_flt }})

  return(filtered_data)
}



# Function to generate color palette
generate_cb_palette <- function(n) {
  # Setting up pallete with color blind safe colors
  # Define the base palette
  cbpal_15 <- c(
    "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", "#490092", "#006ddb", "#b66dff", "#6db6ff",
    "#b6dbff", "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d"
  )

  # If n is less than or equal to 15, return the first n colors
  if (n <= length(cbpal_15)) {
    cat("Palette generated!\n")
    return(cbpal_15[1:n])
  }


  cat("Palette generated\n")
  return(colorRampPalette(cbpal_15)(n))
}


#' Extract Study IDs and Filter Data for Photic Zone
#'
#' This helper function performs two main tasks:
#' 1. Extracts unique study IDs from the input data frame.
#' 2. Filters the data frame for entries in the photic zone.
#'
#' @param data A data frame containing at least 'studyID' and 'plt_flt' columns.
#' @param photic_filter The value used to filter for photic zone in the 'plt_flt' column.
#'
#' @return A list containing two elements:
#'   - study_ids: A vector of unique study IDs.
#'   - photic_data: A data frame filtered for the photic zone.
#'
#' @examples
#' result <- extract_study_ids_and_filter_photic(cmap_coloc, photic = "photic")
#' print(result$study_ids)
#' head(result$photic_data)
extract_study_ids_and_filter_photic <- function(data, photic_filter) {
  # Extract unique study IDs
  study_ids <- data %>%
    distinct(studyID) %>%
    pull(studyID)

  cat("The study IDs are:", paste(study_ids, collapse = ", "), "\n")

  # Filter data for photic zone
  # photic_data <- data %>%
  #   filter(plt_flt == photic_filter)

  photic_data <- filter_df(x = data, plt_flt = photic)

  cat("Filtered data for photic zone. Rows remaining:", nrow(photic_data), "\n")

  return(list(
    study_ids = study_ids,
    photic_data = photic_data
  ))
}

#' Create and Print World Map Plots
#'
#' This function creates a world map plot and a faceted version by ocean,
#' then prints both plots and a summary of the parameters used.
#'
#' @param data A data frame containing the plot data
#' @param x_var The variable to use for x-axis (longitude)
#' @param y_var The variable to use for y-axis (latitude)
#' @param color_var The variable to use for color coding points
#' @param legend_nrow Number of rows in the legend
#' @param fill_palette Color palette to use
#' @param legend_position Position of the legend
#' @param facet_var Variable to use for faceting (default is "ocean")
#'
#' @return A list containing the main plot and the faceted plot
#'
#' @examples
#' result <- create_and_print_world_maps(
#'   data = df_photic,
#'   x_var = lon,
#'   y_var = lat,
#'   color_var = studyID,
#'   legend_nrow = 3,
#'   fill_palette = cb_palette,
#'   legend_position = "bottom"
#' )
create_and_print_world_maps <- function(data, x_var, y_var, color_var, legend_nrow, fill_palette, legend_position, facet_var = "ocean") {
  cat("Creating plots...\n")
  # Load world map data
  world <- map_data("world")

  # Create the main world map plot
  world_map_plot <- create_custom_worldmap_plot(
    data = data,
    x_var = {{ x_var }},
    y_var = {{ y_var }},
    color_var = {{ color_var }},
    legend_nrow = legend_nrow,
    fill_palette = fill_palette,
    legend_position = legend_position,
    world = world
  ) %>% suppressWarnings()

  # Create the faceted plot
  world_map_plot_faceted <- world_map_plot +
    facet_wrap(as.formula(paste("~", facet_var))) +
    theme(legend.position = legend_position)

  # Print the plots
  print(world_map_plot)
  print(world_map_plot_faceted)

  # Print a statement with the parameters
  cat("Custom plots created with the following parameters:\n")
  cat("x_var:", deparse(substitute(x_var)), "\n")
  cat("y_var:", deparse(substitute(y_var)), "\n")
  cat("color_var:", deparse(substitute(color_var)), "\n")
  cat("facet_var:", facet_var, "\n")
  cat("legend_nrow:", legend_nrow, "\n")
  cat("legend_position:", legend_position, "\n")

  # Return the plots
  return(list(main_plot = world_map_plot, faceted_plot = world_map_plot_faceted))
}


main <- function(
    files_to_read,
    files_in_path,
    plot_var,
    files_out_path,
    plot_ext,
    plot_device) {
  # Load the data
  data_list <- load_files(files_to_read, files_in_path)

  df_photic <- extract_study_ids_and_filter_photic(
    data = data_list$cmap_coloc,
    photic_filter = "photic"
  )

  # Generate colorblind safe pallette
  cb_palette <- generate_cb_palette(
    nrow(df_photic$photic_data %>% distinct({{ plot_var }}))
  )

  plot_results <- create_and_print_world_maps(
    data = df_photic$photic_data,
    x_var = lon,
    y_var = lat,
    color_var = {{ plot_var }},
    legend_nrow = 3,
    fill_palette = cb_palette,
    legend_position = "bottom"
  )


  if (length(plot_results) == 2) {
    save_custom_plot(
      plot = plot_results$main_plot,
      output_dir = files_out_path,
      filename = "worldStnMap_nifhdb_allstudies",
      file_ext = plot_ext,
      height = 8.5,
      width = 14,
      dpi = 300,
      device = plot_device
    )

    save_custom_plot(
      plot = plot_results$faceted_plot,
      output_dir = files_out_path,
      filename = "worldStnMap_nifhdb_allstudies_fct_ocean",
      file_ext = plot_ext,
      height = 8.5,
      width = 14,
      dpi = 300,
      device = plot_device
    )
  } else {
    message("Plots did not get saved")
  }
}


# Run if the script is being called directly
if (sys.nframe() == 0 && !interactive()) {
  source_file(files_to_source)

  parser <- setup_parser()
  args <- parse_arg(parser)

  validate_parsed_args(args)

  # Create output directory if it doesn't exist
  create_dir(args$files_out_path)

  main(
    files_to_read = args$files_to_read,
    files_in_path = args$files_in_path,
    plot_var = studyID,
    files_out_path = args$files_out_path,
    plot_ext = args$plot_ext,
    plot_device = args$plot_device
  )
}
