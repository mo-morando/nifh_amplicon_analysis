#!/usr/bin/env Rscript

#' @title World Map Plot Generator for nifH Database
#' @description This script processes nifH amplicon data and generates comprehensive world map visualizations for each study in the nifH database. It provides a flexible and robust approach to data visualization with error handling and informative logging.
#'
#' @details The script performs the following main steps:
#' * Loads required libraries (tidyverse, rlang, patchwork, argparser)
#' * Sources necessary external utility functions
#' * Sets up and parses command-line arguments for input/output paths and plot settings
#' * Defines custom themes and plotting functions for world maps
#' * Loads and processes input data, filtering for photic zone samples
#' * Generates color-blind safe palettes for plot aesthetics
#' * Creates world map plots for all studies and faceted by ocean
#' * Saves generated plots in specified formats
#'
#' Key functions include:
#' * create_custom_worldmap_plot(): Creates a customized world map plot
#' * extract_study_ids_and_filter_photic(): Processes input data for plotting
#' * create_and_print_world_maps(): Generates main and faceted world map plots
#' * main(): Orchestrates the entire data processing and visualization workflow
#'
#' @usage Rscript wrldMaps_nifhDB.R [--files FILES] [--input_path PATH] [--output_path PATH] [--plot_ext EXT] [--plot_device DEVICE]
#'
#' @param --files Comma-separated list of input files (default: cmap_coloc)
#' @param --input_path Directory path for input files (default: ../analysis/out_files)
#' @param --output_path Directory path for output plots (default: ../analysis/plots)
#' @param --plot_ext File extension for saved plots (default: .jpeg)
#' @param --plot_device Device to use for plot generation (default: jpeg)
#'
#' @author Michael (Mo) Morando
#'  2025-02-02
#'
#' @note This script requires the following R packages: tidyverse, rlang, patchwork, argparser
#'
#' @examples
#' Rscript wrldMaps_nifhDB.R --files cmap_coloc --input_path ../data/processed --output_path ../results/plots --plot_ext .png --plot_device png
#'
#' @export



# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(rlang)
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
  tryCatch(
    {
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
    },
    error = function(e) {
      stop("Error in source_file:", conditionMessage(e))
    } # ,
    # warning = function(w) {
    #   warning(paste("Warning in source_file:", conditionMessage(w)))
    # }
  )
}


# Source needed files
files_to_source <- c(
  "functions.R",
  "basic_plotting.R"
)


#' Set up the argument parser
#'
#' This function creates and configures an argument parser for generating world plots
#'
#' @return A configured argument parser object
#' @importFrom argparser arg_parser add_argument
#' @export
#'
#' @examples
#' parser <- setup_parser()
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
  tryCatch(
    {
      # Parse the arguments
      argv <- parse_args(parser)

      return(list(
      # Convert the comma-separated string to a vector
      files_to_read = strsplit(argv$files, ",")[[1]],
      files_in_path = argv$input_path,
      files_out_path = argv$output_path,
      plot_ext = argv$plot_ext,
      plot_device = argv$plot_device
      ))
      # return(list(
      #   files_to_read = files_to_read,
      #   files_in_path = files_in_path,
      #   files_out_path = files_out_path,
      #   plot_ext = plot_ext,
      #   plot_device = plot_device
      # ))
    },
    # warning = function(w) {
    #   cat("Warning occurred:\n")
    #   cat("Warning message:", conditionMessage(w), "\n")
    #   cat("Warning call:", deparse(conditionCall(w)), "\n")
    # },
    error = function(e) {
      cat("Error call in:", deparse(conditionCall(e)), "\n")
      stop("Error parsing arguments due to:\n", conditionMessage(e))
    }
  )
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


#' Create a custom world map plot
#'
#' @param data Data frame containing plot data
#' @param x_var Variable for x-axis (longitude)
#' @param y_var Variable for y-axis (latitude)
#' @param color_var Variable for color coding points
#' @param legend_nrow Number of rows in the legend
#' @param world World map data
#' @param fill_palette Color palette to use
#' @param legend_position Position of the legend
#' @return A ggplot object
create_custom_worldmap_plot <- function(
    data,
    x_var,
    y_var,
    color_var,
    legend_nrow,
    fill_palette = NULL,
    legend_position = "bottom") {
  tryCatch(
    {
      tryCatch(
        {
          # Load world map data
          world <- map_data("world")
        },
        error = function(e) {
          cat("Error call in:", deparse(conditionCall(e)), "\n")
          stop(
            "World map from map_data could not load due to:\n",
            conditionMessage(e)
          )
        }
      )

      # withCallingHandlers({
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
    },
    #   warning = function(w) {
    #     cat("Warning occurred during create_custom_worldmap_plot:\n")
    #     cat("Warning message:", conditionMessage(w), "\n")
    #     cat("Warning call:", deparse(conditionCall(w)), "\n")
    #     invokeRestart("muffleWarning")
    #   }
    # )
    # },
    error = function(e) {
      cat("Error call in:\n", deparse(conditionCall(e)), "\n")
      stop("No plots created due to:\n", conditionMessage(e), "\n")
    }
  )
}


### write function where you do not need quotes to enter information
filter_df <- function(x, plt_flt) {
  tryCatch(
    {
      filtered_data <- x %>%
        filter({{ plt_flt }})

      cat("Filtered df for photic zone samples only. Rows remaining:", nrow(filtered_data), "\n")
      return(filtered_data)
    },
    error = function(e) {
      cat("Error call in:", deparse(conditionCall(e)), "\n")
      stop("Filtered df was not created due to:\n", conditionMessage(e), "\n")
    }
  )
}


#' Generate a color-blind safe palette
#'
#' @param n Number of colors to generate
#' @return Vector of color codes
generate_cb_palette <- function(n) {
  tryCatch(
    {
      # Setting up pallete with color blind safe colors
      # Define the base palette
      cbpal_15 <- c(
        "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", "#490092", "#006ddb", "#b66dff", "#6db6ff",
        "#b6dbff", "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d"
      )

      # If n is less than or equal to 15, return the first n colors
      if (n <= length(cbpal_15)) {
        cat("Palette generated with", n, "colors.\n")
        return(cbpal_15[1:n])
      }


      cat("Palette generated with", n, "colors using color ramp.\n")
      return(colorRampPalette(cbpal_15)(n))
    },
    # warning = function(w) {
    #   cat("Warning occurred in generate_cb_palette:\n")
    #   cat("Warning message:", conditionMessage(w), "\n")
    #   cat("Warning call:", deparse(conditionCall(w)), "\n")
    # },
    error = function(e) {
      cat("Error call in: '", deparse(conditionCall(e)), "'\n")
      stop(
        "Error caused no color palette to be generated:\n",
        conditionMessage(e), "\n"
      )
    }
  )
}


#' Extract Study IDs and Filter Data for Photic Zone
#'
#' This helper function performs two main tasks:
#' 1. Extracts unique study IDs from the input data frame.
#' 2. Filters the data frame for entries in the photic zone.
#'
#' @param data A data frame containing at least 'studyID' and 'plt_flt' columns
#' @param photic_filter The value used to filter for photic zone in the 'plt_flt' column
#' @return A list containing study_ids and photic_data
#' @examples
#' result <- extract_study_ids_and_filter_photic(cmap_coloc, photic = "photic")
extract_study_ids_and_filter_photic <- function(data, photic_filter) {
  tryCatch(
    {
      # Extract unique study IDs
      study_ids <- data %>%
        distinct(studyID) %>%
        pull(studyID)
      cat("There are", length(study_ids), "study IDs:", paste(study_ids, collapse = ", "), "\n")

      # Filter data for photic zone
      photic_data <- filter_df(x = data, plt_flt = photic)

      return(list(
        study_ids = study_ids,
        photic_data = photic_data
      ))
    },
    error = function(e) {
      cat("Error call in:", deparse(conditionCall(e)), "\n")
      stop("No filtered df produced:\n", conditionMessage(e), "\n")
    }
  )
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
  cat("Creating world plots...\n")
  tryCatch(
    {
      # Create the main world map plot
      world_map_plot <- create_custom_worldmap_plot(
        data = data,
        x_var = {{ x_var }},
        y_var = {{ y_var }},
        color_var = {{ color_var }},
        legend_nrow = legend_nrow,
        fill_palette = fill_palette,
        legend_position = legend_position
      ) %>% suppressWarnings()

      # Create the faceted plot
      world_map_plot_faceted <- world_map_plot +
        facet_wrap(as.formula(paste("~", facet_var))) +
        theme(legend.position = legend_position)

      # Print the plots
      print(world_map_plot)
      print(world_map_plot_faceted)

      if (inherits(world_map_plot, "ggplot") && inherits(world_map_plot_faceted, "ggplot")) {
        # Print a statement with the parameters
        cat("Custom plots created with the following parameters:\n")
        cat("x_var:", deparse(substitute(x_var)), "\n")
        cat("y_var:", deparse(substitute(y_var)), "\n")
        color_var_expr <- enquo(color_var)
        cat("color_var:", as_name(color_var_expr), "\n")
        cat("facet_var:", facet_var, "\n")
        cat("legend_nrow:", legend_nrow, "\n")
        cat("legend_position:", legend_position, "\n")
      }

      # Return the plots
      return(list(main_plot = world_map_plot, faceted_plot = world_map_plot_faceted))
    },
    error = function(e) {
      cat("Error call in:", deparse(conditionCall(e)), "\n")
      stop("Plots were not created due to:\n", conditionMessage(e), "\n")
    }
  )
}


#' Main function to process data and generate plots
#'
#' @param files_to_read List of files to read
#' @param files_in_path Input directory path
#' @param plot_var Variable to use for plotting
#' @param files_out_path Output directory path
#' @param plot_ext Plot file extension
#' @param plot_device Plot device
main <- function(
    files_to_read,
    files_in_path,
    plot_var,
    files_out_path,
    plot_ext,
    plot_device) {
  cat("Starting main function...\n")

  tryCatch(
    {
      # Load the data
      data_list <- load_files(files_to_read, files_in_path)
      cat("Data loaded successfully.\n")

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
        message("Error: Plots did not get saved")
      }
    },
    error = function(e) {
      cat("Error in call:", deparse(conditionCall(e)), "\n")
      stop("Main function did not execute properly resulting in no plots being generated or saved due to:\n", conditionMessage(e), "\n")
    }
  )
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
