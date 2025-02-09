#!/usr/bin/env Rscript

#' @title NifH Amplicon Data Pipeline Basic Plots
#' @description This script processes NifH amplicon data and generates a comprehensive set of visualizations to analyze sample distributions, oceanographic parameters, and other relevant metrics. It provides a modular and flexible approach to data visualization with robust error handling.
#'
#' @details The pipeline executes the following main steps:
#' * Loads and validates input files from specified paths
#' * Generates sample distribution plots (e.g., by type, nucleic acid, study ID, ocean, hemisphere)
#' * Creates oceanographic data plots (e.g., SST, PO4, Fe, PP, CHL)
#' * Produces a combined plot for key metrics (latitude, ocean, SST, PO4)
#' * Saves all generated plots in specified formats
#'
#' Key functions include:
#' * generate_sample_distribution_plots(): Creates various sample distribution visualizations
#' * generate_oceanographic_plots(): Produces plots for oceanographic parameters
#' * create_combined_plot(): Assembles a composite plot from individual visualizations
#' * main(): Orchestrates the entire visualization workflow
#'
#' @usage Rscript plots_out.R [--files FILES] [--input_path PATH] [--output_path PATH] [--tables_out_path PATH] [--plot_ext EXT] [--plot_device DEVICE]
#'
#' @param --files Comma-separated list of input files
#' @param --input_path Directory path for input files
#' @param --output_path Directory path for output plots
#' @param --tables_out_path Directory path for output tables
#' @param --plot_ext File extension for saved plots
#' @param --plot_device Device to use for plot generation
#'
#' @author Michael Morando
#' @date 2025-02-02
#'
#' @note This script requires the following R packages: tidyverse, patchwork, argparser
#'
#' @examples
#' Rscript plots_out.R --files cmap_coloc,samples_per_nucleicAcidType,samples_per_nucleicAcidType_studyid,samples_per_photic,samples_per_photic_nucacid,sample_type,unique_sample_id_key --input_path ../data/processed --output_path ../results/plots --tables_out_path ../results/tables --plot_ext .png --plot_device png
#'
#' @export


# Load required libraries with error handling
suppressPackageStartupMessages({
  tryCatch(
    {
      library(tidyverse)
      library(patchwork)
      library(argparser)
    },
    error = function(e) {
      cat("Error loading required packages:", conditionMessage(e), "\n")
      cat("Error call in:", deparse(conditionCall(e)), "\n")
    }
  )
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
      stop(paste("Error in source_file:", conditionMessage(e)))
    }
  )
}

# Source needed files
files_to_source <- c(
  "functions.R",
  "basic_plotting.R"
)

#' Set up the argument parser
#'
#' This function creates and configures an argument parser for processing workspace files.
#'
#' @return A configured argument parser object
#' @importFrom argparser arg_parser add_argument
#' @export
#'
#' @examples
#' parser <- setup_parser()
setup_parser <- function() {
  parser <- arg_parser("Process the files needed to make primary plots")
  parser <- add_argument(parser, "--files",
    help = "CVS list of files to read in",
    default = "cmap_coloc,samples_per_nucleicAcidType,samples_per_nucleicAcidType_studyid,samples_per_photic,samples_per_photic_nucacid,sample_type,unique_sample_id_key"
  )
  parser <- add_argument(parser,
    arg = "--input_path",
    help = "Input directory path",
    default = "../analysis/out_files"
  )
  parser <- add_argument(parser,
    arg = "--output_path",
    help = "Output directory path",
    default = "../analysis/plots"
  )

  parser <- add_argument(parser,
    arg = "--tables_out_path",
    help = "Output directory for tables",
    default = "../analysis/out_files/"
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
#'    \item{plot_ext}{Extension to add to each plot}
#' @importFrom argparser parse_args
#' @export
#'
#' @examples
#' parser <- setup_parser()
#' args <- parse_arg(parser)
parse_arg <- function(parser) {
  # Parse the arguments
  argv <- parse_args(parser)


  return(list(
    files_to_read = strsplit(argv$files, ",")[[1]],
    files_in_path = argv$input_path,
    files_out_path = argv$output_path,
    plot_ext = argv$plot_ext,
    plot_device = argv$plot_device,
    tables_out_path = argv$tables_out_path
  ))
}


#' Generate Sample Distribution Plots
#'
#' This function creates various plots to visualize the distribution of samples
#' across different categories such as sample type, nucleic acid type, study ID,
#' ocean, hemisphere, month, latitude, season, and depth.
#'
#' @param query_df A data frame containing the query results
#' @param files_out_path The output directory for saving plots
#' @param plot_ext The file extension for saved plots
#' @param plot_device The device to use for plotting
#'
#' @return NULL
#'
#' @details
#' The function creates the following plots:
#' - Sample type distribution
#' - Combined nucleic acid type, light level, and sample type distribution
#' - Samples per study ID
#' - Samples per ocean and hemisphere
#' - Samples per month
#' - Samples per month by study ID
#' - Samples distribution by absolute latitude
#' - Samples distribution by season
#' - Samples distribution by depth (above and below 150m)
#'
#' @export
generate_sample_distribution_plots <- function(
    query_df, samples_per_nucleicAcidType, samples_per_nucleicAcidType_studyid,
    samples_per_photic, samples_per_photic_nucacid,
    sample_type, files_out_path, plot_ext, plot_device) {
  tryCatch(
    {
      # Sample type plot
      tryCatch(
        {
          # sample_type <- count_and_arrange(query_df, c("sample_type", "new_group"))
          viridis_color_pallete <- get_viridis_colors(sample_type, new_group, "H", -1, 0)

          sample_type_plot <- bar_plot(
            df = sample_type,
            x = percentage,
            y = "sample type",
            fill = new_group,
            fill_pallete = viridis_color_pallete,
            y_lab = NULL,
            x_lab = "% of total",
            legend_position = "right",
            x_axis_angle = TRUE,
            n_row = 10
          )

          save_custom_plot(
            plot = sample_type_plot,
            output_dir = files_out_path,
            filename = "sample_type_distribution",
            file_ext = plot_ext,
            height = 10.5,
            width = 14,
            dpi = 300,
            device = plot_device
          )

          cat("Sample type plot created successfully.\n")
        },
        error = function(e) {
          cat("Error in creating sample type plot:", conditionMessage(e), "\n")
        }
      )

      # Combined nucleic acid type, light level, and sample type plot
      tryCatch(
        {
          combined_nuca_phtc_tibble <- samples_per_nucleicAcidType %>%
            mutate(tibble_id = "nucleicAcidType") %>%
            rename(new_group = nucleicAcidType) %>%
            bind_rows(samples_per_photic %>%
              mutate(
                tibble_id = "light level",
                photic = if_else(photic %in% "TRUE", "photic", "aphotic")
              ) %>%
              rename(new_group = photic))

          combined_nuca_phtc_smptp_tibble <- sample_type %>%
            bind_rows(combined_nuca_phtc_tibble)

          viridis_color_pallete <- get_viridis_colors(combined_nuca_phtc_smptp_tibble, new_group, "H", -1, 0)

          combined_plot <- suppressMessages(bar_plot(
            df = combined_nuca_phtc_smptp_tibble,
            y = tibble_id,
            x = percentage,
            fill = new_group,
            fill_pallete = viridis_color_pallete,
            y_lab = NULL,
            x_lab = "% of total",
            legend_position = "right",
            n_row = 10,
            x_axis_angle = TRUE
          ))

          save_custom_plot(
            plot = combined_plot,
            output_dir = files_out_path,
            filename = "samp_nucleic_acid_type_light_bar_cmb_plot",
            file_ext = plot_ext,
            height = 10.5,
            width = 14,
            dpi = 300,
            device = plot_device
          )

          cat("Combined nucleic acid type, light level, and sample type plot created successfully.\n")
        },
        error = function(e) {
          cat("Error in creating combined plot:", conditionMessage(e), "\n")
        }
      )

      # Study ID plot
      tryCatch(
        {
          samples_per_studyid <- count_and_arrange(query_df, c("studyID", "nucleicAcidType", "ocean"))
          viridis_color_pallete <- get_viridis_colors(samples_per_studyid, ocean, "magma", -1, 0)

          samples_per_studyid_plot <- bar_plot(
            df = samples_per_studyid,
            y = studyID,
            x = n,
            fill = ocean,
            fill_pallete = viridis_color_pallete,
            y_lab = "Study ID",
            x_lab = "Number of samples",
            legend_position = "right",
            legend_direction = "vertical",
            n_row = 5,
            print_out = TRUE
          ) + facet_wrap(~nucleicAcidType)

          save_custom_plot(
            plot = samples_per_studyid_plot,
            output_dir = files_out_path,
            filename = "Samples_per_studyID_ocean_fill_bar_facet_nucelic",
            file_ext = plot_ext,
            height = 8.5,
            width = 14,
            dpi = 300,
            device = plot_device
          )

          cat("Study ID plot created successfully.\n")
        },
        error = function(e) {
          cat("Error in creating Study ID plot:", conditionMessage(e), "\n")
        }
      )

      # Ocean and hemisphere plot
      tryCatch(
        {
          samples_per_ocean <- count_and_arrange(query_df, c("ocean", "hemi"))
          samples_per_ocean <- add_percentage(samples_per_ocean, n, percentage, grouping_by = NULL, remove_columns = "total")

          viridis_color_pallete <- get_viridis_colors(samples_per_ocean, hemi, "magma", -1, 0)

          samples_per_ocean_plot <- suppressMessages(bar_plot(
            df = samples_per_ocean,
            y = ocean,
            x = n,
            fill = hemi,
            fill_pallete = viridis_color_pallete,
            y_lab = "Ocean",
            x_lab = "Number of samples",
            fill_lab = "Hemisphere",
            legend_position = "right",
            legend_direction = "vertical",
            print_out = TRUE
          ) +
            scale_fill_manual(
              values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
              label = c("northern", "southern")
            ))

          save_custom_plot(
            plot = samples_per_ocean_plot,
            output_dir = files_out_path,
            filename = "Samples_per_ocean_hemi_fill_bar",
            file_ext = plot_ext,
            height = 8.5,
            width = 14,
            dpi = 300,
            device = plot_device
          )

          cat("Ocean and hemisphere plot created successfully.\n")
        },
        error = function(e) {
          cat("Error in creating Ocean and hemisphere plot:", conditionMessage(e), "\n")
        }
      )

      # Month plot
      tryCatch(
        {
          samples_per_month <- count_and_arrange(query_df, c("month"))

          samples_per_month_plot <- suppressMessages(bar_plot(
            df = samples_per_month,
            y = month,
            x = n,
            fill = NULL,
            y_lab = "Month",
            x_lab = "Number of samples",
            legend_position = "none",
            legend_direction = "vertical"
          ) +
            scale_y_discrete(
              limits = rev(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")),
              labels = rev(c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")),
            ))

          save_custom_plot(
            plot = samples_per_month_plot,
            output_dir = files_out_path,
            filename = "Samples_per_month_bar",
            file_ext = plot_ext,
            height = 8.5,
            width = 14,
            dpi = 300,
            device = plot_device
          )

          cat("Month plot created successfully.\n")
        },
        error = function(e) {
          cat("Error in creating Month plot:", conditionMessage(e), "\n")
        }
      )

      tryCatch(
        {
          ### month by study ID
          (samples_per_month_study_id <- count_and_arrange(query_df, c("month", "studyID")))

          viridis_color_pallete <- get_viridis_colors(samples_per_month_study_id, studyID, "H", -1, 0)

          samples_per_month_study_id_plot <- bar_plot(
            df = samples_per_month_study_id,
            y = month,
            x = n,
            fill = studyID,
            fill_pallete = viridis_color_pallete,
            y_lab = "Month",
            x_lab = "Number of samples",
            legend_position = "bottom"
          ) +
            # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
            scale_y_discrete(
              limits =
                rev(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")),
              labels =
                rev(c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")),
            )


          save_custom_plot(
            plot = samples_per_month_study_id_plot,
            output_dir = files_out_path,
            filename = "Samples_per_month_study_id_bar",
            file_ext = plot_ext,
            height = 8.5,
            width = 14,
            dpi = 300,
            device = plot_device
          )

          cat("Month plot by study ID created successfully.\n")
        },
        error = function(e) {
          cat("Error in creating Study ID Monthly plot:", conditionMessage(e), "\n")
        }
      )

      # Latitude plot
      tryCatch(
        {
          samples_per_lat_abs <- count_and_arrange(query_df, c("lat_abs", "hemi"))

          samples_per_lat_abs_plot <- suppressMessages(histogram_plot_x_or_y(samples_per_lat_abs, "x", lat_abs,
            fill = hemi, binwidth = 2,
            x_lab = "Absolute latitude"
          ) +
            theme(legend.position = "bottom") +
            labs(fill = "Hemisphere") +
            scale_fill_manual(
              values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
              label = c("northern", "southern")
            ))

          save_custom_plot(
            plot = samples_per_lat_abs_plot,
            output_dir = files_out_path,
            filename = "Samples_per_abslat_hemi_hist",
            file_ext = plot_ext,
            height = 8.5,
            width = 14,
            dpi = 300,
            device = plot_device
          )

          cat("Latitude plot created successfully.\n")
        },
        error = function(e) {
          cat("Error in creating Latitude plot:", conditionMessage(e), "\n")
        }
      )

      # Season plot
      tryCatch(
        {
          samples_per_season <- count_and_arrange(query_df, c("season", "hemi"))

          samples_per_season_plot <- suppressMessages(bar_plot(
            df = samples_per_season,
            x = n,
            y = season,
            fill = hemi,
            y_lab = "Season",
            x_lab = "Number of samples"
          ) +
            labs(fill = "Hemisphere") +
            scale_fill_manual(
              values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
              label = c("northern", "southern")
            ))

          save_custom_plot(
            plot = samples_per_season_plot,
            output_dir = files_out_path,
            filename = "Samples_per_season_hemi_bar",
            file_ext = plot_ext,
            height = 8.5,
            width = 14,
            dpi = 300,
            device = plot_device
          )

          cat("Season plot created successfully.\n")
        },
        error = function(e) {
          cat("Error in creating Season plot:", conditionMessage(e), "\n")
        }
      )

      # Depth plots
      tryCatch(
        {
          samples_per_depth <- count_and_arrange(query_df, c("depth", "hemi")) %>%
            mutate(less_than_150 = depth <= 155)

          # Depths above 150m
          depths_above_150m_bar_plt <- suppressMessages(bar_plot(
            df = samples_per_depth %>% filter(less_than_150 == TRUE),
            y = n,
            x = depth,
            fill = hemi,
            x_lab = "Depth (m)",
            y_lab = "Number of samples",
            width = 1,
            legend_position = "none"
          ) +
            scale_fill_manual(
              values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
              label = c("northern", "southern")
            ) +
            scale_y_continuous(expand = c(0, 0)) +
            scale_x_reverse() +
            coord_flip())

          save_custom_plot(
            plot = depths_above_150m_bar_plt,
            output_dir = files_out_path,
            filename = "Samples_per_depth_hemi_above150m_hist",
            file_ext = plot_ext,
            height = 8.5,
            width = 14,
            dpi = 300,
            device = plot_device
          )

          # Depths below 150m
          depths_below_150m_bar_plt <- suppressMessages(bar_plot(
            df = samples_per_depth %>% filter(less_than_150 == FALSE),
            y = n,
            x = depth,
            fill = hemi,
            x_lab = "Depth (m)",
            width = 30,
            legend_position = "right",
            legend_direction = "vertical"
          ) +
            scale_fill_manual(
              values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
              label = c("northern", "southern")
            ) +
            scale_x_reverse(
              limits = c(3015, 150),
              breaks = c(150, 300, 500, 750, 1000, 2000, 3000)
            ) +
            scale_y_continuous(
              limits = c(0, 10),
              expand = c(0, 0)
            ) +
            coord_flip())

          save_custom_plot(
            plot = depths_below_150m_bar_plt,
            output_dir = files_out_path,
            filename = "Samples_per_depth_hemi_below150_hist",
            file_ext = plot_ext,
            height = 8.5,
            width = 14,
            dpi = 300,
            device = plot_device
          )

          # Combined depth plot
          combined_plot_depths <- depths_above_150m_bar_plt | depths_below_150m_bar_plt

          save_custom_plot(
            plot = combined_plot_depths,
            output_dir = files_out_path,
            filename = "Samples_per_depth_hemi_hist_cmb_plot",
            file_ext = plot_ext,
            height = 8.5,
            width = 14,
            dpi = 300,
            device = plot_device
          )

          cat("Depth plots created successfully.\n")

          return(list(
            samples_per_lat_abs_plot = samples_per_lat_abs_plot,
            samples_per_ocean_plot = samples_per_ocean_plot
          ))
        },
        error = function(e) {
          cat("Error in creating Depth plots:", conditionMessage(e), "\n")
        }
      )


      cat("All plots generated successfully.\n")
    },
    error = function(e) {
      cat("Error in generate_sample_distribution_plots:", conditionMessage(e), "\n")
      stop("Plot generation failed. Please check the error messages above.")
    }
  )
}


#' Generate Oceanographic Data Plots
#'
#' This function creates various plots to visualize oceanographic data
#' distributions across different parameters such as SST, PO4, Fe, PP, and CHL.
#'
#' @param query_df A data frame containing the query results
#' @param files_out_path The output directory for saving plots
#' @param plot_ext The file extension for saved plots
#' @param plot_device The device to use for plotting
#'
#' @return A list containing samples_per_ocean and samples_per_studyid_ocean data frames
#'
#' @export
generate_oceanographic_plots <- function(query_df, files_out_path, plot_ext, plot_device) {
  tryCatch(
    {
      # SST plot
      tryCatch(
        {
          samples_per_sst_tblSST_AVHRR_OI_NRT <- count_and_arrange(query_df, c("sst_tblSST_AVHRR_OI_NRT", "hemi"))

          sst_hist <- suppressMessages(histogram_plot_x_or_y(
            df = samples_per_sst_tblSST_AVHRR_OI_NRT, aes_var = "x",
            x = sst_tblSST_AVHRR_OI_NRT,
            fill = hemi,
            binwidth = 1,
            x_lab = "SST ˚C"
          ) +
            theme(legend.position = "bottom") +
            labs(fill = "Hemisphere") +
            scale_fill_manual(
              values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
              label = c("northern", "southern")
            ))

          save_custom_plot(
            plot = sst_hist,
            output_dir = files_out_path,
            filename = "Samples_per_SST_hemi_hist",
            file_ext = plot_ext,
            height = 8.5,
            width = 14,
            dpi = 300,
            device = plot_device
          )

          cat("SST plot created successfully.\n")
        },
        error = function(e) {
          cat("Error call in:", deparse(conditionCall(e)), "\n")
          stop(cat("SST plot was not created due to:\n", conditionMessage(e), "\n"))
        }
      )

      # PO4 plot
      tryCatch(
        {
          samples_per_PO4_tblPisces_NRT <- count_and_arrange(query_df, c("PO4_tblPisces_NRT", "hemi"))

          po4_hist <- suppressMessages(histogram_plot_x_or_y(
            df = samples_per_PO4_tblPisces_NRT, aes_var = "x",
            x = PO4_tblPisces_NRT,
            fill = hemi,
            binwidth = 0.02,
            x_lab = expression(bold(paste("PO"[4]^"3-" * " (µmol L"^"-1", ")")))
          ) +
            theme(legend.position = "bottom") +
            labs(fill = "Hemisphere") +
            scale_fill_manual(
              values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
              label = c("northern", "southern")
            ) +
            xlim(0, 2.0))

          save_custom_plot(
            plot = po4_hist,
            output_dir = files_out_path,
            filename = "Samples_per_PO4_hemi_hist",
            file_ext = plot_ext,
            height = 8.5,
            width = 14,
            dpi = 300,
            device = plot_device
          )

          cat("PO4 plot created successfully.\n")
        },
        error = function(e) {
          cat("Error call in:", deparse(conditionCall(e)), "\n")
          stop(cat("PO4 plot was not created due to:\n", conditionMessage(e), "\n"))
        }
      )

      # logFe plot
      tryCatch(
        {
          samples_per_logFe <- count_and_arrange(query_df, c("logFe", "hemi"))

          samples_per_logFe_plot <- histogram_plot_x_or_y(
            df = samples_per_logFe, aes_var = "x",
            x = logFe,
            fill = hemi,
            binwidth = 0.1,
            x_lab = expression(bold(paste("log(Fe)")))
          ) +
            theme(legend.position = "bottom") +
            labs(fill = "Hemisphere") +
            scale_fill_manual(
              values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
              label = c("northern", "southern")
            )

          save_custom_plot(
            plot = samples_per_logFe_plot,
            output_dir = files_out_path,
            filename = "Samples_per_logFe_hemi_hist",
            file_ext = plot_ext,
            height = 8.5,
            width = 14,
            dpi = 300,
            device = plot_device
          )

          cat("logFe plot created successfully.\n")
        },
        error = function(e) {
          cat("Error call in:", deparse(conditionCall(e)), "\n")
          stop(cat("logFe plot was not created due to:\n", conditionMessage(e), "\n"))
        }
      )

      # PP plot
      tryCatch(
        {
          samples_per_PP_tblPisces_NRT <- count_and_arrange(query_df, c("PP_tblPisces_NRT", "hemi"))

          samples_per_PP_tblPisces_NRT_plot <- histogram_plot_x_or_y(
            df = samples_per_PP_tblPisces_NRT, aes_var = "x",
            x = PP_tblPisces_NRT,
            fill = hemi,
            binwidth = 0.001,
            x_lab = expression(bold(paste("PP")))
          ) +
            theme(legend.position = "bottom") +
            labs(fill = "Hemisphere") +
            scale_fill_manual(
              values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
              label = c("northern", "southern")
            )

          save_custom_plot(
            plot = samples_per_PP_tblPisces_NRT_plot,
            output_dir = files_out_path,
            filename = "Samples_per_PP_tblPisces_NRT_hemi_hist",
            file_ext = plot_ext,
            height = 8.5,
            width = 14,
            dpi = 300,
            device = plot_device
          )

          cat("PP plot created successfully.\n")
        },
        error = function(e) {
          cat("Error call in:", deparse(conditionCall(e)), "\n")
          stop(cat("PP plot was not created due to:\n", conditionMessage(e), "\n"))
        }
      )

      # CHL plot
      tryCatch(
        {
          samples_per_CHL_tblPisces_NRT <- count_and_arrange(query_df, c("CHL_tblPisces_NRT", "hemi"))

          samples_per_CHL_tblPisces_NRT_plot <- histogram_plot_x_or_y(
            df = samples_per_CHL_tblPisces_NRT, aes_var = "x",
            x = CHL_tblPisces_NRT,
            fill = hemi,
            binwidth = 0.1,
            x_lab = expression(bold(paste("chlorophyll a")))
          ) +
            theme(legend.position = "bottom") +
            labs(fill = "Hemisphere") +
            scale_fill_manual(
              values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
              label = c("northern", "southern")
            )

          save_custom_plot(
            plot = samples_per_CHL_tblPisces_NRT_plot,
            output_dir = files_out_path,
            filename = "Samples_per_CHL_tblPisces_NRT_hemi_hist",
            file_ext = plot_ext,
            height = 8.5,
            width = 14,
            dpi = 300,
            device = plot_device
          )

          cat("CHL plot created successfully.\n")
        },
        error = function(e) {
          cat("Error call in:", deparse(conditionCall(e)), "\n")
          stop(cat("CHL plot was not created due to:\n", conditionMessage(e), "\n"))
        }
      )

      # Calculate and return additional data
      samples_per_ocean <- count_and_arrange(query_df, c("ocean", "hemi"))
      samples_per_studyid_ocean <- count_and_arrange(query_df, c("studyID", "ocean"))

      cat("All oceanographic plots generated successfully.\n")
      return(list(
        sst_hist = sst_hist,
        po4_hist = po4_hist,
        samples_per_ocean = samples_per_ocean,
        samples_per_studyid_ocean = samples_per_studyid_ocean
      ))
    },
    error = function(e) {
      cat("Error in generate_oceanographic_plots:", conditionMessage(e), "\n")
      stop("Oceanographic plot generation failed. Please check the error messages above.")
    }
  )
}


#' Create Combined Plot
#'
#' This function creates a combined plot from four individual plots: latitude, ocean, SST, and PO4.
#'
#' @param samples_per_lat_abs_plot Latitude plot
#' @param samples_per_ocean_plot Ocean plot
#' @param sst_hist SST histogram
#' @param po4_hist PO4 histogram
#' @param files_out_path Output directory for saving plots
#' @param plot_ext File extension for saved plots
#' @param plot_device Device to use for plotting
#'
#' @return NULL
#'
#' @export
create_combined_plot <- function(samples_per_lat_abs_plot, samples_per_ocean_plot, sst_hist, po4_hist, files_out_path, plot_ext, plot_device) {
  tryCatch(
    {
      # Alter plots for combining
      tryCatch(
        {
          sst_hist_sub <- sst_hist + theme(
            legend.position = "none",
            plot.tag = element_text(face = "bold", size = 25)
          ) +
            labs(tag = "(c)")

          po4_hist_sub <- po4_hist + theme(
            legend.position = "none",
            plot.tag = element_text(face = "bold", size = 25)
          ) +
            labs(tag = "(d)")

          samples_per_lat_abs_plot_sub <- samples_per_lat_abs_plot + theme(
            legend.position = "none",
            plot.tag = element_text(face = "bold", size = 25)
          ) +
            labs(tag = "(a)")

          samples_per_ocean_plot_sub <- samples_per_ocean_plot + theme(
            legend.position = "none",
            axis.text.y = element_text(
              angle = 45,
              hjust = 1,
              size = 19,
              face = "bold"
            ),
            plot.tag = element_text(face = "bold", size = 25)
          ) +
            labs(tag = "(b)")

          cat("Individual plots modified successfully.\n")
        },
        error = function(e) {
          cat("Error call in:", deparse(conditionCall(e)), "\n")
          stop(cat("Individual plots modification failed due to:\n", conditionMessage(e), "\n"))
        }
      )

      # Combine plots
      tryCatch(
        {
          combined_plot_fig_5 <- (
            (samples_per_lat_abs_plot_sub + samples_per_ocean_plot_sub) /
              (sst_hist_sub + po4_hist_sub)
          )

          print(combined_plot_fig_5)

          cat("Combined plot created successfully.\n")
        },
        error = function(e) {
          cat("Error call in:", deparse(conditionCall(e)), "\n")
          stop(cat("Combined plot creation failed due to:\n", conditionMessage(e), "\n"))
        }
      )

      # Save combined plot
      tryCatch(
        {
          save_custom_plot(
            plot = combined_plot_fig_5,
            output_dir = files_out_path,
            filename = "f05",
            file_ext = plot_ext,
            height = 8.5,
            width = 17,
            dpi = 300,
            device = plot_device
          )

          cat("Combined plot saved successfully.\n")
        },
        error = function(e) {
          cat("Error call in:", deparse(conditionCall(e)), "\n")
          stop(cat(
            "Combined plot saving failed due to:\n",
            conditionMessage(e), "\n"
          ))
        }
      )

      cat("All operations completed successfully.\n")
    },
    error = function(e) {
      cat("Error in create_combined_plot:", conditionMessage(e), "\n")
      stop("Combined plot creation failed. Please check the error messages above.")
    }
  )
}


#' Main function to execute analysis
main <- function(
    files_out_path,
    plot_ext,
    plot_device,
    files_to_read, files_in_path, tables_out_path) {
  # Load the data
  data_list <- load_files(files_to_read, files_in_path)

  # Generate query dataframe by deduplicating by group
  query_df <- dedup_by_group(
    data_list$cmap_coloc,
    group_id_key = data_list$unique_sample_id_key,
    group_id
  )

  samp_plot_list <- generate_sample_distribution_plots(
    query_df,
    data_list$samples_per_nucleicAcidType,
    data_list$samples_per_nucleicAcidType_studyid,
    data_list$samples_per_photic,
    data_list$samples_per_photic_nucacid,
    data_list$sample_type,
    files_out_path,
    plot_ext,
    plot_device
  )

  ocean_plot_list <- generate_oceanographic_plots(
    query_df,
    files_out_path,
    plot_ext,
    plot_device
  )

  create_combined_plot(
    samp_plot_list$samples_per_lat_abs_plot,
    samp_plot_list$samples_per_ocean_plot,
    ocean_plot_list$sst_hist,
    ocean_plot_list$po4_hist,
    files_out_path,
    plot_ext,
    plot_device
  )

  # Write out plots
  if (!is.null(samp_plot_list) && !is.null(ocean_plot_list)) {
    create_dir(files_out_path)
    write_file_list(
      file_list = samp_plot_list,
      path = tables_out_path
    )
  }

  message("\n\nPlots generated successfully.\n")
}


# Run if the script is being run directly
if (sys.nframe() == 0 && !interactive()) {
  source_file(files_to_source)

  parser <- setup_parser()
  args <- parse_arg(parser)

  validate_parsed_args(parsed_args = args)

  # Create output directory if it doesn't exist
  create_dir(args$files_out_path)

  main(
    files_to_read = args$files_to_read,
    files_in_path = args$files_in_path,
    files_out_path = args$files_out_path,
    plot_ext = args$plot_ext,
    plot_device = args$plot_device,
    tables_out_path = args$tables_out_path
  )
}
