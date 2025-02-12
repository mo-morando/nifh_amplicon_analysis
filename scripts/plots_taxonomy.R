#!/usr/bin/env Rscript

#' @title NifH Amplicon Data Analysis Taxonomic Plots
#' @description This script processes NifH amplicon sequencing data, performs statistical analyses, and generates comprehensive visualizations to explore nifH cluster distributions across different studies and environmental parameters.
#'
#' @details The pipeline executes the following main steps:
#' * Loads and validates input files containing count data, relative abundance, and metadata
#' * Generates a colorblind-safe palette for nifH cluster visualization
#' * Creates pie charts for overall nifH cluster distribution (counts and relative abundance)
#' * Produces bar plots showing nifH cluster percentages across different studies
#' * Generates scatter plots of relative abundance vs. absolute latitude and sea surface temperature
#' * Saves all visualizations in specified formats
#'
#' Key functions include:
#' * pie_charts(): Creates pie charts for nifH cluster percentages
#' * bar_plots(): Generates bar plots for nifH cluster percentages by study
#' * make_scatter(): Produces scatter plots for relative abundance data
#' * main(): Orchestrates the entire analysis and visualization workflow
#'
#' @usage Rscript plots_taxonomy.R [--files FILES] [--input_path PATH] [--output_path PATH] [--plot_ext EXT] [--plot_device DEVICE]
#'
#' @param --files Comma-separated list of input files
#' @param --input_path Directory path for input files
#' @param --output_path Directory path for output plots
#' @param --plot_ext File extension for saved plots (e.g., ".jpeg", ".png")
#' @param --plot_device Device to use for plot generation (e.g., "jpeg", "png", "cairo_pdf")
#'
#' @author Michael (Mo) Morando
#'  2025-02-02
#'
#' @note This script requires the following R packages: tidyverse, patchwork, argparser
#'
#' @examples
#' Rscript plots_taxonomy.R --files nifhdb_all_counts_AUID_dedup_clean,nifhdb_all_rel_abund_AUID_dedup_clean,nifhdb_all_counts_AUID_dedup_total_study_id_clean,nifhdb_all_rel_abund_AUID_dedup_total_study_id_clean,nifhdb_all_rel_abund_AUID_dedup_study_id_total_clean,nifhdb_all_counts_AUID_dedup_study_id_total_clean,nifhdb_all_counts_AUID_dedup_study_id_total,nifhdb_all_rel_abund_AUID_dedup_study_id_total,RA_df_T_lng_mean_RA_AUID_deduped,cmap_coloc,annoNifHDB_updt --input_path ../analysis/out_files --output_path ../analysis/plots --plot_ext .jpeg --plot_device jpeg
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
# files_to_source <- c(
#   "functions.R",
#   "basic_plotting.R"
# )


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
  parser <- arg_parser("Process input files and generate plots for nifH database analysis")
  parser <- add_argument(parser, "--files",
    help = "CVS list of files to read in",
    default = "nifhdb_all_counts_AUID_dedup_clean,nifhdb_all_rel_abund_AUID_dedup_clean,nifhdb_all_counts_AUID_dedup_total_study_id_clean,nifhdb_all_rel_abund_AUID_dedup_total_study_id_clean,nifhdb_all_rel_abund_AUID_dedup_study_id_total_clean,nifhdb_all_counts_AUID_dedup_study_id_total_clean,nifhdb_all_counts_AUID_dedup_study_id_total,nifhdb_all_rel_abund_AUID_dedup_study_id_total,RA_df_T_lng_mean_RA_AUID_deduped,cmap_coloc,annoNifHDB_updt"
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
}



#' Create a Pie Chart for nifH Cluster Percentages
#'
#' This function creates a pie chart for nifH cluster percentages.
#'
#' @param data Data frame containing the plot data.
#' @param colors Named vector of colors for the nifH clusters.
#' @return A ggplot2 object representing the pie chart.
#'
create_pie_chart <- function(data, colors) {
  tryCatch({
    ggplot(data, aes(x = 1, y = percentage_nifH_cluster, fill = nifH_cluster_modified)) +
      geom_bar(stat = "identity") +
      coord_polar(theta = "y") +
      theme_void() +
      labs(fill = "nifH Cluster") +
      geom_text(aes(label = ifelse(percentage_nifH_cluster >= 5, paste0(percentage_nifH_cluster, "%"), ""), x = 1.15),
        position = position_stack(vjust = 0.5), hjust = 0.5, size = 12, show.legend = FALSE
      ) +
      geom_text(aes(label = ifelse(percentage_nifH_cluster < 5, paste0(percentage_nifH_cluster, "%"), ""), x = 1.5),
        position = position_stack(vjust = 0.5), hjust = 0.5, size = 5, show.legend = FALSE
      ) +
      scale_fill_manual(values = colors) +
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 17, face = "bold")
      ) +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "top", title.hjust = 0.5))
  },
  error = function(e) {
    cat("Error in:", deparse(conditionCall(e)), "\n")
    stop("Error creating pie charts due to:\n", conditionMessage(e))
  }
  )
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
  cat("Generating nifH palette..")
  tryCatch({
    # Generate palette
    pal <- c(RColorBrewer::brewer.pal(12, "Paired")[c(8, 2, 11, 9, 5, 4, 10, 3)], "#777777")

    # Assign colors to nifH cluster names
    names(pal) <- names(nifh_cluster_colours)
    return(pal)
  },
  error = function(e) {
    cat("Error in:", deparse(conditionCall(e)), "\n")
    stop("Error generating nifH palette due to:\n", conditionMessage(e))
  }
  )
}



#' Generate Pie Charts for nifH Cluster Percentages
#'
#' This function creates pie charts for nifH cluster percentages based on counts and relative abundance.
#' It saves the generated plots as image files and returns them as a list.
#'
#' @param counts_data Data frame containing counts data for nifH clusters.
#' @param rel_abund_data Data frame containing relative abundance data for nifH clusters.
#' @param output_dir Directory where the plots will be saved.
#' @param colors Named vector of colors for the nifH clusters.
#' @param plot_ext File extension for saving plots (e.g., ".jpeg", ".png").
#' @param plot_device Plot device used for writing plots (e.g., "cairo_pdf", "svg").
#' @return A list containing three ggplot objects:
#'   \itemize{
#'     \item counts: Pie chart of counts data
#'     \item rel_abund: Pie chart of relative abundance data
#'     \item combined: Combined plot of both pie charts
#'   }
#'   Returns NULL if any error occurs during plot creation.
#' @examples
#' plots <- pie_charts(
#'   nifhdb_all_counts_AUID_dedup_clean,
#'   nifhdb_all_rel_abund_AUID_dedup_clean,
#'   "/path/to/output",
#'   nifh_cluster_colours_colbldsafe,
#'   ".jpeg"
#' )
pie_charts <- function(counts_data, rel_abund_data, output_dir, colors, plot_ext = ".jpeg", plot_device = "jpeg") {
  # Check if the output directory exists
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist.")
  }

  # Create and save pie chart for counts
  cat("Creating pie chart for counts data...\n")
  tryCatch(
    {
      pie_chart_counts <- create_pie_chart(counts_data, colors)

      save_custom_plot(
        plot = pie_chart_counts,
        output_dir = output_dir,
        filename = "DNA_perc_cnts_clst_pie",
        file_ext = plot_ext,
        height = 8.5,
        width = 14,
        dpi = 300,
        device = plot_device
      )
    },
    # warning = function(w) {
    #   cat("Warning occurred:\n")
    #   cat("Warning message:", conditionMessage(w), "\n")
    #   cat("Warning call:", deparse(conditionCall(w)), "\n")
    # },
    error = function(e) {
      cat("An error occurred in pie_charts() while creating the counts pie chart:\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Error call:", deparse(conditionCall(e)), "\n")
    }
  )


  # Create and save pie chart for relative abundance
  cat("Creating pie chart for relative abundance data...\n")
  tryCatch(
    {
      pie_chart_rel_abund <- create_pie_chart(rel_abund_data, colors)

      save_custom_plot(
        plot = pie_chart_rel_abund,
        output_dir = output_dir,
        filename = "DNA_perc_rel_abund_clst_pie",
        file_ext = plot_ext,
        height = 8.5,
        width = 14,
        dpi = 300,
        device = plot_device
      )
    },
    # warning = function(w) {
    #   cat("Warning occurred:\n")
    #   cat("Warning message:", conditionMessage(w), "\n")
    #   cat("Warning call:", deparse(conditionCall(w)), "\n")
    # },
    error = function(e) {
      cat("An error occurred in pie_charts() while creating the relative abundance pie chart:\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Error call:", deparse(conditionCall(e)), "\n")
    }
  )


  # Combine and save combined plot
  cat("Combining pie charts using patchwork...\n")
  tryCatch(
    {
      if (
        all(
          map_lgl(c("pie_chart_counts", "pie_chart_rel_abund"), ~exists(.x))
        ) &&
          all(!map_lgl(c(pie_chart_counts, pie_chart_rel_abund), is.null))
      ) {
        combined_plot <- (pie_chart_counts | pie_chart_rel_abund)

        save_custom_plot(
          plot = combined_plot,
          output_dir = output_dir,
          filename = "DNA_perc_combined_clst_pie",
          file_ext = plot_ext,
          height = 8.5,
          width = 16,
          dpi = 300,
          device = plot_device
        )

      } else {
        stop("One of the pie charts is missing or NULL. Combined plot not created.")
      }
    },
    # warning = function(w) {
    #   cat("Warning occurred:\n")
    #   cat("Warning message:", conditionMessage(w), "\n")
    #   cat("Warning call:", deparse(conditionCall(w)), "\n")
    # },
    error = function(e) {
      cat("An error occurred in pie_charts() while combining the pie charts:\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Error call:", deparse(conditionCall(e)), "\n")
      return(NULL)
    }
  )


  # Define plot list with actual objects
  pie_chart_results <- list(
    counts = if (exists("pie_chart_counts")) pie_chart_counts else NULL,
    rel_abund = if (exists("pie_chart_rel_abund"))
      pie_chart_rel_abund else NULL,
    combined = if (exists("combined_plot")) combined_plot else NULL
  )

  # Check existence of all plots
  cat(ifelse(
    all(!map_lgl(pie_chart_results, is.null)),
    "All pie chart plots exist and will be returned\n\n",
    "At least one pie chart plot does not exist and will not be saved or returned.
    Please see above error(s).\n\n"
  ))


  return(pie_chart_results)
}


#' Generate Bar Plots for nifH Cluster Percentages
#'
#' This function creates bar plots for nifH cluster percentages based on counts and relative abundance.
#' It saves the generated plots as image files with the specified extension.
#'
#' @param counts_data_st_id_clean Cleaned counts data with study ID.
#' @param rel_abund_data_st_id_clean Cleaned relative abundance data with study ID.
#' @param counts_data Raw counts data.
#' @param counts_data_st_id Counts data with study ID.
#' @param rel_abund_data Raw relative abundance data.
#' @param rel_abund_data_st_id Relative abundance data with study ID.
#' @param output_dir Directory where the plots will be saved.
#' @param colors Named vector of colors for the nifH clusters.
#' @param plot_ext File extension for saving plots.
#' @param plot_device Plot device used for writing plots.
#' @return Logical indicating success or failure.
#' @examples
#' bar_plots(nifhdb_all_counts_AUID_dedup_study_id_total, nifhdb_all_rel_abund_AUID_dedup_study_id_total, "/path/to/output", nifh_cluster_colours_colbldsafe, ".jpeg")
bar_plots <- function(
    counts_data_st_id_clean,
    rel_abund_data_st_id_clean,
    counts_data,
    counts_data_st_id,
    rel_abund_data,
    rel_abund_data_st_id,
    output_dir,
    colors,
    plot_ext = ".jpeg",
    plot_device = "jpeg") {

  # Check if the output directory exists
  cat("Checking if the output directory exists...\n")
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist.")
  }


  # Create and save bar plot for counts
  cat("Creating bar plot for counts data...\n")
  tryCatch(
    {
      bar_plot_counts <- bar_plot(
        df = counts_data_st_id_clean,
        x = studyID,
        y = percentage_total_counts_nifH_cluster_study_id_total_clean,
        fill = nifH_cluster_modified,
        fill_pallete = colors,
        x_lab = "Study ID",
        y_lab = bquote(bold("%" ~ total ~ reads)),
        legend_position = "bottom",
        x_axis_angle = TRUE,
        print_out = TRUE
      )

      save_custom_plot(
        plot = bar_plot_counts,
        output_dir = output_dir,
        filename = "DNA_perc_cnts_clst_bar",
        file_ext = plot_ext,
        height = 8.5,
        width = 14,
        dpi = 300,
        device = plot_device
      )
    },
    # warning = function(w) {
    #   cat("Warning occurred:\n")
    #   cat("Warning message:", conditionMessage(w), "\n")
    #   cat("Warning call:", deparse(conditionCall(w)), "\n")
    # },
    error = function(e) {
      cat("An error occurred in bar_plots() while creating the counts bar plot:\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Error call:", deparse(conditionCall(e)), "\n")
    }
  )


  # Create and save bar plot for relative abundance
  cat("Creating bar plot for relative abundance data...\n")
  tryCatch(
    {
      bar_plot_rel_abund <- bar_plot(
        df = rel_abund_data_st_id_clean,
        x = studyID,
        y = percentage_total_rel_abund_nifH_cluster_study_id_total_clean,
        fill = nifH_cluster_modified,
        fill_pallete = colors,
        x_lab = "Study ID",
        y_lab = bquote(bold("%" ~ relative ~ abundance)),
        legend_position = "bottom",
        x_axis_angle = TRUE,
        print_out = TRUE
      )

      save_custom_plot(
        plot = bar_plot_rel_abund,
        output_dir = output_dir,
        filename = "DNA_perc_rel_abund_clst_bar",
        file_ext = plot_ext,
        height = 8.5,
        width = 14,
        dpi = 300,
        device = plot_device
      )
    },
    # warning = function(w) {
    #   cat("Warning occurred:\n")
    #   cat("Warning message:", conditionMessage(w), "\n")
    #   cat("Warning call:", deparse(conditionCall(w)), "\n")
    # },
    error = function(e) {
      cat("An error occurred in bar_plots() while creating the relative abundance bar plot:\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Error call:", deparse(conditionCall(e)), "\n")
    }
  )


  # Create a plot with all the pooled data as a column for better comparison
  # First, counts data
  custom_order_w_total <- c(
    "pooled data",
    "AK2HI",
    "BentzonTilia_2015",
    "Ding_2021",
    "Gradoville_2020_G1",
    "Gradoville_2020_G2",
    "Hallstrom_2021",
    "Hallstrom_2022",
    "Harding_2018",
    "Mulholland_2018",
    "NEMO",
    "Raes_2020",
    "Sato_2021",
    "Selden_2021",
    "Shiozaki_2017",
    "Shiozaki_2018GBC",
    "Shiozaki_2018LNO",
    "Shiozaki_2020",
    "Tang_2020",
    "Wu_2019",
    "Wu_2021",
    "TurkKubo_2021"
  )

  cat("Combining pooled data with study ID data for counts...\n")
  tryCatch(
    {
      combined_data_counts <- counts_data_st_id %>%
        bind_rows(counts_data %>%
          rename(
            # nifH_cluster = nifH_cluster_modified,
            percentage_total_counts_nifH_cluster_study_id_total = percentage_nifH_cluster
          ) %>%
          mutate(
            studyID = "pooled data", #* # make new studyId for pooled data
            cluster_stats = nifH_cluster_modified #* # this is needed for plotting
          )) %>%
        mutate(
          studyID = factor(studyID, levels = custom_order_w_total) #* # make new plotting order
        )
      cat("Successfully combined data...\n")
    },
    # warning = function(w) {
    #   cat("Warning occurred while combining study data:\n")
    #   cat("Warning message:", conditionMessage(w), "\n")
    #   cat("Warning call:", deparse(conditionCall(w)), "\n")
      # return(combined_data_counts)
    # },
    error = function(e) {
      cat("An error occurred while combining study data:\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Error call:", deparse(conditionCall(e)), "\n")
      return(NULL)
    }
  )


  # Create bar plot using the combined data
  cat("Creating combined bar plot of pooled and study ID data...\n")
  tryCatch({
    if (exists("combined_data_counts") && !is.null(combined_data_counts)) {
      bar_plot_counts_combined <- bar_plot(
        df = combined_data_counts,
        x = studyID,
        y = percentage_total_counts_nifH_cluster_study_id_total,
        fill = cluster_stats,
        fill_pallete = colors,
        fill_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
        x_lab = "Study ID",
        y_lab = bquote(bold("%" ~ total ~ reads)),
        legend_position = "bottom",
        x_axis_angle = TRUE,
        print_out = TRUE
      )

      save_custom_plot(
        plot = bar_plot_counts_combined,
        output_dir = output_dir,
        filename = "DNA_perc_cnts_clst_all_bar",
        file_ext = plot_ext,
        height = 8.5,
        width = 14,
        dpi = 300,
        device = plot_device
      )
    } else {
      stop("There is a problem wih combined_data_counts.\n")
    }
  },
      # warning = function(w) {
      #   cat("Warning occurred while creating bar plot:\n")
      #   cat("Warning message:", conditionMessage(w), "\n")
      #   cat("Warning call:", deparse(conditionCall(w)), "\n")
      # },
      error = function(e) {
        cat("An error occurred while creating bar plot:\n")
        cat("Error message:", conditionMessage(e), "\n")
        cat("Error call:", deparse(conditionCall(e)), "\n")
        return(NULL)
      }
  )


  # Next, relative abundance data
  cat("Combining pooled data with study ID data...\n")
  tryCatch(
    {
      combined_data_rel_abund <- rel_abund_data_st_id %>%
        bind_rows(rel_abund_data %>%
          rename(
            percentage_total_rel_abund_nifH_cluster_study_id_total = percentage_nifH_cluster
          ) %>%
          mutate(
            studyID = "pooled data", #* # make new studyId for pooled data
            cluster_stats = nifH_cluster_modified #* # this is needed for plotting
          )) %>%
        mutate(
          studyID = factor(studyID, levels = custom_order_w_total) #* # make new plotting order
        )
      cat("Successfully combined data...\n")
    },
    # warning = function(w) {
    #   cat("Warning occurred while combining study data:\n")
    #   cat("Warning message:", conditionMessage(w), "\n")
    #   cat("Warning call:", deparse(conditionCall(w)), "\n")
      # return(combined_data_rel_abund)
    # },
    error = function(e) {
      cat("An error occurred while combining study data:\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Error call:", deparse(conditionCall(e)), "\n")
      return(NULL)
    }
  )

  # Create bar plot using the combined data
  cat("Creating combined bar plot of pooled and study ID data...\n")
  tryCatch({
    if (
      exists("combined_data_rel_abund") && !is.null(combined_data_rel_abund)) {
        bar_plot_rel_abund_combined <- bar_plot(
          df = combined_data_rel_abund,
          x = studyID,
          y = percentage_total_rel_abund_nifH_cluster_study_id_total,
          fill = cluster_stats,
          fill_pallete = colors,
          fill_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
          x_lab = "Study ID",
          y_lab = bquote(bold("%" ~ relative ~ abundance)),
          legend_position = "bottom",
          x_axis_angle = TRUE,
          print_out = TRUE
        )

        save_custom_plot(
        plot = bar_plot_rel_abund_combined,
        output_dir = output_dir,
        # filename = "DNA_perc_rel_abund_clst_all_bar",
        filename = "f07",
        file_ext = plot_ext,
        height = 8.5,
        width = 14,
        dpi = 300,
        device = plot_device
      )
    } else {
      stop("There was a problem with combined_data_rel_abund.")
    }
      },
      # warning = function(w) {
      #   cat("Warning occurred while creating bar plot:\n")
      #   cat("Warning message:", conditionMessage(w), "\n")
      #   cat("Warning call:", deparse(conditionCall(w)), "\n")
      # },
      error = function(e) {
        cat("An error occurred while creating combined bar plot of pooled and study ID data:\n")
        cat("Error message:", conditionMessage(e), "\n")
        cat("Error call:", deparse(conditionCall(e)), "\n")
        return(NULL)
      }
  )


  # Define plot list with actual objects
  bar_plot_results <- list(
    bar_counts = if (exists("bar_plot_counts")) bar_plot_counts else NULL,
    bar_rel_abund = if (exists("bar_plot_rel_abund")) bar_plot_rel_abund else NULL,
    bar_counts_combined = if (exists("bar_plot_counts_combined")) bar_plot_counts_combined else NULL,
    bar_rel_abund_combined = if (exists("bar_plot_rel_abund_combined")) bar_plot_rel_abund_combined else NULL
  )

  # Check existence of all plots
  cat(
    ifelse(
      all(!map_lgl(bar_plot_results, is.null)),
      "All bar plots exist and will be returned\n\n",
      "At least one bar plot does not exist and will not be saved or returned.\nPlease see above error(s).\n\n"
    )
  )


  return(bar_plot_results)
}



#' Create Scatter Plots for Relative Abundance Data
#'
#' This function generates scatter plots for relative abundance data against absolute latitude 
#' and sea surface temperature (SST).
#'
#' @param rel_abund_df Data frame containing relative abundance data.
#' @param anno_table Data frame containing annotation information.
#' @param meta_table Data frame containing metadata information.
#' @param colors Named vector of colors for plotting.
#' @param plot_ext File extension for saving plots (e.g., ".jpeg").
#' @param plot_device Plot device used for writing plots (e.g., "cairo_pdf", "svg").
#' @param output_dir Directory where the plots will be saved.
#' @return A list containing:
#'   \itemize{
#'     \item data: Processed data frame used for plotting
#'     \item plots: List of generated plots (latitude, sst, combined)
#'   }
#'   Returns NULL if any error occurs during execution.
make_scatter <- function(
    rel_abund_df,
    anno_table,
    meta_table,
    colors,
    plot_ext,
    plot_device,
    output_dir) {

  # Check if the output directory exists
  cat("Checking if the output directory exists...\n")
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist")
  }


  cat("Merging relative abundance data with annotations and metadata...\n")
  tryCatch({
    query_df <- rel_abund_df %>%
      merge_annotations(annotation_table = anno_table) %>%
      merge_cmap(cmap = meta_table, by_join = c("SAMPLEID"))
      cat("Successfully merged files and generated query_df\n")

    # Sum each sample by nifH cluster and select relevant columns
    cat("Summing each sample by nifH cluster and selecting relevant columns...\n")
    query_df_cluster <- query_df %>%
      group_by(SAMPLEID, nifH_cluster) %>%
      mutate(RA = sum(RA, na.rm = T)) %>%
      select(SAMPLEID, studyID, time, depth, lat, lon, nifH_cluster, RA, CyanoCON, lat_abs, hemi, everything()) %>%
      ungroup() %>%
      distinct(SAMPLEID, nifH_cluster, .keep_all = T) %>%
      ungroup()

    # Filter for specific nifH clusters of interest
    cat("Filtering for specific nifH clusters of interest...\n")
    query_df_cluster %>%
      filter(nifH_cluster %in% c("1J/1K", "3", "1A", "1B", "1G", "1O/1P"))
    # cat(xxx)
  }, error = function(e) {
    cat("Error processing query_df:", conditionMessage(e), "\n")
  })


  cat("Creating scatter plots...\n")
  tryCatch({
    if (!exists("query_df_cluster") || is.null(query_df_cluster)) {
      stop("There is a problem with query_df_cluster\n")
    } else {
      cat("Creating scatter plot: Relative Abundance vs Absolute Latitude...\n")
      # Create scatter plot: Relative Abundance vs Absolute Latitude
      tryCatch({
        lat_abs <- scatter_line_plot(
        df = query_df_cluster,
        y = RA,
        x = lat_abs,
        colour = nifH_cluster,
        group = nifH_cluster,
        colour_pallete = colors,
        colour_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
        fill_lab = NULL,
        fill = nifH_cluster,
        fill_pallete = colors,
        title_lab = NULL,
        subtitle_lab = NULL,
        x_lab = "absolute latitude",
        y_lab = "% relative abundance",
        legend_position = "bottom",
        legend_direction = "horizontal"
      )
      }, error = function(e) {
        cat("Error creating latitude plot:", conditionMessage(e), "\n")
      })


    # Create scatter plot: Relative Abundance vs SST
    cat("Creating scatter plot: Relative Abundance vs SST...\n")
    ## RA versus SST
    tryCatch({
      sst <- scatter_line_plot(
      df = query_df_cluster,
      y = RA,
      x = sst_tblSST_AVHRR_OI_NRT,
      colour = nifH_cluster,
      group = nifH_cluster,
      colour_pallete = colors,
      colour_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
      fill_lab = NULL,
      fill = nifH_cluster,
      fill_pallete = colors,
      title_lab = NULL,
      subtitle_lab = NULL,
      x_lab = "SST ˚C",
      y_lab = "% relative abundance",
      legend_position = "bottom",
      legend_direction = "horizontal"
    ) +
      scale_x_reverse()
    }, error = function(e) {
      cat("Error creating SST plot:", conditionMessage(e), "\n")
    })

    # Combine plots using patchwork
    cat("Combining plots using patchwork...\n")
    tryCatch({
      if (
        all(map_lgl(c("lat_abs", "sst"), ~exists(.x))) &&
        all(!map_lgl(c(lat_abs, sst), is.null))
        ) {
      combined_sct_plot <- (lat_abs + theme(legend.position = "none")) / sst

      save_custom_plot(
          plot = combined_sct_plot,
          output_dir = output_dir,
          filename = "NifHsumCld_scat_bi_RAv_lat_SST_noLegend",
          file_ext = plot_ext,
          height = 12.5,
          width = 21,
          dpi = 300,
          device = plot_device
        )
      } else {
        stop("One of the bar plots required to genereate combined plot is missing or NULL")
      }
  },
  # warning = function(w) {
  #   cat("A warning occured \n")
  #   cat("Warning messege:", conditionMessage(w), "\n")
  #   cat("Warning call:", deparse(conditionCall(w)), "\n")
  # },
  error = function(e) {
    cat("An error occurred making combined scatter plots...\n")
    cat("Error message:", conditionMessage(e), "\n")
    cat("Error call:", deparse(conditionCall(e)), "\n")
  })
  }},
  error = function(e) {
    cat("An error occurred while making scatter plots...\n")
    cat("Error message:", conditionMessage(e), "\n")
    cat("Error call:", deparse(conditionCall(e)), "\n")
  })


  # Define plot mapping with actual objects
  scatter_plot_results <- list(
    data = if (exists("query_df_cluster")) query_df_cluster else NULL,
    plots = list(
      lat_abs = if (exists("lat_abs")) lat_abs else NULL,
      sst = if (exists("sst")) sst else NULL,
      combined = if (exists("combined_sct_plot")) combined_sct_plot else NULL
    )
  )

  # Check existence of all plot elements
  cat(
    ifelse(
      all(!map_lgl(c(scatter_plot_results$data, scatter_plot_results$plots), is.null)),
      "All scatter plots exist and will be returned\n\n",
      "At least one scatter plot does not exist and will not be saved or returned.\nPlease see above error(s).\n\n"
    )
  )


  # Return the plots and data
  return(scatter_plot_results)

}



#' Main execution function for the nifH Amplicon Data Analysis Script
#'
#' This function orchestrates the entire workflow of the script. It sets up
#' the argument parser, parses command-line arguments, sources necessary files,
#' reads in data files, and calls functions to generate various plots such as
#' pie charts, bar plots, and scatter plots based on the parsed input.
#'
#' The main steps performed by this function include:
#' 1. Setting up the argument parser to handle command-line inputs.
#' 2. Parsing the command-line arguments to retrieve file paths and options.
#' 3. Sourcing required R scripts that contain necessary functions.
#' 4. Reading in data files specified by the user through command-line arguments.
#' 5. Generating visualizations based on the processed data.
#'
#' @return None
#' @examples
#' # To run the script from the command line:
#' Rscript your_script_name.R --files "file1.csv,file2.csv" --input_path "data/" --output_path "plots/"
main <- function(files_to_source, files_to_read, files_in_path, files_out_path, plot_ext, plot_device, nifh_cluster_colours) {
  # Load the data
  data_list <- load_files(files_to_read, files_in_path)

  # Generate the color blind safe paletter for nifH clusters
  cluster_palette <- generate_nifh_palette(nifh_cluster_colours = nifh_cluster_colours)

  # Generate and save plots
  pie_plot_list <- pie_charts(
    counts_data =
      data_list$nifhdb_all_counts_AUID_dedup_clean,
    rel_abund_data =
      data_list$nifhdb_all_rel_abund_AUID_dedup_clean,
    output_dir = files_out_path,
    colors = cluster_palette,
    plot_ext = plot_ext,
    plot_device = plot_device
  )

  bar_plot_list <- bar_plots(
    counts_data_st_id_clean =
      data_list$nifhdb_all_counts_AUID_dedup_study_id_total_clean,
    rel_abund_data_st_id_clean =
      data_list$nifhdb_all_rel_abund_AUID_dedup_study_id_total_clean,
    counts_data_st_id = data_list$nifhdb_all_counts_AUID_dedup_study_id_total,
    rel_abund_data_st_id =
      data_list$nifhdb_all_rel_abund_AUID_dedup_study_id_total,
    counts_data = data_list$nifhdb_all_counts_AUID_dedup_clean,
    rel_abund_data = data_list$nifhdb_all_rel_abund_AUID_dedup_clean,
    output_dir = files_out_path,
    colors = cluster_palette,
    plot_ext = plot_ext,
    plot_device = plot_device
  )

  scatter_results_list <- make_scatter(
    rel_abund_df = data_list$RA_df_T_lng_mean_RA_AUID_deduped,
    anno_table = data_list$annoNifHDB_updt,
    meta_table = data_list$cmap_coloc,
    output_dir = files_out_path,
    colors = cluster_palette,
    plot_ext = plot_ext,
    plot_device = plot_device
  )

  if (!any(map_lgl(c(pie_plot_list, bar_plot_list, scatter_results_list$data, scatter_results_list$plots), is.null))) {
    # cat("\033[32m[SUCCESS] Pipeline completed successfully\033[0m\n")
    message("\nAll plots generated successfully.\n")
  } else {
    message("\nThere was an error and not all plots were generated\n")
  }
}



# Run if the script is being run directly
if (sys.nframe() == 0 && !interactive()) {
  source_file(files_to_source)

  parser <- setup_parser()
  args <- parse_arg(parser)

  validate_parsed_args(args)

  # Create output directory if it doesn't exist
  create_dir(args$files_out_path)

  final_results <- main(
    files_to_source = files_to_source,
    files_to_read = args$files_to_read,
    files_in_path = args$files_in_path,
    files_out_path = args$files_out_path,
    plot_ext = args$plot_ext,
    plot_device = args$plot_device,
    nifh_cluster_colours = nifh_cluster_colours
  )
}
