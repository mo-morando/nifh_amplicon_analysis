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
  parser <- arg_parser("Process the workspace files")
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

  # Convert the comma-separated string to a vector
  files_to_read <- strsplit(argv$files, ",")[[1]]
  files_in_path <- argv$input_path
  files_out_path <- argv$output_path
  plot_ext <- argv$plot_ext

  return(list(
    files_to_read = files_to_read,
    files_in_path = files_in_path,
    files_out_path = files_out_path,
    plot_ext = plot_ext
  ))
}



# #' Validate parsed arguments
# #'
# #' Validate the objects returned by the parse_arg() function
# #'
# #' @param parsed_args A list containing the parsed arguments returned by parse_arg()
# #'
# #' @return NULL
# #'
# #' @example
# #' parser <- setup_parser()
# #' args <- parse_arg(parser)
# #' validate_parsed_args(args)
# validate_parsed_args <- function(parsed_args) {
#   # Check if files_to_read is a character vector
#   if (!is.character(parsed_args$files_to_read)) {
#     stop("files_to_read must be a character vector")
#   }

#   # Check if files_in_path is a valid directory path
#   if (!dir.exists(parsed_args$files_in_path)) {
#     stop("files_in_path: '", parsed_args$files_in_path, "' must be a valid directory path")
#   }

#   # Check if files_out_path is a valid directory path
#   if (!dir.exists(parsed_args$files_out_path)) {
#     warning("files_in_path: '", parsed_args$files_in_path, "' must be a valid directory path")
#   }

#   # Check if each file in files_to_read exists in files_in_path
#   missing_files <- character(0)
#   for (file in parsed_args$files_to_read) {
#     # if (!file.exists(file.path(parsed_args$files_in_path, file))) {
#     if (!any(startsWith(list.files(parsed_args$files_in_path), paste0(file, ".")) )) {
#       missing_files <- c(missing_files, file)
#     }
#   }
#   if (length(missing_files > 0)) {
#     stop("The following files are missing in files_in_path:\n", paste("\t", missing_files, collapse = "\n"))
#   }

#   # All validations passed
#   cat("All parsed arguments are valid.\n")
# }





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











#' Create a Pie Chart for nifH Cluster Percentages
#'
#' This function creates a pie chart for nifH cluster percentages.
#'
#' @param data Data frame containing the plot data.
#' @param colors Named vector of colors for the nifH clusters.
#' @return A ggplot2 object representing the pie chart.
#'
create_pie_chart <- function(data, colors) {
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



#' Generate Pie Charts for nifH Cluster Percentages
#'
#' This function creates pie charts for nifH cluster percentages based on counts and relative abundance.
#' It saves the generated plots as image files with the specified extension.
#'
#' @param counts_data Data frame containing counts data for nifH clusters.
#' @param rel_abund_data Data frame containing relative abundance data for nifH clusters.
#' @param output_dir Directory where the plots will be saved.
#' @param colors Named vector of colors for the nifH clusters.
#' @param plot_ext File extension for saving plots (e.g., ".jpeg", ".png").
#' @return None
#' @examples
#' pie_charts(nifhdb_all_counts_AUID_dedup_clean, nifhdb_all_rel_abund_AUID_dedup_clean, "/path/to/output", nifh_cluster_colours_colbldsafe, ".jpeg")
pie_charts <- function(counts_data, rel_abund_data, output_dir, colors, plot_ext = ".jpeg") {
  # Check if the output directory exists
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist.")
  }

  # Setup success flag for if all plots are successful
  success <- TRUE


  # Create and save pie chart for counts
  tryCatch(
    {
      cat("Creating pie chart for counts data...\n")
      pie_chart_counts <- create_pie_chart(counts_data, colors)
      print(pie_chart_counts)
      ggsave(filename = file.path(output_dir, paste0("DNA_perc_cnts_clst_pie", plot_ext)), plot = pie_chart_counts, height = 8.5, width = 14, units = "in", dpi = 300)
      cat("Pie chart for counts data saved successfully.\n")
    },
    warning = function(w) {
      cat("Warning occurred:\n")
      cat("Warning message:", conditionMessage(w), "\n")
      cat("Warning call:", deparse(conditionCall(w)), "\n")
    },
    error = function(e) {
      cat("An error occurred in pie_charts() while creating the counts pie chart:\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Error call:", deparse(conditionCall(e)), "\n")
      # return(NULL)
      success <<- FALSE
    }
  )


  # Create and save pie chart for relative abundance
  tryCatch(
    {
      cat("Creating pie chart for relative abundance data...\n")
      pie_chart_rel_abund <- create_pie_chart(rel_abund_data, colors)
      print(pie_chart_rel_abund)
      ggsave(filename = file.path(output_dir, paste0("DNA_perc_rel_abund_clst_pie", plot_ext)), plot = pie_chart_rel_abund, height = 8.5, width = 14, units = "in", dpi = 300)
      cat("Pie chart for relative abundance data saved successfully.\n")
    },
    warning = function(w) {
      cat("Warning occurred:\n")
      cat("Warning message:", conditionMessage(w), "\n")
      cat("Warning call:", deparse(conditionCall(w)), "\n")
    },
    error = function(e) {
      cat("An error occurred in pie_charts() while creating the relative abundance pie chart:\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Error call:", deparse(conditionCall(e)), "\n")
      # return(NULL)
      success <<- FALSE
    }
  )


  # Combine and save combined plot
  tryCatch(
    {
      cat("Combining pie charts...\n")
      combined_plot <- (pie_chart_counts | pie_chart_rel_abund)
      print(combined_plot)
      ggsave(filename = file.path(output_dir, paste0("DNA_perc_combined_clst_pie", plot_ext)), plot = combined_plot, height = 8.5, width = 16, units = "in", dpi = 300)
      cat("Combined pie charts saved successfully.\n")
    },
    warning = function(w) {
      cat("Warning occurred:\n")
      cat("Warning message:", conditionMessage(w), "\n")
      cat("Warning call:", deparse(conditionCall(w)), "\n")
    },
    error = function(e) {
      cat("An error occurred in pie_charts() while combining the pie charts:\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Error call:", deparse(conditionCall(e)), "\n")
      # return(NULL)
      success <<- FALSE
    }
  )


  return(success)
}


#' Generate Bar Plots for nifH Cluster Percentages
#'
#' This function creates bar plots for nifH cluster percentages based on counts and relative abundance.
#' It saves the generated plots as image files with the specified extension.
#'
#' @param counts_data Data frame containing counts data for nifH clusters.
#' @param rel_abund_data Data frame containing relative abundance data for nifH clusters.
#' @param output_dir Directory where the plots will be saved.
#' @param colors Named vector of colors for the nifH clusters.
#' @param plot_ext File extension for saving plots (e.g., ".jpeg", ".png").
#' @return None
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
    plot_ext = ".jpeg") {
  
  # Check if the output directory exists
  cat("Checking if the output directory exists...\n")
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist.")
  }

  # Setup success flag for if all plots work
  success <- TRUE


  # Create and save bar plot for counts
  tryCatch(
    {
      cat("Creating bar plot for counts data...\n")
      bar_plot_counts <- bar_plot(
        df = counts_data_st_id_clean,
        x = studyID,
        y = percentage_total_counts_nifH_cluster_study_id_total_clean,
        # fill = cluster_stats,
        fill = nifH_cluster_modified,
        fill_pallete = colors,
        x_lab = "Study ID",
        y_lab = bquote(bold("%" ~ total ~ reads)),
        legend_position = "bottom",
        x_axis_angle = TRUE,
        print_out = TRUE
      )
      ggsave(file.path(output_dir, paste0("DNA_perc_cnts_clst_bar", plot_ext)), plot = bar_plot_counts, height = 8.5, width = 14, units = "in", dpi = 300)
      cat("Bar plot for counts data saved successfully.\n")
    },
    warning = function(w) {
      cat("Warning occurred:\n")
      cat("Warning message:", conditionMessage(w), "\n")
      cat("Warning call:", deparse(conditionCall(w)), "\n")
    },
    error = function(e) {
      cat("An error occurred in bar_plots() while creating the counts bar plot:\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Error call:", deparse(conditionCall(e)), "\n")
      success <<- FALSE
    }
  )

  # Create and save bar plot for relative abundance
  tryCatch(
    {
      cat("Creating bar plot for relative abundance data...\n")
      bar_plot_rel_abund <- bar_plot(
        df = rel_abund_data_st_id_clean,
        x = studyID,
        y = percentage_total_rel_abund_nifH_cluster_study_id_total_clean,
        # fill = cluster_stats,
        fill = nifH_cluster_modified,
        fill_pallete = colors,
        x_lab = "Study ID",
        y_lab = bquote(bold("%" ~ relative ~ abundance)),
        legend_position = "bottom",
        x_axis_angle = TRUE,
        print_out = TRUE
      )
      ggsave(file.path(output_dir, paste0("DNA_perc_rel_abund_clst_bar", plot_ext)), plot = bar_plot_rel_abund, height = 8.5, width = 14, units = "in", dpi = 300)
      cat("Bar plot for relative abundance data saved successfully.\n")
    },
    warning = function(w) {
      cat("Warning occurred:\n")
      cat("Warning message:", conditionMessage(w), "\n")
      cat("Warning call:", deparse(conditionCall(w)), "\n")
    },
    error = function(e) {
      cat("An error occurred in bar_plots() while creating the relative abundance bar plot:\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Error call:", deparse(conditionCall(e)), "\n")
      success <<- FALSE
    }
  )


  # Create a plot with all the pooled data as a column for better comparison
  # First, counts data
  tryCatch(
    {
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

      #* # make Tibble that joins the pooled data with study ID data
      #* # rename columns to be the same for binding
      #* # make new studyId for pooled data
      cat("Combining pooled data with study Id data...\n")
      combined_data <- counts_data_st_id %>%
        bind_rows(counts_data %>%
          rename(
            # nifH_cluster = nifH_cluster_modified,
            percentage_total_counts_nifH_cluster_study_id_total = percentage_nifH_cluster
          ) %>%
          mutate(
            studyID = "pooled data", #* # make new studyId for pooled data
            # cluster_stats = nifH_cluster #* # this is needed for plotting
            cluster_stats = nifH_cluster_modified #* # this is needed for plotting
          )) %>%
        mutate(
          studyID = factor(studyID, levels = custom_order_w_total) #* # make new plotting order
        )
      cat("Successfully combined data...\n")

      # return(combined_data)
    },
    warning = function(w) {
      cat("Warning occurred while combining study data:\n")
      cat("Warning message:", conditionMessage(w), "\n")
      cat("Warning call:", deparse(conditionCall(w)), "\n")
      # return(combined_data)
    },
    error = function(e) {
      cat("An error occurred while combining study data:\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Error call:", deparse(conditionCall(e)), "\n")
      # return(NULL)
      success <<- FALSE
    }
  )

  # Create bar plot using the combined data
  if (!is.null(combined_data)) {
    cat("Creating combined bar plot of pooled and study ID data...\n")
    tryCatch(
      {
        bar_plot_obj <- bar_plot(
          df = combined_data,
          x = studyID,
          # y = percentage_total_counts_nifH_cluster_study_id_total,
          y = percentage_total_counts_nifH_cluster_study_id_total,
          # fill = nifH_cluster_modified,
          fill = cluster_stats,
          fill_pallete = colors,
          # fill_lab = expression(italic("nifH") "cluster"),
          fill_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
          x_lab = "Study ID",
          # y_lab = "% of total nifH cluster counts per study ID",
          y_lab = bquote(bold("%" ~ total ~ reads)),
          legend_position = "bottom",
          x_axis_angle = TRUE,
          print_out = TRUE
        )

        ggsave(file.path(output_dir, paste0("DNA_perc_cnts_clst_all_bar", plot_ext)), plot = bar_plot_obj, height = 8.5, width = 14, units = "in", dpi = 300)
        cat("Combined bar plot of pooled and study ID data saved successfully.\n")
      },
      warning = function(w) {
        cat("Warning occurred while creating bar plot:\n")
        cat("Warning message:", conditionMessage(w), "\n")
        cat("Warning call:", deparse(conditionCall(w)), "\n")
      },
      error = function(e) {
        cat("An error occurred while creating bar plot:\n")
        cat("Error message:", conditionMessage(e), "\n")
        cat("Error call:", deparse(conditionCall(e)), "\n")
        success <<- FALSE
      }
    )
  }


  # Next, relative abundance data
  tryCatch(
    {
      #* # make Tibble that joins the pooled data with study ID data
      #* # rename columns to be the same for binding
      #* # make new studyId for pooled data
      cat("Combining pooled data with study Id data...\n")
      combined_data <- rel_abund_data_st_id %>%
        bind_rows(rel_abund_data %>%
          rename(
            # nifH_cluster = nifH_cluster_modified,
            percentage_total_rel_abund_nifH_cluster_study_id_total = percentage_nifH_cluster
          ) %>%
          mutate(
            studyID = "pooled data", #* # make new studyId for pooled data
            # cluster_stats = nifH_cluster #* # this is needed for plotting
            cluster_stats = nifH_cluster_modified #* # this is needed for plotting
          )) %>%
        mutate(
          studyID = factor(studyID, levels = custom_order_w_total) #* # make new plotting order
        )
      cat("Successfully combined data...\n")

      # return(combined_data)
    },
    warning = function(w) {
      cat("Warning occurred while combining study data:\n")
      cat("Warning message:", conditionMessage(w), "\n")
      cat("Warning call:", deparse(conditionCall(w)), "\n")
      # return(combined_data)
    },
    error = function(e) {
      cat("An error occurred while combining study data:\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Error call:", deparse(conditionCall(e)), "\n")
      # return(NULL)
      success <<- FALSE
    }
  )

  # Create bar plot using the combined data
  if (!is.null(combined_data)) {
    cat("Creating combined bar plot of pooled and study ID data...\n")
    tryCatch(
      {
        bar_plot_obj <- bar_plot(
          df = combined_data,
          x = studyID,
          # y = percentage_total_rel_abund_nifH_cluster_study_id_total,
          y = percentage_total_rel_abund_nifH_cluster_study_id_total,
          # fill = nifH_cluster_modified,
          fill = cluster_stats,
          fill_pallete = colors,
          # fill_lab = expression(italic("nifH") "cluster"),
          fill_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
          x_lab = "Study ID",
          # y_lab = "% of total nifH cluster counts per study ID",
          y_lab = bquote(bold("%" ~ relative ~ abundance)),
          legend_position = "bottom",
          x_axis_angle = TRUE,
          print_out = TRUE
        )

        ggsave(file.path(output_dir, paste0("DNA_perc_rel_abund_clst_all_bar", plot_ext)), plot = bar_plot_obj, height = 8.5, width = 14, units = "in", dpi = 300)
        cat("Combined bar plot of pooled and study ID data saved successfully.\n")
      },
      warning = function(w) {
        cat("Warning occurred while creating bar plot:\n")
        cat("Warning message:", conditionMessage(w), "\n")
        cat("Warning call:", deparse(conditionCall(w)), "\n")
      },
      error = function(e) {
        cat("An error occurred while creating bar plot:\n")
        cat("Error message:", conditionMessage(e), "\n")
        cat("Error call:", deparse(conditionCall(e)), "\n")
        success <<- FALSE
      }
    )
  }

  return(success)
}



#' Create Scatter Plots for Relative Abundance Data
#'
#' This function generates scatter plots for relative abundance data against absolute latitude and sea surface temperature (SST).
#'
#' @param rel_abund_df Data frame containing relative abundance data.
#' @param anno_table Data frame containing annotation information.
#' @param meta_table Data frame containing metadata information.
#' @param colors Named vector of colors for plotting.
#' @param plot_ext File extension for saving plots (e.g., ".jpeg").
#' @param output_dir Directory where the plots will be saved.
#'
#' @return Logical indicating success (TRUE) or failure (FALSE).
# make_scatter <- function(
#     rel_abund_df,
#     anno_table,
#     meta_table,
#     colors,
#     plot_ext,
#     output_dir) {

#   # Check if the output directory exists
#   cat("Checking if the output directory exists...\n")
#   if (!dir.exists(output_dir)) {
#     stop("Output directory does not exist")
#   }

#   # # Setup succes flag
#   # success <- TRUE
#   success <- FALSE

#   tryCatch({
#     # Merge relative abundance data with annotations and metadata
#     cat("Merging relative abundance data with annotations and metadata...\n")
#     query_df <- rel_abund_df %>%
#       #   add_group_id()  %>%
#       merge_annotations(annotation_table = anno_table) %>%
#       #   merge_cmap(by_join = c("SAMPLEID", "group_id"))
#       merge_cmap(cmap = meta_table, by_join = c("SAMPLEID"))
#       cat("Successfully merge files and generated query_df\n")

#     # Sum each sample by nifH cluster and select relevant columns
#     cat("Summing each sample by nifH cluster and selecting relevant columns...\n")
#     query_df_cluster <- query_df %>%
#       group_by(SAMPLEID, nifH_cluster) %>%
#       mutate(RA = sum(RA, na.rm = T)) %>%
#       select(SAMPLEID, studyID, time, depth, lat, lon, nifH_cluster, RA, CyanoCON, lat_abs, hemi, everything()) %>%
#       ungroup() %>%
#       distinct(SAMPLEID, nifH_cluster, .keep_all = T) %>%
#       ungroup()

#     # Filter for specific nifH clusters of interest
#     cat("Filtering for specific nifH clusters of interest...\n")
#     query_df_cluster_sub <- query_df_cluster %>%
#       filter(nifH_cluster %in% c("1J/1K", "3", "1A", "1B", "1G", "1O/1P"))

#     # Create scatter plot: Relative Abundance vs Absolute Latitude
#     cat("Creating scatter plot: Relative Abundance vs Absolute Latitude...\n")
#     ra_v_abslat_plot <- scatter_line_plot(
#       df = query_df_cluster_sub,
#       y = RA,
#       x = lat_abs,
#       colour = nifH_cluster,
#       group = nifH_cluster,
#       colour_pallete = colors,
#       colour_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
#       fill_lab = NULL,
#       fill = nifH_cluster,
#       fill_pallete = colors,
#       title_lab = NULL,
#       subtitle_lab = NULL,
#       x_lab = "absolute latitude",
#       y_lab = "% relative abundance",
#       legend_position = "bottom",
#       legend_direction = "horizontal"
#     )


#     # Create scatter plot: Relative Abundance vs SST
#     cat("Creating scatter plot: Relative Abundance vs SST...\n")
#     ## RA versus SST
#     ra_v_sst_plot <- scatter_line_plot(
#       df = query_df_cluster_sub,
#       y = RA,
#       x = sst_tblSST_AVHRR_OI_NRT,
#       colour = nifH_cluster,
#       group = nifH_cluster,
#       colour_pallete = colors,
#       colour_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
#       fill_lab = NULL,
#       fill = nifH_cluster,
#       fill_pallete = colors,
#       title_lab = NULL,
#       subtitle_lab = NULL,
#       x_lab = "SST ˚C",
#       y_lab = "% relative abundance",
#       legend_position = "bottom",
#       legend_direction = "horizontal"
#     ) +
#       scale_x_reverse()


#     # Combine plots using patchwork
#     cat("Combining plots using patchwork...\n")
#     ra_v_abslat_plot_sub <- ra_v_abslat_plot +
#       theme(legend.position = "none")
#     ra_v_sst_plot_sub <- ra_v_sst_plot
#     combined_sct_plot <- ra_v_abslat_plot_sub / ra_v_sst_plot_sub

#   # cat("Saving the combined scatter plot temp...\n")
#   #   ggsave(
#   #     filename = file.path(output_dir, paste0("plot1", plot_ext)),
#   #   plot = ra_v_sst_plot, height = 8.5, width = 14, units = "in", dpi = 300)

#     # Save the combined scatter plot
#     cat("Saving the combined scatter plot...\n")
#     ggsave(file.path(output_dir, paste0("NifHsumCld_scat_bi_RAv_lat_SST_noLegend", plot_ext)),
#       plot = combined_sct_plot, height = 12.5, width = 21, units = "in",
#       dpi = 300
#     )

#     # Setup succes flag
#     success <- TRUE

#   },
#   warning = function(w) {
#     cat("A warning occured \n")
#     cat("Warning messege:", conditionMessage(w), "\n")
#     cat("Warning call:", deparse(conditionCall(w)), "\n")
#   },
#   error = function(e) {
#     cat("An error occurred running make_scatter()")
#     cat("Error message:", conditionMessage(e), "\n")
#     cat("Error call:", deparse(conditionCall(e)), "\n")
#     success <<- FALSE
#   })


#   return(success)
# }












#' Main function to execute plots
main <- function(files_to_source, files_to_read, files_in_path, files_out_path, plot_ext, nifh_cluster_colours) {
  # Load the data
  data_list <- load_files(files_to_read, files_in_path)

  # Generate the color blind safe paletter for nifH clusters
  cluster_palette <- generate_nifh_palette(nifh_cluster_colours = nifh_cluster_colours)

  # Generate and save plots
  success_pie <- pie_charts(
    counts_data =
      data_list$nifhdb_all_counts_AUID_dedup_clean,
    rel_abund_data =
      data_list$nifhdb_all_rel_abund_AUID_dedup_clean,
    output_dir = files_out_path,
    colors = cluster_palette,
    plot_ext = plot_ext
  )

  success_bar <- bar_plots(
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
    plot_ext = plot_ext
  )

  # success_scatter <- make_scatter(
  #   rel_abund_df = data_list$RA_df_T_lng_mean_RA_AUID_deduped,
  #   anno_table = data_list$annoNifHDB_updt,
  #   meta_table = data_list$cmap_coloc,
  #   output_dir = files_out_path,
  #   colors = cluster_palette,
  #   plot_ext = plot_ext
  # )


  # if (success_pie && success_bar && success_scatter) {
  if (success_pie && success_bar) {
    message("\n\nAll plots generated successfully.\n")
  } else {
    message("\n\nThere was an error and not all plots were generated
    successfully\n")
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
    nifh_cluster_colours = nifh_cluster_colours
  )
}















# ### - pie chart of nifh cluster percentages

# #### by counts
# #### just DNA, deduplicated by AUID and group_id
# pie_chart_AUID_depud_DNA_counts <- ggplot(nifhdb_all_counts_AUID_dedup_clean, aes(
#   x = 1,
#   # x = studyID,
#   y = percentage_nifH_cluster,
#   fill = nifH_cluster_modified,
#   # colour = studyID,
#   # group = nifH_cluster
# )) +
#   geom_bar(
#     stat = "identity",
#     # position = position_stack(reverse = FALSE)
#   ) +
#   # geom_histogram(binwidth = 1) +
#   coord_polar(theta = "y") +
#   theme_void() +
#   labs(
#     # title = "nifH dabtbase - % of total counts for each clusters",
#     # subtitle = "DNA only with replicates averaged",
#     fill = "nifH Cluster"
#   ) +

#   # Add labels with lines for small slices
#   geom_text(
#     # aes(label = ifelse(percentage_nifH_cluster >= 5, paste(nifH_cluster, "\n", percentage_nifH_cluster, "%"), ""), x = 1.25),
#     aes(label = ifelse(percentage_nifH_cluster >= 5, paste0(percentage_nifH_cluster, "%"), ""), x = 1.15),
#     position = position_stack(vjust = 0.5),
#     hjust = 0.5,
#     size = 12,
#     show.legend = FALSE
#   ) +
#   geom_text(
#     # aes(label = ifelse(percentage_nifH_cluster >= 5, paste(nifH_cluster, "\n", percentage_nifH_cluster, "%"), ""), x = 1.25),
#     aes(
#       label = ifelse(percentage_nifH_cluster < 5, paste0(percentage_nifH_cluster, "%"), ""),
#       x = 1.5
#     ),
#     position = position_stack(vjust = 0.5),
#     hjust = 0.5,
#     size = 5,
#     show.legend = FALSE,
#     # nudge_y = 0.05
#   ) +
#   scale_fill_manual(values = nifh_cluster_colours_colbldsafe) +
#   theme(
#     legend.position = "bottom",
#     legend.title = element_text(size = 20, face = "bold"),
#     legend.text = element_text(size = 17, face = "bold"),
#   ) +
#   guides(fill = guide_legend(
#     nrow = 1,
#     byrow = TRUE,
#     title.position = "top",
#     title.hjust = 0.5
#   ))

# # pie_chart_facet <- pie_chart +
# #   facet_wrap("studyID")

# # Display the pie chart
# print(pie_chart_AUID_depud_DNA_counts)

# # print(pie_chart_facet)

# ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/pie_chart_nifhDF_DNA_dedup_perc_tot_cnts_nifH_cluster.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)

# annotation_table %>%
#   filter(nifH_cluster == "1O/1P") %>%
#   distinct(CON)

# #### by relative abundance
# #### just DNA, deduplicated by AUID and group_id
# pie_chart_AUID_depud_DNA_rel_abund <- ggplot(
#   nifhdb_all_rel_abund_AUID_dedup_clean,
#   aes(
#     x = 1,
#     # x = studyID,
#     y = percentage_nifH_cluster,
#     fill = nifH_cluster_modified,
#     # colour = studyID,
#     # group = nifH_cluster
#   )
# ) +
#   geom_bar(
#     stat = "identity",
#     # position = position_stack(reverse = FALSE)
#   ) +
#   scale_fill_manual(values = nifh_cluster_colours_colbldsafe) +
#   # geom_histogram(binwidth = 1) +
#   coord_polar(theta = "y") +
#   theme_void() +
#   labs(
#     # title = "nifH dabtbase - % of total relative abundance for each clusters",
#     # subtitle = "DNA only with replicates averaged",
#     fill = "nifH Cluster"
#   ) +

#   # Add labels with lines for small slices
#   geom_text(
#     # aes(label = ifelse(percentage_nifH_cluster >= 5, paste(nifH_cluster, "\n", percentage_nifH_cluster, "%"), ""), x = 1.25),
#     aes(label = ifelse(percentage_nifH_cluster >= 5, paste0(percentage_nifH_cluster, "%"), ""), x = 1.15),
#     position = position_stack(vjust = 0.5),
#     hjust = 0.5,
#     size = 12,
#     show.legend = FALSE
#   ) +
#   geom_text(
#     # aes(label = ifelse(percentage_nifH_cluster >= 5, paste(nifH_cluster, "\n", percentage_nifH_cluster, "%"), ""), x = 1.25),
#     aes(
#       label = ifelse(percentage_nifH_cluster < 5, paste0(percentage_nifH_cluster, "%"), ""),
#       x = 1.5
#     ),
#     position = position_stack(vjust = 0.5),
#     hjust = 0.5,
#     size = 5,
#     show.legend = FALSE,
#     # nudge_y = 0.05
#   )

# # pie_chart_facet <- pie_chart +
# #   facet_wrap("studyID")

# # Display the pie chart
# print(pie_chart_AUID_depud_DNA_rel_abund)
# # print(pie_chart_facet)

# ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/pie_chart_nifhDF_DNA_dedup_perc_tot_rel_abund_nifH_cluster.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)


# ### - combine plots

# ## * make new plots that are friendly to combining
# pie1 <- pie_chart_AUID_depud_DNA_counts +
#   theme(
#     legend.position = "bottom", # Set the legend position to "bottom"
#     legend.justification = c(1, 0), # Adjust the justification to place it at the bottom-right
#     legend.box.just = "right"
#   )

# pie2 <- pie_chart_AUID_depud_DNA_rel_abund +
#   theme(
#     legend.position = "none"
#   )

# ## * combine plots with patchwork
# # combined_plot_pie_nifH_cluster_abundance <- pie1 + pie2
# combined_plot_pie_nifH_cluster_abundance <- (pie1 | pie2)

# print(combined_plot_pie_nifH_cluster_abundance)

# ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/pie_chart_nifhDF_DNA_dedup_perc_abund_nifH_cluster_cmb_plot.jpeg", height = 8.5, width = 16, units = "in", dpi = 300)






# ### - bar chart of all samples
# ### - bar chart of all samples
# ### - bar chart of all samples

# #### * by counts
# #### * just DNA, deduplicated by AUID and group_id
# #### * percentages based on total samples

# # ## ! If I use this tibble, I can show the individual contributions of each study ID, but if I want solid bars that don't have the lines within each bar that show each study ID, this doesn't work. Instead use futher cleaned up tibble below
# # to_plot <- nifhdb_all_counts_AUID_dedup_total_study_id
# # ## ! Instead use futher cleaned up tibble just below here:

# nifhdb_all_counts_AUID_dedup_total_study_id_clean_nolines <- nifhdb_all_counts_AUID_dedup_total_study_id_clean %>%
#   group_by(studyID) %>%
#   summarise("percentage_total_counts_nifH_cluster_total_study_id_clean" = sum(percentage_total_counts_nifH_cluster_total_study_id_clean))
# to_plot <- nifhdb_all_counts_AUID_dedup_total_study_id_clean_nolines

# # viridis_color_pallete <- get_viridis_colors(to_plot, nifH_cluster, "H", -1, 0)

# (nifhdb_all_counts_AUID_dedup_total_study_id_bar_plot <- bar_plot(
#   df = to_plot,
#   x = studyID,
#   # y = percentage_total_counts_nifH_cluster_total_study_id_clean,
#   y = percentage_total_counts_nifH_cluster_total_study_id_clean,
#   fill = NULL,
#   x_lab = "Study ID",
#   # y_lab = "% of total nifH",
#   # y_lab = bquote(bold("%"~total ~ bold(italic(nifH)))),
#   y_lab = bquote(bold("%" ~ total ~ reads)),
#   legend_position = "none"
# ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), legend.direction = "vertical"))

# ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/bar_chart_nifhDB_DNA_dedup_perc_tot_cnts_nifH_cluster_total_studyID.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)


# ####* by relative abundance
# ####* just DNA, deduplicated by AUID and group_id
# ####* percentages based on total samples

# # ## ! If I use this tibble, I can show the individual contributions of each study ID, but if I want solid bars that don't have the lines within each bar that show each study ID, this doesn't work. Instead use futher cleaned up tibble below
# # to_plot <- nifhdb_all_rel_abund_AUID_dedup_total_study_id
# # ## ! Instead use futher cleaned up tibble just below here:

# nifhdb_all_rel_abund_AUID_dedup_total_study_id_clean_nolines <- nifhdb_all_rel_abund_AUID_dedup_total_study_id_clean %>%
#   group_by(studyID) %>%
#   summarise("percentage_total_rel_abund_total_nifH_cluster_study_id_clean" = sum(percentage_total_rel_abund_total_nifH_cluster_study_id_clean))

# to_plot <- nifhdb_all_rel_abund_AUID_dedup_total_study_id_clean_nolines
# # viridis_color_pallete <- get_viridis_colors(to_plot, nifH_cluster, "H", -1, 0)

# (nifhdb_all_rel_abund_AUID_dedup_total_study_id_bar_plot <- bar_plot(
#   df = to_plot,
#   x = studyID,
#   # y = percentage_total_rel_abund_total_nifH_cluster_study_id,
#   y = percentage_total_rel_abund_total_nifH_cluster_study_id_clean,
#   # fill = nifH_cluster,
#   # fill_pallete = viridis_color_pallete,
#   # fill = nifH_cluster_modified,
#   # fill_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
#   x_lab = "Study ID",
#   y_lab = bquote(bold("%" ~ relative ~ abundance)),
#   # y_lab = "% of total nifH cluster relative abundance per study ID",
#   legend_position = "none"
# ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), legend.direction = "vertical"))

# ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/bar_chart_nifhDB_DNA_dedup_perc_tot_RA_nifH_cluster_total_studyID.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)


# #### by counts
# #### just DNA, deduplicated by AUID and group_id
# #### percentages recalculated based on individual studyID samples totals

# # to_plot <- nifhdb_all_counts_AUID_dedup_study_id_total_clean
# custom_order_w_total <- c(
#   "pooled data",
#   "AK2HI",
#   "BentzonTilia_2015",
#   "Ding_2021",
#   "Gradoville_2020_G1",
#   "Gradoville_2020_G2",
#   "Hallstrom_2021",
#   "Hallstrom_2022",
#   "Harding_2018",
#   "Mulholland_2018",
#   "NEMO",
#   "Raes_2020",
#   "Sato_2021",
#   "Selden_2021",
#   "Shiozaki_2017",
#   "Shiozaki_2018GBC",
#   "Shiozaki_2018LNO",
#   "Shiozaki_2020",
#   "Tang_2020",
#   "Wu_2019",
#   "Wu_2021",
#   "TurkKubo_2021"
# )

# #* # make tibble that joins the pooled data with study ID data
# #* # rename columns to be the same for binding
# #* # make new studyId for pooled data
# temp <- nifhdb_all_counts_AUID_dedup_study_id_total %>%
#   bind_rows(nifhdb_all_counts_AUID_dedup_clean %>%
#     rename(
#       nifH_cluster = nifH_cluster_modified,
#       percentage_total_counts_nifH_cluster_study_id_total = percentage_nifH_cluster
#     ) %>%
#     mutate(
#       studyID = "pooled data", #* # make new studyId for pooled data
#       cluster_stats = nifH_cluster #* # this is needed for plotting
#     )) %>%
#   mutate(
#     studyID = factor(studyID, levels = custom_order_w_total) #* # make new plotting order
#   )

# to_plot <- temp




# # viridis_color_pallete <- get_viridis_colors(to_plot, nifH_cluster, "H", -1, 0)

# nifhdb_all_counts_AUID_dedup_study_id_total_bar_plot <- bar_plot(
#   df = to_plot,
#   x = studyID,
#   # y = percentage_total_counts_nifH_cluster_study_id_total,
#   y = percentage_total_counts_nifH_cluster_study_id_total,
#   # fill = nifH_cluster_modified,
#   fill = cluster_stats,
#   fill_pallete = nifh_cluster_colours_colbldsafe,
#   # fill_lab = expression(italic("nifH") "cluster"),
#   fill_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
#   x_lab = "Study ID",
#   # y_lab = "% of total nifH cluster counts per study ID",
#   y_lab = bquote(bold("%" ~ total ~ reads)),
#   legend_position = "bottom",
#   x_axis_angle = TRUE,
#   print_out = TRUE
# )

# ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/bar_chart_nifhDB_DNA_dedup_perc_tot_cnts_nifH_cluster_studyID_total.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)

# #### by relative abundance
# #### just DNA, deduplicated by AUID and group_id
# #### percentages recalculated based on individual studyID samples totals

# temp <- nifhdb_all_rel_abund_AUID_dedup_study_id_total %>%
#   left_join(
#     metatable %>%
#       distinct(studyID, study_ocean)
#   ) %>%
#   bind_rows(nifhdb_all_rel_abund_AUID_dedup_clean %>%
#     rename(
#       nifH_cluster = nifH_cluster_modified,
#       percentage_total_rel_abund_nifH_cluster_study_id_total = percentage_nifH_cluster
#     ) %>%
#     mutate(
#       studyID = "pooled data", #* # make new studyId for pooled data
#       cluster_stats = nifH_cluster #* # this is needed for plotting
#     )) %>%
#   mutate(
#     studyID = factor(studyID, levels = custom_order_w_total) #* # make new plotting order
#   )

# to_plot <- temp

# # viridis_color_pallete <- get_viridis_colors(to_plot, nifH_cluster, "H", -1, 0)

#  nifhdb_all_rel_abund_AUID_dedup_study_id_total_bar_plot <- bar_plot(
#   df = to_plot,
#   x = studyID,
#   y = percentage_total_rel_abund_nifH_cluster_study_id_total,
#   # y = percentage_total_rel_abund_total_nifH_cluster_study_id,
#   # fill = nifH_cluster_modified,
#   fill = cluster_stats,
#   fill_pallete = nifh_cluster_colours_colbldsafe,
#   # fill_lab = expression(italic("nifH") "cluster"),
#   fill_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
#   x_lab = "Study ID",
#   y_lab = bquote(bold("%" ~ relative ~ abundance)),
#   # y_lab = "% of total nifH cluster relative abundance per study ID",
#   legend_position = "bottom",
#   legend_direction = "horizontal",
#   n_row = 1,
#   x_axis_angle = TRUE,
#   print_out = TRUE
# )
# #+
# #   geom_hline(yintercept = 20, color = "red") +
# #   geom_hline(yintercept = 28, color = "orange") +
# #   geom_hline(yintercept = 2.4, color = "magenta") +
# #   geom_hline(yintercept = 5, , color = "brown")
# # +
# #   guides(fill = guide_legend(
# #     nrow = 1,
# #     byrow = TRUE,
# #     title.position = "top",
# #     title.hjust = 0.5
# # )



# ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/bar_chart_nifhDB_DNA_dedup_perc_tot_RA_nifH_cluster_studyID_total.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)

# nifhdb_all_rel_abund_AUID_dedup_study_id_total_bar_plot +
#   facet_wrap(~study_ocean)
