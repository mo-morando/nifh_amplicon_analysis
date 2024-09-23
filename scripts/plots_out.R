#!/usr/bin/env Rscript

#' @title Plots Output Script
#' @description This script processes workspace files and generates various plots for nifH amplicon data analysis.
#' @details The script performs the following steps:
#' 1. Loads necessary libraries and sources required files
#' 2. Sets up command-line argument parsing
#' 3. Loads and processes input data
#' 4. Generates and saves various plots including:
#'    - Sample type plots
#'    - Ocean and hemisphere plots
#'    - Depth distribution plots
#'    - Oceanographic data plots (SST, PO4, Fe, etc.)
#' @author Michael Morando
#' @date 2023-09-10

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
#     warning("files_out_path: '", parsed_args$files_out_path, "' must be a valid directory path")
#   }

#   # Check if each file in files_to_read exists in files_in_path
#   missing_files <- character(0)
#   for (file in parsed_args$files_to_read) {
#     # if (!file.exists(file.path(parsed_args$files_in_path, file))) {
#     if (!any(startsWith(list.files(parsed_args$files_in_path), file))) {
#       missing_files <- c(missing_files, file)
#     }
#   }
#   if (length(missing_files) > 0) {
#     stop("The following files are missing in files_in_path:\n", paste("\t",missing_files, collapse = "\n"))
#   }

#   # All validations passed
#   cat("All parsed arguments are valid.\n")
# }



# #' Load and Assign CSV Files
# #'
# #' @param file_list Character vector of file names to load
# #' @param path Directory path for input files
# #' @param verbose Logical. Whether to print status messages
# #' @return List of loaded data frames or character vectors
# #' @importFrom readr read_csv
# #' @examples
# #' \dontrun{
# #' files <- c("data1", "data2", "single_line")
# #' data_list <- load_files(files, "path/to/csv/files")
# #' }
# load_files <- function(file_list, path, verbose = TRUE) {
#   data_list <- list()

#   for (file in file_list) {
#     file_path <- file.path(path, paste0(file, ".csv"))
#     tryCatch(
#       {
#         if (file.exists(file_path)) {
#           if (verbose) cat("Loading file:", file_path, "\n")

#           # Check if the file has only one line
#           lines <- suppressWarnings(readLines(file_path, n = 2)) # Read only the first two lines
#           if (length(lines) == 1) {
#             # File contains a single line, read as character vector
#             data_list[[file]] <- read_csv(file_path, show_col_types = FALSE, col_names = FALSE)
#             if (verbose) cat("  File loaded as a data frame but without column
#           headers.\n")
#           } else {
#             # File contains multiple lines, read as a data frame
#             data_list[[file]] <- read_csv(file_path, show_col_types = FALSE)
#             if (verbose) {
#               cat(
#                 "  File loaded as a data frame with",
#                 nrow(data_list[[file]]), "rows and",
#                 ncol(data_list[[file]]), "columns.\n"
#               )
#             }
#           }
#         } else {
#           warning(paste("File not found:", file_path))
#         }
#       },
#       error = function(e) {
#         cat("Error loading file", file_path, ":", conditionMessage(e), "\n")
#       }
#     )
#   }
#   if (verbose) cat("Finished loading", length(data_list), "file.\n")
#   return(data_list)
# }






#' Generate and save plots
#'
#' @param data_list List of loaded data frames
#' @param files_out_path Path to save the plots
#' @return NULL
generate_plots <- function(
    cmap_coloc,
    samples_per_nucleicAcidType, samples_per_nucleicAcidType_studyid,
    samples_per_photic, samples_per_photic_nucacid,
    sample_type,
    unique_sample_id_key,
    files_out_path,
    plot_ext) {
  # remake query df
  query_df <- dedup_by_group(
    cmap_coloc,
    group_id_key = unique_sample_id_key,
    group_id
  )

  ### * Plot
  to_plot <- sample_type

  viridis_color_pallete <- get_viridis_colors(to_plot, new_group, "H", -1, 0)

  bar_plot(
    df = to_plot,
    x = percentage,
    y = "sample type",
    # y = percentage_total_rel_abund_total_nifH_cluster_study_id_clean,
    fill = new_group,
    fill_pallete = viridis_color_pallete,
    # fill = nifH_cluster_modified,
    # fill_lab = expression(italic("nifH") "cluster"),
    # fill_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
    y_lab = NULL,
    x_lab = "% of total",
    legend_position = "right",
    x_axis_angle = TRUE,
    n_row = 10,
  )

  ### * combine nucleic acid type, light level, and sample type so you can plot them all together

  ##
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


  to_plot <- combined_nuca_phtc_smptp_tibble

  viridis_color_pallete <- get_viridis_colors(to_plot, new_group, "H", -1, 0)

  (combined_nuca_phtc_smptp_tibble_plot <- suppressMessages(bar_plot(
    df = to_plot,
    y = tibble_id,
    x = percentage,
    # y = percentage_total_rel_abund_total_nifH_cluster_study_id_clean,
    fill = new_group,
    fill_pallete = viridis_color_pallete,
    # fill = nifH_cluster_modified,
    # fill_lab = bquote((bold(italic(nifH)) ~ cluster)),
    y_lab = NULL,
    x_lab = "% of total",
    legend_position = "right",
    n_row = 10,
    x_axis_angle = TRUE
  )))

  ggsave(file.path(files_out_path, paste0("samp_nucleic_acid_type_light_bar_cmb_plot", plot_ext)), height = 10.5, width = 14, units = "in", dpi = 300)



  ### Study ID
  (samples_per_studyid <- count_and_arrange(
    query_df, c("studyID", "nucleicAcidType", "ocean")
  ))

  to_plot <- samples_per_studyid

  viridis_color_pallete <- get_viridis_colors(to_plot, ocean, "magma", -1, 0)

  (samples_per_studyid_plot <- bar_plot(
    df = to_plot,
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
  ) #+
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15))
  )

  samples_per_studyid_plot +
    facet_wrap(~nucleicAcidType)

  ggsave(file.path(files_out_path, paste0("Samples_per_studyID_ocean_fill_bar_facet_nucelic", plot_ext)), height = 8.5, width = 14, units = "in", dpi = 300)


  # ### * combine plot

  # samples_per_studyid_plot

  # combined_nuca_phtc_smptp_tibble_plot

  # (combine_test <- samples_per_studyid_plot | combined_nuca_phtc_smptp_tibble_plot +
  #   plot_layout(widths = c(1.5, 0))
  # #  plot_layout(widths = c(2, 1), heights = unit(c(5, 1), c('cm', 'null')))
  # )

  # ggsave(file.path(files_out_path, paste0("bar_chart_nifhDB_DNA_dedup_perc_tot_nifH_cluster_studyID_samples_per_cmb_plot", plot_ext)), height = 8.5, width = 14, units = "in", dpi = 300)



  ### -  Ocean and hemisphere
  samples_per_ocean <- count_and_arrange(
    query_df, c("ocean", "hemi")
  ) 

  # samples_per_ocean <- query_df  %>% 
  # count_and_arrange( c("ocean", "hemi")) %>%
  #   add_percentage(n,
  #     percentage,
  #     grouping_by = NULL,
  #     remove_columns = "total"
  #   )

  samples_per_ocean <- add_percentage(samples_per_ocean,
    n,
    percentage,
    grouping_by = NULL,
    remove_columns = "total"
    )

  # samples_per_ocean %>%
  #   left_join(
  #     count_and_arrange(query_df, c("ocean"), "n_ocean", "n_ocean")
  #   ) %>%
  #   arrange(desc(ocean))  
  # # print(samples_per_studyid, n = 50)

  # 299 / 429 * 100

  to_plot <- samples_per_ocean

  viridis_color_pallete <- get_viridis_colors(to_plot, hemi, "magma", -1, 0)

  (samples_per_ocean_plot <- suppressMessages(bar_plot(
    df = to_plot,
    y = ocean,
    x = n,
    fill = hemi,
    fill_pallete = viridis_color_pallete,
    y_lab = "Oceans",
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
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15))
  )


  ggsave(file.path(files_out_path, paste0("Samples_per_ocean_hemi_fill_bar", plot_ext)), height = 8.5, width = 14, units = "in", dpi = 300)


  (samples_per_studyid_ocean <- count_and_arrange(
    query_df, c("studyID", "ocean")) %>% 
  add_percentage(n,
            percentage,
            grouping_by = NULL,
            remove_columns = "total")
  )


  files_write_out_path <- "../analysis/out_files"

  write_csv(samples_per_ocean, file.path(files_write_out_path, "samples_per_ocean.csv"))

  write_csv(samples_per_studyid_ocean, file.path(files_write_out_path, "samples_per_studyid_ocean.csv"))


  ### month
  (samples_per_month <- count_and_arrange(query_df, c("month")))

  # query_df  %>%
  # filter(grepl("-12-", x = time))

  to_plot <- samples_per_month

  viridis_color_pallete <- get_viridis_colors(to_plot, month, "H", -1, 0)

  (samples_per_month_plot <- suppressMessages(bar_plot(
    df = to_plot,
    y = month,
    x = n,
    fill = NULL,
    # fill_pallete = viridis_color_pallete,
    y_lab = "Month",
    x_lab = "Number of samples",
    legend_position = "none",
    legend_direction = "vertical"
  ) +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
    scale_y_discrete(
      limits =
        rev(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")),
      labels =
      # c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
        rev(c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")),
    )))


  ggsave(file.path(files_out_path, paste0("Samples_per_month_bar", plot_ext)), height = 8.5, width = 14, units = "in", dpi = 300)

  ### month by study ID
  (samples_per_month_study_id <- count_and_arrange(query_df, c("month", "studyID")))

  # query_df  %>%
  # filter(grepl("-12-", x = time))

  to_plot <- samples_per_month_study_id

  viridis_color_pallete <- get_viridis_colors(to_plot, studyID, "H", -1, 0)

  bar_plot(
    df = to_plot,
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
      # c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
        rev(c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")),
    )


  ggsave(file.path(files_out_path, paste0("Samples_per_month_study_id_bar", plot_ext)), height = 8.5, width = 14, units = "in", dpi = 300)




  ### - Histograms

  ### * lat_abs
  samples_per_lat_abs <- count_and_arrange(query_df, c("lat_abs", "hemi"))

  (samples_per_lat_abs_plot <- suppressMessages(histogram_plot_x_or_y(samples_per_lat_abs, "x", lat_abs,
    fill = hemi, binwidth = 2,
    x_lab = "absolute latitude"
  ) +
    theme(
      legend.position = "bottom"
    ) +
    labs(
      fill = "Hemisphere"
    ) +
    scale_fill_manual(
      values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
      label = c("northern", "southern")
    )))

  ggsave(file.path(files_out_path, paste0("Samples_per_abslat_hemi_hist", plot_ext)), height = 8.5, width = 14, units = "in", dpi = 300)


  ### season
  samples_per_season <- count_and_arrange(
    query_df, c("season", "hemi")
  )

  suppressMessages((bar_plot(
    df = samples_per_season,
    x = n,
    y = season,
    fill = hemi,
    y_lab = "season",
    x_lab = "Number of samples"
  ) +
    # theme(
    #     legend.position = "bottom"
    # ) +
    labs(
      fill = "Hemisphere"
    ) +
    scale_fill_manual(
      values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
      label = c("northern", "southern")
    )))

  ggsave(file.path(files_out_path, paste0("Samples_per_season_hemi_bar", plot_ext)), height = 8.5, width = 14, units = "in", dpi = 300)


  ### depth
  samples_per_depth <- count_and_arrange(
    query_df, c("depth", "hemi")
  ) %>%
    mutate(less_than_150 = depth <= 155)

  to_plot <- samples_per_depth

  suppressMessages(bar_plot(
    df = to_plot,
    y = n,
    x = depth,
    fill = hemi,
    x_lab = "depth (m)",
    y_lab = "Number of samples",
    width = 1,
    # legend_position = "none"
  ) +
    scale_fill_manual(
      values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
      label = c("northern", "southern")
    ) +
    coord_flip() +
    scale_x_reverse() +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(
      ~less_than_150,
      ncol = 1,
      scales = "free"
    ))


  ### * just depths above 150 m
  (depths_above_150m_bar_plt <- suppressMessages(bar_plot(
    df = samples_per_depth %>%
      filter(less_than_150 == T),
    y = n,
    x = depth,
    fill = hemi,
    x_lab = "depth (m)",
    y_lab = "Number of samples",
    # legend_position = "none",
    width = 1,
    legend_position = "none",
    # legend_direction = "horizontal"
  ) +
    scale_fill_manual(
      values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
      label = c("northern", "southern")
    ) +
    # theme(
    #   axis.title.x = element_blank()
    # ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_reverse() +
    coord_flip()))

  ggsave(file.path(files_out_path, paste0("Samples_per_depth_hemi_above150m_hist", plot_ext)), height = 8.5, width = 14, units = "in", dpi = 300)

  ### * just depths below 150 m
  (depths_below_150m_bar_plt <- suppressMessages(bar_plot(
    df = samples_per_depth %>%
      filter(less_than_150 == FALSE),
    y = n,
    x = depth,
    fill = hemi,
    x_lab = "depth (m)",
    # legend_position = "none",
    width = 30,
    legend_position = "right",
    legend_direction = "vertical"
  ) +
    scale_fill_manual(
      values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
      label = c("northern", "southern")
    ) +
    # scale_x_reverse(
    #   limits = c(150, 3015),
    #   breaks = c(150, 300, 500, 750, 1000, 2000, 3000)
    # ) +
    scale_x_reverse(
      limits = c(3015, 150),
      breaks = c(150, 300, 500, 750, 1000, 2000, 3000)
    ) +
    # scale_x_continuous(
    #   limits = rev(c(150, 3015)),
    #   breaks = rev(c(150, 300, 500, 750, 1000, 2000, 3000))
    # ) +
    scale_y_continuous(
      limits = c(0, 10),
      expand = c(0, 0)
    ) +
    coord_flip()))

  ggsave(file.path(files_out_path, paste0("Samples_per_depth_hemi_below150_hist", plot_ext)), height = 8.5, width = 14, units = "in", dpi = 300)

  ### * combine plots

  combined_plot_depths <- depths_above_150m_bar_plt | depths_below_150m_bar_plt

  # print(combined_plot_depths)

  ggsave(file.path(files_out_path, paste0("Samples_per_depth_hemi_hist_cmb_plot", plot_ext)), height = 8.5, width = 14, units = "in", dpi = 300)

  ### _ oceanographic data
  ### _ oceanographic data
  ### _ oceanographic data

  ###* sst_tblSST_AVHRR_OI_NRT
  samples_per_sst_tblSST_AVHRR_OI_NRT <- count_and_arrange(
    query_df, c("sst_tblSST_AVHRR_OI_NRT", "hemi")
  )

  sst_hist <- suppressMessages(histogram_plot_x_or_y(
    df = samples_per_sst_tblSST_AVHRR_OI_NRT, aes_var = "x",
    x = sst_tblSST_AVHRR_OI_NRT,
    fill = hemi,
    binwidth = 1,
    x_lab = "SST ˚C"
  ) +
    theme(
      legend.position = "bottom"
    ) +
    labs(
      fill = "Hemisphere"
    ) +
    scale_fill_manual(
      values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
      label = c("northern", "southern")
    ))

  print(sst_hist)

  ggsave(file.path(files_out_path, paste0("Samples_per_SST_hemi_hist", plot_ext)), height = 8.5, width = 14, units = "in", dpi = 300)

  ### * PO4_tblPisces_NRT
  samples_per_PO4_tblPisces_NRT <- count_and_arrange(
    query_df, c("PO4_tblPisces_NRT", "hemi")
  )

  po4_hist <- suppressMessages(histogram_plot_x_or_y(
    df = samples_per_PO4_tblPisces_NRT, aes_var = "x",
    x = PO4_tblPisces_NRT,
    fill = hemi,
    binwidth = 0.02,
    x_lab = expression(bold(paste(
      "PO"[4]^"3-" * " (µmol L"^"-1", ")"
    )))
  ) +
    theme(
      legend.position = "bottom"
    ) +
    labs(
      fill = "Hemisphere"
    ) +
    scale_fill_manual(
      values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
      label = c("northern", "southern")
    ))

  print(po4_hist)

  ggsave(file.path(files_out_path, paste0("Samples_per_PO4_hemi_hist", plot_ext)), height = 8.5, width = 14, units = "in", dpi = 300)



  ## - Make combined plot
  #* # absolute lat ---> samples_per_lat_abs_plot
  #* # ocean  ---> samples_per_ocean_plot
  #* # SST ---> sst_hist
  #* # PO4 ---> po4_hist

  #* # alter plots for combining

  sst_hist_sub <- sst_hist + theme(legend.position = "none")
  po4_hist_sub <- po4_hist + theme(legend.position = "none")
  samples_per_lat_abs_plot_sub <- samples_per_lat_abs_plot + theme(legend.position = "none")
  samples_per_ocean_plot_sub <- samples_per_ocean_plot + theme(
    legend.position = "none",
    axis.text.y = element_text(
      angle = 45,
      hjust = 1,
      size = 19,
      face = "bold"
    )
  )

  #* # combine with patchwork

  combined_plot_fig_5 <- (
    (samples_per_lat_abs_plot_sub + samples_per_ocean_plot_sub) / 
    (sst_hist_sub + po4_hist_sub)
  ) #+
    # plot_layout(guides = "collect") &
    # theme(
    #   legend.position = "bottom",
    #   legend.box.just = "center",
    #   legend.direction = "horizontal"
    # )

  print(combined_plot_fig_5)

  ggsave(file.path(files_out_path, paste0("Samples_per_abslat_ocean_sst_po4_hemi_cmb_hist", plot_ext)), height = 8.5, width = 17, units = "in", dpi = 300)


  ### * logFe
  samples_per_logFe <- count_and_arrange(
    query_df, c("logFe", "hemi")
  )

  # print(samples_per_logFe, n = 1000)

  histogram_plot_x_or_y(
    df = samples_per_logFe, aes_var = "x",
    x = logFe,
    fill = hemi,
    binwidth = 0.1,
    x_lab = expression(bold(paste(
      "log(Fe)"
    )))
  ) +
    theme(
      legend.position = "bottom"
    ) +
    labs(
      fill = "Hemisphere"
    ) +
    scale_fill_manual(
      values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
      label = c("northern", "southern")
    )

  ggsave(file.path(files_out_path, paste0("Samples_per_logFe_hemi_hist", plot_ext)), height = 8.5, width = 14, units = "in", dpi = 300)


  ### * PP_tblPisces_NRT
  samples_per_PP_tblPisces_NRT <- count_and_arrange(
    query_df, c("PP_tblPisces_NRT", "hemi")
  )

  # print(samples_per_PP_tblPisces_NRT, n = 1000)

  histogram_plot_x_or_y(
    df = samples_per_PP_tblPisces_NRT, aes_var = "x",
    x = PP_tblPisces_NRT,
    fill = hemi,
    binwidth = 0.001,
    x_lab = expression(bold(paste(
      "PP"
    )))
  ) +
    theme(
      legend.position = "bottom"
    ) +
    labs(
      fill = "Hemisphere"
    ) +
    scale_fill_manual(
      values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
      label = c("northern", "southern")
    )

  ggsave(file.path(files_out_path, paste0("Samples_per_PP_tblPisces_NRT_hemi_hist", plot_ext)), height = 8.5, width = 14, units = "in", dpi = 300)


  ### * CHL_tblPisces_NRT
  samples_per_CHL_tblPisces_NRT <- count_and_arrange(
    query_df, c("CHL_tblPisces_NRT", "hemi")
  )

  # print(samples_per_CHL_tblPisces_NRT, n = 1000)

  histogram_plot_x_or_y(
    df = samples_per_CHL_tblPisces_NRT, aes_var = "x",
    x = CHL_tblPisces_NRT,
    fill = hemi,
    binwidth = 0.1,
    x_lab = expression(bold(paste(
      "chlorophyll a"
    )))
  ) +
    theme(
      legend.position = "bottom"
    ) +
    labs(
      fill = "Hemisphere"
    ) +
    scale_fill_manual(
      values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
      label = c("northern", "southern")
    )

  ggsave(file.path(files_out_path, paste0("Samples_per_CHL_tblPisces_NRT_hemi_hist", plot_ext)), height = 8.5, width = 14, units = "in", dpi = 300)
}


#' Main function to execute analysis
main <- function(files_to_read, files_in_path, files_out_path, plot_ext) {

  # Load the data
  data_list <- load_files(files_to_read, files_in_path)

  # Generate and save plots
  generate_plots(
    cmap_coloc = data_list$cmap_coloc,
    samples_per_nucleicAcidType = data_list$samples_per_nucleicAcidType,
    samples_per_nucleicAcidType_studyid =
      data_list$samples_per_nucleicAcidType_studyid,
    samples_per_photic = data_list$samples_per_photic,
    samples_per_photic_nucacid = data_list$samples_per_photic_nucacid,
    sample_type = data_list$sample_type,
    unique_sample_id_key = data_list$unique_sample_id_key,
    plot_ext = plot_ext,
    files_out_path = files_out_path
  )

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

  final_results <- main(
    # files_to_source,
    args$files_to_read,
    args$files_in_path,
    args$files_out_path,
    args$plot_ext
  )
}
