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
        default = "annoNifHDB_updt,cmap_coloc,RA_df_T_lng_mean_RA_AUID_deduped"
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



# temp <- RA_df_T_lng_mean_RA_AUID_deduped %>%
# #   add_group_id()  %>%
#   merge_annotations(annotation_table = annoNifHDB_updt) %>%
# #   merge_cmap(by_join = c("SAMPLEID", "group_id"))
#   merge_cmap(cmap = cmap_coloc, by_join = c("SAMPLEID"))


#' Create Scatter Plots for Relative Abundance Data
#'
#' This function generates scatter plots for relative abundance data against absolute latitude and sea surface temperature (SST).
#'
#' @param rel_abund_df Data frame containing relative abundance data.
#' @param anno_table Data frame containing annotation information.
#' @param meta_table Data frame containing metadata information.
#' @param colors Named vector of colors for plotting.
#' @param plot_ext File extension for saving plots (e.g., ".jpeg").
#' @param plot_device Plot device used for writing plots (e.g., "cairo_pdf",
#' "svg").
#' @param output_dir Directory where the plots will be saved.
#'
#' @return Logical indicating success (TRUE) or failure (FALSE).
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

    # # Setup succes flag
    # success <- TRUE
    success <- FALSE

    # tryCatch(
    #     {
            # Merge relative abundance data with annotations and metadata
            cat("Merging relative abundance data with annotations and metadata...\n")
            query_df <- rel_abund_df %>%
                left_join(anno_table, by = "AUID") %>%
                left_join(meta_table, by = "SAMPLEID")
            #   add_group_id()  %>%
            # merge_annotations(annotation_table = anno_table) %>%
                #   merge_cmap(by_join = c("SAMPLEID", "group_id"))
                # merge_cmap(cmap = meta_table, by_join = c("SAMPLEID"))
            cat("Successfully merge files and generated query_df\n")

            # Sum each sample by nifH cluster and select relevant columns
            cat("Summing each sample by nifH cluster and selecting relevant columns...\n")
            query_df_cluster <- query_df %>%
                # group_by(SAMPLEID, nifH_cluster) %>%
                # mutate(RA = sum(RA, na.rm = T)) %>%
                # select(SAMPLEID, studyID, time, depth, lat, lon, nifH_cluster, RA, CyanoCON, lat_abs, hemi, everything()) %>%
                # ungroup() %>%
                # distinct(SAMPLEID, nifH_cluster, .keep_all = T) %>%
                # ungroup()
  group_by(SAMPLEID, nifH_cluster) %>%
  mutate(RA = sum(RA, na.rm = TRUE)) %>%
  select(SAMPLEID, studyID, time, depth, lat, lon, nifH_cluster, RA, CyanoCON, lat_abs, hemi, everything()) %>%
  ungroup() %>%
  distinct(SAMPLEID, nifH_cluster, .keep_all = TRUE)

            print(query_df_cluster)
            names(query_df_cluster)
            cat("\n")

            # Filter for specific nifH clusters of interest
            cat("Filtering for specific nifH clusters of interest...\n")
            query_df_cluster_sub <- query_df_cluster %>%
                filter(nifH_cluster %in% c("1J/1K", "3", "1A", "1B", "1G", "1O/1P"))
            print(query_df_cluster_sub)
            names(query_df_cluster_sub)
            cat("\n")

            # Create scatter plot: Relative Abundance vs Absolute Latitude
            cat("Creating scatter plot: Relative Abundance vs Absolute Latitude...\n")
            ra_v_abslat_plot <- scatter_line_plot(
                df = query_df_cluster_sub,
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

            # Create a basic ggplot scatter plot
            basic_plot <- ggplot(query_df_cluster_sub, aes(x = lat_abs, y = RA, color = nifH_cluster)) +
                geom_point() +
                labs(title = "Basic Scatter Plot", x = "Absolute Latitude", y = "% Relative Abundance") +
                theme_minimal()

            # Create scatter plot: Relative Abundance vs SST
            cat("Creating scatter plot: Relative Abundance vs SST...\n")
            ## RA versus SST
            ra_v_sst_plot <- scatter_line_plot(
                df = query_df_cluster_sub,
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


            # Combine plots using patchwork
            cat("Combining plots using patchwork...\n")
            ra_v_abslat_plot_sub <- ra_v_abslat_plot +
                theme(legend.position = "none")
            ra_v_sst_plot_sub <- ra_v_sst_plot
            combined_sct_plot <- ra_v_abslat_plot_sub / ra_v_sst_plot_sub

            # Save the combined scatter plot
            cat("Saving the combined scatter plot...\n")
            cat("Unique nifH_clusters:", unique(query_df_cluster_sub$nifH_cluster), "\n")
            cat("Number of colors in palette:", length(colors), "\n")
            cat("Checking for NA or non-finite values...\n")
            print(sum(is.na(query_df_cluster_sub$RA)))
            print(sum(!is.finite(query_df_cluster_sub$RA)))
            print(sum(is.na(query_df_cluster_sub$lat_abs)))
            print(sum(!is.finite(query_df_cluster_sub$lat_abs)))


            # Setup succes flag
            success <- TRUE
    #     },
    #     warning = function(w) {
    #         cat("A warning occured \n")
    #         cat("Warning messege:", conditionMessage(w), "\n")
    #         cat("Warning call:", deparse(conditionCall(w)), "\n")
    #     },
    #     error = function(e) {
    #         cat("An error occurred running make_scatter()")
    #         cat("Error message:", conditionMessage(e), "\n")
    #         cat("Error call:", deparse(conditionCall(e)), "\n")
    #         success <<- FALSE
    #     }
    # )


    return(success)
}



#' Main function to execute plots
main <- function(files_to_source, files_to_read, files_in_path, files_out_path, plot_ext, plot_device, nifh_cluster_colours) {
    # Load the data
    data_list <- load_files(files_to_read, files_in_path)

    # # Generate color palette
    # nifh_cluster_colours <- setNames(RColorBrewer::brewer.pal(8, 'Paired')[c(8, 2, 11, 9, 5, 4)], c("1J/1K", "3", "1A", "1B", "1G", "1O/1P"))

    # Generate the color blind safe paletter for nifH clusters
    cluster_palette <- generate_nifh_palette(nifh_cluster_colours = nifh_cluster_colours)

    print(cluster_palette)

    # Create sample data
    set.seed(123) # For reproducibility
    sample_data <- data.frame(
        x = rnorm(100), # 100 random normal values for x-axis
        y = rnorm(100) # 100 random normal values for y-axis
    )

    # Define a simple color palette
    simple_colors <- c("blue", "red")

    # Create a scatter plot
    scatter_plot <- ggplot(sample_data, aes(x = x, y = y)) +
        geom_point(color = simple_colors[1]) + # Use the first color for points
        labs(
            title = "Simple Scatter Plot",
            x = "X-axis label",
            y = "Y-axis label"
        ) +
        theme_minimal()

    # Print the plot to the console
    print(scatter_plot)

    # Generate and save plots
    success_scatter <- make_scatter(
        rel_abund_df = data_list$RA_df_T_lng_mean_RA_AUID_deduped,
        anno_table = data_list$annoNifHDB_updt,
        meta_table = data_list$cmap_coloc,
        output_dir = files_out_path,
        colors = cluster_palette,
        plot_ext = plot_ext,
        plot_device = plot_device
    )


    # if (success_pie && success_bar && success_scatter) {
    if (success_scatter) {
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
        plot_device = args$plot_device,
        nifh_cluster_colours = nifh_cluster_colours
    )
}
