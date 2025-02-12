#!/usr/bin/env Rscript

#' @title NifH Cluster Analysis Pipeline Sample and Taxonomy Breakdown
#' @description This module processes and analyzes nifH cluster data from amplicon sequencing, generating comprehensive count tables and relative abundance metrics for both DNA and RNA samples. It performs various statistical analyses, data transformations, and produces cleaned and threshold-filtered results for downstream analysis.
#'
#' @details The pipeline executes the following main steps:
#' * Loads input files containing raw count data, relative abundance data, and metadata
#' * Processes raw count and abundance data, removing RNA samples and unknown AUIDs
#' * Generates cluster-specific statistics for both pooled and study-specific analyses
#' * Calculates percentages by nifH cluster and study ID
#' * Cleans up tables by grouping low-abundance clusters into an 'other' category
#' * Produces both count-based and relative abundance metrics
#'
#' Key functions include:
#' * count_tables(): Generates count tables for nifH clusters
#' * rel_abund_tables(): Generates relative abundance tables for nifH clusters
#' * main(): Orchestrates the entire analysis workflow
#'
#' @usage Rscript samp_n_tax_breakdown.R [--files FILES] [--input_path PATH] [--output_path PATH]
#'
#' @param --files Comma-separated list of input files (default: counts_df_T_lng_AUID_deduped,RA_df_T_lng_mean_RA_AUID_deduped,cmap_coloc,annoNifHDB_updt,unique_sample_id_key,photic_samples_key,DNA_samples_key)
#' @param --input_path Directory path for input files (default: ../analysis/out_files)
#' @param --output_path Directory path for output files (default: ../analysis/out_files)
#'
#' @author Michael Morando
#'  2025-02-02
#'
#' @note This script requires the following R packages: tidyverse, patchwork, argparser
#'
#' @examples
#' Rscript samp_n_tax_breakdown.R --files counts_df_T_lng_AUID_deduped,RA_df_T_lng_mean_RA_AUID_deduped,cmap_coloc,annoNifHDB_updt,unique_sample_id_key,photic_samples_key,DNA_samples_key --input_path ../data/processed --output_path ../results/analysis
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
  parser <- arg_parser("Process the workspace files")
  parser <- add_argument(parser, "--files",
    help = "CVS list of files to read in",
    default = "counts_df_T_lng_AUID_deduped,RA_df_T_lng_mean_RA_AUID_deduped,cmap_coloc,annoNifHDB_updt,unique_sample_id_key,photic_samples_key,DNA_samples_key"
  )
  parser <- add_argument(parser,
    arg = "--input_path",
    help = "Input directory path",
    default = "../analysis/out_files"
  )
  parser <- add_argument(parser,
    arg = "--output_path",
    help = "Output directory path",
    default = "../analysis/out_files"
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
    files_out_path = argv$output_path
  ))
}


#### - count data

#' Generate count tables for nifH clusters
#'
#' This function generates count tables for nifH clusters based on the provided abundance table and annotation table.
#' It calculates percentages by nifH cluster and cleans up the tables by grouping percentages below a threshold into an 'other' category.
#'
#' @param abundance_table A data frame containing the abundance data with columns: SAMPLEID, AUID, and counts.
#' @param annotation_table A data frame containing the annotation information for nifH clusters.
#' @param threshold A numeric value specifying the threshold for grouping percentages into the 'other' category (default: 1.0).
#'
#' @return A list containing the following data frames:
#'   - nifhdb_all_counts_AUID_dedup_clean: Cleaned count table with percentages by nifH cluster.
#'   - nifhdb_all_counts_AUID_dedup_total_study_id_clean: Cleaned count table with percentages by nifH cluster and study ID (total counts).
#'   - nifhdb_all_counts_AUID_dedup_study_id_total_clean: Cleaned count table with percentages by nifH cluster and study ID (study-specific counts).
#'
#' @examples
#' count_tables(abundance_table, annotation_table, threshold = 1.0)
#'
count_tables <- function(abundance_table, annotation_table, metatable, DNA_samples_key, threshold = 1.0) {
  tryCatch(
    {
      ### * make query data frame where you remove RNA samples
      ## * add studyID column
      cat("Generating query data frame...\n")
      query_df <- abundance_table %>%
        replace_na(list(. = 0)) %>%
        remove_samples_nucleic_acid(
          nucleic_acid_type = "DNA",
          nucleic_acid_samples_key = DNA_samples_key
        ) %>%
        # remove_unknown_auids() %>%
        left_join(metatable %>%
          select(SAMPLEID, studyID))

      num_samples <- query_df %>%
        distinct(SAMPLEID) %>%
        nrow()
      cat("Number of distinct samples in query data frame:", num_samples, "\n")


      cat("Calculating percentages by nifH cluster...\n")
      (nifhdb_all_counts_AUID_dedup <- sum_tax(
        abundance_table = query_df,
        annotation_table = annotation_table,
        joining_by_col = "AUID",
        grouping_by_1 = nifH_cluster,
        grouping_by_2 = NULL,
        sum_by = counts,
        total_counts_id = total_counts_nifH_cluster,
        grouping_for_percentage = NULL,
        #   sum_by_percent = total_counts,
        percentage_id = percentage_nifH_cluster,
        nifH_cluster
      ))


      cat("Cleaning up table by grouping percentages below the threshold '", threshold, "' into group 'other'...\n")
      (nifhdb_all_counts_AUID_dedup_clean <- clean_percentages(
        df = nifhdb_all_counts_AUID_dedup,
        grouping_by = NULL,
        clean_var = percentage_nifH_cluster,
        threshold = 1.0,
        # sum_by = percentage_nifH_cluster_study_id,
        percentage_id = percentage_nifH_cluster,
        column_var_mutate = nifH_cluster_modified,
        column_var = nifH_cluster
      ))

      #### - By study ID

      cat("Calculating percentages by nifH cluster and study ID (total counts)...\n")
      (nifhdb_all_counts_AUID_dedup_total_study_id <- sum_tax(
        abundance_table = query_df,
        annotation_table = annotation_table,
        joining_by_col = "AUID",
        sum_by = counts,
        grouping_by_1 = nifH_cluster,
        grouping_by_2 = studyID,
        grouping_for_percentage = NULL, ### make this NULL to get percentages of total instead of percentages by studyID
        total_counts_id = total_counts_nifH_cluster_total_study_id,
        #   sum_by_percent = total_counts,
        percentage_id = percentage_total_counts_nifH_cluster_total_study_id
      ))


      cat("Cleaning up table by grouping percentages threshold '", threshold, "' into group 'other' (total counts)...\n")
      (nifhdb_all_counts_AUID_dedup_total_study_id_clean <- clean_percentages(
        df = nifhdb_all_counts_AUID_dedup_total_study_id,
        grouping_by = studyID,
        clean_var = percentage_total_counts_nifH_cluster_total_study_id,
        # threshold = 2.5,
        threshold = 1.0,
        # sum_by = percentage_nifH_cluster_study_id,
        percentage_id = percentage_total_counts_nifH_cluster_total_study_id_clean,
        column_var_mutate = nifH_cluster_modified,
        column_var = nifH_cluster
      ))
      ### ! FIXME: the percentages are being calculated by totals and not totals within each study ID


      cat("Calculating percentages by nifH cluster and study ID (study-specific counts)...\n")
      (nifhdb_all_counts_AUID_dedup_study_id_total <- sum_tax(
        abundance_table = query_df,
        annotation_table = annotation_table,
        joining_by_col = "AUID",
        sum_by = counts,
        grouping_by_1 = nifH_cluster,
        grouping_by_2 = studyID,
        grouping_for_percentage = studyID, ### make this NULL to get percentages of total instead of percentages by studyID
        total_counts_id = total_counts_nifH_cluster_study_id_total,
        #   sum_by_percent = total_counts,
        percentage_id = percentage_total_counts_nifH_cluster_study_id_total,
        nifH_cluster, studyID
      ))


      cat("Cleaning up table by grouping percentages below the threshold into 'other' (study-specific counts)...\n")
      (nifhdb_all_counts_AUID_dedup_study_id_total_clean <- clean_percentages(
        df = nifhdb_all_counts_AUID_dedup_study_id_total,
        grouping_by = studyID,
        clean_var = percentage_total_counts_nifH_cluster_study_id_total,
        # threshold = 2.5,
        threshold = 1.0,
        # sum_by = percentage_nifH_cluster_study_id,
        percentage_id = percentage_total_counts_nifH_cluster_study_id_total_clean,
        column_var_mutate = nifH_cluster_modified,
        column_var = nifH_cluster
      ))


      cat("Returning the cleaned count tables as a list...\n")
      # Return the cleaned count tables
      return(list(
        nifhdb_all_counts_AUID_dedup_study_id_total = nifhdb_all_counts_AUID_dedup_study_id_total,
        nifhdb_all_counts_AUID_dedup_clean = nifhdb_all_counts_AUID_dedup_clean,
        nifhdb_all_counts_AUID_dedup_total_study_id_clean = nifhdb_all_counts_AUID_dedup_total_study_id_clean,
        nifhdb_all_counts_AUID_dedup_study_id_total_clean = nifhdb_all_counts_AUID_dedup_study_id_total_clean
      ))
    },
    # warning = function(w) {
    #     cat("Warning occurred:\n")
    #     cat("Warning message:", conditionMessage(w), "\n")
    #     cat("Warning call:", deparse(conditionCall(w)), "\n")
    #   },
    error = function(e) {
      cat("An error occurred in running count_tables():\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Error call:", deparse(conditionCall(e)), "\n")
      return(NULL)
    }
  )
}




#### - relative abundance data

#' Generate relative abundance tables for nifH clusters
#'
#' This function generates relative abundance tables for nifH clusters based on the provided abundance table and annotation table.
#' It calculates percentages by nifH cluster and cleans up the tables by grouping percentages below a threshold into an 'other' category.
#'
#' @param abundance_table A data frame containing the abundance data with columns: SAMPLEID, AUID, and RA.
#' @param annotation_table A data frame containing the annotation information for nifH clusters.
#' @param threshold A numeric value specifying the threshold for grouping percentages into the 'other' category (default: 1.0).
#'
#' @return A list containing the following data frames:
#'   - nifhdb_all_rel_abund_AUID_dedup_clean: Cleaned relative abundance table with percentages by nifH cluster.
#'   - nifhdb_all_rel_abund_AUID_dedup_total_study_id_clean: Cleaned relative abundance table with percentages by nifH cluster and study ID (total relative abundance).
#'   - nifhdb_all_rel_abund_AUID_dedup_study_id_total_clean: Cleaned relative abundance table with percentages by nifH cluster and study ID (study-specific relative abundance).
#'
#' @examples
#' relative_abundance_tables(abundance_table, annotation_table, threshold = 1.0)
#'
rel_abund_tables <- function(abundance_table, annotation_table, metatable, DNA_samples_key, threshold = 1.0) {
  tryCatch(
    {
      ### * make query data frame where you remove RNA samples
      ## * add studyID column
      query_df <- abundance_table %>%
        replace_na(list(. = 0)) %>%
        remove_samples_nucleic_acid(
          nucleic_acid_type = "DNA",
          nucleic_acid_samples_key = DNA_samples_key
        ) %>%
        # remove_unknown_auids() %>%
        left_join(metatable %>%
          select(SAMPLEID, studyID))

      ## * calculate percentages by nifH cluster
      (nifhdb_all_rel_abund_AUID_dedup <- sum_tax(
        abundance_table = query_df,
        annotation_table = annotation_table,
        joining_by_col = "AUID",
        grouping_by_1 = nifH_cluster,
        grouping_by_2 = NULL,
        sum_by = average_value_AUID,
        total_counts_id = total_rel_abund_nifH_cluster,
        grouping_for_percentage = NULL,
        #   sum_by_percent = total_counts,
        percentage_id = percentage_nifH_cluster,
        nifH_cluster
      ))


      ## * Clean up table so that everything below given threshold is summed together as 'other'
      (nifhdb_all_rel_abund_AUID_dedup_clean <-
        clean_percentages(
          df = nifhdb_all_rel_abund_AUID_dedup,
          grouping_by = NULL,
          clean_var = percentage_nifH_cluster,
          threshold = 1.0,
          # sum_by = percentage_nifH_cluster_study_id,
          percentage_id = percentage_nifH_cluster,
          column_var_mutate = nifH_cluster_modified,
          column_var = nifH_cluster
        ))

      #### - By study ID

      ## * calculate over pooled data, so each study ID is a percentage of the total
      (nifhdb_all_rel_abund_AUID_dedup_total_study_id <- sum_tax(
        abundance_table = query_df,
        annotation_table = annotation_table,
        joining_by_col = "AUID",
        grouping_by_1 = nifH_cluster,
        grouping_by_2 = studyID,
        sum_by = average_value_AUID,
        grouping_for_percentage = NULL, ### make this NULL to get percentages of total instead of percentages by studyID
        total_counts_id = total_rel_abund_total_nifH_cluster,
        #   sum_by_percent = total_counts,
        percentage_id = percentage_total_rel_abund_total_nifH_cluster_study_id,
        nifH_cluster, studyID
      ))


      ## * Clean up table so that everything below a certain threshold is summed together as 'other'
      nifhdb_all_rel_abund_AUID_dedup_total_study_id_clean <- suppressMessages(nifhdb_all_rel_abund_AUID_dedup_total_study_id %>%
        clean_percentages(
          grouping_by = studyID,
          clean_var = percentage_total_rel_abund_total_nifH_cluster_study_id,
          # threshold = 2.5,
          threshold = 1.0,
          # sum_by = percentage_nifH_cluster_study_id,
          percentage_id = percentage_total_rel_abund_total_nifH_cluster_study_id_clean,
          column_var_mutate = nifH_cluster_modified,
          column_var = nifH_cluster
        ))


      ## * calculate over each study ID, so each percentage is based on its own study.
      (nifhdb_all_rel_abund_AUID_dedup_study_id_total <- sum_tax(
        abundance_table = query_df,
        annotation_table = annotation_table,
        joining_by_col = "AUID",
        sum_by = average_value_AUID,
        grouping_by_1 = nifH_cluster,
        grouping_by_2 = studyID,
        grouping_for_percentage = studyID, ### make this NULL to get percentages of total instead of percentages by studyID
        total_counts_id = total_rel_abund_nifH_cluster_study_id_total,
        #   sum_by_percent = total_counts,
        percentage_id = percentage_total_rel_abund_nifH_cluster_study_id_total,
        nifH_cluster, studyID
      ))


      ## * Clean up table so that everything below a certain threshold is summed together as 'other'
      (nifhdb_all_rel_abund_AUID_dedup_study_id_total_clean <- suppressMessages(clean_percentages(
        df = nifhdb_all_rel_abund_AUID_dedup_study_id_total,
        grouping_by = studyID,
        clean_var = percentage_total_rel_abund_nifH_cluster_study_id_total,
        # threshold = 2.5,
        threshold = 1.0,
        # sum_by = percentage_nifH_cluster_study_id,
        percentage_id = percentage_total_rel_abund_nifH_cluster_study_id_total_clean,
        column_var_mutate = nifH_cluster_modified,
        column_var = nifH_cluster
      )))

      cat("Returning the cleaned relative abundance tables as a list...\n")
      # Return the cleaned count tables
      list(
        nifhdb_all_rel_abund_AUID_dedup_study_id_total = nifhdb_all_rel_abund_AUID_dedup_study_id_total,
        nifhdb_all_rel_abund_AUID_dedup_clean = nifhdb_all_rel_abund_AUID_dedup_clean,
        nifhdb_all_rel_abund_AUID_dedup_total_study_id_clean = nifhdb_all_rel_abund_AUID_dedup_total_study_id_clean,
        nifhdb_all_rel_abund_AUID_dedup_study_id_total_clean = nifhdb_all_rel_abund_AUID_dedup_study_id_total_clean
      )
    },
    error = function(e) {
      cat("An error occurred in running rel_abund_tables():\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Error call:", deparse(conditionCall(e)), "\n")
      return(NULL)
    }
  )
}


### ! \TODO: Starting to develop AUID stats
## ! needed for Ecology manuscript analysis


# (query_df <- abundance_table %>%
#   replace_na(list(. = 0)) %>%
#   remove_samples_nucleic_acid(nucleic_acid_type = "DNA",
#   nucleic_acid_samples_key = DNA_samples_key) %>%
#   remove_aphotic_samples() %>%
#   # remove_unknown_auids() %>%
#   # filter(AUID %in% unknownnan_auids) %>%
#   left_join(metatable %>%
#     select(SAMPLEID, studyID)) %>%
#   arrange(desc(counts))
# )

# ## * calculate percentages by AUID
# (temp <- sum_tax(
#   abundance_table = query_df,
#   annotation_table = annotation_table,
#   joining_by_col = "AUID",
#   grouping_by_1 = AUID,
#   grouping_by_2 = NULL,
#   sum_by = counts,
#   total_counts_id = total_counts_AUID,
#   grouping_for_percentage = NULL,
#   #   sum_by_percent = total_counts,
#   percentage_id = percentage_AUID,
#   nifH_cluster
# ))

## lets look a bit at the unknown
# temp_unknown <- temp %>%
#   filter_df(plt_flt = nifH_cluster == "unknown")


main <- function(files_to_source, files_to_read, files_in_path, files_out_path) {
  #  Load other scripts and files

  # Load data
  data_list <- load_files(files_to_read, files_in_path)

  # Make tables
  count_table_list <- count_tables(
    abundance_table = data_list$counts_df_T_lng_AUID_deduped,
    annotation_table = data_list$annoNifHDB_updt,
    metatable = data_list$cmap_coloc,
    DNA_samples_key = data_list$DNA_samples_key,
    threshold = 1.0
  )


  rel_aund_table_list <- rel_abund_tables(
    abundance_table = data_list$RA_df_T_lng_mean_RA_AUID_deduped,
    annotation_table = data_list$annoNifHDB_updt,
    metatable = data_list$cmap_coloc,
    DNA_samples_key = data_list$DNA_samples_key,
    threshold = 1.0
  )

  # Concatenate the two lists into a single list
  final_results <- c(count_table_list, rel_aund_table_list)

  if (!is.null(final_results)) {
    # Create output directory if it doesn't exist
    create_dir(files_out_path)
    write_file_list(final_results, files_out_path)
  }
}


# Execute script if being called from command-line
if (sys.nframe() == 0 && !interactive()) {
  source_file(files_to_source)

  parser <- setup_parser()
  args <- parse_arg(parser)

  validate_parsed_args(parsed_args = args)

  # Disable summarise grouped output messages globally
  options(dplyr.summarise.inform = FALSE)

  main(
    files_to_source = files_to_source,
    files_to_read = args$files_to_read,
    files_in_path = args$files_in_path,
    files_out_path = args$files_out_path
  )
}
