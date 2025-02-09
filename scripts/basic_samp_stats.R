#!/usr/bin/env Rscript

#' @title NIFH Database Analysis Pipeline Basic Sampling Stats
#' @description This script performs comprehensive analysis on NIFH database samples and read counts, providing detailed statistics about sample distributions, read counts, and sample characteristics. It processes input data, calculates various metrics, and outputs results for further analysis.
#'
#' @details The pipeline executes the following main steps:
#' * Loads and validates input files from specified paths
#' * Counts total reads in the database (raw and deduplicated)
#' * Calculates sample numbers across different categories (e.g., photic/aphotic, DNA/RNA)
#' * Analyzes sample statistics by various groupings (e.g., nucleic acid type, study ID)
#' * Generates and outputs result tables for further analysis
#'
#' Key functions include:
#' * count_total_reads(): Calculates total read counts for various sample subsets
#' * calculate_sample_numbers(): Computes sample counts for different categories
#' * analyze_sample_statistics(): Performs detailed statistical analysis on samples
#' * main(): Orchestrates the entire analysis workflow
#'
#' @usage Rscript basic_samp_stats.R [--files FILES] [--input_path PATH] [--output_path PATH]
#'
#' @param --files Comma-separated list of input files (default: counts_df_T_lng,counts_df_T_lng_AUID_deduped,cmap_coloc,unique_sample_id_key,photic_samples_key,DNA_samples_key)
#' @param --input_path Directory path for input files (default: ../analysis/out_files)
#' @param --output_path Directory path for output files (default: ../analysis/out_files)
#'
#' @author Michael Morando
#' @date 2025-02-02
#'
#' @note This script requires the following R packages: tidyverse, argparser
#'
#' @examples
#' Rscript basic_samp_stats.R --files counts_df_T_lng,counts_df_T_lng_AUID_deduped,cmap_coloc,unique_sample_id_key,photic_samples_key,DNA_samples_key --input_path ../data/processed --output_path ../results/analysis
#'
#' @export


# Load required libraries with error handling
suppressPackageStartupMessages({
  tryCatch(
    {
      library(tidyverse)
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
          } # ,
          # warning = function(w) {
          #   cat("Warning while sourcing : ", file_path, ":", conditionMessage(w), "\n")
          # }
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
    default = "counts_df_T_lng,counts_df_T_lng_AUID_deduped,cmap_coloc,unique_sample_id_key,photic_samples_key,DNA_samples_key"
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
#'    \item{files_in_pathd}{The input directory path}
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


#' Count total reads in the nifH database
#'
#' @param counts_df_T_lng Dataframe with read counts
#' @param counts_df_T_lng_AUID_deduped Dataframe with deduplicated read counts
#' @param photic_samples_key Key for photic samples
#' @importFrom dplyr select
count_total_reads <- function(counts_df_T_lng, counts_df_T_lng_AUID_deduped, photic_samples_key) {
  tryCatch(
    {
      if (nrow(counts_df_T_lng) == 0 || is.null(counts_df_T_lng)) {
        stop("The provided tibble: '", deparse(substitute(counts_df_T_lng)), "' is empty\n")
      } else {
        # cat(photic_samples_key)
        # break()
        # Total reads in nifh database
        tot_cnts <- counts_df_T_lng %>%
          select(counts) %>%
          sum()
        cat("Total reads in nifh database:", tot_cnts, "\n")

        # Total reads in nifH database for photic samples
        tot_cnts_photic <- counts_df_T_lng %>%
          remove_aphotic_samples(photic_samples_key) %>%
          select(counts) %>%
          sum()
        cat("Total reads in nifH database, just photic samples:", tot_cnts_photic, "\n")
      }

      if (nrow(counts_df_T_lng_AUID_deduped) == 0 || is.null(counts_df_T_lng_AUID_deduped)) {
        stop("The provided tibble: '", deparse(substitute(counts_df_T_lng_AUID_deduped)), "' is empty\n")
      } else {
        # Count total deduplicated reads in nifH database
        tot_cnts_deduped <- counts_df_T_lng_AUID_deduped %>%
          select(counts) %>%
          sum()
        cat("Total reads in nifH database, deduplicated:", tot_cnts_deduped, "\n")

        # Count total deduplicated reads in nifH database for photic samples
        tot_cnts_deduped_photic <- counts_df_T_lng_AUID_deduped %>%
          remove_aphotic_samples(photic_samples_key) %>%
          select(counts) %>%
          sum()
        cat("Total reads in nifH database, deduplicated, just photic samples:", tot_cnts_deduped_photic, "\n")
      }
    },
    error = function(e) {
      cat("Error:", conditionMessage(e), "\n")
    } # ,
    # warning = function(w) {
    #   cat("Warning:", conditionMessage(w), "\n")
    # }
  )
}


#' Calculate sample numbers in the nifH database
#'
#' @param cmap_coloc Dataframe with sample information
#' @param DNA_samples_key Key for DNA samples
#' @param unique_sample_id_key Key for unique sample IDs
#' @param photic_samples_key Key for photic samples
#' @importFrom dplyr filter count pull
calculate_sample_numbers <- function(cmap_coloc, DNA_samples_key, unique_sample_id_key, photic_samples_key) {
  tryCatch(
    {
      ## -# photic samples
      number_of_aphotic_samples <- nrow(cmap_coloc) - nrow(remove_aphotic_samples(cmap_coloc, photic_samples_key))

      number_of_photic_samples <- nrow(cmap_coloc) -
        number_of_aphotic_samples

      cat("Total number of samples: ", nrow(cmap_coloc), "\nNumber of aphotic samples:", number_of_aphotic_samples, "\nNumber of photic samples:", number_of_photic_samples, "\n\n",
        sep = " "
      )
    },
    error = function(e) {
      cat("Error calculating sample numbers:", conditionMessage(e), "\n")
    } # ,
    # warning = function(w) {
    #   cat("Warning encountered while calculating sample numbers:", conditionMessage(w), "\n")
    # }
  )

  tryCatch(
    {
      ## -# DNA/RNA samples
      number_of_DNA_samples <- nrow(cmap_coloc) - nrow(remove_samples_nucleic_acid(cmap_coloc, "RNA", DNA_samples_key))

      number_of_RNA_samples <- nrow(cmap_coloc) - nrow(remove_samples_nucleic_acid(cmap_coloc, "DNA", DNA_samples_key))

      cat("Total number of samples: ", nrow(cmap_coloc), "\nTotal number of DNA samples:", number_of_DNA_samples, "\nTotal number of RNA samples:", number_of_RNA_samples, "\n\n",
        sep = " "
      )
    },
    error = function(e) {
      cat("Error calculating DNA/RNA samples:", conditionMessage(e), "\n")
    } # ,
    # warning = function(w) {
    #   cat("Warning encountered while calculating DNA/RNA samples:", conditionMessage(w), "\n")
    # }
  )

  tryCatch(
    {
      ## -# replicate samples
      number_of_replicate_samples <- unique_sample_id_key %>%
        filter(replicate_flag %in% "Replicate_Sample") %>%
        count() %>%
        pull()

      cat("Number of replicate samples:", number_of_replicate_samples, "\n\n",
        sep = " "
      )
    },
    error = function(e) {
      cat("Error calculating replicate samples:", conditionMessage(e), "\n")
    } # ,
    # warning = function(w) {
    #   cat("Warning encountered while calculating replicate samples:", conditionMessage(w), "\n")
    # }
  )

  tryCatch(
    {
      #-# Size_fraction
      nsf <- count_and_arrange(unique_sample_id_key, c("size_frac_flag"))

      one_sf <- paste(nsf %>% filter(size_frac_flag == "One_Size_Fractions"), collapse = ": ", sep = ":")
      two_sf <- paste(nsf %>% filter(size_frac_flag == "Two_Size_Fractions"), collapse = ": ", sep = ":")

      cat("Number of", one_sf, "\n")
      cat("Number of", two_sf, "\n\n")
    },
    error = function(e) {
      cat("Error calculating size fractions:", conditionMessage(e), "\n")
    } # ,
    # warning = function(w) {
    #   cat("Warning encountered while calculating size fractions:", conditionMessage(w), "\n")
    # }
  )
}


#' Analyze sample statistics in the NIFH database
#'
#' @param cmap_coloc Dataframe with sample information
#' @param unique_sample_id_key Key for unique sample IDs
#' @importFrom dplyr group_by summarize mutate select
analyze_sample_statistics <- function(cmap_coloc, unique_sample_id_key) {
  cat("Analyzing sample stats...\n")
  tryCatch(
    {
      tryCatch(
        {
          query_df <- dedup_by_group(cmap_coloc, group_id_key = unique_sample_id_key, group_id)
        },
        error = function(e) {
          cat("Error making query_df:", conditionMessage(e), "\n")
        }
      )

      ## - nucleicAcidType
      samples_per_nucleicAcidType <- count_and_arrange(query_df, "nucleicAcidType") %>%
        add_percentage(n,
          percentage,
          grouping_by = NULL,
          remove_columns = "total"
        ) %>%
        add_total_row(
          column_name = "nucleicAcidType",
          summary_column = n,
          # pull_name = n,
          all_columns = TRUE
        )

      print(samples_per_nucleicAcidType, n = 1000)

      ## * by study id
      samples_per_nucleicAcidType_studyid <- query_df %>%
        count_and_arrange(c("studyID", "nucleicAcidType")) %>%
        add_percentage(n,
          percentage,
          grouping_by = NULL,
          remove_columns = "total"
        ) %>%
        add_total_row(
          column_name = "studyID",
          summary_column = n,
          all_columns = FALSE
        )

      print(samples_per_nucleicAcidType_studyid, n = 50)
    },
    error = function(e) {
      cat("Error analyzing nucleic acid types by study id:", conditionMessage(e), "\n")
    }
    # warning = function(w) {
    #   cat("Warning encountered while analyzing nucleic acid types by study id:", conditionMessage(w), "\n")
    # }
  )

  tryCatch(
    {
      query_df <- dedup_by_group(cmap_coloc, unique_sample_id_key, group_id)

      ### -  photic
      samples_per_photic <- count_and_arrange(query_df, "photic") %>%
        add_percentage(n,
          percentage,
          grouping_by = NULL,
          remove_columns = "total"
        ) # %>%
      # add_total_row(
      #   column_name = "photic",
      #   summary_column = n,
      #   all_columns = FALSE
      # )

      print(samples_per_photic, n = 1000)
    },
    error = function(e) {
      cat("Error analyzing photic samples:", conditionMessage(e), "\n")
    } # ,
    # warning = function(w) {
    #   cat("Warning encountered while analyzing photic samples:", conditionMessage(w), "\n")
    # }
  )

  tryCatch(
    {
      samples_per_photic_nucacid <- count_and_arrange(cmap_coloc, c("photic", "nucleicAcidType"))


      print(samples_per_photic_nucacid, n = 100)
    },
    error = function(e) {
      cat("Error analyzing photic samples and nucleic acid types:", conditionMessage(e), "\n")
    } # ,
    # warning = function(w) {
    #   cat("Warning encountered while analyzing photic samples and nucleic acid types:", conditionMessage(w), "\n")
    # }
  )

  tryCatch(
    {
      ## - sample type
      (sample_type <- count_and_arrange(unique_sample_id_key, c("replicate_flag", "nucleicAcidType", "size_frac_flag"), replicate_flag) %>%
        mutate(
          new_group = paste(replicate_flag, nucleicAcidType, size_frac_flag, sep = "_"),
          tibble_id = "sample_id"
        ) %>%
        add_percentage(n,
          percentage,
          grouping_by = NULL,
          remove_columns = "total"
        ))

      print(sample_type, n = 100)
    },
    error = function(e) {
      cat("Error analyzing sample types:", conditionMessage(e), "\n")
    },
    # warning = function(w) {
    #   cat("Warning encountered while analyzing sample types:", conditionMessage(w), "\n")
    # },
    error = function(e) {
      cat("Error analyzing sample statistics:", conditionMessage(e), "\n")
    }
  )

  return(list(
    samples_per_nucleicAcidType = samples_per_nucleicAcidType,
    samples_per_nucleicAcidType_studyid = samples_per_nucleicAcidType_studyid,
    samples_per_photic = samples_per_photic,
    samples_per_photic_nucacid = samples_per_photic_nucacid,
    sample_type = sample_type
  ))
}



#' Main function to execute analysis
main <- function(files_to_read, files_in_path, files_out_path) {

  # Load the data
  data_list <- load_files(files_to_read, files_in_path)

  # Run the full script
  count_total_reads(
    data_list$counts_df_T_lng,
    data_list$counts_df_T_lng_AUID_deduped,
    data_list$photic_samples_key
  )
  calculate_sample_numbers(
    data_list$cmap_coloc,
    data_list$DNA_samples_key,
    data_list$unique_sample_id_key,
    data_list$photic_samples_key
  )
  stats_list <- analyze_sample_statistics(
    data_list$cmap_coloc,
    data_list$unique_sample_id_key
  )

  # return(stats_list)
  if (!is.null(stats_list)) {
    # Create output directory if it doesn't exist
    create_dir(args$files_out_path)
    write_file_list(stats_list, args$files_out_path)
  }
}


# Run if the script is being run directly
if (sys.nframe() == 0 && !interactive()) {
  source_file(files_to_source)

  parser <- setup_parser()
  args <- parse_arg(parser)

  validate_parsed_args(parsed_args = args)

  final_results <- main(
    args$files_to_read,
    args$files_in_path,
    args$files_out_path
  )
}
