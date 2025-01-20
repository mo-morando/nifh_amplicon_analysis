#!/usr/bin/env Rscript

#' Summarizing and Formatting Tables
#'
#' This script processes and summarizes data from the DADA2 nifH pipeline and
#' workflow stages.
#' It performs the following steps:
#' 1. Loads required libraries
#' 2. Sources necessary files
#' 3. Sets up argument parsing
#' 4. Reads and processes input CSV files
#' 5. Calculates percentage of retained read pairs
#' 6. Summarizes and formats tables
#' 7. Renames columns for better readability
#' 8. Writes output CSV files
#'
#' @author Michael (Mo) Morando
#' @date 2024-12-04

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
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




# # Source needed files
# files_to_source <- c(
#   "/Users/mo/Projects/nifH_amp_project/myWork/scripts/functions.R",
#   "/Users/mo/Projects/nifH_amp_project/myWork/scripts/basic_plotting.R"
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
  parser <- arg_parser("Process the workspace files")
  parser <- add_argument(parser, "--files",
    help = "CVS list of files to read in",
    default = "Table_SreadsAtEachStage_samples,workflowTable"
  )
  parser <- add_argument(parser,
    arg = "--input_path",
    help = "Input directory path",
    default = "../data/csvs_for_tables"
  )
  parser <- add_argument(parser,
    arg = "--output_path",
    help = "Output directory path",
    default = "../analysis/tables"
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
    files_to_read = strsplit(argv$files, ",")[[1]],
    files_in_path = argv$input_path,
    files_out_path = argv$output_path
  ))
}


#' Summarize workflow stages and produce summary tables
#'
#' @param df Data frame containing workflow stage data
#' @param grp_by Column to group by (unquoted)
#' @param rnd_to Number of decimal places to round to
#' @param initial_col Initial column for percentage calculation (unquoted)
#' @param final_col Final column for percentage calculation (unquoted)
#'
#' @return A summarized data frame
#-# function to summarize these files and produce Tables 4 and 5
summarise_workflow_stages_table <- function(
    df,
    grp_by,
    rnd_to,
    initial_col,
    final_col) {
  message("Summarising workflow stage: ", deparse(substitute(df)))

  summarised_table <- df %>%
    group_by({{ grp_by }}) %>%
    reframe(across(where(is.numeric), ~ mean(.))) %>%
    ungroup() %>%
    bind_rows(
      tibble(
        {{ grp_by }} := "mean",
        df %>%
          summarise(across(where(is.numeric), ~ mean(.))) #* This gives the mean over all the reads, not by study ID
      )
    ) %>%
    bind_rows(
      tibble(
        {{ grp_by }} := "median",
        df %>%
          summarise(across(where(is.numeric), ~ median(.))) #* This gives the median over all the reads, not by study ID
      )
    ) %>%
    bind_rows(
      tibble(
        {{ grp_by }} := "sum",
        df %>%
          summarise(across(where(is.numeric), ~ sum(.))) #* this gives the total over all the reads
      )
    ) %>%
    mutate(
      PctReadPairsRetained =
        ifelse(test = studyID == "sum",
          yes = round((1 - ({{ initial_col }} - {{ final_col }}) / {{ initial_col }}), 1) * 100,
          no = PctReadPairsRetained
        )
    ) %>%
    mutate(across(where(is.numeric), ~ round(., rnd_to)))

  cat("Summarized workflow stages for '", deparse(substitute(df)), "' by group: '", deparse(substitute(grp_by)), "'\n", sep = "")

  return(summarised_table)
}


#' Format table for scientific notation
#'
#' @param table Input table to format
#'
#' @return A formatted table with scientific notation
format_table <- function(table) {
  formatted_table <- table %>%
    mutate(across(!contains("studyID") & !contains("PctReadPairsRetained"), ~ format(., scientific = TRUE, digits = 2)))

  cat("Formatted table with scientific notation.\n")
  return(formatted_table)
}


process_data <- function(data_list) {
  cat("Processing tables...\n")
  Table_SreadsAtEachStage_samples <- data_list$Table_SreadsAtEachStage_samples %>%
    rename(
      studyID = Study,
      SAMPLEID = Sample,
      InNonChimericASVs_final = InNonChimericASVs
    )

  workflowTable <- data_list$workflowTable %>% rename(
    studyID = Study,
    ReadsFilterAuids.Length_final = ReadsFilterAuids.Length
  )

  ### add percentage column
  Table_SreadsAtEachStage_samples <- Table_SreadsAtEachStage_samples %>%
    mutate(
      PctReadPairsRetained = round(((1 - (Initial - InNonChimericASVs_final) / Initial) * 100), 1)
      # test = round(((1 - (Initial - contains("final")) / Initial) * 100), 1)
    )

  workflowTable <- workflowTable %>%
    mutate(
      PctReadPairsRetained = ifelse(
        ReadsFilterAuids.Length_final == 0,
        0,
        round(((1 - (ReadsPipeline -
          ReadsFilterAuids.Length_final) / ReadsPipeline) * 100), 1)
      )
    )


  #* # summarise DADA2 nifH pipeline
  Table_SreadsAtEachStage_samples_summary <- summarise_workflow_stages_table(
    df = Table_SreadsAtEachStage_samples,
    grp_by = studyID,
    rnd_to = 1,
    initial_col = Initial,
    final_col = InNonChimericASVs_final
  )

  #* ## summarise workflow
  workflowTable_summary <- summarise_workflow_stages_table(
    df = workflowTable,
    grp_by = studyID,
    rnd_to = 1,
    initial_col = ReadsPipeline,
    final_col = ReadsFilterAuids.Length_final
  )

  ##* format table
  Table_SreadsAtEachStage_samples_summary_format <- Table_SreadsAtEachStage_samples_summary %>%
    format_table()

  workflowTable_summary_format <- workflowTable_summary %>%
    format_table()

  #-##  fix column names

  (Table_SreadsAtEachStage_samples_summary_format_rename <- Table_SreadsAtEachStage_samples_summary_format %>%
    rename(
      "Study ID" = studyID,
      # Initial = Initial,
      Trimmed4 = Trimmed,
      Filtered4 = Filtered,
      Merged7 = Merged,
      "non-Bimera9" = InNonChimericASVs_final,
      "Retained (%)" = PctReadPairsRetained
    ))


  (workflowTable_summary_format_rename <- workflowTable_summary_format %>%
    # rename("DADA2 nifH pipeline" = ReadsPipeline
    rename(
      "Study ID" = studyID,
      "DADA2 nifH pipeline" = ReadsPipeline,
      GatherAsvs = ReadsGatherAsvs,
      "Small samples" = ReadsFilterAuids.SmallSamp,
      Rare = ReadsFilterAuids.Rare,
      NonNifH = ReadsFilterAuids.NonNifH,
      Length = ReadsFilterAuids.Length_final,
      "Retained (%)" = PctReadPairsRetained
    ))
  # mutate(across(where(is.numeric), ~sum(.)))

  cat("Processed data and prepared summaries.\n")

  return(list(
    "Table 4" = Table_SreadsAtEachStage_samples_summary_format_rename,
    "Table 5" = workflowTable_summary_format_rename
  ))
}


main <- function(
    files_to_read,
    files_in_path,
    files_out_path) {
  # Load data from CSV files
  data_list <- load_files(files_to_read, files_in_path)

  # Process data to summarize and format tables
  final_results <- process_data(data_list)

  if (!is.null(final_results)) {
    # Create output directory if it doesn't exist
    create_dir(files_out_path)
    write_file_list(final_results, files_out_path)
  }
  # return(final_results)
}


if (sys.nframe() == 0 && !interactive()) {
  # Source needed files
  files_to_source <- c(
    "/Users/mo/Projects/nifH_amp_project/myWork/scripts/functions.R",
    "/Users/mo/Projects/nifH_amp_project/myWork/scripts/basic_plotting.R"
  )
  source_file(files_to_source)

  # Set up and parse arguments
  parser <- setup_parser()
  args <- parse_arg(parser)

  validate_parsed_args(args)

  main(
    files_to_read = args$files_to_read,
    files_in_path = args$files_in_path,
    files_out_path = args$files_out_path
  )
}
