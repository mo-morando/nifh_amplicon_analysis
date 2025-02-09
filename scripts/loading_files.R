#!/usr/bin/env Rscript

#' @title NifH Amplicon Data Loading
#' @description This script processes NifH amplicon data by loading workspace files,
#'              transforming data objects, and saving processed results. It provides
#'              a robust and modular approach to handle NifH amplicon sequencing data
#'              analysis, with error handling and informative output.
#'
#' @details The pipeline performs the following steps:
#' * Sources required utility functions and plotting scripts
#' * Parses command-line arguments for workspace and output paths
#' * Loads workspace data with error handling
#' * Transforms workspace objects to a standardized format
#' * Writes processed data to output files
#'
#' The script includes several key functions:
#' * source_file(): Sources external R scripts with error handling
#' * setup_parser(): Sets up command-line argument parsing
#' * parse_arg(): Parses command-line arguments
#' * load_workspace(): Loads R workspace with error handling
#' * transform_workspace_objects(): Standardizes data object names and formats
#' * main(): Orchestrates the entire data processing workflow
#'
#' @usage Rscript script_name.R [--workspace PATH] [--output_path PATH]
#'
#' @param --workspace Path to R workspace file (.Rdata)
#' @param --output_path Output directory path
#'
#' @author Michael Morando
#' @date 2023
#'
#' @note This script requires the following R packages: tidyverse, argparser
#'
#' @examples
#' Rscript script_name.R --workspace ../data/workspace/Sep18/workspace.RData --output_path ../analysis/out_files
#'
#' @export


# Load required libraries
tryCatch(
  {
    suppressPackageStartupMessages({
      library(tidyverse)
      library(argparser)
    })
  },
  error = function(e) {
    cat("Error call in:", deparse(conditionCall(e)), "\n")
    stop(cat("Required packages could not be loaded due to:\n", conditionMessage(e), "\n"))
  }
)

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
  parser <- arg_parser("Load and process nifH amplicon data from R workspace")
  parser <- add_argument(
    parser = parser,
    arg = "--workspace",
    help = "Path to R workspace file (.Rdata)",
    default = "../data/workspace/Sep18/workspace.RData"
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
  tryCatch(
    {
      # Parse the arguments
      argv <- parse_args(parser)

      return(list(
        workspace = argv$workspace,
        files_out_path = argv$output_path
      ))
    },
    # warning = function(w) {
    #   cat("Warning occurred:\n")
    #   cat("Warning message:", conditionMessage(w), "\n")
    #   cat("Warning call:", deparse(conditionCall(w)), "\n")
    # },
    error = function(e) {
      cat("Error call in:", deparse(conditionCall(e)), "\n")
      stop(paste("Error parsing arguments due to:\n", conditionMessage(e)))
    }
  )
}


#' Load R workspace with error handling
#'
#' @param workspace_path Character string specifying path to R workspace file
#' @return Invisible NULL
#' @examples
#' load_workspace("path/to/workspace.RData")
#' @export
load_workspace <- function(workspace_path) {
  tryCatch(
    {
      if (!file.exists(workspace_path)) {
        stop("Workspace file does not exist: ", workspace_path)
      }

      load(workspace_path, envir = .GlobalEnv)

      cat("Workspace loaded successfully from:", workspace_path, "\n")

      return(invisible(NULL))
    },
    error = function(e) {
      cat("Error call in:", deparse(conditionCall(e)), "\n")
      stop(cat("Workspace loading failed due to:\n", conditionMessage(e), "\n"))
    }
  )
}


#' Transform workspace objects
#'
#' This helper function takes workspace objects and transforms them to conform
#'  analysis pipeline
#'
#' @return A list containing transformed objects
#' @export
#'
#' @importFrom dplyr rename
transform_workspace_objects <- function() {
  # Define old and new names
  old_names <- c("annotTab", "cmapTab", "abundTab", "relabundTab", "metaTab")
  new_names <- c("annoNifHDB_updt", "cmap_coloc", "nifhDB_cnts", "nifhDB_RA", "metaTab")

  transformed_objects <- list()

  for (i in seq_along(old_names)) {
    old_name <- old_names[i]
    new_name <- new_names[i]

    if (exists(old_name, envir = .GlobalEnv)) {
      transformed_objects[[new_name]] <- get(old_name)

      # Special case for cmapTab
      if (old_name == "cmapTab") {
        transformed_objects[[new_name]] <- rename(transformed_objects[[new_name]], studyID = StudyID)
      }

      cat("Converting file ", old_name, " ---> ", new_name, "...\n")
    } else {
      warning(cat("Input object", old_name, "not found. Skipping.\n"))
    }
  }

  return(transformed_objects)
}


#' Main function to process the data
#'
#' @param workspace_path Path to the R workspace file
#' @param files_out_path Output directory path
#' @return NULL
main <- function(workspace_path, files_out_path) {
  load_workspace(workspace_path)

  transformed_objects <- transform_workspace_objects()

  # Create output directory if it doesn't exist
  create_dir(files_out_path)

  if (!is.null(transformed_objects)) {
    write_file_list(
      file_list = transformed_objects,
      path = files_out_path
    )
  }
  cat("\nDone loading R workspace.\n\n")
}


# Run the script if it's being executed directly
if (sys.nframe() == 0 && !interactive()) {
  source_file(files_to_source)

  parser <- setup_parser()
  args <- parse_arg(parser)

  # validate_parsed_args(parsed_args = args)

  main(
    workspace_path = args$workspace,
    files_out_path = args$files_out_path
  )
}
