#!/usr/bin/env Rscript

## load libraries
library(tidyverse)
library(argparser)

# load in functions and plotting parameter
source("/Users/mo/Projects/nifH_amp_project/myWork/scripts/functions.R")

source("/Users/mo/Projects/nifH_amp_project/myWork/scripts/basic_plotting.R")

# set working directory
# setwd("/Users/mo/Projects/nifH_amp_project/myWork")

### _ Loading in the data
#  Print out what script is running 
script_name <- basename(commandArgs(trailingOnly = FALSE)[4])
cat("Running script:", script_name, "\n")
cat("Load in the data\n")


#' Set up argument parser
#' 
#' @return A configured argument parser object
setup_parser <- function() {
  parser <- arg_parser("Load and process nifH amplicon data from R workspace")
  parser <- add_argument(parser = parser,
  arg = "--workspace",
  help = "Path to R workspace file (.Rdata)",
  default = "data/workspace/NifH_ASV_DB_March28_after_drop_19_empty_samples/workspace.RData")
  parser <- add_argument(parser,
    arg = "--output_path",
    help = "Output directory path",
    default = "analysis/out_files"
  )

  return(parser)
}


#' Load R workspace
#' 
#' @param workspace_path Path to R workspace file
#' @return Invisible NULL (objects are loaded into the global environment)
load_workspace <- function(workspace_path) {
  if ( !file.exists(workspace_path)) {
    stop("Workspace file does not exist: ", workspace_path)
  }

  load(workspace_path, envir = .GlobalEnv)

  cat("Workspace loaded into global environment.\n")
  
  invisible(NULL)
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

  for ( i in seq_along(old_names)) {
    old_name <- old_names[i]
    new_name <- new_names[i]

    if ( exists(old_name, envir = .GlobalEnv)) {
      transformed_objects[[new_name]] <- get(old_name)

      # Special case for cmapTab
      if (old_name == "cmapTab") {
        transformed_objects[[new_name]] <- rename(transformed_objects[[new_name]], studyID = StudyID)
    }

    cat("Converting file ", old_name ," ---> ", new_name, "...\n")
    } else {
      warning(cat("Input object",  old_name, "not found. Skipping.\n"))
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

  write_file_list(
  file_list = transformed_objects,
  path = files_out_path
)
  cat("\nDone loading script!!!\n")

  cat("Woooooooohooooooo!!!\n\n")
}


# Run the script if it's being executed directly
if (sys.nframe() == 0 && !interactive()) {
  parser <- setup_parser()
  args <- parse_args(parser)
  
  main(args$workspace, args$output_path)
}
