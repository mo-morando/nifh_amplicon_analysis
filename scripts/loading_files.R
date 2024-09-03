#!/usr/bin/env Rscript

## load libraries
library(tidyverse)

source("/Users/mo/Projects/nifH_amp_project/myWork/scripts/functions.R")

source("/Users/mo/Projects/nifH_amp_project/myWork/scripts/basic_plotting.R")

### set working directory
setwd("/Users/mo/Projects/nifH_amp_project/myWork")

### _ Loading in the data
script_name <- basename(commandArgs(trailingOnly = FALSE)[4])
cat("Running script:", script_name, "\n")
cat("Load in the data\n")

#-# Load in R workspace from nifh workflow
load("data/workspace/NifH_ASV_DB_March28_after_drop_19_empty_samples/workspace.RData")



# # Loop through each object and print its dimensions
# for (obj in list) {
#   if (!is.function(get(obj))) { # Exclude functions
#     cat("Object:", obj, "\n")
#     cat("Dimensions:", paste(dim(get(obj)), collapse = "x"), "\n\n")
#   }
# }

#' Print Dimensions of Non-Function Objects in a List
#' 
#' This function iterates through each object in a list, excluding functions,
#' and prints the object's name and dimensions.
#' 
#' @param list A character vector or list containing names of objects in the workspace.
#' @return NULL
#' @export
#' 
#' @examples
#' print_object_dimensions(list_of_objects)
#' 
#' @seealso
#' /code{\link{print}}
#' 
print_object_dimensions <- function(list) {
  # Loop through each object in the list
  for (obj in list) {
    # Check if the object is not a function
    if (!is.function(get(obj))) {
      # Print the object's name
      cat("Object:", obj, "\n")
      # Print the dimensions of the object
      if (!is.null(dim(get(obj)))) {
        cat("Dimensions:", paste(dim(get(obj)), collapse = "x"), "\n\n")
      }
      cat("Dimensions: None\n\n")
    }
  }
}

# Get the names of all objects in the workspace
workspace_objects <- ls()

print_object_dimensions(workspace_objects)


#-# Object: annotTab
# Name I use in my scripts for this file
cat("Converting file annotTab ---> annoNifHDB_updt...\n")
annoNifHDB_updt <- annotTab


### - Object: cmapTab
# Name I use in my scripts for this file
cat("Converting file cmapTab ---> cmap_coloc...\n")
cmap_coloc <- cmapTab %>%
  ### columns to rename to fit with my other scripts
  rename(studyID = StudyID)


### - Object: metaTab
# this file now exists but we don't change the name



### - count/abundance tables
cat("Converting file nifhDB_cnts ---> nifhDB_cnts...\n")
nifhDB_cnts <- abundTab
cat("Converting file relabundTab ---> nifhDB_RA...\n")
nifhDB_RA <- relabundTab

# ## check to verify that this worked
# ## Values should all be 1
# rel_abun_check <- nifhDB_RA %>%
#   summarise(across(where(is.numeric), sum))

# ## Verify all columns sum to ~1, and all cols are type double.
# stopifnot(rel_abun_check - 1 < 1e-9)
# stopifnot(sapply(nifhDB_RA[, -1], class) == "numeric")



# Get the names of all objects in the workspace
workspace_objects_final <- ls()

for (obj in workspace_objects_final) {
   cat( obj, class(get(obj)), "\n", sep = '\t')
}

## Sort a list of all files in the environment for the files to write out
## These are tibbles/data frames and lists
object_list <- list()
non_object_list <- list()
for (obj in workspace_objects_final) {
  if (!is.function(get(obj))) { # Exclude functions
#   if (is.data.frame(get(obj)) || is_tibble(get(obj)) || is_list(get(obj))) { # Exclude functions

    cat(
      "Adding file:", obj, "to list",
      deparse(substitute(object_list)), "\n"
    )

    object_list[[obj]] <- get(obj)
  } else {
    non_object_list <- c(non_object_list, obj)
  }
}
# Print out the entire non_data_frame_list with each value on a separate line
cat("Non-data frame or list objects that were not added to list:", paste("\t", non_object_list), sep = "\n")

# cat("Non-data frame or list objects that were not added to list:\n")
# for (obj in non_object_list) {
#   cat("\t", obj, "\n")
# }



#-# Define directory to write files out to
files_out_path <- "analysis/out_files"

create_dir(files_out_path)

#' Write files list
#'
#' This function writes objects from a list to files with specified names, file paths and extentions.
#'
#' @param file_list A list of tibbles to be written to CSV files
#' @param out_ext Optional. The file extension to be added to the file names
#' @param path Optional. The directory where the files will be written. Defaults to the current directory
#'
#' @return NULL
#'
#' @export
#'
#' @examples
#' write_file_list(
#'  file_list = object_list,
#'  path = files_out_path,
#'  out_ext = ".csv"
#')
#' write_files(list(dataframe1, dataframe2), names(list(dataframe1 = dataframe1, dataframe2 = dataframe2)), ".csv")
#'
#' @seealso
#' /code{\link{main}}
#'
#' @author Michael (Mo) Morando
#'
#'
write_file_list <- function(file_list, path = ".", out_ext = ".csv") {
   cat("Writing files from list to path:", path, "\n")
   
   # Ensure file_list is a list and it has names
   if (!is.list(file_list) && !is.null(names(file_list))) {
	  stop("file_list must be a list, and have associated names.")
   }

   # Extract names to list
   file_list_names <- names(file_list)
   
   # Iterate over the file_list and save files
   for (i in seq_along(file_list)) {
	file_path <- file.path(path, paste0(file_list_names[i], out_ext))

	# Determine the type and write accordingly
	file_data <- file_list[[i]]
	if (inherits(file_data, c("tbl_df", "data.frame"))) {
		write_csv(file_data, file_path)
	}

	else if (is.list(file_data)) {
	   write_csv(tibble( "temp_col_id" = unlist(file_data)), file = file_path)
	}

	else if (is.character(file_data)) {
	   if(!is.null(names(file_data)) && all(names(file_data) != "")) {
		write_csv(tibble("key" = names(file_data), "value" = file_data), 
		file = file_path)
	   } else {
		  writeLines(file_data, con = file_path, sep = ",")
	   }
	
	}else {
	   cat("Skipping file:", file_list_names[i], "- unsupported type.\n")
	   next
	}

    cat("Wrote file:", file_list_names[i], "\n")
	}
}


write_file_list(
  file_list = object_list,
  path = files_out_path
)




### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data

cat("Done loading script!!!
")

cat("Woooooooohooooooo!!!")
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
