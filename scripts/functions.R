#' Module Docstring:
#'
#' This script contains a collection of functions designed to perform various data processing tasks on dataframes using the tidyverse package in R. Below are the main functionalities of the script:
#'
#' - Merging Functions:
#'     - `merge_cmap`: Merge dataframe with Colocalization Map data.
#'     - `merge_annotations`: Merge dataframe with annotation table.
#'
#' - Data Transformation:
#'     - `transform_data_lng`: Transform data from wide to long format.
#'
#' - Data Manipulation:
#'     - `add_total_row`: Add a total row to the dataframe.
#'     - `add_group_id`: Add group ID to the dataframe.
#'     - `dedup_by_group`: Remove duplicates based on group ID.
#'     - `main_average_ra_dedup_by_group`: Calculate average and remove duplicates based on group ID.
#'     - `remove_samples_nucleic_acid`: Remove samples based on nucleic acid type.
#'     - `count_and_arrange`: Count occurrences and arrange dataframe.
#'     - `add_percentage`: Add a percentage column to the dataframe.
#'     - `sum_tax`: Summarize taxonomic data.
#'     - `clean_percentages`: Clean percentage data.
#'     - `remove_aphotic_samples`: Remove aphotic samples from the dataframe.
#'     - `add_size_frac_key`: Add size fraction column to dataframe.
#'     - `add_rep_flag`: Add replicate flag column to identify replicate samples.
#'
#' These functions can be used individually or in combination to preprocess and analyze data efficiently.
#'
#' Load required libraries
cat("Loading the tidyverse library...\n")
library(tidyverse)


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



#' Define function to create a new directory
cat("Defining function to create a new directory...\n")

#' Create directory if it doesn't exist
#'
#' This function checks if a directory exists at the specified path.
#' If the directory does not exist, it creates the directory.
#'
#' @param dir_path The path to the directory to be checked/created.
#' @return None
#' @export
#' @examples
#' \dontrun{
#' create_dir("path/to/your/directory")
#' }
create_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path)
    cat(paste("Directory: '", dir_path, "' created\n", sep = ""))
  } else {
    cat(paste("Directory: '", dir_path, "' already exists\n", sep = ""))
  }
}


#' Validate parsed arguments
#'
#' Validate the objects returned by the parse_arg() function
#'
#' @param parsed_args A list containing the parsed arguments returned by parse_arg()
#'
#' @return NULL
#'
#' @example
#' parser <- setup_parser()
#' args <- parse_arg(parser)
#' validate_parsed_args(args)
validate_parsed_args <- function(parsed_args) {
  # Check if files_to_read is a character vector
  if (!is.character(parsed_args$files_to_read)) {
    stop("files_to_read must be a character vector")
  }

  # Check if files_in_path is a valid directory path
  if (!dir.exists(parsed_args$files_in_path)) {
    stop("files_in_path: '", parsed_args$files_in_path, "' must be a valid directory path")
  }

  # Check if files_out_path is a valid directory path
  if (!dir.exists(parsed_args$files_out_path)) {
    warning("files_in_path: '", parsed_args$files_in_path, "' must be a valid directory path")
  }

  # Check if both plot_ext and plot_device are provided
  if (!is.null(parsed_args$plot_ext) && !is.null(parsed_args$plot_device)) {
    plot_ext_check <- gsub("\\.", "", parsed_args$plot_ext)
    if (plot_ext_check %in% "pdf") {
      parsed_args$plot_device <- "pdf"
    }
    if (plot_ext_check != parsed_args$plot_device) {
      stop("plot_ext: '", parsed_args$plot_ext, "' must match plot_device: 
      '", parsed_args$plot_device, "'")
    }
  }

  # Check if each file in files_to_read exists in files_in_path
  missing_files <- character(0)
  for (file in parsed_args$files_to_read) {
    # if (!file.exists(file.path(parsed_args$files_in_path, file))) {
    if (!any(startsWith(list.files(parsed_args$files_in_path), file))) {
      missing_files <- c(missing_files, file)
    }
  }
  if (length(missing_files > 0)) {
    stop("The following files are missing in files_in_path:\n", paste("\t",missing_files, collapse = "\n"))
  }

  # All validations passed
  cat("All parsed arguments are valid.\n")
}





#' Load and Assign CSV Files
#'
#' @param file_list Character vector of file names to load
#' @param path Directory path for input files
#' @param verbose Logical. Whether to print status messages
#' @return List of loaded data frames or character vectors
#' @importFrom readr read_csv
#' @examples
#' \dontrun{
#' files <- c("data1", "data2", "single_line")
#' data_list <- load_files(files, "path/to/csv/files")
#' }
load_files <- function(file_list, path, verbose = TRUE) {
  data_list <- list()

  for (file in file_list) {
    file_path <- file.path(path, paste0(file, ".csv"))
    tryCatch(
      {
        if (file.exists(file_path)) {
          if (verbose) cat("Loading file:", file_path, "\n")

          # Check if the file has only one line
          lines <- suppressWarnings(readLines(file_path, n = 2)) # Read only the first two lines
          if (length(lines) == 1) {
            # File contains a single line, read as character vector
            data_list[[file]] <- read_csv(file_path, show_col_types = FALSE, col_names = FALSE)
            if (verbose) cat("  File loaded as a data frame but without column
          headers.\n")
          } else {
            # File contains multiple lines, read as a data frame
            data_list[[file]] <- read_csv(file_path, show_col_types = FALSE)
            if (verbose) {
              cat(
                "  File loaded as a data frame with",
                nrow(data_list[[file]]), "rows and",
                ncol(data_list[[file]]), "columns.\n"
              )
            }
          }
        } else {
          warning(paste("File not found:", file_path))
        }
      },
      error = function(e) {
        cat("Error loading file", file_path, ":", conditionMessage(e), "\n")
      }
    )
  }
  if (verbose) cat("Finished loading", length(data_list), "file.\n\n")
  return(data_list)
}


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



cat("Defining function to write files from a list...\n")
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
      # file_path <- file.path(path, paste0(file_list_names[i], ".tsv"))
		  writeLines(file_data, con = file_path, sep = ",")
	   }
	
	}else {
	   cat("Skipping file:", file_list_names[i], "- unsupported type.\n")
	   next
	}

    cat("Wrote file:", file_list_names[i], "\n")
	}
}



#' Define merging functions
cat("Defining merging functions...\n")

#' Merge dataframe with Colocalization Map data from CMAP.
#'
#' This function merges the input dataframe with the Colocalization Map data. It performs a left join based on the specified column.
#'
#' @param df Input dataframe.
#' @param cmap Colocalization Map dataframe, default is cmap_coloc.
#' @param by_join Column name to join on, default value is "SAMPLEID".
#'
#' @return Merged dataframe.
#'
merge_cmap <- function(df, cmap = cmap_coloc, by_join = "SAMPLEID") {
  cat("Merging dataframe with Colocalization Map data...\n")
  if ("group_id" %in% names(cmap)) {
    cmap <- cmap %>%
      select(-group_id)
  }
  merged_df <- df %>%
    left_join(cmap, by = by_join)
  return(merged_df)
}

#' Merge dataframe with annotation table.
#'
#' This function merges the input dataframe with the annotation table. It performs a left join based on the specified column.
#'
#' @param df Input dataframe.
#' @param annotation_table Annotation table dataframe, default is annoNifHDB_updt.
#' @param by_join Column name to join on, default is "AUID".
#'
#' @return Merged dataframe.
#'
merge_annotations <- function(df, annotation_table = annoNifHDB_updt, by_join = "AUID") {
  cat("Merging dataframe with annotation table...\n")
  merged_df <- df %>%
    left_join(annotation_table, by = by_join)
  return(merged_df)
}

#' Define data transformation function
cat("Defining data transformation function...\n")

#' Transform data from wide to long format.
#'
#' This function transforms the input dataframe from wide to long format using pivot_longer function from tidyverse package.
#'
#' @param input_df Input dataframe.
#' @param starts_with_col Prefix for columns to pivot.
#' @param names_to_col Column name for pivoted column names.
#' @param values_to_col Column name for pivoted values.
#'
#' @return Transformed dataframe.
#'
transform_data_lng <- function(input_df, starts_with_col, names_to_col, values_to_col) {
  cat("Transforming data from wide to long format...\n")
  result <- input_df %>%
    {
      if ("group_id" %in% names(input_df)) select(., -group_id) else .
    } %>%
    pivot_longer(
      cols = starts_with(starts_with_col),
      names_to = if ("group_id" %in% names(input_df)) starts_with_col else names_to_col,
      values_to = values_to_col
    )
  return(result)
}

# transform_data_lng <- function(input_df, starts_with_col, names_to_col, values_to_col) {
#   cat("Transforming data from wide to long format...\n")
#   if ("group_id" %in% names(input_df)) {
#     result <- input_df %>%
#       select(-group_id) %>%
#       pivot_longer(cols = starts_with(starts_with_col),
#                    names_to = starts_with_col,
#                    values_to = values_to_col)
#   } else {
#     result <- input_df %>%
#       pivot_longer(cols = starts_with(starts_with_col),
#                    names_to = names_to_col,
#                    values_to = values_to_col)
#   }
#   return(result)
# }


#' Define function to remove aphotic samples from a tibble
cat("Defining function to remove aphotic samples from a tibble...\n")

#' Remove photic samples from the dataframe, returning just the photic samples
#'
#' This function removes aphotic samples from the dataframe and returns only the photic samples.
#'
#' @param df Input dataframe.
#' @param photic_samples Photic sample key used to filter dataframe, default is photic_samples_key.
#'
#' @return Dataframe containing only photic samples.
#'
remove_aphotic_samples <- function(df, photic_samples = photic_samples_key) {
  subbed_df <- df %>%
    filter(SAMPLEID %in% photic_samples)

  return(subbed_df)
}


#' Define function to add a total row
cat("Defining function to add a total row...\n")

#' Add total row to the dataframe.
#'
#' This function adds a total row to the dataframe, summarizing the specified column or all numeric columns if not specified.
#'
#' @param df Input dataframe.
#' @param summary_column Column to summarize.
#' @param column_name Name of the new total row.
#' @param pull_name Name of the column to pull for summary.
#' @param all_columns All columns to include in the total row.
#'
#' @return Dataframe with total row added.
#'
add_total_row <- function(df, summary_column = NULL, column_name, pull_name = NULL, all_columns) {
  cat("Adding total row to the dataframe...\n")
  if (all_columns == FALSE) {
    df_with_total_count <- df %>%
      bind_rows(tibble(
        !!column_name := "Total",
        {{ summary_column }} := sum(df %>% pull({{ summary_column }}))
      ))
  } else {
    df_with_total_count <- df %>%
      bind_rows(tibble(
        !!column_name := "Total",
        df %>% summarise(across(where(is.numeric), ~ sum(.)))
      ))
  }
  return(df_with_total_count)
}

#' Define function to add a size fraction row
cat("Defining function to add a size fraction row...\n")
#' Add size fraction row to dataframe
#'
#' This function adds a size fraction column to the dataframe based on the specified join key.
#'
#' @param df Input dataframe.
#' @param join_key Colocalization Map dataframe.
#' @param remove_col Column name to join on.
#'
#' @return Dataframe with size fraction column added
#'
add_size_frac_key <- function(
    df,
    join_key = size_fraction_key,
    remove_col = num_distinct_size_fractions) {
  cat("Add '", deparse(substitute(join_key)), "' column to dataframe: ", deparse(substitute(df)), " ...\n", sep = "")
  df_size_frac <- df %>%
    left_join(join_key) %>%
    suppressMessages() %>%
    select(-{{ remove_col }})

  return(df_size_frac)
}


# Now I want to identify replicate samples
# To do this, I think I have to separate nucleic acids into their own dfs
# The make a replicate flag to identify reps


#' Merge dataframe with Colocalization Map data from CMAP.
#'
#' This function merges the input dataframe with the Colocalization Map data. It performs a left join based on the specified column.
#'
#' @param df Input dataframe.
#' @param cmap Colocalization Map dataframe, default is cmap_coloc.
#' @param by_join Column name to join on, default value is "SAMPLEID".
#'
#' @return Merged dataframe.
#'
merge_cmap <- function(df, cmap = cmap_coloc, by_join = "SAMPLEID") {
  cat("Merging dataframe with Colocalization Map data...\n")
  if ("group_id" %in% names(cmap)) {
    cmap <- cmap %>%
      select(-group_id)
  }
  merged_df <- df %>%
    left_join(cmap, by = by_join)
  return(merged_df)
}

#' Merge dataframe with annotation table.
#'
#' This function merges the input dataframe with the annotation table. It performs a left join based on the specified column.
#'
#' @param df Input dataframe.
#' @param annotation_table Annotation table dataframe, default is annoNifHDB_updt.
#' @param by_join Column name to join on, default is "AUID".
#'
#' @return Merged dataframe.
#'
merge_annotations <- function(df, annotation_table = annoNifHDB_updt, by_join = "AUID") {
  cat("Merging dataframe with annotation table...\n")
  merged_df <- df %>%
    left_join(annotation_table, by = by_join)
  return(merged_df)
}

#' Transform data from wide to long format.
#'
#' This function transforms the input dataframe from wide to long format using pivot_longer function from tidyverse package.
#'
#' @param input_df Input dataframe.
#' @param starts_with_col Prefix for columns to pivot.
#' @param names_to_col Column name for pivoted column names.
#' @param values_to_col Column name for pivoted values.
#'
#' @return Transformed dataframe.
#'
transform_data_lng <- function(input_df, starts_with_col, names_to_col, values_to_col) {
  cat("Transforming data from wide to long format...\n")
  result <- input_df %>%
    {
      if ("group_id" %in% names(input_df)) select(., -group_id) else .
    } %>%
    pivot_longer(
      cols = starts_with(starts_with_col),
      names_to = if ("group_id" %in% names(input_df)) starts_with_col else names_to_col,
      values_to = values_to_col
    )
  return(result)
}

#' Add replicate flag column to dataframe to identify replicate samples
#'
#' This function adds a replicate flag column to the dataframe to identify which samples are replicates.
#'
#' @param df Input dataframe.
#'
#' @return Dataframe with replicate flag column added
#'
add_rep_flag <- function(df) {
  cat("Add replicate flag column to dataframe: ", deparse(substitute(df)), " ...\n", sep = "")
  df_w_rep <- df %>%
    group_by(sample_point) %>%
    mutate(
      replicate_flag = case_when(
        (num_dist_samp_pnts == 1) ~ "Single_Sample",
        # (num_dist_samp_pnts > 2) ~ "Replicate_Sample",
        (n_distinct(nucleicAcidType) == 1 & size_frac_flag == "Two_Size_Fractions") ~ "Single_Sample",
        (n_distinct(nucleicAcidType) == 1 & size_frac_flag == "Two_Size_Fractions" & num_dist_samp_pnts > 3) ~ "Replicate_Sample",
        (n_distinct(nucleicAcidType) == 1 & size_frac_flag == "One_Size_Fractions") ~ "Replicate_Sample",
      )
    ) %>%
    ungroup()

  return(df_w_rep)
}




#' Define function to add group ID
cat("Defining function to add group ID...\n")

#' Add group ID to the dataframe.
#'
#' This function adds a group ID column to the dataframe based on the specified column.
#'
#' @param df Input dataframe.
#' @param var_select Column to select.
#' @param by_var Column name to join on.
#'
#' @return Dataframe with group ID added.
#'
add_group_id <- function(df, group_id_key = unique_sample_id_key, var_select = SAMPLEID, by_var = "SAMPLEID") {
  cat("Adding group ID to the dataframe...\n")
  df_with_group_id <- df %>%
    left_join(group_id_key %>%
      select({{ var_select }}, group_id), by = by_var)
  return(df_with_group_id)
}

#' Define function to remove duplicates based on group ID
cat("Defining function to remove duplicates based on group ID...\n")

#' Remove duplicates based on group ID.
#'
#' This function removes duplicates based on the group ID column in the dataframe.
#'
#' @param df Input dataframe.
#' @param ... Additional arguments to pass to distinct().
#'
#' @return Dataframe with duplicates removed.
#'
dedup_by_group <- function(df, group_id_key = unique_sample_id_key, ...) {
  cat("Removing duplicates based on group ID...\n")
  contains_group_id <- "group_id" %in% names(df)
  if (contains_group_id) {
    df_deup <- df %>%
      distinct(..., .keep_all = TRUE)
  } else {
    df_deup <- df %>%
      add_group_id(group_id_key) %>%
      distinct(..., .keep_all = TRUE)
  }
  rows_removed <- nrow(df) - nrow(df_deup)
  cat("Number of rows removed:", rows_removed, "\n")
  return(df_deup)
}

#' Define function to calculate average and remove duplicates based on group ID in the relative abundance table
cat("Defining function to calculate average and remove duplicates based on group ID in the relative abundance table...\n")

#' Calculate average and remove duplicates based on group ID in the relative abundance table.
#'
#' This function calculates the average and removes duplicates based on the group ID in the relative abundance table column in the dataframe.
#'
#' @param df_lng Input dataframe.
#' @param ... Additional arguments to pass to group_by() and dedup_by_group().
#' @param mean_by Column to compute mean.
#'
#' @return Dataframe with averages calculated and duplicates removed.
#'
main_average_ra_dedup_by_group <- function(df_lng, ..., grp_key, mean_by) {
  cat("Calculating average and removing duplicates based on group ID in the relative abundance table...\n")
  cat("This is going to take awhile...\n")
  df_lng_mean_ra_depup_by_group <- df_lng %>%
    add_group_id(group_id_key = {{ grp_key }}) %>%
    group_by(...) %>%
    mutate(average_value_AUID = mean({{ mean_by }}, na.rm = TRUE)) %>%
    ungroup() %>%
    dedup_by_group(group_id_key = unique_sample_id_key, ...) %>%
    ungroup()

    
    return(df_lng_mean_ra_depup_by_group)
}


#' Define function to remove duplicates based on group ID in count table
cat("Defining function to remove duplicates based on group ID in count table...\n")

#' remove duplicates based on group ID in count table.
#'
#' This function calculates the average and removes duplicates based on the group ID in count table column.
#'
#' @param df_lng Input dataframe (count table).
#' @param ... Additional arguments to pass to dedup_by_group().
#'
#' @return Dataframe (count table) with duplicates removed.
#'
main_cnts_dedup_by_group <- function(df_lng, ..., grp_key) {
  cat("Removing duplicates based on group ID in count table...\n")
  df_lng_mean_ra_depup_by_group <- df_lng %>%
    add_group_id(group_id_key = {{ grp_key }}) %>%
    dedup_by_group(group_id_key = unique_sample_id_key, ...) %>%
    ungroup() %>%
    return(df_lng_mean_ra_depup_by_group)
}


#' Define function to remove samples based on nucleic acid type
cat("Defining function to remove samples based on nucleic acid type...\n")

#' Remove samples based on nucleic acid type.
#'
#' This function removes samples from the dataframe based on the specified nucleic acid type.
#'
#' @param df Input dataframe.
#' @param nucleic_acid_type Type of nucleic acid to remove.
#' @param nucleic_acid_samples_key Samples key.
#'
#' @return Dataframe with samples removed based on nucleic acid type.
#'
remove_samples_nucleic_acid <- function(df, nucleic_acid_type, nucleic_acid_samples_key = DNA_samples_key) {
  nucleic_acid_flag <- c("DNA", "RNA") %in% nucleic_acid_type
  if (nucleic_acid_flag[1]) {
    na_removed <- "RNA"
    nucleic_acid <- df %>%
      filter(SAMPLEID %in% {{ nucleic_acid_samples_key }})
  } else if (nucleic_acid_flag[2]) {
    na_removed <- "DNA"
    nucleic_acid <- df %>%
      filter(!SAMPLEID %in% {{ nucleic_acid_samples_key }})
  } else {
    cat("You supplied '", nucleic_acid_type, ",' which is neither DNA nor RNA. Please make sure it is in ALL CAPS.\n\n")
    nucleic_acid <- NULL
  }
  if (!is.null(nucleic_acid)) {
    removed_samples <- nrow(df) - nrow(nucleic_acid)
    cat("Samples of type nucleic acid type '", nucleic_acid_type, "' were retained.\nNumber of ", na_removed, " samples removed: ", removed_samples, "\n\n", sep = "")
  }
  return(nucleic_acid)
}

#' Define function to count occurrences and arrange dataframe
cat("Defining function to count occurrences and arrange dataframe...\n")

#' Count occurrences and arrange dataframe.
#'
#' This function counts occurrences and arranges the dataframe based on the specified variable.
#'
#' @param data Input dataframe.
#' @param group_vars Variables to group by.
#' @param arrange_var Variable to arrange by.
#' @param count_col_name Name of the count column. Default value is "n"
#'
#' @return Dataframe with counts and arranged by specified variable.
#'
count_and_arrange <- function(
    data,
    group_vars,
    arrange_var = n,
    count_col_name = "n") {
  cat_text <- paste(group_vars, collapse = ", ")
  # cat("Counting occurrences and arranging dataframe by: '", deparse(substitute(group_vars)), "'...\n", sep = "")
  cat("Counting occurrences and arranging dataframe by: '", cat_text, "'...\n", sep = "")
  df_count_and_arrange <- data %>%
    count(across(all_of(group_vars)), name = count_col_name) %>%
    arrange(desc({{ arrange_var }}))
  return(df_count_and_arrange)
}

#' Define function to add percentage column to the dataframe
cat("Defining function to add percentage column to the dataframe...\n")

#' Add percentage column to the dataframe.
#'
#' This function adds a percentage column to the dataframe, which is calculated based on the specified column and grouping variable.
#'
#' @param df Input dataframe.
#' @param sum_by_percent Column to sum by.
#' @param percentage_id Name of the percentage column.
#' @param grouping_by Variable to group by. Default is NULL.
#' @param remove_columns Columns to remove. Default is NULL.
#'
#' @return Dataframe with percentage column added.
#'
add_percentage <- function(df, sum_by_percent, percentage_id, rnd_pct = 1, grouping_by = NULL, remove_columns = NULL) {
  cat("Adding percentage column to the dataframe that is summed by '", deparse(substitute(sum_by_percent)), "' and grouped by '", grouping_by, "'...\n", sep = "")
  df_with_total <- df %>%
    group_by({{ grouping_by }}) %>%
    mutate(
      total = sum({{ sum_by_percent }}),
      {{ percentage_id }} := round({{ sum_by_percent }} / total * 100, rnd_pct)
    ) %>%
    ungroup()
  if (!is.null({{ remove_columns }})) {
    df_with_total <- df_with_total %>%
      select(-{{ remove_columns }})
  }
  return(df_with_total)
}

#' Define function to summarize taxonomic data
cat("Defining function to summarize taxonomic data...\n")

#' Summarize taxonomic data.
#'
#' This function summarizes taxonomic data by joining abundance and annotation tables, calculating total counts and percentage, and arranging the dataframe.
#'
#' @param abundance_table Abundance table dataframe.
#' @param annotation_table Annotation table dataframe.
#' @param joining_by_col Column name to join on.
#' @param total_counts_id Name of the total counts column.
#' @param grouping_for_percentage Variable for percentage grouping. Default is NULL.
#' @param percentage_id Name of the percentage column.
#' @param grouping_by_1 First variable to group by. Default is NULL.
#' @param grouping_by_2 Second variable to group by. Default is NULL.
#' @param sum_by Column to sum by.
#' @param ... Additional arguments.
#'
#' @return Summarized dataframe.
#'
sum_tax <- function(abundance_table, annotation_table, joining_by_col, total_counts_id, grouping_for_percentage = NULL,
                    percentage_id, grouping_by_1 = NULL, grouping_by_2 = NULL, sum_by, ...) {
  cat("Summarizing taxonomic data...\n")
  result <- abundance_table %>%
    left_join(annotation_table, joining_by_col) %>%
    group_by({{ grouping_by_1 }}, {{ grouping_by_2 }}) %>%
    mutate({{ total_counts_id }} := sum({{ sum_by }}, na.rm = TRUE)) %>%
    ungroup() %>%
    distinct({{ total_counts_id }}, {{ grouping_by_1 }}, {{ grouping_by_2 }}, .keep_all = TRUE) %>%
    arrange(desc({{ total_counts_id }})) %>%
    group_by({{ grouping_for_percentage }}) %>%
    mutate(
      total = sum({{ total_counts_id }}),
      {{ percentage_id }} := {{ total_counts_id }} / total * 100
    ) %>%
    ungroup() %>%
    select(..., {{ total_counts_id }}, total, {{ percentage_id }}, everything())
  return(result)
}

#' Define function to clean percentage data
cat("Defining function to clean percentage data...\n")

#' Clean percentage data.
#'
#' This function cleans percentage data in the dataframe based on the specified threshold.
#'
#' @param df Input dataframe.
#' @param grouping_by Variable to group by.
#' @param clean_var Variable to clean.
#' @param threshold Threshold for cleaning. Default value is 1%.
#' @param percentage_id Name of the percentage column.
#' @param column_var_mutate Column to mutate.
#' @param column_var Original column.
#'
#' @return Cleaned dataframe.
#'
clean_percentages <- function(df, grouping_by, clean_var, threshold = 1.0, percentage_id,
                              column_var_mutate, column_var) {
  cat("Cleaning percentage data...\n")
  cat("Using clean_var:", deparse(substitute(clean_var)), "\n")
  cat("Threshold required: ", threshold, "%\n", sep = "")

  cleaned_df <- df %>%
    mutate({{ column_var_mutate }} := ifelse({{ clean_var }} < threshold, "other", {{ column_var }})) %>%
    group_by({{ column_var_mutate }}, {{ grouping_by }}) %>%
    summarise({{ percentage_id }} := sum({{ clean_var }})) %>%
    mutate({{ percentage_id }} := round({{ percentage_id }}, 1)) %>%
    ungroup() %>%
    arrange(desc({{ percentage_id }}))

  if (nrow(cleaned_df) < nrow(df)) {
    removed_rows <- df %>%
      filter({{ clean_var }} < threshold) %>%
      #   distinct({{ column_var }})  %>%
      pull({{ column_var }})
    unique_removed_rows <- unique(removed_rows)
    cat("Unique column headers that did not reach the threshold and were converted to 'other':\n")
    cat(paste(unique_removed_rows, collapse = ", "), "\n")
    cat("Number of rows removed:", length(removed_rows), "\n\n")
  }

  return(cleaned_df)
}

#' Provide usage instructions
cat("\nSCRIPT LOADED SUCCESSFULLY!\n\n")
cat("\nUSAGE:\n\n")
cat("merge_cmap(df, cmap = cmap_coloc, by_join = 'SAMPLEID'): Merge dataframe with Colocalization Map data.\n\n")
cat("merge_annotations(df, annotation_table = annoNifHDB_updt, by_join = 'AUID'): Merge dataframe with annotation table.\n\n")
cat("transform_data_lng(input_df, starts_with_col, names_to_col, values_to_col): Transform data from wide to long format.\n\n")
cat("add_total_row(df, summary_column = NULL, column_name, pull_name = NULL, all_columns): Add total row to the dataframe.\n\n")
cat("add_group_id(df, var_select = SAMPLEID, by_var = 'SAMPLEID'): Add group ID to the dataframe.\n\n")
cat("dedup_by_group(df, ...): Remove duplicates based on group ID.\n\n")
cat("main_average_ra_dedup_by_group(df_lng, ..., mean_by): Calculate average and remove duplicates based on group ID in the relative abundance table.\n\n")
cat("main_cnt_dedup_by_group(df_lng, ...): Calculate average and remove duplicates based on group ID in the count table.\n\n")
cat("remove_samples_nucleic_acid(df, nucleic_acid_type, nucleic_acid_samples_key = DNA_samples_key): Remove samples based on nucleic acid type.\n\n")
cat("count_and_arrange(data, group_vars, arrange_var = n, count_col_name = 'n'): Count occurrences and arrange dataframe.\n\n")
cat("add_percentage(df, sum_by_percent, percentage_id, grouping_by = NULL, remove_columns = NULL): Add percentage column to the dataframe.\n\n")
cat("sum_tax(abundance_table, annotation_table, joining_by_col, total_counts_id, grouping_for_percentage = NULL, percentage_id, grouping_by_1 = NULL, grouping_by_2 = NULL, sum_by, ...): Summarize taxonomic data.\n\n")
cat("clean_percentages(df, grouping_by, clean_var, threshold = 1.0, percentage_id, column_var_mutate, column_var): Clean percentage data.\n\n")
cat("remove_aphotic_samples(df, photic_samples = photic_samples_key): Remove aphotic samples from the dataframe.\n\n")
cat("add_size_frac_key(df, join_key, remove_col): Add size fraction column to dataframe.\n\n")
cat("add_rep_flag(df): Add replicate flag column to identify replicate samples.\n\n")



### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
