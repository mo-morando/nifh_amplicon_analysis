#!/usr/bin/env Rscript

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
    }#,
    # warning = function(w) {
    #   cat("Warning while sourcing : ", file_path, ":", conditionMessage(w), "\n")
    # }
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

  # Convert the comma-separated string to a vector
  files_to_read <- strsplit(argv$files, ",")[[1]]
  files_in_path <- argv$input_path
  files_out_path <- argv$output_path

  return(list(
    files_to_read = files_to_read,
    files_in_path = files_in_path,
    files_out_path = files_out_path
  ))
}

#' Format and Print Table
#'
#' This function formats a table with custom options and prints it.
#'
#' @param data A data frame or tibble to be formatted and printed.
#' @param zero_print Character to use for printing zero values. Default is ".".
#' @param max_rows Maximum number of rows to print. Default is 1000.
#' @param ... Additional arguments to be passed to format() function.
#'
#' @return Invisibly returns the formatted table if successful, or NULL if an error occurs.
#'
#' @examples
#' # Assuming 'my_data' is your data frame or tibble
#' format_and_print_table(my_data)
#' format_and_print_table(my_data, zero_print = "0", max_rows = 100)
#'
#' @export
format_and_print_table <- function(data, zero_print = ".", max_rows = 1000, ...) {
  # Check if data is a data frame or tibble
  if (!is.data.frame(data)) {
    stop("Input must be a data frame or tibble.")
  }
  
  tryCatch({
    # Format the table
    formatted_table <- tryCatch({
      format(data, zero.print = zero_print, ...)
    }, error = function(e) {
      cat("Error call in: '", deparse(conditionCall(e)), "'\n")
      stop(cat(
        "Error caused no formatted table to be created:\n",
        conditionMessage(e), "\n"
      ))
    })
    
    # Print the formatted table
    tryCatch({
      print(formatted_table, n = max_rows)
    }, error = function(e) {
      cat("Error call in: '", deparse(conditionCall(e)), "'\n")
      warning(cat(
        "Error caused formatted table not to be printed:\n",
        conditionMessage(e), "\n"
      ))
    })
    
    # Return the formatted table invisibly
    invisible(formatted_table)
  }, error = function(e) {
    cat("Error call in: '", deparse(conditionCall(e)), "'\n")
    cat("Unexpected error:\n", conditionMessage(e), "\n")
    invisible(NULL)
  })
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

      if ( nrow(counts_df_T_lng_AUID_deduped) == 0 || is.null(counts_df_T_lng_AUID_deduped)) {
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
    }#,
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
    }#,
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
    }#,
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
    }#,
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
    }#,
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

  # tryCatch(
  #   {
  #     query_df <- dedup_by_group(cmap_coloc, group_id)},
  #       error = function(e) {
  #         cat("Error making query_df:", conditionMessage(e), "\n")
  #       }
  # )

  # cat(query_df)

  tryCatch(
    {
    query_df <- dedup_by_group(cmap_coloc, group_id_key = unique_sample_id_key, group_id)

  #   cat("Removing duplicates based on group ID...\n")
  #   query_df <- cmap_coloc %>%
  #     add_group_id(unique_sample_id_key) %>%
  #     distinct(group_id, .keep_all = TRUE)
  # # }
  #   rows_removed <- nrow(cmap_coloc) - nrow(query_df)
    
  #   cat("Number of rows removed:", rows_removed, "\n")

    # cat(query_df)

    ## - nucleicAcidType
    samples_per_nucleicAcidType <- count_and_arrange(query_df, "nucleicAcidType") %>%
          add_percentage(n,
            percentage,
            grouping_by = NULL,
            remove_columns = "total"
          ) %>%
          # mutate(
          #   total = sum(n),
          #   percentage = n / total * 100
          # ) %>%
          # select(nucleicAcidType, n, percentage)
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
      )  #%>% 
      # add_total_row(
      #   column_name = "photic",
      #   summary_column = n,
      #   all_columns = FALSE
      # )

          print(samples_per_photic, n = 1000)
        },
        error = function(e) {
          cat("Error analyzing photic samples:", conditionMessage(e), "\n")
        }#,
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
        }#,
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

  # Load other files and scripts needed
  source_file(files_to_source)

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

  return(stats_list)
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

  # print(final_results)

  if (!is.null(final_results)) {
     # Create output directory if it doesn't exist
    create_dir(args$files_out_path)
    write_file_list(final_results, args$files_out_path)
  }

  # if (!is.null(final_results))
  #   done_flag <-list("BSS script finished!")
  #   names(done_flag) <- "BSS_script_finished"
  #   write_file_list(done_flag, args$files_out_path)

}
