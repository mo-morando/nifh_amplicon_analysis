#!/usr/bin/env Rscript

#' @title NifH Amplicon Data Preprocessing Pipeline
#' @description This script preprocesses NifH amplicon data by loading input files, performing data transformations, and generating processed outputs for downstream analysis. It provides a comprehensive and modular approach to handle complex NifH amplicon sequencing data, with robust error handling and informative logging.
#'
#' @details The pipeline performs the following steps:
#' * Sources required utility functions and plotting scripts
#' * Parses command-line arguments for input files, input path, and output path
#' * Loads and processes annoNifHDB update data
#' * Processes metaTab data and merges with CMAP data
#' * Cleans and transforms CMAP dataset
#' * Processes sample types, including DNA/RNA splitting and deduplication
#' * Constructs new tibbles for counts and relative abundance data
#' * Writes processed data to output files
#'
#' Key functions include:
#' * source_file(): Sources external R scripts with error handling
#' * setup_parser(): Configures command-line argument parsing
#' * parse_arg(): Parses command-line arguments
#' * process_annoNifHDB_updt(): Processes annoNifHDB update data
#' * process_metaTab(): Processes metaTab data and merges with CMAP data
#' * cmap_clean_main(): Cleans and transforms CMAP dataset
#' * process_sample_types(): Processes sample types and creates keys
#' * construct_new_tibbles(): Constructs new tibbles for further analysis
#' * main(): Orchestrates the entire data preprocessing workflow
#'
#' @usage Rscript preprocess_files.R [--files FILES] [--input_path PATH] [--output_path PATH]
#'
#' @param --files Comma-separated list of files to read in
#' @param --input_path Input directory path
#' @param --output_path Output directory path
#'
#' @author Michael Morando
#'  2023
#'
#' @note This script requires the following R packages: tidyverse, argparser
#'
#' @examples
#' Rscript preprocess_files.R --files annoNifHDB_updt,metaTab,cmap_coloc,nifhDB_cnts,nifhDB_RA --input_path ../analysis/out_files --output_path ../analysis/processed_files
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
  parser <- arg_parser("Process the workspace files")
  parser <- add_argument(parser, "--files",
    help = "CVS list of files to read in",
    default = "annoNifHDB_updt,metaTab,cmap_coloc,nifhDB_cnts,nifhDB_RA"
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

# TODO: Add tryCatch error handling to all functions
# Including nested tryCatch for various manipulations within each function
#' Process annoNifHDB Update
#'
#' This function processes the annoNifHDB update data, performing various
#' transformations and cleaning steps to prepare it for further analysis.
#'
#' @param data A data frame containing annoNifHDB update data.
#' @return A processed data frame ready for further analysis.
process_annoNifHDB_updt <- function(data) {
  cat("Processing the annoNifHDB_updt files...\n")

  clusters_of_interest <- c(
    "1A",
    "1J/1K",
    "1O/1P",
    "3",
    "1G",
    "1B",
    "4",
    "2"
  )

  processed_data <- data %>%
    mutate(
      CON = primary_id,
      nifH_cluster = case_when(
        cluster %in% c(3, 4, 2) ~ as.character(cluster),
        is.na(subcluster) ~ as.character(cluster),
        TRUE ~ subcluster
      ),
      nifH_cluster = case_when(
        grepl("UCYN-B|croco", UCYNAoligos.id, ignore.case = TRUE) ~ "1B",
        nifH_cluster %in% c("1J", "1K") ~ "1J/1K",
        nifH_cluster %in% c("1P", "1O") ~ "1O/1P",
        is.na(nifH_cluster) | subcluster == "nan" ~ "unknown",
        TRUE ~ as.character(nifH_cluster)
      ),
      cluster_stats = ifelse(test = nifH_cluster %in% clusters_of_interest,
        yes = nifH_cluster,
        no = "other"
      ),
      group1 = if_else(nifH_cluster %in% c("1A", "3"), "1A C3", nifH_cluster),
      group2 = if_else(nifH_cluster %in% c("1A", "3", "1J/1K"),
        "1A,  3,  &  1J/1K",
        nifH_cluster
      ),
      group3 = if_else(nifH_cluster %in% c("1A", "3"), "1A C3",
        if_else(nifH_cluster %in% c("1G", "1J/1K"), "1B 1J/1K", nifH_cluster)
      ),
      group4 = if_else(nifH_cluster %in% c("1A", "3", "1J/1K"),
        "1A/C3 1J/1K", "1B & 1G"
      ),
      CyanoCON = if_else(nifH_cluster == "1B", CON, nifH_cluster),
      crocoCMB = if_else(
        grepl("croco", CON, ignore.case = TRUE) &
          CON == "Crocosphera_DQ118216_Moisander", "CDQmois",
        if_else(grepl("croco", CON, ignore.case = TRUE), "CrocoCMB", CON)
      ),
      CyanoGroups = if_else(CON == "UCYN-A3" | CON == "UCYN-A4", "A3-A4", CON),
      CyanoGroupsII = if_else(grepl("Tricho", CON), "Trichodesmium sp.", CON),
      CyanoGroupsIII = if_else(grepl("Tricho", CON),
        "Trichodesmium sp.", CyanoGroups
      ),
      CyanoGroupsIV = if_else(CON == "UCYN-A3" | CON == "UCYN-A1", "A1-A3",
        if_else(CON == "UCYN-A2" | CON == "UCYN-A4", "A2-A4", CON)
      )
    )

  cat("process_annoNifHDB_updt is done!\n")


  return(processed_data)
}


#' Process metaTab Data
#'
#' This function processes the metaTab data, including fixing up size fractions
#' before merging with other metadata.
#'
#' @param metaTab A data frame containing metaTab data.
#' @param cmap_coloc A data frame containing CMAP data with sample collection data.
#' @return A merged data frame (cmap_coloc) with updated size fractions.
process_metaTab <- function(metaTab, cmap_coloc) {
  cat("Processing the metaTab files...\n")

  processed_metaTab <- metaTab %>%
    mutate(
      Size_fraction =
        if_else(Size_fraction %in%
          c("whole", "Sterivex", "0.22") | is.na(Size_fraction),
        "0.2",
        Size_fraction
        ),
      coastal_class = if_else(
        Coastal_200km == TRUE,
        true = "coastal",
        false = "open ocean"
      ),
      StudyID = case_when(
        StudyID %in% "Turk_2021" ~ "TurkKubo_2021",
        StudyID %in% "TianjUni_2016" ~ "Wu_2021",
        StudyID %in% "TianjUni_2017" ~ "Wu_2019",
        .default = StudyID
      )
    )


  merged_data <- cmap_coloc %>%
    left_join(processed_metaTab, by = "SAMPLEID")


  cat("Processing the metaTab files is done!\n")


  return(merged_data)
}


#' Clean up CMAP dataset
#'
#' This function cleans up the cmap_coloc dataset by performing various data processing steps.
#'
#' @param cmap_coloc The CMAP dataset to be cleaned.
#' @return A list containing the cleaned CMAP dataset and the photic samples key.
cmap_clean_main <- function(cmap_coloc) {
  cat("Cleaning up cmap_coloc...\n")

  # Function to mutate studyID based on conditions
  fix_studyid <- function(data) {
    data %>%
      mutate(
        studyID = case_when(
          studyID == "Gradoville_2020" & grepl("2017", date)
          ~ "Gradoville_2020_G2",
          studyID == "Gradoville_2020" & grepl("2016", date)
          ~ "Gradoville_2020_G1",
          studyID %in% "Turk_2021" ~ "TurkKubo_2021",
          studyID %in% "TianjUni_2016" ~ "Wu_2021",
          studyID %in% "TianjUni_2017" ~ "Wu_2019",
          .default = studyID
        )
      )
  }

  # Function to replace depth zero values with a surface depth
  replace_zero_depth <- function(depth) {
    ifelse(depth == 0, 2, depth)
  }

  # Function to add DNA/RNA column
  add_nucleic_acid_type <- function(LibrarySource) {
    if_else(LibrarySource == "METATRANSCRIPTOMIC", "RNA", "DNA")
  }

  # Function to standardize Size_fraction
  standardize_size_fraction <- function(Size_fraction) {
    if_else(Size_fraction %in% c("whole", "Sterivex", "0.22") | is.na(Size_fraction), "0.2", Size_fraction)
  }

  # Function to determine hemisphere
  determine_hemisphere <- function(lat) {
    hemi <- cut(lat,
      breaks = c(-Inf, 0, Inf),
      labels = c("southernHemi", "northernHemi")
    )
    factor(hemi, c("northernHemi", "southernHemi"))
  }

  # Function to determine season
  determine_season <- function(month, hemi) {
    case_when(
      (month %in% c("03", "04", "05") & hemi == "northernHemi") ~ "spring",
      (month %in% c("06", "07", "08") & hemi == "northernHemi") ~ "summer",
      (month %in% c("09", "10", "11") & hemi == "northernHemi") ~ "fall",
      (month %in% c("12", "01", "02") & hemi == "northernHemi") ~ "winter",
      (month %in% c("03", "04", "05") & hemi == "southernHemi") ~ "fall",
      (month %in% c("06", "07", "08") & hemi == "southernHemi") ~ "winter",
      (month %in% c("09", "10", "11") & hemi == "southernHemi") ~ "spring",
      (month %in% c("12", "01", "02") & hemi == "southernHemi") ~ "summer",
      TRUE ~ NA
    )
  }

  # Function to determine geographical region
  determine_geo_region <- function(lat_abs) {
    geoRegion <- cut(lat_abs,
      breaks = c(-1, 23, 35, 66, 100),
      labels = c("Eq/Trop", "Sub-Trop", "Temp", "Poles")
    )
    factor(geoRegion, c("Eq/Trop", "Sub-Trop", "Temp", "Poles"))
  }

  # Function to rename specific columns in the CMAP tibble
  rename_cmap_columns <- function(column_name) {
    column_name <- str_remove(column_name, "CMAP_")
    column_name <- coalesce(str_extract(column_name, "^[A-Za-z0-9]+_darwin"), column_name)
    column_name <- str_replace(column_name, "clim_tblWOA_2018_qrtdeg_Climatology", "WOA_2018_qrtdeg")
    column_name <- str_replace(column_name, "clim_tblWOA_2018_1deg_Climatology", "WOA_2018_1deg")
    column_name <- str_replace(column_name, "WOA_clim_tblWOA_Climatology", "WOA_clim")
    column_name <- str_replace(column_name, "clim_tblWOA_2018_MLD_qrtdeg_Climatology", "WOA_2018_MLD_qrtdeg")
    column_name <- str_replace(column_name, "^C_", "conductivity_")
    column_name <- str_replace(column_name, "^i_", "sigma_")
    column_name <- str_replace(column_name, "^t_", "temp_")
    column_name <- str_replace(column_name, "^A_", "AOU_")
    column_name <- str_replace(column_name, "^O_", "O2_sat_")
    column_name <- str_replace(column_name, "^n_", "NO3_")
    column_name <- str_replace(column_name, "^p_", "phosphate_")
    column_name <- str_replace(column_name, "^si_", "silica_")
    column_name <- str_replace(column_name, "^M_", "MLD_")
    column_name <- str_replace(column_name, "^s_", "salinity_")
    column_name
  }

  # # Read coastal/open ocean identification data
  # coastal_open_ids <- read_csv("../data/workspace/coastal_openocean_ids/sample_classifications_in_nifH_ASV_DB.csv", show_col_types = FALSE) %>%
  #   rename(
  #     coastal_class = classification,
  #     studyID = StudyID
  #   )

  studyid_regions <- read_csv("../data/misc/studyid_regions.csv", show_col_types = FALSE)

  # Define function to add study region to a tibble
  add_study_region <- function(tibble, key = studyid_regions) {
    tb_sr <- tibble %>%
      left_join(key) %>%
      mutate(
        ocean = if_else(lat_abs >= 60 & hemi == "southernHemi",
          "Southern", study_ocean
        )
      )


    return(tb_sr)
  }


  # Define function to add study region to a tibble
  photic_type <- function(coastal_class, depth) {
    return(case_when(
      (coastal_class %in% "open ocean" & depth <= 100) ~ TRUE,
      (coastal_class %in% "coastal" & depth <= 50) ~ TRUE,
      TRUE ~ FALSE
    ))
  }


  # Perform data cleaning and transformation steps
  cmap_coloc <- cmap_coloc %>%
    # left_join(coastal_open_ids) %>%
    rename_with(rename_cmap_columns) %>%
    mutate(
      depth = replace_zero_depth(depth),
      hemi = determine_hemisphere(lat),
      lat_abs = abs(lat),
      date = str_remove(time, " .+$"),
      month = str_remove_all(str_extract(date, "-.+-"), "-"),
      year = str_extract(date, "^\\d{4}"),
      season = determine_season(month, hemi),
      season = factor(season, c("winter", "spring", "summer", "fall")),
      geoRegion = determine_geo_region(lat_abs),
      logFe = log(Fe_tblPisces_NRT)
    ) %>%
    fix_studyid() %>% # ! date is required for fix_studyid function so must
    # ! go after this column is made
    # Add study regions to CMAP
    # Currently, only ocean region
    add_study_region() %>% # ! fixed studyids are required
    mutate(
      Size_fraction = standardize_size_fraction(Size_fraction),
      nucleicAcidType = add_nucleic_acid_type(LibrarySource),
      photic = photic_type(coastal_class, depth)
      # photic = case_when(
      #   (coastal_class %in% "open ocean" & depth <= 100) ~ TRUE,
      #   (coastal_class %in% "coastal" & depth <= 50) ~ TRUE,
      #   TRUE ~ FALSE
      # )
    )

  # Define photic and aphotic samples
  ## Make photic sample key
  photic_samples_key <- cmap_coloc %>%
    filter(photic == TRUE) %>%
    pull(SAMPLEID)


  cat("Cleaning up cmap_coloc is done!\n")

  return(list(
    cmap_coloc = cmap_coloc,
    photic_samples_key = photic_samples_key
  ))
}


### - We have to have a tibble that identifies the different sample types and nucleic acid types
#### - To do this we need to split the tibble into a DNA and RNA
#### - Then identify the samples types
#### - Average by sample type and nucleic acid type
#### - Deduplicate by group

#' Process Sample Types
#'
#' This function processes the sample types in the CMAP data by splitting them
#' into DNA and RNA, identifying sample types, averaging by sample type and
#' nucleic acid type, and deduplicating by group.
#'
#' @param cmap_coloc A data frame containing CMAP data with sample collection data.
#' @return A list containing processed data frames and keys.
process_sample_types <- function(cmap_coloc) {
  cat("Define some keys and flags to calculate some basic sample stats and process some of the tibbles for downstream analysis...\n")
  # Define columns representing a sampling point
  sampling_point <- c(
    "studyID",
    "lat",
    "lon",
    "time",
    "depth"
  )

  # Assign a unique sample point ID and count distinct sample points
  cmap_coloc <- cmap_coloc %>%
    group_by(across(all_of(sampling_point))) %>%
    mutate(
      sample_point = cur_group_id(), # Assign a unique ID for each sampling point
      num_dist_samp_pnts = n() # Count distinct sampling points
    )

  # Identify size fractions and create a flag for one or two size fractions
  size_fraction_key <- cmap_coloc %>%
    group_by(across(all_of(sampling_point))) %>%
    summarise(num_distinct_size_fractions = n_distinct(Size_fraction)) %>%
    mutate(
      size_frac_flag = case_when(
        (num_distinct_size_fractions == 1) ~ "One_Size_Fractions",
        (num_distinct_size_fractions > 1) ~ "Two_Size_Fractions",
      )
    ) %>%
    ungroup() %>%
    suppressMessages()

  # Add size fraction flag to the main data frame
  cmap_coloc <- add_size_frac_key(
    cmap_coloc,
    size_fraction_key,
    num_distinct_size_fractions
  )

  # Identify DNA samples
  DNA_samples_key <- cmap_coloc %>%
    filter(nucleicAcidType == "DNA") %>%
    pull(SAMPLEID)

  # Split DNA samples
  sample_types_DNA <- remove_samples_nucleic_acid(
    cmap_coloc,
    "DNA",
    DNA_samples_key
  ) %>%
    add_rep_flag() %>%
    select(sample_point, num_dist_samp_pnts, nucleicAcidType, size_frac_flag, Size_fraction, replicate_flag, everything()) %>%
    ungroup()

  # Identify RNA samples
  sample_types_RNA <- remove_samples_nucleic_acid(
    cmap_coloc,
    "RNA",
    DNA_samples_key
  ) %>%
    add_rep_flag() %>%
    select(sample_point, num_dist_samp_pnts, nucleicAcidType, size_frac_flag, Size_fraction, replicate_flag, everything()) %>%
    ungroup()

  # Combine DNA and RNA samples
  sample_types_all <- sample_types_DNA %>%
    bind_rows(sample_types_RNA)

  # Define columns for creating unique sample IDs
  unique_sample_column_ids <- c(
    "studyID",
    "lat",
    "lon",
    "time",
    "depth",
    "replicate_flag",
    "size_frac_flag",
    "nucleicAcidType"
  )

  # Assign a unique ID for each sample group
  cmap_coloc <- sample_types_all %>%
    group_by(across(all_of(unique_sample_column_ids))) %>%
    mutate(group_id = as.character(cur_group_id())) %>%
    ungroup()

  # Create a key for unique sample IDs
  unique_sample_id_key <- cmap_coloc %>%
    select(all_of(unique_sample_column_ids), SAMPLEID, group_id) %>%
    distinct(SAMPLEID, group_id, .keep_all = TRUE)

  cat("Done defining keys and some calculations\n")

  return(list(
    size_fraction_key = size_fraction_key,
    cmap_coloc = cmap_coloc,
    DNA_samples_key = DNA_samples_key,
    sample_types_all = sample_types_all,
    unique_sample_column_ids = unique_sample_column_ids,
    unique_sample_id_key = unique_sample_id_key
  ))
}


#' Construct New Tibbles for Other Types of Analysis
#'
#' This function constructs new tibbles for other types of analysis by transforming and processing the data.
#'
#' @param nifhDB_cnts A data frame containing nifH count data.
#' @param nifhDB_RA A data frame containing nifH relative abundance data.
#' @param unique_sample_id_key A data frame containing unique sample ID information.
#' @return A list containing the constructed tibbles.
construct_new_tibbles <- function(nifhDB_cnts, nifhDB_RA, unique_sample_id_key) {
  cat("Constructing new tibbles for other types of analysis...\n")

  # Internal function to process data
  process_nifhDB_data <- function(data, is_cnts, grp_key) {
    # Transform data
    data_T <- data %>%
      pivot_longer(cols = -AUID, names_to = "SAMPLEID", values_to = "Value") %>%
      pivot_wider(names_from = AUID, values_from = Value)

    value_col <- if (is_cnts) "counts" else "RA"

    # Transform data to long format for averaging and deduplication
    data_T_lng <- transform_data_lng(
      input_df = data_T,
      starts_with_col = "AUID",
      names_to_col = "AUID",
      values_to_col = value_col
    )


    if (is_cnts) {
      df_T_lng_deduped <- main_cnts_dedup_by_group(
        df_lng = data_T_lng,
        AUID, group_id,
        grp_key = grp_key
      )
      #-## now pivot wide #-## now pivot wide
      df_T_deduped <- df_T_lng_deduped %>%
        select(-all_of(c("group_id"))) %>%
        pivot_wider(
          names_from = AUID,
          values_from = counts
        )
      #-## now pivot wide and transform
      df_deduped <- df_T_lng_deduped %>%
        select(-all_of(c("group_id"))) %>%
        pivot_wider(
          names_from = SAMPLEID,
          values_from = counts
        )

      return(list(
        data_T_lng = data_T_lng,
        df_T_lng_deduped = df_T_lng_deduped,
        df_T_deduped = df_T_deduped,
        df_deduped = df_deduped
      ))
    } else {
      # Average over groups and deduplicate data if RA data
      df_T_lng_mean_deduped <- main_average_ra_dedup_by_group(
        df_lng = data_T_lng,
        AUID, group_id,
        grp_key = grp_key,
        mean_by = !!sym(value_col)
      )

      #-## now pivot wide #-## now pivot wide
      df_T_mean_deduped <- df_T_lng_mean_deduped %>%
        select(-any_of(c("RA", "group_id"))) %>%
        pivot_wider(
          names_from = AUID,
          values_from = average_value_AUID
        )
      #-## now pivot wide and transform
      df_mean_deduped <- df_T_lng_mean_deduped %>%
        select(-all_of(c("RA", "group_id"))) %>%
        pivot_wider(
          names_from = SAMPLEID,
          values_from = average_value_AUID
        )

      return(list(
        data_T_lng = data_T_lng,
        df_T_lng_mean_deduped = df_T_lng_mean_deduped,
        df_T_mean_deduped = df_T_mean_deduped,
        df_mean_deduped = df_mean_deduped
      ))
    }
  }

  # Process nifH count data
  nifhDB_cnts_results <- process_nifhDB_data(nifhDB_cnts, is_cnts = TRUE, grp_key = unique_sample_id_key)
  cat("Counts done!!\n\n")

  # Process nifH relative abundance data
  nifhDB_RA_results <- process_nifhDB_data(nifhDB_RA, is_cnts = FALSE, grp_key = unique_sample_id_key)

  cat("Relative abundace done!!\n\n")

  # Return the constructed tibbles as a list
  return(c(nifhDB_cnts_results, nifhDB_RA_results))
}

#' Main function to process the data
#'
#' @param files_to_read List of files to process
#' @param files_in_path Input file path
#' @param files_out_path Output file path
#' @return A list of processed data frames
main <- function(files_to_read, files_in_path, files_out_path) {
  cat("preprocess_files.R script executed!!\n")

  # Load the data
  data_list <- load_files(files_to_read, files_in_path)
  # data_list <- load_files(c("annoNifHDB_updt","metaTab","cmap_coloc","nifhDB_cnts","nifhDB_RA"), "analysis/out_files")

  # Process data

  annoNifHDB_updt <- process_annoNifHDB_updt(data_list$annoNifHDB_updt)

  cmap_coloc <- process_metaTab(data_list$metaTab, data_list$cmap_coloc)

  result <- cmap_clean_main(cmap_coloc)
  cmap_coloc <- result$cmap_coloc
  photic_samples_key <- result$photic_samples_key

  result_list <- process_sample_types(cmap_coloc)

  new_tibbles <- construct_new_tibbles(data_list$nifhDB_cnts, data_list$nifhDB_RA, result_list$unique_sample_id_key)

  # Prepare the final results
  final_results <- list(
    annoNifHDB_updt = annoNifHDB_updt,
    cmap_coloc = cmap_coloc,
    photic_samples_key = photic_samples_key,
    size_fraction_key = result_list$size_fraction_key,
    DNA_samples_key = result_list$DNA_samples_key,
    sample_types_all = result_list$sample_types_all,
    unique_sample_column_ids = result_list$unique_sample_column_ids,
    unique_sample_id_key = result_list$unique_sample_id_key,
    counts_df_T_lng = new_tibbles[[1]],
    counts_df_T_lng_AUID_deduped = new_tibbles$df_T_lng_deduped,
    counts_df_T_AUID_deduped = new_tibbles$df_T_deduped,
    counts_df_AUID_deduped = new_tibbles$df_deduped,
    RA_df_T_lng = new_tibbles[[5]],
    RA_df_T_lng_mean_RA_AUID_deduped = new_tibbles$df_T_lng_mean_deduped,
    RA_df_T_mean_RA_AUID_deduped = new_tibbles$df_T_mean_deduped,
    RA_df_mean_RA_AUID_deduped = new_tibbles$df_mean_deduped
  )

  if (!is.null(final_results)) {
    create_dir(files_out_path)
    write_file_list(
      file_list = final_results,
      path = files_out_path
    )
  }

  cat("\nPreprocessing done.\n\n")
  # return(final_results)
}


# Run if the script is being run directly
if (sys.nframe() == 0 && !interactive()) {
  source_file(files_to_source)


  parser <- setup_parser()
  args <- parse_arg(parser)


  validate_parsed_args(parsed_args = args)

  main(
    args$files_to_read,
    args$files_in_path,
    args$files_out_path
  )
}
