# Check if the argparser library is available and install it if necessary
if (!requireNamespace("argparser", quietly = TRUE)) {
  install.packages("argparser")
}

library(tidyverse)
library(argparser)

source("/Users/mo/Projects/nifH_amp_project/myWork/scripts/functions.R")
source("/Users/mo/Projects/nifH_amp_project/myWork/scripts/basic_plotting.R")

## - functions
cat("

Functions needed:
transform_data_lng
add_group_id
dedupe_by_group
remove_samples_nucleic_acid
count_and_arrange

")

#### -

# _# Processing the files of the workspace
cat("\nProcessing the files of the workspace\n")


files_in_path <- "analysis/out_files"

files_out_path <- "analysis/out_files"

files_to_read <- c(
  "annoNifHDB_updt",
"metaTab",
"cmap_coloc",
"nifhDB_cnts",
"nifhDB_RA")

# for (file in files_to_read) {
#   cat("Loading file:", file.path(files_in_path, paste0(file, ".csv")), "\n")
#   assign(file, read_csv(file.path(files_in_path, paste0(file, ".csv"))))
# }

load_files <- function(file_list, path) {
   for (file in file_list) {
    cat("Loading file:", file.path(path, paste0(file, ".csv")), "\n")
    assign(file, read_csv(file.path(path, paste0(file, ".csv"))), envir = .GlobalEnv)
    }
}

load_files(files_to_read, files_in_path)

cat(ls())
print(annoNifHDB_updt)

#-# Processing annoNifHDB_updt



#' Process annoNifHDB Update
#'
#' This function processes the annoNifHDB update data, performing various transformations
#' and cleaning steps to prepare it for further analysis.
#'
#' @param data A data frame containing annoNifHDB update data.
#' @return A processed data frame ready for further analysis.
#' @export
#' @examples
#' # Load the data
#' data("annoNifHDB_updt")
#'
#' # Process the data
#' processed_data <- process_annoNifHDB_updt(annoNifHDB_updt)
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
      CON = consensus_id,
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

  return(processed_data)
}


# _##########################################################
#-# Processing metaTab
#' Process metaTab Data
#'
#' This function processes the metaTab data, including fixing up size fractions
#' before merging with other metadata.
#'
#' @param metaTab A data frame containing metaTab data.
#' @param cmap_coloc A data frame containing CMAP data with sample collection data.
#' @return A merged data frame (cmap_coloc) with updated size fractions.
#' @export
#' @examples
#' # Load the data
#' data("metaTab")
#' data("cmap_coloc")
#'
#' # Process the data
#' processed_data <- process_metaTab(metaTab, cmap_coloc)
process_metaTab <- function(metaTab, cmap_coloc) {
  cat("Processing the metaTab files...\n")

  processed_metaTab <- metaTab %>%
    mutate(
      Size_fraction =
        if_else(Size_fraction %in%
          c("whole", "Sterivex", "0.22") | is.na(Size_fraction),
        "0.2",
        Size_fraction
        )
    )

  merged_data <- cmap_coloc %>%
    left_join(processed_metaTab, by = "SAMPLEID")

  return(merged_data)
}
# _##########################################################


# _##########################################################
#-# Processing cmap_coloc
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
          studyID == "Gradoville_2020" & grepl("2017", date) ~ "Gradoville_2020_G2",
          studyID == "Gradoville_2020" & grepl("2016", date) ~ "Gradoville_2020_G1",
          TRUE ~ studyID
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

  # Read coastal/open ocean identification data
  coastal_open_ids <- read_csv("data/workspace/coastal_openocean_ids/sample_classifications_in_nifH_ASV_DB.csv") %>%
    rename(
      coastal_class = classification,
      studyID = StudyID
    )

  #-# add Ocean regions for each study ID
  studyid_regions <- read_csv("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/studyid_regions.csv")

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
    left_join(coastal_open_ids) %>%
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

  return(list(cmap_coloc = cmap_coloc, photic_samples_key = photic_samples_key))
}


### _ Define some keys and process some of the tibbles for downstream analysis
# cat("Load in the keys 🔑🗝️ ")
cat("Defining keys...")


### - We have to have a tibble that identifies the different sample types and nucleic acid types
#### - To do this we need to split the tibble into a DNA and RNA
#### - Then identify the samples types
#### - Average by sample type and nucleic acid type
#### - Deduplicate by group

#' Process Sample Types
#'
#' This function processes the sample types in the CMAP data by splitting them into DNA and RNA,
#' identifying sample types, averaging by sample type and nucleic acid type, and deduplicating by group.
#'
#' @param cmap_coloc A data frame containing CMAP data with sample collection data.
#' @return A processed data frame with sample types identified, averaged, and deduplicated.
#' @export
#' @examples
#' # Load the data
#' data("cmap_coloc")
#'
#' # Process the data
#' processed_data <- process_sample_types(cmap_coloc)
process_sample_types <- function(cmap_coloc) {
  cat("Define some keys and flags to calculate some basic sample stats and process some of the tibbles for downstream analysis...")
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
    ungroup()

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

  return(list(
    size_fraction_key = size_fraction_key,
    cmap_coloc = cmap_coloc,
    DNA_samples_key = DNA_samples_key,
    sample_types_all = sample_types_all,
    unique_sample_column_ids = unique_sample_column_ids,
    unique_sample_id_key = unique_sample_id_key
  ))
}



# _# Construct some new tibbles for other types of analysis
#' Construct New Tibbles for Other Types of Analysis
#'
#' This function constructs new tibbles for other types of analysis by transforming and processing the data.
#'
#' @param nifhDB_cnts A data frame containing nifH count data.
#' @param nifhDB_RA A data frame containing nifH relative abundance data.
#' @return A list containing the constructed tibbles.
#' @export
#' @examples
#' # Load the data
#' data("nifhDB_cnts")
#' data("nifhDB_RA")
#'
#' # Construct new tibbles
#' new_tibbles <- construct_new_tibbles(nifhDB_cnts, nifhDB_RA)
construct_new_tibbles <- function(nifhDB_cnts, nifhDB_RA, unique_sample_id_key) {
  cat("Constructing new tibbles for other types of analysis...\n")
  
  # # Initialize the variables
  # df_T_lng_deduped <- NULL
  # df_T_deduped <- NULL
  # df_deduped <- NULL
  # df_T_lng_mean_deduped <- NULL
  # df_T_mean_deduped <- NULL
  # df_mean_deduped <- NULL

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
      print(head(df_T_lng_mean_deduped))
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

  # # Initialize the variables
  # df_T_lng_deduped <- NULL
  # df_T_deduped <- NULL
  # df_deduped <- NULL
  # df_T_lng_mean_deduped <- NULL
  # df_T_mean_deduped <- NULL
  # df_mean_deduped <- NULL

  # Process nifH count data
  nifhDB_cnts_results <- process_nifhDB_data(nifhDB_cnts, is_cnts = TRUE, grp_key = unique_sample_id_key)
  cat("Counts done!!\n\n")

  # Process nifH relative abundance data
  nifhDB_RA_results <- process_nifhDB_data(nifhDB_RA, is_cnts = FALSE, grp_key = unique_sample_id_key)

  cat("Relative abundace done!!\n\n")

  # Return the constructed tibbles as a list
  return(c(nifhDB_cnts_results, nifhDB_RA_results))
}


# Get the names of all objects in the workspace
# workspace_objects_final <- ls()

# files_out_path <- "analysis/files/"

# write_workspace_to_csv(workspace_objects_final, )


main <- function() {
  print("preprocess_files.R script executed!!")
  annoNifHDB_updt <- process_annoNifHDB_updt(annoNifHDB_updt)
  print("process_annoNifHDB_updt is done")

  cmap_coloc <- process_metaTab(metaTab, cmap_coloc)
  print("process_metaTab is done")

  result <- cmap_clean_main(cmap_coloc)
  cmap_coloc <- result$cmap_coloc
  photic_samples_key <- result$photic_samples_key
  print("cmap_clean_main is done")

  # Call the process_sample_types function
  result_list <- process_sample_types(cmap_coloc)
  # Access elements from the returned list by their names
  size_fraction_key <- result_list$size_fraction_key
  cmap_coloc <- result_list$cmap_coloc
  DNA_samples_key <- result_list$DNA_samples_key
  sample_types_all <- result_list$sample_types_all
  unique_sample_column_ids <- result_list$unique_sample_column_ids
  unique_sample_id_key <- result_list$unique_sample_id_key

  new_tibbles <- construct_new_tibbles(nifhDB_cnts, nifhDB_RA, unique_sample_id_key)
  # Extract the constructed tibbles
  # count samples
  counts_df_T_lng <- new_tibbles[[1]]
  counts_df_T_lng_AUID_deduped <- new_tibbles$df_T_lng_deduped
  counts_df_T_AUID_deduped <- new_tibbles$df_T_deduped
  counts_df_AUID_deduped <- new_tibbles$df_deduped
  # RA samples
  RA_df_T_lng <- new_tibbles[[5]]
  RA_df_T_lng_mean_RA_AUID_deduped <- new_tibbles$df_T_lng_mean_deduped
  RA_df_T_mean_RA_AUID_deduped <- new_tibbles$df_T_mean_deduped
  RA_df_mean_RA_AUID_deduped <- new_tibbles$df_mean_deduped
  print("construct_new_tibbles is done")
  
  print("Script is done")

  # # Return the variables you want to access outside of this function
  # return(list(
  #   counts_df_T_lng_AUID_deduped = counts_df_T_lng_AUID_deduped,
  #   counts_df_T_AUID_deduped = counts_df_T_AUID_deduped,
  #   counts_df_AUID_deduped = counts_df_AUID_deduped,
  #   RA_df_T_lng_mean_RA_AUID_deduped = RA_df_T_lng_mean_RA_AUID_deduped,
  #   RA_df_T_mean_RA_AUID_deduped = RA_df_T_mean_RA_AUID_deduped,
  #   RA_df_mean_RA_AUID_deduped = RA_df_mean_RA_AUID_deduped,
  # ))
}


main <- function() {
  print("preprocess_files.R script executed!!")
  annoNifHDB_updt <- process_annoNifHDB_updt(annoNifHDB_updt)
  print("process_annoNifHDB_updt is done")

  cmap_coloc <- process_metaTab(metaTab, cmap_coloc)
  print("process_metaTab is done")

  result <- cmap_clean_main(cmap_coloc)
  cmap_coloc <- result$cmap_coloc
  photic_samples_key <- result$photic_samples_key
  print("cmap_clean_main is done")

  # Call the process_sample_types function
  result_list <- process_sample_types(cmap_coloc)


  # Call the process_sample_types function
  result_list <- process_sample_types(cmap_coloc)
  # Access elements from the returned list by their names
  size_fraction_key <- result_list$size_fraction_key
  cmap_coloc <- result_list$cmap_coloc
  DNA_samples_key <- result_list$DNA_samples_key
  sample_types_all <- result_list$sample_types_all
  unique_sample_column_ids <- result_list$unique_sample_column_ids
  unique_sample_id_key <- result_list$unique_sample_id_key

  new_tibbles <- construct_new_tibbles(nifhDB_cnts, nifhDB_RA, unique_sample_id_key)

  # print(new_tibbles)

  # Access elements from the returned list by their names
  return(list(
    annoNifHDB_updt = annoNifHDB_updt,
    cmap_coloc = cmap_coloc,
    photic_samples_key = photic_samples_key,
    size_fraction_key = result_list$size_fraction_key,
    DNA_samples_key = result_list$DNA_samples_key,
    sample_types_all = result_list$sample_types_all,
    unique_sample_column_ids = result_list$unique_sample_column_ids,
    unique_sample_id_key = unique_sample_id_key,
    counts_df_T_lng = new_tibbles[[1]],
    counts_df_T_lng_AUID_deduped = new_tibbles$df_T_lng_deduped,
    counts_df_T_AUID_deduped = new_tibbles$df_T_deduped,
    counts_df_AUID_deduped = new_tibbles$df_deduped,
    RA_df_T_lng = new_tibbles[[5]],
    RA_df_T_lng_mean_RA_AUID_deduped = new_tibbles$df_T_lng_mean_deduped,
    RA_df_T_mean_RA_AUID_deduped = new_tibbles$df_T_mean_deduped,
    RA_df_mean_RA_AUID_deduped = new_tibbles$df_mean_deduped
  ))
}



main_output = main()

class(main_output)
main_output$annoNifHDB_updt
main_output$photic_samples_key
view(main_output$photic_samples_key)
class(main_output$photic_samples_key)


# Access the elements by their names
annoNifHDB_updt <- main_output$annoNifHDB_updt
cmap_coloc <- main_output$cmap_coloc
photic_samples_key <- main_output$photic_samples_key
size_fraction_key <- main_output$size_fraction_key
DNA_samples_key <- main_output$DNA_samples_key
sample_types_all <- main_output$sample_types_all
unique_sample_column_ids <- main_output$unique_sample_column_ids
unique_sample_id_key <- main_output$unique_sample_id_key
counts_df_T_lng <- main_output$counts_df_T_lng
counts_df_T_lng_AUID_deduped <- main_output$counts_df_T_lng_AUID_deduped
counts_df_T_AUID_deduped <- main_output$counts_df_T_AUID_deduped
counts_df_AUID_deduped <- main_output$counts_df_AUID_deduped
RA_df_T_lng <- main_output$RA_df_T_lng
RA_df_T_lng_mean_RA_AUID_deduped <- main_output$RA_df_T_lng_mean_RA_AUID_deduped
RA_df_T_mean_RA_AUID_deduped <- main_output$RA_df_T_mean_RA_AUID_deduped
RA_df_mean_RA_AUID_deduped <- main_output$RA_df_mean_RA_AUID_deduped





# main_output = list(
#     annoNifHDB_updt = annoNifHDB_updt,
#     cmap_coloc = cmap_coloc,
#     photic_samples_key = photic_samples_key,
#     size_fraction_key = result_list$size_fraction_key,
#     DNA_samples_key = result_list$DNA_samples_key,
#     sample_types_all = result_list$sample_types_all,
#     unique_sample_column_ids = result_list$unique_sample_column_ids,
#     unique_sample_id_key = unique_sample_id_key,
#     counts_df_T_lng = new_tibbles[[1]],
#     counts_df_T_lng_AUID_deduped = new_tibbles$df_T_lng_deduped,
#     counts_df_T_AUID_deduped = new_tibbles$df_T_deduped,
#     counts_df_AUID_deduped = new_tibbles$df_deduped,
#     RA_df_T_lng = new_tibbles[[5]],
#     RA_df_T_lng_mean_RA_AUID_deduped = new_tibbles$df_T_lng_mean_deduped,
#     RA_df_T_mean_RA_AUID_deduped = new_tibbles$df_T_mean_deduped,
#     RA_df_mean_RA_AUID_deduped = new_tibbles$df_mean_deduped
#   )

# annoNifHDB_updt = annoNifHDB_updt
#     cmap_coloc = cmap_coloc
#     photic_samples_key = photic_samples_key
#     size_fraction_key = result_list$size_fraction_key
#     DNA_samples_key = result_list$DNA_samples_key
#     sample_types_all = result_list$sample_types_all
#     unique_sample_column_ids = result_list$unique_sample_column_ids
#     unique_sample_id_key = unique_sample_id_key
#     counts_df_T_lng = new_tibbles[[1]]
#     counts_df_T_lng_AUID_deduped = new_tibbles$df_T_lng_deduped
#     counts_df_T_AUID_deduped = new_tibbles$df_T_deduped
#     counts_df_AUID_deduped = new_tibbles$df_deduped
#     RA_df_T_lng = new_tibbles[[5]]
#     RA_df_T_lng_mean_RA_AUID_deduped = new_tibbles$df_T_lng_mean_deduped
#     RA_df_T_mean_RA_AUID_deduped = new_tibbles$df_T_mean_deduped
#     RA_df_mean_RA_AUID_deduped = new_tibbles$df_mean_deduped

### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
cat("

Done loading script!!!

")
cat("

Woooooooohooooooo

!!!")

### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data




# # Replace depth zero values to a surface depth
# cmap_coloc <- cmap_coloc %>%
#   mutate(
#     depth =
#       ifelse(depth == 0,
#         2,
#         depth
#       )
#   )

# #-# to add DNA/RNA column
# cmap_coloc <- cmap_coloc %>%
#   mutate(nucleicAcidType = if_else(condition = LibrarySource == "METATRANSCRIPTOMIC", true = "RNA", false = "DNA"))


# #-# Add coastal/open ocean identification column
# coastal_open_ids <- read_csv("data/workspace/coastal_openocean_ids/sample_classifications_in_nifH_ASV_DB.csv") %>%
#   rename(
#     coastal_class = classification,
#     studyID = StudyID
#   )


# #-# Define photic and aphotic samples
# cmap_coloc <- cmap_coloc %>%
#   left_join(coastal_open_ids) %>%
#   mutate(
#     photic = case_when(
#       (coastal_class %in% "open ocean" & depth <= 100) ~ TRUE,
#       (coastal_class %in% "coastal" & depth <= 50) ~ TRUE,
#       TRUE ~ FALSE
#     )
#   )

# ## Make photic sample key
# photic_samples_key <- cmap_coloc %>%
#   filter(photic == TRUE) %>%
#   pull(SAMPLEID)


# cmap_coloc <- cmap_coloc %>%
#   mutate(
#     # Determine hemisphere
#     hemi = cut(.$lat,
#       breaks = c(-Inf, 0, Inf),
#       labels = c("southernHemi", "northernHemi")
#     ),
#     hemi = factor(hemi, c("northernHemi", "southernHemi")),
#     lat_abs = abs(lat),
#     # Standardize Size_fraction
#     Size_fraction = if_else(
#       Size_fraction %in% c("whole", "Sterivex", "0.22") | is.na(Size_fraction),
#       "0.2",
#       Size_fraction
#     ),
#     # Extract month and determine season
#     date = str_remove(time, " .+$"),
#     month = str_remove_all(str_extract(date, "-.+-"), "-"),
#     year = str_extract(date, "^\\d{4}"),
#     season = case_when(
#       (month %in% c("03", "04", "05") & hemi == "northernHemi") ~ "spring",
#       (month %in% c("06", "07", "08") & hemi == "northernHemi") ~ "summer",
#       (month %in% c("09", "10", "11") & hemi == "northernHemi") ~ "fall",
#       (month %in% c("12", "01", "02") & hemi == "northernHemi") ~ "winter",
#       (month %in% c("03", "04", "05") & hemi == "southernHemi") ~ "fall",
#       (month %in% c("06", "07", "08") & hemi == "southernHemi") ~ "winter",
#       (month %in% c("09", "10", "11") & hemi == "southernHemi") ~ "spring",
#       (month %in% c("12", "01", "02") & hemi == "southernHemi") ~ "summer",
#       TRUE ~ NA
#     ),
#     season = factor(season, c("winter", "spring", "summer", "fall")),
#     # Renaming columns
#     CMAP_NP_darwin = CMAP_NO3_darwin_clim_tblDarwin_Nutrient_Climatology / CMAP_PO4_darwin_clim_tblDarwin_Nutrient_Climatology,
#     CMAP_NP_Pisces_NRT = CMAP_NO3_tblPisces_NRT / CMAP_PO4_tblPisces_NRT,
#     CMAP_NP_WOA_clim = CMAP_nitrate_WOA_clim_tblWOA_Climatology / CMAP_phosphate_WOA_clim_tblWOA_Climatology,
#     CMAP_NP_WOA_clim = CMAP_n_an_clim_tblWOA_2018_1deg_Climatology / CMAP_p_an_clim_tblWOA_2018_1deg_Climatology
#   ) %>%
#   # Renaming columns with patterns
#   rename_with(~ str_remove(string = ., pattern = "CMAP_")) %>%
#   rename_with(~ coalesce(str_extract(string = ., pattern = "^[A-Za-z0-9]+_darwin"), .)) %>%
#   rename_with(~ str_replace(string = ., pattern = "clim_tblWOA_2018_qrtdeg_Climatology", replacement = "WOA_2018_qrtdeg")) %>%
#   rename_with(~ str_replace(string = ., pattern = "clim_tblWOA_2018_1deg_Climatology", replacement = "WOA_2018_1deg")) %>%
#   rename_with(~ str_replace(string = ., pattern = "WOA_clim_tblWOA_Climatology", replacement = "WOA_clim")) %>%
#   rename_with(~ str_replace(string = ., pattern = "clim_tblWOA_2018_MLD_qrtdeg_Climatology", replacement = "WOA_2018_MLD_qrtdeg")) %>%
#   rename_with(~ str_replace(string = ., pattern = "^C_", replacement = "conductivity_")) %>%
#   rename_with(~ str_replace(string = ., pattern = "^i_", replacement = "sigma_")) %>%
#   rename_with(~ str_replace(string = ., pattern = "^t_", replacement = "temp_")) %>%
#   rename_with(~ str_replace(string = ., pattern = "^A_", replacement = "AOU_")) %>%
#   rename_with(~ str_replace(string = ., pattern = "^O_", replacement = "O2_sat_")) %>%
#   rename_with(~ str_replace(string = ., pattern = "^n_", replacement = "NO3_")) %>%
#   rename_with(~ str_replace(string = ., pattern = "^p_", replacement = "phosphate_")) %>%
#   rename_with(~ str_replace(string = ., pattern = "^si_", replacement = "silica_")) %>%
#   rename_with(~ str_replace(string = ., pattern = "^M_", replacement = "MLD_")) %>%
#   rename_with(~ str_replace(string = ., pattern = "^s_", replacement = "salinity_")) %>%
#   # Determine geographical region
#   mutate(
#     geoRegion = cut(.$lat_abs,
#       breaks = c(-1, 23, 35, 66, 100),
#       labels = c("Eq/Trop", "Sub-Trop", "Temp", "Poles")
#     ),
#     geoRegion = factor(geoRegion, c("Eq/Trop", "Sub-Trop", "Temp", "Poles")),
#     logFe = log(Fe_tblPisces_NRT)
#   ) %>%
#   # ! FIXME: Fix Gradient crusies
#   mutate(
#     studyID = if_else(studyID == "Gradoville_2020" & grepl("2017", date), "Gradoville_2020_G2",
#       if_else(studyID == "Gradoville_2020" & grepl("2016", date), "Gradoville_2020_G1", studyID)
#     )
#   )

# ### add Ocean regions for each study ID
# studyid_regions <- read_csv("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/studyid_regions.csv")

# ### identy studies from Southern Ocean
# cmap_coloc <- cmap_coloc %>%
#   left_join(studyid_regions) %>%
#   # select(photic, SAMPLEID, studyID, study_ocean, everything()) %>%
#   mutate(
#     ocean = if_else(lat_abs >= 60 & hemi == "southernHemi", "Southern", study_ocean)
#   )





# #* Run functions to clean up CMAP
# cmap_coloc <- cmap_coloc %>%
#   rename_with(rename_cmap_columns) %>%
#   mutate(
#     hemi = determine_hemisphere(lat),
#     lat_abs = abs(lat),
#     Size_fraction = standardize_size_fraction(Size_fraction),
#     date = str_remove(time, " .+$"),
#     month = str_remove_all(str_extract(date, "-.+-"), "-"),
#     year = str_extract(date, "^\\d{4}"),
#     season = determine_season(month, hemi),
#     season = factor(season, c("winter", "spring", "summer", "fall")),
#     depth = replace_zero_depth(depth), # Replace depth zero values with a surface depth
#     logFe = log(Fe_tblPisces_NRT),
#     nucleicAcidType = add_nucleic_acid_type(LibrarySource), # Add DNA/RNA column
#     geoRegion = determine_geo_region(lat_abs)
#   ) %>%
#   fix_studyid() %>%
#   left_join(coastal_open_ids) %>%
#   mutate(
#     photic = case_when(
#       (coastal_class %in% "open ocean" & depth <= 100) ~ TRUE,
#       (coastal_class %in% "coastal" & depth <= 50) ~ TRUE,
#       TRUE ~ FALSE
#     )
#   )
