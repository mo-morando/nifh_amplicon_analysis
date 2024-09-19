library(tidyverse)

### tables

### make Table 1
## ! requires manual work so don't really want to mess with this too much

### - putting together Table 1
all_studies <- read_csv("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/alllStudiesAssessedMain.csv")

all_together <- read_csv("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/AllTogether.csv")

studyid_regions <- read_csv("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/studyid_regions.csv")



# all_together_way <- read_csv("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/AllTogetherWay.csv")

# all_studiesStages <- read_csv("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/all_studiesAssessed_Stages.csv")

## * make subsets to containing the columns i want to keep
(all_together_sub <- all_together %>%
  mutate(
    samples = ifelse(studyID == "Rahav_2016", NA, samples)
  ) %>%
  select(
    studyID,
    # Location,
    SRAinformation,
    samples,
    bytes_gb,
    bases_g,
    references,
    DOI
  )
)

all_studies_sub <- all_studies %>%
  select(studyID, reason, Refs, comments, primers_used, used, reason_excluded)

## * merge them into Supplementary Table 1
(supp_table_1 <- all_together_sub %>%
  full_join(all_studies_sub) %>%
  filter(!studyID %in% c("Henke_2018", "Fran_2019")) %>%
  select(-c(bytes_gb, bases_g, reason, Refs, comments))
# %>%
# full_join(all_together_way_sub)
)
# view()

## * write out Supplementary_Table_1
write_csv(supp_table_1, "~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/Supplementary_Table_1.csv")

### * Filter columns to produced Table 1
## We need less information for a tighter table
# names(supp_table_1)

table_1_names_key <- c(
  "studyID",
  "samples",
  "SRAinformation",
  "references",
  "DOI"
)

table_1 <- supp_table_1 %>%
  filter(used == "y") %>%
  select(all_of(table_1_names_key)) %>%
  left_join(studyid_regions) %>%
  select(studyID, study_ocean, everything()) %>%
  rename(
    "Study ID" = studyID,
    "Ocean" = study_ocean,
    "SRA Information" = SRAinformation,
    "References" = references
  )

## * write out Table 1
write_csv(table_1, "~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/Table_1.csv")

### - Supplemental Table 2

##* functino to extract columns from our CMAP file to use to merge with the CMAP cataloge (which I will download from CMAP) that has all the information all the CMAP data we pulled so I can make a supp table

extract_cmap_info_from_names <- function(cmap_names) {
  extracted_cmap_info <- cmap_names %>%
    extract(
      CMAP_ID,
      into = c("Variable", "Table_Name"),
      regex = "CMAP_(.*?)_(tbl.*)",
      remove = FALSE
    )

  return(extracted_cmap_info)
}

cmap_names <- tibble("CMAP_ID" = names(cmap_coloc)) %>%
  filter(grepl("CMAP", CMAP_ID))

cmap_variables <- extract_cmap_info_from_names(cmap_names)


# cmap_variables <- tibble("CMAP_ID" = names(CMAP_20230719)) %>%
#   filter(grepl("CMAP", CMAP_ID)) %>%
#   extract(CMAP_ID, into = c("Variable", "Table"), regex = "CMAP_(.*?)_(tbl.*)", remove = FALSE)

tibble(str_match(names(cmap_coloc), "CMAP_(.*?)_(tbl.*)"))

### _ write out table that will be used to query the cmap catalog and extract only the varibles used in this study
## _ the script to do this next step can be found
# _ /Users/mo/Projects/nifH_amp_project/myWork/data/databases/CMAP/scripts/get_cmap_catalog.R
## _ must use conda environmnet --> cmap4r in order to use this R package
write_csv(cmap_variables, "~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/cmap_variables.csv")

##* Read in cmap catalog with nifh database specific
cmap_catalog <- read_csv("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/merged_metadata/nifHDB/firstAttempt/20230719_colocalized_nifH_CMAP_catalog.csv")


# !# there are temporarily NAs that I need to deal with
## ! use the bit for now
(cmap_var_cat_na <- cmap_variables %>%
  left_join(cmap_catalog) %>%
  filter(!is.na(Long_Name)))

cmap_catalog %>%
  filter(
    # Table_Name == "tblWind_NRT_hourly" #| Table_Name == "tblArgoMerge_REP"
    Table_Name == "tblMercator_MLD_NRT" #| Table_Name == "tblArgoMerge_REP"
    # filter(grepl("argo",Table_Name, ignore.case = TRUE) #| Table_Name == "tblWOA_2018_1deg_Climatology"
  ) %>%
  view()
  print(n = 30)


cmap_var_cat_replace <- tibble(
  CMAP_ID = c(
    "CMAP_sla_nrt_tblAltimetry_NRT_Signal",
    "CMAP_wind_stress_curl_tblWind_NRT_hourly",
    "CMAP_wind_stress_tblWind_NRT_hourly",
    "CMAP_i_an_clim_tblWOA_2018_qrtdeg_Climatology",
    "CMAP_i_mn_clim_tblWOA_2018_qrtdeg_Climatology",
    "CMAP_i_an_clim_tblWOA_2018_1deg_Climatology",
    "CMAP_i_mn_clim_tblWOA_2018_1deg_Climatology",
    "CMAP_mld_da_median_argo_climmld_dt_median_argo_clim_tblArgo_MLD_Climatology"
  ),
  CMAP_ID_new = c(
    "CMAP_sla_tblAltimetry_NRT_Signal",
    "CMAP_wind_curl_tblWind_NRT_hourly",
    "CMAP_stress_curl_tblWind_NRT_hourly",
    "CMAP_I_an_clim_tblWOA_2018_qrtdeg_Climatology",
    "CMAP_I_mn_clim_tblWOA_2018_qrtdeg_Climatology",
    "CMAP_I_an_clim_tblWOA_2018_1deg_Climatology",
    "CMAP_I_mn_clim_tblWOA_2018_1deg_Climatology",
    "CMAP_mld_da_median_argo_clim_tblArgo_MLD_Climatology"
  )
)

cmap_var_cat_na_replace <- cmap_var_cat_na %>%
  left_join(cmap_var_cat_replace, by = "CMAP_ID") %>%
  mutate(CMAP_ID = coalesce(CMAP_ID_new, CMAP_ID)) %>%
  select(-CMAP_ID_new)

temp_cat <- extract_cmap_info_from_names(
  cmap_var_cat_na_replace
  %>% select(CMAP_ID)
)

temp_cat %>%
  # filter(!Variable == "mld_da_median_argo_climmld_dt_median_argo_clim") %>%
  left_join(cmap_catalog, by = c("Variable", "Table_Name")) # %>%
# filter(is.na(Long_Name)) %>%
view()

## ! use the bit for now


supp_table_2 <- cmap_var_cat_na %>%
  select(-(Variable_25th:Unstructured_Variable_Metadata))



write_csv(supp_table_2, "~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/SupplementaryTable_2.csv")



