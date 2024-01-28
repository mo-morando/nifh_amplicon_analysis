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
    Location,
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
  select(studyID, ocean, everything())

## * write out Table 1
write_csv(table_1, "~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/Table_1.csv")

### - Table 2

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

cmap_names <- tibble("CMAP_ID" = names(CMAP_20230719)) %>%
  filter(grepl("CMAP", CMAP_ID))

cmap_variables <- extract_cmap_info_from_names(cmap_names)


# cmap_variables <- tibble("CMAP_ID" = names(CMAP_20230719)) %>%
#   filter(grepl("CMAP", CMAP_ID)) %>%
#   extract(CMAP_ID, into = c("Variable", "Table"), regex = "CMAP_(.*?)_(tbl.*)", remove = FALSE)

tibble(str_match(names(CMAP_20230719), "CMAP_(.*?)_(tbl.*)"))

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
  filter(is.na(Long_Name)))

cmap_catalog %>%
  filter(
    Table_Name == "tblWOA_Climatology" #| Table_Name == "tblWOA_2018_1deg_Climatology"
    # filter(grepl("argo",Table_Name, ignore.case = TRUE) #| Table_Name == "tblWOA_2018_1deg_Climatology"
  ) %>%
  view()


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


supp_table_2 <- cmap_catalog %>%
  select(-(Variable_25th:Unstructured_Variable_Metadata))



write_csv(supp_table_2, "~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/SupplementaryTable_2.csv")




### - Table 3
### - Table 3

# read_csv("~/mmorando@ucsc.edu_gDrive/My Drive/data/amplicon_review/all_studies/Jmag/Tables_and_plots/Pipeline_stages/Table_SreadsAtEachStage_samples.csv")

Table_SreadsAtEachStage_samples <- read_csv("analysis/Jmags/tables/Table_SreadsAtEachStage_samples.csv") %>%
  rename(
    studyID = Study,
    SAMPLEID = Sample,
    InNonChimericASVs_final = InNonChimericASVs
  )

### - Table 4
# read_csv("~/mmorando@ucsc.edu_gDrive/My Drive/data/amplicon_review/all_studies/Jmag/Tables_and_plots/Workflow_stages/July25_has_updated_params_for_FilterAuids/workflowTable.csv")
workflowTable <- read_csv("analysis/Jmags/tables/workflowTable.csv") %>% rename(
  studyID = Study,
  SAMPLEID = Sample,
  ReadsFilterAuids.Length_final = ReadsFilterAuids.Length
)

### add percentage column
Table_SreadsAtEachStage_samples <- Table_SreadsAtEachStage_samples %>%
  mutate(
    PctReadPairsRetained = round(((1 - (Initial - InNonChimericASVs_final) / Initial) * 100), 1)
    # test = round(((1 - (Initial - contains("final")) / Initial) * 100), 1)
  )

##
workflowTable <- workflowTable %>%
  mutate(
    PctReadPairsRetained = ifelse(
      ReadsFilterAuids.Length_final == 0,
      0,
      round(((1 - (ReadsPipeline -
        ReadsFilterAuids.Length_final) / ReadsPipeline) * 100), 1)
    )
  )

#-# function to summarize these files and produce Tables 3 and 4
summarise_workflow_stages_table <- function(
    df,
    grp_by,
    rnd_to,
    initial_col,
    final_col) {
  summarised_table <- df %>%
    group_by({{ grp_by }}) %>%
    reframe(across(where(is.numeric), ~ mean(.))) %>%
    ungroup() %>%
    bind_rows(
      tibble(
        {{ grp_by }} := "Total mean",
        df %>%
          summarise(across(where(is.numeric), ~ mean(.))) #* This gives the mean over all the reads, not by study ID
      )
    ) %>%
    bind_rows(
      tibble(
        {{ grp_by }} := "Total median",
        df %>%
          summarise(across(where(is.numeric), ~ median(.))) #* This gives the median over all the reads, not by study ID
      )
    ) %>%
    bind_rows(
      tibble(
        {{ grp_by }} := "Total sum",
        df %>%
          summarise(across(where(is.numeric), ~ sum(.))) #* this gives the total over all the reads
      )
    ) %>%
    mutate(
      PctReadPairsRetained =
        ifelse(test = studyID == "Total sum",
          yes = round((1 - ({{ initial_col }} - {{ final_col }}) / {{ initial_col }}), 1) * 100,
          no = PctReadPairsRetained
        )
    ) %>%
    mutate(across(where(is.numeric), ~ round(., rnd_to)))

  return(summarised_table)
}

format_table <- function(table) {
  formatted_table <- table %>%
    mutate(across(!contains("studyID") & !contains("PctReadPairsRetained"), ~ format(., scientific = TRUE, digits = 2)))

  return(formatted_table)
}

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
    Rare = ReadsFilterAuids.Rare,
    NonNifH = ReadsFilterAuids.NonNifH,
    Length = ReadsFilterAuids.Length_final,
    "Retained (%)" = PctReadPairsRetained
  ))
# mutate(across(where(is.numeric), ~sum(.)))

#-# write out files

write_csv(Table_SreadsAtEachStage_samples_summary_format_rename, "/Users/mo/Projects/nifH_amp_project/myWork/analysis/tables/Table_SreadsAtEachStage_samples_summary_format_rename.csv")

write_csv(Table_SreadsAtEachStage_samples_summary_format_rename, "~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/Table_3.csv")


write_csv(workflowTable_summary_format_rename, "/Users/mo/Projects/nifH_amp_project/myWork/analysis/tables/workflowTable_summary_format_rename.csv")

write_csv(workflowTable_summary_format_rename, "~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/Table_4.csv")
