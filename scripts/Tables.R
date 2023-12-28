### tables

### make Table 1
## ! requires manual work so don't really want to mess with this too much

### - putting together Table 1
all_studies <- read_csv("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/alllStudiesAssessedMain.csv")

all_together <- read_csv("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/AllTogether.csv")

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
  select(all_of(table_1_names_key))

## * write out Table 1
write_csv(table_1, "~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/Table_1.csv")

### - Table 2

cmap_variables <- CMAP_20230719 %>%
  # select(-c(SAMPLEID:geo_loc_name)) %>%
  select(contains("CMAP")) %>%
  names()

table_2 <- tibble(cmap_variables, "brief_description" = NA)

write_csv(table_2, "~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/Table_2.csv")


### - Table 3
### - Table 3

Table_SreadsAtEachStage_samples <- read_csv("~/mmorando@ucsc.edu_gDrive/My Drive/data/amplicon_review/all_studies/Jmag/Tables_and_plots/Pipeline_stages/Table_SreadsAtEachStage_samples.csv") %>%
  rename(
    studyID = Study,
    SAMPLEID = Sample,
    InNonChimericASVs_final = InNonChimericASVs
  )

### - Table 4
workflowTable <- read_csv("~/mmorando@ucsc.edu_gDrive/My Drive/data/amplicon_review/all_studies/Jmag/Tables_and_plots/Workflow_stages/July25_has_updated_params_for_FilterAuids/workflowTable.csv") %>% rename(
  studyID = Study,
  SAMPLEID = Sample,
  ReadsFilterAuids.Length_final = ReadsFilterAuids.Length
)

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
  ) %>%
  view()


summarise_table <- function(
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
        {{ grp_by }} := "Total sum",
        df %>%
          summarise(across(where(is.numeric), ~ sum(.)))
      ) %>%
        bind_rows(
          tibble(
            {{ grp_by }} := "Study mean",
            df %>%
              summarise(across(where(is.numeric), ~ mean(.)))
          )
        ) %>%
        bind_rows(
          tibble(
            {{ grp_by }} := "Study median",
            df %>%
              summarise(across(where(is.numeric), ~ median(.)))
          )
        )
    ) %>%
    ##  FIXME:
    ## ! Right now, workflow file doesn't have a PctReadPairsRetained column yet
    ## ! So we can't do this part below and instead need to do it outside this function in the next steps below
    mutate(
      PctReadPairsRetained = ifelse(studyID == "Total_sum",
        ({{ initial_col }} - {{ final_col }}) / {{ initial_col }} * 100,
        PctReadPairsRetained
      )
    ) %>%
    mutate(across(where(is.numeric), ~ round(., rnd_to)))

  return(summarised_table)
}

Table_SreadsAtEachStage_samples_summary <- summarise_table(
  df = Table_SreadsAtEachStage_samples,
  grp_by = studyID,
  rnd_to = 1,
  initial_col = Initial,
  final_col = InNonChimericASVs_final
)

# ##  FIXME:
# ## ! Right now, workflow file doesn't have a PctReadPairsRetained column yet
# ## ! So we can't do this part below and instead need to do it outside this function in the next steps below
# Table_SreadsAtEachStage_samples_summary %>% mutate(
#   PctReadPairsRetained = ifelse(studyID == "Total_sum",
#     round((Initial - InNonChimericASVs) / Initial * 100, 1),
#     PctReadPairsRetained
#   )
# )

workflowTable_summary <- summarise_table(
  df = workflowTable,
  grp_by = studyID,
  rnd_to = 1,
  initial_col = ReadsPipeline,
  final_col = ReadsFilterAuids.Length_final
)

write_csv(Table_SreadsAtEachStage_samples_summary, "/Users/mo/Projects/nifH_amp_project/myWork/analysis/tables/Table_SreadsAtEachStage_samples_summary.csv")

write_csv(Table_SreadsAtEachStage_samples_summary, "~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/Table_3.csv")


write_csv(workflowTable_summary, "/Users/mo/Projects/nifH_amp_project/myWork/analysis/tables/workflow_summary.csv")

write_csv(workflowTable_summary, "~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/Table_4.csv")
