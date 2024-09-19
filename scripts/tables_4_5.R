#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
})


### - Table 3
### - Table 3



Table_SreadsAtEachStage_samples <- read_csv("analysis/Jmags/tables/Csvs_for_Tables_4_and_5/Table_SreadsAtEachStage_samples.csv") %>%
  rename(
    studyID = Study,
    SAMPLEID = Sample,
    InNonChimericASVs_final = InNonChimericASVs
  )

### - Table 4
# read_csv("~/mmorando@ucsc.edu_gDrive/My Drive/data/amplicon_review/all_studies/Jmag/Tables_and_plots/Workflow_stages/July25_has_updated_params_for_FilterAuids/workflowTable.csv")
workflowTable <- read_csv("analysis/Jmags/tables/Csvs_for_Tables_4_and_5/table_5_for_manuscript.csv") %>% rename(
  studyID = Study,
#   SAMPLEID = Sample,
  ReadsFilterAuids.Length_final = ReadsFilterAuids.Length
)  %>%
  select(-PctRetained)

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
