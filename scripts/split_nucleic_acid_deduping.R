### ! right now, some of the data is wrong
### ! I believe that I averaged over both duplicates and size fractions
## ! should be ok for now but must
# \FIXME:
## ! also, using lat and lon as is is likely not working do to precisoin
## ! there are lots of decimal points and this is likely making some duplicate
## ! samples into single samples bc the lat or lon is is .0001 off but the time is the same


### _ Loading in functions
cat("Load in the functions")

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

### _ Loading in keys
cat("Load in the keys 🔑🗝️ ")

## make DNA samples_key
DNA_samples_key <- CMAP_coloc %>%
  filter(nucleicAcidType == "DNA") %>%
  pull(SAMPLEID)

# make key to remove all unknown AUIDs and those without a nifh clusters identfied through annotations
(auids_w_annotations <- annoNifHDB_updt %>%
filter(!is.na(AUID)) %>% 
  pull(AUID)) 

### make a function to remove these from a dataframe
remove_unknown_auids <- function(
    df,
    annotation_key = auids_w_annotations) {
  df_filt <- df %>%
    filter(AUID %in% {{ annotation_key }})

  return(df_filt)
}



### - We have to have a df that identifies the different sample types and nucleic acid types
#### - To do this we need to split the df into a DNA and RNA
#### - Then identify the samples types
#### - Average by sample type and nucleic acid type
#### - Dedup

# Group by studyID, lat, lon, time, and depth
sample_types_DNA <- remove_samples_nucleic_acid(
  CMAP_coloc,
  "DNA",
  DNA_samples_key
) %>%
  group_by(studyID, lat, lon, time, depth) %>%
  mutate(sample_type = case_when(
    n() > 1 & n_distinct(Size_fraction) > 1 ~ "Two_Size_Fractions",
    n() > 1 & n_distinct(Size_fraction) == 1 ~ "Duplicate_Samples",
    n() == 1 ~ "Single_Sample"
  )) %>%
  ungroup() %>%
  select(SAMPLEID, sample_type, nucleicAcidType)

# sample_types_DNA_counts <- count_and_arrange(sample_types_DNA, "sample_type")

# add_total_row(sample_types_DNA_counts)
#

sample_types_RNA <- remove_samples_nucleic_acid(
  CMAP_coloc,
  "RNA",
  DNA_samples_key
) %>%
  group_by(studyID, lat, lon, time, depth) %>%
  mutate(sample_type = case_when(
    n() > 1 & n_distinct(Size_fraction) > 1 ~ "Two_Size_Fractions",
    n() > 1 & n_distinct(Size_fraction) == 1 ~ "Duplicate_Samples",
    n() == 1 ~ "Single_Sample"
  )) %>%
  ungroup() %>%
  select(SAMPLEID, sample_type, nucleicAcidType)

# sample_types_RNA_counts <- count_and_arrange(sample_types_RNA, "sample_type")

# add_total_row(sample_types_RNA_counts)


# sample_types_all_counts <- sample_types_RNA_counts %>%
#     bind_rows(sample_types_DNA_counts) %>%
#     add_total_row()

sample_types_all <- sample_types_DNA %>%
  bind_rows(sample_types_RNA)


### make unique sample id key

## i think this desiginates a unique sample
unique_sample_column_ids <- c(
  "studyID",
  "lat",
  "lon",
  "time",
  "depth",
  "sample_type",
  "nucleicAcidType"
)

unique_sample_id_key <- sample_types_all %>%
  left_join(
    CMAP_coloc
    # %>%
    # select(SAMPLEID, studyID, lat, lon, time, depth)
  ) %>%
  group_by(across(all_of(unique_sample_column_ids))) %>%
  # mutate(group_id = ifelse(n() > 1, as.character(cur_group_id()), NA)) %>%
  mutate(group_id = as.character(cur_group_id())) %>%
  ungroup() %>%
  select(SAMPLEID, sample_type, nucleicAcidType, group_id)

### we need a df to average over
## long form df work best so transform it
RA_df_T_lng <- transform_data_lng(
  input_df = nifhDB_RA_T,
  starts_with_col = "AUID",
  names_to_col = "AUID",
  values_to_col = "RA"
)


RA_df_T_lng_mean_RA_AUID_deduped <- main_average_ra_dedup_by_group(
  df_lng = RA_df_T_lng,
  AUID, group_id,
  mean_by = RA
)

# test <- RA_df_T_lng %>%
#   add_group_id() %>%
#   group_by(AUID, group_id) %>%
#   mutate(
#     # Calculate averages for the relevant columns
#     average_value_AUID = mean(RA, na.rm = T)
#     # across(starts_with("AUID"), ~ mean(., na.rm = T))
#   ) %>%
#   ungroup() %>% # print(n=100)
#   dedup_by_group(AUID, group_id)


# all.equal(RA_df_T_lng_mean_RA_AUID_deduped, test)

### for the counts data
### we need a df to average over
## long form df work best so transform it
counts_df_T_lng <- transform_data_lng(
  input_df = nifhDB_cnts_T,
  starts_with_col = "AUID",
  names_to_col = "AUID",
  values_to_col = "counts"
)

#  _ TODO: Need to add main function for just counts since you don't average by mean
### _ I need to used Jonathan's code that averages but keeps the total counts the same...
# RA_df_T_lng_mean_RA_AUID_deduped <- main_average_ra_dedup_by_group(
#   df_lng =
#     RA_df_T_lng, mean_by = RA,
#   AUID, group_id
# )


counts_df_T_lng_AUID_deduped <- counts_df_T_lng %>%
  add_group_id() %>%
  group_by(AUID, group_id) %>%
  # mutate(
  #   # Calculate averages for the relevant columns
  #   average_value_AUID = mean(counts, na.rm = T)
  #   # across(starts_with("AUID"), ~ mean(., na.rm = T))
  # ) %>%
  ungroup() %>%
  dedup_by_group(AUID, group_id)


# write_csv(RA_df_T_lng_mean_RA_AUID_deduped, "RA_df_T_lng_mean_RA_AUID_deduped.csv")



### _ now pivot wide ### _ now pivot wide
RA_df_T_mean_RA_AUID_deduped <- RA_df_T_lng_mean_RA_AUID_deduped %>%
  select(-all_of(c("RA", "group_id"))) %>%
  pivot_wider(
    names_from = AUID,
    values_from = average_value_AUID
  )

RA_df_mean_RA_AUID_deduped <- RA_df_T_lng_mean_RA_AUID_deduped %>%
  select(-all_of(c("RA", "group_id"))) %>%
  pivot_wider(
    names_from = SAMPLEID,
    values_from = average_value_AUID
  )


counts_df_T_AUID_deduped <- counts_df_T_lng_AUID_deduped %>%
  select(-all_of(c("group_id"))) %>%
  pivot_wider(
    names_from = AUID,
    values_from = counts
  )

counts_df_AUID_deduped <- counts_df_T_lng_AUID_deduped %>%
  select(-all_of(c("group_id"))) %>%
  pivot_wider(
    names_from = SAMPLEID,
    values_from = counts
  )


# write_csv(RA_df_T_mean_RA_AUID_deduped, "RA_df_T_mean_RA_AUID_deduped.csv")



### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
cat("

Done loading script!!!

")
cat("

Woooooooohooooooo

!!!")
cat("

Woooooooohooooooo

!!!")
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
