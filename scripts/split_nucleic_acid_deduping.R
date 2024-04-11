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

# Count the number of distinct size fractions for each sampling point
size_fraction_key <- CMAP_coloc %>%
  group_by(studyID, lat, lon, time, depth) %>%
  summarise(num_distinct_size_fractions = n_distinct(Size_fraction)) %>%
  mutate(
    size_frac_flag = case_when(
      (num_distinct_size_fractions == 1) ~ "One_Size_Fractions",
      (num_distinct_size_fractions > 1) ~ "Two_Size_Fractions",
    )
  ) %>%
  ungroup()




CMAP_coloc <- add_size_frac(CMAP_coloc, size_fraction_key, num_distinct_size_fractions)

### - We have to have a df that identifies the different sample types and nucleic acid types
#### - To do this we need to split the df into a DNA and RNA
#### - Then identify the samples types
#### - Average by sample type and nucleic acid type
#### - Dedup

# Define a sampling point for each by studyID, lat, lon, time, and depth
# This represents an individual sampling point
# These points can have multiple size fractions, nucleic acid types, or replicates
# Sampling can happen here at multiple times and depth but are considered different sampling points
sampling_point <- c("studyID", "lat", "lon", "time", "depth")

# Make this a column in CMAP df
CMAP_coloc <- CMAP_coloc %>%
  group_by(across(all_of(sampling_point))) %>%
  mutate(
    sample_point = cur_group_id(),
    num_dist_samp_pnts = n()
  )



sample_types_DNA <- remove_samples_nucleic_acid(
  CMAP_coloc,
  "DNA",
  DNA_samples_key
) %>%
  add_rep_flag() %>%
  select(sample_point, num_dist_samp_pnts, nucleicAcidType, size_frac_flag, Size_fraction, replicate_flag, everything()) %>% # view()
  ungroup() # %>%
# count(replicate_flag)


sample_types_RNA <- remove_samples_nucleic_acid(
  CMAP_coloc,
  "RNA",
  DNA_samples_key
) %>%
  add_rep_flag() %>%
  select(sample_point, num_dist_samp_pnts, nucleicAcidType, size_frac_flag, Size_fraction, replicate_flag, everything()) %>% # view()
  ungroup() # %>%
# count(replicate_flag)

sample_types_all <- sample_types_DNA %>%
  bind_rows(sample_types_RNA)

write_csv(x = sample_types_all, file = "sample_types_all.csv")

### make unique sample id key

## i think this desiginates a unique sample
unique_sample_column_ids <- c(
  "studyID",
  "lat",
  "lon",
  "time",
  "depth",
  # "sample_type",
  "replicate_flag",
  "size_frac_flag",
  "nucleicAcidType"
)

CMAP_coloc <- sample_types_all %>%
  group_by(across(all_of(unique_sample_column_ids))) %>%
  # mutate(group_id = ifelse(n() > 1, as.character(cur_group_id()), NA)) %>%
  mutate(group_id = as.character(cur_group_id())) %>%
  ungroup() # %>%
# select(SAMPLEID, sample_type, nucleicAcidType, group_id)
# distinct(group_id)

count_and_arrange(CMAP_coloc, c("replicate_flag", "nucleicAcidType", "size_frac_flag"), replicate_flag)

# Create unique id key for all sample ids
unique_sample_id_key <- CMAP_coloc %>%
  select(all_of(unique_sample_column_ids), SAMPLEID, group_id) %>%
  distinct(SAMPLEID, group_id, .keep_all = TRUE)









# sample_types_DNA <- remove_samples_nucleic_acid(
#   CMAP_coloc,
#   "DNA",
#   DNA_samples_key
# ) %>%
#   group_by(studyID, lat, lon, time, depth) %>%
#   mutate(sample_type = case_when(
#     n() > 1 & n_distinct(Size_fraction) > 1 ~ "Two_Size_Fractions",
#     n() > 1 & n_distinct(Size_fraction) == 1 ~ "Duplicate_Samples",
#     n() == 1 ~ "Single_Sample"
#   )) %>%
#   ungroup() %>%
#   select(SAMPLEID, sample_type, nucleicAcidType)

# # sample_types_DNA_counts <- count_and_arrange(sample_types_DNA, "sample_type")

# # add_total_row(sample_types_DNA_counts)
# #

# sample_types_RNA <- remove_samples_nucleic_acid(
#   CMAP_coloc,
#   "RNA",
#   DNA_samples_key
# ) %>%
#   group_by(studyID, lat, lon, time, depth) %>%
#   mutate(sample_type = case_when(
#     n() > 1 & n_distinct(Size_fraction) > 1 ~ "Two_Size_Fractions",
#     n() > 1 & n_distinct(Size_fraction) == 1 ~ "Duplicate_Samples",
#     n() == 1 ~ "Single_Sample"
#   )) %>%
#   ungroup() %>%
#   select(SAMPLEID, sample_type, nucleicAcidType)

# # sample_types_RNA_counts <- count_and_arrange(sample_types_RNA, "sample_type")

# # add_total_row(sample_types_RNA_counts)


# # sample_types_all_counts <- sample_types_RNA_counts %>%
# #     bind_rows(sample_types_DNA_counts) %>%
# #     add_total_row()

# sample_types_all <- sample_types_DNA %>%
#   bind_rows(sample_types_RNA)

# write_csv(x = sample_types_all, file = "sample_types_all.csv")


# ### make unique sample id key

# ## i think this desiginates a unique sample
# unique_sample_column_ids <- c(
#   "studyID",
#   "lat",
#   "lon",
#   "time",
#   "depth",
#   # "sample_type",
#   "replicate_flag",
#   "size_frac_flag",
#   "nucleicAcidType"
# )

# unique_sample_id_key <- sample_types_all %>%
#   left_join(
#     CMAP_coloc
#     # %>%
#     # select(SAMPLEID, studyID, lat, lon, time, depth)
#   ) %>%
#   group_by(across(all_of(unique_sample_column_ids))) %>%
#   # mutate(group_id = ifelse(n() > 1, as.character(cur_group_id()), NA)) %>%
#   mutate(group_id = as.character(cur_group_id())) %>%
#   ungroup() %>%
#   select(SAMPLEID, sample_type, nucleicAcidType, group_id)



# temp_bind  %>%
#   group_by(across(all_of(unique_sample_column_ids))) %>%
#   # mutate(group_id = ifelse(n() > 1, as.character(cur_group_id()), NA)) %>%
#   mutate(group_id = as.character(cur_group_id())) %>%
#   ungroup() %>%
#   # select(SAMPLEID, sample_type, nucleicAcidType, group_id)  %>%
# distinct(group_id)




# CMAP_coloc %>%
#   group_by(studyID, lat, lon, time, depth) %>%
#   mutate(
#     sample_point = cur_group_id(),
#     num_dist_samp_pnts = n()
#   ) %>%
#   mutate() %>%
#   mutate(
#     replicate_flag = case_when(
#       (num_dist_samp_pnts == 1) ~ "Single_Sample",
#       (num_dist_samp_pnts > 2) ~ "Replicate_Sample",
#       (n_distinct(nucleicAcidType) > 1 & n_distinct(size_frac_flag) == 1) ~ "Replicate_Sample",
#       (n_distinct(nucleicAcidType) == 1 & n_distinct(size_frac_flag) == 1) ~ "Replicate_Sample",
#     )
#   ) %>%
#   select(sample_point, num_dist_samp_pnts, nucleicAcidType, size_frac_flag, Size_fraction, replicate_flag, everything()) %>%
#   view()
# ungroup() %>%
#   count(replicate_flag)

# view()



# count_and_arrange(sample_types_all, c("sample_type", "nucleicAcidType"), sample_type)


# CMAP_coloc %>%
#   group_by(studyID, lat, lon, time, depth) %>%
#   mutate(n_distinct(nucleicAcidType)) %>%
#   view()








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

# now for the counts
# ! Do not apply means here because they are counts and mean counts would not work without normalization, the total number of reads would need to be maintained or else the information of proportions would be lost
# Instead, just one replicate is used
# Would not be too different from just collecting a single sample
counts_df_T_lng_AUID_deduped <- counts_df_T_lng %>%
  add_group_id() %>%
  # group_by(AUID, group_id) %>%
  # mutate(
  #   # Calculate averages for the relevant columns
  #   average_value_AUID = mean(counts, na.rm = T)
  #   # across(starts_with("AUID"), ~ mean(., na.rm = T))
  # ) %>%
  ungroup() %>%
  dedup_by_group(AUID, group_id) %>%
  ungroup()



# Total reads in nifh database
counts_df_T_lng %>%
  # remove_aphotic_samples()  %>%
  select(counts) %>%
  sum()

# Total photic reads in nifh database
counts_df_T_lng %>%
  remove_aphotic_samples() %>%
  select(counts) %>%
  sum()

# write_csv(RA_df_T_lng_mean_RA_AUID_deduped, "RA_df_T_lng_mean_RA_AUID_deduped.csv")

# Total deduped
counts_df_T_lng_AUID_deduped %>%
  # remove_aphotic_samples()  %>%
  select(counts) %>%
  sum()

# Total deduped photic
counts_df_T_lng_AUID_deduped %>%
  remove_aphotic_samples() %>%
  select(counts) %>%
  sum()


### _ now pivot wide ### _ now pivot wide
RA_df_T_mean_RA_AUID_deduped <- RA_df_T_lng_mean_RA_AUID_deduped %>%
  select(-any_of(c("RA", "group_id"))) %>%
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

### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
