# # Transform data to long format for averaging and deduplication
# counts_df_T_lng <- transform_data_lng(
#     input_df = nifhDB_cnts,
#     starts_with_col = "AUID",
#     names_to_col = "AUID",
#     values_to_col = "counts"
# )

# Total reads in nifh database
tot_cnts <- counts_df_T_lng %>%
    # remove_aphotic_samples()  %>%
    select(counts) %>%
    sum()
cat("Total reads in nifh database:", tot_cnts, "\n")

# Total photic reads in nifh database
tot_cnts_photic <- counts_df_T_lng %>%
    remove_aphotic_samples() %>%
    select(counts) %>%
    sum()
cat("Total reads in nifh database, just photic samples:", tot_cnts_photic, "\n")

# write_csv(RA_df_T_lng_mean_RA_AUID_deduped, "RA_df_T_lng_mean_RA_AUID_deduped.csv")

# Total deduped
tot_cnts_deduped <- counts_df_T_lng_AUID_deduped %>%
    # remove_aphotic_samples()  %>%
    select(counts) %>%
    sum()
cat("Total reads in nifh database, deplicated:", tot_cnts_deduped, "\n")

# Total deduped photic
tot_cnts_deduped_photic <- counts_df_T_lng_AUID_deduped %>%
    remove_aphotic_samples() %>%
    select(counts) %>%
    sum()
cat("Total reads in nifh database, deplicated, photic:", tot_cnts_deduped_photic, "\n")




## -# photic samples
number_of_aphotic_samples <- nrow(cmap_coloc) - nrow(remove_aphotic_samples
(cmap_coloc))

number_of_photic_samples <- nrow(cmap_coloc) -
    number_of_aphotic_samples

cat("Total number of samples: ", nrow(cmap_coloc), "
Number of aphotic samples:", number_of_aphotic_samples, "
Number of photic samples:", number_of_photic_samples, "\n\n",
    sep = " "
)


## -# DNA/RNA samples
count_and_arrange(unique_sample_id_key, c("nucleicAcidType"))

number_of_DNA_samples <- nrow(cmapTab) - nrow(remove_samples_nucleic_acid(cmap_coloc, "RNA", DNA_samples_key))

number_of_RNA_samples <- nrow(cmapTab) - nrow(remove_samples_nucleic_acid(cmap_coloc, "DNA", DNA_samples_key))

cat("Total number of samples: ", nrow(cmap_coloc), "
Total number of DNA samples:", number_of_DNA_samples, "
Total number of RNA samples:", number_of_RNA_samples, "\n\n",
    sep = " "
)


## -# replicate samples
count_and_arrange(unique_sample_id_key, c("replicate_flag"))

number_of_replicate_samples <- unique_sample_id_key %>%
    filter(replicate_flag %in% "Replicate_Sample") %>%
    count() %>%
    pull()


cat("Number of replicate samples:", number_of_replicate_samples, "\n\n",
    sep = " "
)


#-# Size_fraction
nsf <- count_and_arrange(unique_sample_id_key, c("size_frac_flag"))

one_sf <- paste(nsf %>% filter(size_frac_flag == "One_Size_Fractions"), collapse = ": ", sep = ":")
two_sf <- paste(nsf %>% filter(size_frac_flag == "Two_Size_Fractions"), collapse = ": ", sep = ":")

cat("Number of", one_sf, "\n")
cat("Number of", two_sf, "\n")

# Stats on samples types
count_and_arrange(unique_sample_id_key, c("replicate_flag", "nucleicAcidType", "size_frac_flag"), replicate_flag)


## - nucleicAcidType

query_df <- dedup_by_group(cmap_coloc, group_id)

samples_per_nucleicAcidType <- count_and_arrange(query_df, "nucleicAcidType") %>%
    mutate(
        total = sum(n),
        percentage = n / total * 100
    ) %>%
    select(nucleicAcidType, n, percentage)

print(samples_per_nucleicAcidType, n = 1000)

## * by study id
samples_per_nucleicAcidType_studyid <- dedup_by_group(cmap_coloc, group_id) %>%
    # distinct(studyID, nucleicAcidType) %>%
    count_and_arrange(c("studyID", "nucleicAcidType")) %>%
    add_total_row(
        column_name = "studyID",
        summary_column = n,
        # pull_name = n,
        all_columns = FALSE
    )

print(samples_per_nucleicAcidType_studyid, n = 50)

### -  photic
samples_per_photic <- count_and_arrange(query_df, "photic") %>%
    mutate(
        total = sum(n),
        percentage = n / total * 100
    ) %>%
    select(photic, n, percentage)

print(samples_per_photic, n = 1000)

count_and_arrange(cmap_coloc, c("photic", "nucleicAcidType"))


# count_and_arrange(sample_types_all, c("sample_type", "SAMPLEID"), sample_type)  %>% add_total_row

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

# (sample_type <- count_and_arrange(sample_types_all, c("sample_type", "nucleicAcidType"), sample_type) %>%
#   # add_total_row(n, "sample_type", all_columns = FALSE) %>%
#   mutate(
#     new_group = paste(sample_type, nucleicAcidType, sep = "_"),
#     tibble_id = "sample_id"
#   ) %>%
#   add_percentage(n,
#     percentage,
#     grouping_by = NULL,
#     remove_columns = "total"
#   ))
