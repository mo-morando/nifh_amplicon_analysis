# nifhDB_cnts
# annoNifHDB_updt

library(patchwork)

### _ Loading in functions
cat("Load in the functions")

## - functions
cat("Functions needed:
sum_tax
remove_samples_nucleic_acid")

# cluster_stats <- sum_tax(
#   df1 = nifhDB_cnts,
#   df2 = annoNifHDB_updt,
#   sw_col = "AUID",
#   nm_col = "SAMPLEID",
#   vls_col = "cnts",
#   by_col = "AUID",
#   gp1 = nifH_cluster,
#   value = cnts
# )

# print(cluster_stats, n = 50)


#### - count data

### * make query data frame where you remove RNA samples
## * add studyID column
query_df <- counts_df_T_lng_AUID_deduped %>%
  replace_na(list(. = 0)) %>%
  remove_samples_nucleic_acid(nucleic_acid_type = "DNA") %>%
  left_join(CMAP_coloc %>%
    select(SAMPLEID, studyID))

## * calculate percentages by nifH cluster
(nifhdb_all_counts_AUID_dedup <- sum_tax(
  abundance_table = query_df,
  annotation_table = annoNifHDB_updt,
  joining_by_col = "AUID",
  grouping_by_1 = nifH_cluster,
  grouping_by_2 = NULL,
  sum_by = counts,
  total_counts_id = total_counts_nifH_cluster,
  grouping_for_percentage = NULL,
  #   sum_by_percent = total_counts,
  percentage_id = percentage_nifH_cluster,
  nifH_cluster
))

## * Clean up table so that everything below given threshold is summed together as 'other'
(nifhdb_all_counts_AUID_dedup_clean <- clean_percentages(
  df = nifhdb_all_counts_AUID_dedup,
  grouping_by = NULL,
  clean_var = percentage_nifH_cluster,
  threshold = 1.0,
  # sum_by = percentage_nifH_cluster_study_id,
  percentage_id = percentage_nifH_cluster,
  column_var_mutate = nifH_cluster_modified,
  column_var = nifH_cluster
))

#### - By study ID

## * calculate over total data, so each study ID is a percentage of the total
(nifhdb_all_counts_AUID_dedup_total_study_id <- sum_tax(
  abundance_table = query_df,
  annotation_table = annoNifHDB_updt,
  joining_by_col = "AUID",
  sum_by = counts,
  grouping_by_1 = nifH_cluster,
  grouping_by_2 = studyID,
  grouping_for_percentage = NULL, ### make this NULL to get percentages of total instead of percentages by studyID
  total_counts_id = total_counts_nifH_cluster_total_study_id,
  #   sum_by_percent = total_counts,
  percentage_id = percentage_total_counts_nifH_cluster_total_study_id
))


## * Clean up table so that everything below a certain threshold is summed together as 'other'
(nifhdb_all_counts_AUID_dedup_total_study_id_clean <- clean_percentages(
  df = nifhdb_all_counts_AUID_dedup_total_study_id,
  grouping_by = studyID,
  clean_var = percentage_total_counts_nifH_cluster_total_study_id,
  # threshold = 2.5,
  threshold = 1.0,
  # sum_by = percentage_nifH_cluster_study_id,
  percentage_id = percentage_total_counts_nifH_cluster_total_study_id_clean,
  column_var_mutate = nifH_cluster_modified,
  column_var = nifH_cluster
))
### ! FIXME: the percentages are being calculated by totals and not totals within each study ID


## * calculate over each study ID, so each percentage is based on its own study.
(nifhdb_all_counts_AUID_dedup_study_id_total <- sum_tax(
  abundance_table = query_df,
  annotation_table = annoNifHDB_updt,
  joining_by_col = "AUID",
  sum_by = counts,
  grouping_by_1 = nifH_cluster,
  grouping_by_2 = studyID,
  grouping_for_percentage = studyID, ### make this NULL to get percentages of total instead of percentages by studyID
  total_counts_id = total_counts_nifH_cluster_study_id_total,
  #   sum_by_percent = total_counts,
  percentage_id = percentage_total_counts_nifH_cluster_study_id_total,
  nifH_cluster, studyID
))


## * Clean up table so that everything below a certain threshold is summed together as 'other'
(nifhdb_all_counts_AUID_dedup_study_id_total_clean <- clean_percentages(
  df = nifhdb_all_counts_AUID_dedup_study_id_total,
  grouping_by = studyID,
  clean_var = percentage_total_counts_nifH_cluster_study_id_total,
  # threshold = 2.5,
  threshold = 1.0,
  # sum_by = percentage_nifH_cluster_study_id,
  percentage_id = percentage_total_counts_nifH_cluster_study_id_total_clean,
  column_var_mutate = nifH_cluster_modified,
  column_var = nifH_cluster
))

# nifhdb_all_counts_AUID_dedup_study_id_total %>%
#   select(studyID, nifH_cluster, percentage_total_counts_nifH_cluster_study_id_total) %>%
#   arrange(studyID) %>%
#   print(n = 300)


#### - relative abundance data

### * make query data frame where you remove RNA samples
## * add studyID column
query_df <- RA_df_T_lng_mean_RA_AUID_deduped %>%
  replace_na(list(. = 0)) %>%
  remove_samples_nucleic_acid(nucleic_acid_type = "DNA") %>%
  left_join(CMAP_coloc %>%
    select(SAMPLEID, studyID))

## * calculate percentages by nifH cluster
(nifhdb_all_RA_AUID_dedup <- sum_tax(
  abundance_table = query_df,
  annotation_table = annoNifHDB_updt,
  joining_by_col = "AUID",
  grouping_by_1 = nifH_cluster,
  grouping_by_2 = NULL,
  sum_by = average_value_AUID,
  total_counts_id = total_rel_abund_nifH_cluster,
  grouping_for_percentage = NULL,
  #   sum_by_percent = total_counts,
  percentage_id = percentage_nifH_cluster,
  nifH_cluster
))


## * Clean up table so that everything below given threshold is summed together as 'other'
(nifhdb_all_RA_AUID_dedup_clean <-
  clean_percentages(
    df = nifhdb_all_RA_AUID_dedup,
    grouping_by = NULL,
    clean_var = percentage_nifH_cluster,
    threshold = 1.0,
    # sum_by = percentage_nifH_cluster_study_id,
    percentage_id = percentage_nifH_cluster,
    column_var_mutate = nifH_cluster_modified,
    column_var = nifH_cluster
  ))

#### - By study ID

## * calculate over total data, so each study ID is a percentage of the total
(nifhdb_all_RA_AUID_dedup_total_study_id <- sum_tax(
  abundance_table = query_df,
  annotation_table = annoNifHDB_updt,
  joining_by_col = "AUID",
  grouping_by_1 = nifH_cluster,
  grouping_by_2 = studyID,
  sum_by = average_value_AUID,
  grouping_for_percentage = NULL, ### make this NULL to get percentages of total instead of percentages by studyID
  total_counts_id = total_rel_abund_total_nifH_cluster,
  #   sum_by_percent = total_counts,
  percentage_id = percentage_total_rel_abund_total_nifH_cluster_study_id,
  nifH_cluster, studyID
))

# nifhdb_all_RA_AUID_dedup_total_study_id %>%
#   select(nifH_cluster, percentage_total_rel_abund_total_nifH_cluster_study_id) %>%
#   filter(percentage_total_rel_abund_total_nifH_cluster_study_id > 0) %>%
#   view()


## * Clean up table so that everything below a certain threshold is summed together as 'other'
nifhdb_all_RA_AUID_dedup_total_study_id_clean <- nifhdb_all_RA_AUID_dedup_total_study_id %>%
  clean_percentages(
    grouping_by = studyID,
    clean_var = percentage_total_rel_abund_total_nifH_cluster_study_id,
    # threshold = 2.5,
    threshold = 1.0,
    # sum_by = percentage_nifH_cluster_study_id,
    percentage_id = percentage_total_rel_abund_total_nifH_cluster_study_id_clean,
    column_var_mutate = nifH_cluster_modified,
    column_var = nifH_cluster
  )


## * calculate over each study ID, so each percentage is based on its own study.
(nifhdb_all_rel_abund_AUID_dedup_study_id_total <- sum_tax(
  abundance_table = query_df,
  annotation_table = annoNifHDB_updt,
  joining_by_col = "AUID",
  sum_by = average_value_AUID,
  grouping_by_1 = nifH_cluster,
  grouping_by_2 = studyID,
  grouping_for_percentage = studyID, ### make this NULL to get percentages of total instead of percentages by studyID
  total_counts_id = total_rel_abund_nifH_cluster_study_id_total,
  #   sum_by_percent = total_counts,
  percentage_id = percentage_total_rel_abund_nifH_cluster_study_id_total,
  nifH_cluster, studyID
))


## * Clean up table so that everything below a certain threshold is summed together as 'other'
(nifhdb_all_rel_abund_AUID_dedup_study_id_total_clean <- clean_percentages(
  df = nifhdb_all_rel_abund_AUID_dedup_study_id_total,
  grouping_by = studyID,
  clean_var = percentage_total_rel_abund_nifH_cluster_study_id_total,
  # threshold = 2.5,
  threshold = 1.0,
  # sum_by = percentage_nifH_cluster_study_id,
  percentage_id = percentage_total_rel_abund_nifH_cluster_study_id_total_clean,
  column_var_mutate = nifH_cluster_modified,
  column_var = nifH_cluster
))

# nifhdb_all_rel_abund_AUID_dedup_study_id_total %>%
#   select(studyID, nifH_cluster, percentage_total_rel_abund_nifH_cluster_study_id_total) %>%
#   arrange(studyID) %>%
#   print(n = 300)



### ! \TODO: Starting to develop AUID stats
## ! needed for Ecology manuscript analysis


(query_df <- counts_df_T_lng_AUID_deduped %>%
  replace_na(list(. = 0)) %>%
  remove_samples_nucleic_acid(nucleic_acid_type = "DNA") %>%
  remove_aphotic_samples() %>%
  # filter(AUID %in% unknownnan_auids) %>%
  left_join(CMAP_coloc %>%
    select(SAMPLEID, studyID)) %>%
  arrange(desc(counts))
)

## * calculate percentages by AUID
(temp <- sum_tax(
  abundance_table = query_df,
  annotation_table = annoNifHDB_updt,
  joining_by_col = "AUID",
  grouping_by_1 = AUID,
  grouping_by_2 = NULL,
  sum_by = counts,
  total_counts_id = total_counts_AUID,
  grouping_for_percentage = NULL,
  #   sum_by_percent = total_counts,
  percentage_id = percentage_AUID,
  nifH_cluster
))

#### _ PLOTS#### _ PLOTS#### _ PLOTS#### _ PLOTS#### _ PLOTS
#### _ PLOTS#### _ PLOTS#### _ PLOTS#### _ PLOTS#### _ PLOTS

### - all the df just made
# List all variables in the environment
all_vars <- ls()

# Find variables containing "cluster" in their names and remove them
cat("

Tibbles created:", all_vars[grep("nifhdb", all_vars)],
  sep = "\n"
)


# [1] "nifhdb_all_counts_AUID_dedup"
# [2] "nifhdb_all_counts_AUID_dedup_clean"
# [3] "nifhdb_all_RA_AUID_dedup"
# [4] "nifhdb_all_RA_AUID_dedup_clean"
# [5] "nifhdb_all_RA_AUID_dedup_total_study_id"
# [6] "nifhdb_all_RA_AUID_dedup_total_study_id_clean"
# [7] "nifhdb_all_rel_abund_AUID_dedup_study_id_total"
# [8] "nifhdb_all_rel_abund_AUID_dedup_study_id_total_clean"

### - pie chart of nifh cluster percentages

#### by counts
#### just DNA, deduplicated by AUID and group_id
pie_chart_AUID_depud_DNA_counts <- ggplot(nifhdb_all_counts_AUID_dedup_clean, aes(
  x = 1,
  # x = studyID,
  y = percentage_nifH_cluster,
  fill = nifH_cluster_modified,
  # colour = studyID,
  # group = nifH_cluster
)) +
  geom_bar(
    stat = "identity",
    # position = position_stack(reverse = FALSE)
  ) +
  # geom_histogram(binwidth = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(
    # title = "nifH dabtbase - % of total counts for each clusters",
    # subtitle = "DNA only with replicates averaged",
    fill = "nifH Cluster"
  ) +

  # Add labels with lines for small slices
  geom_text(
    # aes(label = ifelse(percentage_nifH_cluster >= 5, paste(nifH_cluster, "\n", percentage_nifH_cluster, "%"), ""), x = 1.25),
    aes(label = ifelse(percentage_nifH_cluster >= 5, paste0(percentage_nifH_cluster, "%"), ""), x = 1.15),
    position = position_stack(vjust = 0.5),
    hjust = 0.5,
    size = 12,
    show.legend = FALSE
  ) +
  geom_text(
    # aes(label = ifelse(percentage_nifH_cluster >= 5, paste(nifH_cluster, "\n", percentage_nifH_cluster, "%"), ""), x = 1.25),
    aes(
      label = ifelse(percentage_nifH_cluster < 5, paste0(percentage_nifH_cluster, "%"), ""),
      x = 1.5
    ),
    position = position_stack(vjust = 0.5),
    hjust = 0.5,
    size = 5,
    show.legend = FALSE,
    # nudge_y = 0.05
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 17, face = "bold"),
  ) +
  guides(fill = guide_legend(
    nrow = 1,
    byrow = TRUE,
    title.position = "top",
    title.hjust = 0.5
  ))

# pie_chart_facet <- pie_chart +
#   facet_wrap("studyID")

# Display the pie chart
print(pie_chart_AUID_depud_DNA_counts)

# print(pie_chart_facet)

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/pie_chart_nifhDF_DNA_dedup_perc_tot_cnts_nifH_cluster.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)

annoNifHDB_updt %>%
  filter(nifH_cluster == "1O/1P") %>%
  distinct(CON)

#### by relative abundance
#### just DNA, deduplicated by AUID and group_id
pie_chart_AUID_depud_DNA_rel_abund <- ggplot(nifhdb_all_RA_AUID_dedup_clean, aes(
  x = 1,
  # x = studyID,
  y = percentage_nifH_cluster,
  fill = nifH_cluster_modified,
  # colour = studyID,
  # group = nifH_cluster
)) +
  geom_bar(
    stat = "identity",
    # position = position_stack(reverse = FALSE)
  ) +
  # geom_histogram(binwidth = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(
    # title = "nifH dabtbase - % of total relative abundance for each clusters",
    # subtitle = "DNA only with replicates averaged",
    fill = "nifH Cluster"
  ) +

  # Add labels with lines for small slices
  geom_text(
    # aes(label = ifelse(percentage_nifH_cluster >= 5, paste(nifH_cluster, "\n", percentage_nifH_cluster, "%"), ""), x = 1.25),
    aes(label = ifelse(percentage_nifH_cluster >= 5, paste0(percentage_nifH_cluster, "%"), ""), x = 1.15),
    position = position_stack(vjust = 0.5),
    hjust = 0.5,
    size = 12,
    show.legend = FALSE
  ) +
  geom_text(
    # aes(label = ifelse(percentage_nifH_cluster >= 5, paste(nifH_cluster, "\n", percentage_nifH_cluster, "%"), ""), x = 1.25),
    aes(
      label = ifelse(percentage_nifH_cluster < 5, paste0(percentage_nifH_cluster, "%"), ""),
      x = 1.5
    ),
    position = position_stack(vjust = 0.5),
    hjust = 0.5,
    size = 5,
    show.legend = FALSE,
    # nudge_y = 0.05
  )

# pie_chart_facet <- pie_chart +
#   facet_wrap("studyID")

# Display the pie chart
print(pie_chart_AUID_depud_DNA_rel_abund)
# print(pie_chart_facet)

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/pie_chart_nifhDF_DNA_dedup_perc_tot_rel_abund_nifH_cluster.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)

annoNifHDB_updt %>%
  filter(nifH_cluster == "1P") %>%
  distinct(CON)

### - combine plots

## * make new plots that are friendly to combining
pie1 <- pie_chart_AUID_depud_DNA_counts +
  theme(
    legend.position = "bottom", # Set the legend position to "bottom"
    legend.justification = c(1, 0), # Adjust the justification to place it at the bottom-right
    legend.box.just = "right"
  )

pie2 <- pie_chart_AUID_depud_DNA_rel_abund +
  theme(
    legend.position = "none"
  )

## * combine plots with patchwork
# combined_plot_pie_nifH_cluster_abundance <- pie1 + pie2
combined_plot_pie_nifH_cluster_abundance <- (pie1 | pie2)

print(combined_plot_pie_nifH_cluster_abundance)

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/pie_chart_nifhDF_DNA_dedup_perc_abund_nifH_cluster_cmb_plot.jpeg", height = 8.5, width = 16, units = "in", dpi = 300)


### - bar chart of all samples
### - bar chart of all samples
### - bar chart of all samples

#### * by counts
#### * just DNA, deduplicated by AUID and group_id
#### * percentages based on total samples

to_plot <- nifhdb_all_counts_AUID_dedup_total_study_id
# to_plot <- nifhdb_all_counts_AUID_dedup_total_study_id_clean

viridis_color_pallete <- get_viridis_colors(to_plot, nifH_cluster, "H", -1, 0)

(nifhdb_all_counts_AUID_dedup_total_study_id_bar_plot <- bar_plot(
  df = to_plot,
  x = studyID,
  y = percentage_total_counts_nifH_cluster_total_study_id,
  # fill = nifH_cluster_modified,
  fill = nifH_cluster,
  # fill_lab = expression(italic("nifH") "cluster"),
  fill_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
  x_lab = "Study ID",
  y_lab = "% of total nifH cluster counts per study ID",
  legend_position = "right"
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), legend.direction = "vertical"))

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/bar_chart_nifhDB_DNA_dedup_perc_tot_cnts_nifH_cluster_total_studyID.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)


####* by relative abundance
####* just DNA, deduplicated by AUID and group_id
####* percentages based on total samples

to_plot <- nifhdb_all_RA_AUID_dedup_total_study_id

viridis_color_pallete <- get_viridis_colors(to_plot, nifH_cluster, "H", -1, 0)

(nifhdb_all_RA_AUID_dedup_total_study_id_bar_plot <- bar_plot(
  df = to_plot,
  x = studyID,
  y = percentage_total_rel_abund_total_nifH_cluster_study_id,
  # y = percentage_total_rel_abund_total_nifH_cluster_study_id_clean,
  fill = nifH_cluster,
  # fill = nifH_cluster_modified,
  # fill_lab = expression(italic("nifH") "cluster"),
  fill_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
  x_lab = "Study ID",
  y_lab = "% of total nifH cluster relative abundance per study ID",
  legend_position = "right"
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), legend.direction = "vertical"))

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/bar_chart_nifhDB_DNA_dedup_perc_tot_RA_nifH_cluster_total_studyID.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)


#### by counts
#### just DNA, deduplicated by AUID and group_id
#### percentages recalculated based on individual studyID samples totals

to_plot <- nifhdb_all_counts_AUID_dedup_study_id_total

viridis_color_pallete <- get_viridis_colors(to_plot, nifH_cluster, "H", -1, 0)

(nifhdb_all_counts_AUID_dedup_study_id_total_bar_plot <- bar_plot(
  df = to_plot,
  x = studyID,
  y = percentage_total_counts_nifH_cluster_study_id_total,
  # y = percentage_total_rel_abund_total_nifH_cluster_study_id,
  fill = nifH_cluster,
  # fill_lab = expression(italic("nifH") "cluster"),
  fill_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
  x_lab = "Study ID",
  y_lab = "% of total nifH cluster counts per study ID",
  legend_position = "right"
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), legend.direction = "vertical"))

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/bar_chart_nifhDB_DNA_dedup_perc_tot_cnts_nifH_cluster_studyID_total.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)

#### by relative abundance
#### just DNA, deduplicated by AUID and group_id
#### percentages recalculated based on individual studyID samples totals

to_plot <- nifhdb_all_rel_abund_AUID_dedup_study_id_total

viridis_color_pallete <- get_viridis_colors(to_plot, nifH_cluster, "H", -1, 0)

(nifhdb_all_rel_abund_AUID_dedup_study_id_total_bar_plot <- bar_plot(
  df = to_plot,
  x = studyID,
  y = percentage_total_rel_abund_nifH_cluster_study_id_total,
  # y = percentage_total_rel_abund_total_nifH_cluster_study_id,
  fill = nifH_cluster,
  # fill_lab = expression(italic("nifH") "cluster"),
  fill_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
  x_lab = "Study ID",
  y_lab = "% of total nifH cluster relative abundance per study ID",
  legend_position = "right",
  legend_direction = "vertical"
)
+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), legend.direction = "vertical")
# +
#   guides(fill = guide_legend(
#     nrow = 1,
#     byrow = TRUE,
#     title.position = "top",
#     title.hjust = 0.5
# )
)

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/bar_chart_nifhDB_DNA_dedup_perc_tot_RA_nifH_cluster_studyID_total.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)


### - Make 4 panel plot
##* left side is counts
##* rights side is relative abundance

#* List all variables in the environment
all_vars <- ls()

#* Find variables containing "cluster" in their names and remove them
cat(all_vars[grep("_bar_plot", all_vars)], sep = "\n")


p1 <- nifhdb_all_counts_AUID_dedup_total_study_id_bar_plot +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 15, face = "bold")
  )

p2 <- nifhdb_all_RA_AUID_dedup_total_study_id_bar_plot +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 11.5, face = "bold")
  )

p3 <- nifhdb_all_counts_AUID_dedup_study_id_total_bar_plot +
  theme(
    legend.position = "none",
    # axis.text.x = element_blank(),
    axis.title.y = element_text(size = 15, face = "bold")
  )

p4 <- nifhdb_all_rel_abund_AUID_dedup_study_id_total_bar_plot +
  theme(
    legend.position = "none",
    # axis.text.x = element_blank(),
    axis.title.y = element_text(size = 11.5, face = "bold")
  )



(combined_plot_nifh_cluster_by_study_plot <- (p1 |
  p2) / (p3 | p4))

### * save plot
ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/bar_chart_nifhDB_DNA_dedup_perc_tot_nifH_cluster_studyID_cmb_plot.jpeg", height = 12.5, width = 14, units = "in", dpi = 300)

# combined_plot_nifh_cluster_by_study_plot <- (nifhdb_all_counts_AUID_dedup_total_study_id_bar_plot |
#   nifhdb_all_RA_AUID_dedup_total_study_id_bar_plot) / (nifhdb_all_counts_AUID_dedup_study_id_total_bar_plot | nifhdb_all_rel_abund_AUID_dedup_study_id_total_bar_plot)


# print(combined_plot_nifh_cluster_by_study_plot)
