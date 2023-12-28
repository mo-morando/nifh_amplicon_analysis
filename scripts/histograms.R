### Not sure how to display this data
## But we need to know how many or the percentage of certain parts of our data


library(patchwork)



### photic samples
number_of_aphotic_samples <- nrow(CMAP_coloc) - nrow(remove_aphotic_samples
(CMAP_coloc))

number_of_aphotic_samples <- nrow(CMAP_20230719) - nrow(remove_aphotic_samples
(CMAP_coloc))

number_of_photic_samples <- nrow(CMAP_coloc) -
  number_of_aphotic_samples

cat("

    Total number of samples: ", nrow(CMAP_coloc), "
    Number of aphotic samples:", number_of_aphotic_samples, "
    Number of photic samples:", number_of_photic_samples,
  sep = " "
)
### DNA/RNA samples

number_of_DNA_samples <- nrow(CMAP_20230719) - nrow(remove_samples_nucleic_acid(CMAP_coloc, "RNA", DNA_samples_key))

cat("

    Number of DNA samples:", number_of_DNA_samples,
  sep = " "
)

number_of_RNA_samples <- nrow(CMAP_20230719) - nrow(remove_samples_nucleic_acid(CMAP_coloc, "DNA", DNA_samples_key))

cat("

Number of RNA samples:", number_of_RNA_samples,
  sep = " "
)

### duplicate samples
number_of_duplicate_samples <- unique_sample_id_key %>%
  filter(sample_type %in% "Duplicate_Samples") %>%
  count() %>%
  pull()


cat("Number of duplicate samples:", number_of_duplicate_samples,
  sep = " "
)





## - nucleicAcidType

query_df <- dedup_by_group(CMAP_coloc, group_id)

samples_per_nucleicAcidType <- count_and_arrange(query_df, "nucleicAcidType") %>%
  mutate(
    total = sum(n),
    percentage = n / total * 100
  ) %>%
  select(nucleicAcidType, n, percentage)

print(samples_per_nucleicAcidType, n = 1000)

## * by study id
dedup_by_group(CMAP_coloc, group_id) %>%
  # distinct(studyID, nucleicAcidType) %>%
  count_and_arrange(c("studyID", "nucleicAcidType")) %>%
  add_total_row(
    column_name = "studyID",
    summary_column = n,
    # pull_name = n,
    all_columns = FALSE
  )

### -  photic
samples_per_photic <- count_and_arrange(query_df, "photic") %>%
  mutate(
    total = sum(n),
    percentage = n / total * 100
  ) %>%
  select(photic, n, percentage)

print(samples_per_photic, n = 1000)

count_and_arrange(CMAP_coloc, c("photic", "nucleicAcidType"))


## - sample type
(sample_type <- count_and_arrange(sample_types_all, c("sample_type", "nucleicAcidType"), sample_type) %>%
  # add_total_row(n, "sample_type", all_columns = FALSE) %>%
  mutate(
    new_group = paste(sample_type, nucleicAcidType, sep = "_"),
    tibble_id = "sample_id"
  ) %>%
  add_percentage(n,
    percentage,
    grouping_by = NULL,
    remove_columns = "total"
  ))

### * Plot
to_plot <- sample_type


viridis_color_pallete <- get_viridis_colors(to_plot, new_group, "H", -1, 0)

bar_plot(
  df = to_plot,
  x = percentage,
  y = "sample type",
  # y = percentage_total_rel_abund_total_nifH_cluster_study_id_clean,
  fill = new_group,
  # fill = nifH_cluster_modified,
  # fill_lab = expression(italic("nifH") "cluster"),
  # fill_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
  y_lab = NULL,
  x_lab = "% of total",
  legend_position = "right"
) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 19,
      face = "bold"
    ),
    legend.direction = "vertical"
  )

### * combine nucleic acid type, light level, and sample type so you can plot them all together

##
combined_nuca_phtc_tibble <- samples_per_nucleicAcidType %>%
  mutate(tibble_id = "nucleicAcidType") %>%
  rename(new_group = nucleicAcidType) %>%
  bind_rows(samples_per_photic %>%
    mutate(
      tibble_id = "light level",
      photic = if_else(photic %in% "TRUE", "photic", "aphotic")
    ) %>%
    rename(new_group = photic))


combined_nuca_phtc_smptp_tibble <- sample_type %>%
  # mutate(var_type = "sample_id") %>%
  # rename(tibble_id = sample_type) %>%
  bind_rows(combined_nuca_phtc_tibble)


to_plot <- combined_nuca_phtc_smptp_tibble

viridis_color_pallete <- get_viridis_colors(to_plot, new_group, "H", -1, 0)

(combined_nuca_phtc_smptp_tibble_plot <- bar_plot(
  df = to_plot,
  y = tibble_id,
  x = percentage,
  # y = percentage_total_rel_abund_total_nifH_cluster_study_id_clean,
  fill = new_group,
  # fill = nifH_cluster_modified,
  # fill_lab = expression(italic("nifH") "cluster"),
  # fill_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
  y_lab = NULL,
  x_lab = "% of total",
  legend_position = "right"
) +
  theme(
    axis.text.y = element_text(
      angle = 45,
      hjust = 1,
      size = 19,
      face = "bold"
    ),
    legend.direction = "vertical"
  )
)


ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/samp_nucleic_acid_type_light_bar_cmb_plot.jpeg", height = 10.5, width = 14, units = "in", dpi = 300)

### * combine plot

samples_per_studyid_plot

combined_nuca_phtc_smptp_tibble_plot

(combine_test <- samples_per_studyid_plot | combined_nuca_phtc_smptp_tibble_plot +
  plot_layout(widths = c(1.5, 0))
#  plot_layout(widths = c(2, 1), heights = unit(c(5, 1), c('cm', 'null')))
)

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/bar_chart_nifhDB_DNA_dedup_perc_tot_nifH_cluster_studyID_samples_per_cmb_plot.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)



# TODO:
### information for barplots and histograms


### Study ID
# samples_per_studyid <- query_df %>%
#     # distinct(SAMPLEID, .keep_all = T) %>%
#     group_by(studyID, hemi) %>%
#     # group_by(lat_abs, hemi) %>%
#     # group_by(lat) %>%
#     # count() %>%
#     summarize(sample_count = n()) %>%
#     ungroup() %>%
#     arrange(desc(n))

(samples_per_studyid <- count_and_arrange(
  query_df, c("studyID")
))

# print(samples_per_studyid, n = 50)

to_plot <- samples_per_studyid

viridis_color_pallete <- get_viridis_colors(to_plot, studyID, "H", -1, 0)

(samples_per_studyid_plot <- bar_plot(
  df = to_plot,
  y = studyID,
  x = n,
  fill = NULL,
  y_lab = "Study ID",
  x_lab = "Number of samples",
  legend_position = "none"
) #+
# theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15))
)

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_studyID_bar.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)

### month
(samples_per_month <- count_and_arrange(query_df, c("month")))

# query_df  %>%
# filter(grepl("-12-", x = time))

to_plot <- samples_per_month

viridis_color_pallete <- get_viridis_colors(to_plot, month, "H", -1, 0)

bar_plot(
  df = to_plot,
  y = month,
  x = n,
  fill = month,
  y_lab = "Month",
  x_lab = "Number of samples",
  legend_position = "none",
  legend_direction = "vertical"
) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
  scale_y_discrete(
    limits =
      rev(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")),
    labels =
    # c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
      rev(c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")),
  )


ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_month_bar.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)

### month by study ID
(samples_per_month_study_id <- count_and_arrange(query_df, c("month", "studyID")))

# query_df  %>%
# filter(grepl("-12-", x = time))

to_plot <- samples_per_month_study_id

viridis_color_pallete <- get_viridis_colors(to_plot, studyID, "H", -1, 0)

bar_plot(
  df = to_plot,
  y = month,
  x = n,
  fill = studyID,
  y_lab = "Month",
  x_lab = "Number of samples",
  legend_position = "bottom"
) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
  scale_y_discrete(
    limits =
      rev(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")),
    labels =
    # c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
      rev(c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")),
  )


ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_month_study_id_bar.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)




### - Histograms

### * lat_abs
samples_per_lat_abs <- count_and_arrange(query_df, c("lat_abs", "hemi"))

print(samples_per_lat_abs, n = 1000)

# histogram_plot(samples_per_lat_abs, lat_abs,
#     hemi, 1,
#     y_lab = "absolute latitude"
# )

histogram_plot_x_or_y(samples_per_lat_abs, "x", lat_abs,
  fill = hemi, binwidth = 1,
  x_lab = "absolute latitude"
) +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    fill = "Hemisphere"
  ) +
  scale_fill_manual(
    values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
    label = c("northern", "southern")
  )

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_abslat_hemi_hist.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)


### season
samples_per_season <- count_and_arrange(
  query_df, c("season", "hemi")
)

# samples_per_season <- count_and_arrange(
#     CMAP_coloc, c("season", "hemi")
# )

# print(samples_per_season, n = 1000)

bar_plot(
  df = samples_per_season,
  x = n,
  y = season,
  fill = hemi,
  y_lab = "season",
  x_lab = "Number of samples"
) +
  # theme(
  #     legend.position = "bottom"
  # ) +
  labs(
    fill = "Hemisphere"
  ) +
  scale_fill_manual(
    values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
    label = c("northern", "southern")
  )

## !x value has to be continous so "season" doesn't work
## !could change to number c(1,2,3,4) and then label them by their actual season

# histogram_plot_x_or_y(
#     df = samples_per_season, aes_var = "x",
#     x = season,
#     fill = hemi,
#     binwidth = 1,
#     x_lab = "season"
# ) +
#     theme(
#         legend.position = "bottom"
#     ) +
#     labs(
#         fill = "Hemisphere"
#     ) +
#     scale_fill_manual(
#         values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
#         label = c("northern", "southern")
#     )

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_season_hemi_bar.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)


### depth
samples_per_depth <- count_and_arrange(
  query_df, c("depth", "hemi")
) %>%
  mutate(less_than_150 = depth <= 155)


# samples_per_depth <- count_and_arrange(
#     CMAP_coloc, c("depth", "hemi")
# )

# print(samples_per_depth, n = 1000)
# print(samples_per_depth %>%
#   arrange(depth), n = 1000)

# histogram_plot_x_or_y(
#     df = samples_per_depth, aes_var = "x",
#     x = depth,
#     fill = hemi,
#     binwidth = 100,
#     x_lab = "depth (m)"
# ) +
#     theme(
#         legend.position = "bottom"
#     ) +
#     labs(
#         fill = "Hemisphere"
#     ) +
#     scale_fill_manual(
#         values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
#         label = c("northern", "southern")
#     )

to_plot <- samples_per_depth

# viridis_color_pallete <- get_viridis_colors(to_plot, studyID, "A", 1, 0)

bar_plot(
  df = to_plot,
  y = n,
  x = depth,
  fill = hemi,
  x_lab = "depth (m)",
  y_lab = "Number of samples",
  width = 1,
  # legend_position = "none"
) +
  scale_fill_manual(
    values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
    label = c("northern", "southern")
  ) +
  coord_flip() +
  scale_x_reverse() +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(
    ~less_than_150,
    ncol = 1,
    scales = "free"
  )


### * just depths above 150 m
(depths_above_150m_bar_plt <- bar_plot(
  df = samples_per_depth %>%
    filter(less_than_150 == T),
  y = n,
  x = depth,
  fill = hemi,
  x_lab = "depth (m)",
  y_lab = "Number of samples",
  # legend_position = "none",
  width = 1,
  legend_position = "none",
  # legend_direction = "horizontal"
) +
  scale_fill_manual(
    values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
    label = c("northern", "southern")
  ) +
  # theme(
  #   axis.title.x = element_blank()
  # ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_reverse() +
  coord_flip()
)

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_depth_hemi_above150m_hist.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)

### * just depths below 150 m
(depths_below_150m_bar_plt <- bar_plot(
  df = samples_per_depth %>%
    filter(less_than_150 == FALSE),
  y = n,
  x = depth,
  fill = hemi,
  x_lab = "depth (m)",
  # legend_position = "none",
  width = 30,
  legend_position = "right",
  legend_direction = "vertical"
) +
  scale_fill_manual(
    values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
    label = c("northern", "southern")
  ) +
  # scale_x_reverse(
  #   limits = c(150, 3015),
  #   breaks = c(150, 300, 500, 750, 1000, 2000, 3000)
  # ) +
  scale_x_reverse(
    limits = c(3015, 150),
    breaks = c(150, 300, 500, 750, 1000, 2000, 3000)
  ) +
  # scale_x_continuous(
  #   limits = rev(c(150, 3015)),
  #   breaks = rev(c(150, 300, 500, 750, 1000, 2000, 3000))
  # ) +
  scale_y_continuous(
    limits = c(0, 10),
    expand = c(0, 0)
  ) +
  coord_flip()
)

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_depth_hemi_below150_hist.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)

### * combine plots

combined_plot_depths <- depths_above_150m_bar_plt | depths_below_150m_bar_plt

print(combined_plot_depths)

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_depth_hemi_hist_cmb_plot.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)



### _ oceanographic data
### _ oceanographic data
### _ oceanographic data

###* sst_tblSST_AVHRR_OI_NRT
samples_per_sst_tblSST_AVHRR_OI_NRT <- count_and_arrange(
  query_df, c("sst_tblSST_AVHRR_OI_NRT", "hemi")
)

print(samples_per_sst_tblSST_AVHRR_OI_NRT, n = 1000)

# histogram_plot(
#     df = samples_per_sst_tblSST_AVHRR_OI_NRT,
#     y = sst_tblSST_AVHRR_OI_NRT,
#     fill = hemi,
#     binwith = 1,
#     y_lab = "SST ˚C"
# )

# histogram_plot_x_or_y(
#   df = samples_per_sst_tblSST_AVHRR_OI_NRT, aes_var = "x",
#   x = sst_tblSST_AVHRR_OI_NRT,
#   fill = hemi,
#   binwidth = 1,
#   x_lab = "SST ˚C"
# ) +
#   theme(
#     legend.position = "bottom"
#   ) +
#   labs(
#     fill = "Hemisphere"
#   ) +
#   scale_fill_manual(
#     values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
#     label = c("northern", "southern")
#   )

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_SST_hemi_hist.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)


### * PO4_tblPisces_NRT
samples_per_PO4_tblPisces_NRT <- count_and_arrange(
  query_df, c("PO4_tblPisces_NRT", "hemi")
)

print(samples_per_PO4_tblPisces_NRT, n = 1000)

# histogram_plot(
#     df = samples_per_PO4_tblPisces_NRT,
#     y = PO4_tblPisces_NRT,
#     fill = hemi,
#     binwith = 1,
#     y_lab = "SST ˚C"
# )

histogram_plot_x_or_y(
  df = samples_per_PO4_tblPisces_NRT, aes_var = "x",
  x = PO4_tblPisces_NRT,
  fill = hemi,
  binwidth = 0.01,
  x_lab = expression(bold(paste(
    "PO"[4]^"3-" * " (µmol L"^"-1", ")"
  )))
) +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    fill = "Hemisphere"
  ) +
  scale_fill_manual(
    values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
    label = c("northern", "southern")
  )

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_PO4_hemi_hist.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)


### * logFe
samples_per_logFe <- count_and_arrange(
  query_df, c("logFe", "hemi")
)

print(samples_per_logFe, n = 1000)

histogram_plot_x_or_y(
  df = samples_per_logFe, aes_var = "x",
  x = logFe,
  fill = hemi,
  binwidth = 0.1,
  x_lab = expression(bold(paste(
    "log(Fe)"
  )))
) +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    fill = "Hemisphere"
  ) +
  scale_fill_manual(
    values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
    label = c("northern", "southern")
  )

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_logFe_hemi_hist.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)


### * PP_tblPisces_NRT
samples_per_PP_tblPisces_NRT <- count_and_arrange(
  query_df, c("PP_tblPisces_NRT", "hemi")
)

print(samples_per_PP_tblPisces_NRT, n = 1000)

histogram_plot_x_or_y(
  df = samples_per_PP_tblPisces_NRT, aes_var = "x",
  x = PP_tblPisces_NRT,
  fill = hemi,
  binwidth = 0.001,
  x_lab = expression(bold(paste(
    "PP"
  )))
) +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    fill = "Hemisphere"
  ) +
  scale_fill_manual(
    values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
    label = c("northern", "southern")
  )

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_PP_tblPisces_NRT_hemi_hist.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)


### * CHL_tblPisces_NRT
samples_per_CHL_tblPisces_NRT <- count_and_arrange(
  query_df, c("CHL_tblPisces_NRT", "hemi")
)

print(samples_per_CHL_tblPisces_NRT, n = 1000)

histogram_plot_x_or_y(
  df = samples_per_CHL_tblPisces_NRT, aes_var = "x",
  x = CHL_tblPisces_NRT,
  fill = hemi,
  binwidth = 0.1,
  x_lab = expression(bold(paste(
    "chlorophyll a"
  )))
) +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    fill = "Hemisphere"
  ) +
  scale_fill_manual(
    values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
    label = c("northern", "southern")
  )

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_CHL_tblPisces_NRT_hemi_hist.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)
