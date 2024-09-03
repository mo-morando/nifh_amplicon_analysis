### Not sure how to display this data
## But we need to know how many or the percentage of certain parts of our data


library(patchwork)

# remake query df
query_df <- dedup_by_group(cmap_coloc, group_id)


### * Plot
to_plot <- sample_type


viridis_color_pallete <- get_viridis_colors(to_plot, new_group, "H", -1, 0)

bar_plot(
  df = to_plot,
  x = percentage,
  y = "sample type",
  # y = percentage_total_rel_abund_total_nifH_cluster_study_id_clean,
  fill = new_group,
  fill_pallete = viridis_color_pallete,
  # fill = nifH_cluster_modified,
  # fill_lab = expression(italic("nifH") "cluster"),
  # fill_lab = bquote(bold(bold(italic(nifH)) ~ cluster)),
  y_lab = NULL,
  x_lab = "% of total",
  legend_position = "right",
  x_axis_angle = TRUE,
  n_row = 10,
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

(combined_nuca_phtc_smptp_tibble_plot <- suppressMessages(bar_plot(
  df = to_plot,
  y = tibble_id,
  x = percentage,
  # y = percentage_total_rel_abund_total_nifH_cluster_study_id_clean,
  fill = new_group,
  fill_pallete = viridis_color_pallete,
  # fill = nifH_cluster_modified,
  # fill_lab = bquote((bold(italic(nifH)) ~ cluster)),
  y_lab = NULL,
  x_lab = "% of total",
  legend_position = "right",
  n_row = 10,
  x_axis_angle = TRUE
) ))


ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/samp_nucleic_acid_type_light_bar_cmb_plot.jpeg", height = 10.5, width = 14, units = "in", dpi = 300)



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
  query_df, c("studyID", "nucleicAcidType", "ocean")
))

# print(samples_per_studyid, n = 50)

to_plot <- samples_per_studyid

viridis_color_pallete <- get_viridis_colors(to_plot, ocean, "magma", -1, 0)

(samples_per_studyid_plot <- bar_plot(
  df = to_plot,
  y = studyID,
  x = n,
  fill = ocean,
  fill_pallete = viridis_color_pallete,
  y_lab = "Study ID",
  x_lab = "Number of samples",
  legend_position = "right",
  legend_direction = "vertical",
  n_row = 5,
  print_out = TRUE
) #+
# theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15))
)

samples_per_studyid_plot +
  facet_wrap(~nucleicAcidType)
# facet_wrap(~hemi)

# ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_studyID_nucacid_fill_bar.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_studyID_ocean_fill_bar_facet_nucelic.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)


# ### * combine plot

# samples_per_studyid_plot

# combined_nuca_phtc_smptp_tibble_plot

# (combine_test <- samples_per_studyid_plot | combined_nuca_phtc_smptp_tibble_plot +
#   plot_layout(widths = c(1.5, 0))
# #  plot_layout(widths = c(2, 1), heights = unit(c(5, 1), c('cm', 'null')))
# )

# ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/bar_chart_nifhDB_DNA_dedup_perc_tot_nifH_cluster_studyID_samples_per_cmb_plot.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)



### -  Ocean and hemisphere
(samples_per_ocean <- count_and_arrange(
  #   query_df, c("studyID", "nucleicAcidType")
  #   query_df, c("ocean")
  query_df, c("ocean", "hemi")
))
samples_per_ocean %>%
  left_join(
    count_and_arrange(query_df, c("ocean"), "n_ocean", "n_ocean")
  ) %>%
  arrange(desc(ocean))
# print(samples_per_studyid, n = 50)

299 / 429 * 100

to_plot <- samples_per_ocean

viridis_color_pallete <- get_viridis_colors(to_plot, hemi, "magma", -1, 0)

(samples_per_ocean_plot <- suppressMessages(bar_plot(
  df = to_plot,
  y = ocean,
  x = n,
  fill = hemi,
  fill_pallete = viridis_color_pallete,
  y_lab = "Oceans",
  x_lab = "Number of samples",
  fill_lab = "Hemisphere",
  legend_position = "right",
  legend_direction = "vertical",
  print_out = TRUE
) +
  scale_fill_manual(
    values = c("northernHemi" = "darkorange", "southernHemi" = "green"),
    label = c("northern", "southern")
  ))
# theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15))
)

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_ocean_hemi_fill_bar.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)

### month
(samples_per_month <- count_and_arrange(query_df, c("month")))

# query_df  %>%
# filter(grepl("-12-", x = time))

to_plot <- samples_per_month

viridis_color_pallete <- get_viridis_colors(to_plot, month, "H", -1, 0)

(samples_per_month_plot <- suppressMessages(bar_plot(
  df = to_plot,
  y = month,
  x = n,
  fill = NULL,
  # fill_pallete = viridis_color_pallete,
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
  )))


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
  fill_pallete = viridis_color_pallete,
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

# print(samples_per_lat_abs, n = 1000)

# histogram_plot(samples_per_lat_abs, lat_abs,
#     hemi, 1,
#     y_lab = "absolute latitude"
# )

(samples_per_lat_abs_plot <- suppressMessages(histogram_plot_x_or_y(samples_per_lat_abs, "x", lat_abs,
  fill = hemi, binwidth = 2,
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
  )))

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_abslat_hemi_hist.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)


### season
samples_per_season <- count_and_arrange(
  query_df, c("season", "hemi")
)

# samples_per_season <- count_and_arrange(
#     cmap_coloc, c("season", "hemi")
# )

# print(samples_per_season, n = 1000)

suppressMessages((bar_plot(
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
  )))

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
#     cmap_coloc, c("depth", "hemi")
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

suppressMessages(bar_plot(
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
  ))


### * just depths above 150 m
(depths_above_150m_bar_plt <- suppressMessages(bar_plot(
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
  coord_flip()))

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_depth_hemi_above150m_hist.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)

### * just depths below 150 m
(depths_below_150m_bar_plt <- suppressMessages(bar_plot(
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
  coord_flip()))

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_depth_hemi_below150_hist.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)

### * combine plots

combined_plot_depths <- depths_above_150m_bar_plt | depths_below_150m_bar_plt

# print(combined_plot_depths)

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_depth_hemi_hist_cmb_plot.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)



### _ oceanographic data
### _ oceanographic data
### _ oceanographic data

###* sst_tblSST_AVHRR_OI_NRT
samples_per_sst_tblSST_AVHRR_OI_NRT <- count_and_arrange(
  query_df, c("sst_tblSST_AVHRR_OI_NRT", "hemi")
)

# print(samples_per_sst_tblSST_AVHRR_OI_NRT, n = 1000)

# histogram_plot(
#     df = samples_per_sst_tblSST_AVHRR_OI_NRT,
#     y = sst_tblSST_AVHRR_OI_NRT,
#     fill = hemi,
#     binwith = 1,
#     y_lab = "SST ˚C"
# )

sst_hist <- suppressMessages(histogram_plot_x_or_y(
  df = samples_per_sst_tblSST_AVHRR_OI_NRT, aes_var = "x",
  x = sst_tblSST_AVHRR_OI_NRT,
  fill = hemi,
  binwidth = 1,
  x_lab = "SST ˚C"
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
  ))

print(sst_hist)

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_SST_hemi_hist.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)

### * PO4_tblPisces_NRT
samples_per_PO4_tblPisces_NRT <- count_and_arrange(
  query_df, c("PO4_tblPisces_NRT", "hemi")
)

# print(samples_per_PO4_tblPisces_NRT, n = 1000)

# histogram_plot(
#     df = samples_per_PO4_tblPisces_NRT,
#     y = PO4_tblPisces_NRT,
#     fill = hemi,
#     binwith = 1,
#     y_lab = "SST ˚C"
# )

po4_hist <- suppressMessages(histogram_plot_x_or_y(
  df = samples_per_PO4_tblPisces_NRT, aes_var = "x",
  x = PO4_tblPisces_NRT,
  fill = hemi,
  binwidth = 0.02,
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
  ))

print(po4_hist)

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_PO4_hemi_hist.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)



## - Make combined plot
#* # absolute lat ---> samples_per_lat_abs_plot
#* # ocean  ---> samples_per_ocean_plot
#* # SST ---> sst_hist
#* # PO4 ---> po4_hist

#* # alter plots for combining

sst_hist_sub <- sst_hist + theme(legend.position = "none")
po4_hist_sub <- po4_hist + theme(legend.position = "none")
samples_per_lat_abs_plot_sub <- samples_per_lat_abs_plot + theme(legend.position = "none")
samples_per_ocean_plot_sub <- samples_per_ocean_plot + theme(
  legend.position = "none",
  axis.text.y = element_text(
    angle = 45,
    hjust = 1,
    size = 19,
    face = "bold"
  )
)

#* # combine with patchwork
combined_plot_fig_5 <- (samples_per_lat_abs_plot_sub / sst_hist_sub) | (samples_per_ocean_plot_sub / po4_hist_sub)

print(combined_plot_fig_5)

# (combined_plot_fig_5 <- (samples_per_lat_abs_plot_sub | samples_per_ocean_plot_sub) / (sst_hist_sub | po4_hist_sub))

# (combined_plot_fig_5 <- samples_per_lat_abs_plot_sub / samples_per_ocean_plot_sub + sst_hist_sub / po4_hist_sub)

# combined_plot_fig_5 + plot_layout(
#   ncol = 2,
#   heights = c(1, 2), # This adjusts the height ratio between the top and bottom plots
#   # widths = c(1, 5)  # This adjusts the height ratio between the top and bottom plots
# )

ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_abslat_ocean_sst_po4_hemi_cmb_hist.svg", height = 8.5, width = 17, units = "in", dpi = 300, device = "svg")


# ## combine plots   OLD PLOT HERE

# sst_hist <- sst_hist + theme(legend.position = "none")

# combine_sst_po4_hist <- (sst_hist / po4_hist)

# print(combine_sst_po4_hist)

# ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/Samples_per_PO4_N_SST_hemi_cmb_hist.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)




### * logFe
samples_per_logFe <- count_and_arrange(
  query_df, c("logFe", "hemi")
)

# print(samples_per_logFe, n = 1000)

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

# print(samples_per_PP_tblPisces_NRT, n = 1000)

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

# print(samples_per_CHL_tblPisces_NRT, n = 1000)

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
