### * Make 4 panel plot
## left side is counts
## rights side is relative abundance

# List all variables in the environment
all_vars <- ls()

# Find variables containing "cluster" in their names and remove them
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
ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/bar_chart_nifhDB_DNA_dedup_perc_tot_cnts_nifH_cluster_studyID_total.jpeg", height = 8.5, width = 14, units = "in", dpi = 300)

combined_plot_nifh_cluster_by_study_plot <- (nifhdb_all_counts_AUID_dedup_total_study_id_bar_plot |
    nifhdb_all_RA_AUID_dedup_total_study_id_bar_plot) / (nifhdb_all_counts_AUID_dedup_study_id_total_bar_plot | nifhdb_all_rel_abund_AUID_dedup_study_id_total_bar_plot)


print(combined_plot_nifh_cluster_by_study_plot)
