# Function to define the custom theme
eb <- element_blank()

theme_custom <- function() {
  theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.major = eb,
      panel.grid.minor = eb,
      strip.background = eb,
      strip.text.x = element_text(size = 16, face = "bold"),
      panel.border = element_rect(colour = "black"),
      axis.title = element_text(size = 15, face = "bold"),
      axis.title.x = element_text(
        size = 23,
        face = "bold",
        colour = "black"
      ),
      axis.title.y = element_text(
        size = 21,
        face = "bold",
        colour = "black"
      ),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 17, face = "bold"),
      legend.box = "horizontal",
      legend.direction = "horizontal",
      axis.text.y = element_text(size = 21, face = "bold"),
      axis.text.x = element_text(size = 21, face = "bold")
    )
}

### custom colors
get_viridis_colors <- function(
    df_to_plot,
    plot_by,
    pallette_type = "D",
    color_direction,
    begin_point) {
  # Define the number of samples
  number_of_samples <- df_to_plot %>%
    distinct({{ plot_by }}) %>%
    nrow()
  # Create a custom color palette based on an existing one (e.g., viridis)
  custom_palette <- viridisLite::viridis(
    number_of_samples,
    option = pallette_type,
    direction = color_direction,
    begin = begin_point
  )

  return(custom_palette)
}

## colours specific to nifh clusters
# nifh_cluster_colours <- c(
#   "1A" = "steelblue",
#   "1J/1K" = "chocolate4",
#   "1O/1P" = "magenta",
#   # "3" = "goldenrod2",
#   "3" = "red",
#   "1G" = "chocolate1",
#   "1B" = "green2",
#   "other" = "lightgrey",
#   "unknown" = "darkgrey"
# )

nifh_cluster_colours <- c(
  "1A" = "steelblue",
  "1J/1K" = "chocolate4",
  "1O/1P" = "magenta",
  # "3" = "goldenrod2",
  "3" = "red",
  "1G" = "chocolate1",
  "1B" = "green2",
  "other" = "lightgrey",
  "unknown" = "darkgrey",
  "4" = "black",
  "2" = "blue"
)

# Function to create a custom dot plot
create_custom_plot_point <- function(
    data1, data2, x1, y1, x2, y2,
    colour, group, size, label, legend_position) {
  ggplot() +
    geom_point(
      data = data1,
      aes(x = {{ x1 }}, y = {{ y1 }}, colour = {{ colour }}, group = {{ group }}, size = {{ size }}),
      na.rm = TRUE,
      shape = 1,
      stroke = 1
    ) +
    geom_text(
      data = data2,
      aes(x = {{ x2 }}, y = {{ y2 }}, label = {{ label }}),
      size = 10,
      fontface = "bold",
      vjust = -1
    ) +
    labs(
      y = expression(bold(paste(
        "PO"[4]^"3-" * " (µmol L"^"-1", ")"
      ))),
      x = expression(bold(paste("SST (˚C)"))),
      title = "",
      shape = "hemisphere",
      fill = "Taxa"
    ) +
    scale_size_continuous(
      breaks = c(0, 0.1, 0.33, 0.66, 1.00),
      range = c(.1, 10),
      limits = c(0, 1)
    ) +
    scale_color_fermenter(
      type = "div",
      palette = "Spectral",
      direction = -1,
      breaks = c(-11, -9, -7)
    ) +
    ylim(-0.05, 1.95) +
    xlim(-2, 30.4) +
    theme_custom() +
    theme(legend.position = legend_position)
}

# Function to create a simple custom bar plot
bar_plot <- function(
    df, x, y,
    fill = NULL,
    fill_pallete = NULL,
    width = NULL,
    fill_lab = NULL,
    colour_lab = NULL,
    title_lab = NULL,
    subtitle_lab = NULL,
    x_lab, y_lab = "Number of samples",
    legend_position = "bottom",
    legend_direction = "horizontal") {
  gg <- ggplot({{ df }}, aes(x = {{ x }}, y = {{ y }}, fill = {{ fill }})) +
    geom_bar(stat = "identity", colour = "#3a3838", width = width) +
    labs(
      title = title_lab,
      subtitle = subtitle_lab,
      fill = fill_lab,
      colour = colour_lab,
      x = x_lab,
      y = y_lab
    ) +
    # theme_bw() +
    # theme_minimal() +
    scale_fill_manual(values = fill_pallete) +
    # scale_fill_manual(values = viridis_color_pallete) +
    theme_custom() +
    theme(
      legend.position = legend_position,
      legend.direction = legend_direction
    )
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
  # guides(fill = "none")
  return(gg)
}

# Function to create a simple custom histogram plot
# histogram_plot <- function(
#     df, y, fill = NULL, binwith = 1,
#     x_lab = "Number of samples", y_lab) {
#   # FIXME: gg <- ggplot({{ df }}, aes(y = {{ y }}, fill = as.factor({{ fill }}))) +
#   gg <- ggplot({{ df }}, aes(y = {{ y }}, fill = {{ fill }})) +
#     geom_histogram(binwidth = binwith, colour = "black") +
#     labs(
#       x = x_lab,
#       y = y_lab
#     ) +
#     # theme_bw() +
#     # theme_minimal() +
#     theme_custom() +
#     # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
#     guides(fill = "none")
#   return(gg)
# }

histogram_plot <- function(
    df, y, fill = NULL, binwith = 1,
    x_lab = "Number of samples", y_lab) {
  # FIXME: gg <- ggplot({{ df }}, aes(y = {{ y }}, fill = as.factor({{ fill }}))) +
  gg <- ggplot({{ df }}, aes(y = {{ y }}, fill = {{ fill }})) +
    geom_histogram(binwidth = binwith, colour = "black") +
    labs(
      x = x_lab,
      y = y_lab
    ) +
    # theme_bw() +
    # theme_minimal() +
    theme_custom() #+
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
  # guides(fill = "none")
  return(gg)
}

histogram_plot_x_or_y <- function(
    df, aes_var,
    x = NULL, y = NULL, fill = NULL, binwidth = 1,
    x_lab = "Number of samples", y_lab = "Number of samples") {
  aes_mapping <- if ({{ aes_var }} == "x") {
    aes(x = {{ x }}, fill = {{ fill }})
  } else if ({{ aes_var }} == "y") {
    aes(y = {{ y }}, fill = {{ fill }})
  } else {
    stop("Invalid value for aes_var. Use 'x' or 'y'.")
  }

  gg <- ggplot({{ df }}, aes_mapping) +
    geom_histogram(binwidth = binwidth, colour = "black") +
    labs(
      x = x_lab,
      y = y_lab
    ) +
    theme_custom() #+
  # guides(fill = "none")

  return(gg)
}




# # Function to create a simple custom histogram plot
# pie_chart <- ggplot(nifhdb_all_counts_AUID_dedup_clean, aes(
#   x = 1,
#   # x = studyID,
#   y = percentage_nifH_cluster,
#   fill = nifH_cluster,
#   # colour = studyID,
#   # group = nifH_cluster
# )) +
#   geom_bar(
#     stat = "identity",
#     # position = position_stack(reverse = FALSE)
#   ) +
#   # geom_histogram(binwidth = 1) +
#   coord_polar(theta = "y") +
#   theme_void() +
#   labs(title = "nifH dabtbase - % of total counts for each clusters", subtitle = "DNA only with replicates averaged", fill = "nifH Cluster") +

#   # Add labels with lines for small slices
#   geom_text(
#     # aes(label = ifelse(percentage_nifH_cluster >= 5, paste(nifH_cluster, "\n", percentage_nifH_cluster, "%"), ""), x = 1.25),
#     aes(label = ifelse(percentage_nifH_cluster >= 5, paste0(percentage_nifH_cluster, "%"), ""), x = 1.15),
#     position = position_stack(vjust = 0.5),
#     hjust = 0.5,
#     size = 12,
#     show.legend = FALSE
#   ) +
#   geom_text(
#     # aes(label = ifelse(percentage_nifH_cluster >= 5, paste(nifH_cluster, "\n", percentage_nifH_cluster, "%"), ""), x = 1.25),
#     aes(
#       label = ifelse(percentage_nifH_cluster < 5, paste0(percentage_nifH_cluster, "%"), ""),
#       x = 1.5
#     ),
#     position = position_stack(vjust = 0.5),
#     hjust = 0.5,
#     size = 5,
#     show.legend = FALSE,
#     # nudge_y = 0.05
#   )


# # * Usage

### create df to plot with features needed

# str(RA_df_T_lng)

# # Create the plot
# p <- create_custom_plot_point(
#   # data1 = nifHDB_lng_Clst,
#   # data1 = nifhDB_cnts_T_DNA_dedup_lng,
#   data1 = RA_df_T_lng  %>%
#   left_join(annoNifHDB_updt) %>%
#   left_join(CMAP_coloc),
#   data2 = plt_rgns,
#   x1 = sst_tblSST_AVHRR_OI_NRT,
#   y1 = PO4_tblPisces_NRT,
#   x2 = SSTs,
#   y2 = PO4s,
#   colour = logFe,
#   group = CON,
#   size = RA,
#   label = regions,
#   legend_position = "none"
# )

# print(p)


# getwd()

# # Add the date stamp to a string
# out_dir <- "~/mmorando@ucsc.edu_gDrive/My Drive/data/amplicon_review/all_studies/figures/"
# out_sub <- "station_plots/All/"
# out_nm <- "test"
# out_type <- ".jpeg"

# out_obj <- "_DNAstudies"
# #
# string_with_date <- paste(out_dir, out_sub, out_nm, out_obj, "_", date_stamp, out_type, sep = "")
# string_with_date



### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
cat("Done loading script!!!")
cat("Woooooooohooooooo!!!")
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
