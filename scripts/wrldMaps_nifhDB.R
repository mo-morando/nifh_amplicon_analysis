#!/usr/bin/env Rscript

library(tidyverse)

cmap_coloc <- read_csv("/Users/mo/Projects/nifH_amp_project/myWork/analysis/out_files/cmap_coloc.csv")

## - Load in the world map as a variable for plotting
world <- map_data("world")

#' Define a custom theme for ggplot2 plots
#'
#' This function defines a custom theme for ggplot2 plots, including modifications to plot titles, axes, grids, and legends.
#'
#' @return A ggplot2 theme object.
#'
theme_custom_maps <- function() {
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 13, face = "bold"),
    panel.border = element_rect(colour = "black"),
    axis.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 25, face = "bold", colour = "black"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 13, face = "bold"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold")
  )
}


### - load function for plotting world maps
create_custom_worldmap_plot <- function(
    data,
    x_var,
    y_var,
    color_var,
    legend_nrow,
    fill_pallete = NULL,
    legend_position = "bottom") {
  plot <- ggplot(data) +
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "black", fill = "black", size = 0.2
    ) +
    geom_point(
      aes(
        x = {{ x_var }}, y = {{ y_var }}, color = {{ color_var }}
      ),
      data = data # %>% mutate(studyID = str_replace_all(studyID, "_", " "))
      ,
      alpha = 0.15,
      size = 4, shape = 1, stroke = 3
    ) +
    # scale_color_manual(values = getPalette(palette_count)) +
    scale_fill_manual(values = fill_pallete) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      y = "Latitude",
      x = "Longitude",
      colour = "Study ID"
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank()) +
    theme_custom_maps() +
    theme(legend.position = legend_position) +
    guides(colour = guide_legend(nrow = legend_nrow, byrow = TRUE, override.aes = list(alpha = 1))) +
    coord_quickmap()

  return(plot)
}



# ### write function where you enter information in quotes
# testFF <- function(x, plt_flt) {
#   x %>%
#     filter(!!sym(plt_flt))
# }

### write function where you do not need quotes to enter information
filter_df <- function(x, plt_flt) {
  cat("Input Variables:\n")
  cat("x:", deparse(substitute(x)), "\n")
  cat("plt_flt:", deparse(substitute(plt_flt)), "\n")

  filtered_data <- x %>%
    filter({{ plt_flt }})

  return(filtered_data)
}

####* cmap_coloc has 871 stations so its the nifHDB level, after FilterAUIDs

## * only photic
data_sub <- cmap_coloc
data_sub <- filter_df(x = cmap_coloc, plt_flt = photic)

## * only aphotic
# data_sub <- testFF(x = cmap_coloc, plt_flt = !photic)

## * look at specific study ids
cmap_coloc %>%
  distinct(studyID) %>%
  view()

# data_sub <- filter_df(x = cmap_coloc, plt_flt = studyID == "Shiozaki_2018LNO")
# data_sub <- filter_df(x = cmap_coloc, plt_flt = studyID == "TurkKubo_2021")
# data_sub <- filter_df(x = cmap_coloc, plt_flt = studyID == "Shiozaki_2017")
# data_sub <- filter_df(x = cmap_coloc, plt_flt = studyID == "Shiozaki_2020")
# data_sub <- filter_df(x = cmap_coloc, plt_flt = studyID == "Shiozaki_2018GBC")
# data_sub <- filter_df(x = cmap_coloc, plt_flt = studyID == "NEMO")
# data_sub <- filter_df(x = cmap_coloc, plt_flt = studyID == "Mulholland_2018")
# data_sub <- filter_df(x = cmap_coloc, plt_flt = ocean == "Atlantic")





# Setting up pallete with color blind safe colors
# Define the base palette
cbpal_15 <- c(
  "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", "#490092", "#006ddb", "#b66dff", "#6db6ff",
  "#b6dbff", "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d"
)

# Function to generate color palette
generate_cb_palette <- function(n) {
  # If n is less than or equal to 15, return the first n colors
  if (n <= length(cbpal_15)) {
    return(cbpal_15[1:n])
  }
  # If n is greater than 15, use colorRampPalette to interpolate and create more colors
  return(colorRampPalette(cbpal_15)(n))
}


cb_pallette <- generate_cb_palette(nrow(data_sub  %>% distinct(studyID)))



# Call the function with your desired variables and parameters
(world_map_plot <- create_custom_worldmap_plot(
  # data = mapping_df_sub,
  # data = mapping_df_sub %>% distinct(SAMPLEID, .keep_all = TRUE),
  data = data_sub,
  x_var = lon,
  y_var = lat,
  color_var = studyID,
  legend_nrow = 3,
  fill_pallete = cb_pallette,
  legend_position = "bottom"
) %>%
  suppressWarnings())


world_map_plot +
  # facet_wrap(~studyID) +
  facet_wrap(~ocean) +
  theme(legend.position = "bottom")

# Print the plot
print(world_map_plot)

# Print a statement with the parameters
cat("Custom plot created with the following parameters:\n")
cat("x_var:", "lon\n")
cat("y_var:", "lat\n")
cat("color_var:", "studyID\n")
cat("palette_count:", colourCount, "\n")


ggsave("/Users/mo/Projects/nifH_amp_project/myWork/analysis/plots/worldStnMap_nifhdb_allstudies.svg", height = 10.4, width = 17, units = "in", dpi = 300, device = "svg")



# ggsave("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/figures/worldMaps/metadata/worldStnMap_nifhdb_allstudies.jpeg", height = 8.5, width = 14)
# ggsave("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/figures/worldMaps/metadata/worldStnMap_nifhdb_DNAstudies.jpeg", height = 8.5, width = 14)
# ggsave("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/figures/worldMaps/metadata/worldStnMap_nifhdb_RNAstudies.jpeg", height = 8.5, width = 14)
# ggsave("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/figures/worldMaps/metadata/worldStnMap_nifhdb_photicAllstudies.jpeg", height = 8.5, width = 14)
# ggsave("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/figures/worldMaps/metadata/worldStnMap_nifhdb_aphoticAllstudies.jpeg", height = 8.5, width = 14)






### Data to look at

mapping_df <- nifhDB_RA

## searches

## 1G
# clust_1G = "nifH_cluster == '1G' & Values > 0"
# alpha_beta = "nifH_cluster == '1G' & nifH_cluster == '1G' &Values > 0"

mapping_df_sub <- mapping_df %>%
  pivot_longer(cols = -AUID, names_to = "SAMPLEID", values_to = "Values") %>%
  merge_annotations() %>%
  filter_df(plt_flt = nifH_cluster == "1O/1P" & Values > 0) %>%
  # pivot_wider(names_from = AUID, values_from = "Values") %>%
  merge_cmap()

mapping_df_sub %>%
  distinct(ocean)

## if using viridis to select colors
viridis_color_pallete <- get_viridis_colors(mapping_df_sub, studyID, "H", -1, 0)
# viridis_color_pallete <- get_viridis_colors(cmap_coloc, studyID, "H", -1, 0)

# Call the function with your desired variables and parameters
world_map_plot <- create_custom_worldmap_plot(
  # data = mapping_df_sub,
  data = mapping_df_sub %>% distinct(SAMPLEID, .keep_all = TRUE),
  x_var = lon,
  y_var = lat,
  color_var = ocean,
  legend_position = "bottom"
)

# Assuming your data frame is named 'your_data_frame_name'
na_instances <- cmap_coloc %>%
  select(where(~ any(is.na(.)))) ## checks by columns for any NAs <- this is way more useful
filter_all(any_vars(is.na(.))) ## checks by row

# 'na_instances' will contain all rows where at least one column has an NA value.
