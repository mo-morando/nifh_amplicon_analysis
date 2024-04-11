### updated nifH dataDB

## !! I am not sure if this is the correct file I should use or if I should use
## !! /Users/mo/Library/CloudStorage/GoogleDrive-mmorando@ucsc.edu/My Drive/scripts/amp_review_Specific/stn_plots/STTPFe/SSTPFe_scatter.R

#### no faceted
data <- tibble(
    regions = c("poles", "high latitidue", "subpolar gyre", "subtropical gyre", "transition zone", "equator"),
    SSTs = c(1, 6, 12, 22, 16, 26),
    PO4s = c(0.75, 1.52, 0.625, 0.45, 0.125, 1.125)
)

plt_rgns <- tibble(
    regions = c("poles", "high latitidue", "subpolar gyre", "subtropical gyre", "transition zone", "equator"),
    SSTs = c(2, 6, 12, 24, 14.5, 26),
    PO4s = c(0.75, 1.52, 0.625, 0.35, -0.05, 1.0)
)

#### facetted
data_GammaAfct <- tibble(
    regions = c("poles", "high latitidue", "subpolar gyre", "subtropical gyre", "transition zone", "equator"),
    SSTs = c(2, 6, 12, 24, 14.5, 26),
    PO4s = c(0.75, 1.52, 0.625, 0.35, -0.05, 1.0)
)



### functions

plt_thm <- function() {
    theme(
        plot.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 16, face = "bold"),
        # titles of each facet
        panel.border = element_rect(colour = "black"),
        # axis.ticks.x = element_blank(),
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
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 17, face = "bold"),
        legend.position = "none",
        legend.box = "horizontal",
        legend.direction = "horizontal",
        # legend.position = 'none',
        axis.text.y = element_text(size = 21, face = "bold"),
        # actual text of tick marks
        axis.text.x = element_text(size = 21, face = "bold")
        # axis.text.x = element_text(size = 3, angle = 45
        # colour = COLORS[x_axis$key]
        # )
    )
}




#####



p <- ggplot() +
    geom_point(
        data = nifHDB_lng_Clst %>%
            mutate(
                CyanoGroupsIII = if_else(grepl("Tricho", CON), "Trichodesmium spp.", CyanoGroups),
                Clst3Grps = if_else(grepl("Cluster_III", CON), "Cluster_III_L", false = CON)
            ),
        aes(
            # y =RA,
            y = PO4_tblPisces_NRT,
            x = sst_tblSST_AVHRR_OI_NRT,
            group = CON,
            colour = logFe,
            size = RA
        ),
        na.rm = T,
        shape = 1,
        # size=1,
        stroke = 1
    ) +
    geom_text(
        data = plt_rgns,
        aes(
            y = PO4s,
            x = SSTs,
            label = regions
        ),
        size = 10,
        fontface = "bold",
        vjust = -1
    ) +
    labs(
        y = expression(bold(paste(
            "PO"[4]^"3-" * " (µmol L"^"-1", ")"
        ))),
        x = expression(bold(paste("SST (˚C)"))),
        # y = expression(bold(paste('', h^-1))),
        title = "",
        shape = "hemisphere",
        fill = "Taxa"
    ) +
    # guides(#colour = guide_legend(nrow = 2, byrow = TRUE),) +
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
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 16, face = "bold"),
        # titles of each facet
        panel.border = element_rect(colour = "black"),
        # axis.ticks.x = element_blank(),
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
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 17, face = "bold"),
        legend.position = "none",
        legend.box = "horizontal",
        legend.direction = "horizontal",
        # legend.position = 'none',
        axis.text.y = element_text(size = 21, face = "bold"),
        # actual text of tick marks
        axis.text.x = element_text(size = 21, face = "bold")
        # axis.text.x = element_text(size = 3, angle = 45
        # colour = COLORS[x_axis$key]
        # )
    )
print(p)


p +
    facet_wrap(~nifH_cluster,
        # ncol = 2,
        nrow = 2,
        # labeller=labeller(RA = studyID
        #                          ),
        scales = "free_y"
    )

p +
    facet_wrap(~group1,
        # ncol = 2,
        # nrow = 2,
        # labeller=labeller(RA = studyID
        #                          ),
        scales = "fixed"
    )

p +
    facet_wrap(~CON,
        # ncol = 2,
        # nrow = 2,
        # labeller=labeller(RA = studyID
        #                          ),
        scales = "fixed"
    )


p +
    facet_wrap(~studyID,
        # ncol = 2,
        # nrow = 1,
        # labeller=labeller(RA = studyID
        #                          ),
        scales = "fixed"
    )



p +
    facet_wrap(~Clst3Grps,
        # ncol = 2,
        # nrow = 2,
        # labeller=labeller(RA = studyID
        #                          ),
        scales = "fixed"
    )



p +
    facet_wrap(~CON,
        # ncol = 2,
        # ncol = 2,
        # labeller = ,
        # labeller=labeller(RA = studyID
        #                          ),
        scales = "fixed"
    )

p +
    facet_wrap(~CyanoGroupsIII,
        # ncol = 2,
        # nrow = 3,
        # labeller=labeller(RA = studyID
        #                          ),
        scales = "fixed"
    )
# p +
#   facet_wrap(~group1,
#              # ncol = 2,
#              nrow = 3,
#              # labeller=labeller(RA = studyID
#              #                          ),
#              scales = 'free_y'
#   )

p +
    facet_wrap(~group1,
        # ncol = 2,
        # nrow = 2,
        # labeller=labeller(RA = studyID
        #                          ),
        scales = "free_y"
    )

p +
    facet_wrap(~group2,
        # ncol = 2,
        # nrow = 3,
        # labeller=labeller(RA = studyID
        #                          ),
        scales = "free_y"
    )

p +
    facet_wrap(~group4,
        # ncol = 2,
        # nrow = 3,
        # labeller=labeller(RA = studyID
        #                          ),
        scales = "free_y"
    )

p +
    facet_wrap(~AUID,
        # ncol = 2,
        # nrow = 3,
        # labeller=labeller(RA = studyID
        #                          ),
        scales = "fixed"
    )


### updated nifH dataDB






## write out plots

ggsave()

















###### old stuff from ealier iterations


# p= ggplot(
#
# ) +
#   geom_point(
#     data = nifHDB_lng_Clst %>%
#       # mutate(
#       #   logFe = log(Fe_tblPisces_NRT)
#       # ) %>%
#       # filter(AUID %in% tempFilt) %>%
#       # %>%
#       # filter(!CON=='oligo6_A1_ATCTCGCTTTTTT') %>%
#       # filter(RA>0) %>%
#       # filter(logFe>=-9) %>%
#       # filter(CON %in% nonPCRclades) %>%
#       # filter(str_detect(pattern = 'gamma A', string = CON, negate = F)) %>%      ### this block here
#       # filter(nifH_cluster %in% c('1G')) %>%                                      ### is for gamma A specific filtering
#     # filter(!grepl('unknown', CON)) %>%                                         ### I did for ASLO 2023
#     # filter(AUID %in% tempFilt) %>%                                             ### #####
#     # %>%
#     # filter(!CON=='oligo6_A1_ATCTCGCTTTTTT') %>%
#     # filter(CON %in% c('Trichodesmium erythraeum', 'gamma 4', 'Pseudomonas stutzeri', 'UCYN-A1', 'gamma A')) %>%
#     # filter(grepl('Hall', studyID)) %>%
#     # filter(!CON %in% c('ALV82265.1' , 'CAL79071.1'))
#     # %>%
#     mutate(
#       CyanoGroupsIII = if_else( grepl('Tricho', CON) , 'Trichodesmium spp.', CyanoGroups),
#       Clst3Grps = if_else(grepl('Cluster_III', CON),  'Cluster_III_L', false =  CON)
#     )
#     , aes( #y =RA,
#       y = PO4_tblPisces_NRT,
#       # y = NO3_tblPisces_NRT,
#       # y = NP_Pisces_NRT,
#       # y = log(Fe_tblPisces_NRT),
#       # y = log(POFe_darwin),
#       # y =log(SiO2_darwin),
#       # x = lat_abs,
#       x = sst_tblSST_AVHRR_OI_NRT,
#       group = CON,
#       colour = logFe,
#       # colour = nifH_cluster,
#       # group = nifH_cluster,
#       # shape = AUID,
#       size = RA
#     ),
#     na.rm = T,
#     shape= 1, #size=1,
#     stroke = 1) +
#   geom_text(
#     data = data_GammaAfct,
#     aes(
#       # y = yy,
#       #   x = xx,
#       # label = labels,
#       y = PO4s,
#       x = SSTs,
#       label = regions), size = 10,fontface = "bold", vjust = -1) +
#   # stat_summary(aes(y = RA,
#   #                  # x = lat,
#   #                  x = lat_abs,
#   #                  # colour = genome879_nifh_c,
#   #                  colour = genome879_nifh_c,
#   #                  # colour = top30,
#   #                  group = genome879_nifh_c,
#   #                  # shape=df_type
#   # ),fun.data= mean_cl_normal) +
#   # geom_smooth(data = testGamma, aes( #y =RA,
#   #     y = sst_tblSST_AVHRR_OI_NRT / 30,
# #     # x = lat,
# #     x = lat_abs,
# #     # x = sst_tblSST_AVHRR_OI_NRT,
# #     # x = wind_stress_tblWind_NRT,
# #     # colour = CyanoGammaCON,
# #     # group = CyanoGammaCON,
# #     # colour = CyanoCON,
# #     # group = CyanoCON,
# #     # colour = CON,
# #     # group = CON,
# #     # colour = nifH_cluster,
# #     # group = nifH_cluster,
# #     shape = hemi
# # ), method='loess', formula= y~x, se = T, show.legend = F) +
# # geom_smooth(data = testGamma, aes( #y =RA,
# #   y = PO4_tblPisces_NRT / 2,
# #   # x = lat,
# #   x = lat_abs,
# #   # x = sst_tblSST_AVHRR_OI_NRT,
# #   # x = wind_stress_tblWind_NRT,
# #   # colour = CyanoGammaCON,
# #   # group = CyanoGammaCON,
# #   # colour = CyanoCON,
# #   # group = CyanoCON,
# #   # colour = CON,
# #   # group = CON,
# #   # colour = nifH_cluster,
# #   # group = nifH_cluster,
# #   shape = hemi
# # ), method='loess', formula= y~x, se = T, show.legend = F) +
# # stat_poly_line(formula = y ~ poly(x, 2, raw = TRUE)) +
# # stat_poly_line(aes(group = nifH_cluster)) +
# # aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
# # label.x.npc = "right", label.y.npc = 0.15,
# # formula = formula, parse = TRUE, size = 3) +
# # stat_poly_eq(formula = y ~ poly(x, 2, raw = TRUE),
# #              aes(label = after_stat(rr.label)),
# #              label.x.npc = "right",
# #              label.y.npc = "top",) +
# # stat_poly_eq(formula = y ~ x,
# #              aes(label = after_stat(rr.label)),
# #              label.x.npc = "left",
# #              label.y.npc = "top") +
# # geom_bar(aes(y = RA, x = Latitude,
# #              fill = phylo,
# #              group = hemi
# # ), na.rm = T,
# # position = 'fill',
# # stat = 'identity',
# # width = 0.2
# # ) +
# labs(
#   # x = expression(bold(paste('absolute latitude'))),
#   # y = expression(bold(paste( 'relative abundance/PO4' ))),
#   y = expression(bold(paste("PO"[4]^"3-"*" (µmol L"^"-1", ")" ))),
#   x = expression(bold(paste( 'SST (˚C)' ))),
#   #y = expression(bold(paste('', h^-1))),
#   title = '',
#   # title = 'ASV relative abundance of Top 30 Unknown ASVs all depths',
#   # title = 'Alpha/Beta (1J/1K) photic depths',
#   # title = 'Cyano clades photic depths',
#   # title = 'Major nifH clusters from photic depths',
#   # title = 'Major nifH clades from photic depths',
#   # title = 'Major nifH clades from photic depths',
#   colour = expression(bold(paste('Fe log'))),
#   shape = 'hemisphere',
#   fill = 'Taxa'
# ) +
#   guides(#colour = guide_legend(nrow = 2, byrow = TRUE),
#     # shape = guide_legend(nrow = 2, byrow = TRUE),
#     # colour = 'none',
#     # size = 'none'
#   ) +
#   # guides(colour = guide_legend(nrow = 3, byrow = TRUE),
#   #        shape = 'none') +
#   # scale_size_binned(n.breaks = 10) +
#   # scale_size_continuous( breaks = c( 0.1,0.3, 0.6, 1.00), range = c(.1,10), limits = c(0.01,1))+
#   scale_size_continuous(breaks = c(0,0.1,0.33,0.66, 1.00), range = c(.1,10), limits = c(0,1)) +
#   # scale_shape_binned() +
#   # scale_fill_fermenter(type = 'div',palette =  "RdBu", direction = -1)+
#   # scale_color_distiller(type = 'div',palette =  "Spectral", direction = -1 ) +
#   scale_color_fermenter(type = 'div',palette =  "Spectral", direction = -1, breaks = c(-11,-9,-7))+
#   # scale_colour_stepsn(breaks = c(0, 23, 35, 66), colours = getPalette(colourCount) #, labels = c(  'Eq/Trop',  'Sub-Trop',   'Temp',  'Poles')
#   #                     ) +
#   # scale_shape_binned(breaks = c(0, 23, 35, 66) #, colours = getPalette(colourCount) #, labels = c(  'Eq/Trop',  'Sub-Trop',   'Temp',  'Poles')
#   # ) +
#   # scale_colour_stepsn(breaks = c(10, 40, 100), colours = getPalette(colourCount)) +
#   # scale_color_fermenter(type = 'div',palette =  "RdYlBu", direction = -1)+
#   # scale_color_fermenter(type = 'div',palette =  "RdBu", direction = -1)+
#   # scale_color_gradientn(colours = rainbow(5)) +
#   # scale_color_manual(values = UCYN_Colors) +
#   # scale_colour_manual( values = c('UCYN-A1' = 'chartreuse4',
#   #                                 'UCYN-A2' = 'deeppink',
# #                                 'UCYN-A3' = 'blue',
# #                                 'UCYN-A4' = 'brown1',
# #                                 'Trichodesmium erythraeum' = 'darkorchid1',
# #                                 'unknown1B' = 'grey') )+
# # labels = c('oligo1_A1_ATCTCGCTTCTTT' = 'UCYN-A1',
# #            'oligo3_A2_ATTCTATTTTCTT' = 'UCYN-A2',
# #            'oligo2_A3_GCTCTATCCTTCA' = 'UCYN-A3',
# #            'oligo4_A4_ACTCTATTTCCCT' = 'UCYN-A4') )+ ## this lets you change the legend labels
# # scale_fill_manual(values = UCYN_Colors) + ## this lets you change the legend labels
# # scale_fill_manual(values = getPalette(colourCount)) + ## this lets you change the legend labels
# # scale_color_manual(values = c("1A" = "purple",
# #                               "1J/1K"="orange",
# #                               "1G"="steelblue", c('UCYN-A1' = 'green')
# #                               "1B" = "green2")) +
# # scale_fill_manual(values = c("1A" = "purple",
# #                               "1J/1K"="orange",
# #                               "1G"="steelblue",
# #                               "1B" = "green2")) +
# # scale_color_manual(values = c("1A" = "black",
# #                               "1J/1K"="steelblue",
# #                               "3" = "goldenrod2",
# #                               "1G"="purple",
# #                               "1B" = "green2")) +
# # scale_fill_manual(values = c("1A" = "black",
# #                              "1J/1K"="chocolate4",
# #                              "3" = "goldenrod2",
# #                              "1G"="steelblue",
# #                              "1B" = "green2")) +
# # scale_shape_manual(values = c('northernHemi'= 1,
# #                               'southernHemi' = 2),
# #                    labels = c('northernHemi'= 'North',
# #                               'southernHemi' = 'South')) +
# # scale_y_continuous(sec.axis = sec_axis(trans = ~. *30)) +
# # scale_fill_manual(values = getPalette(colourCount), labels =toPlot$phylo) + ## this lets you change the legend labels
# # scale_fill_manual(values = cbp2,  labels = p_labs) +
# # scale_shape_discrete(labels = c('Tricho', 'UCYN-A', 'Croco')) +
# # scale_color_discrete(labels = c('Tricho', 'UCYN-A', 'Croco')) +
# ylim(-0.05,1.95)+
#   xlim(-2,30.4)+
#   theme_bw() +
#   theme(plot.title = element_text(face = 'bold'),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_blank(),
#         strip.text.x = element_text(size = 16, face = 'bold'), # titles of each facet
#         panel.border = element_rect(colour = "black"),
#         #axis.ticks.x = element_blank(),
#         axis.title = element_text(size = 15, face = 'bold'),
#         axis.title.x = element_text(size = 23, face = 'bold', colour = 'black'),
#         axis.title.y = element_text(size = 21, face = 'bold', colour = 'black'),
#         legend.title = element_text(size = 13, face = 'bold'),
#         legend.text = element_text(size = 17, face = 'bold'),
#         legend.position="none",
#         legend.box = "horizontal",
#         legend.direction = "horizontal",
#         # legend.position = 'none',
#         axis.text.y = element_text(size = 21, face = 'bold'), # actual text of tick marks
#         axis.text.x = element_text(size = 21, face = 'bold')
#         # axis.text.x = element_text(size = 3, angle = 45
#         # colour = COLORS[x_axis$key]
#         # )
#
#   );p
#
# p +
#   facet_wrap(~CON,
#              # ncol = 2,
#              # nrow = 2,
#              # labeller=labeller(RA = studyID
#              #                          ),
#              scales = 'fixed'
#   )
#
# p +
#   facet_wrap(~CON,
#              # ncol = 2,
#              # nrow = 2,
#              # labeller=labeller(RA = studyID
#              #                          ),
#              scales = 'fixed'
#   )
#
# p +
#   facet_wrap(~studyID,
#              # ncol = 2,
#              # nrow = 1,
#              # labeller=labeller(RA = studyID
#              #                          ),
#              scales = 'fixed'
#   )
#
#
#
# p +
#   facet_wrap(~Clst3Grps,
#              # ncol = 2,
#              # nrow = 2,
#              # labeller=labeller(RA = studyID
#              #                          ),
#              scales = 'fixed'
#   )
#
# p +
#   facet_wrap(~nifH_cluster,
#              # ncol = 2,
#              nrow = 2,
#              # labeller=labeller(RA = studyID
#              #                          ),
#              scales = 'free_y'
#   )
#
# p +
#   facet_wrap(~CON,
#              # ncol = 2,
#              # ncol = 2,
#              # labeller = ,
#              # labeller=labeller(RA = studyID
#              #                          ),
#
#              scales = 'fixed'
#   )
#
# p +
#   facet_wrap(~CyanoGroupsIII,
#              # ncol = 2,
#              # nrow = 3,
#              # labeller=labeller(RA = studyID
#              #                          ),
#              scales = 'fixed'
#   )
# # p +
# #   facet_wrap(~group1,
# #              # ncol = 2,
# #              nrow = 3,
# #              # labeller=labeller(RA = studyID
# #              #                          ),
# #              scales = 'free_y'
# #   )
#
# p +
#   facet_wrap(~group2,
#              # ncol = 2,
#              nrow = 2,
#              # labeller=labeller(RA = studyID
#              #                          ),
#              scales = 'free_y'
#   )
#
# p +
#   facet_wrap(~group2,
#              # ncol = 2,
#              nrow = 3,
#              # labeller=labeller(RA = studyID
#              #                          ),
#              scales = 'free_y'
#   )
#
# p +
#   facet_wrap(~group3,
#              # ncol = 2,
#              nrow = 3,
#              # labeller=labeller(RA = studyID
#              #                          ),
#              scales = 'free_y'
#   )
#
# p +
#   facet_wrap(~AUID,
#              # ncol = 2,
#              # nrow = 3,
#              # labeller=labeller(RA = studyID
#              #                          ),
#              scales = 'fixed'
#   )
