What i have left to do

### plots for kenda
#Bar chart (or the like): describing # samples per depth
Bar chart (or the like): describing # samples from each study from cDNA vs. DNA
**This can be part of a figure that show this and other information similar to this.
** As percentages of total counts
# two different bar plots, side by side 
-nucleic acids
-photic
# one bar plot 
-duplicate samples - 155
-single samples - 578
-size fractions - 138
Bar chart (or the like): describing relative abundances of nifH clusters from each study (total dataset)

## NEW PLOT ###
month vs season 
-on same axis
-month is veritcal axis

### tables
Table 1
Table 2 
Table 3
Table 4



### look into this
## summarise the counts and relataive abundace data to see there importances
## may need to also add CMAP file so regions can be explored as well.
annoNifHDB_updt %>%
  filter(nifH_cluster == "1P") %>%
  distinct(CON)
