## To generate the top 20 rows of Table 4.
tab4 <- read.csv('Table_SreadsAtEachStage_samples.csv')
aggregate( . ~ Study, tab4[,-2], function(v) mean(v[v>0]))

## Calc'd by my script that makes Fig. 4.  Not sure if my summary rows (last 3 rows)
## are calc'd as in the manuscript.
tab5 <- read.csv('table_5_for_manuscript.csv')

