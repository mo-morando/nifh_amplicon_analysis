library(tidyverse)

# ### i think this is what comes 
# auid.abundancesRA = read_tsv('/Users/mo/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/master_ASVht/FilterAuids/auid.abundances.filteredRA.tsv') %>% 
#   rename(AUID = ...1) %>% 
#   column_to_rownames('AUID')
# 
# dim(auid.abundancesRA)
# rm(auid.abundancesRA)

### get a data time stamp

# Get the current date
current_date <- Sys.Date()

# Format the date string
date_stamp <- format(current_date, "%Y%m%d")

# # Add the date stamp to a string
# original_string <- "Example string"
# string_with_date <- paste0(original_string, "_", date_stamp)

# # Print the result
# print(string_with_date)
cat("Today's date is",print(date_stamp))

paste('all_studies/master_fasta/nifHcatalog/nifHcatalog_CntPer',date_stamp, '.csv', sep = '')



#### read in data


fastaFile_DB = Biostrings::readDNAStringSet("all_studies/master_fasta/nifHDB/auid.filtered.nifHDB.fasta")
seq_name = names(fastaFile_DB)
sequence = paste(fastaFile_DB)
nifHDB.fa_df <- data.frame(seq_name, sequence)


dfExpnd.fa = nifHDB.fa_df %>% 
  mutate (
    AUID = str_extract('AUID\\.\\d+',string = nifHDB.fa_df$seq_name),
    studies = str_remove('AUID\\.\\d+',string = nifHDB.fa_df$seq_name)
  )


dfExpnd.fa_nifHDB = dfExpnd.fa %>%
#   filter(AUID %in% TargetAUIDs ) %>%
  select(AUID, sequence)


### convert into actual fasta file
seqs =  Biostrings::DNAStringSet(dfExpnd.fa_nifHDB$sequence)

# Assign AUID values as names to the seqs vector
names(seqs) <- dfExpnd.fa_nifHDB$AUID

Biostrings::writeXStringSet(x = seqs,  
                            filepath = paste('~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/master_fasta/nifHDB/nifHDB_sequences_CFH',date_stamp, '.fa', sep = ''), 
                            format = 'fasta')


# dir.create('all_studies/master_fasta/nifHcatalog/')

# write_csv(A2KHI.df, 'all_studies/master_fasta/A2KHI_AUID.csv')
write_csv(dfExpnd_nifHDB, paste('all_studies/master_fasta/nifHDB/nifHDB_sequences',date_stamp, '.csv', sep = ''))

write_csv(dfExpnd_nifHDB, paste('all_studies/master_fasta/nifHDB/nifHDB_sequences',date_stamp, '.csv', sep = ''))













### abundance table
nifhDB_cnts =read.table("all_studies/master_fasta/nifHDB/auid.abundances.filtered.nifHDB.tsv", header = TRUE, sep = "\t", row.names = 1) %>% 
  rownames_to_column('AUID') %>%
  rename_all(~str_remove(., pattern = '.*___'))

dim(nifhDB_cnts)
names(nifhDB_cnts)

### remove samples from NEMO that I have no idea what they are 

nifhDB_cnts = nifhDB_cnts %>% 
  select(-matches('Turk\\d+\\.e')) %>%
    select(-matches("Harding229\\.66705_S229|Harding230\\.66706_S230|Harding231\\.66709_S231")) %>%
  view()


##### do some calcluations

nifhDB_RA = nifhDB_cnts %>% 
  # column_to_rownames('AUID') %>% 
  # group_by(SAMPLEID) %>%
  mutate(across(where(is.numeric), ~./sum(., na.rm = T)))
# mutate(relative_abundance = count / sum(count)) %>%
# pivot_wider(names_from = ASV, values_from = relative_abundance)

nifhDB_RA %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  # filter(.!=1) %>%  
  view()


# Assuming your dataframe is named 'data'
nifhDB_cnts_T <- nifhDB_cnts %>%
  pivot_longer(cols = -AUID, names_to = "SAMPLEID", values_to = "Value") %>%
  pivot_wider(names_from = AUID, values_from = Value) #%>% 
  # mutate(
  #   SAMPLEID = str_remove(pattern = '.*___', SAMPLEID)
  # )

nifhDB_RA_T <- nifhDB_RA %>%
  pivot_longer(cols = -AUID, names_to = "SAMPLEID", values_to = "Value") %>%
  pivot_wider(names_from = AUID, values_from = Value) #%>% 
  # mutate(
  #   SAMPLEID = str_remove(pattern = '.*___', SAMPLEID)
  # )


#### read in annotations table so you can make catalog
annotations

annotationsNifHDB = annotations %>% 
  filter(AUID %in% nifhDB_cnts$AUID) %>% 
  left_join(phyloAll, by = 'AUID')%>% 
  mutate(
    CON = if_else(gamma_cluster=='unk' | is.na(gamma_cluster), ConsClassKTK, gamma_cluster)
  )%>% 
  mutate(
    group1 = if_else(nifH_cluster %in% c('1A', '3'), '1A C3', nifH_cluster), 
    # group2 = if_else(nifH_cluster %in% c('1A', '3', '1J/1K'), '1A C3 1J/1K', nifH_cluster),
    group2 = if_else(nifH_cluster %in% c('1A', '3', '1J/1K'), '1A,  3,  &  1J/1K', nifH_cluster),
    group3 = if_else(nifH_cluster %in% c('1A', '3'), '1A C3',
                     if_else(nifH_cluster %in% c('1G', '1J/1K'), '1B 1J/1K', nifH_cluster)),
    group4 = if_else(nifH_cluster %in% c('1A', '3', '1J/1K'), '1A/C3 1J/1K', '1B & 1G'),
    CyanoCON = if_else(nifH_cluster=='1B', CON, nifH_cluster),
    JKalphaCON = if_else(nifH_cluster=='1J/1K' & ConsClassKTK %in% JKalphaSlct, ConsClassKTK, nifH_cluster),
    crocoCMB = if_else(grepl('croco', CON, ignore.case = T) & CON=='Crocosphera_DQ118216_Moisander', 'CDQmois', 
                       if_else(grepl('croco', CON, ignore.case = T), 'CrocoCMB', CON)),
    CyanoGroups = if_else( CON=='UCYN-A3' | CON=='UCYN-A4' , 'A3-A4', CON),
    CyanoGroupsII = if_else( grepl('Tricho', CON) , 'Trichodesmium sp.', CON),
    CyanoGroupsIII = if_else( grepl('Tricho', CON) , 'Trichodesmium sp.', CyanoGroups),
    CyanoGroupsIV = if_else( CON=='UCYN-A3' | CON=='UCYN-A1' , 'A1-A3', 
                             if_else(CON=='UCYN-A2' | CON=='UCYN-A4', 'A2-A4',CON ))
    # CyanoGammaCON = if_else(nifH_cluster=='1B' | nifH_cluster=='1G', CON, nifH_cluster),
    # CyanoGammaCON = if_else(CON %in% CyanosToRemove$CON | CON %in% GammasToRemove, nifH_cluster, CyanoGammaCON)
  ) 


annotationsNifHDB %>% 
  filter(is.na(CON)) %>% 
  view()

#### bring in the metadata
# 
#### these files were before we keep all the data like RNA so they are not helpful with nifhDB
# CMAP_full = read_tsv('all_studies/merged_metadata/CMAP/GatherMetadata/metadata.cmap.tsv') %>% 
#   rename(studyID = StudyID)
# 
# CMAP_full = read_csv('all_studies/merged_metadata/CMAP/20220831CMAP_colocalized_nifH.csv') %>% 
#   # mutate(
#   #   StudyID = if_else(StudyID=='Gradoville_2020' & grepl('2017', DateTimeUTC), 'Gradoville_2020_G2', 
#   #                     if_else(StudyID=='Gradoville_2020' & grepl('2016', DateTimeUTC), 'Gradoville_2020_G1', StudyID))
#   # ) %>% 
#   rename(studyID = StudyID)


CMAP_full = read_tsv('~/mmorando@ucsc.edu_gDrive/My Drive/data/amplicon_review/all_studies/merged_metadata/nifHDB/firstAttempt/MO_here_are_tables_with_also_transcriptomic/metadata.cmap.tsv') %>% 
    rename(studyID = StudyID,
           lat = Lat,
           lon = Lon,
           depth = Depth) %>% 
  mutate(
    transcipt_flag = ifelse(grepl("_transcriptomic", SAMPLEID), T, F),
    nucleicAcidType = if_else(transcipt_flag==T, 'RNA', 'DNA'),
    SAMPLEID = str_remove(SAMPLEID, "_transcriptomic"),
    studyID = if_else(studyID=='Gradoville_2020' & grepl('2017', DateTimeUTC), 'Gradoville_2020_G2', 
                                            if_else(studyID=='Gradoville_2020' & grepl('2016', DateTimeUTC), 'Gradoville_2020_G1', studyID)
    ),
    time = str_replace(DateTimeUTC, ' ', 'T'),
    time = str_replace(time, "T.+$", "T12:00:00"),
      depth = str_replace(depth, 'surface', '3'),
      depth = str_replace(depth, '^0$', '3'),
    depth = str_replace(depth, '<2\\\\, depending on the tide', '2'),
    depth = str_replace(depth, '~2', '2')
    # depth = as.numeric(depth)
  ) %>%
     select(SAMPLEID, transcipt_flag, nucleicAcidType, time,  everything())

CMAP_full %>% 
  # filter(studyID=='Harding_2018') %>% 
  # filter(is.na(DateTimeUTC)) %>% 
  filter(is.na(lat)) %>% 
  view()

raes_id = CMAP_20230719  %>% 
filter(grepl("Raes", studyID)) %>% 
pull(SAMPLEID)

CMAP_full  %>% 
# distinct(studyID)  %>% print(n=50)
filter(!studyID %in% "Raes_2020")


filter(SAMPLEID %in% raes_id)

# 
# CMAP_full%>% 
#   # filter(studyID=='BentzonTilia_2015') %>% 
#   left_join(CMAPtempII %>% 
#               rename(studyID = StudyID) ) %>% 
#   view()


# CMAPtempII %>% 
#   rename(studyID = StudyID) %>% 
#   select(SAMPLEID, studyID, DateTimeUTC)

### identify the columns i want
CMAP_cols = c("SAMPLEID"    ,     "transcipt_flag"  , "studyID" ,"LibrarySelection", "LibrarySource"  ,"Collection_Date",  "time"  ,"depth"    ,        "Station"  ,  "Size_fraction", "lat", "lon" , "DateTimeUTC" )
CMAP_full %>% 
  select(CMAP_cols)

  
CMAP_20230719 = read_csv('all_studies/merged_metadata/nifHDB/firstAttempt/20230719_colocalized_nifH.csv')

CMAP_coloc = CMAP_20230719


# dim(CMAP_full_old)
dim(CMAP_full)
names(CMAP_full)
dim(CMAP_coloc)
names(CMAP_coloc)

CMAP_full %>% 
  distinct(time)

# all.equal(CMAP_full$SAMPLEID, CMAP_full_old$SAMPLEID)


### this a parsing problem with this file right now so we have to remove some bad columns

# CMAP_full = CMAP_full %>% 
#   # filter(grepl(x = studyID, pattern =("Harding229\\.66705_S229|Harding230\\.66706_S230|Harding231\\.66709_S231"))) %>% 
#   filter((studyID=='BentzonTilia_2015' & is.na(DateTimeUTC)))  



CMAP_pnts = CMAP_coloc %>% 
  distinct(
    # depth, 
    lat, lon, time,studyID) %>% mutate(
      # time = str_extract(pattern = "\\d{4}-\\d{2}-\\d{2}", time),
      time = format(time, "%Y-%m-%d")
    ) %>% view()



### we need distance from shore to make coast and photic columns
### we need distance from shore to make coast and photic columns
### we need distance from shore to make coast and photic columns
### we need distance from shore to make coast and photic columns

read_csv("~/mmorando@ucsc.edu_gDrive/My Drive/data/amplicon_review/all_studies/Thal_RwrkSpc/smetaTab.csv")

smetaTab_filt = smetaTab %>% 
  filter(!grepl(x = SAMPLEID, pattern ='Turk\\d+\\.e' )) %>%
  filter(!grepl(x = SAMPLEID,pattern ="Harding229\\.66705_S229|Harding230\\.66706_S230|Harding231\\.66709_S231")) %>% 
  rename(
    depth = Depth,
    lat = Lat,
    lon = Lon,
    time = DateTimeUTC,
    studyID = StudyID) %>% 
  mutate(
    depth = str_replace(depth, 'surface', '3'),
    depth = str_replace(depth, '^0$', '3'),
    depth = replace_na(depth, '5'),
    depth = as.numeric(depth),
    studyID = if_else(studyID=='Gradoville_2020' & grepl('2017', time), 'Gradoville_2020_G2', 
                      if_else(studyID=='Gradoville_2020' & grepl('2016', time), 'Gradoville_2020_G1', studyID))
    # photic = if_else(condition = kmToCoast>=10 & depth<=100, true =  T, 
    #                  false = if_else(kmToCoast<=10 & kmToCoast>=1 & depth<=50, T,
    #                                  if_else(kmToCoast<1 & depth<=20, T, F ))
    # )
  ) %>% 
  select(everything()
  ) %>%
  select(
    # SAMPLEID,
    # depth, 
    lat, lon, time,studyID, kmToCoast, nearestLand.Lon,  nearestLand.Lat) %>%
  distinct() %>%
  mutate(
    # time = str_extract(pattern = "\\d{4}-\\d{2}-\\d{2}", time),
    time = format(time, "%Y-%m-%d")
  )


### replace missing values
bentz <- list(
  studyID = 'BentzonTilia_2015',
  # photic = T,
  kmToCoast = 0,
  nearestLand.Lon = 12.0861,
  nearestLand.Lat = 55.6985
)

Mulh <- list(
  studyID = 'Mulholland_2018',
  # photic = T,
  kmToCoast = 188,
  nearestLand.Lon = -69.9247,
  nearestLand.Lat = 41.8381
)

shi2017 <- list(
  studyID = 'Shiozaki_2017',
  # photic = T,
  kmToCoast = 1123,
  nearestLand.Lon = -167.9999,
  nearestLand.Lat = 25.0164
)

# Combine the lists into a single data frame
replacement_df <- bind_rows(bentz, Mulh, shi2017)

# Left join 'temp' with 'replacement_df' based on 'studyID', then use coalesce to fill NAs
library(lubridate)

CMAP_pnts_smeta = CMAP_pnts %>% 
  left_join(
    smetaTab_filt
  ) %>%
  # filter(is.na(kmToCoast)) %>%
  left_join(replacement_df, by = "studyID") %>% 
  mutate(
    # photic = coalesce(photic.x, photic.y),
    kmToCoast = coalesce(kmToCoast.x, kmToCoast.y),
    nearestLand.Lon = coalesce(nearestLand.Lon.x, nearestLand.Lon.y),
    nearestLand.Lat = coalesce(nearestLand.Lat.x, nearestLand.Lat.y)
  ) %>%
  select(-ends_with(".x"), -ends_with(".y")) %>% 
  mutate(
    time = paste(time, "12:00:00 UTC"),
    time = lubridate::as_datetime(time, tz = "UTC")
  ) %>%  
  # filter(is.na(kmToCoast)) %>%
  view()

CMAP_pnts_smeta %>% 
  left_join(CMAP_coloc %>% 
              distinct(
                # depth, 
                SAMPLEID, lat, lon, time,studyID) %>% mutate(
                  # time = str_extract(pattern = "\\d{4}-\\d{2}-\\d{2}", time),
                  time = format(time, "%Y-%m-%d"),
                  time = paste(time, "12:00:00 UTC"),
                  time = as_datetime(time, tz = "UTC")
                )) %>%  left_join(replacement_df, by = "studyID") %>%
  mutate(
    kmToCoast = coalesce(kmToCoast.x, kmToCoast.y),
    nearestLand.Lon = coalesce(nearestLand.Lon.x, nearestLand.Lon.y),
    nearestLand.Lat = coalesce(nearestLand.Lat.x, nearestLand.Lat.y)
  ) %>%
  select(-ends_with(".x"), -ends_with(".y")) %>% view()


CMAP_coloc_WF = CMAP_coloc %>%
  left_join(CMAP_pnts_smeta) %>% 
  
  mutate(#lat = as.numeric(Lat),
    photic = if_else(condition = kmToCoast>=10 & depth<=100, true =  T, 
                     false = if_else(kmToCoast<=10 & kmToCoast>=1 & depth<=50, T,
                                     if_else(kmToCoast<1 & depth<=20, T, F ))
    ),
    hemi = cut(.$lat,
                    breaks=c(-Inf, 0, Inf),
                    labels=c("southernHemi","northernHemi")),
         hemi = factor(hemi, c("northernHemi","southernHemi")),
         lat_abs = abs(lat),
         date = str_remove(time, ' .+$'),
         season = str_remove_all(str_extract(date, '-.+-'), '-'),
         season = if_else(season %in% c('03','04','05') & hemi=='northernHemi', 'spring',
                          if_else(season %in% c('09','10','11') & hemi=='northernHemi', 'fall',
                                  if_else(season %in% c('06','07','08') & hemi=='northernHemi', 'summer',
                                          if_else(season %in% c('12','01','02') & hemi=='northernHemi', 'winter',
                                                  if_else(season %in% c('03','04','05') & hemi=='southernHemi', 'fall',
                                                          if_else(season %in% c('09','10','11') & hemi=='southernHemi', 'spring',
                                                                  if_else(season %in% c('06','07','08') & hemi=='southernHemi', 'winter',
                                                                          if_else(season %in% c('12','01','02') & hemi=='southernHemi', 'summer', 'XXX')
                                                                  ))))))),
         season = factor( season, c('winter', 'spring','summer', 'fall')),
         CMAP_NP_darwin = CMAP_NO3_darwin_clim_tblDarwin_Nutrient_Climatology / CMAP_PO4_darwin_clim_tblDarwin_Nutrient_Climatology,
         CMAP_NP_Pisces_NRT = CMAP_NO3_tblPisces_NRT / CMAP_PO4_tblPisces_NRT,
         CMAP_NP_WOA_clim = CMAP_nitrate_WOA_clim_tblWOA_Climatology / CMAP_phosphate_WOA_clim_tblWOA_Climatology,
         CMAP_NP_WOA_clim = CMAP_n_an_clim_tblWOA_2018_1deg_Climatology / CMAP_p_an_clim_tblWOA_2018_1deg_Climatology
  ) %>%
  # select(-c(time, Region, geo_loc_name, Size_fraction, nearestLand.Lon, nearestLand.Lat, lat)) %>%
  # rename_with(~str_remove(string = ., pattern = 'CMAP')) %>%
  rename_with(~str_remove(string = ., pattern = 'CMAP_')) %>%
  # rename_with(~str_extract(string = ., pattern = '^.+darwin'), contains('darwin')) %>%
  rename_with(~coalesce(str_extract(string = ., pattern = "^[A-Za-z0-9]+_darwin"), .)) %>%
  rename_with(~str_replace(string = ., pattern = 'clim_tblWOA_2018_qrtdeg_Climatology',replacement = 'WOA_2018_qrtdeg')) %>%
  rename_with(~str_replace(string = ., pattern = 'clim_tblWOA_2018_1deg_Climatology',replacement = 'WOA_2018_1deg')) %>%
  rename_with(~str_replace(string = ., pattern = 'WOA_clim_tblWOA_Climatology',replacement = 'WOA_clim')) %>%
  rename_with(~str_replace(string = ., pattern = 'clim_tblWOA_2018_MLD_qrtdeg_Climatology',replacement = 'WOA_2018_MLD_qrtdeg')) %>%
  rename_with(~str_replace(string = ., pattern = '^C_',replacement = 'conductivity_')) %>%
  rename_with(~str_replace(string = ., pattern = '^i_',replacement = 'sigma_')) %>%
  rename_with(~str_replace(string = ., pattern = '^t_',replacement = 'temp_')) %>%
  rename_with(~str_replace(string = ., pattern = '^A_',replacement = 'AOU_')) %>%
  rename_with(~str_replace(string = ., pattern = '^O_',replacement = 'O2_sat_')) %>%
  rename_with(~str_replace(string = ., pattern = '^n_',replacement = 'NO3_')) %>%
  rename_with(~str_replace(string = ., pattern = '^p_',replacement = 'phosphate_')) %>%
  rename_with(~str_replace(string = ., pattern = '^si_',replacement = 'silica_')) %>%
  rename_with(~str_replace(string = ., pattern = '^M_',replacement = 'MLD_')) %>%
  rename_with(~str_replace(string = ., pattern = '^s_',replacement = 'salinity_')) %>%
  mutate(
    geoRegion = cut(.$lat_abs,
                    breaks = c(-1, 23, 35, 66, 100),
                    labels = c(  'Eq/Trop',  'Sub-Trop',   'Temp',  'Poles')),
    geoRegion = factor(geoRegion, c(  'Eq/Trop',  'Sub-Trop',   'Temp',  'Poles')),
    logFe = log(Fe_tblPisces_NRT)
  )
#%>%
  # select(names(CMAPmeta))

names(CMAP_coloc_WF)
dim(CMAP_coloc_WF)



CMAP_coloc_WF %>% 
  filter(!is.na(kmToCoast)) %>% 
  view()












CMAP_coloc

abundTab =nifhDB_cnts

metaCols = c("lat", "lon", "time", "studyID")


metaTab = read.table("CMAP_coloc_20230810.cvs", header = T, sep = ",", colClasses = "character")
names(metaTab)

if (!all(metaCols %in% colnames(metaTab))) {
  cat(
    "The following columns are missing in the metadata file: ",
    setdiff(metaCols, colnames(metaTab)), "\n"
  )
  stop("Aborting because cannot partition samples using fields asked.")
}

apply(metaTab[, metaCols], 1, paste, collapse = ",")

SAMPLEID

meta2acol <- lapply(metaTab[, 1], grep, x = colnames(abundTab))
names(meta2acol) <- metaTab[, 1]


x <- apply(metaTab[, metaCols], 1, paste, collapse = ",")
grpInfo <- data.frame(metaTab[, metaCols], GroupNum = sapply(x, function(v) match(v, x)))
grpSz <- table(grpInfo$GroupNum)
grpInfo <- data.frame(grpInfo, GroupSize = as.vector(grpSz[as.character(grpInfo$GroupNum)]))

x <- sapply(meta2acol, length)
stopifnot(x <= 1)
meta2acol <- unlist(meta2acol[x == 1]) # Just the metadata needed for the abundance table.
grpInfo <- grpInfo[x == 1, ] # Similiary restrict grpInfo
## For each sample, add its matched abundanced column.
grpInfo$abundCol <- colnames(abundTab)[as.numeric(meta2acol)]
cat("done.\n")
stopifnot(nrow(grpInfo) <= ncol(abundTab))

setdiff(metaCols, colnames(metaTab))

setdiff(c("SAMPLEID, lat, lon, time, studyID"), colnames(metaTab))

setdiff("SAMPLEID", colnames(metaTab))
normalizeSampSizes <- F



## Helper funcs that used grpInfo and the abundTab
combineSampsWithNormalization <- function(gn) {
  ## Get in-group samples that are not 0 for every ASV.
  idx <- grpInfo[grpInfo$GroupNum == gn, "abundCol"]
  stopifnot(length(idx) > 0)
  gtab <- abundTab[, idx, drop = F]
  idx <- which(colSums(gtab) > 0)
  asvCounts <- rep(0, nrow(gtab)) # Assume no ASVs
  if (length(idx) >= 1) {
    gtab <- gtab[, idx, drop = F]
    m <- max(colSums(gtab))
    asvCounts <- rowSums(scale(gtab, center = F, scale = colSums(gtab) / m))
    asvCounts <- (asvCounts * sum(gtab) / ncol(gtab)) / m
  }
  ## Check that normalization did not change the total reads, allowing
  ## for tiny diff due to integer vs. float.
  tots <- c(sum(gtab), sum(asvCounts))
  stopifnot(abs(tots[1] - tots[2]) < 1E-9)
  asvCounts
}

combineSampsBySimplyAdding <- function(gn) {
  idx <- grpInfo[grpInfo$GroupNum == gn, "abundCol"]
  rowSums(abundTab[, idx, drop = F]) # okay if just 1 sample in the group
}


grpNums <- unique(grpInfo$GroupNum)
if (normalizeSampSizes) {
  cat(
    "Normalizing samples within each group to have same count, without changing\n",
    "the total reads in the group.\n"
  )
  newAbundTab <- sapply(grpNums, combineSampsWithNormalization)
} else {
  cat("Simply adding up the reads (ASV by ASV) across samples within each group.\n")
  newAbundTab <- sapply(grpNums, combineSampsBySimplyAdding)
  ## Two ways to check that each ASV's total is as before, w.r.t. samples with metadata
  stopifnot(rowSums(newAbundTab) == rowSums(abundTab[, grpInfo$abundCol]))
  stopifnot(rowSums(newAbundTab) == rowSums(abundTab[, -noMetaSampsIdx]))
}
stopifnot(nrow(newAbundTab) == nrow(abundTab))

idx <- match(grpNums, grpInfo$GroupNum)
stopifnot(!is.na(idx))
colnames(newAbundTab) <- grpInfo[idx, "abundCol"]
stopifnot(colnames(newAbundTab) %in% colnames(abundTab))

noMetaSampsIdx <- setdiff(1:ncol(abundTab), as.numeric(meta2acol))

## Tack on the metadata-less samples.
stopifnot(rownames(newAbundTab) == rownames(abundTab)) # cbind should check but...
x <- cbind(newAbundTab, AUID = abundTab[, noMetaSampsIdx])
newAbundTab <- x
length(newAbundTab)
nrow(newAbundTab)


ti <- sum(abundTab[-1])
tf <- sum(apply(newAbundTab[,-ncol(newAbundTab)], 2, as.numeric) )
cat(ti, " Total reads in input abundance table\n")
cat(tf, " Total reads after combining samples by group\n")
cat(round(100 * (tf - ti) / ti, 1), " % change\n")
rm(ti, tf)


names(newAbundTab)
str(newAbundTab)








nifhDB_cnts_T %>% 
  left_join(CMAP_coloc_WF)%>% 
  # left_join(CMAP_full_sub %>% 
  #             select(SAMPLEID, studyID))%>% 
  filter(is.na(studyID)) %>% 
  # filter(grepl('SRR1', SAMPLEID)) %>% 
  # mutate(
  #   SAMPLEID =  str_extract(pattern = '......', SAMPLEID)
  # ) %>% 
  distinct(SAMPLEID) %>% 
  # group_by(SAMPLEID) %>% 
  # count() %>% 
  # ungroup()
  view()

nifhDB_cnts_T_CMAP = nifhDB_cnts_T %>% 
  # left_join(CMAP_full %>% 
  #             select( any_of(CMAP_cols))) %>% 
  left_join(CMAP_coloc_WF) %>% 
  select( -contains('AUID'), everything(), contains('AUID')) %>% 
  mutate(
    depth = as.numeric(depth),
    depth = if_else(studyID=='Raes_2020', 5, depth)
  ) 

nifhDB_RA_T_CMAP = nifhDB_RA_T %>% 
  # left_join(CMAP_full %>% 
  #             select( any_of(CMAP_cols))) %>% 
  left_join(CMAP_coloc_WF) %>% 
  select( -contains('AUID'), everything(), contains('AUID')) %>% 
  mutate(
    depth = as.numeric(depth),
    depth = if_else(studyID=='Raes_2020', 5, depth)
  )

nifhDB_cnts_T_CMAP %>% 
  filter(studyID=='Hallstrom_2021') %>% 
  nrow()

# Raes_2020_SAMPLEID = nifhDB_cnts_T_CMAP %>% 
#   filter(studyID=='Raes_2020') %>% 
#   distinct(SAMPLEID) %>% 
#   pull()


# CMAP_full_sub = CMAP_full %>% 
#   # distinct(depth) %>% view()
#   filter(SAMPLEID %in% Raes_2020_SAMPLEID) %>%
#   mutate(
#     depth = if_else(studyID=='Raes_2020', '5', depth),
#     depth = as.numeric(depth),
#   ) 





# #### make new collection file to retrive CMAP data
# CMAP_20230718 = nifhDB_RA_T_CMAP %>% 
#   select(SAMPLEID, depth,  lat, lon , time,studyID, nucleicAcidType, Size_fraction, geo_loc_name ) %>% 
#   view()
# 
# 
# write_csv(CMAP_20230718, 'all_studies/merged_metadata/CMAP/20230718CMAP_colocalized_nifH.csv')
# write_csv(CMAP_20230718, '/Users/mo/mikemojr@gmail.com - Google Drive/My Drive/Colab Notebooks/Untitled Folder/20230718CMAP_colocalized_nifH.csv')








#### do some filtering 



nifhDB_cnts_T_CMAP_lng = nifhDB_cnts_T_CMAP %>%
  pivot_longer(cols = starts_with('AUID'), names_to = 'AUID', values_to = 'cnt') %>%
  filter(cnt>0)

nifhDB_RA_T_CMAP_lng = nifhDB_RA_T_CMAP %>%
  pivot_longer(cols = starts_with('AUID'), names_to = 'AUID', values_to = 'RA') %>%
  filter(RA>0)

# # Gather the data into a long format
# df_long <- nifhDB_RA %>%
#   pivot_longer(-AUID, names_to = "station", values_to = "RA") %>% 
#   mutate(
#     RA = replace_na(RA,0))

df = nifhDB_RA_T_CMAP
df_long = nifhDB_RA_T_CMAP_lng
df_long %>% 
  distinct(SAMPLEID) %>% 
  nrow()



####

df_long %>% names()


### nucleic acid type

# # Group by AUID and check if it has both DNA and RNA
# df_long <- df_long %>%
#   group_by(AUID) %>%
#   mutate(nucleicAcidTypeSummary = ifelse("DNA" %in% nucleicAcidType & "RNA" %in% nucleicAcidType, "BOTH",
#                                          ifelse("DNA" %in% nucleicAcidType, "DNA",
#                                                 ifelse("RNA" %in% nucleicAcidType, "RNA", NA_character_))))

# # Print the modified data frame
# NA_type = df_long %>% 
#         distinct(AUID, nucleicAcidTypeSummary) 

NA_type2 = nifhDB_cnts_T_CMAP_lng %>% 
  distinct(AUID, nucleicAcidType) %>%
  group_by(AUID) %>% 
  # pivot_wider(names_from = nucleicAcidType, values_from = AUID, values_fn = length, names_prefix = "has_") %>%
  mutate(nucleicAcidTypeSummary = case_when(
    sum(nucleicAcidType == "DNA") > 0 & sum(nucleicAcidType == "RNA") > 0 ~ "BOTH",
    sum(nucleicAcidType == "DNA") > 0 ~ "DNA",
    sum(nucleicAcidType == "RNA") > 0 ~ "RNA",
    TRUE ~ "NONE"
  )) %>%
  select(-c(nucleicAcidType)) %>% 
  distinct(AUID, .keep_all = T) %>% 
  view()

# # Compare the nucleicAcidTypeSummary values between the two data frames
# result <- inner_join(NA_type, NA_type2, by = "AUID") %>%
#   mutate(summary_equal = nucleicAcidTypeSummary.x == nucleicAcidTypeSummary.y)
# 
# result %>% 
#   filter(summary_equal==T)

NA_type = NA_type2

### number of 
NA_type.sum = NA_type %>% 
  group_by(nucleicAcidTypeSummary) %>%
  # summarise(n = n())
  count() #%>% 
  # arrange(match(light, c("photic", "BOTH", "aphotic")))



#### d3pth bins

df_SSTbins = nifhDB_cnts_T_CMAP%>% 
  filter(!is.na(sst_tblSST_AVHRR_OI_NRT)) %>% 
  mutate(
    SSTbins = cut(sst_tblSST_AVHRR_OI_NRT, breaks = c(-5,5,10,15,20,25,30,35), labels = c('<5', '5-10', '10-15', '15-20', '20-25', '25-30', '30-35'), )
  ) %>%  
  select(SSTbins, starts_with('AUID')) %>% 
  group_by(SSTbins) %>% 
  summarise(across(contains('AUID'), ~sum(., na.rm = T))) 

# Using pivot_longer and pivot_wider
df_SSTbins_T <- df_SSTbins %>%
  pivot_longer(cols = starts_with("AUID"), names_to = "AUID", values_to = "Value") %>%
  pivot_wider(names_from = SSTbins, values_from = Value)



##### photic/aphotic/both

nifhDB_cnts_T_CMAP_lng %>% distinct(AUID,photic) %>% 
  # group_by(AUID) %>% 
  mutate(
    dup = duplicated(AUID)
  ) %>%
  distinct(AUID, .keep_all = T) %>% 
  view()



lghtLvl = nifhDB_cnts_T_CMAP_lng %>% 
  distinct(AUID, photic) %>%
  group_by(AUID) %>% 
  # pivot_wider(names_from = nucleicAcidType, values_from = AUID, values_fn = length, names_prefix = "has_") %>%
  mutate(light = case_when(
    sum(photic == T) > 0 & sum(photic == F) > 0 ~ "BOTH",
    sum(photic == T) > 0 ~ "photic",
    sum(photic == F) > 0 ~ "aphotic",
    TRUE ~ "NONE"
  )) %>% 
  # mutate(light = case_when(
  #   photic == T & photic == F ~ "BOTH",
  #   photic == T  ~ 'T',
  #   photic == F  ~ 'F',
  #   TRUE ~ "NONE"
  # )) %>%
  # select(-c(has_DNA, has_RNA))
  distinct(AUID, light) 


### number of 
lghtLvl.sum = lghtLvl %>% 
  group_by(light) %>%
  # summarise(n = n())
  count() %>% 
  arrange(match(light, c("photic", "BOTH", "aphotic")))
  


# Print the ordered data frame with the statement
cat("The number of counts for each group:\n", 
    lghtLvl.sum$light, "\n ")
print(lghtLvl.sum)

# Assuming 'lghtLvl' is your data frame
lghtLvl_count <- lghtLvl %>%
  group_by(light) %>%
  count() %>% 
  arrange(match(light, c("photic", "BOTH", "aphotic")))

# # Use rowwise to iterate through each group and create the statement
# lghtLvl_count %>%
#   rowwise() %>%
#   mutate(statement = paste("The number of counts for", light, "is", n)) %>%
#   select(-n) %>%  # Remove the 'n' column, if needed
#   arrange(match(light, c("photic", "BOTH", "aphotic")))

# Loop through each row and print the statement using cat
for (i in seq_len(nrow(lghtLvl_count))) {
  cat("The number of counts for", lghtLvl_count$light[i], "is", lghtLvl_count$n[i], "\n")
}





# Distance from shore- coastal/non-coastal 

T1 = nifhDB_cnts_T_CMAP_lng %>% 
  # distinct(SAMPLEID, AUID,photic, kmToCoast, cnt) %>% 
  mutate(
    cnts_k2C = cut(kmToCoast, breaks = c(-5,200,Inf), labels = c('<200 km', '>200 km') ),
    # perc_k2C = 
  ) %>% 
  group_by(AUID) %>% 
  count(cnts_k2C)
  
  group_by(AUID) %>%
  mutate(cntGThzero =sum(cnt>0)
        ) #%>% ### number of samples for each cluster that has at least one hit in samples that are less than or equal to 200 km off the coast) #%>%  ### number of samples for each cluster that has at least one hit
  # mutate(
  #   dup = duplicated(AUID)
  # ) %>%
  # distinct(AUID, .keep_all = T) %>% 
  # view()


T1 = RAtableSummed_lng_Clst %>% 
  group_by(nifH_cluster) %>% 
  # summarise(sum(RA>0)) %>% 
  summarise(cntGThzero =sum(RA>0))%>%  ### number of samples for each cluster that has at least one hit
  ungroup() %>% 
  left_join(RAtableSummed_lng_Clst %>% 
              filter(kmToCoast<=200) %>%
              group_by(nifH_cluster) %>% 
              summarise(cntGThzero201 =sum(RA>0))) %>% ### number of samples for each cluster that has at least one hit in samples that are less than or equal to 200 km off the coast
  ungroup() %>% 
  group_by(nifH_cluster) %>% 
  mutate(fracCntGThzero201 =cntGThzero201/cntGThzero) %>% ### fraction of samples that have at least one hit and are located at least <=200 km off of shore
  ungroup() %>% 
  arrange(desc(fracCntGThzero201)) %>% 
  view()



T2 = RAtableSummed_lng_Clst %>% 
  group_by(nifH_cluster) %>% 
  # summarise(sum(RA>0)) %>% 
  summarise(RAsum =sum(RA))%>%  ### sum the RA over all samples
  ungroup() %>% 
  left_join(RAtableSummed_lng_Clst %>% 
              filter(kmToCoast<=200) %>%
              group_by(nifH_cluster) %>% 
              summarise(RAgThzero201 =sum(RA))) %>% ### sum the RA over all samples that are less than or equal to 200 km off the coast
  ungroup() %>% 
  group_by(nifH_cluster) %>% 
  mutate(fracRAgThzero201 =RAgThzero201/RAsum) %>% # fraction of summed RA of all samples located less than or equal to 200 km off the coast 
  ungroup() %>% 
  arrange(desc(fracRAgThzero201)) %>% 
  view()

T3 = T1 %>% 
  bind_cols(T2)  %>% 
  select(-nifH_cluster...5) %>% 
  rename(nifH_cluster = nifH_cluster...1) %>% 
  view()



lghtLvl = nifhDB_cnts_T_CMAP_lng %>% 
  distinct(AUID, photic) %>%
  group_by(AUID) %>% 
  # pivot_wider(names_from = nucleicAcidType, values_from = AUID, values_fn = length, names_prefix = "has_") %>%
  mutate(light = case_when(
    sum(photic == T) > 0 & sum(photic == F) > 0 ~ "BOTH",
    sum(photic == T) > 0 ~ "photic",
    sum(photic == F) > 0 ~ "aphotic",
    TRUE ~ "NONE"
  )) %>% 
  # mutate(light = case_when(
  #   photic == T & photic == F ~ "BOTH",
  #   photic == T  ~ 'T',
  #   photic == F  ~ 'F',
  #   TRUE ~ "NONE"
  # )) %>%
  # select(-c(has_DNA, has_RNA))
  distinct(AUID, light) 




### identify all ASVs that reach a threshold in X number of samples

threshold = 0.01
num_samp_filt =  5

# Group by AUID and count the number of samples each ASV appears in
asv_counts_samples <- df_long %>%
  group_by(AUID) %>%
  summarise(num_samples = sum(RA >= threshold))

# Filter ASVs that meet the criteria
filtered_asvs <- asv_counts_samples %>%
  filter(num_samples >=num_samp_filt)

# Print the selected ASVs
selected_asvs =filtered_asvs$AUID
num_selected_asvs <- length(selected_asvs)

selected_asvs_RAperSamp = selected_asvs

cat( length(selected_asvs_RAperSamp), "ASV reach relative abundances of at least", threshold, "in at least", num_samp_filt, "samples\n")
cat("Selected ASVs:", selected_asvs, "\n")






threshold_studyID = 0.01
num_samp_filt =  5

# Group by studyID and AUID and count the number of samples each ASV appears in each studyID
asv_counts_stdID <- df_long %>%
  group_by(studyID, AUID) %>%
  summarise(num_samples = sum(RA > threshold_studyID)) %>% 
  ungroup()

# Filter ASVs that meet the criteria
filtered_asvs <- asv_counts_stdID %>%
  filter(num_samples >=num_samp_filt)


# Print the selected ASVs
selected_asvs =filtered_asvs$AUID
num_selected_asvs <- length(selected_asvs)

selected_asvs_SampperStudy = selected_asvs
cat( length(selected_asvs_SampperStudy), "ASV reach relative abundances >", threshold_studyID, "in at least", num_samp_filt, "samples\n")
cat("Selected ASVs:", selected_asvs, "\n")


#### identify each study each ASV is in and make a column with the list as well as separate columns for each 

asv_inStdID <- df_long %>%
  group_by(AUID, studyID) %>%
  summarise(num_samples = RA >= 0) %>% 
  ungroup() %>% 
  # summarise()
  # filter(AUID=='AUID.23') %>% 
  distinct()

asv_inStdID = asv_inStdID %>%
  group_by(AUID) %>%
  summarize(total_num_studies = n()) %>% 
  distinct(AUID, .keep_all = T) %>% 
  view()
  


asv_inStdID_wide = asv_inStdID %>% 
  pivot_wider(id_cols = AUID, names_from = studyID, values_from = num_samples) %>% 
  # replace_na(., F ) %>% 
  view()


###

asv_inSamp <- df_long %>%
  group_by(AUID, SAMPLEID) %>%
  summarise(num_samples = RA >= 0) %>% 
  ungroup() %>% 
  # summarise()
  # filter(AUID=='AUID.23') %>% 
  distinct()

asv_inSamp = asv_inSamp %>%
  group_by(AUID) %>%
  summarize(total_num_studies = n()) %>% 
  distinct(AUID, .keep_all = T) %>% 
  view()


AUID_cnts = nifhDB_cnts_T_CMAP_lng %>% 
  group_by(AUID) %>% 
  summarise(totalCnts = sum(cnt)) %>%  ##summed all counts for each clade
  ungroup() %>% 
  mutate(percentage = totalCnts / sum(totalCnts) * 100) %>% 
  view()
  


NA_type

df_SSTbins_T

lghtLvl

asv_inStdID

asv_counts_samples

AUID_cnts

test = c('NA_type', 'df_SSTbins_T')

lghtLvl

asv_inStdID

asv_counts_samples

AUID_cnts

bigSUB = dfExpnd.fa_nifHDB %>% 
  left_join(AUID_cnts) %>% 
  
  left_join(NA_type) %>% 
  left_join(AUID_cnts) %>% 
  left_join(lghtLvl) %>% 
  left_join(asv_inStdID) %>% 
  left_join(asv_inStdID) %>% 
  left_join(df_SSTbins_T) %>% 
  
  view()


# Assuming 'annotations' is your data frame
for (df in test) {
  annotationsSUB <- left_join(annotations, df, by = "AUID")
}

# Find values that are the same in both columns
same_values <- intersect(selected_asvs_RAperSamp, selected_asvs_SampperStudy)

length(same_values)

# Find values that differ between the two columns
differing_values <- union(setdiff(selected_asvs_RAperSamp, selected_asvs_SampperStudy), setdiff(selected_asvs_SampperStudy, selected_asvs_RAperSamp))

length(union(setdiff(selected_asvs_RAperSamp, selected_asvs_SampperStudy), setdiff(selected_asvs_SampperStudy, selected_asvs_RAperSamp)))


### pull AUIDs out of main file 
dfExpnd.fa_nifHDB_filtRA = dfExpnd.fa_nifHDB %>%
    filter(AUID %in% selected_asvs_RAperSamp ) %>%
  select(AUID, sequence)
### pull AUIDs out of main file 
dfExpnd.fa_nifHDB_filtStdID = dfExpnd.fa_nifHDB %>%
  filter(AUID %in% selected_asvs_SampperStudy ) %>%
  select(AUID, sequence)


bigSUB

### pull AUIDs out of main file 
dfExpnd.fa_nifHDB_filtRA = bigSUB %>%
  filter(AUID %in% selected_asvs_RAperSamp ) 

### pull AUIDs out of main file 
dfExpnd.fa_nifHDB_filtStdID = bigSUB %>%
  filter(AUID %in% selected_asvs_SampperStudy ) 




# dir.create('all_studies/master_fasta/nifHcatalog/')

# write_csv(A2KHI.df, 'all_studies/master_fasta/A2KHI_AUID.csv')
write_csv(dfExpnd_nifHDB, paste('all_studies/master_fasta/nifHDB/nifHDB_sequences',date_stamp, '.csv', sep = ''))



write_csv(dfExpnd.fa_nifHDB_filtRA, paste('all_studies/master_fasta/nifHDB/nifHDB_filtRA_sequences',date_stamp, '.csv', sep = ''))

write_csv(dfExpnd.fa_nifHDB_filtStdID, paste('all_studies/master_fasta/nifHDB/nifHDB_filtStdID_sequences',date_stamp, '.csv', sep = ''))

nifhDB_RA_T_CMAP_Isaac = nifhDB_RA_T_CMAP %>% 
  filter(nucleicAcidType=='DNA') %>% 
  select(-(time:logFe))
  

write_csv(nifhDB_RA_T_CMAP_Isaac, paste('all_studies/master_fasta/nifHDB/nifHDB_sequences',date_stamp, '_Isaac.csv', sep = ''))
names(nifhDB_RA_T_CMAP)

# 
#   
# tempParse = nifhDB_RA_T_CMAP_lng
# 
# # Can one of you easily create: 
# #   
# #   1) an AUID table with AUIDs that are common to more than 1 study? 
# #   
# #   2) an AUID table with AUIDs representing AUIDs that reach more than 1% relative abundance in greater than 2 samples? (Or some filtering criteria that makes more sense in the context of the paper). 
# # 
# # Hopefully one of these approaches will reduce the size of the tree, while still making the point we are trying to make. 
# 
# names(tempParse)
# 
# 
# ## 1 an AUID table with AUIDs that are common to more than 1 study? 
# tempParseStdCntAUID = tempParse %>% 
#   select(SAMPLEID, studyID, AUID, RA) %>% 
#   filter(RA>0.005) %>% 
#   group_by(AUID) %>%
#   summarize(studyID_count = n_distinct(studyID)) 
# 
# 
# TargAUID_stdID2 = tempParseStdCntAUID %>% 
#   filter(studyID_count >= 2) %>% 
#   arrange(desc(studyID_count)) 
# 
# length(TargAUID_stdID2$AUID)
# 
# # 2) an AUID table with AUIDs representing AUIDs that reach more than 1% relative abundance in greater than 2 samples? (Or some filtering criteria that makes more sense in the context of the paper). 
# 
# tempParseStdCnt0.1AUID = tempParse %>% 
#   select(SAMPLEID, studyID, AUID, RA) %>% 
#   filter(RA>0.01) %>% 
#   group_by(AUID) %>%
#   summarize(SAMPLEID_count = n_distinct(SAMPLEID)) %>% 
#   arrange(desc(SAMPLEID_count)) 
# 
# TargAUID_SAMP2 = tempParseStdCnt0.1AUID %>% 
#   filter(SAMPLEID_count > 2) %>% 
#   arrange(desc(SAMPLEID_count)) 
# 
# 
# # Find values that are the same in both columns
# same_values <- intersect(TargAUID_stdID2$AUID, TargAUID_SAMP2$AUID)
# 
# length(same_values)
# 
# # Find values that differ between the two columns
# differing_values <- union(setdiff(TargAUID_stdID2$AUID, TargAUID_SAMP2$AUID), setdiff(TargAUID_SAMP2$AUID, TargAUID_stdID2$AUID))
# 
# length(union(setdiff(TargAUID_stdID2$AUID, TargAUID_SAMP2$AUID), setdiff(TargAUID_SAMP2$AUID, TargAUID_stdID2$AUID)))
# 
# ### pull AUIDs out of main file 
# dfExpnd_stdID2_Catalog = dfExpnd_Catalog %>% 
#   filter(AUID %in% TargAUID_stdID2$AUID)
# 
# ### pull AUIDs out of main file 
# dfExpnd_SAMP2_Catalog = dfExpnd_Catalog %>% 
#   filter(AUID %in% TargAUID_SAMP2$AUID)
# 
# 
# write_csv(dfExpnd_stdID2_Catalog, paste('all_studies/master_fasta/nifHcatalog/nifHcatalog_2studies_sequences',date_stamp, '.csv', sep = ''))
# 
# write_csv(dfExpnd_SAMP2_Catalog, paste('all_studies/master_fasta/nifHcatalog/nifHcatalog_0.1in3samps_sequences',date_stamp, '.csv', sep = ''))

x = nifhDB_cnts %>% 
  pull(AUID)


# Assuming your dataframe is named 'data'
nifhDB_cnts_T <- nifhDB_cnts %>%
  pivot_longer(cols = -AUID, names_to = "SAMPLEID", values_to = "Value") %>%
  pivot_wider(names_from = AUID, values_from = Value) #%>% 
# mutate(
#   SAMPLEID = str_remove(pattern = '.*___', SAMPLEID)
# )

y = hTmain_Catalog %>% 
  select(contains('AUID')) %>% 
  names()

length(x)
length(y)


# Find values that are the same in both columns
same_values <- intersect(x, y)

length(same_values)

# Find values that differ between the two columns
differing_values <- union(setdiff(x, y), setdiff(y, x))


dif_value = setdiff(y, x)

nifhDB_cnts %>% 
  # filter(AUID %in% dif_value) %>% 
  rowwise() %>% 
  summarise(sum(c_across(where(is.numeric)), na.rm = T)
            ) %>% 
  view()


hT_phtSansUnkwn.satCMAPmeta_merg_dedup %>% 
  select( dif_value) %>%
  pivot_longer(names_to = 'AUID', values_to = 'cnts', cols = everything()) %>% 
  filter(cnts>0) %>% 
  # count(AUID)
  group_by(AUID) %>% 
  mutate(
    samples = n()
  ) %>% 
  # rowwise() %>% 
  # summarise(sum(c_across(where(is.numeric)), na.rm = T)
  # summarise(across(everything(), ~sum(., na.rm = TRUE), .names = "sum_{.col}")
  # ) %>% 
  view() 

df <- df %>%
  mutate(across(everything(), ~sum(., na.rm = TRUE), .names = "sum_{.col}"))



hTmain_Catalog$AUID.100
sum(nifhDB_cnts$AUID=='AUID.100')

length(union(setdiff(x, y), setdiff(y, x)))















tempPht = nifhDB_cnts_T_CMAP %>% 
  filter(photic==T)

tempApht = nifhDB_cnts_T_CMAP %>% 
  filter(!photic==T)

tempApht %>% 
  summarise(across(where(contains("AUID")), ~sum(.)))

nifhDB_cnts_T_CMAP_lng %>% 
  # summarise(across(where(is.numeric), sum))
  summarise(across(contains("AUID"), ~sum(.)))

tempPht%>% 
  # summarise(across(where(is.numeric), sum))
  summarise(across(contains("AUID"), ~sum(.)))

tempApht %>% 
  # summarise(across(where(is.numeric), sum))
  summarise(across(contains("AUID"), ~sum(.))) %>% 
  filter(.>0)

YOU ARE HERERHERHERHEHRERH

x = tempPht %>% 
  select(contains('AUID')) %>% 
  names()


# Assuming your dataframe is named 'data'
nifhDB_cnts_T <- nifhDB_cnts %>%
  pivot_longer(cols = -AUID, names_to = "SAMPLEID", values_to = "Value") %>%
  pivot_wider(names_from = AUID, values_from = Value) #%>% 
# mutate(
#   SAMPLEID = str_remove(pattern = '.*___', SAMPLEID)
# )

y = tempApht %>% 
  select(contains('AUID')) %>% 
  names()

length(x)
length(y)


# Find values that are the same in both columns
same_values <- intersect(x, y)

length(same_values)

# Find values that differ between the two columns
differing_values <- union(setdiff(x, y), setdiff(y, x))


dif_value = setdiff(y, x)

nifhDB_cnts %>% 
  # filter(AUID %in% dif_value) %>% 
  rowwise() %>% 
  summarise(sum(c_across(where(is.numeric)), na.rm = T)
  ) %>% 
  view()

