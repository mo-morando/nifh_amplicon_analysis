### set working directory
setwd("/Users/mo/Projects/nifH_amp_project/myWork")

## load libraries
library("tidyverse")

### _ Loading in the data
cat("Load in the data")

### annotation file with consensus ID already intact
#### read in annotations table
# ## - new annotation file
# _# updated by me 20240219
# _# updated by me 20240219

annotations_new_auids <- read_csv("data/annotations/nifhDB/merged_database/annotation_tab.csv")

annoNifHDB_updt <- annotations_new_auids %>%
  rename(AUID = AUID_OLD)

# _# updated by me 20240219
# _# updated by me 20240219


#### ! Testing my scripts/workflow with an annotation file produced by the nifH
## ! worflow from thallasophile
## ! BUT they do not have the UCYN-A oligos
# annotations_new_auids <- read_tsv("data/annotations/nifhDB/jonathan/FOR_MO_ESSD/auids.annot.tsv")
# annotations_new_auids <- read_csv("annotationsall_updatedauids_nifhdb.csv")
# annotations_new_auids <- read_csv("annotationsall_updatedauids_nifhdb_consensus_tax_annotation.csv")

# annoNifHDB_updt <- annotations_new_auids %>%
#   select(-AUID) %>%
#   rename(AUID = AUID_OLD)
# ## !! Remove down to here


## - make CON column as the consesus id column since this is what the rest of
## then make the nifH_cluster
## All used in the subsequent scripts
annoNifHDB_updt <- annoNifHDB_updt %>%
  mutate(
    CON = consensus_id,
    nifH_cluster = if_else(
      cluster == 3 | cluster == 4 | cluster == 2,
      as.character(cluster),
      subcluster
    ),
    nifH_cluster = if_else(
      grepl("UCYN-B|croco", ncd_cyano_id, ignore.case = T),
      "1B",
      nifH_cluster
    ),
    nifH_cluster = if_else(
      nifH_cluster == "1J" | nifH_cluster == "1K",
      "1J/1K",
      nifH_cluster
    ),
    nifH_cluster = if_else(
      nifH_cluster == "1P" | nifH_cluster == "1O",
      "1O/1P",
      nifH_cluster
    ),
    nifH_cluster = if_else(
      is.na(nifH_cluster) | subcluster == "nan",
      "unknown",
      nifH_cluster
    )
  )

# ! TODO: The consensus ID/CON need to be cleaned up
## _ MAGs need to be identified
## _ TARA prefixes should be removed or MAG### should come first
## _ Lots of these need to be made more broad or another column needs to identify things like this is a gamma A and this is a Psuedomonas
# annoNifHDB_updt %>%
#   select(consecutive_id, subcluster, nifH_cluster) %>%
#   distinct(CON, .keep_all = T) %>%
#   view()

## add a column for parsing nifH clusters of interest
## variable with clusters of interest
clusters_of_interest <- c(
  "1A",
  "1J/1K",
  "1O/1P",
  "3",
  "1G",
  "1B",
  "4",
  "2"
)
## add new column
annoNifHDB_updt <- annoNifHDB_updt %>%
  mutate(
    cluster_stats = ifelse(test = nifH_cluster %in% clusters_of_interest,
      yes = nifH_cluster,
      no = "other"
    )
  )


# _##  remove ASVs that cannot be resolved and Syne-like
# _##  remove ASVs that cannot be resolved and Syne-like
#### Synechococcus like nifh sequences make it through to this point but the
## annotaions help identify these.

## - make Synechococcus auid key
grep_text <- "Synechococcus"

syne_auid_key <- annoNifHDB_updt %>%
  filter(grepl(grep_text, genome879_nifh_sseqid, ignore.case = TRUE)) %>%
  pull(AUID)

## - make unknown auid key
unknown_con <- "unknownnan"
unknown_cluster <- "unknown"

unknownnan_auid_key <- annoNifHDB_updt %>%
  filter(CON %in% unknown_con | nifH_cluster %in% unknown_cluster) %>%
  pull(AUID)

## - remove these from the annotation table
annoNifHDB_updt <- annoNifHDB_updt %>%
  filter(!AUID %in% syne_auid_key) %>%
  filter(!AUID %in% unknownnan_auid_key)

dim(annoNifHDB_updt)

# annoNifHDB_updt %>%
#     filter(is.na(nifH_cluster) | nifH_cluster == "nan") %>%
#     view()

cat("Load in annotation file 'annoNifHDB_updt' the data")

## - Add groups for stats and plotting
annoNifHDB_updt <- annoNifHDB_updt %>%
  mutate(
    group1 = if_else(nifH_cluster %in% c("1A", "3"), "1A C3", nifH_cluster),
    group2 = if_else(nifH_cluster %in% c("1A", "3", "1J/1K"), "1A,  3,  &  1J/1K", nifH_cluster),
    group3 = if_else(nifH_cluster %in% c("1A", "3"), "1A C3",
      if_else(nifH_cluster %in% c("1G", "1J/1K"), "1B 1J/1K", nifH_cluster)
    ),
    group4 = if_else(nifH_cluster %in% c("1A", "3", "1J/1K"), "1A/C3 1J/1K", "1B & 1G"),
    CyanoCON = if_else(nifH_cluster == "1B", CON, nifH_cluster),
    crocoCMB = if_else(grepl("croco", CON, ignore.case = T) & CON == "Crocosphera_DQ118216_Moisander", "CDQmois",
      if_else(grepl("croco", CON, ignore.case = T), "CrocoCMB", CON)
    ),
    CyanoGroups = if_else(CON == "UCYN-A3" | CON == "UCYN-A4", "A3-A4", CON),
    CyanoGroupsII = if_else(grepl("Tricho", CON), "Trichodesmium sp.", CON),
    CyanoGroupsIII = if_else(grepl("Tricho", CON), "Trichodesmium sp.", CyanoGroups),
    CyanoGroupsIV = if_else(CON == "UCYN-A3" | CON == "UCYN-A1", "A1-A3",
      if_else(CON == "UCYN-A2" | CON == "UCYN-A4", "A2-A4", CON)
    )
  )

## _### CMAP environmental data #####
## _### CMAP environmental data #####

### - CMAP
CMAP_20230719 <- read_csv("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/merged_metadata/nifHDB/firstAttempt/20230719_colocalized_nifH.csv")

CMAP_coloc <- CMAP_20230719

# _ FIXME: this needs to be clean up and moved to other helper script to process all these files and make keys for photic, nucleic acids, etc.
#### ! read this old file in for now so you can fix/add some extra data
## ! eventually you will get all this made by the main script and this can be removed
temp <- read_csv("/Users/mo/mmorando@ucsc.edu_gDrive/My Drive/data/amplicon_review/all_studies/merged_metadata/RAtabCMAPmeta_merg20230126.csv")

kmToCoast_key <- temp %>%
  distinct(lat, lon, .keep_all = T) %>%
  select(lat, lon, kmToCoast)

photic_key <- temp %>%
  distinct(lat, lon, .keep_all = T) %>%
  select(lat, lon, photic)

rm(temp)

CMAP_coloc <- CMAP_coloc %>%
  left_join(kmToCoast_key)

CMAP_coloc <- CMAP_coloc %>%
  left_join(photic_key) %>%
  select(SAMPLEID, studyID, lat, lon, photic, everything()) # %>%
# filter(is.na(photic))

photic_samples_key <- CMAP_coloc %>%
  filter(photic == TRUE) %>%
  pull(SAMPLEID)

# remove_aphotic_samples <- function(df, photic_samples = photic_samples_key) {
#   subbed_df <- df %>%
#     filter(SAMPLEID %in% photic_samples)

#   return(subbed_df)
# }


CMAP_coloc <- CMAP_coloc %>%
  # filter((is.na(photic) & studyID %in% c("Raes_2020", "Mulholland_2018"))) %>%
  # filter((is.na(photic) & studyID %in% c("Shiozaki_2017")) & depth==104) #%>%
  mutate(photic = if_else((is.na(photic) & studyID %in% c("Shiozaki_2017") & depth == 104), FALSE,
    if_else((is.na(photic) & studyID %in% c("Shiozaki_2017") & depth == 3), TRUE, photic)
  )) %>%
  mutate(photic = if_else((is.na(photic) & studyID %in% c("Raes_2020", "Mulholland_2018")), TRUE, photic)) %>%
  select(photic, everything()) # %>%
# filter(is.na(photic))
#




## _  FIXME: There is a major problem with the distance to shore (kmToCoast)
## _ Many distances are much shorter than they should be due to random landmasses that we would not include, e.g., super small islands
## _ TODO: Script determining kmToCoast needs to be fixed
## _ "Ras_2020", "Mulholland_2018", "Shiozaki_2017" definitely need to be fixed
# merged <- merged %>%
#     # filter(lat == 40.4755) %>%
#     # filter(is.na(kmToCoast)) %>%
#     filter(studyID %in% c("Ras_2020", "Mulholland_2018", "Shiozaki_2017")e) %>%
#     select(SAMPLEID, studyID, lat, lon, kmToCoast, everything())

# write_csv(merged, "test_merge.csv")

# write_csv(merged_photic, "test_merge_photic.csv")


## Once CMAP file is all fixed and clean, add the extra columns for parsing

CMAP_coloc <- CMAP_coloc %>%
  mutate(
    hemi = cut(.$lat,
      breaks = c(-Inf, 0, Inf),
      labels = c("southernHemi", "northernHemi")
    ),
    hemi = factor(hemi, c("northernHemi", "southernHemi")),
    lat_abs = abs(lat),
    Size_fraction = if_else(
      Size_fraction %in% c("whole", "Sterivex", "0.22") | is.na(Size_fraction),
      "0.2",
      Size_fraction
    ),
    date = str_remove(time, " .+$"),
    month = str_remove_all(str_extract(date, "-.+-"), "-"),
    season = if_else(month %in% c("03", "04", "05") & hemi == "northernHemi", "spring",
      if_else(month %in% c("09", "10", "11") & hemi == "northernHemi", "fall",
        if_else(month %in% c("06", "07", "08") & hemi == "northernHemi", "summer",
          if_else(month %in% c("12", "01", "02") & hemi == "northernHemi", "winter",
            if_else(month %in% c("03", "04", "05") & hemi == "southernHemi", "fall",
              if_else(month %in% c("09", "10", "11") & hemi == "southernHemi", "spring",
                if_else(month %in% c("06", "07", "08") & hemi == "southernHemi", "winter",
                  if_else(month %in% c("12", "01", "02") & hemi == "southernHemi", "summer", "XXX")
                )
              )
            )
          )
        )
      )
    ),
    season = factor(season, c("winter", "spring", "summer", "fall")),
    CMAP_NP_darwin = CMAP_NO3_darwin_clim_tblDarwin_Nutrient_Climatology / CMAP_PO4_darwin_clim_tblDarwin_Nutrient_Climatology,
    CMAP_NP_Pisces_NRT = CMAP_NO3_tblPisces_NRT / CMAP_PO4_tblPisces_NRT,
    CMAP_NP_WOA_clim = CMAP_nitrate_WOA_clim_tblWOA_Climatology / CMAP_phosphate_WOA_clim_tblWOA_Climatology,
    CMAP_NP_WOA_clim = CMAP_n_an_clim_tblWOA_2018_1deg_Climatology / CMAP_p_an_clim_tblWOA_2018_1deg_Climatology
  ) %>%
  rename_with(~ str_remove(string = ., pattern = "CMAP_")) %>%
  rename_with(~ coalesce(str_extract(string = ., pattern = "^[A-Za-z0-9]+_darwin"), .)) %>%
  rename_with(~ str_replace(string = ., pattern = "clim_tblWOA_2018_qrtdeg_Climatology", replacement = "WOA_2018_qrtdeg")) %>%
  rename_with(~ str_replace(string = ., pattern = "clim_tblWOA_2018_1deg_Climatology", replacement = "WOA_2018_1deg")) %>%
  rename_with(~ str_replace(string = ., pattern = "WOA_clim_tblWOA_Climatology", replacement = "WOA_clim")) %>%
  rename_with(~ str_replace(string = ., pattern = "clim_tblWOA_2018_MLD_qrtdeg_Climatology", replacement = "WOA_2018_MLD_qrtdeg")) %>%
  rename_with(~ str_replace(string = ., pattern = "^C_", replacement = "conductivity_")) %>%
  rename_with(~ str_replace(string = ., pattern = "^i_", replacement = "sigma_")) %>%
  rename_with(~ str_replace(string = ., pattern = "^t_", replacement = "temp_")) %>%
  rename_with(~ str_replace(string = ., pattern = "^A_", replacement = "AOU_")) %>%
  rename_with(~ str_replace(string = ., pattern = "^O_", replacement = "O2_sat_")) %>%
  rename_with(~ str_replace(string = ., pattern = "^n_", replacement = "NO3_")) %>%
  rename_with(~ str_replace(string = ., pattern = "^p_", replacement = "phosphate_")) %>%
  rename_with(~ str_replace(string = ., pattern = "^si_", replacement = "silica_")) %>%
  rename_with(~ str_replace(string = ., pattern = "^M_", replacement = "MLD_")) %>%
  rename_with(~ str_replace(string = ., pattern = "^s_", replacement = "salinity_")) %>%
  mutate(
    geoRegion = cut(.$lat_abs,
      breaks = c(-1, 23, 35, 66, 100),
      labels = c("Eq/Trop", "Sub-Trop", "Temp", "Poles")
    ),
    geoRegion = factor(geoRegion, c("Eq/Trop", "Sub-Trop", "Temp", "Poles")),
    logFe = log(Fe_tblPisces_NRT)
  )

### add Ocean regions for each study ID
studyid_regions <- read_csv("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/tables/Table1/studyid_regions.csv")

### identy studies from Southern Ocean
CMAP_coloc <- CMAP_coloc %>%
  left_join(studyid_regions) %>%
  select(photic, SAMPLEID, studyID, study_ocean, everything()) %>%
  mutate(
    ocean = if_else(lat_abs >= 60 & hemi == "southernHemi", "Southern", study_ocean)
  )

# ### Size fractions
# CMAP_coloc <- CMAP_coloc %>%
#     mutate(
#         Size_fraction = if_else(Size_fraction %in% c("whole", "Sterivex", "0.22") | is.na(Size_fraction), "0.2", Size_fraction)
#     ) # %>%
# # select(studyID, Size_fraction, everything()) %>%
# # view()


# _### FASTA FILE ####
# _### FASTA FILE ####
# _### FASTA FILE ####


### - nifh database fasta
fastaFile_DB <- Biostrings::readDNAStringSet("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/master_fasta/nifHDB/auid.filtered.nifHDB.fasta")
seq_name <- names(fastaFile_DB)
sequence <- paste(fastaFile_DB)
nifHDB.fa_df <- data.frame(seq_name, sequence)


dfExpnd.fa <- nifHDB.fa_df %>%
  mutate(
    AUID = str_extract("AUID\\.\\d+", string = nifHDB.fa_df$seq_name),
    studies = str_remove("AUID\\.\\d+", string = nifHDB.fa_df$seq_name)
  ) %>%
  ## remove unknown and syne-like sequences
  filter(!AUID %in% syne_auid_key) %>%
  filter(!AUID %in% unknownnan_auid_key)

dfExpnd.fa_nifHDB <- dfExpnd.fa %>%
  #   filter(AUID %in% TargetAUIDs ) %>%
  select(AUID, sequence)


### convert into actual fasta file
seqs <- Biostrings::DNAStringSet(dfExpnd.fa_nifHDB$sequence)

# Assign AUID values as names to the seqs vector
names(seqs) <- dfExpnd.fa_nifHDB$AUID

## write out as fasta
# Biostrings::writeXStringSet(
#     x = seqs,
#     filepath = paste("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/master_fasta/nifHDB/nifHDB_filtered", date_stamp, ".fa", sep = ""),
#     format = "fasta"
# )

Biostrings::writeXStringSet(
  x = seqs,
  filepath = paste("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/master_fasta/nifHDB/nifHDB_filtered.fa", sep = ""),
  format = "fasta"
)

Biostrings::writeXStringSet(
  x = seqs,
  filepath = paste("../analysis/files/nifHDB_filtered.fa", sep = ""),
  format = "fasta"
)

### - fasta
# dfExpnd.fa_nifHDB
# str(dfExpnd.fa_nifHDB)

# dir.create('all_studies/master_fasta/nifHcatalog/')

# ##* write out as dataframe (.csv)
# # write_csv(A2KHI.df, '~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies//master_fasta/A2KHI_AUID.csv')
# write_csv(dfExpnd.fa_nifHDB, paste("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies//master_fasta/nifHDB/nifHDB_sequences", date_stamp, ".csv", sep = ""))


### - count/abundance tables
nifhDB_cnts <- read.table("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/master_fasta/nifHDB/auid.abundances.filtered.nifHDB.tsv", header = TRUE, sep = "\t", row.names = 1) %>%
  rownames_to_column("AUID") %>%
  rename_all(~ str_remove(., pattern = ".*___")) %>%
  as_tibble()

dim(nifhDB_cnts)
names(nifhDB_cnts)

## -# remove unknown and syne-like sequences
#* # remove unknown and syne-like sequences
nifhDB_cnts <- nifhDB_cnts %>%
  filter(!AUID %in% syne_auid_key) %>%
  filter(!AUID %in% unknownnan_auid_key)

### remove samples from NEMO that do not fit the criteria for the nifh database
(nifhDB_cnts <- nifhDB_cnts %>%
  select(-matches("Turk\\d+\\.e")) %>%
  select(-matches("Harding229\\.66705_S229|Harding230\\.66706_S230|Harding231\\.66709_S231")))


## -# Create a relative abundance table
##### calculate relative abundance
nifhDB_RA <- nifhDB_cnts %>%
  mutate(across(where(is.numeric), ~ . / sum(., na.rm = T)))

## check to verify that this worked
## Values should all be 1
nifhDB_RA %>%
  summarise(across(where(is.numeric), sum)) %>%
  view()


# transform data
nifhDB_cnts_T <- nifhDB_cnts %>%
  pivot_longer(cols = -AUID, names_to = "SAMPLEID", values_to = "Value") %>%
  pivot_wider(names_from = AUID, values_from = Value) # %>%
# mutate(
#   SAMPLEID = str_remove(pattern = '.*___', SAMPLEID)
# )

nifhDB_RA_T <- nifhDB_RA %>%
  pivot_longer(cols = -AUID, names_to = "SAMPLEID", values_to = "Value") %>%
  pivot_wider(names_from = AUID, values_from = Value) # %>%
# mutate(
#   SAMPLEID = str_remove(pattern = '.*___', SAMPLEID)
# )


### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
cat("
Files loaded:
dfExpnd.fa_nifHDB
nifhDB_cnts
nifhDB_cnts_T
nifhDB_RA
nifhDB_RA_T
annoNifHDB_updt
CMAP_coloc

")
cat("Done loading script!!!
")

cat("Woooooooohooooooo!!!")
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data

### collection

### These are the files
# dfExpnd.fa_nifHDB

# nifhDB_cnts
# nifhDB_cnts_T

# nifhDB_RA
# nifhDB_RA_T

# annoNifHDB_updt

# CMAP_coloc
