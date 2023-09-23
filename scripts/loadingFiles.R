library("tidyverse")

### _ Loading in the data
cat("Load in the data")


### - nifh database fasta
fastaFile_DB <- Biostrings::readDNAStringSet("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/master_fasta/nifHDB/auid.filtered.nifHDB.fasta")
seq_name <- names(fastaFile_DB)
sequence <- paste(fastaFile_DB)
nifHDB.fa_df <- data.frame(seq_name, sequence)


dfExpnd.fa <- nifHDB.fa_df %>%
    mutate(
        AUID = str_extract("AUID\\.\\d+", string = nifHDB.fa_df$seq_name),
        studies = str_remove("AUID\\.\\d+", string = nifHDB.fa_df$seq_name)
    )


dfExpnd.fa_nifHDB <- dfExpnd.fa %>%
    #   filter(AUID %in% TargetAUIDs ) %>%
    select(AUID, sequence)


### convert into actual fasta file
seqs <- Biostrings::DNAStringSet(dfExpnd.fa_nifHDB$sequence)

# Assign AUID values as names to the seqs vector
names(seqs) <- dfExpnd.fa_nifHDB$AUID

### write out as fasta
# Biostrings::writeXStringSet(
#     x = seqs,
#     filepath = paste("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/master_fasta/nifHDB/nifHDB_sequences_CFH", date_stamp, ".fa", sep = ""),
#     format = "fasta"
# )

### - fasta
dfExpnd.fa_nifHDB

# dir.create('all_studies/master_fasta/nifHcatalog/')

# ##* write out as dataframe (.csv)
# # write_csv(A2KHI.df, '~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies//master_fasta/A2KHI_AUID.csv')
# write_csv(dfExpnd.fa_nifHDB, paste("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies//master_fasta/nifHDB/nifHDB_sequences", date_stamp, ".csv", sep = ""))


### - count/abundance tables
nifhDB_cnts <- read.table("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/master_fasta/nifHDB/auid.abundances.filtered.nifHDB.tsv", header = TRUE, sep = "\t", row.names = 1) %>%
    rownames_to_column("AUID") %>%
    rename_all(~ str_remove(., pattern = ".*___"))

dim(nifhDB_cnts)
names(nifhDB_cnts)

### remove samples from NEMO that I have no idea what they are

nifhDB_cnts <- nifhDB_cnts %>%
    select(-matches("Turk\\d+\\.e")) %>%
    select(-matches("Harding229\\.66705_S229|Harding230\\.66706_S230|Harding231\\.66709_S231")) %>%
    view()


##### calculate relative abundance
nifhDB_RA <- nifhDB_cnts %>%
    # column_to_rownames('AUID') %>%
    # group_by(SAMPLEID) %>%
    mutate(across(where(is.numeric), ~ . / sum(., na.rm = T)))
# mutate(relative_abundance = count / sum(count)) %>%
# pivot_wider(names_from = ASV, values_from = relative_abundance)

nifhDB_RA %>%
    summarise(across(where(is.numeric), sum)) %>%
    # filter(.!=1) %>%
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

### - count files
nifhDB_cnts
nifhDB_cnts_T
nifhDB_cnts_T_lng


### annotation files
#### read in annotations table
## - new annotation file
annoNifHDB_updt <- read_tsv("~/mmorando@ucsc.edu_gDrive/My Drive/data/amplicon_review/all_studies/master_annotation/filtered_annotations/nifHDB/auids.annot_nifhDB.tsv")

cat("Load in annotation file", annoNifHDB_updt, "the data")


### - CMAP
CMAP_20230719 <- read_csv("~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/merged_metadata/nifHDB/firstAttempt/20230719_colocalized_nifH.csv")
CMAP_coloc <- CMAP_20230719

### - CMAP
CMAP_coloc


### collection

### These are the files
dfExpnd.fa_nifHDB

nifhDB_cnts
nifhDB_cnts_T

nifhDB_RA
nifhDB_RA_T

annoNifHDB_updt

CMAP_coloc
