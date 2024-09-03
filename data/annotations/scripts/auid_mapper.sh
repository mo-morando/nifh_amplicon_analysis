#!/bin/bash
### makes key specific to old version of nifHDB
# analysis/Jmags/mappingOLD2NEWauids/fastaMapper.sh current.asv2auid.fasta auid.filteredRA.fasta >fasta_id_map2.tsv
# shell realpath"data/annotations/scripts/fastaMapper.sh"
# "~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/master_fasta/nifHDB/auid.filtered.nifHDB.fasta"
### makes key specific to nifHDB
# fastaMapper.sh current.asv2auid.fasta auid.filtered.nifHDB.fasta >nifHDB_old2new_auids_map.tsv

/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/scripts/fastaMapper.sh ~/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/master_fasta/nifHDB/auid.filtered.nifHDB.fasta
