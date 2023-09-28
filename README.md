I am trying to put together a full annotation file for the nifH DB. Currently, I am using different AUIDs than Jonathan so they need to be matched up properply. 

My files are:
Fasta : input_fasta_file="/Users/mo/Projects/nifH_amp_project/myWork/analysis/Jmags/mappingOLD2NEWauids/auid.filtered.nifHDB.fasta"
Annotations: **Must add** UCYN-A oligo reps
UCYN-A oligo rep fasta: "/Users/mo/Projects/nifH_amp_project/myWork/UCYN-A_oligoreps.fasta"


Jonathans files:
Annotations for ARB/genomes879/CART/NCD/cyano catalog: /Users/mo/Library/CloudStorage/GoogleDrive-mmorando@ucsc.edu/My Drive/data/amplicon_review/all_studies/master_annotation/filtered_annotations/nifHDB/auids.annot_nifHDB_newAUID.tsv

### makes key specific to nifHDB
fastaMapper.sh current.asv2auid.fasta auid.filtered.nifHDB.fasta > nifHDB_old2new_auids_map.tsv
This file nifHDB_old2new_auids_map.tsv can be used to map AUIDs between old AUIDs (my files) and new AUIDs (Jonathan files)
/Users/mo/Projects/nifH_amp_project/myWork/analysis/Jmags/mappingOLD2NEWauids/nifHDB_old2new_auids_map.tsv

In the end, I want AUID that reflect the old convention. So I will need to:
1. Make new annotation file that contains