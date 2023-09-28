#### Jonathan changed all the AUIDs when he updated the chimera checking. Now most files, including all my working files, figures, and tables have the old AUIDs. He wrote a scripts, fastaMapper.R, to match of the sequences between the files and provide a key to match old and new AUIDs.

### using current.asv2auid.fasta which is the full output from DADA2 pipeline and has the updated AUIDs (Jonthan's annotation files)

### using auid.filteredRA.fasta has all the old AUIDs

### makes key specific to old version of nifHDB
fastaMapper.sh current.asv2auid.fasta auid.filteredRA.fasta > fasta_id_map2.tsv

### makes key specific to nifHDB
fastaMapper.sh current.asv2auid.fasta auid.filtered.nifHDB.fasta > nifHDB_old2new_auids_map.tsv

This file nifHDB_old2new_auids_map.tsv can be used to map AUIDs between old AUIDs (my files) and new AUIDs (Jonathan files)