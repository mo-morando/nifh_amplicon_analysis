Making annotation pipeline

## Making annotation files
## make blast databases
### I think people could use this option and its simple to add
data/annotations/scripts/make_blastdb_fun.py

## make 2 seperate blasts and CART
### Will use Jonathans existing scripts
## plus 2 more blasts of NCD and oligos
data/annotations/scripts/blastn.sh

## sorting each blast
data/annotations/scripts/temp_filter_blast_output_pandas.py

## Merging all outputs together into single file
data/annotations/scripts/parsing_annotations.py
data/annotations/scripts/temp_parsingAUIDs_function.py
## consensus ID
data/annotations/scripts/make_consensus_taxonomy.py

## cleaning up everything
### Remove the many columns we no longer want
### Rename columns to what we want

### clean up name
### consensus ID/CON
    ## need a better, cleaner name for NCDs and oligos
### nifh_clusters 
    ## combing them 
## Removing unknown
## Removing syn


## Fasta file to use of nifH database:
/Users/mo/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/master_fasta/nifHDB/auid.filtered.nifHDB.fasta


## libraries used
snakemake
python
blast
libraries
    argparse
    subprocess
    sys
    pandas