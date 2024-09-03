#!/bin/bash

original_headers=$(head -n 1 /Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/UCYNA_oligos/filt_final_output_csv)
# original_string="sseqid length qcovs pident mismatch gapopen qstart qend sstart send evalue bitscore qseqid"
header_prefix="uncyna_oligos_"
# modified_string=$(echo "$original_headers" | awk '{ for (i=1; i<=NF; i++) $i="'${header_prefix}'" $i } 1')
modified_string=$(echo "$original_headers" | awk -v prefix="$header_prefix" 'BEGIN {FS=OFS=","} { for (i=1; i<=NF; i++) $i= prefix $i } 1')
# modified_string=$(echo "$original_headers" | awk 'BEGIN {FS=OFS=","} { for (i=1; i<=NF; i++) $i="'${header_prefix}'" $i } 1')
echo "$modified_string"
