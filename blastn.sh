#!/usr/bin/bash

# Example usage
# input_fasta_file="/Users/mo/Projects/nifH_amp_project/myWork/test.fasta"
input_fasta_file="/Users/mo/Projects/nifH_amp_project/myWork/analysis/Jmags/mappingOLD2NEWauids/auid.filtered.nifHDB.fasta"
# blast_database = "./UCYN_A_oligos/UCYN-A_oligoreps.fasta"
blast_database="/Users/mo/Projects/nifH_amp_project/myWork/data/databases/UCYN-A_oligoreps_db"
# identity_threshold=97
# blast_type="blastn"
num_threads=4
output_csv=auid.filtered.nifHDB_file.csv

# Define column headers
# column_headers="qseqid,sseqid,alignment_length,qcovs,PID,mismatches,gap_openings,q_start,q_end,s_start,s_end,evalue,bit_score"
# column_headers="sseqid,length,qcovs,pident,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qseqid"

column_headers_blast="sseqid length qcovs pident mismatch gapopen qstart qend sstart send evalue bitscore qseqid"
column_headers_tsv=$(echo "$column_headers_blast" | tr ' ' '\t')
column_headers_csv=$(echo "$column_headers_blast" | tr ' ' ',')

# column_headers_tsv="sseqid\tlength\tqcovs\tpident\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqseqid"
# column_headers_csv="qseqid,sseqid,alignment_length,qcovs,PID,mismatches,gap_openings,q_start,q_end,s_start,s_end,evalue,bit_score"

# Run BLAST with format 6 and include headers
blastn -db "$blast_database" -query "$input_fasta_file" -max_target_seqs 1 -outfmt "6 $column_headers_blast" -num_threads "$num_threads" -out output_file.tsv

# Run BLAST with format 10 and include headers
blastn -db "$blast_database" -query "$input_fasta_file" -max_target_seqs 1 -outfmt "10 $column_headers_blast" -out "$output_csv" -num_threads "$num_threads"

# # Run BLAST with format 6 and include headers
# blastn -db "$blast_database" -query "$input_fasta_file" -max_target_seqs 1 -outfmt "6 qseqid sseqid length qcovs pident mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads "$num_threads" -out output_file.tsv

# # Run BLAST with format 10 and include headers
# blastn -db "$blast_database" -query "$input_fasta_file" -max_target_seqs 1 -outfmt "10 qseqid sseqid length qcovs pident mismatch gapopen qstart qend sstart send evalue bitscore" -out "$output_csv" -num_threads "$num_threads"

# Add column headers to the output file
echo "10 $column_headers_blast"
echo "$column_headers_csv" >temp_output_file.csv
cat $output_csv >>temp_output_file.csv
mv temp_output_file.csv $output_csv

# Add column headers to the output file
echo -e "6 $column_headers_blast"
echo -e "$column_headers_tsv" >temp_output_file.tsv
cat output_file.tsv >>temp_output_file.tsv
mv temp_output_file.tsv output_file.tsv

## sort out significant matches to a filtered table

# # Add column headers to the output file
# echo "10 sseqid length qcovs pident mismatch gapopen qstart qend sstart send evalue bitscore qseqid"
# echo "$column_headers" >temp_output_file.csv
# cat $output_csv >>temp_output_file.csv
# mv temp_output_file.csv $output_csv

# # Add column headers to the output file
# echo -e "6 sseqid length qcovs pident mismatch gapopen qstart qend sstart send evalue bitscore qseqid"
# echo -e "$column_headers_t" >temp_output_file.tsv
# cat output_file.tsv >>temp_output_file.tsv
# mv temp_output_file.tsv output_file.tsv
