#!/bin/bash

#Define column headers
column_headers_blast="sseqid length qcovs pident mismatch gapopen qstart qend sstart send evalue bitscore qseqid"

# Function to run BLAST
# Parameters:
#   - input_fasta_file: Path to input FASTA file
#   - blast_database: Path to BLAST database
#   - blast_type: Type of BLAST to run
#   - num_threads: Number of threads for BLAST
#   - output_csv: Path to output CSV file
run_blast() {
    input_fasta_file="$1"
    blast_database="$2"
    blast_type="$3"
    num_threads="$4"
    output_csv="$5"

    # Run BLAST with format 10 and include headers
    "$blast_type" -db "$blast_database" -query "$input_fasta_file" -max_target_seqs 1 -outfmt "10 $column_headers_blast" -out "$output_csv" -num_threads "$num_threads"

    # Check if BLAST execution was successful
    if [ $? -eq 0 ]; then
        echo "BLAST execution was successful! :)"
    else
        echo "BLAST execution fail. :("
        exit 1
    fi

}

# Function to add column headers to CSV file
# Parameters:
#   - input_csv: Path to input CSV file
#   - output_csv: Path to output CSV file
add_headers_csv() {
    input_csv="$1"
    output_csv="$2"
    # header_prefix="$3"

    # # Append prefix specific to database being generated
    # column_headers_blast_sub=$(echo "$column_headers_blast" | awk '{ for (i=1; i<=NF; i++) $i="'${header_prefix}'" $i } 1')

    # Replace spaces with commas to create CSV headers
    column_headers_csv=$(echo "$column_headers_blast" | tr ' ' ',')

    echo "$column_headers_csv" >temp_output_file.csv
    cat "$input_csv" >>temp_output_file.csv
    mv temp_output_file.csv "$output_csv"
    if [ $? -eq 0 ]; then
        echo "Added headers to csv."
        rm "$input_csv"
        if [ $? -eq 0 ]; then
            echo "Removed input csv."
        else
            echo "Failed to remove input csv."
            exit 1
        fi
    else
        echo "Failed to add headers to csv."
        exit 1
    fi
}

## Check if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <input_fasts> <blast_database> <blast_type> <num_threads> <output_csv> <final_output_csv>"
    exit 1
fi

## Call the functions with provided arguments
run_blast "$1" "$2" "$3" "$4" "$5"
add_headers_csv "$5" "$6" #"$7"

# bash ../../../scripts/blast_function.sh "/Users/mo/Projects/nifH_amp_project/myWork/analysis/Jmags/mappingOLD2NEWauids/auid.filtered.nifHDB.fasta" "/Users/mo/Projects/nifH_amp_project/myWork/data/databases/UCYNA_oligoreps/ucyna_oligoreps_db/ucyna_oligoreps_db"   "blastn" 4 "output_blast.csv"
