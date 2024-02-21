#!/bin/bash

#-# Function to add prefixes to each headers to a csv file
add_headers_prefix_csv() {
    # input csv file
    input_csv="$1"
    # output csv file
    header_output_csv="$2"
    # Prefix to be added to each header
    header_prefix="$3"
    # Headers to keep

    ### Some checks

    # check if input csv file exists
    if [ ! -f "$input_csv" ]; then
        echo "Error: Input csv file '$input_csv' not found"
        exit 1
    fi

    # check if input file is recongnized as a csv
    file_type=$(file -b --mime-type "$input_csv")
    if [[ "$file_type" != "text/csv" ]]; then
        echo "Error: Input file is not recognized as a csv file. Please make sure you are using the proper format."
        exit 1
    fi

    # Grab the headers from the input csv
    column_headers_blast=$(head -n 1 "$input_csv")

    # Append prefix specific to database being generated to each header
    column_headers_blast_sub=$(echo "$column_headers_blast" | awk -v prefix="$header_prefix" 'BEGIN {FS=OFS=","} { for (i=1; i<=NF; i++) $i= prefix $i } 1')
    # Replace spaces with commas to create CSV headers
    # column_headers_csv=$(echo "$column_headers_blast_sub" | tr ' ' ',')

    # Add modified headers to the output file
    echo "$column_headers_blast_sub" >temp_output_file.csv
    # Append data from input csv but exclude the headers
    # This will allow you to replace them with the modified headers
    tail -n +2 "$input_csv" >>temp_output_file.csv
    # Rename temp file to the specified output csv file
    mv temp_output_file.csv "$header_output_csv"
    # rm "$input_csv"

    # Make sure the output now exists
    if [ -f "$header_output_csv" ]; then
        echo """
Csv file with modified headers was successfully created!
File is fould here: '$header_output_csv'
"""
    else
        echo "Error: Failed to create modified CSV file with headers"
        exit 1
    fi
}

# Usage statement
print_usage() {
    echo "Usage: $0 <input_csv> <header_output_csv> <header_prefix>"
    echo "Example: $0 input.csv output.csv prefix_"
}

## Check if the correct number of arguments is provided or if --help or -h option is passed
if [ "$#" -ne 3 ] || [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    print_usage
    exit 1
fi

## Call the functions with provided arguments
add_headers_prefix_csv "$1" "$2" "$3"
