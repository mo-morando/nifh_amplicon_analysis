#!/bin/bash

# Define function
## to filter csv file by headers
filter_csv() {
    input_csv="$1"
    filtered_csv="$2"
    headers_to_keep="$3"

    # Some checks

    # Check if csv file exists
    if [ ! -f "$input_csv" ]; then
        echo "Error: csv file '$input_csv' not found"
        return 1
    fi

    # Were headers provided?
    if [ -z "$headers_to_keep" ]; then
        echo "Error: No '$headers_to_keep' provided"
        return 1
    fi

    # # Extract headers from input CSV
    # headers=$(head -n 1 "$input_csv")

    # Execute filteration via awk
    # awk -v cols="$headers_to_keep" 'BEGIN { FS=OFS","; split(cols, arr, " "); }
    # NR==1 {
    #     for (i=1; i<=NF; i++) {
    #         if (i == arr) {
    #             printf "%s%s", $i, (i<(NF) ? OFS : ORS);
    #         }
    #     }
    # }' "$input_csv" >"$filtered_csv"
    # awk -v cols="$headers_to_keep" 'BEGIN { FS=OFS=","; print "Cols before split: " cols; split(cols, arr, ","); print "Arr1: " arr[1]; print "arr[3]"; }

    ######

    # awk -v cols="$headers_to_keep" 'BEGIN { FS=OFS=","; split(cols, arr, ","); }
    # NR==1 {
    #     for (i=1; i<=NF; i++) {
    #         for (j=1; j<=length(arr); j++)
    #             if ($i == arr[j]) {
    #                 printf "%s%s", $i, (j<(length(arr)) ? OFS : ORS);
    #         }
    #     }
    # }' "$input_csv" >"$filtered_csv"

    #####

    awk -v headers="$headers_to_keep" 'BEGIN { FS=OFS=","; split(headers, arr, ","); }
    NR==1 {
        for (i=1; i<=NF; i++) {
            for (j=1; j<=length(arr); j++) {
                if ($i == arr[j]) {
                    printf "%s%s", $i, (j<(length(arr)) ? OFS : ORS);
                }
            }
        }
    }' "$input_csv" >"$filtered_csv"

    #    awk -v cols="$headers_to_keep" '
    # BEGIN { FS=OFS=","; split(cols, arr, ","); }
    # {
    #     if (NR == 1) {
    #         for (i=1; i<=NF; i++) {
    #             if ($i in arr) {
    #                 header_index[$i] = i;
    #             }
    #         }
    #         print "Header indices:"
    #         for (hdr in header_index) {
    #             print hdr, header_index[hdr]
    #         }
    #     } else {
    #         for (i in header_index) {
    #             printf "%s%s", $(header_index[i]), (i < length(header_index) ? OFS : ORS);
    #         }
    #     }
    # }' "$input_csv" >"$filtered_csv"

    # Make sure output file exits
    if [ -f "$filtered_csv" ]; then
        echo """
Csv file with modified headers was successfully created!
File is fould here: '$filtered_csv'
"""
    else
        echo "Error: Failed to create filtered CSV file with headers"
        exit 1
    fi
}

# Usage statement
print_usage() {
    echo "Usage: $0 <input_csv> <filtered_csv> <headers_to_keep>"
    echo "Example usages: $0 input.csv output.csv sseqid length pident evalue"
}

# Call function:
filter_csv "$1" "$2" "$3"
