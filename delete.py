import csv


def filter_csv(input_csv_path, output_csv_path, min_qcovs, min_pident):
    """Take blast output and filters by query coverage and percent identity

    Args:
        input_csv_path (csv): blast output as csv
        output_csv_path (csv): filtered blast output
        min_qcovs (float): minimum query coverage threshold
        min_pident (float): minimum percent identity threshold
    """
    # Initialize a list to store the filtered rows
    filtered_rows = []

    # Open the input CSV file
    with open(input_csv_path, "r") as input_csv_file:
        # Create a CSV reader
        csv_reader = csv.reader(input_csv_file, delimiter=",")

        # Read the header row
        column_headers_blast = next(csv_reader)

        # Append the header row to the filtered rows
        filtered_rows.append(column_headers_blast)

        # Find the indices of the columns you need
        qcovs_index = column_headers_blast.index("qcovs")
        pident_index = column_headers_blast.index("pident")

        # Iterate through rows in the input CSV and filter based on criteria
        for row in csv_reader:
            qcovs = float(row[qcovs_index])
            pident = float(row[pident_index])

            # Check if the row meets the filter criteria
            if qcovs >= min_qcovs and pident >= min_pident:
                filtered_rows.append(row)

    # Write the filtered rows to the output CSV file
    with open(output_csv_path, "w", newline="") as output_csv_file:
        csv_writer = csv.writer(output_csv_file, delimiter=",")
        csv_writer.writerows(filtered_rows)

    print("Filtered data has been written to", output_csv_path)


# Usage:
filter_csv(
    input_csv_path="auid.filtered.nifHDB_UCYNAoligos.csv",
    output_csv_path="filtered_output.csv",
    min_qcovs=70,
    min_pident=97,
)
