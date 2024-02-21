import os
import argparse
import pandas as pd

# Get the full path to the currently running script
script_path = __file__

# Extract the name of the script from the full path
script_name = os.path.basename(script_path)

print("The name of the currently running script is:", script_name)


def filter_blast_output(
    input_blast_csv, filtered_blast_output, min_qcovs=70, min_pident=97
) -> None:
    """
    Filter blast output based on specified criteria.

    Parameters:
    - input_blast_csv (str): Path to the input blast output CSV file.
    - filtered_blast_output (str): Path to the filtered blast output CSV file.
    - min_qcovs (float): Minimum qcovs value threshold (default: 70).
    - min_pident (float): Minimum pident value threshold (default: 97).

    Returns:
    - None
    """

    success = False  # Flag to track success

    try:
        # Read the CSV file into a dataframe
        df: DataFrame = pd.read_csv(filepath_or_buffer=input_blast_csv)

        # Filter rows based on criteria
        filtered_df: DataFrame = df[
            (df["qcovs"] >= min_qcovs) & (df["pident"] >= min_pident)
        ]
        print("Filtering completed successfully!")

        # Write the filtered DataFrame to the output CSV file
        filtered_df.to_csv(filtered_blast_output, index=False)
        # Print statement for feedback
        print(
            f"""
Filtered blast table has been written to:
{filtered_blast_output}
"""
        )

        success = True  # Set success flag to True

    # Define some errors:
    except FileNotFoundError:
        print("Error: File not found.")
    except pd.errors.EmptyDataError:
        print("Error: One of the DataFrames is empty.")
    except pd.errors.ParserError:
        print("Error: There was an error parsing the CSV file.")
    except Exception as e:
        print(f"Error: {e}")
        print("Filtering failed :(")

    # Define succes block with print statments
    if success:
        print(f"Script '{script_name}' executed successfully!!")
    else:
        print(
            f"""Script '{script_name}' exited with an error. Filtation failed 
and '{filtered_blast_output}' was not written. Please see log for details :("""
        )


## Define main function to call from command line
def main():
    """
    Main function to handle command-line arguments and call the filtering function.

    Returns:
    - None
    """
    parser = argparse.ArgumentParser(
        description="Filter blast output based on specific criteria"
    )
    parser.add_argument("--input_blast_csv", help="Path to blast file to parse")
    parser.add_argument(
        "--filtered_blast_output", help="Path to filter blast file output"
    )
    parser.add_argument(
        "--min_qcovs",
        type=float,
        default=70.0,
        help="Minimum qcovs value threshold (default: 70)",
    )
    parser.add_argument(
        "--min_pident",
        type=float,
        default=97.0,
        help="Minimum pident value threshold (default: 97)",
    )

    args = parser.parse_args()

    if args.input_blast_csv is None or args.filtered_blast_output is None:
        parser.error(
            "Please provide paths for both input_blast_csv and filtered_blast_output files. :)"
        )

    filter_blast_output(
        input_blast_csv=args.input_blast_csv,
        filtered_blast_output=args.filtered_blast_output,
        min_qcovs=args.min_qcovs,
        min_pident=args.min_pident,
    )


if __name__ == "__main__":
    # try:
    main()
#     print("Script executed successfully.")
# except SystemExit:
#     print("Script exited with an error.")
# except Exception as e:
#     print(f"Unexpected error: {e}")
