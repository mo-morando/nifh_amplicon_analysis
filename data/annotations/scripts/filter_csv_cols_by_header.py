import os
import argparse
import pandas as pd

# Get the full path of the currently running script
script_path: str = __file__

# Extract the name of the script for later use
script_name: str = os.path.basename(p=script_path)

# Print the name of the script to signify start
print(f"{script_name} is currently running...")


# define function to read in and filter csv file
def filter_csv(
    input_csv,
    filtered_csv,
    headers_to_keep,
) -> None:
    """
    Reads in a CSV file, filters it based on specified headers, and writes the
    filtered data to a new CSV file.

    Args:
        input_csv (str): Path to the input CSV file.
        filtered_csv (str): Path to write the filtered CSV file to.
        headers_to_keep (list): Headers to retain in the filtered CSV file.

    Returns:
        None
    """

    # Define flag to track success in the script
    success = False

    # Start try block
    try:
        # Read in csv file
        csv: pd.DataFrame = pd.read_csv(filepath_or_buffer=input_csv)

        # Filter df based on specified input headers
        filt_csv = csv[headers_to_keep]

        # Write file out as a csv
        filt_csv.to_csv(filtered_csv, index=False)

        # Set success flag to true
        success = True

    # Define some errors
    except FileExistsError:
        print(f"Error: File not found at '{input_csv}'")
    except pd.errors.EmptyDataError:
        print(f"Error: DataFrame: '{input_csv}' is empty")
    except pd.errors.ParserError:
        print(
            f"""Error: There was an error parsing the input csv file:
'{input_csv}'"""
        )
    except KeyError as e:
        print(
            f"""Error: Column '{e.args[0]}' not found in the input DataFrame.
If error states 'None' is missing, an argument was not pass. Please see usage under --help and check your values.
"""
        )
    except Exception as e:
        print(f"Error: {e}")
    finally:
        # success bloc with print statements
        if success:
            print(f"Script '{script_name}' executed successfully!!")
        else:
            print(
                f"""Script '{script_name}' exited with an error. No filtration
was carried out and '{filtered_csv}' was not written! See logs for details...
"""
            )


def main() -> None:
    """
    Main function to parse command-line arguments and initiate CSV filtering.
    """
    # Define argparse
    parser = argparse.ArgumentParser(
        description="Filter CSV file based on specified headers supplied."
    )
    parser.add_argument(
        "--input_csv",
        type=str,
        help="Path to the input csv file.",
    )
    parser.add_argument(
        "--filtered_csv",
        type=str,
        help="Path to write the filtered csv file to.",
    )
    parser.add_argument(
        "--headers_to_keep",
        nargs="+",
        type=str,
        help="Headers to retain in the filtered csv file.",
    )

    # Parse the arguments
    args = parser.parse_args()

    # Call the function to filter the csv file
    filter_csv(
        input_csv=args.input_csv,
        filtered_csv=args.filtered_csv,
        headers_to_keep=args.headers_to_keep,
    )


if __name__ == "__main__":
    main()
