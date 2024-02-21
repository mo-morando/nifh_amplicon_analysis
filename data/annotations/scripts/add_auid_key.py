import os
import argparse
import pandas as pd

# Get the full path to the currently running script
script_path: str = __file__

# Extract the name of the script from the full path
script_name: str = os.path.basename(script_path)

# print(f"The name of the currently running script is:", script_name)
print(f"{script_name} is currently running...")


# Define function
def add_auid_key(
    annotation_table,
    output_path,
    left_on,
    right_on,
    how,
    auid_key=None,
    columns_to_select=None,
):
    success = False  # Flag to track success

    # Set default value for columns_to_select if not provided
    # if columns_to_select is None:
    #     columns_to_select = ["AUID_updt", "AUID_OLD"]

    # # Set default value for columns_to_select if not provided
    # if auid_key is None:
    #     auid_key = "nifHDB_old2new_auids_map.tsv"

    try:
        # Read in annotation_table csv
        annotation_table: pd.DataFrame = pd.read_csv(
            filepath_or_buffer=annotation_table
        )
        # except FileNotFoundError:
        #     print(f"Error: File not found at '{annotation_table}'!")
        # Read in auid_key csv
        auid_key: pd.DataFrame = pd.read_csv(
            filepath_or_buffer=auid_key, delimiter="\t"
        )
        # except FileNotFoundError:
        #     print(f"Error: File not found at '{auid_key}'!")
        # Merge annotation_table and auid_key based on specified parameters

        # print(annotation_table.columns)
        # print(auid_key.columns)

        merged_df: pd.DataFrame = pd.merge(
            left=annotation_table,
            right=auid_key[columns_to_select],
            # right=auid_key,
            left_on=left_on,
            right_on=right_on,
            how=how,
        )

        # Provide feeback on mergin operation
        print(
            f"""
Merging files annotation_table by {left_on} and auid_key by {right_on} via a {how} join
"""
        )

        # Save the merged DataFrame to the specified output path
        merged_df.to_csv(output_path, index=False)
        print(
            f"""
Merged file merged_df written out to:
{output_path}
"""
        )

        success = True  # Set success flag to True

    # Define some errors
    except pd.errors.EmptyDataError:
        print("Error: One of the DataFrames is empty.")
    except pd.errors.ParserError:
        print("Error: There was an error parsing the CSV file.")
    except KeyError as e:
        print(
            f"""Error: Column '{e.args[0]}' not found in one of the DataFrames.
If error states 'None' is missing, an arugment was not pass. Please see usage under --help and check your values.
"""
        )
    except Exception as e:
        print(f"Error: {e}")
        print("Merging failed :(")

    # Success block with print statements
    if success:
        print(f"Script '{script_name}' executed successfully!!")
        # return merged_df
    else:
        print(
            f"""Script '{script_name}' exited with an error. Merged failed and '{output_path}' was not written! See logs for details :( """
        )
        # return None


## Define maine script to call from command line
def main():
    """
    Main function to handle command-line arguments and call the filtering function.

    Returns:
    - None
    """
    parser = argparse.ArgumentParser(
        description="Merge two DataFrames based on specified columns and merge parameters."
    )
    parser.add_argument("--left_csv", help="Path to left side csv")
    parser.add_argument(
        "--right_csv",
        default="../../databases/auid_key/nifHDB_old2new_auids_map.tsv",
        help="Path to right side csv",
    )
    parser.add_argument("--merged_csv", help="Path to merged csv")
    parser.add_argument(
        "--columns_to_select",
        nargs="+",
        default=["AUID_updt", "AUID_OLD"],
        # default="AUID_updt, AUID_OLD",
        help="Columns to select and keep from right side csv",
    )
    parser.add_argument(
        "--left_on",
        type=str,
        help="Column as the left join key to merge left side csv",
    )
    parser.add_argument(
        "--right_on",
        type=str,
        help="Column as the right join key to merge right side csv",
    )
    parser.add_argument(
        "--how",
        type=str,
        choices=["left", "right", "inner", "outer"],
        help="The type of merge to perform",
    )

    args = parser.parse_args()

    # Ensure required arguments are provided

    if args.left_csv is None or args.merged_csv is None:
        parser.error(
            "Please provide paths for both left input file as well as an output path for the new csv file. :)"
        )

    add_auid_key(
        annotation_table=args.left_csv,
        auid_key=args.right_csv,
        output_path=args.merged_csv,
        columns_to_select=args.columns_to_select,
        left_on=args.left_on,
        right_on=args.right_on,
        how=args.how,
    )


if __name__ == "__main__":
    # try:
    main()
