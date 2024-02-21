import os
import argparse
import pandas as pd

# Get the full path to the currently running script
script_path: str = __file__

# Extract the name of the script from the full path
script_name: str = os.path.basename(script_path)

# print(f"The name of the currently running script is:", script_name)
print(f"{script_name} is currently running...")


# Define function to merge cvs files
def merge_dataframes(
    df1,
    df2,
    output_path,
    # columns_to_select,
    left_on,
    right_on,
    how,
) -> None:
    """
    Merge two DataFrames based on specified columns and merge parameters.

    Args:
        df1 (pd.DataFrame): The first DataFrame.
        df2 (pd.DataFrame): The second DataFrame.
        output_path (str): The path to save the merged DataFrame as a CSV file.
        columns_to_select (list): The columns to select from df2 for merging.
        left_on (str): The column in df1 to use as the left join key.
        right_on (str): The column in df2 to use as the right join key.
        how (str): The type of merge to perform (e.g., 'left', 'inner', 'outer').

    Returns:
        None

    Notes:
        This function merges the specified DataFrames and writes the result to a CSV file at the specified output path. It does not return any value.
    """

    # Initialize merged_df with None
    # merged_df = None  # Initialize merged_df to None

    success = False  # Flag to track success

    try:
        # Read in df1 csv
        df1: pd.DataFrame = pd.read_csv(filepath_or_buffer=df1)
        # except FileNotFoundError:
        #     print(f"Error: File not found at '{df1}'!")
        # Read in df2 csv
        df2: pd.DataFrame = pd.read_csv(filepath_or_buffer=df2)
        # except FileNotFoundError:
        #     print(f"Error: File not found at '{df2}'!")
        # Merge df1 and df2 based on specified parameters
        merged_df: pd.DataFrame = pd.merge(
            left=df1,
            # right=df2[columns_to_select],
            right=df2,
            left_on=left_on,
            right_on=right_on,
            how=how,
        )

        # Provide feeback on mergin operation
        print(
            f"""
Merging files df1 by {left_on} and df2 by {right_on} through via {how} join
"""
        )

        # * Coalesce conflicting columns due to differences in AUIDs that were
        # annotated by each step

        # Iterate over each column
        for column in merged_df.columns:
            # Check for columns ending in '_x' or '_y'
            if column.endswith("_x"):
                # Get the base column name without the suffix
                base_column_name = column[:-2]
                # Fill missing values in '_x' with values from '_y' when they
                # exist
                merged_df[column].fillna(
                    value=merged_df[f"{base_column_name}_y"], inplace=True
                )
                # Coalesce the conflicting columns into a single column with
                # '_x' suffix
                merged_df.rename(columns={column: base_column_name}, inplace=True)
                merged_df.drop(columns=[f"{base_column_name}_y"], inplace=True)

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
    parser.add_argument("--right_csv", help="Path to right side csv")
    parser.add_argument("--merged_csv", help="Path to merged csv")
    # parser.add_argument(
    #     "--columns_to_select",
    #     nargs="+",
    #     help="Columns to select and keep from right side csv",
    # )
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

    if args.left_csv is None or args.right_csv is None or args.merged_csv is None:
        parser.error(
            "Please provide paths for both left and right input files as well as merged csv output path. :)"
        )

    merge_dataframes(
        df1=args.left_csv,
        df2=args.right_csv,
        output_path=args.merged_csv,
        # columns_to_select=args.columns_to_select,
        left_on=args.left_on,
        right_on=args.right_on,
        how=args.how,
    )


if __name__ == "__main__":
    # try:
    main()
    #     print("Script executed successfully.")
    # except SystemExit:
    #     print("Script exited with an error.")
    # except Exception as e:
    #     print(f"Unexpected error: {e}")

# # Usage example

# python3 ../scripts/temp_parsingAUIDs_function.py --left_csv --right_csv UCYNA_oligos/filt_cols_filt_blast_output --merged_csv --columns_to_select --left_on --right_on --how

# merge_params: dict[str, list[str] | str] = {
#     "columns_to_select": ["AUID_updt", "AUID_OLD"],
#     "left_on": "AUID_OLD",
#     "right_on": "AUID_updt",
#     # "how": "outer",
#     "how": "left",
# }


# merged_df = merge_dataframes(
#     df1="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/jonathan/FOR_MO_ESSD/auids.annot.tsv",
#     # df2="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/",
#     df2="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/nifHDB_old2new_auids_map.tsv",
#     # output_path="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/annotationsAll_mergeAll_updatedAUIDs.csv",
#     output_path="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/annotationsall_updatedauids_nifhdb.csv",
#     **merge_params,
# )


# merged_df = merge_dataframes(
#     df1="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/temp_merged_output_file.csv",
#     df2="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/",
#     # df2="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/nifHDB_old2new_auids_map.tsv",
#     # output_path="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/annotationsAll_mergeAll_updatedAUIDs.csv",
#     output_path="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/TEST_annotationsAll_mergeSub_updatedAUIDs.csv",
#     **merge_params,
# )
