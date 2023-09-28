import pandas as pd


def merge_dataframes(
    df1, df2, output_path, columns_to_select, left_on, right_on, how
) -> pd.DataFrame:
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
        pd.DataFrame: The merged DataFrame.
    """
    # Usage example
    df1: pd.DataFrame = pd.read_csv(filepath_or_buffer=df1)
    df2: pd.DataFrame = pd.read_csv(
        filepath_or_buffer=df2,
        delimiter="\t",
    )
    print("Files imported")

    # Merge df1 and df2 based on specified parameters
    merged_df: pd.DataFrame = pd.merge(
        left=df1,
        right=df2[columns_to_select],
        left_on=left_on,
        right_on=right_on,
        how=how,
    )

    print(
        f"Merging files df1 by {left_on} and df2 by {right_on} through via {how} join"
    )

    # Save the merged DataFrame to the specified output path
    merged_df.to_csv(output_path, index=False)
    print(f"Merged file merged_df written out to: {output_path}")

    return merged_df


# Usage example

merge_params: dict[str, list[str] | str] = {
    "columns_to_select": ["AUID_updt", "AUID_OLD"],
    "left_on": "AUID_OLD",
    "right_on": "AUID_updt",
    # "how": "outer",
    "how": "left",
}

merged_df = merge_dataframes(
    df1="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/temp_merged_output_file.csv",
    df2="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/",
    # df2="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/nifHDB_old2new_auids_map.tsv",
    # output_path="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/annotationsAll_mergeAll_updatedAUIDs.csv",
    output_path="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/TEST_annotationsAll_mergeSub_updatedAUIDs.csv",
    **merge_params,
)
