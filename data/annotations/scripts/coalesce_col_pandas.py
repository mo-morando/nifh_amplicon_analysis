import pandas as pd


def coalesce_columns(file_path):
    """
    Coalesce columns ending with '_x' or '_y' into single columns with '_x' suffix.
    Fill missing values in '_x' with values from '_y' when they exist.

    Args:
    file_path (str): The file path to the CSV file.

    Returns:
    None. The edited DataFrame is directly overwritten to the input file.
    """
    try:
        # Read in the CSV file
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        print("File not found.")
        return
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return

    try:
        # Iterate over columns
        for column in df.columns:
            # Check if column ends with '_x' or '_y'
            if column.endswith("_x"):
                # Get the base column name without suffix
                base_column_name = column[:-2]
                # Fill missing values in '_x' with values from '_y' when they exist
                df[column].fillna(df[f"{base_column_name}_y"], inplace=True)
                # Coalesce the conflicting columns into a single column with '_x' suffix
                df.rename(columns={column: base_column_name}, inplace=True)
                # Drop the corresponding '_y' column
                df.drop(columns=[f"{base_column_name}_y"], inplace=True)
    except Exception as e:
        print(f"An error occurred: {e}")
        return

    # Write out the edited DataFrame to the same CSV file
    try:
        df.to_csv(file_path, index=False)
        print("File successfully overwritten.")
    except Exception as e:
        print(f"An error occurred while writing the file: {e}")


# Example usage:
coalesce_columns(
    "/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/merged_database/arb_genome879_merge_copy.csv"
)
