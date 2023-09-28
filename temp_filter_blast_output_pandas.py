import importlib
from logging import Filterer

# Check if pandas is installed
try:
    importlib.import_module(name="pandas")
except ImportError:
    print("pandas is not installed. Installing...")
    import subprocess

    subprocess.check_call(args="pip3", bufsize="install", executable="pandas")

# Now import pandas
import pandas as pd

print("pandas imported!")


def filter_blast_output(
    input_blast_output_path, output_blast_path, min_qcovs=70, min_pident=97
) -> None:
    # Read the CSV file into a dataframe
    df: DataFrame = pd.read_csv(filepath_or_buffer=input_blast_output_path)

    # Filter rows based on criteria
    filtered_df: DataFrame = df[
        (df["qcovs"] >= min_qcovs) & (df["pident"] >= min_pident)
    ]

    # Write the filtered DataFrame to the output CSV file
    filtered_df.to_csv(output_blast_path, index=False)

    print(f"Filtered blast table has been writte to {output_blast_path}")


## Usage
filter_blast_output(
    input_blast_output_path="auid.filtered.nifHDB_UCYNAoligos.csv",
    output_blast_path="filtered_blast_UCYNA_oligoreps.csv",
    min_qcovs=70,
    min_pident=97,
)
# filter_blast_output( \
#     input_blast_output_path="auid.filtered.nifHDB_UCYNAoligos.csv",\
#     output_blast_path="filtered_blast_UCYNA_oligoreps.csv",\
#     min_qcovs=70,\
#     min_pident=97,\
# )
