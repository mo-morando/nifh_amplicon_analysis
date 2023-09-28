import string
import pandas as pd


def assign_consensus_id(df, output_path, min_pid=97):
    # Read in CSV
    df = pd.read_csv(filepath_or_buffer=df)

    # Create a new column with the specified name and data type
    df.loc[df["Genomes879.pctId"] >= min_pid, "Genoes879.pid_flag"] = 1

    # Convert cluster column to string
    # print(df["cluster"])
    # print(df["cluster"].dtype)
    df["cluster"] = (
        df["cluster"]
        .fillna("unresolved")
        .astype(dtype=str)
        .str.replace(pat=r"\.0*", repl="")
    )
    # print(df["cluster"])
    # print(df["cluster"].dtype)

    # print(df["cluster"].drop_duplicates())

    # Create new column "MarineDiazo.ID". append subject to description separated by ";" and place value here
    df["MarineDiazo.ID"] = ""
    df["MarineDiazo.ID"] = (
        df["MarineDiazo.description"] + ";" + df["MarineDiazo.subject"]
    )
    # Initialize new column with default value "unknown"
    df["consensus_id"] = "unknown"

    # Define the condition for updating "consensus_id"
    condition_pid_flag = df["Genoes879.pid_flag"] == 1

    # Conditionally fill "consensus_id" based on my criteria
    print(df.dtypes)

    # condition_cluster3_4_flag = (df["cluster"] == "3.0") | (df["cluster"] == "4.0")
    # condition_cluster3_4_flag = row["cluster"] == "3.0" or row["cluster"] == "4.0"
    # print(condition_cluster3_4_flag)

    df["consensus_id"] = (
        df["UCYNA_oligos"]
        .fillna(df["MarineDiazo.ID"])
        .fillna(df["Genomes879.id"].where(condition_pid_flag))
        # .fillna("unknown" + df["cluster"].where(condition_cluster3_4_flag).astype(str))
        # .fillna(df["subcluster"])
        .fillna(
            df.apply(
                lambda row: "unknown" + row["cluster"]
                if (row["cluster"] == "3.0" or row["cluster"] == "4.0")
                else "unknown" + str(row["subcluster"]),
                axis=1,
            )
        )
        # .fillna("unknown" + df["subcluster"])
    )

    print(df["consensus_id"].head())

    # Fix subcluster
    # combine 1J and 1K
    # remove letters after 3 and 4

    # p = df["consensus_id"].str.contains("unknown").drop_duplicates()
    # print(p)

    # Write out file as CSV
    df.to_csv(output_path)


### Usage
assign_consensus_id(
    # df=data/annotations/nifhDB/temp_annotation_merge_allAUID.csv
    df="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/temp_annotation_merge_allAUID.csv",
    output_path="test_consensus_tax_annotation.csv",
)
