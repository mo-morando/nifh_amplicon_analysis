import pandas as pd


def assign_consensus_id(df, output_path, min_pid=97):
    # Read in CSV
    df = pd.read_csv(filepath_or_buffer=df)

    # Create a new column with the specified name and data type
    df.loc[df["Genomes879.pctId"] >= min_pid, "Genoes879.pid_flag"] = 1

    # Convert cluster and subcluster columns to string
    df["subcluster"] = df["subcluster"].astype(str)
    df["cluster"] = df["cluster"].apply(lambda x: str(int(x)) if not pd.isna(x) else x)

    # Fix subcluster
    # combine 1J and 1K

    ### Clean up oligo names
    # df["UCYNA_oligos"] = df["UCYNA_oligos"].replace(replacement_dict)
    # df["UCYNA_oligos"] = "UCYN-" + df["UCYNA_oligos"].str.split("_").str[1]

    # Create new column "MarineDiazo.ID". append subject to description separated by ";" and place value here
    df["MarineDiazo.ID"] = ""
    df["MarineDiazo.ID"] = (
        df["MarineDiazo.description"] + ";" + df["MarineDiazo.subject"]
    )

    # Initialize new column with default value "unknown"
    df["consensus_id"] = "unknown"

    # Define the condition for updating "consensus_id"
    condition_pid_flag = df["Genoes879.pid_flag"] == 1
    # condition_cluster3_4_flag = df["cluster"] == "3" or df["cluster"] == "4"

    # Create column
    df["consensus_id"] = (
        # df["UCYNA_oligos"]
        df["MarineDiazo.ID"]
        # .fillna(df["MarineDiazo.ID"])
        # .fillna(df["MarineDiazo.ID"])
        .fillna(df["Genomes879.id"].where(condition_pid_flag))
        # .fillna("unknown" + df["cluster"].where(condition_cluster3_4_flag).astype(str))
        # .fillna(df["subcluster"])
        .fillna(
            df.apply(
                lambda row: "unknown" + row["cluster"]
                if (row["cluster"] == "3" or row["cluster"] == "4")
                else "unknown" + row["subcluster"],
                axis=1,
            )
            # df.apply(
            #     lambda row: "unknown" + str(int(row["cluster"]))
            #     if (row["cluster"] == 3.0 or row["cluster"] == 4.0)
            #     else "unknown" + str(row["subcluster"]),
            #     axis=1,
            # )
        )
        # .fillna("unknown" + df["subcluster"])
    )

    # Write out file as CSV
    df.to_csv(output_path, index=False)
    print(f"File written out as '{output_path}'")


# ### Usage
# assign_consensus_id(
#     # df=data/annotations/nifhDB/temp_annotation_merge_allAUID.csv
#     df="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/temp_annotation_merge_allAUID.csv",
#     output_path="test_consensus_tax_annotation.csv",
# )
### Usage
# assign_consensus_id(
#     # df=data/annotations/nifhDB/temp_annotation_merge_allAUID.csv
#     df="/Users/mo/Projects/nifH_amp_project/myWork/analysis/files/annoNifHDB_updt.csv",
#     output_path="test_annoNifHDB_updt_consensus_tax_annotation.csv",
# )
### Usage
assign_consensus_id(
    # df=data/annotations/nifhDB/temp_annotation_merge_allAUID.csv
    df="annotationsall_updatedauids_nifhdb.csv",
    output_path="annotationsall_updatedauids_nifhdb_consensus_tax_annotation.csv",
)
