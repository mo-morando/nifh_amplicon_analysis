import pandas as pd


# Read in dataframes
df1 = pd.read_csv(
    # "../data/annotations/nifhDB/auids.annot_nifHDB_newAUID.tsv",
    "/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/auids.annot_nifHDB_newAUID.tsv",
    delimiter="\t",
)
df2 = pd.read_csv(
    "/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/nifHDB_old2new_auids_map.tsv",
    delimiter="\t",
)

print("Files imported")
print(df1.columns)
print(df2.columns)

# Merge df1 and df2 on the 'AUID' column to get the corresponding 'AUID_OLD' values
merged_df = pd.merge(
    df1,
    df2[["AUID_updt", "AUID_OLD"]],
    left_on="AUID",
    right_on="AUID_updt",
    how="left",
)

# # Replace the 'AUID' column in df1 with the 'AUID_OLD' values from merged_df
# df1["AUID"] = merged_df["AUID_OLD"]

# # Drop the 'AUID_updt' and 'AUID_OLD' columns from merged_df
# merged_df = merged_df.drop(["AUID_updt", "AUID_OLD"], axis=1)

# Write the modified df1 to a new CSV file
merged_df.to_csv("temp_merged_output_file.csv", index=False)
