import pandas as pd

# merge dfs
df_reps = pd.read_csv(
    filepath_or_buffer="/Users/mo/Projects/nifH_amp_project/myWork/filtered_blast_UCYNA_oligoreps.csv"
)
# convert df_annotations into CVS
# df_annotations = pd.read_csv(
#     "/Users/mo/Library/CloudStorage/GoogleDrive-mmorando@ucsc.edu/My Drive/data/amplicon_review/all_studies/master_annotation/filtered_annotations/nifHDB/auids.annot_nifHDB_newAUID.tsv",
#     delimiter="\t",
# )
df_annotations: pd.DataFrame = pd.read_csv(
    # "../data/annotations/nifhDB/auids.annot_nifHDB_newAUID.tsv",
    filepath_or_buffer="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/temp_merged_output_file.csv",
    # delimiter="\t",
)

# Fix column IDs
df_reps: pd.DataFrame = df_reps.rename(
    columns={"qseqid": "AUID", "sseqid": "UCYNA_oligos"}
)
df_annotations: pd.DataFrame = df_annotations.rename(
    columns={"AUID": "AUID_NEW", "AUID_OLD": "AUID"}
)

# Only grab the columns I want
# subset columns of interest
df_reps = df_reps.loc[
    :, ["AUID", "UCYNA_oligos", "qcovs", "pident", "evalue", "bitscore"]
]
# Remove column I don't want
del df_annotations["AUID_updt"]

# Merge based on a common column
merged_df: pd.DataFrame = pd.merge(
    left=df_annotations,
    right=df_reps,
    on="AUID",
    # right_on="AUID",
    # left_on="AUID_OLD",
    how="outer",  # I ended up with the same results despite using outer or left
)

# Write out file to CSV
merged_df.to_csv(
    path_or_buf="/Users/mo/Projects/nifH_amp_project/myWork/data/annotations/nifhDB/temp_annotation_merge_allAUID.csv"
)
print("Merged file merged_df written out as: 'temp_annotation_merge_allAUID.csv'")

print("Done")
