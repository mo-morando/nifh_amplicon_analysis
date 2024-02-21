import os
import argparse
import pandas as pd

# Get the full path to teh currecctly running script
script_path: str = __file__

# Extract the name of the script for use later
script_name: str = os.path.basename(p=script_path)

print(f"{script_name} is currently running...")


# Define function to make consensus ID
def assign_consensus_id(
    annotation_table,
    output_table,
    min_pid_genomes879=97.0,
) -> None:
    """
    Assigns consensus taxonomy ID for an annotation file.

    Args:
        annotation_table (str): Path to the input annotation table.
        output_table (str): Path to the output annotation table with consensus ID added.
        min_pid_genomes879 (float, optional): Minimum threshold percentage identity to consider Genomes879.id in consensus ID. Defaults to 97.0.
    """

    # Flag to track success in script
    success = False

    try:
        # Read input CSV
        annotation_table: pd.DataFrame = pd.read_csv(
            filepath_or_buffer=annotation_table,
        )

        # Convert cluster and subcluster columns to string
        annotation_table["subcluster"] = annotation_table["subcluster"].astype(str)
        annotation_table["cluster"] = annotation_table["cluster"].apply(
            lambda x: str(object=int(x)) if not pd.isna(x) else x
        )  # -# converting to interger and then strings gets rid of decimals

        ### Clean up oligo names
        annotation_table["ucyna_oligos_sseqid"] = (
            "UCYN-" + annotation_table["ucyna_oligos_sseqid"].str.split("_").str[1]
        )

        ### Clean up ncd and cyano columns

        # Make dictionary of columns to rename
        new_columns_ncd_cyano: dict[str, str] = {
            "seqID": "ncd_cyano_subject",
            "description": "ncd_cyano_description",
        }

        # rename them
        annotation_table = annotation_table.rename(columns=new_columns_ncd_cyano)

        # remove unneeded column
        annotation_table = annotation_table.drop(columns="seqID2")

        # Create new column "ncd_cyano_id" by concatenating description to subject
        annotation_table["ncd_cyano_id"] = ""
        annotation_table["ncd_cyano_id"] = (
            annotation_table["ncd_cyano_description"]
            + ";"
            + annotation_table["ncd_cyano_subject"]
        )

        ### Clean up ncd and cyano columns

        # Rename column
        annotation_table = annotation_table.rename(
            columns={"genus_species": "genome879_tax"}
        )

        # Remove unneeded column
        annotation_table = annotation_table.drop(columns="nifHCluster")

        # - Get things together for the consensus id
        # Create a new column 'genome879_nifh_pid_flag' based on minimum PID
        annotation_table.loc[
            annotation_table["genome879_nifh_pident"] >= min_pid_genomes879,
            "genome879_nifh_pid_flag",
        ] = 1
        # annotation_table.loc[
        #     annotation_table["Genomes879.pctId"] >= min_pid_genomes879,
        #     "Genomes879.pid_flag",
        # ] = 1

        # Define a condition for updating "consensus_id"
        condition_pid_flag: Series[_bool] = (
            annotation_table["genome879_nifh_pid_flag"] == 1
        )

        # Initialize new column "consensus_id" with default value "unknown"
        annotation_table["consensus_id"] = "unknown"

        # _# Assign values to 'consensus_id' based on conditions
        annotation_table["consensus_id"] = (
            annotation_table["ucyna_oligos_sseqid"]
            .fillna(annotation_table["ncd_cyano_id"])
            .fillna(annotation_table["genome879_tax"].where(condition_pid_flag))
            # .fillna("unknown" + annotation_table["cluster"].where(condition_cluster3_4_flag).astype(str))
            # .fillna(annotation_table["subcluster"])
            .fillna(
                annotation_table.apply(
                    lambda row: "unknown" + row["cluster"]
                    if (row["cluster"] == "3" or row["cluster"] == "4")
                    else "unknown" + row["subcluster"],
                    axis=1,
                )
                # annotation_table.apply(
                #     lambda row: "unknown" + str(int(row["cluster"]))
                #     if (row["cluster"] == 3.0 or row["cluster"] == 4.0)
                #     else "unknown" + str(row["subcluster"]),
                #     axis=1,
                # )
            )
            # .fillna("unknown" + annotation_table["subcluster"])
        )

        # Write out file as CSV
        annotation_table.to_csv(output_table, index=False)
        print(f"Output file '{output_table}' written successfully.")
        # print(f"Script '{script_name}' executed successfully!!")

        # set success flag to True
        success = True

    # Define some errors
    except FileNotFoundError:
        print(f"Error: File not fount at '{annotation_table}'!")
    except pd.errors.EmptyDataError:
        print(f"Error: DataFrame: '{annotation_table}' is empty.")
    except pd.errors.ParserError:
        print(f"Error: There was an error parsing the csv file: '{annotation_table}'.")
    except KeyError as e:
        print(
            f"""Error: Column '{e.args[0]}' not found in the input DataFrame.
If error states "'None' is missing", an argument was not pass. Please see usage under --help and check your values.
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
                f"""Script '{script_name}' exited with an error. No consensus ID
was made and '{output_table}' was not written! See logs for details...
"""
            )


## Define main function to call from command line
def main() -> None:
    """
    Parses command-line arguments and calls assign_consensus_id function.
    """
    parser = argparse.ArgumentParser(
        description="Assign consensus taxonomy ID for annotation file"
    )
    parser.add_argument(
        "--annotation_table",
        help="Path to the input annotation table",
    )
    parser.add_argument(
        "--output_table",
        help="Path to output annotation table with consensus ID added",
    )
    parser.add_argument(
        "--min_pid_genomes879",
        type=float,
        default=97.0,
        help="Minimum threshold pid to consider Genomes879.id in consensus ID. If no value is give, 97.0 pid is used by default",
    )

    args: Namespace = parser.parse_args()

    assign_consensus_id(
        annotation_table=args.annotation_table,
        output_table=args.output_table,
        min_pid_genomes879=args.min_pid_genomes879,
    )


if __name__ == "__main__":
    main()
