import argparse
import subprocess

# import sys


def create_blast_database(
    input_fasta_file,
    database_type,
    blast_database,
    # blast_database_dir,
):
    # Define makeblastdb executable
    makeblastdb_executable = "makeblastdb"
    try:
        subprocess.run(
            [
                makeblastdb_executable,
                "-in",
                input_fasta_file,
                "-dbtype",
                database_type,
                "-out",
                blast_database,
                # "-title",
                # blast_database_dir,
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        print(f"BLAST database: '{blast_database}' was created successfully!")
        return True  # This indicates success!
    except subprocess.CalledProcessError as e:
        print(f"Error creating BLAST database '{blast_database}': {e}")
        return False  # This means it didn't work :(


def parse_arguments():
    """
    Parse command-line (CL) arguments using argparse

    Returns:
        argparse.Namespace: Parsed arugments
    """
    parser = argparse.ArgumentParser(
        description="Create BLAST database from a FASTA file."
    )
    parser.add_argument("--input_fasta_file", help="Path to the input FASTA file.")
    parser.add_argument("--database_type", help="Type of database (nucl or prot).")
    parser.add_argument("--blast_database", help="Path to the output BLAST database")
    # parser.add_argument(
    #     "--blast_database_dir", help="Path to the output BLAST database directory"
    # )
    return parser.parse_args()


# * Used to ensure that the following code is only executed when run directly as # a script and not when its imported as a module
if __name__ == "__main__":
    # Calls this function to obtain the CL arguments
    args = parse_arguments()

    ## Calls the `create_blast_database` function with the parsed CL arguments
    success: bool = create_blast_database(
        input_fasta_file=args.input_fasta_file,
        database_type=args.database_type,
        blast_database=args.blast_database,
        # blast_database_dir=args.blast_database_dir,
    )

    if success:
        print("Database creation successful!")
    else:
        print("Database creation failed. Oh no!")
