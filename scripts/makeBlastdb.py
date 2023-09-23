import subprocess


def create_blast_database(input_fasta_file, database_type, blast_database):
    # Specify the full path to makeblastdb executable
    makeblastdb_executable = "makeblastdb"

    # Create the BLAST database
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
            ],
            check=True,  # Check for errors
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,  # Interpret output as text (for Python 3.7+)
        )
        print(f"BLAST database '{blast_database}' created successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error creating BLAST database: {e}")


# Example usage to create a BLAST database from a FASTA file
input_fasta_file = "/Users/mo/Projects/nifH_amp_project/myWork/UCYN-A_oligoreps.fasta"
# input_fasta_file = "UCYN_A_oligos/UCYN-A_oligoreps.fasta"
# input_fasta_file = "/Users/mo/Library/CloudStorage/GoogleDrive-mmorando@ucsc.edu/My Drive/bioinfo/databases/UCYN_A_oligos/UCYN-A_oligoreps.fasta"
# input_fasta_file = "../../../Library/CloudStorage/GoogleDrive-mmorando@ucsc.edu/My Drive/bioinfo/databases/UCYN_A_oligos/UCYN-A_oligoreps.fasta"
database_type = "nucl"  # Change to "prot" if it's a protein database
blast_database = (
    "/Users/mo/Projects/nifH_amp_project/myWork/data/databases/UCYN-A_oligoreps_db"
)

create_blast_database(input_fasta_file, database_type, blast_database)
