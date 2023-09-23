import concurrent.futures
import subprocess

# def creat_blast_database(input_fasta_file, database_type, blast_database)
#     # Specify the full path to blastn executable
#     blastn_executable = "blast"

#     #Convert input fasta to blast_database
#     makeblastdb -in input_fasta_file -dbtype dbtype -out blast_database


def perform_blast_search(
    query_sequence, blast_database, identity_threshold=97, blast_type="blastn"
):
    # Specify the full path to blastn executable
    blastn_executable = "blastn"

    # Create the BLAST command using the full path
    blast_cmd = [
        blastn_executable,
        "-db",
        blast_database,
        "-query",
        "/dev/stdin",
        "-outfmt",
        "6",
        "-max_target_seqs",
        "1",
        "-evalue",
        "10",
    ]

    # Rest of the function remains the same...

    # Run BLAST and capture the output
    blast_process = subprocess.run(
        blast_cmd,
        input=query_sequence,  # Encode query sequence as bytes
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,  # Interpret output as text (for Python 3.7+)
    )

    # Check for errors
    if blast_process.returncode != 0:
        print(f"Error running BLAST:")
        print(blast_process.stderr)
    else:
        # Parse the BLAST results (assuming tabular output format)
        blast_result = blast_process.stdout.strip().split("\n")
        top_hit_data = blast_result[0].split("\t")

        # Extract relevant information
        top_hit_auid = query_sequence.split("\n")[0]
        top_hit_pid = float(top_hit_data[2])  # Percent identity
        top_hit_description = top_hit_data[1]  # Hit description
        top_hit_evalue = float(top_hit_data[10])  # E-value

        # Return the top hit information
        return top_hit_auid, top_hit_pid, top_hit_description, top_hit_evalue


def main(
    input_fasta_file,
    blast_database,
    identity_threshold=97,
    blast_type="blastn",
    num_threads=4,
):
    # Read the input FASTA file
    with open(input_fasta_file, "r") as fasta_file:
        queries = fasta_file.read().split(">")[1:]

    results = []

    # Create a thread pool with concurrent.futures
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Submit BLAST jobs in parallel
        futures = [
            executor.submit(
                perform_blast_search,
                query,
                blast_database,
                identity_threshold,
                blast_type,
            )
            for query in queries
        ]

        # Collect results as they become available
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                results.append(result)

    # Print or process the results
    for i, (auid, pid, description, evalue) in enumerate(results, start=1):
        print(f"Query {i}:")
        print(f"AUID: {auid}")
        print(f"Percent Identity: {pid}%")
        print(f"Top Hit Description: {description}")
        print(f"E-Value: {evalue}")
        print()


##
# Example usage
input_fasta_file = "/Users/mo/Projects/nifH_amp_project/myWork/test.fasta"
# blast_database = "./UCYN_A_oligos/UCYN-A_oligoreps.fasta"
blast_database = (
    "/Users/mo/Projects/nifH_amp_project/myWork/data/databases/UCYN-A_oligoreps_db"
)
identity_threshold = 97
blast_type = "blastn"
num_threads = 4

main(input_fasta_file, blast_database, identity_threshold, blast_type, num_threads)
