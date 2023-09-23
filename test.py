from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


def perform_blast_search(
    input_fasta_file,
    blast_database,
    identity_threshold=97,
    blast_type="blastn",
    nucleotide_type="nucl",
):
    # Initialize a list to store results for each query
    results = []

    # Read the input FASTA file
    with open(input_fasta_file, "r") as fasta_file:
        queries = fasta_file.read().split(">")[1:]

    for query in queries:
        query_lines = query.split("\n")
        query_id = query_lines[0]
        query_sequence = "\n".join(query_lines[1:])

        # Perform BLAST search for each query sequence
        result_handle = NCBIWWW.qblast(
            blast_type,
            blast_database,
            query_sequence,
            megablast=True,
            expect=10.0,
            format_type="XML",
            matrix_name=nucleotide_type,
        )

        # Parse the BLAST results
        blast_records = NCBIXML.parse(result_handle)

        # Initialize variables to store top hit information for this query
        top_hit_auid = None
        top_hit_pid = 0
        top_hit_description = None
        top_hit_evalue = float("inf")

        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    percent_identity = (hsp.identities / hsp.align_length) * 100

                    if (
                        percent_identity >= identity_threshold
                        and hsp.expect < top_hit_evalue
                    ):
                        top_hit_auid = blast_record.query_id
                        top_hit_pid = percent_identity
                        top_hit_description = alignment.title
                        top_hit_evalue = hsp.expect

        # Append the top hit information for this query to the results list
        results.append((top_hit_auid, top_hit_pid, top_hit_description, top_hit_evalue))

        # Close the result handle for this query
        result_handle.close()

    # Return a list of results, one for each query
    return results


# perform_blast_search()
# Specify the input FASTA file and other parameters
input_fasta_file = "/Users/mo/Projects/nifH_amp_project/myWork/test.fasta"
blast_database = "/Users/mo/Library/CloudStorage/GoogleDrive-mmorando@ucsc.edu/My Drive/bioinfo/databases/UCYN_A_oligos/UCYN-A_oligoreps.fasta"
identity_threshold = 97
blast_type = "blastn"
nucleotide_type = "nucl"

# Perform the BLAST search
results = perform_blast_search(
    input_fasta_file, blast_database, identity_threshold, blast_type, nucleotide_type
)

# Print the results for each query sequence
for i, (auid, pid, description, evalue) in enumerate(results, start=1):
    print(f"Query {i}:")
    print(f"AUID: {auid}")
    print(f"Percent Identity: {pid}%")
    print(f"Top Hit Description: {description}")
    print(f"E-Value: {evalue}")
    print()
