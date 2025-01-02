import argparse
import glob
import os
import subprocess

import pandas as pd
from Bio import SeqIO

# Set up argument parsing
parser = argparse.ArgumentParser(
    description="Cluster sequences based on identity threshold."
)
parser.add_argument(
    "fasta_file", help="The path to the input FASTA file containing the sequences."
)
parser.add_argument(
    "db_name", help="The name (path) for the BLAST database to be created or used."
)
parser.add_argument(
    "output_file", help="Path to the file where clusters will be saved."
)

# Parse arguments
args = parser.parse_args()
# Step 1: Read sequences from a multi-FASTA file
fasta_file = args.fasta_file
sequences = list(SeqIO.parse(fasta_file, "fasta"))

# Initialize a pandas DataFrame to store percentage of identical matches (pident)
pident_matrix = pd.DataFrame(
    index=[seq.id for seq in sequences],
    columns=[seq.id for seq in sequences],
    data=None,
)

# BLAST parameters
db_name = args.db_name
makeblastdb_command = [
    "makeblastdb",
    "-in",
    fasta_file,
    "-dbtype",
    "nucl",
    "-out",
    db_name,
]

subprocess.run(makeblastdb_command)
blast_parameters = [
    "blastn",
    "-task",
    "blastn",
    "-db",
    db_name,
    "-outfmt",
    "6",
    "-max_target_seqs",
    "10000",
    "-max_hsps",
    "1",
    "-evalue",
    "10",
    "-reward",
    "1",
    "-penalty",
    "-2",
    "-gapopen",
    "1",
    "-gapextend",
    "1",
]

for i, query_seq in enumerate(sequences):
    # Save the query sequence to a temporary file
    query_file = f"temp_query_{db_name}_{i}.fasta"
    SeqIO.write(query_seq, query_file, "fasta")

    # Run BLASTn using subprocess and capture output
    blast_command = blast_parameters + ["-query", query_file]
    result = subprocess.run(blast_command, capture_output=True, text=True)
    # Process the BLAST output directly from memory
    for line in result.stdout.strip().split("\n"):
        parts = line.strip().split()
        query_id, subject_id, pident, align_length = (
            parts[0],
            parts[1],
            float(parts[2]),
            int(parts[3]),
        )
        query_len = len(query_seq.seq)

        # Ensure alignment length is greater than 80% of the query length
        if align_length >= 0.8 * query_len:
            pident_matrix.at[query_id, subject_id] = pident
        else:
            pident_matrix.at[query_id, subject_id] = (
                None  # Fill with None if not meeting criteria
            )

    # Cleanup: remove temporary query file
    os.remove(query_file)

# Write the matrix to a CSV file
output_matrix_file = f"pident_matrix_{db_name}.csv"
pident_matrix.to_csv(output_matrix_file)

# Create the pattern to match all files starting with db_name
pattern = f"{db_name}*"

# Find all files matching the pattern
files_to_delete = glob.glob(pattern)

# Loop through the files and delete them
for file_path in files_to_delete:
    try:
        os.remove(file_path)
        print(f"Deleted {file_path}")
    except Exception as e:
        print(f"Error deleting {file_path}: {e}")


print(f"Pairwise identity matrix saved to '{output_matrix_file}'.")

# Prepare a dictionary to hold the clusters
clusters = {}


# Function to find the cluster for a sequence
def find_cluster(seq_id):
    for cluster_id, members in clusters.items():
        if seq_id in members:
            return cluster_id
    return None


# Iterate over the matrix to cluster sequences
for seq_id in pident_matrix.columns:
    # Skip if the sequence is already in a cluster
    if find_cluster(seq_id) is not None:
        continue

    # Create a new cluster for this sequence
    cluster_id = len(clusters) + 1
    clusters[cluster_id] = [seq_id]

    # Check against all other sequences
    for other_seq_id in pident_matrix.columns:
        # Skip comparison with itself
        if seq_id == other_seq_id:
            continue

        # Check if the identity is above 90%
        if pident_matrix.at[seq_id, other_seq_id] > 90:
            # Add to the same cluster if not already in another cluster
            if find_cluster(other_seq_id) is None:
                clusters[cluster_id].append(other_seq_id)

# Output the clusters
with open(args.output_file, "w") as f:
    for cluster_id, members in clusters.items():
        f.write(f"Cluster {cluster_id}: {', '.join(members)}\n")

print(f"Clusters saved to '{args.output_file}'.")
