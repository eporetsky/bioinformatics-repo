import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd

def split_and_save_non_redundant_sequences(input_folder, output_folder, output_tsv, sequences_per_file=500):
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Collect all non-redundant sequences and their original IDs
    non_redundant_sequences = {}
    all_gene_ids = []

    # Process each FASTA file in the input folder
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".fasta") or file_name.endswith(".fa"):
            file_path = os.path.join(input_folder, file_name)
            print(f"Processing file: {file_name}")
            genome_id = os.path.basename(file_name).replace(".fa", "")
            for record in SeqIO.parse(file_path, "fasta"):
                seq_str = str(record.seq).strip()  # Ensure no leading/trailing whitespace
                if seq_str and all(c.isalnum() or c in "-.*" for c in seq_str):  # Validate sequence
                    if seq_str not in non_redundant_sequences:
                        # Assign a new numbered ID to the sequence
                        new_id = f"seq_{len(non_redundant_sequences) + 1}"
                        non_redundant_sequences[seq_str] = new_id
                    # Save mapping of original gene ID to the reference ID
                    all_gene_ids.append({"geneID": record.id, "genome_id": genome_id,"sequence_id": non_redundant_sequences[seq_str]})
                else:
                    print(f"Skipped invalid sequence in {file_name}: {record.id}")

    # Save non-redundant sequences into multiple FASTA files
    sequences = list(non_redundant_sequences.items())
    for i in range(0, len(sequences), sequences_per_file):
        split_sequences = sequences[i:i + sequences_per_file]
        output_file = os.path.join(output_folder, f"split_{i // sequences_per_file + 1}.fasta")
        records = [SeqRecord(Seq(seq), id=seq_id, description="") for seq, seq_id in split_sequences]
        SeqIO.write(records, output_file, "fasta")
        print(f"Saved {len(records)} sequences to {output_file}")

    # Save mapping to a TSV file
    df = pd.DataFrame(all_gene_ids)
    df.to_csv(output_tsv, sep="\t", index=False)
    print(f"Mapping saved to {output_tsv}")

# Define paths and parameters
input_folder = "clean/"  # Folder containing input FASTA files
output_folder = "splits/"  # Folder to save split FASTA files
output_tsv = "sequence_mapping.tsv"  # TSV file to save gene ID mapping
split_and_save_non_redundant_sequences(input_folder, output_folder, output_tsv)
