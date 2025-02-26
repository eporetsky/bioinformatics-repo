import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def extract_sequences_from_fasta(diffbind_bed, genome_fasta, output_fasta):
    """
    Extract sequences from a FASTA file based on genomic coordinates in a BED file.

    Parameters:
    - diffbind_bed: Path to the BED file containing genomic coordinates (DiffBind output).
    - genome_fasta: Path to the FASTA file of the genome.
    - output_fasta: Path to the output FASTA file to save extracted sequences.
    """
    # Load the genome FASTA file into a dictionary for quick access
    print("Loading genome FASTA file...")
    genome_dict = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    # Open the BED file and parse the coordinates
    print("Reading BED file and extracting sequences...")
    extracted_sequences = []
    with open(diffbind_bed, "r") as bed_file:
        for line in bed_file:
            if line.startswith("#") or line.strip() == "":
                continue  # Skip comments or empty lines
            
            # Parse BED fields: chrom, start, end, name (optional), score (optional), strand (optional)
            fields = line.strip().split("\t")
            chrom = fields[0]
            start = int(fields[1])  # BED is 0-based, inclusive start
            end = int(fields[2])    # BED is 0-based, exclusive end
            name = fields[3] if len(fields) > 3 else f"{chrom}:{start}-{end}"
            strand = fields[5] if len(fields) > 5 else "+"

            # Extract the sequence from the genome
            if chrom in genome_dict:
                seq = genome_dict[chrom].seq[start:end]  # Extract sequence
                if strand == "-":  # Reverse complement if on the negative strand
                    seq = seq.reverse_complement()
                
                # Create a SeqRecord for the extracted sequence
                record = SeqRecord(seq, id=name, description=f"{chrom}:{start}-{end}({strand})")
                extracted_sequences.append(record)
            else:
                print(f"Warning: Chromosome {chrom} not found in genome FASTA.")

    # Write the extracted sequences to the output FASTA file
    print(f"Writing {len(extracted_sequences)} sequences to {output_fasta}...")
    with open(output_fasta, "w") as output_file:
        SeqIO.write(extracted_sequences, output_file, "fasta")

    print("Sequence extraction complete!")

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Extract sequences from a genome FASTA file based on BED coordinates.")
    parser.add_argument("-b", "--diffbind_bed", required=True, help="Path to the BED file containing genomic coordinates (DiffBind output).")
    parser.add_argument("-f", "--genome_fasta", required=True, help="Path to the FASTA file of the genome.")
    parser.add_argument("-o", "--output_fasta", required=True, help="Path to the output FASTA file to save extracted sequences.")

    # Parse the arguments
    args = parser.parse_args()

    # Call the function with the provided arguments
    extract_sequences_from_fasta(args.diffbind_bed, args.genome_fasta, args.output_fasta)