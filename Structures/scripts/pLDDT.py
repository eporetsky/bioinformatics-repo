import os
import argparse
from concurrent.futures import ProcessPoolExecutor
from Bio import PDB
import multiprocessing

def calculate_average_plddt(pdb_file):
    """
    Calculate the average pLDDT score from a PDB file.

    Parameters:
    pdb_file (str): Path to the PDB file.

    Returns:
    float: The average pLDDT score.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    plddt_scores = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    plddt_scores.append(atom.bfactor)

    if plddt_scores:
        return sum(plddt_scores) / len(plddt_scores)
    return None

def process_single_file(pdb_file_path):
    """
    Process a single PDB file to calculate its average pLDDT score.

    Parameters:
    pdb_file_path (str): Path to the PDB file.

    Returns:
    tuple: Filename and its average pLDDT score.
    """
    avg_plddt = calculate_average_plddt(pdb_file_path)
    return os.path.basename(pdb_file_path), avg_plddt

def process_pdb_files(input_path, output_tsv, num_threads):
    """
    Process PDB files from a single file or a folder in parallel and save the average pLDDT scores to a TSV file.

    Parameters:
    input_path (str): Path to the input PDB file or folder containing PDB files.
    output_tsv (str): Path to the output TSV file.
    num_threads (int): Number of threads to use for parallel processing.
    """
    if os.path.isdir(input_path):
        pdb_files = [os.path.join(input_path, f) for f in os.listdir(input_path) if f.endswith('.pdb')]
    elif os.path.isfile(input_path) and input_path.endswith('.pdb'):
        pdb_files = [input_path]
    else:
        print("Error: Input must be a PDB file or a directory containing PDB files.")
        return

    results = []

    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        for pdb_file, avg_plddt in executor.map(process_single_file, pdb_files):
            if avg_plddt is not None:
                results.append((pdb_file, avg_plddt))
            else:
                print(f"Warning: No pLDDT scores found in {pdb_file}")

    with open(output_tsv, 'w') as out_file:
        out_file.write("Filename\tAverage_pLDDT\n")
        for pdb_file, avg_plddt in results:
            out_file.write(f"{pdb_file}\t{avg_plddt:.2f}\n")

    print(f"Results saved to {output_tsv}")

def main():
    """
    Main function to parse arguments and execute the script.
    """
    parser = argparse.ArgumentParser(description="Calculate average pLDDT scores for PDB files from a single file or a folder.")
    parser.add_argument('--input', '-i', type=str, required=True, help='Path to a PDB file or a folder containing PDB files.')
    parser.add_argument('--output', '-o', type=str, required=True, help='Path to the output TSV file.')
    parser.add_argument('--threads', '-t', type=int, default=multiprocessing.cpu_count(), help='Number of threads to use (default: all available cores).')

    args = parser.parse_args()

    process_pdb_files(args.input, args.output, args.threads)

if __name__ == "__main__":
    main()