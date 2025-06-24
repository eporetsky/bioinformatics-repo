import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import os


def pdb_to_sdf(pdb_file, sdf_file):
    """
    Converts a PDB file to an SDF file using RDKit. Overwrites sdf_file if it exists.
    """
    with open(pdb_file, 'r') as f:
        pdb_block = f.read()
    mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not parse PDB file: {pdb_file}")
    # Generate coordinates if needed
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol)
    writer = Chem.SDWriter(sdf_file)
    writer.write(mol)
    writer.close()
    return sdf_file


def merge_sdf_files(sdf_files, merged_sdf):
    with open(merged_sdf, 'w') as outfile:
        for sdf_file in sdf_files:
            with open(sdf_file, 'r') as infile:
                outfile.write(infile.read())


def main():
    parser = argparse.ArgumentParser(description="Merge multiple ligand SDFs and/or metal PDBs into a single SDF file.")
    parser.add_argument('--input_sdf', nargs='*', default=[], help='Input ligand SDF file(s)')
    parser.add_argument('--input_pdb', nargs='*', default=[], help='Input metal PDB file(s)')
    parser.add_argument('--output_sdf', required=True, help='Output merged SDF file')
    args = parser.parse_args()

    sdf_files = list(args.input_sdf)

    # Convert all PDBs to SDFs and add to list
    for pdb_path in args.input_pdb:
        sdf_path = os.path.splitext(pdb_path)[0] + '.sdf'
        pdb_to_sdf(pdb_path, sdf_path)
        sdf_files.append(sdf_path)

    if len(sdf_files) < 2:
        raise ValueError("At least two molecules (from SDF and/or PDB) are required to merge.")

    # Merge SDFs by concatenation
    merge_sdf_files(sdf_files, args.output_sdf)
    print(f"Merged SDF written to {args.output_sdf}")


if __name__ == "__main__":
    main()

