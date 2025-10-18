import pandas as pd
import numpy as np
import gzip
import os

# File paths
GENO_CONV_PATH = 'name_conversion.tsv'
SAMPLE_ANN_PATH = 'PRJNA383416.samples.tsv'
CPM_PATH = 'PRJNA383416.cpm.tsv.gz'
OUTPUT_DIR = 'data/'  # Saves {condition}.tsv here

def main():
    """Process CPM expression matrix per condition keeping only allowed genotypes, renaming, averaging and log2-transforming."""
    # 1. Load genotype conversion
    geno_conv = pd.read_csv(GENO_CONV_PATH, sep='\t')
    geno_dict = dict(zip(geno_conv['original'], geno_conv['convert']))

    # 2. Load sample/condition annotation
    sample_ann = pd.read_csv(SAMPLE_ANN_PATH, sep='\t')
    
    # 3. Parse annotation to build {sample: (converted_geno, condition)}
    sample_info = {}
    for _, row in sample_ann.iterrows():
        sample = row['sample']
        name = row['name']
        name_parts = name.split(' ')
        orig_geno = name_parts[0]
        condition = '_'.join([p.replace('(','').replace(')','') for p in name_parts[1:]])
        if orig_geno not in geno_dict:
            continue  # skip, no SNP data
        conv_geno = geno_dict[orig_geno]
        sample_info[sample] = {'genotype': conv_geno, 'condition': condition}

    # Find which columns to keep and their mappings
    used_samples = set(sample_info.keys())
    if not used_samples:
        raise ValueError('No samples left after genotype filtering.')

    # 4. Read CPM matrix (samples are columns)
    # Sniff column names:
    with gzip.open(CPM_PATH, 'rt') as f:
        header = f.readline().strip().split('\t')
    gene_col = header[0]
    exp_samples = header[1:]
    kept = [s for s in exp_samples if s in used_samples]
    
    if not kept:
        raise ValueError('No CPM columns match filtered samples!')

    # 5. Read only needed columns + gene
    usecols = [gene_col] + kept
    expr = pd.read_csv(CPM_PATH, sep='\t', compression='gzip', usecols=usecols)
    expr = expr.set_index(gene_col)

    # 6. For each condition, make new file
    # Map: sample -> (geno, cond)
    sample2geno = {s: sample_info[s]['genotype'] for s in kept}
    sample2cond = {s: sample_info[s]['condition'] for s in kept}

    conditions = sorted(set(sample2cond.values()))
    print(f'Writing {len(conditions)} condition matrix files...')
    for cond in conditions:
        # Samples for this condition
        cond_samples = [s for s in kept if sample2cond[s]==cond]
        if not cond_samples:
            continue
        # Data: subset, rename columns to genotype
        sub = expr[cond_samples].copy()
        # Rename columns: sampleID -> genotype
        sub.columns = [sample2geno[s] for s in cond_samples]
        # Log2(CPM+1)
        #sub = np.log2(sub+1)
        # Group by genotype, then mean (for biological replicates)
        group = sub.T.groupby(by=sub.columns).mean().T
        # Save
        outf = os.path.join(OUTPUT_DIR, f'{cond}.tsv')
        group.to_csv(outf, sep='\t')
        print(f'  - {outf} ({group.shape[1]} genotypes)')

if __name__ == "__main__":
    main()
