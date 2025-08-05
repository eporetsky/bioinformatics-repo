import time
import psutil
import gc
import pigz
import pandas as pd
import numpy as np
from pynetcor.cor import corrcoef
import os

if __name__ == "__main__":
    script_start = time.time()
    base_data_dir = "data"
    mr_threshold = 50
    threads = 72  # Adjust as needed
    output_base = "mr"
    os.makedirs(output_base, exist_ok=True)

    def print_mem(msg):
        process = psutil.Process()
        mem = process.memory_info().rss / 1024 ** 2
        print(f"{msg} | Memory usage: {mem:.2f} MB")

    for species in os.listdir(base_data_dir):
        species_dir = os.path.join(base_data_dir, species)
        cpm_dir = os.path.join(species_dir, "cpm")
        if not os.path.isdir(cpm_dir):
            continue
        print(f"Processing species: {species}")
        output_species_dir = os.path.join(output_base, species)
        os.makedirs(output_species_dir, exist_ok=True)
        for fname in os.listdir(cpm_dir):
            if not fname.endswith(".cpm.tsv.gz"):
                continue
            experiment = fname.replace(".cpm.tsv.gz", "")
            file_path = os.path.join(cpm_dir, fname)
            output_path = os.path.join(output_species_dir, f"{experiment}.mr.tsv.gz")
            if os.path.exists(output_path):
                print(f"  Skipping {species} {experiment}: output already exists.")
                continue
            print(f"  Processing experiment: {experiment}")
            t0 = time.time()
            print_mem("  Reading CPM file")
            df = pd.read_csv(file_path, sep="\t", index_col="geneID")
            # Skip if not enough columns
            if df.shape[1] < 12:
                print(f"  Skipping {species} {experiment}: less than 12 columns.")
                continue
            # Zero variance row filter
            df = df.loc[df.var(axis=1) > 0]
            # log2(cpm+1) transformation
            df = np.log2(df + 1)
            print_mem("  Calculating correlation matrix")
            corr_np = corrcoef(df, threads=threads)
            corr_np = corr_np.astype(np.float32)
            print_mem("  Calculating row/col ranks")
            import concurrent.futures
            def rank_row(row):
                return np.argsort(-row).argsort() + 1
            with concurrent.futures.ThreadPoolExecutor() as executor:
                row_ranks = np.array(list(executor.map(rank_row, corr_np)), dtype=np.uint16)
                col_ranks = np.array(list(executor.map(rank_row, corr_np.T)), dtype=np.uint16).T
            del corr_np
            gc.collect()
            print_mem("  Calculating mutual rank matrix")
            mr_np = np.sqrt(row_ranks * col_ranks).astype(np.float32)
            del row_ranks, col_ranks
            gc.collect()
            print_mem("  Filtering and writing output")
            gene_ids = df.index.to_numpy()
            i, j = np.triu_indices_from(mr_np, k=1)
            mask = mr_np[i, j] <= mr_threshold
            i_sel = i[mask]
            j_sel = j[mask]
            mr_sel = mr_np[i_sel, j_sel]
            with pigz.open(output_path, "wt") as f:
                f.write("Gene1\tGene2\tMR\n")
                for a, b, v in zip(gene_ids[i_sel], gene_ids[j_sel], mr_sel):
                    f.write(f"{a}\t{b}\t{v:.6g}\n")
            del mr_np, i, j, mask, i_sel, j_sel, mr_sel
            gc.collect()
            print_mem(f"  Output written for {experiment}")
            print(f"  Done {species} {experiment} in {time.time() - t0:.1f}s")
    print(f"Total runtime: {time.time() - script_start:.1f}s")