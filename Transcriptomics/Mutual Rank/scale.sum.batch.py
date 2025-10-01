import os
import numpy as np
import pandas as pd
from scipy.sparse import load_npz, csr_matrix
import subprocess
import gzip
import itertools
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# conda install bioconda::mcl (mcl-22.282)

def threshold_and_normalize(mat, threshold=5):
    # Remove values below threshold
    mat = mat.tocsr()
    mat.data[mat.data < threshold] = 0
    mat.eliminate_zeros()
    # Row-normalize (each row sums to 1)
    row_sums = np.array(mat.sum(axis=1)).flatten()
    row_indices, col_indices = mat.nonzero()
    data = mat.data / row_sums[row_indices]
    norm_mat = csr_matrix((data, (row_indices, col_indices)), shape=mat.shape)
    return mat

# --- WGCNA-style diagnostics ---
def row_sum_normalize(mat):
    mat = mat.tocsr()
    row_sums = np.array(mat.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1  # avoid division by zero
    row_indices, col_indices = mat.nonzero()
    data = mat.data / row_sums[row_indices]
    return csr_matrix((data, (row_indices, col_indices)), shape=mat.shape)

def soft_threshold(mat, beta):
    mat = mat.tocsr().copy()
    mat.data = np.power(mat.data, beta)
    return mat

def node_connectivity(mat):
    return np.array(mat.sum(axis=1)).flatten()

def scale_free_fit(connectivity):
    connectivity = connectivity[connectivity > 0]
    if len(connectivity) < 10:
        return 0
    hist, bin_edges = np.histogram(connectivity, bins=30)
    x = np.log10(bin_edges[1:][hist > 0]).reshape(-1, 1)
    y = np.log10(hist[hist > 0])
    if len(x) < 2:
        return 0
    model = LinearRegression().fit(x, y)
    r2 = model.score(x, y)
    return r2

def main():
    consistency_dir = 'consistency'
    mcl_dir = 'mcl'
    data_dir = 'data'

    threshold = 10
    
    os.makedirs(mcl_dir, exist_ok=True)
    for fname in os.listdir(consistency_dir):
        if not fname.endswith('.consistency.npz'):
            continue
        species = fname.replace('.consistency.npz', '')
        if species not in ["AtCol-0"]: #,"ZmB73"
             continue
        print(f'Processing {species}')
        npz_path = os.path.join(consistency_dir, fname)
        cmp_dir = os.path.join(data_dir, species, 'cpm')
        mat = load_npz(npz_path)
        # --- WGCNA-style diagnostics ---
        mat_norm = row_sum_normalize(mat)
        betas = [round(x, 2) for x in np.arange(1.0, 2.01, 0.1)]
        diag_png = os.path.join(mcl_dir, f'{species}.wgcna_diag.png')
        r2s = []
        mean_conns = []
        for beta in betas:
            print(f'  Processing beta: {beta}')
            mat_beta = soft_threshold(mat_norm, beta)
            conn = node_connectivity(mat_beta)
            r2 = scale_free_fit(conn)
            r2s.append(r2)
            mean_conns.append(np.mean(conn))
        fig, axs = plt.subplots(1, 2, figsize=(12, 5))
        axs[0].plot(betas, r2s, marker='o')
        axs[0].set_xlabel('Soft-threshold power (β)')
        axs[0].set_ylabel('Scale-free topology fit (R²)')
        axs[0].set_title(f'Scale-free Topology Fit: {species}')
        axs[1].plot(betas, mean_conns, marker='o')
        axs[1].set_xlabel('Soft-threshold power (β)')
        axs[1].set_ylabel('Mean connectivity')
        axs[1].set_title(f'Mean Connectivity: {species}')
        plt.tight_layout()
        plt.savefig(diag_png)
        plt.close()
        print(f'WGCNA diagnostic plot saved: {diag_png}')
        # Print the lowest beta where R2 >= 0.8
        b_candidates = [b for b, r2 in zip(betas, r2s) if r2 >= 0.8]
        if b_candidates:
            print(f'For {species}, minimum beta (B) with R² >= 0.8: {b_candidates[0]}')
        else:
            print(f'For {species}, no beta (B) found with R² >= 0.8')

if __name__ == '__main__':
    main() 