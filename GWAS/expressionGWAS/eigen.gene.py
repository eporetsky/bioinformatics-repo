import pandas as pd
import numpy as np
import glob
import os
from sklearn.decomposition import PCA

CLS_DIR = 'CLS'
DATA_DIR = 'data'
OUT_DIR_PC1 = 'phe'
OUT_DIR_MEAN = 'mean'

os.makedirs(OUT_DIR_PC1, exist_ok=True)
os.makedirs(OUT_DIR_MEAN, exist_ok=True)

exp_files = glob.glob(os.path.join(DATA_DIR, '*.tsv'))
conditions = [os.path.splitext(os.path.basename(f))[0] for f in exp_files]

for cond in conditions:
    print(f'Processing {cond}...')
    cls_path = os.path.join(CLS_DIR, f'{cond}.cls5.tsv')
    data_path = os.path.join(DATA_DIR, f'{cond}.tsv')
    if not os.path.exists(cls_path) or not os.path.exists(data_path):
        print(f'  Skipping (missing file): {cond}')
        continue

    clusters = pd.read_csv(cls_path, sep='\t', usecols=['clusterID','geneID'])
    clust_groups = clusters.groupby('clusterID')['geneID'].apply(list).to_dict()
    expr = pd.read_csv(data_path, sep='\t', index_col=0)

    # Step 1: log2(CPM+1)
    expr_logged = np.log2(expr + 1)
    # Step 2: gene-wise z-score (center and scale)
    expr_z = expr_logged.sub(expr_logged.mean(axis=1), axis=0)
    expr_z = expr_z.div(expr_logged.std(axis=1), axis=0)

    all_pc1 = {}
    all_means = {}
    taxa = expr.columns.tolist()

    for clust, geneids in clust_groups.items():
        sub = expr_z.loc[expr_z.index.intersection(geneids)]
        if sub.shape[0] < 2:
            continue
        vals = sub.T.dropna(axis=0, how='any')
        if vals.shape[0] < 2:
            continue
        # PC1
        pca = PCA(n_components=1)
        pc1 = pca.fit_transform(vals.values)
        pc1 = pd.Series(pc1.flatten(), index=vals.index)
        key = f'{cond}_{clust}'
        all_pc1[key] = pc1

        # Mean cluster expression for each sample/taxa
        all_means[key] = vals.mean(axis=1)

    if not all_pc1:
        print(f'  No clusters with at least two genes and two taxa for {cond}')
        continue

    # Save PC1
    result = pd.DataFrame(all_pc1)
    result.index.name = 'Taxa'
    result = result.reset_index()
    output_pc1 = os.path.join(OUT_DIR_PC1, f'{cond}.tsv')
    result.to_csv(output_pc1, sep='\t', index=False)
    print(f'  Wrote {output_pc1}')

    # Save mean matrix
    mean_result = pd.DataFrame(all_means)
    mean_result.index.name = 'Taxa'
    mean_result = mean_result.reset_index()
    output_mean = os.path.join(OUT_DIR_MEAN, f'{cond}.tsv')
    mean_result.to_csv(output_mean, sep='\t', index=False)
    print(f'  Wrote {output_mean}')
print('DONE.')
