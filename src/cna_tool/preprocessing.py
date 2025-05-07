import scanpy as sc
import numpy as np
import pandas as pd



def annotation_preprocess(adata, n_neighbors=20, n_pcs=10, plot=True):
    """
    Preprocesses an AnnData object by filtering, normalizing, and computing embeddings for clustering and visualization.

    Parameters:
    - adata (AnnData): Input AnnData object (raw counts).
    - n_neighbors (int): Number of neighbors to use when building the neighborhood graph (default: 20).
    - n_pcs (int): Number of principal components to use for PCA and neighbors (default: 10).
    - plot (bool): Whether to plot UMAP with Leiden clusters (default: True).

    Returns:
    - AnnData: A processed AnnData object with filtering, normalization, PCA, neighbors, clustering, and UMAP embeddings.
    """
    data_processed = adata.copy()
    data_processed.var_names_make_unique()
    print("Before filtering: ", data_processed.shape)

    # filter cells based on MT content (cell quality control already performed)
    data_processed.var['mt'] = data_processed.var_names.str.startswith('MT-')
    ribo_prefix = ("RPS","RPL")
    data_processed.var['ribo'] = data_processed.var_names.str.startswith(ribo_prefix)
    sc.pp.calculate_qc_metrics(data_processed, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)

    # Filtering, keep cells with fewer than 20% mitochondrially encoded gene total UMIs
    data_processed = data_processed[data_processed.obs['pct_counts_mt']<20,:].copy()

    # filter based on total number of genes detected (at least 500)
    sc.pp.filter_cells(data_processed, min_genes=500)

    # filter based on total number of counts
    sc.pp.filter_cells(data_processed, max_counts=30000)

    # keep genes that are detected in at least 3 cells
    sc.pp.filter_genes(data_processed, min_cells=3)
    print("After filtering: ", data_processed.shape)

    # normalize
    data_processed.layers['counts'] = data_processed.X.copy()
    sc.pp.normalize_total(data_processed , target_sum=1e4)
    sc.pp.log1p(data_processed )
    sc.pp.highly_variable_genes(data_processed , min_mean=0.0125, max_mean=6, min_disp=0.25)

    # pca
    sc.tl.pca(data_processed)
    sc.pp.neighbors(data_processed, n_neighbors=n_neighbors, n_pcs=n_pcs)

    #sc.pp.neighbors(data_processed, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.leiden(data_processed,.1)
    sc.tl.paga(data_processed)
    sc.pl.paga(data_processed, plot=False)
    sc.tl.umap(data_processed, 0.25, init_pos='paga')
    if plot:
        sc.pl.umap(data_processed,color=['leiden'], alpha=.75, s=15, legend_loc='on data')
    return data_processed

def simple_preprocess(adata):
    """
    Applies basic preprocessing to an AnnData object:
    - Normalizes total counts per cell to a fixed target (10,000),
    - Applies log1p transformation to stabilize variance.

    Parameters:
    - adata: AnnData object containing raw count data.

    Returns:
    - adata: Preprocessed AnnData object with normalized and log-transformed expression values.
    """
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata
