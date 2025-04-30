import scanpy as sc
import numpy as np
import pandas as pd
from .cna_inference import CNAInferer
from .utils import select_control_mask, fetch_gene_coordinates_if_missing

def infer_cnas_from_scrna(
    adata,
    cna_label_col=None,
    diploid_labels=[''],
    control_mask=None,
    gtf_df=None,
    window=50,
    gain_thr=0.15,
    loss_thr=-0.15,
    norm_method='zscore',
    require_gene_coords=True
):
    """
    Infers DNA copy number alterations (CNAs) from an AnnData object containing scRNA-seq data.

    Parameters
    ----------
    adata : AnnData
        AnnData object with raw gene expression. Should contain gene names in `adata.var['gene_name']`.

    cna_label_col : str or None
        Column in `adata.obs` used to define diploid control cells via `diploid_labels`.
        If None, and `control_mask` is also None, diploid control is chosen heuristically.

    diploid_labels : list of str, default=['']
        Values in `adata.obs[cna_label_col]` that identify diploid cells.

    control_mask : np.ndarray of bool or None
        Optional boolean mask to directly specify diploid control cells.

    gtf_df : pd.DataFrame or None
        Optional dataframe with gene annotations: must include 'gene_name', 'chromosome', 'start', 'end'.

    window : int
        Window size for smoothing gene expression across the genome.

    gain_thr : float
        Threshold above which smoothed z-score is called a 'gain'.

    loss_thr : float
        Threshold below which smoothed z-score is called a 'loss'.

    norm_method : str
        Normalization method: 'zscore' or 'control_median'.

    require_gene_coords : bool
        If True, raises an error if gene coordinates are not available. If False, will attempt to infer or fetch them.

    Returns
    -------
    adata : AnnData
        The modified AnnData object with:
            - adata.uns['cna_segments']: list of detected CNA segments
            - adata.obs['cna_profile']: per-cell CNA genotype string
    """
    # Ensure gene coordinates are available
    if require_gene_coords and not all(col in adata.var.columns for col in ['chromosome', 'start', 'end']):
        if gtf_df is not None:
            print("Mapping gene coordinates from provided GTF...")
            adata = fetch_gene_coordinates_if_missing(adata, gtf_df)
        else:
            raise ValueError("Gene coordinates missing in adata.var and no GTF provided.")

    # Select control vs. test cells
    if control_mask is not None:
        control = adata[control_mask].copy()
        test = adata[~control_mask].copy()
    elif cna_label_col is not None:
        ctrl_mask = select_control_mask(adata, cna_label_col, diploid_labels)
        control = adata[ctrl_mask].copy()
        test = adata[~ctrl_mask].copy()
    else:
        # Fallback heuristic: assume largest cluster is diploid
        largest_group = adata.obs.groupby('cell_type').size().idxmax()
        ctrl_mask = adata.obs['cell_type'] == largest_group
        control = adata[ctrl_mask].copy()
        test = adata[~ctrl_mask].copy()
        print(f"[info] No control mask provided. Using largest cluster ('{largest_group}') as control.")

    # Run inference
    inferer = CNAInferer(
        adata=test,
        control_adata=control,
        gtf_df=gtf_df,
        window=window,
        gain_thr=gain_thr,
        loss_thr=loss_thr,
        norm_method=norm_method
    )
    updated_test = inferer.infer()

    # Attach output to original AnnData
    adata.obs['cna_profile'] = ''
    adata.obs.loc[updated_test.obs_names, 'cna_profile'] = updated_test.obs['cna_profile']
    adata.uns['cna_segments'] = updated_test.uns.get('cna_segments', [])

    return adata
