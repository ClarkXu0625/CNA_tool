from .utils import ensure_gene_coords, normalize_expr, sliding_window_segments
import pandas as pd
import gc
import numpy as np
import scanpy as sc
import anndata as ad

class CNAInferer:
    def __init__(self,
                 adata,
                 control_adata,
                 gtf_df: pd.DataFrame = None,
                 window=50,
                 gain_thr=0.2,
                 loss_thr=-0.2,
                 norm_method='log2_ratio'):
        """
        Initializes the CNAInferer with test and control datasets and preprocessing parameters.

        Parameters:
        - adata: AnnData object representing the test cells to be analyzed.
        - control_adata: AnnData object with control/reference cells.
        - gtf_df: Optional pandas DataFrame containing gene coordinate annotations
                  with columns like 'gene_name', 'chromosome', 'start', 'end'.
        - window: Integer, size of the smoothing window for signal processing.
        - gain_thr: Float, threshold above which a region is considered a copy number gain.
        - loss_thr: Float, threshold below which a region is considered a copy number loss.
        - norm_method: Method used to normalize expression between test and control
                       (e.g., 'log2_ratio').
        """

        # Step 1: Ensure gene coordinate information is present in adata.var
        self.adata = ensure_gene_coords(adata.copy(), gtf_df)

        # Step 2: Align the control to have the same genes and order as the test set
        self.control = control_adata[:, self.adata.var_names].copy()

        # Step 3: Store inference parameters
        self.window      = window       # Smoothing window size
        self.gain_thr    = gain_thr     # Gain threshold for CNA calling
        self.loss_thr    = loss_thr     # Loss threshold for CNA calling
        self.norm_method = norm_method  # Normalization method to apply (e.g., log2 ratio)

    # def infer(self):
    #     # gene-specific normalization
    #     Z = normalize_expr(self.adata, self.control, method=self.norm_method)

    #     # segmentation call
    #     segs = sliding_window_segments(
    #         Z,
    #         self.adata.var[['chromosome','start','end']],
    #         window=self.window,
    #         gain_thr=self.gain_thr,
    #         loss_thr=self.loss_thr
    #     )
    #     self.adata.uns['cna_segments'] = segs

    #     # build per-cell calls
    #     profiles = []
    #     for ci in range(Z.shape[0]):
    #         calls = []
    #         for seg in segs:
    #             mask = (
    #                 (self.adata.var['chromosome']==seg['chrom']) &
    #                 (self.adata.var['start']  >= seg['start']) &
    #                 (self.adata.var['end']    <= seg['end'])
    #             )
    #             subz = Z[ci, mask.values]
    #             if subz.size and ((seg['type']=='gain' and subz.mean()>0) or
    #                               (seg['type']=='loss' and subz.mean()<0)):
    #                 calls.append(f"{seg['chrom']}:{seg['start']}-{seg['end']}({seg['type']})")
    #         profiles.append(";".join(calls))
    #     self.adata.obs['cna_profile'] = profiles
    #     return self.adata
    
    def infer(self):
        """
        Perform copy number alteration (CNA) inference on the input AnnData object.

        Steps:
        - Normalize gene expression between test and control cells.
        - Identify genome segments with potential CNAs using a sliding window.
        - For each cell, determine which CNA segments are active based on mean expression.
        - Store CNA profile strings in `adata.obs['cna_profile']`.
        
        Returns:
        - Updated AnnData object with CNA calls stored in `.obs['cna_profile']`
        and CNA segment definitions in `.uns['cna_segments']`.
        """
        # Normalize expression (row-wise normalization preferred)
        Z = normalize_expr(self.adata, self.control, method=self.norm_method)

        # Convert to float32 for memory savings
        if Z.dtype != np.float32:
            Z = Z.astype(np.float32, copy=False)

        # Segment genome (runs only once)
        segs = sliding_window_segments(
            Z,
            self.adata.var[['chromosome', 'start', 'end']],
            window=self.window,
            gain_thr=self.gain_thr,
            loss_thr=self.loss_thr
        )
        self.adata.uns['cna_segments'] = segs

        # Precompute segment-gene indices only once
        chrom = self.adata.var['chromosome'].to_numpy()
        start = self.adata.var['start'].to_numpy()
        end = self.adata.var['end'].to_numpy()

        segment_gene_indices = []
        for seg in segs:
            mask = (chrom == seg['chrom']) & (start >= seg['start']) & (end <= seg['end'])
            idxs = np.flatnonzero(mask)
            segment_gene_indices.append(idxs)

        # Free variables not needed anymore
        del chrom, start, end
        gc.collect()

        # Write CNA profiles in chunks to avoid high memory usage
        n_cells = Z.shape[0]
        batch_size = 500  # adjust to fit memory constraints
        profile_chunks = []

        for i in range(0, n_cells, batch_size):
            batch_profiles = []
            Z_batch = Z[i:i + batch_size]  # get only current batch

            for ci in range(Z_batch.shape[0]):
                calls = []
                for seg, idxs in zip(segs, segment_gene_indices):
                    subz = Z_batch[ci, idxs]
                    if subz.size and (
                        (seg['type'] == 'gain' and subz.mean() > 0) or
                        (seg['type'] == 'loss' and subz.mean() < 0)
                    ):
                        calls.append(f"{seg['chrom']}:{seg['start']}-{seg['end']}({seg['type']})")
                batch_profiles.append(";".join(calls))

            profile_chunks.extend(batch_profiles)

            # Clean up
            del Z_batch, batch_profiles
            gc.collect()

        # Store the final profile list
        self.adata.obs['cna_profile'] = profile_chunks

        return self.adata



def prepare_control(adata, time_key='timepoint', ctrl_label='D0', min_genes=200, n_subsample=500):
    """
    Extracts and preprocesses control cells from an AnnData object.

    Parameters:
    - adata: AnnData object containing the full dataset.
    - time_key: The column in adata.obs that indicates the timepoint or condition.
    - ctrl_label: The value in time_key that corresponds to the control group (e.g., 'D0').
    - min_genes: Minimum number of genes a cell must express to be kept.
    - n_subsample: Number of control cells to randomly sample, if more than this number exist.

    Returns:
    - control: Preprocessed AnnData object containing a filtered and possibly subsampled control group.
    """
    ctrl_mask = adata.obs[time_key] == ctrl_label
    control = adata[ctrl_mask].copy()
    
    # Optional: filter low-quality cells
    sc.pp.filter_cells(control, min_genes=min_genes)
    
    # Subsample to reduce size
    if control.n_obs > n_subsample:
        idx = np.random.choice(control.n_obs, n_subsample, replace=False)
        control = control[idx]
    
    return control

def infer_cna_by_timepoint(adata, control, gtf_df=None, time_key='timepoint', exclude=['D0'], gain_thr=0.1, loss_thr=-0.5, window=30):
    """
    Infers Copy Number Alterations (CNAs) for each timepoint relative to a control dataset.

    Parameters:
    - adata: AnnData object with test cells.
    - control: AnnData object with control cells (from prepare_control).
    - gtf_df: Optional DataFrame with gene annotation (chromosome, start, end).
    - time_key: The column in adata.obs indicating the timepoint of each cell.
    - exclude: List of timepoints to skip (e.g., control timepoint like 'D0').
    - gain_thr: Threshold above which a region is considered a CNA gain.
    - loss_thr: Threshold below which a region is considered a CNA loss.
    - window: Window size for smoothing CNA signal across genomic regions.

    Returns:
    - combined: AnnData object combining CNA-inferred results across all timepoints.
    """
    results = []
    timepoints = sorted(set(adata.obs[time_key]) - set(exclude))

    for tp in timepoints:
        print(f"[INFO] Inferring CNAs for {tp}...")
        tp_mask = adata.obs[time_key] == tp
        test = adata[tp_mask].copy()
        
        inferer = CNAInferer(test, control, gtf_df=gtf_df, gain_thr=gain_thr, loss_thr=loss_thr, window=window, norm_method='log2_ratio')
        inferred = inferer.infer()
        results.append(inferred)

    # Combine all inferred AnnData
    combined = ad.concat(results, join='outer', index_unique=None)
    return combined


# A sample pipline workflow
def test_pipeline_on_slice(
    adata,
    obs_key_control: str,
    control_values,
    obs_key_test: str,
    test_values,
    chromosome: str,
    gtf_df: pd.DataFrame = None,
    window: int = 20,
    gain_thr: float = 0.2,
    loss_thr: float = -0.2,
    norm_method: str = 'log2_ratio'
    ):
    """
    1) Pick control cells via adata.obs[obs_key_control] ∈ control_values
    2) Pick test  cells via adata.obs[obs_key_test]    ∈  test_values
    3) Subset genes on `chromosome`
    4) Run CNAInferer on that slice
    """
    # build masks
    ctrl_mask = adata.obs[obs_key_control].isin(control_values)
    tst_mask  = adata.obs[obs_key_test]   .isin(test_values)

    # subset
    ctrl = adata[ctrl_mask].copy()
    sub  = adata[tst_mask].copy()
    # keep only genes on the desired chromosome
    sub  = sub[:, sub.var['chromosome'] == chromosome].copy()
    ctrl = ctrl[:, sub.var_names].copy()

    # infer
    inferer = CNAInferer(
        adata        = sub,
        control_adata= ctrl,
        gtf_df       = gtf_df,
        window       = window,
        gain_thr     = gain_thr,
        loss_thr     = loss_thr,
        norm_method  = norm_method
    )
    return inferer.infer()
