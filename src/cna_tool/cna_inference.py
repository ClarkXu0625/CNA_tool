from .utils import ensure_gene_coords, normalize_expr, sliding_window_segments
import pandas as pd
import gc
import numpy as np

class CNAInferer:
    def __init__(self,
                 adata,
                 control_adata,
                 gtf_df: pd.DataFrame=None,
                 window=50,
                 gain_thr=0.2,
                 loss_thr=-0.2,
                 norm_method='log2_ratio'):
        # 1) ensure coords
        self.adata = ensure_gene_coords(adata.copy(), gtf_df)
        # 2) align control
        self.control = control_adata[:, self.adata.var_names].copy()
        # 3) params
        self.window      = window
        self.gain_thr    = gain_thr
        self.loss_thr    = loss_thr
        self.norm_method = norm_method

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
