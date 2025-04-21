from utils import ensure_gene_coords, normalize_expr, sliding_window_segments
import pandas as pd

class CNAInferer:
    def __init__(self,
                 adata,
                 control_adata,
                 gtf_df: pd.DataFrame=None,
                 window=50,
                 gain_thr=0.2,
                 loss_thr=-0.2,
                 norm_method='zscore'):
        # 1) ensure coords
        self.adata = ensure_gene_coords(adata.copy(), gtf_df)
        # 2) align control
        self.control = control_adata[:, self.adata.var_names].copy()
        # 3) params
        self.window      = window
        self.gain_thr    = gain_thr
        self.loss_thr    = loss_thr
        self.norm_method = norm_method

    def infer(self):
        # gene-specific normalization
        Z = normalize_expr(self.adata, self.control, method=self.norm_method)

        # segmentation call
        segs = sliding_window_segments(
            Z,
            self.adata.var[['chromosome','start','end']],
            window=self.window,
            gain_thr=self.gain_thr,
            loss_thr=self.loss_thr
        )
        self.adata.uns['cna_segments'] = segs

        # build per-cell calls
        profiles = []
        for ci in range(Z.shape[0]):
            calls = []
            for seg in segs:
                mask = (
                    (self.adata.var['chromosome']==seg['chrom']) &
                    (self.adata.var['start']  >= seg['start']) &
                    (self.adata.var['end']    <= seg['end'])
                )
                subz = Z[ci, mask.values]
                if subz.size and ((seg['type']=='gain' and subz.mean()>0) or
                                  (seg['type']=='loss' and subz.mean()<0)):
                    calls.append(f"{seg['chrom']}:{seg['start']}-{seg['end']}({seg['type']})")
            profiles.append(";".join(calls))
        self.adata.obs['cna_profile'] = profiles
        return self.adata


# Replace your old test_pipeline_on_slice with this:

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
    norm_method: str = 'zscore'
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
