# CNA_tool
A tool to infer CNAs from scRNA-seq data


## Setup:

    conda create -n cna_tool python=3.11
    conda activate cna_tool
    cd CNA_tool
    pip install scanpy python-igraph leidenalg scipy numpy umap-learn leidenalg matplotlib scikit-learn mygene pandas


A minimal working example of using the tool:


    import scanpy as sc
    from cna_tool import infer_cnas_from_scrna

    adata = sc.read_h5ad("PBMC_simulated_cnas.h5ad")
    adata = infer_cnas_from_scrna(adata)


## code Structure

- cna_tool/
    - \_\_init__.py
    - infer.py
    - preprocessing.py
    - cna_inference.py
    - utils.py
    - tl/


## APIs

| Function/Classes | Parameters | Description | 
| ---------- | ---- | -------------------------- |
| CNAInferer | adata, control_adata |main object |
| CNAInferer.infer | self |Infer CNA from given adata |
| utils.select_control_mask| adata, obs_key: str, control_values | Build a boolean mask selecting control (diploid) cells| 
| utils.ensure_gene_coords | adata, gtf_df: pd.DataFrame=None | If 'chromosome','start','end' are present in adata.var, do nothing. Else merge with user‑supplied gtf_df. Else fetch missing via pybiomart (if installed). Drops any genes still lacking coords.|
| utils.normalize_expr, 
| utils.sliding_window_segments, 
| utils.fetch_gene_coordinates_if_missing, 
| utils.map_gene_coordinates