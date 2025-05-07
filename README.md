# CNA_tool
A tool to infer CNAs from scRNA-seq data. 

*please read through the following content*


## Setup:

    conda create -n cna_tool python=3.11
    conda activate cna_tool
    cd CNA_tool
    pip install scanpy python-igraph leidenalg scipy numpy umap-learn leidenalg matplotlib scikit-learn mygene pandas
    git clone https://github.com/ClarkXu0625/CNA_tool.git
    cd CNA_tool


A minimal working example of using the tool:

    import scanpy as sc
    import sys
    sys.path.append('../src')
    import cna_tool
    from cna_tool import infer_cnas_from_scrna

    adata = sc.read_h5ad("PBMC_simulated_cnas.h5ad")
    adata = infer_cnas_from_scrna(adata)


## Code Structure
CNA_Tool
- src/cna_tool/
    - \_\_init__.py
    - infer.py
    - preprocessing.py
    - cna_inference.py
    - utils.py
    - tl/
      - \_\_init__.py
- notebooks/
  - Task2b.ipynb
  - 
- raw_data_extraction/
  - dataset5.ipynb

## Data

Dataset being used in this projects are dataset [3](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263152), [4](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277604), [5](//www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194214), [6](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131736), and [7](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195467
)
The data could found in [Google Drive](https://drive.google.com/drive/folders/10LGU_CHLHJkABwpyEqT1XuFGGef-xL7m?usp=sharing)

### Data Structure
- processed_data/
  - labeled_data/
    - num3_adata_filtered.h5ad
    - num4_human_adata_filtered.h5ad
    - num6_adata_filtered.h5ad
    - num7_adata_filtered.h5ad
  - PBMC_simulated_cnas_041025.h5ad
  - num3_adata.h5ad
  - num4_adata.h5ad
  - num6_adata.h5ad
  - num7_adata.h5ad

### Description
- num*_adata.h5ad: Raw data (different formats) being converted to adata structure without any preprocessing. The code to obtain the h5ad raw data could be found in /raw_data_extraction/
- PBMC_simulated_cnas_041025.h5ad: Raw labeled data get from the website


## APIs
Import the functions in the following format
    import sys
    sys.path.append('src')
    import cna_tool
    from cna_tool import infer_cnas_from_scrna
    
For the following function/classes, please type help(function) to get detailed description. Sample usage:

    help(infer_cnas_from_scrna)

You may also view [Tutorial](Tutorial.ipynb) for api call

| Function/Classes | Parameters | Description | 
| ---------- | ---- | -------------------------- |
| infer.infer_cnas_from_scrna | adata, cna_label_col=None, diploid_labels=[''], control_mask=None, gtf_df=None, window=50, gain_thr=0.15, loss_thr=-0.15, norm_method='zscore', require_gene_coords=True| Infers DNA copy number alterations (CNAs) from an AnnData object containing scRNA-seq data. |
| CNAInferer | adata, control_adata |main object |
| CNAInferer.infer | self |Infer CNA from given adata |
| prepare_control | adata, time_key, ctrl_label, min_genes, n_subsample | Extracts and preprocesses control cells from an AnnData object |
| infer_cna_by_timepoint | adata, control, gtf_df=None, time_key, exclude, gain_thr, loss_thr, window | Infers Copy Number Alterations (CNAs) for each timepoint relative to a control dataset |
| utils.select_control_mask| adata, obs_key: str, control_values | Build a boolean mask selecting control (diploid) cells| 
| utils.ensure_gene_coords | adata, gtf_df: pd.DataFrame=None | If 'chromosome','start','end' are present in adata.var, do nothing. Else merge with user‑supplied gtf_df. Else fetch missing via pybiomart (if installed). Drops any genes still lacking coords.|
| utils.normalize_expr | adata, control_adata, method='log2_ratio' | Gene‐specific normalization against diploid control, method could be either zscore or log2_ratio |
| utils.sliding_window_segments | Z, var_df, window=50, gain_thr=0.2, loss_thr=-0.2 | Identifies copy number alteration (CNA) segments from smoothed, average z-scored expression across genes using a sliding window approach. |
| utils.fetch_gene_coordinates_if_missing | adata, gtf_df | Map the gene coordinate from provided df generated from gtf file, creating .var columns "chromosome", "start", and "end" |
| utils.map_gene_coordinates | adata, adata | Map the gene coordinate from one adata (labeled) to another, creating .var columns "chromosome", "start", and "end" |
| preprocessing.annotation_preprocess | adata | full preprocess pipeline: filter by counts and mt content; normalize; perform clustering, and umap |
| preprocessing.simple_preprocess | adata | normalize and log transfer the adata |
| tl.evaluate_cna_call |adata, truth_col, pred_col |  Evaluates the performance of CNA prediction by comparing predicted CNA profiles with simulated ground truth labels |
| tl.run_cna_evaluation | adata, params | Runs the full CNA inference and evaluation pipeline using the provided parameters |
