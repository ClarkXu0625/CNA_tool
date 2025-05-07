# CNA_tool
A tool to infer CNAs from scRNA-seq data.

*please read through the following content*


## Table of Contents

1. [Introduction](#cna_tool)
2. [Setup](#setup)
3. [Minimal Working Example](#setup)
4. [Reproducing Results](#want-to-replicate-our-result)
5. [Code Structure](#code-structure)
6. [Deliverables](#deliverables)
7. [Data Used](#data-used)  
   7.1 [Data Structure](#data-structure)  
   7.2 [Description](#description)
8. [APIs](#apis)  
   8.1 [Function and Class Overview](#apis)


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
    sys.path.append('.src')
    import cna_tool
    from cna_tool import infer_cnas_from_scrna

    adata = sc.read_h5ad("PBMC_simulated_cnas.h5ad")
    adata = infer_cnas_from_scrna(adata)


## Want to replicate our result?

You may start with notebooks/sampleworkflow
prerequisit: Open download num3_adata_filtered.h5ad , num5_adata_filtered.h5ad, num7_adata_filtered.h5ad, and PBMC_simulated_cnas_041025.h5ad from [Google drive](https://drive.google.com/drive/folders/10LGU_CHLHJkABwpyEqT1XuFGGef-xL7m), place it under /data. Then, you are good to go!

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
  - Task2b_task4.ipynb
  - Task2b.ipynb
  - sample_workflow.ipynb
- raw_data_extraction/
  - data3.ipynb: concatenate raw dataset3 to a h5ad file
  - data5.ipynb: concatenate raw dataset5 to a h5ad file
  - data7.ipynb: concatenate raw dataset7 to a h5ad file
  - map_gene_coordinates.ipynb: map chromosome coordinate in raw data

## Deliverables
- Task 1: 
  - src/cna_tool: the package code
  - notbooks/sample_workflow.ipynb: section "Task 1"
- Task 2
  - Task 2a: notebooks/Task2a.ipynb
  - Task 2b: notebooks/Task2b_task4b.ipynb: section "Task2b"
- Task 3
  - notebooks/sample_workflow.ipynb: section "Task 3"
- Task 4
  - notebooks/Task2b_task4b.ipynb: section "Task4"

## Data Used

Dataset being used in this projects are dataset [3](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263152), [5](//www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194214), and [7](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195467
)
All data used for analysis could found in [Google Drive](https://drive.google.com/drive/folders/10LGU_CHLHJkABwpyEqT1XuFGGef-xL7m?usp=sharing)


### Data Structure
- processed_data/
  - labeled_data/
    - num3_adata_filtered.h5ad
    - num5_adata_filtered.h5ad
    - num7_adata_filtered.h5ad
  - cna_profile_results/
    - adata3_results.h5ad
    - adata5_results.h5ad
    - adata7_results.h5ad
  - PBMC_simulated_cnas_041025.h5ad
  - num3_adata.h5ad
  - num5_adata.h5ad
  - num7_adata.h5ad

### Description
- num*_adata.h5ad: Raw data (different formats) being converted to adata structure without any preprocessing. The code to obtain the h5ad raw data could be found in /raw_data_extraction/. These datasets are the raw data that needs downstream process (i.e. chromosome coordinate assignment) and will be used in task 3.
- PBMC_simulated_cnas_041025.h5ad: Raw labeled data that is used for task 1 and 2a. You may also see test_adata.h5ad in some jupyter notebooks. They are the same thing. Just change the name from "test_adata.h5ad" to "PBMC_simulated_cnas_041025.h5ad"
- labaled_data/num*_adata_filtered.h5ad: Raw data being processed to add gene coordinate. These data are used to perform task 3.
- cna_profile_results/adata*_results.h5ad: CNA inferred result is generated and stored in .obs['cna_profile']. These datasets are generated from task 3.



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
