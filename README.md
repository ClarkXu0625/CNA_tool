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
    - 
<pre lang="md"> ``` cna_tool/ ├── __init__.py ├── infer.py # Top-level wrapper: infer_cnas_from_scrna ├── cna_inference.py # Core CNAInferer class ├── preprocessing.py # Gene filtering, normalization, GTF loading ├── utils.py # Helper functions: sliding windows, formatting └── data/ └── gencode.v47.annotation.gtf # Optional GTF file tests/ ├── test_inference.py └── example_usage.ipynb README.md setup.py ``` </pre>

# APIs