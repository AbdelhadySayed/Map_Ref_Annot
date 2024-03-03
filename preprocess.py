
import numpy as np
import pandas as pd
import scanpy as sc
import anndata


def preprocess(adata_com, mt=20, ribo=62, min_genes=200, min_cells=3, scale=True):
    # mitochondrial genes
    adata_com.var["mt"] = adata_com.var_names.str.startswith("MT-")
    # ribosomal genes
    adata_com.var["ribo"] = adata_com.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
    adata_com, qc_vars=["mt", "ribo"], inplace=True, percent_top=[20], log1p=False)

    adata2 = adata_com[(adata_com.obs.pct_counts_mt<mt) & (adata_com.obs.pct_counts_ribo<ribo),:]

    sc.pp.filter_cells(adata2, min_genes=min_genes) #get rid of cells with fewer than 200 genes
    sc.pp.filter_genes(adata2,min_cells=min_cells) #get rid of genes that are found in fewer than 3 cells
   
    adata2 = adata2[:, (adata2.var["mt"]==False) & (adata2.var["ribo"]==False)]
    adata2.layers['counts'] = adata2.X.copy()
    sc.pp.normalize_total(adata2, target_sum=1e4)
    sc.pp.log1p(adata2) #change to log counts
    sc.pp.highly_variable_genes(adata2) #these are default values

    adata_processed = adata2[:, adata2.var.highly_variable] #filter highly variable
    if scale:
        sc.pp.scale(adata_processed, max_value=10) #scale each gene to unit variance
        
    return adata_processed