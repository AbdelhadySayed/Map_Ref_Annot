import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import torch
import scvi
import matplotlib.pyplot as plt

def integrate(adata_processed, batch_col="study_id", n_latent=20, max_epochs=400, vis_latent=True):
    adata_scvi = adata_processed.copy()
    scvi.model.SCVI.setup_anndata(adata_scvi, layer="counts", batch_key=batch_col)
    max_epochs_scvi = max_epochs 
    #if max_epochs:
    #    max_epochs_scvi = max_epochs
    #else:
    #    max_epochs_scvi = np.min([round((20000 / data.n_obs) * 400), 400])
    model_scvi = scvi.model.SCVI(adata_scvi, n_latent=n_latent)
    model_scvi.train(max_epochs=max_epochs)
    adata_scvi.obsm["X_scVI"] = model_scvi.get_latent_representation()

    adata_scvi.layers['scvi_normalized'] = model_scvi.get_normalized_expression(library_size = 1e4)
    sc.pp.neighbors(adata_scvi, use_rep="X_scVI")
    sc.tl.umap(adata_scvi)
    sc.tl.leiden(adata_scvi)
    if vis_latent:
        sc.pl.umap(adata_scvi, color = ['leiden', batch_col], frameon = False)
    return adata_scvi