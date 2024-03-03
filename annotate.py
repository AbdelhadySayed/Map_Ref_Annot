import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from sklearn.model_selection import train_test_split
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.neighbors import KNeighborsClassifier
import torch
import scvi
import matplotlib.pyplot as plt


def annotate(adata_scvi, method="cosine", n_neighbors=30, cell_type_col="cell_type_ref"):
    lab_data = adata_scvi[adata_scvi.obs[cell_type_col] != "not defined"]
    unlab_data = adata_scvi[adata_scvi.obs[cell_type_col] == "not defined"]
    # if method =="cosine":
    cos = cosine_similarity(unlab_data.X, lab_data.X, dense_output=True)
    cos_df = pd.DataFrame(cos, index=unlab_data.obs.index, columns=lab_data.obs[cell_type_col])
    #fine cell types 
    cell_types = []
    for cell in cos_df.index:
        i = np.argmax(cos_df.loc[cell])
        cell_type = cos_df.columns[i]
        cell_types.append([cell, cell_type])
    cell_type_df = pd.DataFrame(cell_types, columns=["id", "type"])
    cell_type_df.set_index("id", inplace=True)
    unlab_data.obs["fine_cell_types"] = cell_type_df["type"]
    
    # major cell types
    neigh = KNeighborsClassifier(n_neighbors=n_neighbors, metric=method)
    neigh.fit(lab_data.X.toarray(), lab_data.obs[cell_type_col])
    # transfer cell types to unlab data 
    unlab_data.obs['major_cell_types'] = neigh.predict(unlab_data.X)
    # visualize the unlab data after asigning cell types major and fine
    sc.pl.umap(unlab_data, color = ['major_cell_types'], frameon = False)
    
    sc.pl.umap(unlab_data, color = ['fine_cell_types'], frameon = False)
    return unlab_data, unlab_data.obs.major_cell_types, unlab_data.obs.fine_cell_types