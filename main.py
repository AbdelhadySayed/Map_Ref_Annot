import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.neighbors import KNeighborsClassifier
import torch
import scvi

from concat import con_cat
from preprocess import preprocess
from integrate import integrate
from annotate import annotate
import argparse
import matplotlib.pyplot as plt



def main(args):
    
    # query and ref data expected to be in h5ad format till now
    query_data = sc.read(args.query_data_dir)
    ref_data   = sc.read(args.ref_data_dir)

    adata_com = con_cat(query_data, ref_data, args.cell_type_col, args.batch_col, ref_to_target =args.ref_to_query_percent)
    adata_processed = preprocess(adata_com, args.mt, args.ribo, args.min_genes, args.min_cell, args.scale)
    adata_scvi = integrate(adata_processed, args.batch_col, args.n_latent, args.max_epochs, args.vis_latent)
    adata_labelled, major_cells, fine_cells = annotate(adata_scvi, args.method, args.n_neighbors, args.cell_type_col)
    # transfer cell types to final adata
    final_data = sc.read(args.final_data)
    if adata_labelled.obs.index == final_data.obs.index:
        final_data.obs.major_cell_types = major_cells
        final_data.obs.fine_cell_types = fine_cells
    else:
        "Precessed data has rownames different than final data, assure that you make the same processing steps and parameters for final as used in this method"
    return final_data





if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument(
        "--query-data-dir", 
        type=str,
        default="./Myloma.h5ad",
        metavar="DD",
        help="query data directory need to be annotated",
    )

    parser.add_argument(
        "--ref-data-dir", 
        type=str,
        default="./Bone_marrow.h5ad",
        metavar="DD",
        help="reference data directory need to map to",
    )

    parser.add_argument(
        "--final-data-dir", 
        type=str,
        default="./Myloma.h5ad",
        metavar="DD",
        help="query data directory need to be annotated to be returned annotated finally without modification",
    )
#   
    parser.add_argument(
        "--cell-type-col", 
        type=str,
        default="cell_type",
        metavar="DD",
        help="cell type column name in reference data",
    )

    parser.add_argument(
        "--batch-col", 
        type=str,
        default="study",
        metavar="DD",
        help="batch column name in query data",
    )

    parser.add_argument(
        "--ref-to-query-percent",
        type=float,
        default=0.333,
        metavar="N",
        help="Percent of reference data size to query data (default: 1/3)",
    )


    parser.add_argument(
        "--mt",
        type=int,
        default=20,
        metavar="N",
        help="Threshold of mitochonrial genes percent to remove dead cells based on having more than it",
    )

    parser.add_argument(
        "--ribo",
        type=int,
        default=62,
        metavar="N",
        help="Threshold of ribosomal genes percent to remove cells based on having more than it",
    )

    parser.add_argument(
        "--min-genes",
        type=int,
        default=200,
        metavar="N",
        help="Threshold of number genes should be expressed to filter cells based on it",
    )
    
    parser.add_argument(
        "--min-cells",
        type=int,
        default=3,
        metavar="N",
        help="Threshold of number cells should the gene be expressed in to filter genes based on it",
    )

    parser.add_argument(
        "--scale",
        type=bool,
        default=False,
        metavar="N",
        help="Whether to scale data in preprocessing step or not",
    )

    parser.add_argument(
        "--vis-latent",
        type=bool,
        default=False,
        metavar="N",
        help="Whether to view latent spaces after integration or not",
    )

    parser.add_argument(
        "--n-latent",
        type=int,
        default=30,
        metavar="N",
        help="latent space dimension in integration",
    )

    parser.add_argument(
        "--max-epochs",
        type=int,
        default=400,
        metavar="N",
        help="max epochs used in integration with ScVI",
    )

    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=30,
        metavar="N",
        help="no of neighbors should be used in KNN classifier to make label transfer",
    )

    parser.add_argument(
        "--method",
        type=str,
        default="cosine",
        metavar="N",
        help="Similarity method used in KNN classifier",
    )

    args = parser.parse_args()
    main(args)