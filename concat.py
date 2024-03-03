
import numpy as np
import pandas as pd
import scanpy as sc
import anndata


def con_cat(target_data, ref_data, cell_type_col, batch_col, ref_to_target =1/3):
    if ref_data.shape[0]>= ref_to_target * target_data.shape[0]:
        print("Warning: ref_data shouldnt be larger than target data, we will splice it")
        percent = target_data.shape[0]/ref_data.shape[0]*ref_to_target
        _,X,_2,y = train_test_split(ref_data.X.toarray(), ref_data.obs[cell_type_col], test_size=percent, random_state=42)
        ref_data = anndata.AnnData(X, dict(obs_names= y.index, cell_type_ref=y), dict(var_names= adata_ref.var.index))
    else:
        ref_data.obs.rename(columns={cell_type_col: "cell_type_ref"}, inplace=True)
    print(ref_data)
    ref_data = ref_data[:, ref_data.var.index.isin(target_data.var.index)]
    target_data.obs["cell_type_ref"] = "not defined"
    ref_data.obs[batch_col] = "ref"
    #for col in target_data.obs.columns[:-1]:
    #    ref_data.obs[col] = "ref"
    concat_data = sc.concat([target_data, ref_data])
    return concat_data