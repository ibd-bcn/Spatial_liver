


# prior to this convert seurat obj todas to ann data objet 
import anndata as ad
import scanpy as sc
import liana as li
import numpy as np


# Perform liana rank aggregation for liana results 
# per each condition an andata object can be created and run as the following
liana_res_group = li.mt.rank_aggregate(
    adata,
    groupby='annotation',
    resource_name='consensus',
    expr_prop=0.1,
    verbose=True,
    inplace=False  
)


print(liana_res_group.head())

# tensor building 

li.mt.rank_aggregate.by_sample(
    adata,
    groupby='annotation',
    resource_name='consensus',
    sample_key='sample_key'sample, 
    use_raw=False, 
    verbose=True, 
    n_perms=None, 
    return_all_lrs=True, 
    )
    
tensor = li.multi.to_tensor_c2c(adata,
                                sample_key='sample_key'sample,
                                score_key='magnitude_rank', 
                                how='outer_cells' 
                                )
                                
c2c.io.export_variable_with_pickle(tensor, "tensor.pkl")
