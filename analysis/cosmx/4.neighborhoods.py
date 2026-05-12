#!/usr/bin/env python3
# ==============================================================================
# Script: 4.neighborhoods.py
# Description: Spatial neighborhood composition and clustering using KDTree
# ==============================================================================

import os
import time
import random
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from pathlib import Path
from sklearn.neighbors import KDTree
import warnings

warnings.filterwarnings("ignore")

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
META_PATH = "/path/to/Post_analysis/Neighborhood/meta.csv"
OUT_BASE_DIR = "/path/to/Post_analysis/Neighborhood/Run/"

# Parameters
KNN_LIST = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95]
RESOLUTIONS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
TIER_FINE_COL = 'refined'   # Tier 2/3
TIER_BROAD_COL = 'subset'   # Tier 1
SCALING_FACTOR = 0.5

# Ensure reproducibility
random.seed(46)
np.random.seed(46)

# ==============================================================================
# Helper Functions
# ==============================================================================
def calculate_composition_matrix(neighbor_identities, knn, num_classes):
    """Applies distance-based weights, calculates frequency histograms, and adds noise."""
    weight_vector = np.linspace(knn, 1, num=knn, dtype=int)
    weighted_neighbors = np.repeat(neighbor_identities, weight_vector, axis=1)
    
    def hist_func(x):
        return np.histogram(x, bins=range(num_classes + 1))[0]
    
    comp_matrix = np.apply_along_axis(hist_func, 1, weighted_neighbors)
    comp_matrix = comp_matrix / knn  # Normalize
    
    # Add noise
    noise = np.random.uniform(size=comp_matrix.shape, low=0.01, high=0.1)
    return comp_matrix + noise

# ==============================================================================
# Main Execution
# ==============================================================================
def main():
    print("Loading metadata...")
    meta = pd.read_csv(META_PATH)
    slice_ids = meta['sample'].unique()

    # Extract unique labels for both tiers
    unique_fine_labels = meta[TIER_FINE_COL].unique()
    unique_broad_labels = meta[TIER_BROAD_COL].unique()

    # Create conversion dictionaries
    fine_name_to_idx = {name: idx for idx, name in enumerate(unique_fine_labels)}
    broad_name_to_idx = {name: idx for idx, name in enumerate(unique_broad_labels)}

    # Apply mappings
    meta['fine_index'] = meta[TIER_FINE_COL].apply(lambda x: fine_name_to_idx[x])
    meta['broad_index'] = meta[TIER_BROAD_COL].apply(lambda x: broad_name_to_idx[x])
    
    n_cells = meta.shape[0]

    # Iterate over each kNN parameter
    for knn in KNN_LIST:
        print(f"\n{'='*50}\nStarting Analysis for kNN = {knn}\n{'='*50}")
        out_dir = Path(OUT_BASE_DIR) / str(knn)
        out_dir.mkdir(parents=True, exist_ok=True)

        # Pre-allocate arrays for this kNN
        neighbor_id_fine = np.zeros((n_cells, knn), dtype=int)
        neighbor_id_broad = np.zeros((n_cells, knn), dtype=int)
        
        # We process coordinates slice by slice
        array_pos = 0
        for slice_id in slice_ids:
            local_t = time.time()
            
            meta_slice = meta[meta['sample'] == slice_id]
            num_cells = len(meta_slice)
            
            # Spatial coordinates
            coords = meta_slice[['CenterX_global_px', 'CenterY_global_px']].to_numpy()
            
            # KDTree Spatial Search (Computed ONCE per slice)
            tree = KDTree(coords)
            dist, ind = tree.query(coords, k=(knn + 1))
            ind = ind[:, 1:]  # Drop the first entry (self)
            
            # Map Fine Identities (Tier 2/3)
            local_fine_idx = meta_slice['fine_index'].to_numpy()
            neighbor_id_fine[array_pos:(array_pos + num_cells), :] = local_fine_idx[ind]
            
            # Map Broad Identities (Tier 1)
            local_broad_idx = meta_slice['broad_index'].to_numpy()
            neighbor_id_broad[array_pos:(array_pos + num_cells), :] = local_broad_idx[ind]

            array_pos += num_cells
            print(f"Processed slice {slice_id} (n={num_cells}) in {time.time() - local_t:.2f}s")

        # Compute composition matrices
        print("Calculating weighted composition matrices...")
        comp_matrix_fine = calculate_composition_matrix(neighbor_id_fine, knn, len(unique_fine_labels))
        comp_matrix_broad = calculate_composition_matrix(neighbor_id_broad, knn, len(unique_broad_labels))

        # Convert to DataFrames
        df_fine = pd.DataFrame(comp_matrix_fine, index=meta['cell_names'], columns=unique_fine_labels)
        
        broad_col_names = [f"{label}_Tier1" for label in unique_broad_labels]
        df_broad = pd.DataFrame(comp_matrix_broad, index=meta['cell_names'], columns=broad_col_names)

        # Merge Tier 1 (scaled) and Tier 2/3
        df_merged = pd.concat([(df_broad * SCALING_FACTOR), df_fine], axis=1)
        df_merged.to_csv(out_dir / "weight_mat.csv")

        # ----------------------------------------------------------------------
        # Scanpy / AnnData Processing
        # ----------------------------------------------------------------------
        print("Creating AnnData object and running Scanpy pipeline...")
        adata = ad.AnnData(df_merged)

        local_t = time.time()
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        print(f"PCA, KNN, and UMAP completed in {time.time() - local_t:.2f}s")

        # Leiden Clustering across resolutions
        for res in RESOLUTIONS:
            local_t = time.time()
            sc.tl.leiden(adata, key_added=f"leiden_res{res}", resolution=res, flavor="igraph", n_iterations=2, directed=False)
            print(f"Leiden clustering (res={res}) completed in {time.time() - local_t:.2f}s")

        # ----------------------------------------------------------------------
        # Save Final Outputs
        # ----------------------------------------------------------------------
        print("Saving outputs...")
        
        # Save H5AD
        adata.write(out_dir / 'andata.h5ad')
        
        # Save Dataframe
        dataframe = adata.obs.copy()
        dataframe.to_csv(out_dir / 'dataframe.csv', index=True)
        
        # Save UMAP coordinates
        umap_df = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names, columns=['UMAP1', 'UMAP2'])
        umap_df.to_csv(out_dir / 'umap.csv', index=True)
        
        # Merge with metadata and save
        dataframe['cell_names'] = dataframe.index
        meta_merged = meta.merge(dataframe, on='cell_names', how='left')
        meta_merged.to_csv(out_dir / 'seu_neigh.csv', index=False)

        print(f"Finished kNN = {knn} neighborhood analysis.")

if __name__ == "__main__":
    main()
