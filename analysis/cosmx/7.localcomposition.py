#!/usr/bin/env python3
# ==============================================================================
# Script: 7.localcomposition.py
# Description: Spatial neighborhood enrichment and local composition analysis
# ==============================================================================

import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import squidpy as sq
import seaborn as sns
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
META_PATH = "/path/to/your/meta.csv"
OUT_DIR   = "/path/to/your/output/directory/"
PLOT_OUT_PATH = os.path.join(OUT_DIR, "neighborhood_enrichment_subset.png")

# Processing Parameters
N_NEIGHS = 60

# ------------------------------------------------------------------------------
# Define custom cell groupings and targets here
# ------------------------------------------------------------------------------
# List of specific subtypes you want to group into a broader category
SUBTYPES_TO_GROUP = [
    'Subtype_1', 'Subtype_2', 'Subtype_3', 'Subtype_4'
]
GROUP_NAME = 'Grouped_Major_Type'

# Define which cell types to visualize on the Heatmap (Rows vs Columns)
TARGET_ROWS = ["Cell_Type_A", "Cell_Type_B"]
TARGET_COLS = ["Cell_Type_C", "Cell_Type_D", "Cell_Type_E", "Grouped_Major_Type"]

# Ensure output directory exists
os.makedirs(OUT_DIR, exist_ok=True)

# ==============================================================================
# Main Execution
# ==============================================================================
def main():
    print("Loading metadata...")
    meta = pd.read_csv(META_PATH)

    # Group specific subtypes into a broad category
    meta['refined'] = meta['refined'].replace(SUBTYPES_TO_GROUP, GROUP_NAME)

    # Create dummy matrix for AnnData (Squidpy spatial functions require an X matrix)
    X_dummy = np.zeros((len(meta), 1)) 
    adata = ad.AnnData(X=X_dummy, obs=meta)

    # Assign spatial coordinates
    coords = meta[["CenterX_global_px", "CenterY_global_px"]].to_numpy()
    adata.obsm["spatial"] = coords
    
    # Ensure correct data types for categorical grouping
    adata.obs["refined"] = adata.obs["refined"].astype("category")
    adata.obs["fov"] = adata.obs["fov"].astype("category")

    # Calculate spatial neighbors
    print(f"Calculating spatial neighbors (n_neighs={N_NEIGHS})...")
    sq.gr.spatial_neighbors(
        adata,
        spatial_key="spatial", 
        coord_type="generic", 
        n_neighs=N_NEIGHS,
        library_key="fov"
    )

    # Calculate neighborhood enrichment
    print("Calculating neighborhood enrichment...")
    sq.gr.nhood_enrichment(
        adata, 
        cluster_key="refined"
    )

    # ==============================================================================
    # Visualization
    # ==============================================================================
    print("Generating heatmap...")
    zscores = adata.uns["refined_nhood_enrichment"]["zscore"]
    cell_types = adata.obs["refined"].cat.categories

    # Convert z-scores into a Pandas DataFrame for easy filtering
    enrichment_df = pd.DataFrame(zscores, index=cell_types, columns=cell_types)

    # Subset the DataFrame based on target rows and columns
    subset_df = enrichment_df.loc[TARGET_ROWS, TARGET_COLS]

    # Plot configuration
    plt.figure(figsize=(10, 4))
    sns.heatmap(
        subset_df, 
        cmap="coolwarm", 
        center=0, 
        annot=True, 
        fmt=".2f", 
        cbar_kws={'label': 'Z-score'}
    )

    plt.title("Neighborhood Enrichment Interaction Matrix")
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()

    # Save output
    print(f"Saving plot to {PLOT_OUT_PATH}...")
    plt.savefig(PLOT_OUT_PATH, dpi=600, bbox_inches='tight')
    print("Analysis complete.")

if __name__ == "__main__":
    main()
