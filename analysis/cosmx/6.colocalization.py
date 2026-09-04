#!/usr/bin/env python3
# ==============================================================================
# Script: 6.colocalization.py
# Description: Spatial colocalization enrichment via KDTree and permutations
# ==============================================================================

import os
import gc
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import squidpy as sq
from sklearn.neighbors import KDTree
from scipy.stats import norm
from statsmodels.stats.multitest import fdrcorrection
import warnings

warnings.filterwarnings("ignore")

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
META_PATH = "/path/to/meta.csv"
OUT_DIR   = "/path/to/python_colocalization_dir/"

# Processing parameters
DISTANCES      = range(100, 1600, 50)
N_PERMUTATIONS = 1000
BATCH_SIZE     = 10_000

# Ensure output directory exists
os.makedirs(OUT_DIR, exist_ok=True)

# ==============================================================================
# Main Execution
# ==============================================================================
def main():
    print("Loading metadata...")
    df = pd.read_csv(META_PATH)
    samples = df['slide_fov'].unique()

    for smp in samples:
        print(f"\n{'='*50}\nProcessing sample: {smp}\n{'='*50}")
        start_time = time.time()

        df_sample = df[df["slide_fov"] == smp]
        coords = np.vstack([df_sample["CenterX_global_px"].values, df_sample["CenterY_global_px"].values]).T
        annotation = df_sample["refined"].values
        del df_sample

        tree = KDTree(coords)
        all_bins = []

        for radius in DISTANCES:
            print(f"\nStep 1: Querying neighbors at radius {radius}...")
            neighbors = tree.query_radius(coords, r=radius)
            #Delete self-connections
            neighbors = np.array([n[n != i] for i, n in enumerate(neighbors)], dtype=object)
            n_cells = len(neighbors)
            sparse_chunks = []

            print(f"Step 2: Building sparse matrix in batches (batch size: {BATCH_SIZE})...")
            for i_start in range(0, n_cells, BATCH_SIZE):
                i_end = min(i_start + BATCH_SIZE, n_cells)
                print(f"   - Processing batch: cells {i_start} to {i_end - 1}")

                neighbors_batch = neighbors[i_start:i_end]
                lengths = np.array([len(n) for n in neighbors_batch])
                row_idx = np.repeat(np.arange(i_start, i_end), lengths)
                col_idx = np.concatenate(neighbors_batch)
                data = np.ones_like(col_idx, dtype=np.uint8)

                chunk = sp.coo_matrix((data, (row_idx, col_idx)), shape=(n_cells, n_cells))
                sparse_chunks.append(chunk)

                # Clear processed neighbors to save memory
                for j in range(i_start, i_end):
                    neighbors[j] = None
                
                del neighbors_batch, row_idx, col_idx, data, lengths, chunk
                gc.collect()
            del neighbors

            print("Step 3: Merging batches into CSR matrix...")
            A_sparse = sparse_chunks[0].tocsr()
            sparse_chunks[0] = None
            gc.collect()

            for i in range(1, len(sparse_chunks)):
                print(f"   - Merging chunk {i + 1}/{len(sparse_chunks)}")
                A_sparse += sparse_chunks[i].tocsr()
                sparse_chunks[i] = None
                gc.collect()
            del sparse_chunks
            gc.collect()

            print("Step 4: Creating AnnData object...")
            adata = ad.AnnData(
                X=None,
                obs=pd.DataFrame({"refined": annotation}),
                obsm={"spatial": coords}
            )
            
            adata.obsp["spatial_connectivities"] = A_sparse
            del A_sparse
            adata.obs["refined"] = pd.Categorical(adata.obs["refined"])

            print("Step 5: Computing observed interactions...")
            sq.gr.interaction_matrix(adata, cluster_key="refined")
            anot_interactions = adata.uns['refined_interactions']
            idx = adata.obs["refined"].cat.categories
            
            df_matrix = pd.DataFrame(anot_interactions, index=idx, columns=idx)
            observed = df_matrix.reset_index().melt(id_vars='index')
            observed.columns = ['from', 'to', 'count']
            observed = observed[observed['count'] > 0].reset_index(drop=True)
            observed['bin'] = radius
            del df_matrix

            print(f"Step 6: Running simulations ({N_PERMUTATIONS} permutations)...")
            all_results = []
            for i in range(N_PERMUTATIONS):
                if i % 10 == 0:
                    print(f"   - Simulation {i}/{N_PERMUTATIONS}")
                
                np.random.seed(1234 + i)
                adata.obs['refined_random'] = pd.Categorical(np.random.permutation(adata.obs['refined'].values))
                sq.gr.interaction_matrix(adata, cluster_key="refined_random")

                anot_interactions = adata.uns['refined_random_interactions']
                df_matrix = pd.DataFrame(anot_interactions, index=idx, columns=idx)
                df_long = df_matrix.reset_index().melt(id_vars='index')
                df_long.columns = ['from', 'to', 'count']
                df_long = df_long[df_long['count'] > 0].reset_index(drop=True)
                df_long['bin'] = radius
                df_long['sim'] = i
                all_results.append(df_long)
                
                del df_matrix, df_long

            results = pd.concat(all_results, ignore_index=True)
            del all_results

            print("Step 7: Aggregating simulation results...")
            simulations_agg = (
                results.groupby(['from', 'to', 'bin'])['count']
                .agg(mean_count='mean', std_count='std')
                .reset_index()
            )
            del results

            print("Step 8: Merging observed with simulations...")
            final_df = pd.merge(
                observed,
                simulations_agg,
                on=['from', 'to', 'bin'],
                how="left"
            )
            del observed, simulations_agg

            print("Step 9: Cleanup and enrichment calculation...")
            final_df['mean_count'] = final_df['mean_count'].fillna(0)
            final_df['std_count']  = final_df['std_count'].fillna(0)
            final_df = final_df[~((final_df['count'] == 0) & (final_df['mean_count'] == 0))]
            final_df['enrichment'] = np.log2((final_df['count'] + 1) / (final_df['mean_count'] + 1))
            
            stats_df = final_df[final_df['std_count'] > 0].copy()
            del final_df

            print("Step 10: Computing z-scores and p-values...")
            stats_df['zscore'] = (stats_df['count'] - stats_df['mean_count']) / stats_df['std_count']
            stats_df['pvalue'] = 2 * norm.cdf(-abs(stats_df['zscore']))
            stats_df['padj']   = fdrcorrection(stats_df['pvalue'])[1]
            stats_df['sample'] = smp
            all_bins.append(stats_df)

            del adata, stats_df
            gc.collect()

        print(f"Step 11: Saving results for {smp}...")
        final_sample = pd.concat(all_bins, ignore_index=True)
        final_sample.to_csv(os.path.join(OUT_DIR, f"{smp}.csv"), index=False)
        
        del all_bins, final_sample
        gc.collect()

        total_time = time.time() - start_time
        print(f"Sample {smp} completed in {total_time:.2f} seconds ({total_time/60:.2f} minutes).")

if __name__ == "__main__":
    main()
