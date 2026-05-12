#!/usr/bin/env python3
# ==============================================================================
# Script: scotia_run.py
# Description: Execute SCOTIA Optimal Transport for cell interactions
# ==============================================================================

import os
import sys
import warnings
import numpy as np
import pandas as pd
import scotia
from scipy.spatial import distance_matrix

warnings.filterwarnings("ignore")

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
BASE_DIR = "/path/to/output/SCOTIA/Files/"
LR_PAIR_FILE = os.path.join(BASE_DIR, "lr_pair.csv")

# List of patients to process
SAMPLES = ["Slide_1", "Slide_2", "Slide_3", "Slide_healthy"]

def read_data(file):
    try:
        return pd.read_csv(file, header=0, index_col=None, sep="\t")
    except Exception as e:
        print(f"Error reading {file}: {e}")
        return pd.DataFrame()

def main():
    # Load known Ligand-Receptor pairs
    known_lr_pairs = pd.read_csv(LR_PAIR_FILE, sep=',', index_col=None)
    print(f"Loaded {len(known_lr_pairs)} LR pairs.")

    for patient in SAMPLES:
        print(f"Processing patient: {patient}")
        meta_path = os.path.join(BASE_DIR, f"{patient}_meta_def.csv")
        exp_path = os.path.join(BASE_DIR, f"{patient}_exp_def.csv")

        meta_df = pd.read_csv(meta_path, sep=',', index_col=None)
        exp_df = pd.read_csv(exp_path, sep=',', index_col=None)
        
        fovs = meta_df.iloc[:, 1].unique()

        # Step 1: DBSCAN cell clustering per FOV
        for fov in fovs:
            print(f"  Clustering FOV: {fov}")
            celltype = []
            cell_idx = []
            
            meta_df_fov = meta_df[meta_df['fov'] == fov].copy()
            meta_df_fov['index'] = range(meta_df_fov.shape[0])
            
            # Save FOV subset
            fov_csv_out = os.path.join(BASE_DIR, f"{patient}_fov_{fov}.csv")
            meta_df_fov.to_csv(fov_csv_out, header=True, index=False, sep="\t")

            cell_type_l = list(set(meta_df_fov['annotation']))
            
            for ct in cell_type_l:
                meta_df_sel = meta_df_fov[meta_df_fov['annotation'] == ct]
                X = np.array(meta_df_sel[['x_positions', 'y_positions']])
                
                fov_size_x = meta_df_fov['x_positions'].max() - meta_df_fov['x_positions'].min()
                fov_size_y = meta_df_fov['y_positions'].max() - meta_df_fov['y_positions'].min()
                fov_size = np.max([fov_size_x, fov_size_y])
                
                if X.shape[0] >= 5:
                    idx_l, _ = scotia.dbscan_ff_cell(
                        X, 
                        X_index_arr=np.array(meta_df_sel['index']),
                        min_cluster_size=5, 
                        eps_l=list(range(10, int(fov_size), 1))
                    )
                    
                    if len(idx_l) > 0:
                        celltype += [ct for _ in idx_l]
                        cell_idx += idx_l
                        
            tmp_df = pd.DataFrame([celltype, cell_idx]).T
            cluster_out = os.path.join(BASE_DIR, f"{patient}_fov_{fov}_dbscan.cell.clusters.test.npy")
            np.save(cluster_out, tmp_df)

        # Step 2: Optimal Transport
        exp_df_norm = exp_df.iloc[:, 2:]
        exp_df_norm = exp_df_norm[exp_df_norm > 0]
        df_quantile = exp_df_norm.quantile(q=0.99, axis=0) 

        for fov in fovs:
            print(f"  Calculating OT for FOV: {fov}")
            
            # Load clustering result
            cluster_out = os.path.join(BASE_DIR, f"{patient}_fov_{fov}_dbscan.cell.clusters.test.npy")
            cluster_df = pd.DataFrame(np.load(cluster_out, allow_pickle=True))
            cluster_df.columns = ['cell_type', 'cell_idx']
            
            # Coordinates
            meta_df_fov = meta_df[meta_df['fov'] == fov].copy()
            meta_df_fov.index = range(meta_df_fov.shape[0])
            cell_id_all = np.array(range(meta_df_fov.shape[0]))
            
            coord = np.array(meta_df_fov[['x_positions', 'y_positions']])
            S_all_arr = distance_matrix(coord, coord)
            
            # Expression
            exp_df_fov = exp_df[exp_df['fov'] == fov].iloc[:, 2:]
            exp_df_fov = exp_df_fov / df_quantile
            exp_df_fov[exp_df_fov > 1] = 1
            exp_df_fov.index = cell_id_all

            # Select potentially communicating cell cluster pairs
            S_all_arr_new = scotia.sel_pot_inter_cluster_pairs(S_all_arr, cluster_df)

            # Optimal transport between source and target cells
            ga_df_final = scotia.source_target_ot(S_all_arr_new, exp_df_fov, meta_df_fov, known_lr_pairs)
            
            if ga_df_final.shape[0] > 0:
                ga_df_final.columns = ['source_cell_idx', 'receptor_cell_idx', 'likelihood', 'ligand_recptor', 'source_cell_type', 'target_cell_type']
                ga_df_final['cell_pairs'] = ga_df_final['source_cell_type'] + "_" + ga_df_final['target_cell_type']
                
                ot_out = os.path.join(BASE_DIR, f"{patient}_fov_{fov}.ot.csv")
                ga_df_final.to_csv(ot_out, header=True, index=False, sep="\t")

if __name__ == "__main__":
    main()
