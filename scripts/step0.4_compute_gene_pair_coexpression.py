import os
import re
import argparse
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
from joblib import Parallel, delayed
from sklearn.covariance import LedoitWolf
import matplotlib.pyplot as plt



# Gene filtering
def filter_genes_by_fraction(adata, min_frac=0.1):
    # keep genes expressed in at least min_frac of cells
    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    frac_expr = (X > 0).sum(axis=0) / X.shape[0]
    keep = np.asarray(frac_expr).ravel() >= min_frac
    adata_filtered = adata[:, keep].copy()
    return adata_filtered

# HVG selection (training pool only)
def select_hvgs(adata, n_hvg=1000, donor_subset=None):
    if donor_subset is not None:
        ad = adata[adata.obs['Region_Donor_ID'].isin(donor_subset)].copy()
    else:
        ad = adata.copy()
    expr = ad.X.toarray() if hasattr(ad.X, 'toarray') else ad.X
    gene_variances = np.var(expr, axis=0)
    top_gene_idx = np.argsort(-gene_variances)[:min(n_hvg, ad.n_vars)]
    hvgs = ad.var_names[top_gene_idx].tolist()
    return hvgs


def shrinkage_corr_matrix(Xd):
    """
    X: cells x genes array 
    returns: correlation matrix (genes x genes)
    """
    if Xd.shape[0] < 3:
        return None
    
    #Drop rows (cells) with any NaN
    if isinstance(Xd, np.ndarray):
        mask_valid = ~np.isnan(Xd).any(axis=1)
        Xd = Xd[mask_valid, :]
    
    # Ledoit-Wolf covariance shrinkage
    lw = LedoitWolf().fit(Xd)
    cov = lw.covariance_
    
    diag = np.sqrt(np.diag(cov))
    denom = np.outer(diag, diag)
    with np.errstate(divide='ignore', invalid='ignore'):
        corr = cov / denom
    corr[np.isnan(corr)] = 0.0
    np.fill_diagonal(corr, 1.0)
    return corr


def donor_shrinkage_corr_vector(adata, donor_id, triu_idx):
    """
    adata: AnnData (already filtered to hvgs)
    donor_id: Region_Donor_ID or Donor_ID (match adata.obs)
    triu_idx: np.triu_indices(n_genes, k=1)
    returns: 1D numpy vector length = n_pairs (upper triangle)
    """
    # select donor's cells
    mask = adata.obs['Region_Donor_ID'].values == donor_id
    if mask.sum() == 0:
        return np.full(triu_idx[0].shape[0], np.nan, dtype=np.float32)
    Xd = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    Xd = Xd[mask, :]
    
    corr_mat = shrinkage_corr_matrix(Xd)
    if corr_mat is None:
        return np.full(triu_idx[0].shape[0], np.nan, dtype=np.float32)
    return corr_mat[triu_idx].astype(np.float32)



# --- Load data ---

parser = argparse.ArgumentParser()
parser.add_argument("--Subclass", type=str)
args = parser.parse_args()

chosen_cell_type = args.Subclass
print(f"Running for cell type: {chosen_cell_type}")

celltype_sanitized = re.sub(r'(\d+)/(\d+)', r'\1-\2', chosen_cell_type)  
celltype_sanitized = re.sub(r'\s+', '_', celltype_sanitized)          

celltype_col = 'Subclass'
adata_all = sc.read_h5ad(f'data/SEAAD_A9_RNAseq_DREAM.2025-07-15.h5ad')
adata_a9 = adata_all[adata_all.obs[celltype_col] == chosen_cell_type].copy()
print(f'Cells in A9: {chosen_cell_type}: {adata_a9.n_obs}')

adata_all = sc.read_h5ad(f'data/SEAAD_MTG_RNAseq_DREAM.2025-07-15.h5ad')
adata_mtg = adata_all[adata_all.obs[celltype_col] == chosen_cell_type].copy()
print(f'Cells in MTG: {chosen_cell_type}: {adata_mtg.n_obs}')

adata_combined = sc.concat([adata_a9, adata_mtg], label="Region", keys=["A9", "MTG"])
adata_combined.obs['Region_Donor_ID'] = adata_combined.obs['Region'].astype(str) + "_" + adata_combined.obs['Donor ID'].astype(str)
adata_combined.obs.shape, adata_combined.n_vars


OUTDIR = f"results/pair_features/{celltype_sanitized}"
os.makedirs(OUTDIR, exist_ok=True)



# --- Preprocessing (select genes) ---

# Preprocess: filter cells and genes, select HVGs
print(f"\nBefore filtering: {adata_combined.n_obs} cells, {adata_combined.n_vars} genes.")
adata_filtered = filter_genes_by_fraction(adata_combined, min_frac=0.1)
print(f"After gene filtering: {adata_filtered.n_vars} genes remain.")

hvgs = select_hvgs(adata_filtered, n_hvg=1000) # using all donors and all cells
hvgs_file = os.path.join(OUTDIR, f"HVGs1000.txt")
pd.Series(hvgs).to_csv(hvgs_file, index=False, header=False)

# Select genes and prepare for correlation analysis
adata_selected = adata_filtered[:, hvgs].copy()
n_genes = len(hvgs)
triu_idx = np.triu_indices(n_genes, k=1)
n_pairs = len(triu_idx[0])
print(f"\nUsing {n_genes} genes for correlation analysis")
print(f"Number of gene pairs = {n_pairs:,}")

donors = np.unique(adata_selected.obs['Region_Donor_ID'].values)
print(f"Number of donors = {len(donors)}")

pair_names = [f"{hvgs[i]}_{hvgs[j]}" for i, j in zip(*triu_idx)]



# --- All cells ---
print("\n------ Using all cells -------")

print(f"\nProcessing all cells ...")
print(f"   Number of all cells: {adata_selected.n_obs}")

N_JOBS = 8
results = Parallel(n_jobs=N_JOBS)(
    delayed(donor_shrinkage_corr_vector)(adata_selected, donor_id, triu_idx) 
    for donor_id in donors
)

mmap = np.vstack(results).astype(np.float32)  # donors x pairs
X_donors_pairs = pd.DataFrame(mmap, index=donors, columns=pair_names)
print(f"   Shape of donor gene-pair features (correlation matrix): {X_donors_pairs.shape}")

# Check for donors with NaN values
donor_nan =np.isnan(X_donors_pairs).sum(axis=1)
if (donor_nan != 0).sum():
    print(f"   Donors with NaN values in gene-pair features: {X_donors_pairs.index[donor_nan != 0]}")
X_donors_pairs.fillna(0.0, inplace=True)  # fill NaN with 0

corr_file = os.path.join(OUTDIR, "pair_features_LW_all.csv")
X_donors_pairs.to_csv(corr_file)

corr_pkl = os.path.join(OUTDIR, "pair_features_LW_all.pkl")
with open(corr_pkl, "wb") as f:
    pickle.dump(X_donors_pairs, f)

print(f"   Saved gene-pair features to {corr_pkl}")



# --- 4 cell splits ---

print("\n------ Using 4 cell splits (20%) -------")

CELLID_DIR = "data/splitData_cellIDs/sampledData_supertype_25per"
os.makedirs(os.path.join(OUTDIR, "cell_splits_25per"), exist_ok=True)

N_JOBS = 8
pairs_df_list = []

for i in range(1, 5):

    print(f"\nProcessing split {i} ...")

    cell_splits_a9 = pd.read_csv(os.path.join(CELLID_DIR, f"SEAAD_A9_RNAseq_DREAM.sample25_{i:02d}_cellIDs.csv"), index_col=0)
    cell_splits_mtg = pd.read_csv(os.path.join(CELLID_DIR, f"SEAAD_MTG_RNAseq_DREAM.sample25_{i:02d}_cellIDs.csv"), index_col=0)
    adata_subset = adata_selected[adata_selected.obs.index.isin(cell_splits_a9.index) | adata_selected.obs.index.isin(cell_splits_mtg.index)].copy()
    print(f"   Number of cells in split {i}: {adata_subset.n_obs}")
    
    results = Parallel(n_jobs=N_JOBS)(
        delayed(donor_shrinkage_corr_vector)(adata_subset, donor_id, triu_idx)  
        for donor_id in donors
    )

    mmap = np.vstack(results).astype(np.float32)  # donors x pairs
    X_donors_pairs = pd.DataFrame(mmap, index=donors, columns=pair_names)
    print(f"   Shape of correlation matrix: {X_donors_pairs.shape}")

    # Check for donors with NaN values
    donor_nan =np.isnan(X_donors_pairs).sum(axis=1)
    if (donor_nan != 0).sum():
        print(f"   Donors with NaN values in gene-pair features: {X_donors_pairs.index[donor_nan != 0]}")
    X_donors_pairs.fillna(0.0, inplace=True)  # fill NaN with 0


    corr_file = os.path.join(OUTDIR, f"cell_splits_25per/pair_features_LW_{i:02d}.csv")
    X_donors_pairs.to_csv(corr_file)

    corr_pkl = os.path.join(OUTDIR, f"cell_splits_25per/pair_features_LW_{i:02d}.pkl")
    with open(corr_pkl, "wb") as f:
        pickle.dump(X_donors_pairs, f)

    pairs_df_list.append(X_donors_pairs)
    print(f"   Saved gene-pair features to {corr_pkl}")


# Aggregate results across splits (median)
pairs_concat = np.stack([df.values for df in pairs_df_list], axis=2)  # Stack: donors x pairs x 4
median_corr = np.nanmedian(pairs_concat, axis=2)  # take median across splits: donors x pairs

median_corr_df = pd.DataFrame(
    median_corr, 
    index=pairs_df_list[0].index, 
    columns=pairs_df_list[0].columns
)

median_corr_df.to_csv(os.path.join(OUTDIR, f"pair_features_LW_median.csv"))

corr_pkl = os.path.join(OUTDIR, f"pair_features_LW_median.pkl")
with open(corr_pkl, "wb") as f:
    pickle.dump(median_corr_df, f)

