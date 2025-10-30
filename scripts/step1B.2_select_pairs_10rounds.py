import os
import re
import pickle
import argparse
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


def compute_assoc(pairs_mat, y_train):
    
    rhos, pvals = [], []
    for i in range(pairs_mat.shape[1]):
        col = pairs_mat.iloc[:, i]
        if np.any(np.isnan(col)):
            rhos.append(np.nan)
            pvals.append(np.nan)
            continue
        rho, p = spearmanr(col, y_train)
        rhos.append(rho)
        pvals.append(p)
    
    pvals = np.array(pvals)
    rhos = np.array(rhos)

    # Adjust p-values
    valid_mask = ~np.isnan(pvals)
    qvals_full = np.full_like(pvals, np.nan)
    if valid_mask.sum() > 0:
        _, qvals, _, _ = multipletests(pvals[valid_mask], method="fdr_bh")
        qvals_full[valid_mask] = qvals

    results_df = pd.DataFrame({
        "feature": pairs_mat.columns,
        "gene1": pairs_mat.columns.str.split("_").str[0],
        "gene2": pairs_mat.columns.str.split("_").str[1],
        "rho": rhos,
        "pval": pvals,
        "qval": qvals_full
    })
    return results_df



# ----------- Load data -----------

parser = argparse.ArgumentParser()
parser.add_argument("--Subclass", type=str)
parser.add_argument("--Region", type=str)
args = parser.parse_args()

chosen_cell_type = args.Subclass
celltype_sanitized = re.sub(r'(\d+)/(\d+)', r'\1-\2', chosen_cell_type)  
celltype_sanitized = re.sub(r'\s+', '_', celltype_sanitized)   

chosen_region = args.Region
print(f"##### Running analysis for cell type: {chosen_cell_type}, region: {chosen_region} #####")

DATA_DIR = f"results/v7_pair_features/{celltype_sanitized}/"
os.makedirs(DATA_DIR, exist_ok=True)


# ----------- Seperate for selected region -----------

with open(os.path.join(DATA_DIR, f"pair_features_LW_all.pkl"), "rb") as f:
    pairs_all = pickle.load(f)

pairs_region_all = pairs_all[pairs_all.index.str.startswith(f"{chosen_region}_")]
pairs_region_all.index = pairs_region_all.index.str.replace(f"{chosen_region}_", "")
pairs_region_all.head()


with open(os.path.join(DATA_DIR, f"{chosen_region}_pair_features_LW_all.pkl"), "wb") as f:
    pickle.dump(pairs_region_all, f)


with open(os.path.join(DATA_DIR, f"pair_features_LW_median-25per.pkl"), "rb") as f:
    pairs_median = pickle.load(f)

pairs_region_median = pairs_median[pairs_median.index.str.startswith(f"{chosen_region}_")]
pairs_region_median.index = pairs_region_median.index.str.replace(f"{chosen_region}_", "")

with open(os.path.join(DATA_DIR, f"A9_pair_features_LW_median-25per.pkl"), "wb") as f:
    pickle.dump(pairs_region_median, f)


donors_meta = pd.read_csv(f"data/donors_meta_{chosen_region}.csv", index_col=0)
donors_df = pd.read_csv(f"data/{chosen_region}_donors_split.csv", index_col=0)

os.makedirs(os.path.join(DATA_DIR, f"{chosen_region}_donor_splits"), exist_ok=True)


# -----------  Pair features: median 25 per ------------

print("\n----------- Using pair features: median 25 per (4 cell splits) -----------")


with open(os.path.join(DATA_DIR, f"{chosen_region}_pair_features_LW_median-25per.pkl"), "rb") as f:
     pairs_df = pickle.load(f)

donor_removed = (pairs_df == 0).all(axis=1).sum()
if donor_removed > 0:   
    print(f"{chosen_cell_type}: remove {donor_removed} donors with all-zero (NaN) values.")
    pairs_df = pairs_df.loc[~(pairs_df == 0).all(axis=1)]
    donors_df = donors_df.loc[pairs_df.index]


N_ROUNDS = 10

for r in range(N_ROUNDS):

    print(f"\nProcessing round {r} ...")

    subset_df = donors_df[donors_df[f"train_r{r}"] == 1]                                 
    train_donors = subset_df.index
    print(f"Selected donors: {len(train_donors)}")
    
    for protein in ['6e10', 'AT8']:
        print(f"    Analyzing protein: {protein}")
        corr_mat = pairs_df.loc[train_donors]  # donors × pairs
        y_train = subset_df[f"percent {protein} positive area"]

        results_df = compute_assoc(corr_mat, y_train)
        results_df.to_csv(os.path.join(DATA_DIR, f"{chosen_region}_donor_splits/pairs_LW_median-25per_assoc_{protein}_r{r}.csv"), index=False)

