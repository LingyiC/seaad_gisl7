import anndata as ad
import numpy as np

from sklearn.model_selection import train_test_split


def get_donor_meta(adata):
    order_maps = {
        "ADNC":   {"Not AD": 0, "Low": 1, "Intermediate": 2, "High": 3},
        "Braak":  {f"Braak {i}": idx for idx, i in enumerate(["0","I","II","III","IV","V","VI"])},  
        "Thal":   {f"Thal {i}": i for i in range(6)},                             
        "CERAD":  {"Absent": 0, "Sparse": 1, "Moderate": 2, "Frequent": 3}
    }

    cell_cols = ['library_prep', 'Method', 'Class', 'Subclass', 'Supertype', 'n_cells_aggregated']
    donor_cols = [c for c in adata.obs.columns if c not in cell_cols]

    agg_dict = {c: (lambda x: x.iloc[0]) for c in donor_cols}

    meta_donor = (
        adata.obs[donor_cols]
        .groupby("Donor ID")
        .agg(agg_dict)
    )

    for col, mapper in order_maps.items():
        if col in meta_donor.columns:
            meta_donor[f"{col}_ord"] = meta_donor[col].astype(str).map(mapper).astype(int)

    return meta_donor



def generate_train_val_sets(meta_donors, random_seed, n_rounds=10):
    output_df = meta_donors.copy()
    
    all_rounds = []
    for r in range(n_rounds):

        train_idx, val_idx = train_test_split(
            output_df.index,
            test_size=0.5,
            stratify=output_df["ADNC"],
            random_state=random_seed + r
        )

        print(f"Round {r}: train={len(train_idx)} val={len(val_idx)} donors")

        all_rounds.append({
            "round": r,
            "train_idx": np.array(train_idx),
            "val_idx": np.array(val_idx)
        })

        output_df[f"train_r{r}"] = output_df.index.isin(train_idx).astype(int)
        output_df[f"val_r{r}"] = output_df.index.isin(val_idx).astype(int)

    return output_df



for chosen_region in ['A9', 'MTG']:

    adata = ad.read_h5ad(f"data/SEAAD_{chosen_region}_RNAseq_DREAM.2025-07-15.h5ad")
    meta_donors = get_donor_meta(adata)
    meta_donors.to_csv(f"data/donors_meta_{chosen_region}.csv", index=True)

    #meta_donors = pd.read_csv(f"data/donors_meta_{chosen_region}.csv", index_col=0)
    output_df = generate_train_val_sets(meta_donors, random_seed=42, n_rounds=10)
    output_df.to_csv(f"data/{chosen_region}_donors_split.csv")




