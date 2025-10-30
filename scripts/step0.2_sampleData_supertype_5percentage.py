
"""
Lingyi notes: supertype stratified sampling seems to work better than subclass stratified sampling. 5% much more stable than 25% sampling regarding finding MI/attractor features. 

Each sample contains approximately 5% of the original dataset's cells, 
and together all 20 samples reconstruct the complete original dataset without any overlap.
"""
import os
import numpy as np
import scanpy as sc
import pickle

def create_mutually_exclusive_samples(adata, n_samples=20, stratify_by=None, random_state=42):
    """
    Create mutually exclusive samples that together cover the entire dataset.
    Each sample will contain approximately 1/n_samples fraction of the data.
    """
    rng = np.random.default_rng(random_state)
    obs = adata.obs.copy()
    
    if stratify_by is None:
        # Simple random split
        indices = np.arange(adata.n_obs)
        rng.shuffle(indices)
        
        samples = []
        for i in range(n_samples):
            start_idx = int(i * len(indices) / n_samples)
            end_idx = int((i + 1) * len(indices) / n_samples)
            sample_indices = indices[start_idx:end_idx]
            samples.append(adata[sample_indices, :].copy())
        
        return samples
    
    # Stratified split
    if stratify_by not in obs.columns:
        raise ValueError(f"`{stratify_by}` not found in adata.obs")
    
    # Print statistics about the stratification groups
    group_sizes = obs.groupby(stratify_by).size()
    print(f"Stratifying by {stratify_by}:")
    print(f"  Total groups: {len(group_sizes)}")
    print(f"  Groups with < {n_samples} cells: {sum(group_sizes < n_samples)}")
    print(f"  Smallest group: {group_sizes.min()} cells")
    print(f"  Largest group: {group_sizes.max()} cells")
    print(f"  Median group size: {group_sizes.median():.1f} cells")
    
    samples = [[] for _ in range(n_samples)]
    
    # For each group, split into n_samples parts
    for group, group_df in obs.groupby(stratify_by, sort=False, dropna=False):
        group_indices = group_df.index.values
        rng.shuffle(group_indices)
        
        # If group has fewer cells than n_samples, distribute them across samples
        if len(group_indices) < n_samples:
            # Assign each cell to a different sample in round-robin 
            for idx, cell_idx in enumerate(group_indices):
                sample_idx = idx % n_samples
                samples[sample_idx].append(cell_idx)
        else:
            # Normal splitting for groups with enough cells
            for i in range(n_samples):
                start_idx = int(i * len(group_indices) / n_samples)
                end_idx = int((i + 1) * len(group_indices) / n_samples)
                if start_idx < len(group_indices):
                    sample_indices = group_indices[start_idx:end_idx]
                    samples[i].extend(sample_indices)
    
    # Convert to AnnData objects
    adata_samples = []
    for i, sample_indices in enumerate(samples):
        if sample_indices:  # Only create sample if it has cells
            sample_mask = adata.obs.index.isin(sample_indices)
            adata_samples.append(adata[sample_mask, :].copy())
        else:
            print(f"Warning: Sample {i} is empty")
    
    return adata_samples


# --------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------
data_path = "/burg/anastassiou/projects/seaad/data/"
sampled_data_path = "/burg/anastassiou/projects/seaad/data/sampledData"
# os.makedirs(sampled_data_path, exist_ok=True)

# --------------------------------------------------------------------
# A9 data
# --------------------------------------------------------------------
print("Loading A9 data...")
a9_data = sc.read_h5ad(os.path.join(data_path, "SEAAD_A9_RNAseq_DREAM.2025-07-15.h5ad"))
a9_data.obs["donor_supertype"] = a9_data.obs["Donor ID"].astype(str) + "_" + a9_data.obs["Supertype"].astype(str) # LC: subclass? superclass? superclass! 

print("Creating 20 mutually exclusive A9 samples...")
a9_samples = create_mutually_exclusive_samples(a9_data, n_samples=20, stratify_by="donor_supertype", random_state=42)

for i, sample in enumerate(a9_samples):
    sample.obs = sample.obs.drop(columns=["donor_supertype"])
    print(f"A9 sample {i+1:2d}: {sample.n_obs:6d} cells, {sample.n_vars} genes")
    
    # Save as pickle
    base_name = f"SEAAD_A9_RNAseq_DREAM.sample5_{i+1:02d}"
    pkl_path = os.path.join(sampled_data_path, base_name + ".pkl")
    # h5ad_path = os.path.join(sampled_data_path, base_name + ".h5ad")
    
    with open(pkl_path, "wb") as f:
        pickle.dump(sample, f)
    # sample.write_h5ad(h5ad_path)  LC: it seems that saving as h5ad has the "/" in the filename issue.

print(f"Total A9 cells across all samples: {sum(s.n_obs for s in a9_samples)} (original: {a9_data.n_obs})")

# Free up memory by deleting A9 data
del a9_data, a9_samples
import gc
gc.collect()
print("A9 data removed from memory to free up RAM")

# --------------------------------------------------------------------
# MTG data
# --------------------------------------------------------------------
print("\nLoading MTG data...")
mtg_data = sc.read_h5ad(os.path.join(data_path, "SEAAD_MTG_RNAseq_DREAM.2025-07-15.h5ad"))
mtg_data.obs["donor_supertype"] = mtg_data.obs["Donor ID"].astype(str) + "_" + mtg_data.obs["Supertype"].astype(str)

print("Creating 20 mutually exclusive MTG samples...")
mtg_samples = create_mutually_exclusive_samples(mtg_data, n_samples=20, stratify_by="donor_supertype", random_state=42)

for i, sample in enumerate(mtg_samples):
    sample.obs = sample.obs.drop(columns=["donor_supertype"])
    print(f"MTG sample {i+1:2d}: {sample.n_obs:6d} cells, {sample.n_vars} genes")
    
    # Save as pickle
    base_name = f"SEAAD_MTG_RNAseq_DREAM.sample5_{i+1:02d}"
    pkl_path = os.path.join(sampled_data_path, base_name + ".pkl")
    # h5ad_path = os.path.join(sampled_data_path, base_name + ".h5ad")
    
    with open(pkl_path, "wb") as f:
        pickle.dump(sample, f)
    # sample.write_h5ad(h5ad_path)

print(f"Total MTG cells across all samples: {sum(s.n_obs for s in mtg_samples)} (original: {mtg_data.n_obs})")

print("\nAll samples saved successfully!")
