import pickle
import pandas as pd
import numpy as np
import os
import glob

# ============================================================
# Set up paths
# ============================================================
input_folder = "/burg/anastassiou/projects/seaad/data/sampledData_supertype_5per/pickles"
output_folder = "/burg/anastassiou/projects/seaad/data/sampledData_supertype_5per/csv_fraction"

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Get all pickle files in the folder
pickle_files = glob.glob(os.path.join(input_folder, "*.pkl"))
pickle_files.sort()  # Sort for consistent processing order

print(f"Found {len(pickle_files)} pickle files to process:")
for file in pickle_files:
    print(f"  - {os.path.basename(file)}")

# ============================================================
# Process each pickle file
# ============================================================
for i, file_path in enumerate(pickle_files):
    print(f"\n{'='*60}")
    print(f"Processing file {i+1}/{len(pickle_files)}: {os.path.basename(file_path)}")
    print(f"{'='*60}")

    try:
        # Load the data
        with open(file_path, "rb") as f:
            adata = pickle.load(f)

        print("Data loaded successfully!")
        print(f"AnnData shape: {adata.shape}")
        print(f"Available columns in obs: {list(adata.obs.columns)}")

        # -------------------------------------------------
        # Donor metadata
        # -------------------------------------------------
        donor_columns = ['Donor ID', 'Sex', 'Age at Death', 'ADNC', 'Braak', 'CERAD', 'Cognitive Status']
        donor_columns = [col for col in donor_columns if col in adata.obs.columns]

        donor_info = adata.obs[donor_columns].drop_duplicates()
        print(f"\nNumber of unique donors: {len(donor_info)}")

        # -------------------------------------------------
        # Cell-level donor IDs
        # -------------------------------------------------
        cell_data = pd.DataFrame(adata.obs[['Donor ID']])
        cell_data.index = adata.obs.index
        print(f"\nCell data shape: {cell_data.shape}")

        # -------------------------------------------------
        # Gene expression fractions
        # -------------------------------------------------
        gene_names = adata.var.index.tolist()
        print(f"\nTotal number of genes: {len(gene_names)}")

        print("Calculating gene expression fractions for all genes...")

        # Create binary expression matrix (genes expressed > 0)
        if hasattr(adata.X, 'toarray'):  # sparse matrix
            expression_binary = (adata.X > 0).astype(int)
            try:
                expression_dense = expression_binary.toarray()
            except MemoryError:
                raise RuntimeError("Matrix too large to convert to dense. Consider using sparse methods.")
        else:  # dense numpy array
            expression_dense = (adata.X > 0).astype(int)

        # DataFrame with expression presence/absence
        print("Creating expression DataFrame...")
        expression_df = pd.DataFrame(
            expression_dense,
            index=adata.obs.index,
            columns=adata.var.index
        )

        # Add donor ID
        expression_df['Donor ID'] = adata.obs['Donor ID'].values

        # Fractions per donor
        print("Computing fractions per donor...")
        gene_fractions = expression_df.groupby('Donor ID').mean().reset_index()

        print(f"Gene fractions shape: {gene_fractions.shape}")
        print("Sample of gene fractions (first 5 donors, first 5 genes):")
        print(gene_fractions.iloc[:5, :6])

        # -------------------------------------------------
        # Combine donor metadata + gene fractions
        # -------------------------------------------------
        adf = adata.obs.copy()
        # modified on 10/7, adding'Hispanic/Latino',
        metadata_columns = [
            'Donor ID', 'Sex', 'Age at Death', 'Years of education', 'PMI', 'APOE Genotype', 'Hispanic/Latino',
            'Thal', 'Braak', 'CERAD', 'ADNC', 'Cognitive Status', 'LATE', 'Highest Lewy Body Disease',
            'percent 6e10 positive area', 'percent AT8 positive area',
            'percent NeuN positive area', 'percent GFAP positive area',
            'percent aSyn positive area', 'percent pTDP43 positive area'
        ]
        metadata_columns = [col for col in metadata_columns if col in adf.columns]

        adf_unique = adf[metadata_columns].drop_duplicates()

        combined = pd.merge(adf_unique, gene_fractions, on='Donor ID', how='inner')
        combined_sorted = combined.sort_values('Donor ID').reset_index(drop=True)

        print(f"Combined shape: {combined_sorted.shape}")

        # -------------------------------------------------
        # Save to CSV
        # -------------------------------------------------
        input_filename = os.path.basename(file_path)
        output_filename = input_filename.replace('.pkl', '_fractions.csv')
        output_file = os.path.join(output_folder, output_filename)

        print(f"Saving to: {output_file}")
        combined_sorted.to_csv(output_file, index=False)
        print(f"Successfully saved {output_filename}")

        print("Sample of final combined data:")
        print(combined_sorted.head())

    except Exception as e:
        print(f"Error processing {os.path.basename(file_path)}: {str(e)}")
        continue

print(f"\n{'='*60}")
print("Processing complete!")
print(f"Output files saved to: {output_folder}")




####
# sanity check
####

# import os
# import pickle
# import pandas as pd
# file_path = "/burg/anastassiou/projects/seaad/data/sampledData/pickles/SEAAD_MTG_RNAseq_DREAM.sample5_11.pkl"
# with open(file_path, "rb") as f:
#     adata = pickle.load(f)
# gs = "SPP1"
# adf = adata.obs.copy()
# adf[gs] = adata[:, gs].X.toarray().flatten() > 0
# frac_ = adf.groupby("Donor ID", observed=False)[gs].mean().reset_index() #.rename(columns={"CD163p": "FrCD163p"})


# csv_path = "/Users/lingyi/Documents/sea-ad/csv_fraction/SEAAD_MTG_RNAseq_DREAM.sample5_11_fractions.csv"
# df = pd.read_csv(csv_path)
 
# df.loc[df['Donor ID'] == "H21.33.046", "COL1A1"]


# csv_path = "/burg/anastassiou/projects/seaad/data/sampledData/csv_fraction/SEAAD_MTG_RNAseq_DREAM.sample5_11_fractions.csv"
# df = pd.read_csv(csv_path)
 
# df.loc[df['Donor ID'] == "H21.33.046", "SPP1"]