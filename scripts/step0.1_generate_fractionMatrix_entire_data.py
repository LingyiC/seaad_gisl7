
import pandas as pd
import numpy as np
import os
import anndata
from scipy import sparse

# Lingyi 1007: Added memory-efficient processing for gene fractions, by chunk processing


# ============================================================
# Set up paths
# ============================================================
output_folder = "/burg/anastassiou/projects/seaad/data/csv_fraction_all"
os.makedirs(output_folder, exist_ok=True)

# Only one file (MTG). Change to A9 file for A9 processing
adata_file = "/burg/anastassiou/projects/seaad/data/SEAAD_MTG_RNAseq_DREAM.2025-07-15.h5ad"

print(f"Processing single AnnData file: {os.path.basename(adata_file)}")

# ============================================================
# Process MTG file
# ============================================================
try:
    adata = anndata.read_h5ad(adata_file)
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
    # Gene expression fractions, memory-efficient chunk processing, otherwise, too large for Ginsburg processing 
    # -------------------------------------------------
    gene_names = adata.var.index.tolist()
    print(f"\nTotal number of genes: {len(gene_names)}")

    print("Calculating gene expression fractions for all genes using sparse matrices...")

    # Get unique donors (convert to string to avoid categorical issues)
    unique_donors = adata.obs['Donor ID'].astype(str).unique()
    print(f"Number of unique donors: {len(unique_donors)}")

    # Initialize results dictionary
    gene_fractions_dict = {'Donor ID': unique_donors}

    # Process genes in chunks to manage memory, process 1000 genes at a time
    chunk_size = 1000
    n_genes = len(gene_names)
    n_chunks = (n_genes + chunk_size - 1) // chunk_size
    
    print(f"Processing {n_genes} genes in {n_chunks} chunks of {chunk_size}")
    
    for chunk_idx in range(n_chunks):
        start_idx = chunk_idx * chunk_size
        end_idx = min((chunk_idx + 1) * chunk_size, n_genes)
        chunk_genes = gene_names[start_idx:end_idx]
        
        print(f"Processing chunk {chunk_idx + 1}/{n_chunks}: genes {start_idx} to {end_idx-1}")
        
        # Extract expression data for this chunk
        if hasattr(adata.X, 'toarray'):  # sparse matrix
            chunk_expression = adata.X[:, start_idx:end_idx]
            # Create binary expression matrix (genes expressed > 0)
            chunk_binary = (chunk_expression > 0).astype(np.float32)
        else:  # dense numpy array
            chunk_expression = adata.X[:, start_idx:end_idx]
            chunk_binary = (chunk_expression > 0).astype(np.float32)
        
        # Convert to DataFrame for groupby operation
        if sparse.issparse(chunk_binary):
            chunk_df = pd.DataFrame(
                chunk_binary.toarray(),
                index=adata.obs.index,
                columns=chunk_genes
            )
        else:
            chunk_df = pd.DataFrame(
                chunk_binary,
                index=adata.obs.index,
                columns=chunk_genes
            )
        
        # Add donor ID (convert to string to avoid categorical issues)
        chunk_df['Donor ID'] = adata.obs['Donor ID'].astype(str).values
        
        # Calculate fractions for this chunk
        chunk_fractions = chunk_df.groupby('Donor ID', observed=True).mean().reset_index()
        
        # Add to results dictionary
        for gene in chunk_genes:
            donor_fractions = pd.DataFrame({'Donor ID': unique_donors})
            donor_fractions = donor_fractions.merge(
                chunk_fractions[['Donor ID', gene]], 
                on='Donor ID', 
                how='left'
            ).fillna(0)
            gene_fractions_dict[gene] = donor_fractions[gene].values
        
        # Clean up memory
        del chunk_df, chunk_fractions, chunk_binary
        if 'chunk_expression' in locals():
            del chunk_expression

    # Create final gene fractions DataFrame
    print("Creating final gene fractions DataFrame...")
    gene_fractions = pd.DataFrame(gene_fractions_dict)

    print(f"Gene fractions shape: {gene_fractions.shape}")
    print("Sample of gene fractions (first 5 donors, first 5 genes):")
    print(gene_fractions.iloc[:5, :6])

    # -------------------------------------------------
    # Combine donor metadata + gene fractions
    # -------------------------------------------------
    adf = adata.obs.copy()

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
    input_filename = os.path.basename(adata_file)
    output_filename = input_filename.replace('.h5ad', '_fractions.csv')
    output_file = os.path.join(output_folder, output_filename)

    print(f"Saving to: {output_file}")
    combined_sorted.to_csv(output_file, index=False)
    print(f"Successfully saved {output_filename}")

    print("Sample of final combined data:")
    print(combined_sorted.head())

except Exception as e:
    print(f"Error processing {os.path.basename(adata_file)}: {str(e)}")

print(f"\n{'='*60}")
print("Processing complete!")
print(f"Output file saved to: {output_folder}")
