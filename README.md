# README

This repository contains the workflow and scripts used to generate fraction matrices, perform subsampling, identify MI-ranked features, identify gene pairs, and train predictive models for **6e10** and **AT8** protein in SEA-AD datasets.



## **Workflow Overview**

The analysis pipeline consists of the following steps:

### **Step 0: Data Preparation**

1. **`step0.1_generate_fractionMatrix_entire_data.py`**      
   Generates the fraction matrix for the entire MTG and A9 datasets, respectively.

2. **`step0.2_sampleData_supertype_5percentage.py`**      
   Performs subsampling to create **20 mutually exclusive subsets**, each representing 5% of the original dataset.

   * Each sample preserves the data distribution on supertypes.
   * All 20 samples collectively reconstruct the complete dataset.

3. **`step0.3_generate_fractionMatrix_for_5percentage.py`**       
   Generates a fraction matrix for each of the 20 subsampled datasets.

4. **`step0.4_compute_gene_pair_coexpression.py`**  **`step0.4_compute_gene_pair_coexpression.sh`**      
   Calculates cell-type-specific gene pair co-expression matrics using the full set or 4 subsamples of 25% cells (stabilized co-expression).

### **Step 1A: Feature Selection via Mutual Information (MI)**

5. **`step1.1_csv_fraction_supertype_train_MIs_6e10.R`**      
   Calculates MI-ranked genes for **6e10** protein levels across the 20 subsampled datasets.

6. **`step1.1_csv_fraction_supertype_train_MIs_AT8.R`**      
   Calculates MI-ranked genes for **AT8** protein levels across the 20 subsampled datasets.

7. **`step1.2_find_MI_features.Rmd`**     
   Integrates MI results and selects the final set of MI-based features used for model training.

### **Step 1B: Feature Selection for Gene Pairs**

8. **`step1B.1_prepare_donor_splits.py`**        
   Randomly split the donors into discovery and replication halves for 10 rounds    
   
9. **`step1B.2_select_pairs_10rounds.py`**  **`step1B.2_select_pairs_10rounds.sh`**         
   Select gene pairs based on the spearman's correlation between each gene pair's stabilized co-expression and the pathological target.
    
10. **`step1B.3_select_final_consensus_pairs.ipynb`**      
   Get the final consensus pair feature set.

### **Step 2: Model Training**

11. **`step2_final_model_6e10.ipynb`**       
   Trains the predictive model for **6e10** protein data.

12. **`step2_final_model_AT8.ipynb`**       
   Trains the predictive model for **AT8** protein data.




## **Dependencies**

The workflow uses the following key environments:

* Python ≥ 3.12
* R ≥ 4.3

