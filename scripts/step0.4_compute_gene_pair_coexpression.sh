#!/bin/bash
#SBATCH --job-name=compute_pair_features
#SBATCH --time=1-00:00    
#SBATCH --array=0-8
#SBATCH --mem=300G      
#SBATCH --cpus-per-task=8      
#SBATCH --output=logs/%x_%A_%a.out

module load miniconda3
source activate vcc_env
conda activate vcc_env


CELLTYPES=('OPC' 'Oligodendrocyte' 'Sst' 'L2/3 IT' 'L5/6 NP' 'Pvalb' 'Sncg' 'Vip' 'Lamp5')
CELLTYPE="${CELLTYPES[$SLURM_ARRAY_TASK_ID]}"

python step0.4_compute_gene_pair_coexpression.py --Subclass "$CELLTYPE"
