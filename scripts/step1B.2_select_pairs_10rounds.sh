#!/bin/bash
#SBATCH --job-name=select_pair_features
#SBATCH --time=1-00:00
#SBATCH --array=1-9
#SBATCH --mem=200G      
#SBATCH --cpus-per-task=8      
#SBATCH --output=logs/%x_%A_%a.out

module load miniconda3
source activate vcc_env
conda activate vcc_env


CELLTYPES=('Microglia-PVM' 'Astrocyte' 'OPC' 'Oligodendrocyte' 'Sst' 'L2/3 IT' 'L5/6 NP' 'Pvalb' 'Sncg' 'Vip' 'Lamp5')
CELLTYPE="${CELLTYPES[$SLURM_ARRAY_TASK_ID]}"

REGION='MTG'
#REGION='A9'

python step1.4_select_pair_features_repeated_rounds.py --Subclass "$CELLTYPE" --Region $REGION


#cd results/pair_features
#for ct in */; do cp "${ct}/HVGs1000.txt" "../../docker/HVGs/${ct%/}_HVGs1000.txt"; done
