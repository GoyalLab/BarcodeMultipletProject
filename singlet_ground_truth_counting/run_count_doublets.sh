#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH --job-name count_doublets
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 64G
#SBATCH --time 24:00:00
#SBATCH --output=/projects/p31666/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/p31666/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

module purge
source /home/zzj4347/.bashrc

cd /projects/p31666/zzhang/doublet-bchmk/repo/singlet_ground_truth_counting
conda activate doublet-bchmk
python count_doublets_main.py
