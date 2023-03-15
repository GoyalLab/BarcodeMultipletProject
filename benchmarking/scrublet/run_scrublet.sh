#!/bin/bash
#SBATCH -A p31666
#SBATCH -p normal
#SBATCH --job-name bchmk_scrublet
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 64G
#SBATCH --time 1:00:00
#SBATCH --output=/projects/p31666/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/p31666/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

module purge all
source /home/zzj4347/.bashrc

cd /projects/p31666/zzhang/doublet-bchmk/repo/benchmarking/scrublet
conda activate doublet-bchmk
python scrublet_main.py
