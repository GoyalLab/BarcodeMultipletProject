#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH --job-name gatk
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 150G
#SBATCH --time 24:00:00
#SBATCH --output=/projects/p31666/zzhang/cluster_logs/%j.%N.txt
#SBATCH --error=/projects/p31666/zzhang/cluster_logs/%j.%N.err
#SBATCH --verbose

module purge
source /home/zzj4347/.bashrc

#gatk MarkDuplicatesSpark -I /projects/p31666/zzhang/doublet-bchmk/data/demuxlet/bams/possorted_genome_bam.bam \
# -O /projects/p31666/zzhang/doublet-bchmk/data/demuxlet/bams/possorted_genome_bam_marked_duplicates.bam
# -M /projects/p31666/zzhang/doublet-bchmk/data/demuxlet/bams/marked_dup_metrics.txt

gatk --java-options "-Xmx64g" HaplotypeCaller  \
  -R /projects/p31666/zzhang/resource/genomes/hg19/genome.fa \
  -I /projects/p31666/zzhang/doublet-bchmk/data/demuxlet/bams/possorted_genome_bam.bam \
  -O /projects/p31666/zzhang/doublet-bchmk/data/demuxlet/vcf/output.vcf.gz \
  -debug
