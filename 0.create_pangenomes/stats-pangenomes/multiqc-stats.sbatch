#!/bin/bash
#SBATCH --array=0-19                   # adjust to number of communities
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/odgi_%A_%a.out
#SBATCH --error=logs/odgi_%A_%a.err

module load singularity/4.1.0-slurm
mkdir -p odgi_stats

files=(all_communities/*.smooth.final.og)
f=${files[$SLURM_ARRAY_TASK_ID]}

base=$(basename "$f" .smooth.final.og)
yaml="odgi_stats/${base}.og.stats.yaml"

singularity run odgi:0.9.3--h077b44d_0 odgi stats -i "$f" -m -sgdl > "$yaml"
