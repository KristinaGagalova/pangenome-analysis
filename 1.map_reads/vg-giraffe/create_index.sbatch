#!/bin/bash

#SBATCH --array=0-19                # adjust based on how many .gfa files you have
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --account=pawsey1142
#SBATCH --output=logs/vg_%A_%a.out
#SBATCH --error=logs/vg_%A_%a.err

module load singularity/4.1.0-slurm

# Load all GFA files into an array
files=(*.smooth.final.gfa)

# Pick the one for this array task
f=${files[$SLURM_ARRAY_TASK_ID]}

# Extract community number from filename
num=$(echo "$f" | sed -E 's/.*community\.([0-9]+)\..*/\1/')

# Construct prefix
prefix="bnapus.community.${num}"

echo "Processing community $num: $f"
echo "Output prefix: $prefix"

singularity run vg:1.67.0--h9ee0642_0 \
    vg autoindex -w giraffe \
    --prefix "$prefix" \
    --gfa "$f"
