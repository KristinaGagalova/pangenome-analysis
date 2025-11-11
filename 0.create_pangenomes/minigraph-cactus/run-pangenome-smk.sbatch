#!/bin/bash

#SBATCH --job-name=cactus_pangenome
#SBATCH --output=cactus_pangenome.out
#SBATCH --error=cactus_pangenome.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=4-00:00:00
#SBATCH --partition=long
#SBATCH --account=pawsey1142

# Move to your cactus directory
module load singularity/4.1.0-slurm

# Activate your virtual environment
source /scratch/pawsey1142/kgagalova/pangenomes/src/miniforge/etc/profile.d/conda.sh
eval "$(/scratch/pawsey1142/kgagalova/pangenomes/src/miniforge/bin/mamba shell hook --shell bash)"
mamba activate cactus-env

export PATH=/scratch/pawsey1142/kgagalova/pangenomes/src/miniforge/bin/:$PATH

snakemake -j 24 -e slurm \
	-s /scratch/pawsey1142/kgagalova/pangenomes/minigraph-cactus/cactus-snakemake/cactus_minigraph.smk \
	--configfile bnapus-minigraph-cfg.yaml \
	--default-resources slurm_account=pawsey1142 \
	--rerun-incomplete
