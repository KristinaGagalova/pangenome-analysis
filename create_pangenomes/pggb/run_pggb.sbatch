#!/bin/bash --login

# 1) scheduler resources
#SBATCH --account=pawsey1142
#SBATCH --partition=long
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32         # threads for PGGB
#SBATCH --mem=200G                 # or use an integer MiB value, e.g. 307200M
#SBATCH --array=0-15%4             # edit: 16 jobs total, at most 4 concurrent
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

# 2) software environment (Setonix uses modules)
module load singularity/4.1.0-slurm
# (You don't actually use Nextflow below, so no need to load it here)

# 3) inputs / container
PGGB=/scratch/pawsey1142/kgagalova/pangenomes/pggb/pggb_latest.sif
LIST=/scratch/pawsey1142/kgagalova/pangenomes/pggb/per_chromosome/list_batches.in  # <-- one input per line

# 4) Setonix best practices for multithreaded jobs
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=close
# Workaround recommended by Pawsey when sharing nodes:
export FI_CXI_DEFAULT_VNI=$(od -vAn -N4 -tu < /dev/urandom)

# 5) pick the input for this array task
INPUT=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$LIST")
if [[ -z "$INPUT" ]]; then
  echo "No input at index $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

# 6) derive a clean job-specific output dir
base=$(basename "$INPUT")
stem=${base%.gz}
stem=${stem%.fa}
stem=${stem%.fasta}
OUTDIR="output/${stem}.out"
mkdir -p "$OUTDIR"

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: INPUT=$INPUT  OUTDIR=$OUTDIR"
echo "Using ${SLURM_CPUS_PER_TASK} threads"

# 7) always launch the job step with srun on Setonix
srun -N 1 -n 1 -c "${SLURM_CPUS_PER_TASK}" \
  singularity run "$PGGB" pggb \
    -i "$INPUT" \
    -o "$OUTDIR" \
    -V Darmor:1000 \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --poa-threads "${SLURM_CPUS_PER_TASK}"
