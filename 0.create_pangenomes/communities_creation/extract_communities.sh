#!/usr/bin/env bash
set -euo pipefail

FASTA="$1"     # e.g., bnapus16.fasta or bnapus16.fasta.gz
shift
COMM_FILES=("$@")

# Index FASTA if needed
if [[ ! -f "${FASTA}.fai" ]]; then
  echo "Indexing ${FASTA}"
  samtools faidx "$FASTA"
fi

for f in "${COMM_FILES[@]}"; do
  [[ -f "$f" ]] || { echo "Skipping missing file $f"; continue; }

  list="${f%.txt}.names.list"
  awk '{print $1}' "$f" | sed '/^[[:space:]]*$/d' | sort -u > "$list"

  out="${f%.txt}.fa"
  echo "Extracting sequences for ${f} -> ${out}.gz"

  # Extract and compress in one step
  samtools faidx "$FASTA" -r "$list" | bgzip -c > "${out}.gz"

  # Optionally index the compressed output
  samtools faidx "${out}.gz"
done

echo "All community FASTA files extracted and bgzipped."
