#!/usr/bin/env python3

# Convert probes from - https://snpdb.appliedbioinformatics.com.au/download to fastq files, where each variant has a separate entry
# The Phred score is maximum for each base

import argparse
import csv
import gzip
import os
import re
import sys
from glob import glob

SNP_RE = re.compile(r"\[([ACGT]+)/([ACGT]+)\]")

def phred_char(phred: int) -> str:
    """Convert numeric Phred to FASTQ ASCII (Phred+33)."""
    if phred < 0 or phred > 93:
        raise ValueError("Phred must be in [0, 93] for standard FASTQ ASCII range.")
    return chr(phred + 33)

def expand_snp_sequence(seq: str):
    """
    Expand a sequence containing exactly one SNP like ...[A/C]...
    Returns (allele1_seq, allele2_seq, allele1, allele2).
    """
    m = SNP_RE.search(seq)
    if not m:
        return None
    a1, a2 = m.group(1), m.group(2)
    # allow multi-base alleles too (e.g., [AT/A])
    left = seq[:m.start()]
    right = seq[m.end():]
    return left + a1 + right, left + a2 + right, a1, a2

def infer_chrom_from_filename(path: str) -> str:
    """
    Example: Brassica_napus_Darmor_v81_A01.csv -> A01
    Uses the last '_' chunk before .csv
    """
    base = os.path.basename(path)
    if not base.lower().endswith(".csv"):
        base = os.path.splitext(base)[0]
    stem = os.path.splitext(base)[0]
    parts = stem.split("_")
    return parts[-1] if parts else "UNK"

def csv_to_fastq_gz(csv_path: str, out_dir: str, phred: int = 41) -> str:
    chrom = infer_chrom_from_filename(csv_path)
    out_path = os.path.join(out_dir, f"snpdb_{chrom}.fastq.gz")

    qch = phred_char(phred)

    n_written = 0
    n_skipped_no_snp = 0
    n_skipped_short_row = 0

    with open(csv_path, "r", newline="") as f_in, gzip.open(out_path, "wt", newline="\n") as f_out:
        reader = csv.reader(f_in, delimiter=",")
        for row_idx, row in enumerate(reader, start=1):
            # Need at least 10 columns (1-indexed col10 => row[9])
            if len(row) < 10:
                n_skipped_short_row += 1
                continue

            seq = row[9].strip()  # 10th column
            exp = expand_snp_sequence(seq)
            if exp is None:
                n_skipped_no_snp += 1
                continue

            seq1, seq2, allele1, allele2 = exp

            # Build a header with useful metadata; keep it safe for FASTQ (no spaces preferred)
            # Columns based on your example:
            # 1:id, 2:species, 3:array, 4:chr1, 5:chr2, 6:[A/C], 7:strand, 8:pos1, 9:pos2, 10:seq
            snp_id   = row[0].strip()
            species  = row[1].strip() if len(row) > 1 else ""
            array    = row[2].strip() if len(row) > 2 else ""
            chr_a    = row[3].strip() if len(row) > 3 else ""
            chr_b    = row[4].strip() if len(row) > 4 else ""
            snp_tag  = row[5].strip() if len(row) > 5 else ""
            strand   = row[6].strip() if len(row) > 6 else ""
            pos1     = row[7].strip() if len(row) > 7 else ""
            pos2     = row[8].strip() if len(row) > 8 else ""

            src_file = os.path.basename(csv_path)

            # Write two FASTQ records (one per allele)
            for allele, s in [(allele1, seq1), (allele2, seq2)]:
                qual = qch * len(s)
                header = (
                    f"@{snp_id}|src={src_file}|chrom={chrom}|species={species}|array={array}"
                    f"|chrA={chr_a}|chrB={chr_b}|snp={snp_tag}|strand={strand}|pos1={pos1}|pos2={pos2}"
                    f"|allele={allele}"
                )
                f_out.write(header + "\n")
                f_out.write(s + "\n")
                f_out.write("+\n")
                f_out.write(qual + "\n")
                n_written += 1

    sys.stderr.write(
        f"[{os.path.basename(csv_path)}] wrote {n_written} FASTQ records to {out_path} "
        f"(skipped: no_snp={n_skipped_no_snp}, short_row={n_skipped_short_row})\n"
    )
    return out_path

def main():
    ap = argparse.ArgumentParser(
        description="Convert SNP CSVs (10th column sequence with [A/C]) into gzipped FASTQ (two alleles per row)."
    )
    ap.add_argument("-i", "--input", required=True,
                    help="Input CSV file OR a directory containing CSVs.")
    ap.add_argument("-o", "--outdir", default=".",
                    help="Output directory for snpdb_*.fastq.gz files (default: current dir).")
    ap.add_argument("--glob", default="*.csv",
                    help="If --input is a directory, glob pattern to match files (default: *.csv).")
    ap.add_argument("--phred", type=int, default=41,
                    help="Max Phred score to use for every base (default: 41 => 'J').")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    if os.path.isdir(args.input):
        paths = sorted(glob(os.path.join(args.input, args.glob)))
        if not paths:
            raise SystemExit(f"No files matched {args.glob} in {args.input}")
        for p in paths:
            csv_to_fastq_gz(p, args.outdir, phred=args.phred)
    else:
        csv_to_fastq_gz(args.input, args.outdir, phred=args.phred)

if __name__ == "__main__":
    main()
