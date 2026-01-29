#!/usr/bin/env python3
"""
Batch on/off-target analysis from BAMs.

Adds:
- total_reads (mapped + unmapped)
- unmapped_unique_reads

Expected chromosome is extracted from BAM filename (e.g. bnapus_C09.sorted.bam → C09)
Observed chromosome extracted from RNAME (e.g. Darmor#0#C05 → C05)

All alignments treated equally:
- no MAPQ
- no primary/secondary distinction
- ignore only FLAG 0x4 when counting mapped alignments
"""

import argparse
import glob
import os
import re
import subprocess
from collections import defaultdict
from typing import Dict, Optional, Set, Tuple


def samtools_view_stream(bam_path: str) -> subprocess.Popen:
    return subprocess.Popen(
        ["samtools", "view", "-h", bam_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )


def expected_from_bam_name(bam_filename: str, pattern: str) -> str:
    m = re.search(pattern, bam_filename)
    return m.group(1) if m else ""


def observed_from_rname(rname: str) -> str:
    if "#" in rname:
        return rname.split("#")[-1]
    if "|" in rname:
        return rname.split("|")[-1]
    return rname


def compute_stats_from_bam(
    bam_path: str,
    expected_chrom: str,
    write_sam_path: Optional[str] = None,
) -> Tuple[int, int, int, int, int, int, int, int, Dict[str, int]]:
    """
    Returns:
      total_reads,
      mapped_unique_reads,
      unmapped_unique_reads,
      any_on,
      any_off,
      on_only,
      off_only,
      both,
      total_mapped_alignments,
      off_contig_counts
    """
    all_reads: Set[str] = set()
    mapped_reads: Set[str] = set()
    on: Set[str] = set()
    off: Set[str] = set()
    off_contig: Dict[str, int] = defaultdict(int)
    total_align = 0

    sam_fh = open(write_sam_path, "w") if write_sam_path else None
    p = samtools_view_stream(bam_path)

    try:
        assert p.stdout is not None
        for line in p.stdout:
            if sam_fh:
                sam_fh.write(line)

            if not line or line.startswith("@"):
                continue

            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue

            qname = f[0]
            try:
                flag = int(f[1])
            except ValueError:
                continue

            rname = f[2]
            all_reads.add(qname)

            # unmapped
            if flag & 0x4:
                continue

            obs = observed_from_rname(rname)
            mapped_reads.add(qname)
            total_align += 1

            if obs == expected_chrom:
                on.add(qname)
            else:
                off.add(qname)
                off_contig[rname] += 1

        rc = p.wait()
        if rc != 0:
            raise RuntimeError(p.stderr.read())

    finally:
        if sam_fh:
            sam_fh.close()

    total_reads = len(all_reads)
    mapped_unique_reads = len(mapped_reads)
    unmapped_unique_reads = total_reads - mapped_unique_reads

    any_on = len(on)
    any_off = len(off)

    both = sum((q in on) and (q in off) for q in mapped_reads)
    on_only = sum((q in on) and (q not in off) for q in mapped_reads)
    off_only = sum((q in off) and (q not in on) for q in mapped_reads)

    return (
        total_reads,
        mapped_unique_reads,
        unmapped_unique_reads,
        any_on,
        any_off,
        on_only,
        off_only,
        both,
        total_align,
        off_contig,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam-dir", required=True)
    ap.add_argument("--pattern", default="*.sorted.bam")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--expected-regex", default=r"_(A\d{2}|C\d{2})")
    ap.add_argument("--write-sam", action="store_true")
    args = ap.parse_args()

    bams = sorted(glob.glob(os.path.join(args.bam_dir, args.pattern)))
    os.makedirs(args.outdir, exist_ok=True)
    per_bam = os.path.join(args.outdir, "per_bam")
    os.makedirs(per_bam, exist_ok=True)

    sam_dir = os.path.join(args.outdir, "sam") if args.write_sam else None
    if sam_dir:
        os.makedirs(sam_dir, exist_ok=True)

    summary = os.path.join(args.outdir, "mapping_summary.tsv")
    with open(summary, "w") as summ:
        summ.write(
            "bam\texpected_chrom\ttotal_reads\tmapped_unique_reads\tunmapped_unique_reads\t"
            "reads_with_any_ON\treads_with_any_OFF\ton_only\toff_only\tboth_on_and_off\t"
            "total_mapped_alignments\tstatus\n"
        )

        for bam in bams:
            base = os.path.basename(bam)
            expected = expected_from_bam_name(base, args.expected_regex)

            if not expected:
                summ.write(f"{base}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tSKIPPED\n")
                continue

            outprefix = os.path.join(per_bam, base.replace(".bam", ""))
            sam_path = os.path.join(sam_dir, base.replace(".bam", ".sam")) if sam_dir else None

            try:
                (
                    total_reads,
                    mapped_unique,
                    unmapped_unique,
                    any_on,
                    any_off,
                    on_only,
                    off_only,
                    both,
                    total_align,
                    off_contig,
                ) = compute_stats_from_bam(bam, expected, sam_path)

                with open(outprefix + ".counts.tsv", "w") as fh:
                    fh.write(
                        f"{total_reads}\t{mapped_unique}\t{unmapped_unique}\t"
                        f"{any_on}\t{any_off}\t{on_only}\t{off_only}\t{both}\t{total_align}\n"
                    )

                with open(outprefix + ".offtarget_contigs.tsv", "w") as fh:
                    for c, n in off_contig.items():
                        fh.write(f"{n}\t{c}\n")

                summ.write(
                    f"{base}\t{expected}\t{total_reads}\t{mapped_unique}\t{unmapped_unique}\t"
                    f"{any_on}\t{any_off}\t{on_only}\t{off_only}\t{both}\t{total_align}\tOK\n"
                )

            except Exception as e:
                summ.write(f"{base}\t{expected}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tFAILED:{e}\n")

    print(f"Done. Summary written to {summary}")


if __name__ == "__main__":
    main()
