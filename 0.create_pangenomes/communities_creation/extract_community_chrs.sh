#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./make_community_segments.sh file_list.txt
#
# file_list.txt: one path to a .fa.gz per line, e.g.
#   /path/to/bnapus16....community.0.fa.gz
#   /path/to/bnapus16....community.1.fa.gz
#   ...

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 file_list.txt" >&2
    exit 1
fi

list="$1"
out="community_segments.tsv"

# Output header
echo -e "community_id\tsequence_id\tsample\thap\tchromosome\tbp_length" > "$out"

# Read each filename from the list
while IFS= read -r f; do
    # Skip empty lines
    [[ -z "$f" ]] && continue

    # Warn and skip if file not found
    if [[ ! -f "$f" ]]; then
        echo "Warning: file not found: $f" >&2
        continue
    fi

    # Extract community number from the filename: ...community.<NUM>.fa.gz
    comm=$(echo "$f" | sed -E 's/.*community\.([0-9]+)\.fa\.gz/\1/')

    gunzip -c "$f" | awk -v C="$comm" '
        BEGIN {
            OFS    = "\t"
            seq_len = 0
            seq_id  = ""
            sample  = "NA"
            hap     = "NA"
            chr     = "NA"
        }

        # Header line
        /^>/ {
            # Print previous sequence if present
            if (seq_len > 0) {
                clean_id = seq_id
                gsub(/^>/, "", clean_id)
                print C, clean_id, sample, hap, chr, seq_len
            }

            seq_len = 0
            seq_id  = $0

            # Remove leading ">"
            header = substr($0, 2)

            # Expected: Sample#Hap#Chromosome ...
            n = split(header, a, "#")
            if (n >= 3) {
                sample = a[1]
                hap    = a[2]
                chr    = a[3]
            } else {
                sample = "NA"
                hap    = "NA"
                chr    = "NA"
            }

            next
        }

        # Sequence lines
        {
            gsub(/[ \t\r]/, "", $0)
            if (length($0) > 0) {
                seq_len += length($0)
            }
        }

        END {
            # Print last sequence in file
            if (seq_len > 0) {
                clean_id = seq_id
                gsub(/^>/, "", clean_id)
                print C, clean_id, sample, hap, chr, seq_len
            }
        }
    ' >> "$out"

done < "$list"
