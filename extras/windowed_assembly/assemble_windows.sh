#!/bin/bash
# assemble_windows.sh -- run SPAdes on every window written by
# bin_sam_windows.pl and pool the contigs into one FASTA (experimental).
#
# Usage:
#   assemble_windows.sh <windows_dir> <kmers> [threads] [min_contig_len]
#     windows_dir     output directory of bin_sam_windows.pl
#     kmers           SPAdes -k list, e.g. 55,87,121 (use the value from the
#                     Fast-Plast progress log for the same reads)
#     threads         per-SPAdes thread count [4]
#     min_contig_len  drop pooled contigs shorter than this [500]
#
# Result: <windows_dir>/windowed_contigs.fasta with headers
#   >w007_35000-45000_NODE_1_length_9876_cov_41.2
# which can be handed to Fast-Plast's filter/afin stage in place of
# 3_Spades_Assembly/spades_iter1/contigs.fasta for comparison.
#
# Windows with no reads are skipped. A window whose SPAdes run fails is
# reported and skipped rather than aborting the whole batch.
set -uo pipefail

WDIR="${1:?windows_dir}"
KMERS="${2:?kmers}"
THREADS="${3:-4}"
MINLEN="${4:-500}"

command -v spades.py >/dev/null || { echo "spades.py not on PATH" >&2; exit 1; }

POOL="$WDIR/windowed_contigs.fasta"
: > "$POOL"
ok=0; fail=0; empty=0

for d in "$WDIR"/w*/; do
    d="${d%/}"; w="$(basename "$d")"
    args=()
    [[ -s "$d/R1.fq" && -s "$d/R2.fq" ]] && args+=(-1 "$d/R1.fq" -2 "$d/R2.fq")
    [[ -s "$d/U.fq" ]] && args+=(-s "$d/U.fq")
    if [[ ${#args[@]} -eq 0 ]]; then empty=$((empty+1)); continue; fi

    if spades.py -o "$d/spades" "${args[@]}" --only-assembler -k "$KMERS" -t "$THREADS" > "$d/spades.log" 2>&1 \
       && [[ -s "$d/spades/contigs.fasta" ]]; then
        ok=$((ok+1))
        # prefix headers with the window id and apply the length floor
        awk -v w="$w" -v min="$MINLEN" '
            /^>/ { if (seq != "" && length(seq) >= min) { print hdr; print seq }
                   hdr = ">" w "_" substr($0, 2); seq = ""; next }
                 { seq = seq $0 }
            END  { if (seq != "" && length(seq) >= min) { print hdr; print seq } }
        ' "$d/spades/contigs.fasta" >> "$POOL"
    else
        fail=$((fail+1)); echo "WARNING: SPAdes failed for $w (see $d/spades.log)" >&2
    fi
done

echo "windows assembled: $ok   failed: $fail   empty: $empty"
echo "pooled contigs: $(grep -c '^>' "$POOL")  -> $POOL"
