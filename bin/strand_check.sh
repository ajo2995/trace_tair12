#!/usr/bin/env bash
# Usage: strand_check.sh sample.bam
#
# Measures library strandedness from the alignments themselves. It counts
# read 1 of each pair over exons of Chr1 genes that have no neighbour on the
# opposite strand (see make_strand_beds.py). In a TruSeq Stranded (dUTP)
# library read 1 is antisense to the transcript, which is what StringTie's
# --rf flag assumes.
#
# Prints the two counts and the antisense fraction. Exits non-zero when the
# fraction is below 0.85 or there are too few reads to judge, so the caller
# does not go on to assemble with the wrong strand flag.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
bam=$1
plus=$PROJ/ref/strandcheck.plus.bed
minus=$PROJ/ref/strandcheck.minus.bed

# count BED REQUIRED_FLAGS EXCLUDED_FLAGS
# 0x40 read 1, 0x10 reverse strand, 0x4 unmapped, 0x100 secondary, 0x800 supplementary
count() { samtools view -c -M -L "$1" -f "$2" -F "$3" "$bam"; }

antisense=$(( $(count "$plus" 0x50 0x904) + $(count "$minus" 0x40 0x914) ))
sense=$(( $(count "$plus" 0x40 0x914) + $(count "$minus" 0x50 0x904) ))

awk -v a="$antisense" -v s="$sense" 'BEGIN {
  total = a + s
  fraction = total ? a / total : 0
  printf "read1_antisense_RF\t%d\nread1_sense_FR\t%d\nRF_fraction\t%.4f\n", a, s, fraction
  exit !(total >= 100000 && fraction >= 0.85)
}'
