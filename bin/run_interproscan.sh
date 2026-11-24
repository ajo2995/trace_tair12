#!/usr/bin/env bash
# Usage: run_interproscan.sh [CPUS]
#
# Annotates every TAIR12 protein with Pfam domains, computed locally.
#
#   -appl Pfam  Pfam only, as in the TRaCE paper.
#   -dp         do not ask EBI's precomputed-match service; for this release
#               the service address returns a web page, so every lookup
#               would fail and fall back to local search anyway.
#
# Two tables come out of the raw result:
#   pfam.all.tsv  every Pfam hit. This is what TRaCE is given.
#   pfam.ipr.tsv  only hits whose family is integrated into InterPro, which
#                 is what the 2020 runs effectively counted. Kept for the
#                 comparison in the summary, nothing else.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
cd "$PROJ"

cpus=${1:-16}
ipr=$IPRSCAN_HOME/interproscan.sh
mkdir -p "$TMPDIR/ipr"

# The tarball ships its own search binaries; make sure they run on this OS
# before committing an hour to the real job.
if [ ! -e domains/selftest.ok ]; then
  rm -f domains/selftest.tsv
  "$ipr" -i "$IPRSCAN_HOME/test_all_appl.fasta" -appl Pfam -dp -f tsv \
      -o domains/selftest.tsv -T "$TMPDIR/ipr" -cpu 4 > logs/interproscan.selftest.log 2>&1
  [ -s domains/selftest.tsv ]
  touch domains/selftest.ok
fi

if [ ! -e domains/pfam.raw.ok ]; then
  rm -f domains/pfam.raw.part
  "$ipr" -i domains/TAIR12.pep.fa -appl Pfam -dp -f tsv \
      -o domains/pfam.raw.part -T "$TMPDIR/ipr" -cpu "$cpus" > logs/interproscan.log 2>&1
  mv domains/pfam.raw.part domains/pfam.raw.tsv
  touch domains/pfam.raw.ok
fi

# TRaCE.pl does not look at the database column, so only Pfam rows may reach it.
awk -F'\t' '$4=="Pfam"' domains/pfam.raw.tsv > domains/pfam.all.tsv
awk -F'\t' '$4=="Pfam" && $12 ~ /^IPR[0-9]+$/' domains/pfam.raw.tsv > domains/pfam.ipr.tsv

echo "raw rows:              $(wc -l < domains/pfam.raw.tsv)"
echo "Pfam rows:             $(wc -l < domains/pfam.all.tsv)"
echo "  with an IPR entry:   $(wc -l < domains/pfam.ipr.tsv)"
echo "column counts:         $(awk -F'\t' '{n[NF]++} END{for(k in n) printf "%d rows x %d cols  ", n[k], k}' domains/pfam.raw.tsv)"
echo "column 12 values:      $(awk -F'\t' '{k=($12 ~ /^IPR[0-9]+$/) ? "IPR" : "[" $12 "]"; n[k]++} END{for(k in n) printf "%s=%d  ", k, n[k]}' domains/pfam.raw.tsv)"
echo "proteins with a hit:   $(cut -f1 domains/pfam.all.tsv | sort -u | wc -l) of $(grep -c '>' domains/TAIR12.pep.fa)"
echo "IDs not in the GFF:    $(cut -f1 domains/pfam.all.tsv | sort -u | comm -23 - <(grep '>' domains/TAIR12.pep.fa | sed 's/^>//' | sort) | wc -l)"
