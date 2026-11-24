#!/usr/bin/env bash
# Usage: run_gffcompare.sh GTF_LIST
#
# Compares each StringTie assembly with the TAIR12 coding transcripts. This
# is a sanity check on the assemblies, not an input to TRaCE: it shows how
# much of the annotation each tissue recovers, and how many assembled
# transcripts sit antisense to an annotated gene (class code "x"), which
# would be common if the strand flag were wrong.
#
# gffcompare writes its .tmap and .refmap files next to the query GTF, so
# each GTF is linked into qc/gffcmp first to keep stringtie/ free of extras.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
cd "$PROJ/qc/gffcmp"

list=$PROJ/$1
summary=summary.tsv
printf 'sample\tquery_transcripts\tintron_chain_sensitivity\tintron_chain_precision\tlocus_sensitivity\tlocus_precision\tmatch_complete\tantisense_x\tintergenic_u\n' > "$summary"

while read -r gtf; do
  sample=${gtf%.gtf}
  ln -sf "$PROJ/stringtie/$gtf" "$gtf"
  gffcompare -r "$PROJ/ref/TAIR12.trace.gff" -o "$sample" "$gtf" > "$sample.log" 2>&1
  rm -f "$gtf"
  awk -v s="$sample" -v tmap="$sample.$gtf.tmap" '
    /Query mRNAs/        { q = $5 }
    /Intron chain level/ { ics = $4; icp = $6 }
    /Locus level/        { ls = $3; lp = $5 }
    END {
      while ((getline line < tmap) > 0) { split(line, f, "\t"); code[f[3]]++ }
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\n", s, q, ics, icp, ls, lp, code["="], code["x"], code["u"]
    }' "$sample.stats" >> "$summary"
done < "$list"

column -t -s$'\t' "$summary"
