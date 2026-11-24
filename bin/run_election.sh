#!/usr/bin/env bash
# Usage: run_election.sh GTF_LIST OUT_DIR [DOMAIN_TSV]
#
# Holds the TRaCE election on the StringTie assemblies named in GTF_LIST, then
# runs the stability checks, the exact-arithmetic recount and the summary.
#
#   GTF_LIST    file of GTF names (one per line, sorted), relative to stringtie/
#   OUT_DIR     where results go (results/ for the real run)
#   DOMAIN_TSV  Pfam table for the main run (default domains/pfam.all.tsv:
#               every Pfam hit counts, as the user asked)
#
# TRaCE.pl is run from inside stringtie/ on an explicit list of files. It
# derives each sample's name from the path it is given and treats every
# *.gtf argument as a voter, so a glob or a stray file would change the
# electorate without any warning.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
cd "$PROJ"

list=$1
out=$2
domains=${3:-domains/pfam.all.tsv}
mkdir -p "$out/sensitivity"
tag=$(basename "$out")

sort -c "$list"
samples=$(wc -l < "$list")
while read -r gtf; do
  [ -e "stringtie/${gtf%.gtf}.ok" ] || { echo "assembly not finished: $gtf" >&2; exit 1; }
done < "$list"

# run_trace HASH_SEED DOMAIN_TSV OUTPUT LOG
# A fixed hash seed makes the order in which Perl walks its hashes, and so
# the order in which votes are added up, the same on every run.
run_trace() {
  ( cd stringtie
    PERL_HASH_SEED=$1 PERL_PERTURB_KEYS=0 "$SYS_PERL" "$TRACE_PL" ../ref/TAIR12.trace.gff "../$2" \
        "$MIN_TPM" "$MIN_EVI_COV" "$MIN_REF_COV" "$MAX_AED" $(cat "../$list") > "../$3" 2> "../$4" )
}

echo "election: $samples samples, domains from $domains"
run_trace 0 "$domains" "$out/TAIR12_TRaCE_output.txt" "logs/trace.$tag.log" &
pids=($!)
for seed in 1 2 3; do
  run_trace "$seed" "$domains" "$out/sensitivity/seed$seed.txt" "logs/trace.$tag.seed$seed.log" &
  pids+=($!)
done
if [ -s domains/pfam.ipr.tsv ]; then
  run_trace 0 domains/pfam.ipr.tsv "$out/sensitivity/integrated_only.txt" "logs/trace.$tag.integrated_only.log" &
  pids+=($!)
fi
for p in "${pids[@]}"; do wait "$p"; done

# ---- checks on the main run -------------------------------------------------
main=$out/TAIR12_TRaCE_output.txt
log=logs/trace.$tag.log
echo "== TRaCE.pl messages =="
echo "samples loaded: $(grep -c -x -F -f <(sed 's/\.gtf$//' "$list") "$log") of $samples"
echo "warnings: $(grep -c -i -E 'uninitialized|no exons in|isn.t numeric|division' "$log" || true)"
grep -i -E 'uninitialized|no exons in|isn.t numeric|division' "$log" | sort | uniq -c | sort -k1,1nr | head -5 || true

echo "== output shape =="
awk -F'\t' -v want=$(( 6 + samples )) '
  NF != want { bad++ }
  { genes[$1]++; if ($3 > top[$1]) top[$1] = $3 }
  END {
    for (g in genes) if (genes[g] != top[g]) gaps++
    printf "rows %d  genes %d  rows with a column count other than %d: %d  genes whose ranks are not 1..n: %d\n",
           NR, length(genes), want, bad + 0, gaps + 0
  }' "$main"

echo "== exon and CDS lengths against the source GFF3 =="
zcat ref/TAIR12_1Feb26.gff3.gz \
  | awk -F'\t' '!/^#/ && ($3=="exon" || $3=="CDS") { match($9, /Parent=[^;]+/); p = substr($9, RSTART+7, RLENGTH-7); s[p "\t" $3] += $5-$4+1 }
                END { for (k in s) print k "\t" s[k] }' \
  | awk -F'\t' 'NR==FNR { src[$1 "\t" $2] = $3; next }
                { if (src[$2 "\texon"] != $4) e++; if (src[$2 "\tCDS"] != $5) c++ }
                END { printf "transcripts whose exon bp differ: %d   whose CDS bp differ: %d\n", e + 0, c + 0 }' - "$main"

echo "== domain column against the Pfam table =="
echo "transcripts with domain aa > 0: $(awk -F'\t' '$6 > 0' "$main" | wc -l)   distinct proteins in $domains: $(cut -f1 "$domains" | sort -u | wc -l)"

# ---- stability ---------------------------------------------------------------
echo "== genes whose rank-1 transcript changes =="
for other in "$out"/sensitivity/*.txt; do
  awk -F'\t' 'NR==FNR { if ($3 == 1) top[$1] = $2; next } $3 == 1 && top[$1] != $2 { n++ } END { print n + 0 }' "$main" "$other" \
    | sed "s|^|$(basename "$other" .txt)\t|"
done | tee "$out/sensitivity/rank1_changes.tsv"

# ---- recount and summary -----------------------------------------------------
echo "== exact-arithmetic recount =="
python3 bin/recount_exact.py "$main" "$out/exact_recount" "$MAX_AED" | tee "$out/exact_recount.summary.txt"

echo "== summary tables =="
python3 bin/summarize_trace.py --trace "$main" --samples "$list" --gff ref/TAIR12.trace.gff \
    --stringtie-dir stringtie --out-dir "$out" --recount "$out/exact_recount.tsv" --proj "$PROJ" \
    --max-aed "$MAX_AED" --min-tpm "$MIN_TPM"
