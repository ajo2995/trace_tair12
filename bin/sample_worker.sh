#!/usr/bin/env bash
# Usage: sample_worker.sh SAMPLE RUN READ_PAIRS
#
# Takes one RNA-seq run from download to a StringTie assembly:
#   download .sra (AWS mirror) -> FASTQ -> HISAT2 alignment -> sorted BAM
#   -> strand check -> StringTie GTF
#
# Each stage leaves a ".ok" marker and writes its output under a temporary
# name first, so the script can be rerun after an interruption and will pick
# up at the first unfinished stage. The downloaded reads are deleted once the
# BAM has passed its checks; the BAM is kept.
#
# Neither the aligner nor the assembler is shown the TAIR12 annotation: the
# assembly has to be independent evidence for TRaCE to compare against.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
cd "$PROJ"

sample=$1
run=$2
pairs=$3

say() { printf '[%s] %-26s %s\n' "$(date '+%F %T')" "$sample" "$*"; }

if [ -e "stringtie/$sample.ok" ]; then
  say "already finished"
  exit 0
fi

if [ ! -e "bam/$sample.ok" ]; then
  url=https://sra-pub-run-odp.s3.amazonaws.com/sra/$run/$run

  if [ ! -e "fastq/$run.ok" ]; then
    if [ ! -e "sra/$run.ok" ]; then
      want=$(curl -fsSI "$url" | awk 'tolower($1)=="content-length:"{print $2}' | tr -d '\r')
      have=$(stat -c %s "sra/$run.sra" 2>/dev/null || echo 0)
      if [ "$have" -ne "$want" ]; then
        say "downloading $(( want / 1000000 )) MB"
        curl -fsS --retry 8 --retry-delay 15 -C - -o "sra/$run.sra" "$url"
      fi
      [ "$(stat -c %s "sra/$run.sra")" -eq "$want" ]
      vdb-validate "sra/$run.sra" > "logs/worker/$sample.vdb-validate.log" 2>&1
      touch "sra/$run.ok"
    fi

    say "converting to FASTQ"
    rm -f fastq/"$run"*.fastq
    fasterq-dump "sra/$run.sra" --split-3 -e 8 -t "$TMPDIR" -O fastq \
        > "logs/worker/$sample.fasterq-dump.log" 2>&1
    touch "fastq/$run.ok"
  fi

  say "aligning"
  # --dta: report alignments tailored for transcript assemblers.
  # --max-intronlen: 99.99% of annotated Arabidopsis introns are under 8.1 kb,
  #   and the HISAT2 default of 500 kb invites spurious long-range splices in
  #   a genome this compact.
  # No --rna-strandness: HISAT2 then tags strand only on spliced reads, from
  #   the splice motif, and StringTie's --rf handles the rest.
  hisat2 -p 10 --dta --max-intronlen 10000 --temp-directory "$TMPDIR" \
      -x ref/hisat2/TAIR12_CM -1 "fastq/${run}_1.fastq" -2 "fastq/${run}_2.fastq" \
      --new-summary --summary-file "logs/hisat2/$sample.summary" 2> "logs/hisat2/$sample.log" \
    | samtools sort -@ 2 -m 2G -T "$TMPDIR/$sample.sort" -o "bam/$sample.part.bam" - \
        2> "logs/hisat2/$sample.sort.log"

  # A truncated download or a broken pipe shows up as missing reads. HISAT2
  # writes one primary record per read, aligned or not.
  samtools quickcheck "bam/$sample.part.bam"
  records=$(samtools view -@ 4 -c -F 0x900 "bam/$sample.part.bam")
  if [ "$records" -ne $(( 2 * pairs )) ]; then
    say "BAM holds $records reads but ENA lists $(( 2 * pairs ))"
    exit 1
  fi
  mv "bam/$sample.part.bam" "bam/$sample.bam"
  samtools index -@ 4 "bam/$sample.bam"
  samtools idxstats "bam/$sample.bam" > "qc/idxstats.$sample.txt"
  touch "bam/$sample.ok"
  rm -f "sra/$run.sra" "sra/$run.ok" fastq/"$run"*.fastq "fastq/$run.ok"
fi

say "checking strandedness"
if ! bin/strand_check.sh "bam/$sample.bam" > "qc/strand/$sample.txt"; then
  say "strand check FAILED: $(tr '\n\t' ' =' < "qc/strand/$sample.txt")"
  exit 2
fi

say "assembling"
# --rf: read 1 is antisense to the transcript (TruSeq Stranded), confirmed above.
# --conservative: the assembly mode of the published human run.
# -x: organelle transcripts are not assembled, so they stay out of the TPM
#   denominator and TRaCE's 0.5 TPM cutoff means the same in every tissue.
# No -G: the assembly must not be guided by the annotation being ranked.
stringtie "bam/$sample.bam" --rf --conservative -x ChrC,ChrM -p 8 -l "$sample" \
    -o "stringtie/$sample.part" 2> "logs/stringtie/$sample.log"

transcripts=$(awk -F'\t' '$3=="transcript"' "stringtie/$sample.part" | wc -l)
if [ "$transcripts" -lt 10000 ]; then
  say "only $transcripts transcripts assembled"
  exit 1
fi
mv "stringtie/$sample.part" "stringtie/$sample.gtf"
touch "stringtie/$sample.ok"
say "finished: $transcripts transcripts, RF fraction $(awk '$1=="RF_fraction"{print $2}' "qc/strand/$sample.txt")"
