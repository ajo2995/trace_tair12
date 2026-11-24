# TRaCE ranking of TAIR12 transcripts

This repository holds a TRaCE run on the TAIR12 protein-coding gene models of *Arabidopsis thaliana* Col-0 (Col-CC v2 assembly): the pipeline, the ranking of every transcript, and a comparison with TAIR's representative gene models. The run was done on 2026-11-20.

TRaCE ([Olson & Ware 2022, Bioinformatics 38:261](https://pubmed.ncbi.nlm.nih.gov/34297055/); code at https://github.com/warelab/trace) ranks the transcripts of each gene by ranked-choice voting. Each RNA-seq sample votes for the annotated transcript closest to the dominant transcript assembled from that sample. Three more voters rank the transcripts by Pfam domain coverage, CDS length and transcript length. Rank 1 is the canonical transcript.

## Main results

- All 26,867 protein-coding genes (36,572 transcripts) are ranked. The ranking is a real choice for the 5,699 genes that have more than one transcript.
- At least one of the 30 tissues voted in 5,517 of those 5,699 genes (96.8%).
- The canonical transcript is not the one with the longest CDS in 1,601 of them (28.1%).
- The canonical transcript differs from TAIR's representative gene model, which is the longest spliced transcript, in 2,942 genes. That is 51.6% of multi-transcript genes and 11.0% of all coding genes.
- In those 2,942 genes:
  - the representative is longer by a median of 89 bp;
  - the two encode the same protein in 1,058 (36.0%);
  - the representative reads through an intron that the canonical splices out in 1,147 (39.0%);
  - of the 2,907 with RNA-seq evidence, more tissues match the canonical than the representative in 2,857 (98.3%).
- 54 genes are effectively ties that floating-point rounding in `TRaCE.pl` decides, and 182 multi-transcript genes had no RNA-seq vote.

[`results/summary.md`](results/summary.md) has the full numbers, the checks that were run and the caveats.

## What is here

| Path | Contents |
|---|---|
| `results/summary.md` | Headline numbers, checks and caveats. Start here. |
| `results/` | The ranking, the comparisons and the QC tables (listed below). |
| `bin/` | The pipeline scripts (listed below). |
| `config.sh` | Settings shared by every script: paths, the TRaCE commit, the thresholds. |
| `manifest/` | The 30 RNA-seq samples: run accessions, read counts, tissue labels and download URLs. |
| `logs/conda_env.txt` | The exact packages installed in the conda environment. |

### Results

| File | Contents |
|---|---|
| `results/TAIR12_TRaCE_output.txt` | TRaCE's own output, unmodified, no header. This is the file `set_canonical_transcripts_from_TRaCE.pl` reads. |
| `results/TAIR12_TRaCE_output.columns.txt` | The column names for that file, one per line. |
| `results/TAIR12_canonical_transcripts.tsv` | One row per gene: the canonical transcript and how it was supported. |
| `results/representative_vs_trace.tsv` | For each multi-transcript gene: TAIR's representative gene model, where TRaCE ranked it, and how it differs from TRaCE's canonical (length, protein, Pfam coverage, intron retention, tissue support). `representative_vs_trace.summary.txt` has the totals. |
| `results/ensembl_canonical_vs_trace.tsv` | The same comparison against the canonical transcripts of an Ensembl core database built from the TAIR12 files, as that database was first loaded. Those follow longest CDS first and differ from TAIR's representatives in 870 genes. |
| `results/tair_representative_rules.tsv`, `results/ensembl_canonical_rules.tsv` | For each multi-transcript gene: TAIR's representative (or the database's canonical) next to the pick of each candidate selection rule. The `.canonical.tsv` files list the representative (or canonical) of every coding gene. |
| `results/exact_recount.tsv` | Genes whose ranking changes when the votes are added with exact arithmetic. |
| `results/qc_*.tsv` | Per-sample alignment and evidence numbers, the histogram of voting samples, the genes no sample voted on, and each assembly compared with TAIR12. |
| `results/sensitivity/` | Full TRaCE output of four reruns: three with other Perl hash seeds, one counting only InterPro-integrated Pfam hits. `rank1_changes.tsv` gives the number of genes whose canonical changes in each. |

### Scripts

| Script | What it does |
|---|---|
| `bin/get_interproscan.sh` | Downloads InterProScan in ranged streams, checks its md5 and unpacks it. |
| `bin/make_trace_gff.py` | Rewrites the TAIR12 GFF3 into the form TRaCE's parser needs. No coordinate is changed. |
| `bin/make_strand_beds.py` | Writes the exon sets that `strand_check.sh` counts reads over. |
| `bin/build_manifest.py` | Writes the sample sheet and the lists in `manifest/`. |
| `bin/run_interproscan.sh` | Runs Pfam through InterProScan on the TAIR12 proteins. |
| `bin/sample_worker.sh` | Handles one sample: download, align, check the BAM, measure strandedness, assemble. Safe to rerun. |
| `bin/strand_check.sh` | Measures library strandedness from a BAM. |
| `bin/run_gffcompare.sh` | Compares each assembly with the TAIR12 gene models. |
| `bin/run_election.sh` | Runs `TRaCE.pl`, the four reruns, the exact recount and the summary tables. |
| `bin/recount_exact.py` | Recounts every election from TRaCE's output with exact fractions. |
| `bin/summarize_trace.py` | Writes the canonical table and the QC tables. |
| `bin/compare_representative.py` | Compares TRaCE's ranking with one representative transcript per gene. |
| `bin/compare_ensembl_canonical.py` | Tests which selection rule reproduces a database's canonical transcripts, or a supplied list such as TAIR's. |
| `bin/fetch_ensembl_canonical.pl` | Read-only dump of gene models and canonical flags from an Ensembl core database. |
| `bin/inspect_ensembl_organelles.pl`, `bin/check_ensembl_organelles.py` | Check that the organelle regions of an Ensembl core database carry the right names. |
| `bin/check_ensembl_trace_canonical.py` | Confirms that an Ensembl core database's canonical transcripts are TRaCE's rank 1. |

## What is not here

Inputs, third-party software and intermediate files are left out. The sizes are what the run left on disk.

| Left out | Size | How to get it |
|---|---|---|
| TAIR12 genome and GFF3 | 50 MB | Download from TAIR (see Inputs). |
| TAIR's list of representative gene models | 0.5 MB | Obtain from TAIR. |
| RNA-seq reads | 59 GB as `.sra` | `bin/sample_worker.sh` downloads each run and deletes it once the BAM passes its checks. |
| Alignments (`bam/`) | 57 GB | Made by `bin/sample_worker.sh`. |
| StringTie assemblies (`stringtie/`), 30 GTFs | 850 MB | Made by `bin/sample_worker.sh`. These are TRaCE's evidence input. |
| Pfam hits and proteins (`domains/`) | 44 MB | Made by `bin/run_interproscan.sh`. The Pfam table is TRaCE's domain input. |
| Reformatted GFF, alignment index, organelle sequences (`ref/`) | 530 MB | Made by the reference commands below. |
| gffcompare output and strand counts (`qc/`), logs | 1.1 GB | Made by the scripts. The numbers are collected in `results/qc_*.tsv`. |
| `TRaCE.pl` | 16 kB | https://github.com/warelab/trace at the commit in `config.sh`. |
| InterProScan 5.78-109.0 | 42 GB | `bin/get_interproscan.sh`. |
| conda environment (`env/`) | 540 MB | The `mamba create` command below; `logs/conda_env.txt` lists what it installed. |
| Tables dumped from the Ensembl core database (`ensembl/`) | 19 MB | `bin/fetch_ensembl_canonical.pl`, given access to the database. |

## Inputs

- **Gene models and genome.** `TAIR12_1Feb26.gff3.gz` (revision of 2026-08-28) and `Col-CC_v2_genome.fasta.gz`, fetched from TAIR on 2026-11-18 at `https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR12_genome_release/<file>`. The assembly is GenBank GCA_978657495.1. It has Chr1-Chr5 only.

  ```
  75570fce3df601202e4ed09c7518ffb8d943fd836cf22c8d91ad486b86e97298  Col-CC_v2_genome.fasta.gz
  242de453d7ea3a18468799ffe0dbef9977e3293d94a5c258279753b72386efd0  TAIR12_1Feb26.gff3.gz
  ```

- **RNA-seq.** The 30-tissue core (Samples 1-30) of the Mergner et al. 2020 atlas (Nature 579:409), ArrayExpress E-MTAB-7978 / ENA PRJEB32665. Col-0 wild type, polyA, TruSeq Stranded, paired 2x75. One library, and so one voter, per tissue. See `manifest/samples.core30.tsv`.
- **Domains.** Pfam 38.2 through InterProScan 5.78-109.0, computed locally. All Pfam hits count toward domain coverage, including families not integrated into InterPro.
- **TRaCE.** `TRaCE.pl` at commit `836479f6bf55ca823ad7246134bdcb1aa070b267` (git blob `b316baa5`), unmodified, with thresholds 0.5 / 0.5 / 0.5 / 0.5.
- **TAIR representative gene models.** `TAIR12_representative.txt`, one transcript ID per gene (40,206 lines), obtained from TAIR on 2026-11-23. It is used only for the comparison with TRaCE, not in the election. For all 5,699 multi-transcript coding genes the list follows one rule: the longest spliced transcript, with ties going to the lower isoform number.

## Software

| Tool | Version | Source |
|---|---|---|
| HISAT2 | 2.2.3 | bioconda |
| StringTie | 3.0.3 | bioconda |
| samtools | 1.24 | bioconda |
| gffread / gffcompare | 0.12.9 / 0.12.10 | bioconda |
| sra-tools | 3.4.1 | bioconda |
| Python | 3.12 | conda-forge |
| InterProScan | 5.78-109.0 | EBI FTP |
| Perl (for `TRaCE.pl`) | 5.34.0 | system |

## How it was run

These are the commands as they were run, from the top of this directory. A rerun from a fresh clone has not been tested. It took about two hours end to end on a 64-core machine, with four samples processed at a time.

```bash
source config.sh
mkdir -p tmp sw/trace ref/hisat2 domains sra fastq bam stringtie \
         qc/strand qc/gffcmp logs/hisat2 logs/stringtie logs/worker

# Software
$MAMBA create -y -q -p "$PROJ/env" --strict-channel-priority -c conda-forge -c bioconda \
    hisat2 'stringtie=3.0.3' gffread gffcompare samtools sra-tools python=3.12
curl -fsSL -o sw/trace/TRaCE.pl "https://raw.githubusercontent.com/warelab/trace/$TRACE_COMMIT/TRaCE.pl"
bin/get_interproscan.sh

# Reference: nuclear genome, plus the two organelle genomes for alignment only
cp "$SRC/Col-CC_v2_genome.fasta.gz" "$SRC/TAIR12_1Feb26.gff3.gz" "$SRC/SHA256SUMS" ref/
(cd ref && grep -E 'Col-CC_v2_genome.fasta.gz|TAIR12_1Feb26.gff3.gz' SHA256SUMS | sha256sum -c -)
zcat ref/Col-CC_v2_genome.fasta.gz | awk '/^>/{print $1; next} {print}' > ref/genome.nuc.fa
for pair in AP000423.1:ChrC BK010421.1:ChrM; do
  acc=${pair%%:*}; name=${pair##*:}
  curl -fsS --retry 5 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$acc&rettype=fasta&retmode=text" \
    | awk -v n="$name" 'NR==1{print ">"n; next} NF{print}'
  sleep 1
done > ref/organelles.fa
cat ref/genome.nuc.fa ref/organelles.fa > ref/genome.idx.fa
samtools faidx ref/genome.nuc.fa
samtools faidx ref/genome.idx.fa
hisat2-build -p 16 ref/genome.idx.fa ref/hisat2/TAIR12_CM

# Gene models and proteins
python3 bin/make_trace_gff.py ref/TAIR12_1Feb26.gff3.gz ref/TAIR12.trace.gff
gffread ref/TAIR12.trace.gff -g ref/genome.nuc.fa -y ref/pep.raw.fa
awk '/^>/{print $1; next} {print}' ref/pep.raw.fa > domains/TAIR12.pep.fa
python3 bin/make_strand_beds.py ref/TAIR12_1Feb26.gff3.gz ref/TAIR12.trace.gff ref/strandcheck
python3 bin/build_manifest.py manifest

# Domains, assemblies, election
bin/run_interproscan.sh 16
xargs -P 4 -L 1 bin/sample_worker.sh < manifest/worker_args.core30.txt   # one sample per line
bin/run_gffcompare.sh manifest/gtfs.core30.txt
bin/run_election.sh manifest/gtfs.core30.txt results
```

What each step is for:

- **`make_trace_gff.py`** reformats the GFF3 for TRaCE's parser, which needs `biotype=protein_coding` on gene lines and exons in ascending order. TAIR12 uses `locus_biotype=` and lists most minus-strand exons descending.
- **The alignment reference** adds the chloroplast (AP000423.1, as `ChrC`) and the mitochondrion (BK010421.1, as `ChrM`), so organelle reads have somewhere to go other than their nuclear copies.
- **`sample_worker.sh`** downloads one run from the AWS mirror of SRA, aligns it with `hisat2 --dta --max-intronlen 10000`, checks that the BAM holds every read, measures strandedness, and assembles with `stringtie --rf --conservative -x ChrC,ChrM`. Neither tool is given the annotation.
- **`run_election.sh`** runs `TRaCE.pl` on an explicit, sorted list of GTFs with a fixed Perl hash seed, then the reruns, the recount and the summary tables.

### Comparison with the representative gene models

`compare_ensembl_canonical.py` takes its gene structures from tables dumped from the Ensembl core database, so the dump comes first, also when the list being tested is TAIR's.

```bash
db=arabidopsis_thaliana_core_10_108_12
mkdir -p ensembl/$db
perl bin/fetch_ensembl_canonical.pl REGISTRY.pm $db ensembl/$db

python3 bin/compare_ensembl_canonical.py ensembl/$db ref/TAIR12.trace.gff \
    results/TAIR12_TRaCE_output.txt results/ensembl_canonical_rules
python3 bin/compare_ensembl_canonical.py ensembl/$db ref/TAIR12.trace.gff \
    results/TAIR12_TRaCE_output.txt results/tair_representative_rules TAIR12_representative.txt

python3 bin/compare_representative.py --trace results/TAIR12_TRaCE_output.txt \
    --gff ref/TAIR12.trace.gff --proteins domains/TAIR12.pep.fa --recount results/exact_recount.tsv \
    --representatives results/tair_representative_rules.canonical.tsv \
    --out results/representative_vs_trace.tsv > results/representative_vs_trace.summary.txt
python3 bin/compare_representative.py --trace results/TAIR12_TRaCE_output.txt \
    --gff ref/TAIR12.trace.gff --proteins domains/TAIR12.pep.fa --recount results/exact_recount.tsv \
    --representatives results/ensembl_canonical_rules.canonical.tsv \
    --out results/ensembl_canonical_vs_trace.tsv > results/ensembl_canonical_vs_trace.summary.txt
```

- `REGISTRY.pm` is an Ensembl registry file that defines `$def_host`, `$def_port`, `$def_user` and `$def_pass`. The Perl scripts read the connection settings from it, open a read-only session and issue only `SELECT` and `SHOW` statements. They need `DBI` and `DBD::mysql`.
- **The database comparison can no longer be rerun.** The database's canonical transcripts were replaced with TRaCE's on 2026-11-23, so a fresh dump no longer shows the canonicals it was loaded with. `results/ensembl_canonical_rules.canonical.tsv` is the record of those, and `results/tair_representative_rules.tsv` carries them in one of its columns.
- **TAIR's list is not needed to rebuild `representative_vs_trace.tsv`.** Without `--representatives`, `compare_representative.py` applies the rule itself (longest spliced transcript, ties to the lower isoform number) and writes an identical table.
- The remaining scripts, `inspect_ensembl_organelles.pl`, `check_ensembl_organelles.py` and `check_ensembl_trace_canonical.py`, were used to check that database before and after it was corrected. `results/summary.md` reports what they found.

## Choices made for this run

- **Electorate:** the 30-tissue core of the atlas, one library per tissue. The atlas's developmental-series libraries and long-read data were not used.
- **Stranded libraries only.** TRaCE matches evidence to genes on the same strand, and 3,725 TAIR12 coding genes overlap a coding gene on the opposite strand.
- **Domain coverage counts every Pfam hit,** including families not integrated into InterPro. Counting only integrated hits changes the canonical transcript of 9 genes.
- **StringTie flags:** `--rf` for the library type and `--conservative`, the assembly mode used for the human run in the paper. `-x ChrC,ChrM` keeps organelle transcripts out of the TPM denominator.
- **HISAT2:** introns capped at 10 kb and no read trimming, as in the paper. Five genes have a longer annotated intron.
- **A fixed Perl hash seed** (`PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0`), because the canonical transcript of six genes depends on hash order.
- **`TRaCE.pl` was not modified.** The behaviours found along the way, including the rounding effect, are measured and reported in `results/summary.md`, not corrected.

## Notes for a rerun

- `config.sh` finds the project directory from its own location. `SRC`, `MAMBA` and `SYS_PERL` default to the paths used for this run and can be set in the environment instead.
- `config.sh` points `TMPDIR` at `tmp/` and every tool is given an explicit temporary directory, because the root filesystem of the machine used was nearly full.
- The run left about 100 GB in the working directory, nearly all of it the alignments (57 GB) and InterProScan (42 GB). Neither is needed once the assemblies and the Pfam table exist.
