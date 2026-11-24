# TRaCE on TAIR12: results and checks

Run date 2026-11-20. For inputs, software versions and how to rerun, see `../README.md`.

## Result

All 26,867 TAIR12 protein-coding genes are ranked (36,572 transcripts). The ranking is a real choice for the 5,699 genes that have more than one transcript (15,404 transcripts); the other 21,168 genes have a single transcript, which is rank 1 by default.

- `TAIR12_TRaCE_output.txt` is TRaCE's output as written, ready for `set_canonical_transcripts_from_TRaCE.pl`.
- `TAIR12_canonical_transcripts.tsv` has one row per gene with the support behind the choice.

The electorate was 30 tissues (one stranded library each) plus TRaCE's three length voters. Domain coverage counts every Pfam hit.

## Headline numbers

Shares are of the 5,699 multi-transcript genes.

| | Genes | Share |
|---|---|---|
| At least one tissue voted | 5,517 | 96.8% |
| All 30 tissues voted | 881 | 15.5% |
| Canonical has the lowest AED in most of the tissues that voted | 5,041 | 88.5% |
| Canonical is not the longest CDS | 1,601 | 28.1% |
| Canonical is not the longest transcript | 2,930 | 51.4% |
| Canonical is neither the longest CDS nor the longest transcript | 1,376 | 24.1% |
| Canonical is not the `.1` model | 3,037 | 53.3% |
| Canonical has less Pfam coverage than another isoform | 286 | 5.0% |

- In 1,607 of these genes every transcript has the same CDS length, so CDS length cannot separate them. Among the other 4,092, the canonical is not the longest CDS in 39.1%.
- For comparison, the paper reports TRaCE's pick as "not the longest" in 37% of maize and 27% of human multi-transcript genes.
- Voting tissues per gene: none for 182 genes, 1-9 for 782, 10-19 for 733, 20-29 for 3,121, all 30 for 881 (`qc_voting_samples.tsv`).

## TRaCE against the TAIR12 representative models

The representative model of each gene comes from TAIR's own list (`TAIR12_representative.txt`, copied from TAIR on 2026-11-23). For all 5,699 multi-transcript coding genes the list follows one rule exactly: the longest spliced transcript, with ties going to the lower isoform number.

TRaCE ranks the representative below another transcript in 2,942 genes: 51.6% of multi-transcript genes, 11.0% of all coding genes. The representative is rank 2 in 2,618 of them, rank 3 in 252 and rank 4 or lower in 72.

In those 2,942 genes, comparing TRaCE's canonical with the representative:

| Measure | Result |
|---|---|
| Spliced length | representative longer by a median of 89 bp (quartiles 30 to 178), 5.1% of its length; equal in 12 |
| Protein | identical in 1,058 (36.0%); in the other 1,884 the canonical protein is shorter in 1,297 (median 12 aa), longer in 572 (median 41 aa), equal in length in 15 |
| Pfam coverage | canonical has more in 608, the same in 2,173, less in 161 |
| Structure | a median of 93.6% of exonic bases shared; the representative reads through an intron the canonical splices out in 1,147 (39.0%), the reverse in 44 |
| RNA-seq evidence | of the 2,907 with a voting tissue, the canonical is closer to the assembly in more tissues in 2,857 (98.3%) and in every voting tissue in 1,579; median mean AED 0.0003 against 0.033 |

- 35 of the 2,942 had no voting tissue, so domain coverage and CDS length decided them. 38 are among the rounding-decided genes below.
- Per-gene values are in `representative_vs_trace.tsv`.

### The Ensembl core database's canonical transcripts are a different set

The Ensembl core database `arabidopsis_thaliana_core_10_108_12`, built from the TAIR12 files, holds the same 26,867 nuclear coding genes and 36,572 transcripts as the TAIR12 GFF, with identical exon coordinates and CDS lengths. Its canonical flags are not TAIR's representative models. For all 5,699 multi-transcript genes they follow Ensembl's own rule: longest CDS, then longest spliced transcript, then the transcript ID compared as text. They differ from TAIR's representative in 870 genes (15.3%).

TRaCE ranks the database's canonical below another transcript in 2,671 genes (46.9%). Across the 5,699 multi-transcript genes:

| | Genes | Share |
|---|---|---|
| TAIR, database and TRaCE all name the same transcript | 2,533 | 44.4% |
| TAIR and the database agree, TRaCE differs | 2,296 | 40.3% |
| TRaCE agrees with the database, TAIR differs | 495 | 8.7% |
| TRaCE agrees with TAIR, the database differs | 224 | 3.9% |
| All three differ | 151 | 2.6% |

- `ensembl_canonical_vs_trace.tsv` is the per-gene comparison against the database's canonical.
- `tair_representative_rules.tsv` and `ensembl_canonical_rules.tsv` show, gene by gene, what each candidate rule picks against TAIR's list and against the database.
- The comparison above describes the database as it was loaded. On 2026-11-23 its canonical transcripts were replaced with TRaCE's: `gene.canonical_transcript_id` now equals TRaCE's rank 1 for all 26,867 ranked genes, a change in exactly the 2,671 genes where the two differed. The earlier canonicals are kept in `ensembl_canonical_rules.canonical.tsv`.
- The separate `is_canonical` transcript attribute was resynchronised afterwards. A recheck found 40,220 attributes, one per gene with transcripts, each on the transcript the gene table names as canonical.
- The organelle sequence regions were swapped in that database (`Mt` held the chloroplast, `Pt` the mitochondrion) and were renamed on 2026-11-23. A recheck confirmed `Pt` is the chloroplast with the ATCG genes and `Mt` the mitochondrion with the ATMG genes.

## Checks

| Check | Outcome |
|---|---|
| Inputs | SHA256 of genome and GFF3 match; chromosome lengths match GenBank GCA_978657495.1 |
| Reformatted GFF | 26,867 genes, 36,572 mRNAs; exon and CDS lengths identical to the source for every transcript; no coordinate changed |
| Proteins | 36,572, none with an internal stop, all start with M |
| Pfam | 59,494 hits on 31,503 proteins (86.1%); 1,760 hits (3.0%) are families not integrated into InterPro, and they are counted |
| Each library | read count in the BAM equals ENA's; 96.9-99.3% of reads aligned |
| Strandedness | read 1 antisense in 98.6-99.9% of reads in every library, so `--rf` is right |
| Assemblies | 13,173 (pollen) to 29,186 transcripts; Chr1-5 only; at most 2 unstranded transcripts per sample; 187-528 antisense to an annotated gene |
| TRaCE run | 30 samples loaded, no warnings; 36,572 rows, 36 columns, ranks 1..n in every gene |
| Output columns | exon and CDS lengths match the source GFF3 for all transcripts; 31,503 transcripts have domain coverage, the same number as proteins with a Pfam hit |
| TRaCE itself | on the maize data bundled in its repository, lengths, CDS and domain columns reproduce exactly for all 72,539 transcripts |

Per-library numbers are in `qc_per_sample.tsv`, and the comparison of each assembly with TAIR12 is in `qc_gffcompare.tsv`.

## What to be careful about

### 1. 54 genes are decided by floating-point rounding

Each of the n tissues voting on a gene carries weight 6/n, and the three length voters total 6. When the tissues agree on one transcript and all three length criteria prefer another, the two sides should tie. In floating point, n copies of 6/n add up to slightly under 6 for n = 7, 10, 11, 17, 18, 20, 21, 22, and to slightly over 6 for n = 13, 15, 19, 23, 25-28, 30. The winner then depends on how many tissues voted, not on the evidence.

- Recounting every election with exact fractions changes the canonical transcript of 54 genes (0.9%); 53 of those are this rounding effect.
- In 16 of the 53, n is in the first set: TRaCE named the transcript the length criteria prefer, although the voting tissues prefer the other. In 37, n is in the second set: TRaCE named the transcript the tissues prefer, and exact arithmetic falls through to the alphabetical tie-break.
- These are effectively ties. `exact_recount.tsv` lists them, and the `exact_recount_canonical` column of the canonical table gives the alternative.
- The published output is unchanged; nothing was "corrected".

This affects the paper's maize design more directly, because 10 tissues is in the first set. Rerunning the bundled maize inputs with the current script, exact arithmetic changes the canonical of 99 of 15,162 multi-transcript genes, 83 of them with all 10 tissues voting. In `Zm00001eb033040`, for example, the elected transcript has no support in any tissue, and the one supported by all 10 loses.

### 2. Six genes depend on Perl's hash order

Rerunning with three other hash seeds changed the canonical of AT2G16920, AT2G22260, AT2G27350, AT2G34360, AT5G56850 and AT5G67580. The delivered run used seed 0 (`PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0`). Without a fixed seed, repeated runs of `TRaCE.pl` on the same inputs can differ on these genes.

### 3. 182 multi-transcript genes got no vote

Their ranking rests on domain coverage, CDS length and transcript length alone (`qc_silent_genes.tsv`).

- For 103, no assembled transcript above 0.5 TPM overlaps the gene in any tissue.
- For 79, an assembled transcript overlaps the gene in at least one tissue but never passed the 50% overlap thresholds. TRaCE scores only the highest-TPM transcript of each assembled locus, so a gene that shares a locus with a stronger neighbour can be passed over.

### 4. Two libraries carry a lot of ribosomal RNA

In S06_seed, 58% of mapped reads fall in annotated rRNA loci; in S20_silique_septum, 33%. The root and leaf libraries checked for comparison are at 0.3-0.4%. Col-CC v2 includes the rDNA arrays, so these reads now align. The assembled rRNA takes 32% and 12.5% of those samples' TPM, which makes the 0.5 TPM cutoff somewhat stricter there. Both still vote on about 14-15 thousand genes, against a typical 16 thousand.

### 5. Smaller points

- **Introns over 10 kb.** Alignment was capped at 10 kb introns. Five genes have a longer annotated intron (AT1G98163, AT4G08878, AT5G69529, AT5G69535, AT5G69554) and none got a vote. Only AT5G69554 has two transcripts.
- **Identical models.** AT3G09085, AT3G10455 and AT4G31080 each have two top transcripts that no criterion separates; the alphabetically first one is canonical.
- **Pollen** voted on only 6,471 genes, as expected for that tissue.
- **Pfam choice.** Counting only InterPro-integrated hits would change the canonical of 9 genes (`sensitivity/integrated_only.txt`).
- **Sample labels.** The submission's metadata disagrees on which of S01 and S02 is petal and which is sepal. It does not affect the election.
- **Maize rerun.** AED values differ from the bundled maize output in 34% of cells, mostly lower in the rerun, because that file was made by an earlier script revision (it has an extra score column). Rank 1 agrees for 85.5% of multi-transcript genes.

## Not done

- The pipeline loads nothing into a database. The update of the core database described above was made separately.
- `TRaCE.pl` was not modified.
- Non-coding genes are not ranked; TRaCE skips them.
- The atlas's 26 developmental-series libraries and long-read data were not used.
