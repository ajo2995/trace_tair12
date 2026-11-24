#!/usr/bin/env python3
"""Turn TRaCE's headerless output into the tables people will actually open.

Reads
  --trace         TRaCE.pl output: gene, transcript, rank, exon bp, CDS bp,
                  domain aa, then one AED column per sample
  --samples       GTF file names in the order they were given to TRaCE.pl
                  (TRaCE sorts sample names, so this list must be sorted)
  --gff           the TRaCE-format GFF (for strands and intron lengths)
  --stringtie-dir the per-sample StringTie GTFs
  --recount       optional: the .tsv written by recount_exact.py

Writes into --out-dir
  TAIR12_TRaCE_output.columns.txt   one column name per line
  TAIR12_canonical_transcripts.tsv  one row per gene
  qc_metrics.tsv                    headline numbers, one per line
  qc_voting_samples.tsv             histogram of voting samples per gene
  qc_per_sample.tsv                 per-sample alignment and evidence numbers
  qc_silent_genes.tsv               multi-transcript genes no sample voted on, with
                                    the number of samples in which an assembled
                                    transcript nevertheless overlaps the gene

A gene's "voting samples" are the samples in which at least one of its
transcripts has an AED at or under the cutoff; those are the only samples
that take part in its election.
"""
import argparse
import bisect
import collections
import os
import re
import sys


def read_trace(path, n_samples):
    genes = collections.OrderedDict()
    for line in open(path):
        cols = line.rstrip("\n").split("\t")
        if len(cols) != 6 + n_samples:
            sys.exit("%s: expected %d columns, found %d" % (cols[1], 6 + n_samples, len(cols)))
        genes.setdefault(cols[0], []).append({
            "id": cols[1], "rank": int(cols[2]), "len": int(cols[3]), "cds": int(cols[4]),
            "domain": int(cols[5]), "aed": [float(x) for x in cols[6:]]})
    for gene, transcripts in genes.items():
        transcripts.sort(key=lambda t: t["rank"])
        if [t["rank"] for t in transcripts] != list(range(1, len(transcripts) + 1)):
            sys.exit("%s: ranks are not 1..n" % gene)
    return genes


def read_gff(path, max_intron):
    """Return gene -> (seqid, strand, start, end) and the genes with an intron over max_intron."""
    span, parent, exons = {}, {}, collections.defaultdict(list)
    for line in open(path):
        cols = line.rstrip("\n").split("\t")
        attrs = dict(item.split("=", 1) for item in cols[8].split(";"))
        if cols[2] == "gene":
            span[attrs["ID"]] = (cols[0], cols[6], int(cols[3]), int(cols[4]))
        elif cols[2] == "mRNA":
            parent[attrs["ID"]] = attrs["Parent"]
        elif cols[2] == "exon":
            exons[attrs["Parent"]].append((int(cols[3]), int(cols[4])))
    long_intron = set()
    for transcript, segments in exons.items():
        segments.sort()
        if any(b[0] - a[1] - 1 > max_intron for a, b in zip(segments, segments[1:])):
            long_intron.add(parent[transcript])
    return span, long_intron


class AssembledLoci:
    """Assembled transcripts of one sample, for asking 'does anything overlap this gene?'"""

    def __init__(self, gtf_path, min_tpm):
        self.transcripts = 0
        self.loci = set()
        self.seqids = collections.Counter()
        self.unstranded = 0
        intervals = collections.defaultdict(list)
        for line in open(gtf_path):
            if line.startswith("#"):
                continue
            cols = line.split("\t", 9)
            if cols[2] != "transcript":
                continue
            self.transcripts += 1
            self.seqids[cols[0]] += 1
            if cols[6] == ".":
                self.unstranded += 1
            self.loci.add(re.search(r'gene_id "([^"]+)"', cols[8]).group(1))
            tpm = float(re.search(r'TPM "([^"]+)"', cols[8]).group(1))
            if tpm >= min_tpm:
                intervals[(cols[0], cols[6])].append((int(cols[3]), int(cols[4])))
        self.starts, self.max_end = {}, {}
        for key, items in intervals.items():
            items.sort()
            self.starts[key] = [start for start, _ in items]
            running, ends = 0, []
            for _, end in items:
                running = max(running, end)
                ends.append(running)
            self.max_end[key] = ends

    def overlaps(self, seqid, strand, start, end):
        key = (seqid, strand)
        if key not in self.starts:
            return False
        i = bisect.bisect_right(self.starts[key], end)      # transcripts starting at or before the gene's end
        return i > 0 and self.max_end[key][i - 1] >= start


def read_alignment_qc(proj, sample):
    """Collect the numbers the worker left behind for one sample."""
    row = {}
    summary = os.path.join(proj, "logs/hisat2", sample + ".summary")
    if os.path.exists(summary):
        text = open(summary).read()
        for label, key in (("Total pairs", "read_pairs"), ("Aligned concordantly 1 time", "concordant_unique"),
                           ("Aligned concordantly >1 times", "concordant_multi")):
            match = re.search(re.escape(label) + r":\s*(\d+)", text)
            if match:
                row[key] = int(match.group(1))
        match = re.search(r"Overall alignment rate:\s*([\d.]+)%", text)
        if match:
            row["overall_alignment_pct"] = float(match.group(1))
    strand = os.path.join(proj, "qc/strand", sample + ".txt")
    if os.path.exists(strand):
        for line in open(strand):
            key, value = line.split()
            if key == "RF_fraction":
                row["rf_fraction"] = float(value)
    idxstats = os.path.join(proj, "qc", "idxstats.%s.txt" % sample)
    if os.path.exists(idxstats):
        mapped = {cols[0]: int(cols[2]) for cols in (line.split("\t") for line in open(idxstats))}
        total = sum(mapped.values())
        if total:
            row["organelle_pct"] = 100.0 * (mapped.get("ChrC", 0) + mapped.get("ChrM", 0)) / total
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--trace", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--stringtie-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--recount")
    parser.add_argument("--proj", default=os.environ.get("PROJ", "."))
    parser.add_argument("--max-aed", type=float, default=0.5)
    parser.add_argument("--min-tpm", type=float, default=0.5)
    parser.add_argument("--max-intron", type=int, default=10000)
    args = parser.parse_args()

    samples = [line.strip()[:-len(".gtf")] for line in open(args.samples) if line.strip()]
    if samples != sorted(samples):
        sys.exit("sample list is not sorted, so it does not match TRaCE's column order")
    genes = read_trace(args.trace, len(samples))
    span, long_intron = read_gff(args.gff, args.max_intron)
    if set(genes) != set(span):
        sys.exit("TRaCE output and GFF disagree on the gene set")

    exact_top = {}
    if args.recount:
        with open(args.recount) as handle:
            header = handle.readline().rstrip("\n").split("\t")
            for line in handle:
                row = dict(zip(header, line.rstrip("\n").split("\t")))
                exact_top[row["gene_id"]] = row["exact_ranking"].split(",")[0]

    with open(os.path.join(args.out_dir, "TAIR12_TRaCE_output.columns.txt"), "w") as out:
        out.write("\n".join(["gene_id", "transcript_id", "rank", "exon_bp", "cds_bp", "domain_aa"]
                            + ["AED_" + s for s in samples]) + "\n")

    metrics = collections.OrderedDict()
    voting_hist = collections.Counter()
    per_sample_vote = [collections.Counter() for _ in samples]
    silent = []                         # multi-transcript genes that no sample voted on
    counts = collections.Counter()

    with open(os.path.join(args.out_dir, "TAIR12_canonical_transcripts.tsv"), "w") as out:
        out.write("\t".join(["gene_id", "canonical_transcript", "transcripts", "voting_samples",
                             "samples_supporting_canonical", "samples_where_canonical_is_best",
                             "is_longest_cds", "is_longest_transcript", "is_dot1",
                             "top_two_indistinguishable", "intron_over_10kb",
                             "exact_recount_canonical", "exon_bp", "cds_bp", "domain_aa"]) + "\n")
        for gene, transcripts in genes.items():
            top = transcripts[0]
            multi = len(transcripts) > 1
            voters = [i for i in range(len(samples)) if any(t["aed"][i] <= args.max_aed for t in transcripts)]
            supporting = sum(1 for i in voters if top["aed"][i] <= args.max_aed)
            best = sum(1 for i in voters if top["aed"][i] == min(t["aed"][i] for t in transcripts))
            longest_cds = top["cds"] == max(t["cds"] for t in transcripts)
            longest_len = top["len"] == max(t["len"] for t in transcripts)
            dot1 = top["id"].endswith(".1")
            same = multi and all(top[k] == transcripts[1][k] for k in ("len", "cds", "domain", "aed"))
            out.write("\t".join(str(x) for x in [
                gene, top["id"], len(transcripts), len(voters), supporting, best,
                "yes" if longest_cds else "no", "yes" if longest_len else "no", "yes" if dot1 else "no",
                "yes" if same else "no", "yes" if gene in long_intron else "no",
                exact_top.get(gene, top["id"]), top["len"], top["cds"], top["domain"]]) + "\n")

            for i in voters:
                per_sample_vote[i]["genes"] += 1
                if multi:
                    per_sample_vote[i]["multi"] += 1
            counts["genes"] += 1
            counts["transcripts"] += len(transcripts)
            counts["transcripts_with_domain"] += sum(1 for t in transcripts if t["domain"] > 0)
            if not multi:
                continue
            counts["multi"] += 1
            counts["multi_transcripts"] += len(transcripts)
            voting_hist[len(voters)] += 1
            if not voters:
                silent.append(gene)
            counts["multi_voted"] += bool(voters)
            counts["not_longest_cds"] += not longest_cds
            counts["not_longest_len"] += not longest_len
            counts["not_longest_either"] += (not longest_cds) and (not longest_len)
            counts["not_dot1"] += not dot1
            counts["indistinguishable"] += same
            counts["canonical_unsupported_but_voted"] += bool(voters) and supporting == 0
            counts["canonical_best_in_majority"] += bool(voters) and best * 2 > len(voters)
            counts["multi_long_intron"] += gene in long_intron
            counts["exact_differs"] += exact_top.get(gene, top["id"]) != top["id"]

    # Per-sample assemblies: sizes, and whether anything expressed overlaps the silent genes.
    overlapped = collections.Counter()
    with open(os.path.join(args.out_dir, "qc_per_sample.tsv"), "w") as out:
        columns = ["sample", "read_pairs", "overall_alignment_pct", "concordant_unique_pct", "organelle_pct",
                   "rf_fraction", "assembled_transcripts", "assembled_loci", "unstranded_transcripts",
                   "seqids", "genes_voted_on", "multi_transcript_genes_voted_on"]
        out.write("\t".join(columns) + "\n")
        for i, sample in enumerate(samples):
            loci = AssembledLoci(os.path.join(args.stringtie_dir, sample + ".gtf"), args.min_tpm)
            for gene in silent:
                if loci.overlaps(*span[gene]):
                    overlapped[gene] += 1
            row = read_alignment_qc(args.proj, sample)
            if "read_pairs" in row and "concordant_unique" in row:
                row["concordant_unique_pct"] = "%.2f" % (100.0 * row["concordant_unique"] / row["read_pairs"])
            row.update({"sample": sample, "assembled_transcripts": loci.transcripts,
                        "assembled_loci": len(loci.loci), "unstranded_transcripts": loci.unstranded,
                        "seqids": ",".join(sorted(loci.seqids)),
                        "genes_voted_on": per_sample_vote[i]["genes"],
                        "multi_transcript_genes_voted_on": per_sample_vote[i]["multi"]})
            for key in ("organelle_pct", "rf_fraction"):
                if key in row:
                    row[key] = "%.4f" % row[key]
            out.write("\t".join(str(row.get(c, "NA")) for c in columns) + "\n")

    with open(os.path.join(args.out_dir, "qc_voting_samples.tsv"), "w") as out:
        out.write("voting_samples\tmulti_transcript_genes\n")
        for n in range(len(samples) + 1):
            out.write("%d\t%d\n" % (n, voting_hist[n]))

    with open(os.path.join(args.out_dir, "qc_silent_genes.tsv"), "w") as out:
        out.write("gene_id\tsamples_with_an_overlapping_assembled_transcript\n")
        for gene in silent:
            out.write("%s\t%d\n" % (gene, overlapped[gene]))

    multi = counts["multi"]
    metrics["samples"] = len(samples)
    metrics["genes"] = counts["genes"]
    metrics["transcripts"] = counts["transcripts"]
    metrics["transcripts_with_pfam_domain"] = counts["transcripts_with_domain"]
    metrics["multi_transcript_genes"] = multi
    metrics["transcripts_in_multi_transcript_genes"] = counts["multi_transcripts"]
    metrics["multi_genes_with_at_least_one_voting_sample"] = counts["multi_voted"]
    metrics["multi_genes_with_no_voting_sample"] = len(silent)
    metrics["silent_genes_never_overlapped_by_an_assembled_transcript"] = sum(1 for g in silent if not overlapped[g])
    metrics["silent_genes_overlapped_in_at_least_one_sample"] = sum(1 for g in silent if overlapped[g])
    metrics["multi_genes_with_all_samples_voting"] = voting_hist[len(samples)]
    metrics["canonical_is_not_longest_cds"] = counts["not_longest_cds"]
    metrics["canonical_is_not_longest_transcript"] = counts["not_longest_len"]
    metrics["canonical_is_neither_longest_cds_nor_longest_transcript"] = counts["not_longest_either"]
    metrics["canonical_is_not_dot1"] = counts["not_dot1"]
    metrics["canonical_has_best_aed_in_majority_of_voting_samples"] = counts["canonical_best_in_majority"]
    metrics["canonical_unsupported_although_samples_voted"] = counts["canonical_unsupported_but_voted"]
    metrics["top_two_indistinguishable"] = counts["indistinguishable"]
    metrics["multi_genes_with_intron_over_%d" % args.max_intron] = counts["multi_long_intron"]
    if args.recount:
        metrics["canonical_differs_from_exact_recount"] = counts["exact_differs"]
    with open(os.path.join(args.out_dir, "qc_metrics.tsv"), "w") as out:
        for key, value in metrics.items():
            share = "" if key in ("samples", "genes", "transcripts", "transcripts_with_pfam_domain",
                                  "multi_transcript_genes", "transcripts_in_multi_transcript_genes") \
                else "%.1f%%" % (100.0 * value / multi if multi else 0)
            out.write("%s\t%d\t%s\n" % (key, value, share))
            print("%-62s %7d  %s" % (key, value, share))


if __name__ == "__main__":
    main()
