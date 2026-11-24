#!/usr/bin/env python3
"""Rewrite the TAIR12 GFF3 into the layout TRaCE.pl's parser expects.

TRaCE.pl (parseGFF) keeps a gene only when it carries biotype=protein_coding,
takes each transcript's start and end from its first and last exon line, and
drops an exon that appears before its mRNA line. TAIR12 tags genes with
locus_biotype= instead, and lists most minus-strand exons in descending order.

This script therefore writes, for every protein-coding gene:
    gene  (ID=...;biotype=protein_coding)
      mRNA  (ID=...;Parent=...)
        exons, ascending
        CDS segments, ascending
Coordinates, strands and CDS phases are copied unchanged. Everything else in
the source file (non-coding genes, pseudogenes, repeats) is left out, because
TRaCE ignores it.

Usage: make_trace_gff.py TAIR12.gff3[.gz] out.gff
"""
import collections
import gzip
import sys


def parse_attributes(text):
    """Return the key=value pairs of a GFF3 attribute column as a dict."""
    pairs = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            pairs[key] = value
    return pairs


def main(source_path, out_path):
    gene_ids = []                 # protein-coding genes, in file order
    gene_cols = {}                # gene ID -> first eight GFF columns
    transcripts_of = collections.defaultdict(list)   # gene ID -> mRNA IDs, in file order
    transcript_cols = {}
    segments = collections.defaultdict(lambda: {"exon": [], "CDS": []})
    skipped_mrna = 0

    opener = gzip.open if source_path.endswith(".gz") else open
    with opener(source_path, "rt") as source:
        for line in source:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9 or cols[2] not in ("gene", "mRNA", "exon", "CDS"):
                continue
            attrs = parse_attributes(cols[8])
            if cols[2] == "gene":
                if attrs.get("locus_biotype") == "protein_coding":
                    gene_ids.append(attrs["ID"])
                    gene_cols[attrs["ID"]] = cols[:8]
            elif cols[2] == "mRNA":
                if attrs["Parent"] in gene_cols:
                    transcripts_of[attrs["Parent"]].append(attrs["ID"])
                    transcript_cols[attrs["ID"]] = cols[:8]
                else:
                    skipped_mrna += 1
            else:
                segments[attrs["Parent"]][cols[2]].append(cols[:8])

    counts = collections.Counter()
    reordered = 0
    problems = []
    with open(out_path, "w") as out:
        for gene_id in gene_ids:
            out.write("\t".join(gene_cols[gene_id] + ["ID=%s;biotype=protein_coding" % gene_id]) + "\n")
            counts["gene"] += 1
            if not transcripts_of[gene_id]:
                problems.append("%s has no mRNA" % gene_id)
            for transcript_id in transcripts_of[gene_id]:
                out.write("\t".join(transcript_cols[transcript_id]
                                    + ["ID=%s;Parent=%s" % (transcript_id, gene_id)]) + "\n")
                counts["mRNA"] += 1
                for kind in ("exon", "CDS"):
                    rows = segments[transcript_id][kind]
                    ordered = sorted(rows, key=lambda r: int(r[3]))
                    if kind == "exon" and ordered != rows:
                        reordered += 1
                    if not ordered:
                        problems.append("%s has no %s" % (transcript_id, kind))
                    for previous, current in zip(ordered, ordered[1:]):
                        if int(current[3]) <= int(previous[4]):
                            problems.append("%s has overlapping %s segments" % (transcript_id, kind))
                    for row in ordered:
                        if (row[0], row[6]) != (transcript_cols[transcript_id][0], transcript_cols[transcript_id][6]):
                            problems.append("%s has a %s on another sequence or strand" % (transcript_id, kind))
                        out.write("\t".join(row + ["Parent=%s" % transcript_id]) + "\n")
                        counts[kind] += 1
                        counts[kind + "_bp"] += int(row[4]) - int(row[3]) + 1

    print("genes %d  mRNAs %d  exon lines %d  CDS lines %d" %
          (counts["gene"], counts["mRNA"], counts["exon"], counts["CDS"]), file=sys.stderr)
    print("exon bp %d  CDS bp %d" % (counts["exon_bp"], counts["CDS_bp"]), file=sys.stderr)
    print("transcripts whose exons were re-sorted: %d" % reordered, file=sys.stderr)
    print("mRNAs skipped because their gene is not protein-coding: %d" % skipped_mrna, file=sys.stderr)
    if problems:
        for problem in problems[:20]:
            print("PROBLEM: " + problem, file=sys.stderr)
        sys.exit("%d structural problems; not safe to pass to TRaCE" % len(problems))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
