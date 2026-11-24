#!/usr/bin/env python3
"""Confirm that an Ensembl core database's canonical transcripts are TRaCE's.

Run after fetch_ensembl_canonical.pl has refreshed ENSEMBL_DIR.

The loader (set_canonical_transcripts_from_TRaCE.pl) writes one thing: for
each rank-1 line of the TRaCE output it sets gene.canonical_transcript_id.
This script checks that, and three things around it:

  1. every gene TRaCE ranked now has TRaCE's rank-1 transcript as canonical;
  2. how many genes that changed, against the canonical each had before
     (BEFORE: gene, transcript, saved from the database before the update);
  3. genes TRaCE did not rank (non-coding, organelle) still have one
     canonical transcript of their own;
  4. the separate "is_canonical" transcript attribute, which the loader does
     not touch, still points at the same transcript as the gene table.

Usage: check_ensembl_trace_canonical.py ENSEMBL_DIR TRaCE_output.txt BEFORE.tsv
"""
import collections
import csv
import sys


def main(ensembl_dir, trace_path, before_path):
    trace_top, trace_genes = {}, collections.Counter()
    for line in open(trace_path):
        cols = line.split("\t", 3)
        trace_genes[cols[0]] += 1
        if cols[2] == "1":
            trace_top[cols[0]] = cols[1]
    before = dict(line.split()[:2] for line in open(before_path) if line.strip())

    by_gene = collections.defaultdict(list)
    with open(ensembl_dir + "/transcripts.tsv") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_gene[row["gene_id"]].append(row)

    flagged = {g: [r["transcript_id"] for r in rows if r["is_canonical"] == "1"] for g, rows in by_gene.items()}
    attrib = {g: [r["transcript_id"] for r in rows if r["is_canonical_attrib"] == "1"] for g, rows in by_gene.items()}

    print("== 1. genes ranked by TRaCE: %d ==" % len(trace_top))
    missing = [g for g in trace_top if g not in by_gene]
    match = [g for g in trace_top if flagged.get(g) == [trace_top[g]]]
    wrong = [g for g in trace_top if g in by_gene and flagged[g] != [trace_top[g]]]
    print("  canonical in the database is TRaCE's rank 1:  %d" % len(match))
    print("  canonical is something else:                  %d" % len(wrong))
    print("  gene not in the database:                     %d" % len(missing))
    for g in wrong[:10]:
        print("     %s database %s, TRaCE %s" % (g, ",".join(flagged[g]) or "none", trace_top[g]))
    multi = [g for g in trace_top if trace_genes[g] > 1]
    print("  of the %d with more than one transcript, matching: %d"
          % (len(multi), sum(1 for g in multi if flagged.get(g) == [trace_top[g]])))

    print("\n== 2. what the update changed, against the canonicals saved before it ==")
    known = [g for g in trace_top if g in before and g in by_gene]
    changed = [g for g in known if flagged[g] != [before[g]]]
    expected = [g for g in known if trace_top[g] != before[g]]
    print("  genes with a saved earlier canonical: %d" % len(known))
    print("  canonical differs from the earlier one: %d (TRaCE's rank 1 differed from it in %d)"
          % (len(changed), len(expected)))
    print("  changed although TRaCE agreed with the earlier canonical: %d"
          % len(set(changed) - set(expected)))
    print("  not changed although TRaCE disagreed: %d" % len(set(expected) - set(changed)))

    print("\n== 3. genes TRaCE did not rank ==")
    others = [g for g in by_gene if g not in trace_top]
    one = sum(1 for g in others if len(flagged[g]) == 1)
    print("  genes with transcripts: %d; with exactly one canonical of their own: %d; with none: %d; with several: %d"
          % (len(others), one, sum(1 for g in others if not flagged[g]), sum(1 for g in others if len(flagged[g]) > 1)))
    biotypes = collections.Counter(by_gene[g][0]["gene_biotype"] + (" (organelle)" if by_gene[g][0]["seq_region"] in ("Mt", "Pt") else "")
                                   for g in others)
    print("  by biotype: " + ", ".join("%d %s" % (n, b) for b, n in biotypes.most_common(8)))
    same_as_attrib = sum(1 for g in others if flagged[g] == attrib[g])
    print("  canonical equals the transcript carrying the is_canonical attribute: %d of %d" % (same_as_attrib, len(others)))

    print("\n== 4. the is_canonical transcript attribute ==")
    total_attrib = sum(len(v) for v in attrib.values())
    print("  transcripts carrying it: %d; genes with exactly one such transcript: %d of %d"
          % (total_attrib, sum(1 for v in attrib.values() if len(v) == 1), len(by_gene)))
    stale = [g for g in by_gene if attrib[g] != flagged[g]]
    stale_ranked = [g for g in stale if g in trace_top]
    print("  genes where it marks a different transcript than the gene table: %d (%d of them ranked by TRaCE)"
          % (len(stale), len(stale_ranked)))
    if stale:
        old = sum(1 for g in stale_ranked if g in before and attrib[g] == [before[g]])
        print("  in %d of those the attribute still marks the canonical from before the update" % old)
        for g in stale[:5]:
            print("     %s gene table %s, attribute %s" % (g, ",".join(flagged[g]) or "none", ",".join(attrib[g]) or "none"))

    print("\n== verdict ==")
    if wrong or missing:
        print("  PROBLEM: %d ranked genes do not have TRaCE's canonical" % (len(wrong) + len(missing)))
    else:
        print("  gene.canonical_transcript_id is TRaCE's rank 1 for all %d ranked genes." % len(trace_top))
    if stale:
        print("  The is_canonical attribute disagrees with the gene table for %d genes." % len(stale))
    else:
        print("  The is_canonical attribute agrees with the gene table for every gene.")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
