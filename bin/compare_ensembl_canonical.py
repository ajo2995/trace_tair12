#!/usr/bin/env python3
"""Work out which rule a set of "one transcript per gene" choices follows.

The choices tested are either the canonical transcripts of an Ensembl core
database (the default) or a list supplied with a fifth argument, such as
TAIR's own file of representative gene models. Gene structures come from the
tables written by fetch_ensembl_canonical.pl.

For the nuclear protein-coding genes the script:

  1. checks that the database holds the same genes, transcripts, exons and
     CDS as the TAIR12 GFF, so that lengths can be compared at all;
  2. if a list is given, checks it against the gene set and against the
     database's canonical transcripts;
  3. for every gene with more than one transcript, asks which transcript each
     candidate rule would pick, and how often that is the chosen transcript;
  4. writes one row per multi-transcript gene with the chosen transcript,
     each rule's pick, and where TRaCE ranked the chosen transcript.

Candidate rules (every rule ends with "lower isoform number" unless it says
otherwise, so that it always picks exactly one transcript):
  spliced             longest spliced transcript (sum of exons)
  span                longest genomic span
  cds                 longest CDS
  cds_spliced         longest CDS, then longest spliced transcript
  cds_span            longest CDS, then longest genomic span
  cds_spliced_string  as cds_spliced, but final ties by ID compared as text
                      (so .10 sorts before .2), as Ensembl's own selector does
  spliced_cds         longest spliced transcript, then longest CDS
  spliced_string      longest spliced transcript, final ties by ID as text
  lowest_number       the lowest isoform number

Usage: compare_ensembl_canonical.py ENSEMBL_DIR TRACE.gff TRaCE_output.txt OUT_PREFIX [LIST]
  LIST  one transcript ID per line (the last whitespace-separated field is used)
Writes OUT_PREFIX.tsv and OUT_PREFIX.canonical.tsv (gene, chosen transcript).
"""
import collections
import csv
import sys

NUCLEAR = {"1", "2", "3", "4", "5"}


def isoform_number(transcript_id):
    return int(transcript_id.rsplit(".", 1)[1])


RULES = collections.OrderedDict([
    ("spliced", lambda t: (-t["len"], isoform_number(t["id"]))),
    ("span", lambda t: (-t["span"], isoform_number(t["id"]))),
    ("cds", lambda t: (-t["cds"], isoform_number(t["id"]))),
    ("cds_spliced", lambda t: (-t["cds"], -t["len"], isoform_number(t["id"]))),
    ("cds_span", lambda t: (-t["cds"], -t["span"], isoform_number(t["id"]))),
    ("cds_spliced_string", lambda t: (-t["cds"], -t["len"], t["id"])),
    ("spliced_cds", lambda t: (-t["len"], -t["cds"], isoform_number(t["id"]))),
    ("spliced_string", lambda t: (-t["len"], t["id"])),
    ("lowest_number", lambda t: (isoform_number(t["id"]),)),
])


def read_table(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main(ensembl_dir, gff_path, trace_path, out_prefix, list_path=None):
    # ---- the database ---------------------------------------------------
    db_gene = {}                    # transcript -> gene (nuclear protein-coding)
    all_gene = {}                   # transcript -> gene (everything in the database)
    canonical = {}                  # gene -> the database's canonical transcript
    genes = collections.defaultdict(list)
    for row in read_table(ensembl_dir + "/transcripts.tsv"):
        all_gene[row["transcript_id"]] = row["gene_id"]
        if row["gene_biotype"] != "protein_coding" or row["seq_region"] not in NUCLEAR:
            continue
        db_gene[row["transcript_id"]] = row["gene_id"]
        genes[row["gene_id"]].append(row["transcript_id"])
        if row["is_canonical"] == "1":
            canonical[row["gene_id"]] = row["transcript_id"]

    db_exons = collections.defaultdict(list)       # transcript -> [(rank, dbid, start, end)]
    for row in read_table(ensembl_dir + "/exons.tsv"):
        if row["transcript_id"] in db_gene:
            db_exons[row["transcript_id"]].append(
                (int(row["rank"]), row["exon_dbid"], int(row["exon_start"]), int(row["exon_end"])))
    db_cds = {}
    for row in read_table(ensembl_dir + "/translations.tsv"):
        transcript = row["transcript_id"]
        if transcript not in db_gene:
            continue
        exons = sorted(db_exons[transcript])
        rank_of = {dbid: rank for rank, dbid, _, _ in exons}
        first, last = rank_of[row["start_exon_dbid"]], rank_of[row["end_exon_dbid"]]
        length = {rank: end - start + 1 for rank, _, start, end in exons}
        db_cds[transcript] = (sum(length[r] for r in range(first, last + 1))
                              - (int(row["seq_start"]) - 1) - (length[last] - int(row["seq_end"])))

    # ---- TAIR12 as given to TRaCE ---------------------------------------
    gff_exons = collections.defaultdict(list)
    gff_gene = {}
    for line in open(gff_path):
        cols = line.rstrip("\n").split("\t")
        attrs = dict(item.split("=", 1) for item in cols[8].split(";"))
        if cols[2] == "mRNA":
            gff_gene[attrs["ID"]] = attrs["Parent"]
        elif cols[2] == "exon":
            gff_exons[attrs["Parent"]].append((int(cols[3]), int(cols[4])))
    trace = {}
    for line in open(trace_path):
        cols = line.split("\t", 6)
        trace[cols[1]] = {"gene": cols[0], "rank": int(cols[2]), "len": int(cols[3]), "cds": int(cols[4])}

    print("== does the database hold the same gene models as the TAIR12 GFF? ==")
    print("nuclear protein-coding genes: database %d, GFF %d, shared %d"
          % (len(genes), len(set(gff_gene.values())), len(set(genes) & set(gff_gene.values()))))
    print("their transcripts:            database %d, GFF %d, shared %d"
          % (len(db_gene), len(gff_gene), len(set(db_gene) & set(gff_gene))))
    shared = sorted(set(db_gene) & set(gff_gene))
    same_gene = sum(1 for t in shared if db_gene[t] == gff_gene[t])
    same_exons = sum(1 for t in shared if sorted((s, e) for _, _, s, e in db_exons[t]) == sorted(gff_exons[t]))
    cds_delta = collections.Counter(trace[t]["cds"] - db_cds.get(t, 0) for t in shared)
    print("shared transcripts with the same parent gene: %d; identical exon coordinates: %d of %d"
          % (same_gene, same_exons, len(shared)))
    print("CDS length, GFF minus database (bp -> transcripts): %s" % dict(cds_delta.most_common(4)))
    print("genes with a canonical transcript in the database: %d of %d" % (len(canonical), len(genes)))

    # ---- the choices to test ---------------------------------------------
    chosen, what = canonical, "the database's canonical"
    if list_path:
        what = "the listed transcript"
        listed = [line.split()[-1] for line in open(list_path) if line.strip()]
        per_gene = collections.defaultdict(list)
        unknown = 0
        for transcript in listed:
            if transcript in all_gene:
                per_gene[all_gene[transcript]].append(transcript)
            else:
                unknown += 1
        doubled = sum(1 for ts in per_gene.values() if len(ts) > 1)
        print("\n== the list %s ==" % list_path)
        print("entries: %d; distinct: %d; not a transcript in the database: %d" % (len(listed), len(set(listed)), unknown))
        print("database genes covered: %d of %d with transcripts; genes listed more than once: %d"
              % (len(per_gene), len(set(all_gene.values())), doubled))
        chosen = {g: ts[0] for g, ts in per_gene.items() if g in genes}
        print("nuclear protein-coding genes covered: %d of %d" % (len(chosen), len(genes)))
        both = [g for g in chosen if g in canonical]
        same = sum(1 for g in both if chosen[g] == canonical[g])
        multi_both = [g for g in both if len(genes[g]) > 1]
        same_multi = sum(1 for g in multi_both if chosen[g] == canonical[g])
        print("same transcript as the database's canonical: %d of %d genes; among multi-transcript genes %d of %d (%.2f%%)"
              % (same, len(both), same_multi, len(multi_both), 100.0 * same_multi / len(multi_both)))

    # ---- which rule gives the chosen transcript? -------------------------
    multi = {g: sorted(ts, key=isoform_number) for g, ts in genes.items() if len(ts) > 1 and g in chosen}
    table = {}
    for gene, transcripts in multi.items():
        info = []
        for t in transcripts:
            starts = [s for _, _, s, _ in db_exons[t]]
            ends = [e for _, _, _, e in db_exons[t]]
            info.append({"id": t, "len": sum(e - s + 1 for _, _, s, e in db_exons[t]),
                         "span": max(ends) - min(starts) + 1, "cds": db_cds[t]})
        table[gene] = info

    print("\n== multi-transcript genes: %d ==" % len(multi))
    print("how often each rule picks %s:" % what)
    picks = {}
    for name, key in RULES.items():
        picks[name] = {gene: min(info, key=key)["id"] for gene, info in table.items()}
        agree = sum(1 for gene in multi if picks[name][gene] == chosen[gene])
        print("  %-20s %5d  %6.2f%%" % (name, agree, 100.0 * agree / len(multi)))

    best = max(RULES, key=lambda name: sum(1 for g in multi if picks[name][g] == chosen[g]))
    misses = [g for g in multi if picks[best][g] != chosen[g]]
    print("\nbest rule: %s; genes it does not explain: %d" % (best, len(misses)))
    for gene in misses[:25]:
        print("  %s chosen %s, rule picks %s | %s" % (
            gene, chosen[gene], picks[best][gene],
            "; ".join("%s cds=%d len=%d span=%d" % (t["id"].rsplit(".", 1)[1], t["cds"], t["len"], t["span"])
                      for t in table[gene])))

    # how the chosen transcript compares with the longest values in its gene
    c = collections.Counter()
    for gene, info in table.items():
        pick = next(t for t in info if t["id"] == chosen[gene])
        c["longest_cds"] += pick["cds"] == max(t["cds"] for t in info)
        c["longest_len"] += pick["len"] == max(t["len"] for t in info)
        c["longest_span"] += pick["span"] == max(t["span"] for t in info)
        c["dot1"] += chosen[gene].endswith(".1")
    print("\n%s is (ties counted as yes):" % what)
    for key, label in (("longest_cds", "a longest CDS"), ("longest_len", "a longest spliced transcript"),
                       ("longest_span", "a longest genomic span"), ("dot1", "the .1 model")):
        print("  %-30s %5d  %6.2f%%" % (label, c[key], 100.0 * c[key] / len(multi)))

    # ---- the chosen transcripts against TRaCE ----------------------------
    below = [g for g in multi if trace[chosen[g]]["rank"] > 1]
    print("\n== TRaCE against these choices ==")
    print("TRaCE ranks %s below another transcript: %d of %d multi-transcript genes (%.1f%%)"
          % (what, len(below), len(multi), 100.0 * len(below) / len(multi)))
    ranks = collections.Counter(min(trace[chosen[g]]["rank"], 4) for g in below)
    print("  its TRaCE rank: " + ", ".join("%s: %d" % ("4+" if k == 4 else k, ranks[k]) for k in sorted(ranks)))

    with open(out_prefix + ".tsv", "w") as out:
        out.write("\t".join(["gene_id", "transcripts", "chosen", "ensembl_canonical", "trace_canonical",
                             "trace_rank_of_chosen"] + ["rule_" + name for name in RULES]) + "\n")
        for gene in multi:
            trace_top = next(t for t in multi[gene] if trace[t]["rank"] == 1)
            out.write("\t".join([gene, str(len(multi[gene])), chosen[gene], canonical[gene], trace_top,
                                 str(trace[chosen[gene]]["rank"])]
                                + [picks[name][gene] for name in RULES]) + "\n")
    with open(out_prefix + ".canonical.tsv", "w") as out:
        for gene in sorted(chosen):
            out.write("%s\t%s\n" % (gene, chosen[gene]))


if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        sys.exit(__doc__)
    main(*sys.argv[1:])
