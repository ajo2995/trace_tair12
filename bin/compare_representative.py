#!/usr/bin/env python3
"""Compare TRaCE's ranking with TAIR12's representative gene models.

The representative of each gene is read from a file (--representatives:
gene, transcript). For TAIR12 that file is derived from TAIR's own list,
TAIR12_representative.txt. Without a file the script applies TAIR's rule
itself: the longest transcript and, when two are equally long, the one with
the lower isoform number, where "longest" is spliced length (the exon-length
column TRaCE reports). That rule reproduces TAIR's list for every
multi-transcript coding gene. --by-span uses genomic span instead, which
does not.

The same file option serves other "one transcript per gene" sets, such as the
canonical transcripts of the Ensembl core database built from TAIR12. Those
follow Ensembl's rule (longest CDS first) and are a different set: they
differ from TAIR's representatives in 870 genes.

For every gene the script finds that representative and looks up where TRaCE
ranked it. Where TRaCE put another transcript first, it measures how the two
transcripts differ:
  length      spliced length, CDS length, Pfam coverage
  protein     whether the two encode the same protein sequence
  structure   whether they share one intron chain (so differ only at their
              ends) and what fraction of their exonic bases they share
  evidence    in how many voting tissues each is closer to the assembled
              transcript (lower AED), and each one's mean AED

A gene's voting tissues are those in which at least one of its transcripts
has an AED at or under the cutoff. An unsupported transcript has AED 1.

Usage: compare_representative.py --trace OUT.txt --gff TRACE.gff --proteins PEP.fa
                                 --out TABLE.tsv [--recount RECOUNT.tsv] [--by-span]
"""
import argparse
import collections
import statistics


def isoform_number(transcript_id):
    return int(transcript_id.rsplit(".", 1)[1])


def read_trace(path):
    genes = collections.OrderedDict()
    for line in open(path):
        cols = line.rstrip("\n").split("\t")
        genes.setdefault(cols[0], {})[cols[1]] = {
            "rank": int(cols[2]), "len": int(cols[3]), "cds": int(cols[4]),
            "domain": int(cols[5]), "aed": [float(x) for x in cols[6:]]}
    return genes


def read_structures(path):
    exons = collections.defaultdict(list)
    for line in open(path):
        cols = line.rstrip("\n").split("\t")
        if cols[2] == "exon":
            exons[cols[8][len("Parent="):]].append((int(cols[3]), int(cols[4])))
    for segments in exons.values():
        segments.sort()
    return exons


def read_fasta(path):
    sequences, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            sequences[name] = []
        elif name:
            sequences[name].append(line)
    return {name: "".join(parts) for name, parts in sequences.items()}


def shared_bases(a, b):
    """Number of bases covered by both lists of sorted, non-overlapping intervals."""
    total, i, j = 0, 0, 0
    while i < len(a) and j < len(b):
        low, high = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if low <= high:
            total += high - low + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total


def intron_chain(segments):
    return tuple((left[1], right[0]) for left, right in zip(segments, segments[1:]))


def retains_intron(exons_a, exons_b):
    """True when an intron of transcript b lies wholly inside one exon of transcript a."""
    introns_b = [(left[1] + 1, right[0] - 1) for left, right in zip(exons_b, exons_b[1:])]
    return any(exon[0] <= start and end <= exon[1] for start, end in introns_b for exon in exons_a)


def describe(values):
    values = sorted(values)
    quartiles = statistics.quantiles(values, n=4) if len(values) > 1 else [values[0]] * 3
    return "median %g (quartiles %g to %g, max %g)" % (statistics.median(values), quartiles[0], quartiles[2], values[-1])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--trace", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--recount")
    parser.add_argument("--by-span", action="store_true")
    parser.add_argument("--representatives",
                        help="two-column file (gene, transcript) naming each gene's representative; "
                             "used instead of the longest-transcript rule")
    parser.add_argument("--max-aed", type=float, default=0.5)
    args = parser.parse_args()

    genes = read_trace(args.trace)
    exons = read_structures(args.gff)
    proteins = read_fasta(args.proteins)
    rounding = set()
    if args.recount:
        with open(args.recount) as handle:
            header = handle.readline().rstrip("\n").split("\t")
            for line in handle:
                row = dict(zip(header, line.rstrip("\n").split("\t")))
                if row["rank1_changed"] == "yes":
                    rounding.add(row["gene_id"])

    given = {}
    if args.representatives:
        for line in open(args.representatives):
            gene, transcript = line.split()[:2]
            given[gene] = transcript

    def length(gene, transcript):
        if args.by_span:
            return exons[transcript][-1][1] - exons[transcript][0][0] + 1
        return genes[gene][transcript]["len"]

    columns = ["gene_id", "transcripts", "representative", "representative_trace_rank", "trace_canonical",
               "agree", "voting_tissues", "representative_bp", "canonical_bp", "bp_longer_representative",
               "representative_cds_bp", "canonical_cds_bp", "same_protein", "protein_aa_canonical_minus_representative",
               "representative_domain_aa", "canonical_domain_aa", "same_intron_chain",
               "representative_retains_intron", "canonical_retains_intron", "shared_exonic_fraction",
               "tissues_canonical_closer", "tissues_representative_closer", "tissues_equal",
               "tissues_supporting_canonical", "tissues_supporting_representative",
               "mean_aed_canonical", "mean_aed_representative", "decided_by_rounding"]
    rows = []
    single = 0
    for gene, transcripts in genes.items():
        if len(transcripts) == 1:
            single += 1
            continue
        if given:
            representative = given[gene]
            if representative not in transcripts:
                raise SystemExit("%s: representative %s is not one of its transcripts" % (gene, representative))
        else:
            representative = min(transcripts, key=lambda t: (-length(gene, t), isoform_number(t)))
        canonical = min(transcripts, key=lambda t: transcripts[t]["rank"])
        rep, can = transcripts[representative], transcripts[canonical]
        voters = [i for i in range(len(rep["aed"]))
                  if any(t["aed"][i] <= args.max_aed for t in transcripts.values())]
        union = rep["len"] + can["len"] - shared_bases(exons[representative], exons[canonical])
        row = {
            "gene_id": gene, "transcripts": len(transcripts), "representative": representative,
            "representative_trace_rank": rep["rank"], "trace_canonical": canonical,
            "agree": "yes" if representative == canonical else "no", "voting_tissues": len(voters),
            "representative_bp": rep["len"], "canonical_bp": can["len"],
            "bp_longer_representative": rep["len"] - can["len"],
            "representative_cds_bp": rep["cds"], "canonical_cds_bp": can["cds"],
            "same_protein": "yes" if proteins[representative] == proteins[canonical] else "no",
            "protein_aa_canonical_minus_representative": len(proteins[canonical]) - len(proteins[representative]),
            "representative_domain_aa": rep["domain"], "canonical_domain_aa": can["domain"],
            "same_intron_chain": "yes" if intron_chain(exons[representative]) == intron_chain(exons[canonical]) else "no",
            # An intron of one transcript that the other reads straight through.
            "representative_retains_intron": "yes" if retains_intron(exons[representative], exons[canonical]) else "no",
            "canonical_retains_intron": "yes" if retains_intron(exons[canonical], exons[representative]) else "no",
            "shared_exonic_fraction": "%.4f" % ((rep["len"] + can["len"] - union) / union),
            "tissues_canonical_closer": sum(1 for i in voters if can["aed"][i] < rep["aed"][i]),
            "tissues_representative_closer": sum(1 for i in voters if rep["aed"][i] < can["aed"][i]),
            "tissues_equal": sum(1 for i in voters if rep["aed"][i] == can["aed"][i]),
            "tissues_supporting_canonical": sum(1 for i in voters if can["aed"][i] <= args.max_aed),
            "tissues_supporting_representative": sum(1 for i in voters if rep["aed"][i] <= args.max_aed),
            "mean_aed_canonical": "%.4f" % statistics.mean(can["aed"][i] for i in voters) if voters else "NA",
            "mean_aed_representative": "%.4f" % statistics.mean(rep["aed"][i] for i in voters) if voters else "NA",
            "decided_by_rounding": "yes" if gene in rounding else "no",
        }
        rows.append(row)

    with open(args.out, "w") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(str(row[c]) for c in columns) + "\n")

    differ = [r for r in rows if r["agree"] == "no"]
    multi = len(rows)

    def share(n, of):
        return "%d (%.1f%%)" % (n, 100.0 * n / of if of else 0)

    if given:
        print("representatives taken from: %s" % args.representatives)
    else:
        print("definition of longest: %s" % ("genomic span" if args.by_span else "spliced length (sum of exons)"))
    print("genes: %d; single-transcript (representative is rank 1 by default): %d; multi-transcript: %d"
          % (single + multi, single, multi))
    print("representative is TRaCE rank 1:              %s of multi-transcript genes" % share(multi - len(differ), multi))
    print("representative ranked below another transcript: %s of multi-transcript genes, %.1f%% of all genes"
          % (share(len(differ), multi), 100.0 * len(differ) / (single + multi)))
    last = sum(1 for r in differ if r["representative_trace_rank"] == r["transcripts"])
    ranks = collections.Counter(min(r["representative_trace_rank"], 5) for r in differ)
    print("  its TRaCE rank: " + ", ".join("%s: %d" % ("5+" if k == 5 else k, ranks[k]) for k in sorted(ranks))
          + "; ranked last in its gene: %d" % last)
    two = [r for r in rows if r["transcripts"] == 2]
    print("  in two-transcript genes: %s; in genes with three or more: %s"
          % (share(sum(1 for r in two if r["agree"] == "no"), len(two)),
             share(sum(1 for r in rows if r["transcripts"] > 2 and r["agree"] == "no"), multi - len(two))))

    n = len(differ)
    print("\nwhere they differ (%d genes), canonical against representative:" % n)
    gaps = [r["bp_longer_representative"] for r in differ]
    print("  spliced length, representative minus canonical (bp): " + describe(gaps))
    edges = [(0, 0), (1, 49), (50, 199), (200, 499), (500, 10 ** 9)]
    print("    " + "; ".join("%s: %s" % ("same length" if hi == 0 else ("%d+" % lo if hi > 10 ** 8 else "%d-%d" % (lo, hi)),
                                            share(sum(1 for g in gaps if lo <= g <= hi), n)) for lo, hi in edges))
    relative = [100.0 * r["bp_longer_representative"] / r["representative_bp"] for r in differ]
    print("    as a share of the representative's length (%): " + describe([round(x, 1) for x in relative]))

    same = [r for r in differ if r["same_protein"] == "yes"]
    other = [r for r in differ if r["same_protein"] == "no"]
    print("  same protein: %s; different protein: %s" % (share(len(same), n), share(len(other), n)))
    longer = [r for r in other if r["protein_aa_canonical_minus_representative"] > 0]
    shorter = [r for r in other if r["protein_aa_canonical_minus_representative"] < 0]
    print("    different protein: canonical longer in %d, shorter in %d, same length in %d"
          % (len(longer), len(shorter), len(other) - len(longer) - len(shorter)))
    if longer:
        print("      canonical longer by (aa): " + describe([r["protein_aa_canonical_minus_representative"] for r in longer]))
    if shorter:
        print("      canonical shorter by (aa): " + describe([-r["protein_aa_canonical_minus_representative"] for r in shorter]))
    more = sum(1 for r in differ if r["canonical_domain_aa"] > r["representative_domain_aa"])
    less = sum(1 for r in differ if r["canonical_domain_aa"] < r["representative_domain_aa"])
    print("  Pfam coverage: canonical more in %s, equal in %s, less in %s" % (share(more, n), share(n - more - less, n), share(less, n)))

    chain = sum(1 for r in differ if r["same_intron_chain"] == "yes")
    print("  same intron chain, different ends: %s; different intron chain: %s" % (share(chain, n), share(n - chain, n)))
    for label, subset in (("same protein", same), ("different protein", other)):
        if subset:
            c = sum(1 for r in subset if r["same_intron_chain"] == "yes")
            print("    %s: same intron chain %d, different %d" % (label, c, len(subset) - c))
    rep_keeps = [r for r in differ if r["representative_retains_intron"] == "yes"]
    can_keeps = [r for r in differ if r["canonical_retains_intron"] == "yes"]
    both = sum(1 for r in rep_keeps if r["canonical_retains_intron"] == "yes")
    print("  representative reads through an intron the canonical splices out: %s; the reverse: %s; both: %d"
          % (share(len(rep_keeps), n), share(len(can_keeps), n), both))
    print("    of those %d: same protein %d, different protein %d"
          % (len(rep_keeps), sum(1 for r in rep_keeps if r["same_protein"] == "yes"),
             sum(1 for r in rep_keeps if r["same_protein"] == "no")))
    print("  exonic bases shared (intersection over union): "
          + describe([round(float(r["shared_exonic_fraction"]), 3) for r in differ]))

    voted = [r for r in differ if r["voting_tissues"] > 0]
    print("  evidence, in the %d of these genes with at least one voting tissue:" % len(voted))
    can_wins = sum(1 for r in voted if r["tissues_canonical_closer"] > r["tissues_representative_closer"])
    rep_wins = sum(1 for r in voted if r["tissues_representative_closer"] > r["tissues_canonical_closer"])
    print("    canonical closer in more tissues: %s; representative closer in more: %s; level: %s"
          % (share(can_wins, len(voted)), share(rep_wins, len(voted)), share(len(voted) - can_wins - rep_wins, len(voted))))
    unanimous = sum(1 for r in voted if r["tissues_canonical_closer"] == r["voting_tissues"])
    print("    canonical closer in every voting tissue: %s" % share(unanimous, len(voted)))
    tissue_can = sum(r["tissues_canonical_closer"] for r in voted)
    tissue_rep = sum(r["tissues_representative_closer"] for r in voted)
    tissue_all = sum(r["voting_tissues"] for r in voted)
    print("    over all %d gene-tissue votes: canonical closer %.1f%%, representative closer %.1f%%, equal %.1f%%"
          % (tissue_all, 100.0 * tissue_can / tissue_all, 100.0 * tissue_rep / tissue_all,
             100.0 * (tissue_all - tissue_can - tissue_rep) / tissue_all))
    print("    mean AED per gene, canonical:      " + describe([float(r["mean_aed_canonical"]) for r in voted]))
    print("    mean AED per gene, representative: " + describe([float(r["mean_aed_representative"]) for r in voted]))
    print("    difference, representative minus canonical: "
          + describe([round(float(r["mean_aed_representative"]) - float(r["mean_aed_canonical"]), 4) for r in voted]))
    never = sum(1 for r in voted if r["tissues_supporting_representative"] == 0)
    print("    representative matches no assembled transcript in any tissue: %s" % share(never, len(voted)))
    print("    tissues supporting each, median: canonical %g, representative %g"
          % (statistics.median(r["tissues_supporting_canonical"] for r in voted),
             statistics.median(r["tissues_supporting_representative"] for r in voted)))
    print("  no voting tissue, so decided by domain and CDS length alone: %d" % (n - len(voted)))
    print("  among the rounding-decided genes: %d" % sum(1 for r in differ if r["decided_by_rounding"] == "yes"))


if __name__ == "__main__":
    main()
