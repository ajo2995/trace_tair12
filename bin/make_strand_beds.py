#!/usr/bin/env python3
"""Write the two BED files that strand_check.sh counts reads over.

Each file holds the exons of Chr1 protein-coding genes on one strand. Genes
whose span touches any annotated gene, pseudogene or transposable element on
the opposite strand are left out, so that a read falling in these exons can
only have come from a transcript of known strand.

Usage: make_strand_beds.py TAIR12.gff3.gz TAIR12.trace.gff out_prefix
Writes out_prefix.plus.bed and out_prefix.minus.bed.
"""
import collections
import gzip
import sys

CHROMOSOME = "Chr1"
LOCUS_TYPES = ("gene", "pseudogene", "mobile_element")


def main(source_gff, trace_gff, prefix):
    # Every annotated locus on the chromosome, coding or not.
    loci = []
    with gzip.open(source_gff, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) == 9 and cols[0] == CHROMOSOME and cols[2] in LOCUS_TYPES:
                loci.append((int(cols[3]), int(cols[4]), cols[6]))
    loci.sort()

    gene_span = {}
    gene_of = {}
    exons = collections.defaultdict(list)
    for line in open(trace_gff):
        cols = line.rstrip("\n").split("\t")
        if cols[0] != CHROMOSOME:
            continue
        attrs = dict(item.split("=", 1) for item in cols[8].split(";"))
        if cols[2] == "gene":
            gene_span[attrs["ID"]] = (int(cols[3]), int(cols[4]), cols[6])
        elif cols[2] == "mRNA":
            gene_of[attrs["ID"]] = attrs["Parent"]
        elif cols[2] == "exon":
            exons[gene_of[attrs["Parent"]]].append((int(cols[3]), int(cols[4])))

    kept = collections.Counter()
    out = {"+": open(prefix + ".plus.bed", "w"), "-": open(prefix + ".minus.bed", "w")}
    for gene_id, (start, end, strand) in sorted(gene_span.items(), key=lambda item: item[1]):
        if any(s <= end and e >= start and other != strand for s, e, other in loci):
            continue
        kept[strand] += 1
        # Merge the exons of all isoforms so no base is listed twice.
        merged = []
        for exon_start, exon_end in sorted(exons[gene_id]):
            if merged and exon_start <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], exon_end)
            else:
                merged.append([exon_start, exon_end])
        for exon_start, exon_end in merged:
            out[strand].write("%s\t%d\t%d\t%s\n" % (CHROMOSOME, exon_start - 1, exon_end, gene_id))
    for handle in out.values():
        handle.close()
    print("%s coding genes: %d; kept with no opposite-strand neighbour: %d plus, %d minus"
          % (CHROMOSOME, len(gene_span), kept["+"], kept["-"]), file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
