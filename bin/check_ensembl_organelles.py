#!/usr/bin/env python3
"""Check that the organelle regions of an Ensembl core database are labelled correctly.

Run after fetch_ensembl_canonical.pl and inspect_ensembl_organelles.pl have
refreshed the tables and the two region sequences in ENSEMBL_DIR.

For the regions named Pt (plastid) and Mt (mitochondrion) it reports:
  sequence      which real organelle genome the region's sequence equals
  genes         whether the ATCG / ATMG gene IDs sit on the matching region
  translations  whether the proteins translate without internal stops from
                the region they sit on (standard genetic code)
  labels        whether the karyotype_rank value and the names of the
                chunk-level pieces agree with the region's name

Usage: check_ensembl_organelles.py ENSEMBL_DIR ORGANELLES.fa
  ORGANELLES.fa must hold ChrC (chloroplast AP000423.1) and ChrM
  (mitochondrion BK010421.1).
"""
import collections
import csv
import re
import sys

EXPECTED = {
    "Pt": {"reference": "ChrC", "what": "chloroplast AP000423.1", "prefix": "ATCG"},
    "Mt": {"reference": "ChrM", "what": "mitochondrion BK010421.1", "prefix": "ATMG"},
}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AMINO[16 * i + 4 * j + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}


def read_fasta(path):
    sequences, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            sequences[name] = []
        elif name:
            sequences[name].append(line.upper())
    return {name: "".join(parts) for name, parts in sequences.items()}


def read_table(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def inspection_section(text, title):
    """Rows of one '== title ==' block of the inspection report, as dicts."""
    match = re.search(r"== %s ==\n(.*?)\n\n" % re.escape(title), text, re.S)
    if not match:
        return []
    lines = match.group(1).split("\n")
    if len(lines) < 2 or lines[1] == "(no rows)":
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def main(ensembl_dir, organelles_path):
    real = read_fasta(organelles_path)
    region_seq = {name: read_fasta("%s/region_%s.fa" % (ensembl_dir, name))[name] for name in EXPECTED}
    report = open(ensembl_dir + "/organelle_inspection.txt").read()

    transcripts = {r["transcript_id"]: r for r in read_table(ensembl_dir + "/transcripts.tsv")
                   if r["seq_region"] in EXPECTED}
    exons = collections.defaultdict(list)
    for r in read_table(ensembl_dir + "/exons.tsv"):
        if r["transcript_id"] in transcripts:
            exons[r["transcript_id"]].append((int(r["rank"]), r["exon_dbid"], int(r["exon_start"]), int(r["exon_end"])))
    translations = {r["transcript_id"]: r for r in read_table(ensembl_dir + "/translations.tsv")
                    if r["transcript_id"] in transcripts}

    def protein(transcript):
        info, sequence = transcripts[transcript], region_seq[transcripts[transcript]["seq_region"]]
        pieces, offset, start, end = [], 0, None, None
        for _, dbid, exon_start, exon_end in sorted(exons[transcript]):
            piece = sequence[exon_start - 1:exon_end]
            if info["strand"] == "-1":
                piece = piece.translate(COMPLEMENT)[::-1]
            if dbid == translations[transcript]["start_exon_dbid"]:
                start = offset + int(translations[transcript]["seq_start"])
            if dbid == translations[transcript]["end_exon_dbid"]:
                end = offset + int(translations[transcript]["seq_end"])
            offset += len(piece)
            pieces.append(piece)
        cds = "".join(pieces)[start - 1:end]
        return "".join(CODON.get(cds[i:i + 3], "X") for i in range(0, len(cds) - len(cds) % 3, 3))

    attribs = inspection_section(report, "attributes of Mt, Pt and, for comparison, chromosome 1")
    pieces = inspection_section(report, "every piece of Mt and Pt, by name")
    leftovers = inspection_section(report, "regions with a leftover temporary name")
    problems = []

    for region, want in EXPECTED.items():
        sequence = region_seq[region]
        print("== region %s (should be the %s) ==" % (region, want["what"]))

        matches = [name for name, ref in real.items() if ref == sequence]
        ok = matches == [want["reference"]]
        other = {"ChrC": "chloroplast AP000423.1", "ChrM": "mitochondrion BK010421.1"}
        found = other[matches[0]] if matches else "neither organelle genome"
        print("  sequence      %s  %d bp, identical to %s" % ("OK  " if ok else "FAIL", len(sequence), found))
        if not ok:
            problems.append("%s sequence is %s" % (region, found))

        genes = {r["gene_id"]: r for r in transcripts.values() if r["seq_region"] == region}
        prefixes = collections.Counter(gene[:4] for gene in genes)
        ok = set(prefixes) == {want["prefix"]}
        print("  genes         %s  %d genes: %s" % ("OK  " if ok else "FAIL", len(genes),
                                                    ", ".join("%d %s" % (n, p) for p, n in prefixes.most_common())))
        if not ok:
            problems.append("%s carries %s genes" % (region, "/".join(prefixes)))
        print("                    biotypes: " + ", ".join(
            "%d %s" % (n, b) for b, n in collections.Counter(r["gene_biotype"] for r in genes.values()).most_common()))

        coding = sorted(t for t in translations if transcripts[t]["seq_region"] == region)
        broken = []
        for transcript in coding:
            aa = protein(transcript)
            body = aa[:-1] if aa.endswith("*") else aa
            if "*" in body:
                broken.append("%s (%d stops in %d aa)" % (transcript, body.count("*"), len(body)))
        print("  translations  %s  %d of %d proteins have no internal stop" %
              ("OK  " if not broken else "NOTE", len(coding) - len(broken), len(coding)))
        for item in broken:
            print("                    internal stops: " + item)

        ranks = [a["value"] for a in attribs if a["region"] == region and a["code"] == "karyotype_rank"]
        ok = ranks == [region] or all(v.isdigit() for v in ranks)
        print("  karyotype_rank %s value %s" % ("OK  " if ok else "FAIL", ", ".join(ranks) if ranks else "(none)"))
        if not ok:
            problems.append("%s has karyotype_rank %s" % (region, ", ".join(ranks)))

        names = [p["piece"] for p in pieces if p["region"] == region]
        wrong = [n for n in names if not n.startswith(region + ":")]
        print("  piece names   %s  %s" % ("OK  " if not wrong else "FAIL", ", ".join(names) if names else "(none listed)"))
        if wrong:
            problems.append("%s is built from pieces named %s" % (region, ", ".join(wrong)))

        others = [a for a in attribs if a["region"] == region and a["code"] not in ("karyotype_rank", "toplevel")]
        print("  other attributes   %s" % (", ".join("%s=%s" % (a["code"], a["value"]) for a in others) or "none set"))
        print()

    if leftovers:
        problems.append("regions still carry a temporary name: " + ", ".join(r["name"] for r in leftovers))
    print("== verdict ==")
    if problems:
        for problem in problems:
            print("  PROBLEM: " + problem)
    else:
        print("  Both organelle regions are labelled correctly.")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
