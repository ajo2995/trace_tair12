#!/usr/bin/env python3
"""Build the sample manifest for the 30-tissue core of the Mergner 2020 atlas.

Sources, both fetched live:
  - ENA file report for PRJEB32665 (run accession, read-pair count)
  - ArrayExpress sample sheet for E-MTAB-7978 (tissue, stage, description)

"Sample 1" to "Sample 30" of the submission are the atlas's 30 tissues; the
remaining 26 samples are developmental series and are not used. Each sample
gets a short name (S01_petal ...) that becomes its StringTie label, its GTF
file name and its column in the TRaCE output, so names must sort in sample
order, contain no whitespace, and must not contain "_length" (TRaCE.pl uses
that suffix to recognise its own length voters).

Writes, under manifest/:
  samples.core30.tsv     full table with a header line
  worker_args.core30.txt "name run read_pairs" per line, for sample_worker.sh
  gtfs.core30.txt        GTF file names in TRaCE column order

Usage: build_manifest.py manifest_dir
"""
import csv
import io
import sys
import urllib.request

ENA_REPORT = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB32665"
              "&result=read_run&fields=run_accession,sample_alias,library_layout,read_count,base_count"
              "&format=tsv")
SDRF = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/978/E-MTAB-7978/Files/E-MTAB-7978.sdrf.txt"
SRA_URL = "https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"

# Sample number -> (short label, run accession expected from the 2026-11-20 review).
# The run accessions are checked against ENA so that a change upstream is noticed.
CORE = {
    1: ("petal", "ERR3333388"), 2: ("sepal", "ERR3333399"), 3: ("stamen", "ERR3333410"),
    4: ("carpel", "ERR3333421"), 5: ("silique", "ERR3333432"), 6: ("seed", "ERR3333440"),
    7: ("dry_seed", "ERR3333441"), 8: ("imbibed_seed", "ERR3333442"), 9: ("pollen", "ERR3333443"),
    10: ("stem_node", "ERR3333389"), 11: ("stem_internode", "ERR3333390"),
    12: ("cauline_leaf", "ERR3333391"), 13: ("leaf_distal", "ERR3333392"),
    14: ("leaf_proximal", "ERR3333393"), 15: ("petiole", "ERR3333394"),
    16: ("senescent_leaf", "ERR3333395"), 17: ("root", "ERR3333396"), 18: ("pedicel", "ERR3333397"),
    19: ("flower", "ERR3333398"), 20: ("silique_septum", "ERR3333400"),
    21: ("silique_valves", "ERR3333401"), 22: ("egg_cell_like_callus", "ERR3333402"),
    23: ("callus", "ERR3333403"), 24: ("cell_culture_3d", "ERR3333404"),
    25: ("cell_culture_10d", "ERR3333405"), 26: ("cotyledon", "ERR3333406"),
    27: ("hypocotyl", "ERR3333407"), 28: ("root_tip", "ERR3333408"),
    29: ("root_upper_zone", "ERR3333409"), 30: ("shoot_apex", "ERR3333411"),
}


def fetch_table(url):
    with urllib.request.urlopen(url, timeout=120) as response:
        text = response.read().decode("utf-8", errors="replace")
    return list(csv.reader(io.StringIO(text), delimiter="\t"))


def sample_number(text):
    """'E-MTAB-7978:Sample 12' or 'Sample 12' -> 12"""
    return int(text.rsplit("Sample ", 1)[1])


def main(out_dir):
    ena = fetch_table(ENA_REPORT)
    header, rows = ena[0], ena[1:]
    col = {name: i for i, name in enumerate(header)}
    runs = {}
    for row in rows:
        runs[sample_number(row[col["sample_alias"]])] = {
            "run": row[col["run_accession"]],
            "layout": row[col["library_layout"]],
            "read_pairs": int(row[col["read_count"]]),
            "bases": int(row[col["base_count"]]),
        }

    sdrf = fetch_table(SDRF)
    header, rows = sdrf[0], sdrf[1:]
    first = {}
    for i, name in enumerate(header):      # some column names repeat; the first one is wanted
        first.setdefault(name, i)
    tissue = {}
    for row in rows:                       # two rows per sample, one per FASTQ file
        tissue[sample_number(row[first["Source Name"]])] = {
            "sdrf_run": row[first["Comment[ENA_RUN]"]],
            "organism_part": row[first["Characteristics[organism part]"]],
            "stage": row[first["Characteristics[developmental stage]"]],
            "site": row[first["Characteristics[sampling site]"]],
            "description": row[first["Description"]],
            "ecotype": row[first["Characteristics[ecotype]"]],
            "genotype": row[first["Characteristics[genotype]"]],
        }

    lines, worker, gtfs = [], [], []
    for number in sorted(CORE):
        label, expected_run = CORE[number]
        name = "S%02d_%s" % (number, label)
        assert "_length" not in name and not any(c.isspace() for c in name), name
        run, meta = runs[number], tissue[number]
        assert run["run"] == expected_run == meta["sdrf_run"], (number, run["run"], expected_run, meta["sdrf_run"])
        assert run["layout"] == "PAIRED", (number, run["layout"])
        assert meta["ecotype"] == "Col-0" and meta["genotype"] == "wild type genotype", (number, meta)
        lines.append([name, run["run"], str(run["read_pairs"]), str(run["bases"]),
                      meta["organism_part"], meta["stage"], meta["site"], meta["description"],
                      SRA_URL.format(run=run["run"])])
        worker.append("%s %s %d" % (name, run["run"], run["read_pairs"]))
        gtfs.append(name + ".gtf")
    assert gtfs == sorted(gtfs), "sample names must sort in sample order"

    with open(out_dir + "/samples.core30.tsv", "w") as out:
        out.write("#" + "\t".join(["sample", "run", "read_pairs", "bases", "organism_part", "stage",
                                   "sampling_site", "description", "sra_url"]) + "\n")
        for line in lines:
            out.write("\t".join(line) + "\n")
    with open(out_dir + "/worker_args.core30.txt", "w") as out:
        out.write("\n".join(worker) + "\n")
    with open(out_dir + "/gtfs.core30.txt", "w") as out:
        out.write("\n".join(gtfs) + "\n")

    pairs = sum(int(line[2]) for line in lines)
    print("%d samples, %.1f million read pairs, %.1f Gb" %
          (len(lines), pairs / 1e6, sum(int(line[3]) for line in lines) / 1e9), file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
