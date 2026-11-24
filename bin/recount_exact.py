#!/usr/bin/env python3
"""Recount every TRaCE election from TRaCE's own output, with exact arithmetic.

Why: in TRaCE.pl each of the n samples that vote on a gene carries weight 6/n,
and the three length voters (domain 3, CDS 2, exon length 1) total 6. The two
sides are meant to balance, but 6/n is usually not exactly representable as a
floating-point number, so n copies of it can sum to slightly more or slightly
less than 6. When all voting samples prefer one isoform and all three length
criteria prefer another, that rounding error decides the election, in a
direction that depends only on n.

This script repeats the election logic of TRaCE.pl (sub rank_transcripts,
lines 51-165) twice for every multi-transcript gene:
  exact  vote totals as fractions, so true ties stay ties
  float  vote totals as double-precision floats, as in TRaCE.pl
and compares both with the ranking TRaCE.pl printed.

  exact == TRaCE             the outcome did not depend on rounding
  exact != TRaCE == float    the outcome was decided by rounding of 6/n sums
  exact != TRaCE != float    another order-dependent effect (the order in which
                             Perl happened to add votes, or AED values that
                             differ only beyond the 15 digits TRaCE prints)

Usage: recount_exact.py TRaCE_output.txt out_prefix [max_aed]
Writes out_prefix.tsv (one line per gene whose ranking differs) and prints a
summary to standard output.
"""
import collections
import sys
from fractions import Fraction

LENGTH_WEIGHT = {"domain": 3, "cds": 2, "len": 1}


def dense_ranks(items, key, reverse=False):
    """Rank items 1, 2, 3 ... by key; equal keys share a rank."""
    ranks, rank, previous = {}, 0, None
    for item in sorted(items, key=key, reverse=reverse):
        value = key(item)
        if rank == 0 or value != previous:
            rank += 1
            previous = value
        ranks[item] = rank
    return ranks


def elect(transcripts, sample_names, max_aed, exact):
    """Return (ranking, number of voting samples) for one gene.

    transcripts maps transcript ID -> {"len", "cds", "domain", "aed": [..]}.
    Mirrors TRaCE.pl: voters keep their original ranks as winners are removed,
    each voter's best remaining rank r is counted at preference level r - 1,
    and a tie that survives every level goes to the alphabetically first ID.
    """
    candidates = sorted(transcripts)

    # Sample voters: only transcripts at or under the AED cutoff are ranked.
    ballots = {}
    for i, sample in enumerate(sample_names):
        scored = [t for t in candidates if transcripts[t]["aed"][i] <= max_aed]
        if scored:
            ballots[sample] = dense_ranks(scored, key=lambda t: transcripts[t]["aed"][i])
    voting = len(ballots)
    if voting:
        sample_weight = Fraction(6, voting) if exact else 6 / voting
    else:
        sample_weight = None
    weight = {sample: sample_weight for sample in ballots}

    # Length voters rank every transcript, longest first.
    for metric, metric_weight in LENGTH_WEIGHT.items():
        voter = metric + "_length"
        ballots[voter] = dense_ranks(candidates, key=lambda t: transcripts[t][metric], reverse=True)
        weight[voter] = Fraction(metric_weight) if exact else float(metric_weight)

    winners = []
    remaining = list(candidates)
    while len(remaining) > 1:
        votes = collections.defaultdict(dict)       # preference level -> transcript -> total
        for voter in sorted(ballots):
            ranks = ballots[voter]
            if not ranks:
                continue
            level = min(ranks.values()) - 1
            previous = None
            for transcript in sorted(ranks, key=lambda t: ranks[t]):
                if previous is not None and ranks[transcript] > previous:
                    level += 1
                previous = ranks[transcript]
                votes[level][transcript] = votes[level].get(transcript, 0) + weight[voter]

        running = set(remaining)
        max_rounds = len(votes)                      # as in TRaCE.pl: the number of levels in use
        level = 0
        while len(running) > 1 and level < max_rounds:
            tallies = {t: votes[level][t] for t in running if level in votes and t in votes[level]}
            if tallies:
                best = max(tallies.values())
                running = {t for t, total in tallies.items() if total == best}
            level += 1
        winner = min(running)
        winners.append(winner)
        for ranks in ballots.values():
            ranks.pop(winner, None)
        remaining = sorted(ballots["len_length"])
    winners.extend(remaining)
    return winners, voting


def main(trace_path, out_prefix, max_aed=0.5):
    genes = collections.OrderedDict()
    printed_rank = {}
    n_samples = None
    for line in open(trace_path):
        cols = line.rstrip("\n").split("\t")
        gene, transcript, rank = cols[0], cols[1], int(cols[2])
        aed = [float(x) for x in cols[6:]]
        if n_samples is None:
            n_samples = len(aed)
        assert len(aed) == n_samples, "ragged line for %s" % transcript
        genes.setdefault(gene, {})[transcript] = {
            "len": int(cols[3]), "cds": int(cols[4]), "domain": int(cols[5]), "aed": aed}
        printed_rank[transcript] = rank
    sample_names = ["sample%02d" % i for i in range(n_samples)]

    tally = collections.Counter()
    by_voters = collections.defaultdict(collections.Counter)
    with open(out_prefix + ".tsv", "w") as out:
        out.write("\t".join(["gene_id", "transcripts", "voting_samples", "cause", "rank1_changed",
                             "TRaCE_ranking", "exact_ranking", "float_ranking"]) + "\n")
        for gene, transcripts in genes.items():
            if len(transcripts) < 2:
                continue
            tally["genes"] += 1
            printed = sorted(transcripts, key=lambda t: printed_rank[t])
            exact, voting = elect(transcripts, sample_names, max_aed, exact=True)
            by_voters[voting]["genes"] += 1
            if exact == printed:
                continue
            floated, _ = elect(transcripts, sample_names, max_aed, exact=False)
            cause = "weight_rounding" if floated == printed else "other_order_effect"
            top_changed = exact[0] != printed[0]
            tally["ranking_differs"] += 1
            tally["ranking_differs_" + cause] += 1
            if top_changed:
                tally["rank1_differs"] += 1
                tally["rank1_differs_" + cause] += 1
                by_voters[voting]["rank1_differs"] += 1
            out.write("\t".join([gene, str(len(transcripts)), str(voting), cause,
                                 "yes" if top_changed else "no",
                                 ",".join(printed), ",".join(exact), ",".join(floated)]) + "\n")

    print("multi-transcript genes recounted: %d (%d samples)" % (tally["genes"], n_samples))
    print("full ranking differs from the exact recount: %d" % tally["ranking_differs"])
    print("    explained by rounding of the 6/n vote sums: %d" % tally["ranking_differs_weight_rounding"])
    print("    other order-dependent effects:              %d" % tally["ranking_differs_other_order_effect"])
    print("rank 1 (canonical) differs from the exact recount: %d" % tally["rank1_differs"])
    print("    explained by rounding of the 6/n vote sums: %d" % tally["rank1_differs_weight_rounding"])
    print("    other order-dependent effects:              %d" % tally["rank1_differs_other_order_effect"])
    print("by number of voting samples (genes / rank 1 differs):")
    for voting in sorted(by_voters):
        if by_voters[voting]["rank1_differs"]:
            print("    n=%-3d %6d / %d" % (voting, by_voters[voting]["genes"], by_voters[voting]["rank1_differs"]))


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) == 4 else 0.5)
