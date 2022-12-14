#!/usr/bin/env python3
"""
Compare two SAM/BAM files and optionally determine whether the differences
make the results better or worse compared to a third file representing the
"truth".

Queries in all input files must be in the same order (with identical QNAME).
Record identity is based on reference name and reference start position
(RNAME, POS) only. Differences in the other SAM fields (flags, MAPQ, CIGAR,
RNEXT, PNEXT, TLEN) are ignored.
"""
import sys
from argparse import ArgumentParser
from contextlib import ExitStack
from itertools import repeat, zip_longest

from pysam import AlignmentFile, AlignedSegment


def main():
    parser = ArgumentParser()
    parser.add_argument("before")
    parser.add_argument("after")
    parser.add_argument("--truth")
    parser.add_argument(
        "--limit", metavar="N", type=int, default=10,
        help="Report details for at most N changed records"
    )
    args = parser.parse_args()
    has_truth = bool(args.truth)
    limit = args.limit
    with ExitStack() as stack:
        before = stack.enter_context(AlignmentFile(args.before))
        after = stack.enter_context(AlignmentFile(args.after))
        if has_truth:
            truth = stack.enter_context(AlignmentFile(args.truth))
        else:
            truth = repeat(None)
        single_total = 0
        unmapped_same = 0
        became_unmapped = 0
        became_mapped = 0
        identical = 0
        multimapper_same = 0
        multimapper_better = 0
        multimapper_worse = 0
        changed = 0
        reported = 0

        # The following three are only updated if truth is available
        same = 0
        better = 0
        worse = 0

        for b, a, t in zip_longest(before, after, truth):
            if (b is None and a is None) and not has_truth:
                break
            if b is None or a is None:
                sys.exit("Input files have different lengths")
            assert b.query_name[:-2] == a.query_name[:-2]
            if has_truth:
                assert a.query_name[:-2] == t.query_name
            single_total += 1

            if b.is_unmapped and a.is_unmapped:
                unmapped_same += 1
                continue

            if a.is_unmapped:
                became_unmapped += 1
                continue

            if b.is_unmapped:
                became_mapped += 1
                continue

            assert not a.is_unmapped and not b.is_unmapped

            b_tup = (b.reference_name, b.reference_start)
            a_tup = (a.reference_name, a.reference_start)
            if b_tup == a_tup:
                identical += 1
                continue

            b_score = b.get_tag("AS")
            a_score = a.get_tag("AS")
            if b.mapping_quality == 0 and a.mapping_quality == 0 and b.is_proper_pair == a.is_proper_pair:
                if b_score == a_score:
                    multimapper_same += 1
                    continue

                if b_score < a_score:
                    multimapper_better += 1
                    continue

                if b_score > a_score:
                    multimapper_worse += 1
                    if reported < limit:
                        print_comparison(b, a)
                    reported += 1
                    continue

            if has_truth:
                t_tup = (t.reference_name, t.reference_start)
                if a_tup != t_tup and b_tup != t_tup:
                    same += 1
                    continue

                if a_tup == t_tup:
                    better += 1
                    continue
                assert b_tup == t_tup
                worse += 1
            changed += 1
            if reported < limit:
                print_comparison(b, a)
            reported += 1

    if changed > limit:
        print(
            f"Reporting limit reached, not showing {reported - limit} "
            "additional changed records."
        )
    print()

    def stat(description, value, should_be_zero: bool = True):
        ok = not should_be_zero or value == 0
        print(f"{'*' if not ok else ' '} {value:>9} reads {description}")

    print("Before/after comparisons")
    print()
    stat("were unmapped before and after", unmapped_same, False)
    stat("became mapped", became_mapped)
    stat("became unmapped", became_unmapped)
    stat("were mapped to same locus before and after", identical, False)
    stat("were multimapper before and after, same alignment score (AS)", multimapper_same)
    stat("were multimapper before and after, better alignment score (AS)", multimapper_better)
    stat("were multimapper before and after, worse alignment score (AS)", multimapper_worse)
    if has_truth:
        stat("were incorrect before and after (relative to truth)", same, False)
        stat("became correct (relative to truth)", better)
        stat("became incorrect (relative to truth)", worse)
    else:
        stat("changed in another way", changed)

    print(f"  {single_total:>9} total reads")

    if unmapped_same + identical < single_total:
        sys.exit(1)


def print_comparison(b: AlignedSegment, a: AlignedSegment):
    assert b.query_name[:-2] == a.query_name[:-2]
    print(b.query_name)

    def compare(name, before, after):
        if before == after:
            s = str(before)
        else:
            s = f"{before} -> {after}"
        print("  ", name.rjust(12) + ":", s)

    compare("MAPQ", b.mapping_quality, a.mapping_quality)
    compare("score", b.get_tag("AS"), a.get_tag("AS"))
    compare("NM", b.get_tag("NM"), a.get_tag("NM"))
    compare("TLEN", b.template_length, a.template_length)
    compare("proper pair", yesno(b.is_proper_pair), yesno(a.is_proper_pair))
    compare("CIGAR", b.cigarstring, a.cigarstring)
    compare("ref", b.reference_name, a.reference_name)
    compare("pos", b.reference_start, a.reference_start)
    print()


def yesno(b: bool) -> str:
    return "yes" if b else "no"


if __name__ == "__main__":
    main()
