//! Fast, hardcoded sanity tests for the AVX2 aligner.
//!
//! These are deliberately *not* an oracle suite: there is no scalar reference implementation
//! and nothing sweeps randomized inputs. They are hand-written cases with hand-checked expected
//! values, meant to catch an obvious break quickly.
//!
//! Two invariants are worth more than the individual cases, and every test that produces an
//! alignment routes through them via [`check`]:
//!
//! 1. **The CIGAR re-scores to the reported score.** A wrong traceback under a right score is
//!    this kernel's characteristic failure (see `TB_FROM_DIAG` in `avx2_u8.rs`), and comparing
//!    scores alone cannot see it.
//! 2. **The CIGAR's consumption matches the reported spans**, so the coordinates and the
//!    operations cannot disagree about how much sequence was used.

use super::{AlignmentResult, Scores, SimdAligner};

const SCORES: Scores = Scores {
    match_: 2,
    mismatch: 8,
    gap_open: 12,
    gap_extend: 1,
    end_bonus: 10,
};

fn aligner() -> SimdAligner {
    SimdAligner::new(SCORES)
}

/// The kernel's scoring alphabet: A/C/G/T/U, case-insensitive. Anything else is unknown and
/// never matches, *including another unknown byte*.
fn scores_as_match(q: u8, r: u8) -> bool {
    fn code(b: u8) -> Option<u8> {
        match b.to_ascii_uppercase() {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' | b'U' => Some(3),
            _ => None,
        }
    }
    matches!((code(q), code(r)), (Some(a), Some(b)) if a == b)
}

/// Parse a CIGAR string back into `(len, op)` pairs.
fn parse(cigar: &str) -> Vec<(usize, char)> {
    let mut out = vec![];
    let mut len = 0usize;
    for ch in cigar.chars() {
        if let Some(d) = ch.to_digit(10) {
            len = len * 10 + d as usize;
        } else {
            out.push((len, ch));
            len = 0;
        }
    }
    assert_eq!(
        len, 0,
        "CIGAR ended with a length and no operation: {cigar}"
    );
    out
}

/// Re-score `result`'s CIGAR against the sequences it claims to align, and assert it agrees
/// with the reported score and spans. `end_bonus` is added by the caller for the modes that
/// earn it, since the CIGAR itself cannot show whether the query was spanned.
fn check(result: &AlignmentResult, query: &[u8], reference: &[u8], end_bonus: i32) {
    let (mut qi, mut ri) = (result.query_start, result.ref_start);
    let mut score = 0i32;

    for (len, op) in parse(&result.cigar.to_string()) {
        match op {
            '=' | 'X' | 'M' => {
                for _ in 0..len {
                    let (q, r) = (query[qi], reference[ri]);
                    score += if scores_as_match(q, r) {
                        SCORES.match_ as i32
                    } else {
                        -(SCORES.mismatch as i32)
                    };
                    // `=`/`X` are literal byte identity, decoupled from how the cell scored.
                    if op != 'M' {
                        let identical = q.to_ascii_uppercase() == r.to_ascii_uppercase();
                        assert_eq!(
                            identical,
                            op == '=',
                            "op {op} at query[{qi}]={} reference[{ri}]={} disagrees with byte identity",
                            q as char,
                            r as char
                        );
                    }
                    qi += 1;
                    ri += 1;
                }
            }
            'I' => {
                score -= SCORES.gap_open as i32 + (len as i32 - 1) * SCORES.gap_extend as i32;
                qi += len;
            }
            'D' => {
                score -= SCORES.gap_open as i32 + (len as i32 - 1) * SCORES.gap_extend as i32;
                ri += len;
            }
            other => panic!("unexpected CIGAR operation {other}"),
        }
    }

    assert_eq!(
        qi, result.query_end,
        "CIGAR consumed query up to {qi} but result says {}",
        result.query_end
    );
    assert_eq!(
        ri, result.ref_end,
        "CIGAR consumed reference up to {ri} but result says {}",
        result.ref_end
    );
    assert_eq!(
        score + end_bonus,
        result.score,
        "CIGAR {} re-scores to {} but result says {}",
        result.cigar,
        score + end_bonus,
        result.score
    );
}

/// Assert a global alignment's score and CIGAR, and re-score it.
fn global_case(query: &[u8], reference: &[u8], score: i32, cigar: &str) {
    let result = aligner().global_alignment(query, reference, None);
    check(&result, query, reference, 0);
    assert_eq!(result.score, score);
    assert_eq!(result.cigar.to_string(), cigar);
    assert_eq!((result.query_start, result.query_end), (0, query.len()));
    assert_eq!((result.ref_start, result.ref_end), (0, reference.len()));
}

// ---------------------------------------------------------------------------------
// Global.
// ---------------------------------------------------------------------------------

#[test]
fn global_perfect_match() {
    global_case(b"AAATTT", b"AAATTT", 6 * 2, "6=");
}

#[test]
fn global_complete_mismatch() {
    global_case(b"AAA", b"TTT", 3 * -8, "3X");
}

#[test]
fn global_single_mismatch() {
    global_case(b"AAATAA", b"AAAAAA", 5 * 2 - 8, "3=1X2=");
}

#[test]
fn global_gap_in_query() {
    global_case(b"AAATTT", b"AAAATTTT", 6 * 2 - 12 - 1, "3=2D3=");
}

#[test]
fn global_gap_in_reference() {
    global_case(b"AAAATTTT", b"AAATTT", 6 * 2 - 12 - 1, "3=2I3=");
}

// The four end-gap cases exercise the boundary cells, which the fill writes after the chunk
// loop rather than before - a hazard called out in the kernel's module docs.

#[test]
fn global_gap_at_query_start() {
    global_case(b"AAATTT", b"TTAAATTT", 6 * 2 - 12 - 1, "2D6=");
}

#[test]
fn global_gap_at_query_end() {
    global_case(b"AAATTT", b"AAATTTAA", 6 * 2 - 12 - 1, "6=2D");
}

#[test]
fn global_gap_at_reference_start() {
    global_case(b"TTAAATTT", b"AAATTT", 6 * 2 - 12 - 1, "2I6=");
}

#[test]
fn global_gap_at_reference_end() {
    global_case(b"AAATTTAA", b"AAATTT", 6 * 2 - 12 - 1, "6=2I");
}

#[test]
fn global_multiple_mismatches() {
    global_case(
        b"ATATATATAT",
        b"AAAAAAAAAA",
        5 * 2 - 5 * 8,
        "1=1X1=1X1=1X1=1X1=1X",
    );
}

#[test]
fn global_complex_alignment() {
    global_case(
        b"AAACTTAAACCTT",
        b"AAAATTGAAATT",
        10 * 2 - 8 - 2 * 12 - 1,
        "3=1X2=1D3=2I2=",
    );
}

#[test]
fn global_empty_query_is_all_deletion() {
    let result = aligner().global_alignment(b"", b"ACGT", None);
    check(&result, b"", b"ACGT", 0);
    assert_eq!(result.score, -(12 + 3));
    assert_eq!(result.cigar.to_string(), "4D");
}

#[test]
fn global_both_empty() {
    let result = aligner().global_alignment(b"", b"", None);
    assert_eq!(result.score, 0);
    assert!(result.cigar.is_empty());
}

// ---------------------------------------------------------------------------------
// The local modes the piecewise extension path actually calls.
// ---------------------------------------------------------------------------------

#[test]
fn local_end_stops_before_a_bad_tail() {
    // The reference matches for 6 bases then diverges completely. Extending into the tail
    // costs more than it earns, so the alignment stops at the divergence.
    let (query, reference) = (b"AAATTTGGGGGG", b"AAATTTCCCCCC");
    let result = aligner().local_end_alignment(query, reference, None);
    check(&result, query, reference, 0);
    assert_eq!(result.score, 6 * 2);
    assert_eq!(result.cigar.to_string(), "6=");
    assert_eq!((result.query_start, result.query_end), (0, 6));
    assert_eq!((result.ref_start, result.ref_end), (0, 6));
}

#[test]
fn local_end_spans_the_query_for_the_end_bonus() {
    // A single trailing mismatch is worth crossing: -8 for the mismatch against +10 of
    // end_bonus for reaching the end of the query.
    let (query, reference) = (b"AAATTTG", b"AAATTTC");
    let result = aligner().local_end_alignment(query, reference, None);
    check(&result, query, reference, SCORES.end_bonus as i32);
    assert_eq!(result.score, 6 * 2 - 8 + 10);
    assert_eq!(result.cigar.to_string(), "6=1X");
    assert_eq!(result.query_end, 7);
}

#[test]
fn local_end_never_scores_below_zero() {
    // Nothing is worth aligning, so the empty alignment wins.
    let result = aligner().local_end_alignment(b"GGGG", b"CCCC", None);
    assert_eq!(result.score, 0);
    assert!(result.cigar.is_empty());
    assert_eq!((result.query_end, result.ref_end), (0, 0));
}

#[test]
fn local_start_is_the_mirror_of_local_end() {
    // Same shape as `local_end_stops_before_a_bad_tail`, reversed: the junk is a leading
    // prefix and the alignment begins after it.
    let (query, reference) = (b"GGGGGGAAATTT", b"CCCCCCAAATTT");
    let result = aligner().local_start_alignment(query, reference, None);
    check(&result, query, reference, 0);
    assert_eq!(result.score, 6 * 2);
    assert_eq!(result.cigar.to_string(), "6=");
    assert_eq!((result.query_start, result.query_end), (6, 12));
    assert_eq!((result.ref_start, result.ref_end), (6, 12));
}

#[test]
fn local_reference_end_spans_the_whole_query() {
    // The query is spanned in full; only the reference's trailing tail is free.
    let (query, reference) = (b"AAATTT", b"AAATTTGGGGGGGG");
    let result = aligner().local_reference_end_alignment(query, reference, None);
    check(&result, query, reference, 0);
    assert_eq!(result.score, 6 * 2);
    assert_eq!((result.query_start, result.query_end), (0, 6));
    assert_eq!(result.ref_end, 6);
}

#[test]
fn local_reference_start_spans_the_whole_query() {
    let (query, reference) = (b"AAATTT", b"GGGGGGGGAAATTT");
    let result = aligner().local_reference_start_alignment(query, reference, None);
    check(&result, query, reference, 0);
    assert_eq!(result.score, 6 * 2);
    assert_eq!((result.query_start, result.query_end), (0, 6));
    assert_eq!((result.ref_start, result.ref_end), (8, 14));
}

// ---------------------------------------------------------------------------------
// Banding. `piecewisealigner.rs` always passes `Some(bandwidth)` for the end extensions and
// lets the aligner decide, so at the default `--bw 1024` a short-read extension takes the exact
// path anyway and only long reads are really banded. These pin the two things that make it safe:
// a band wide enough to reach every cell reproduces the exact answer, and a band too narrow for
// the optimum still returns a consistent alignment.
// ---------------------------------------------------------------------------------

#[test]
fn a_band_wider_than_the_problem_matches_the_exact_answer() {
    let query = b"ACGTACGTAAGGTTCCAACGTTGGCCAATTGGCCAAGGTTCCAA";
    let reference = b"ACGTACGTAAGGTACCAACGTTGGCCAATTGCCCAAGGTTCCAA";
    let mut aligner = aligner();

    let exact = aligner.global_alignment(query, reference, None);
    for w in [64, 128, 1024] {
        let banded = aligner.global_alignment(query, reference, Some(w));
        assert_eq!(banded.score, exact.score, "band {w} changed the score");
        assert_eq!(
            banded.cigar.to_string(),
            exact.cigar.to_string(),
            "band {w} changed the CIGAR"
        );
    }
}

#[test]
fn a_narrow_band_still_returns_a_consistent_alignment() {
    // A band too narrow to contain the optimum must still produce a self-consistent result:
    // the score it reports is the score its own CIGAR achieves.
    let query = b"AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
    let reference = b"AAAAAAAAAAGGGGGGGGGGTTTTTTTTTT";
    let result = aligner().global_alignment(query, reference, Some(2));
    check(&result, query, reference, 0);
}

#[test]
fn a_band_is_never_worse_than_the_exact_optimum() {
    let query = b"AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
    let reference = b"AAAAAAAAAAGGGGGGGGGGTTTTTTTTTT";
    let mut aligner = aligner();
    let exact = aligner.global_alignment(query, reference, None);
    let banded = aligner.global_alignment(query, reference, Some(2));
    assert!(
        banded.score <= exact.score,
        "banded score {} beat the exact optimum {}",
        banded.score,
        exact.score
    );
}

// ---------------------------------------------------------------------------------
// Alphabet. Scoring is biological; `=`/`X` is literal byte identity. They disagree on
// purpose, and `check` asserts the byte-identity half on every case above.
// ---------------------------------------------------------------------------------

#[test]
fn u_scores_as_t_but_reports_a_mismatch_operation() {
    let result = aligner().global_alignment(b"ACGU", b"ACGT", None);
    assert_eq!(result.score, 4 * 2, "U should score as T");
    assert_eq!(
        result.cigar.to_string(),
        "3=1X",
        "U vs T is not byte-identical, so it reports X"
    );
}

#[test]
fn case_is_ignored() {
    let result = aligner().global_alignment(b"acgt", b"ACGT", None);
    assert_eq!(result.score, 4 * 2);
}

#[test]
fn unknown_bytes_never_match_even_themselves() {
    // Two identical `N`s: byte-identical, so the CIGAR says `=`, but they score as mismatches.
    let result = aligner().global_alignment(b"ANNT", b"ANNT", None);
    assert_eq!(
        result.score,
        2 * 2 - 2 * 8,
        "N must not match N; got {}",
        result.cigar
    );
    assert_eq!(result.cigar.to_string(), "4=");
}

// ---------------------------------------------------------------------------------
// Split reference: one query across two references with a single jump.
// ---------------------------------------------------------------------------------

#[test]
fn split_reference_partitions_the_query() {
    let query = b"AAAATTTTCCCCGGGG";
    let left = b"CCCCGGGG";
    let right = b"AAAATTTT";
    let result = aligner().split_reference_alignment(query, left, right, None);

    check(&result.right, query, right, 0);
    check(&result.left, query, left, 0);

    assert_eq!(result.score, result.left.score + result.right.score);
    assert_eq!(
        result.right.query_start, 0,
        "the right arm covers the query prefix"
    );
    assert_eq!(
        result.left.query_end,
        query.len(),
        "the left arm covers the query suffix"
    );
    assert_eq!(
        result.right.query_end, result.left.query_start,
        "the two arms must meet at the jump point"
    );
    assert_eq!(result.score, 16 * 2);
}

// ---------------------------------------------------------------------------------
// Scale: enough to cross the 32-lane chunking and the anti-diagonal bookkeeping, still fast.
// ---------------------------------------------------------------------------------

#[test]
fn a_few_hundred_bases_stay_consistent() {
    let unit = b"ACGTTGCAAGGCTTAC";
    let query: Vec<u8> = unit.iter().cycle().take(500).copied().collect();
    let mut reference = query.clone();
    reference[100] = b'A';
    reference[101] = b'A';
    reference.remove(300);

    let result = aligner().global_alignment(&query, &reference, None);
    check(&result, &query, &reference, 0);
    assert_eq!(result.query_end, query.len());
    assert_eq!(result.ref_end, reference.len());
    assert!(
        result.score > 400 * 2,
        "a near-identical 500-base pair should score well; got {}",
        result.score
    );
}

#[test]
fn identical_long_sequences_are_one_run_of_matches() {
    let query: Vec<u8> = b"ACGTTGCAAGGCTTAC"
        .iter()
        .cycle()
        .take(300)
        .copied()
        .collect();
    let result = aligner().global_alignment(&query, &query, None);
    check(&result, &query, &query, 0);
    assert_eq!(result.score, 300 * 2);
    assert_eq!(result.cigar.to_string(), "300=");
}
