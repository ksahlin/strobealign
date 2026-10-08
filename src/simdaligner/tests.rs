//! Fast, hardcoded sanity tests for the public aligner API, against whichever backend
//! `select_backend` picks on the host.
//!
//! These are hand-written cases with hand-checked expected values, meant to catch an obvious
//! break quickly. The randomized sweeps against a textbook DP are in [`oracle`] at the bottom of
//! this file, and the differential testing against the portable kernel lives in `backend/tests.rs`.
//!
//! Two invariants are worth more than the individual cases, and every test that produces an
//! alignment routes through them via [`check`]:
//!
//! 1. **The CIGAR re-scores to the reported score.** A wrong traceback under a right score is
//!    this kernel's characteristic failure (see `TB_FROM_DIAG` in `kernel.rs`), and comparing
//!    scores alone cannot see it.
//! 2. **The CIGAR's consumption matches the reported spans**, so the coordinates and the
//!    operations cannot disagree about how much sequence was used.

use super::{AlignmentResult, Scores, SimdAligner};

const SCORES: Scores = Scores {
    match_: 2,
    mismatch: 8,
    gap_open: [12, 12],
    gap_extend: [1, 1],
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
                        let identical = q.eq_ignore_ascii_case(&r);
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
                score -= SCORES.gap_cost(len);
                qi += len;
            }
            'D' => {
                score -= SCORES.gap_cost(len);
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
// loop rather than before, a hazard called out in the kernel's module docs.

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

/// The two-piece gap cost, and the alignment it exists to prevent: a deleted block with a few
/// bases left matching in the middle of it, which one affine level is cheaper splitting in two
/// than spanning.
mod two_piece {
    use super::super::{Scores, SimdAligner};

    const AFFINE: Scores = Scores {
        match_: 2,
        mismatch: 8,
        gap_open: [12, 12],
        gap_extend: [1, 1],
        end_bonus: 10,
    };

    const TWO_PIECE: Scores = Scores {
        match_: 2,
        mismatch: 8,
        gap_open: [12, 36],
        gap_extend: [2, 1],
        end_bonus: 10,
    };

    /// A reference whose middle 100 bases are missing from the query, except for a run of
    /// `island` bases taken from the centre of the deleted block.
    fn deletion_with_an_island(island: usize) -> (Vec<u8>, Vec<u8>) {
        // A fixed xorshift stream, so the reference is the same bytes every run.
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let mut reference = Vec::new();
        for _ in 0..200 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            reference.push(b"ACGT"[(state >> 33) as usize % 4]);
        }
        let (flank, del) = (50usize, 100usize);
        let mid = flank + del / 2;
        let mut query = Vec::new();
        query.extend_from_slice(&reference[..flank]);
        query.extend_from_slice(&reference[mid..mid + island]);
        query.extend_from_slice(&reference[flank + del..]);
        (query, reference)
    }

    fn deletion_runs(cigar: &str) -> usize {
        cigar.matches('D').count()
    }

    #[test]
    fn one_affine_level_splits_a_deletion_around_a_matching_island() {
        let (query, reference) = deletion_with_an_island(5);
        let cigar = SimdAligner::new(AFFINE)
            .global_alignment(&query, &reference, None)
            .cigar
            .to_string();
        assert_eq!(cigar, "50=50D5=45D50=");
        assert_eq!(deletion_runs(&cigar), 2);
    }

    #[test]
    fn a_second_level_keeps_the_deletion_whole() {
        let (query, reference) = deletion_with_an_island(5);
        let cigar = SimdAligner::new(TWO_PIECE)
            .global_alignment(&query, &reference, None)
            .cigar
            .to_string();
        // One deletion, not two. Which island bases land as an insertion is down to the
        // sequence, so only the run count is asserted.
        assert_eq!(
            deletion_runs(&cigar),
            1,
            "the deletion should stay whole, got {cigar}"
        );
    }

    /// The island has to be small for the merge to be the cheaper alignment: past the point
    /// where the matches outweigh one long opening, splitting is genuinely the better
    /// alignment and the scheme says so.
    #[test]
    fn a_long_enough_island_still_splits() {
        let (query, reference) = deletion_with_an_island(20);
        let cigar = SimdAligner::new(TWO_PIECE)
            .global_alignment(&query, &reference, None)
            .cigar
            .to_string();
        assert_eq!(
            deletion_runs(&cigar),
            2,
            "a 20-base island should still split the deletion, got {cigar}"
        );
    }

    /// Short indels are priced off the short level: below the crossover the long one is never
    /// the cheaper.
    #[test]
    fn short_gaps_are_priced_off_the_short_level() {
        for k in 1..=10 {
            assert_eq!(TWO_PIECE.gap_cost(k), 12 + (k as i32 - 1) * 2, "gap of {k}");
        }
        // 12 + 2*(k-1) meets 36 + (k-1) at k = 25.
        assert_eq!(TWO_PIECE.gap_cost(25), 36 + 24);
        assert_eq!(TWO_PIECE.gap_cost(100), 36 + 99);
    }
}

/// A textbook DP, and every mode checked against it.
///
/// The differential suite in `backend/tests.rs` says the backends agree with each other, not
/// that they are right. [`Dp`] is the two-piece recurrence written out cell by cell with no
/// offsets, no anti-diagonals and no `u8`, so it shares nothing with the kernel but the
/// problem. Every CIGAR is re-scored too, which is what a score comparison cannot see.
mod oracle {
    use crate::simdaligner::{AlignmentResult, Scores, SimdAligner};

    const NEG: i64 = i64::MIN / 4;

    /// The kernel's alphabet rule, spelled out: A/C/G/T/U case-insensitively, and an unknown byte
    /// never matches anything, including another unknown byte.
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

    /// The full score matrix of a two-piece affine alignment, filled the obvious way.
    struct Dp {
        /// `s[i][j] == S[i][j]`, `(qlen + 1) * (rlen + 1)`.
        s: Vec<Vec<i64>>,
    }

    impl Dp {
        fn new(query: &[u8], refseq: &[u8], sc: Scores) -> Dp {
            Dp::banded(query, refseq, sc, usize::MAX)
        }

        /// The same recurrence confined to `|i - j| <= w`: out-of-band cells are `NEG` and so are
        /// never anyone's predecessor, which is what the kernel's seals mean.
        fn banded(query: &[u8], refseq: &[u8], sc: Scores, w: usize) -> Dp {
            let (n, m) = (query.len(), refseq.len());
            let (q1, e1) = (sc.gap_open[0] as i64, sc.gap_extend[0] as i64);
            let (q2, e2) = (sc.gap_open[1] as i64, sc.gap_extend[1] as i64);

            let mut s = vec![vec![NEG; m + 1]; n + 1];
            // E/F are the vertical (query-consuming) and horizontal gap states, one pair per level.
            let mut e = vec![vec![NEG; m + 1]; n + 1];
            let mut f = vec![vec![NEG; m + 1]; n + 1];
            let mut e2m = vec![vec![NEG; m + 1]; n + 1];
            let mut f2m = vec![vec![NEG; m + 1]; n + 1];

            let live = |i: usize, j: usize| i.abs_diff(j) <= w;

            s[0][0] = 0;
            for (j, cell) in s[0].iter_mut().enumerate().take(m + 1).skip(1) {
                if live(0, j) {
                    *cell = -sc.gap_cost(j) as i64;
                }
            }
            for (i, row) in s.iter_mut().enumerate().take(n + 1).skip(1) {
                if live(i, 0) {
                    row[0] = -sc.gap_cost(i) as i64;
                }
            }
            for i in 1..=n {
                for j in 1..=m {
                    if !live(i, j) {
                        continue;
                    }
                    e[i][j] = (s[i - 1][j] - q1).max(e[i - 1][j] - e1);
                    e2m[i][j] = (s[i - 1][j] - q2).max(e2m[i - 1][j] - e2);
                    f[i][j] = (s[i][j - 1] - q1).max(f[i][j - 1] - e1);
                    f2m[i][j] = (s[i][j - 1] - q2).max(f2m[i][j - 1] - e2);
                    let diag = s[i - 1][j - 1]
                        + if scores_as_match(query[i - 1], refseq[j - 1]) {
                            sc.match_ as i64
                        } else {
                            -(sc.mismatch as i64)
                        };
                    s[i][j] = diag.max(e[i][j]).max(f[i][j]).max(e2m[i][j]).max(f2m[i][j]);
                }
            }
            Dp { s }
        }
    }

    /// Re-score `result`'s CIGAR against the sequences it claims to align, and check it agrees with
    /// the reported score and spans.
    fn rescore(
        result: &AlignmentResult,
        query: &[u8],
        refseq: &[u8],
        sc: Scores,
        bonus: i64,
    ) -> i64 {
        let (mut qi, mut ri) = (result.query_start, result.ref_start);
        let mut score = 0i64;
        let text = result.cigar.to_string();
        let mut len = 0usize;
        for ch in text.chars() {
            if let Some(d) = ch.to_digit(10) {
                len = len * 10 + d as usize;
                continue;
            }
            match ch {
                '=' | 'X' | 'M' => {
                    for _ in 0..len {
                        score += if scores_as_match(query[qi], refseq[ri]) {
                            sc.match_ as i64
                        } else {
                            -(sc.mismatch as i64)
                        };
                        qi += 1;
                        ri += 1;
                    }
                }
                'I' => {
                    score -= sc.gap_cost(len) as i64;
                    qi += len;
                }
                'D' => {
                    score -= sc.gap_cost(len) as i64;
                    ri += len;
                }
                other => panic!("unexpected CIGAR operation {other}"),
            }
            len = 0;
        }
        assert_eq!(qi, result.query_end, "CIGAR vs query_end in {text}");
        assert_eq!(ri, result.ref_end, "CIGAR vs ref_end in {text}");
        assert_eq!(
            score + bonus,
            result.score as i64,
            "CIGAR {text} re-scores to {score}+{bonus}, result says {}",
            result.score
        );
        score + bonus
    }

    fn rev(v: &[u8]) -> Vec<u8> {
        v.iter().rev().copied().collect()
    }

    /// Every mode's objective, computed from the plain DP, against what the kernel returned.
    fn check_all(query: &[u8], refseq: &[u8], refseq2: &[u8], sc: Scores) {
        let mut al = SimdAligner::new(sc);
        let (n, m) = (query.len(), refseq.len());
        let dp = Dp::new(query, refseq, sc);
        let bonus = sc.end_bonus as i64;
        let ctx = || {
            format!(
                "q={:?} r={:?} sc={sc:?}",
                String::from_utf8_lossy(query),
                String::from_utf8_lossy(refseq)
            )
        };

        let r = al.global_alignment(query, refseq, None);
        assert_eq!(r.score as i64, dp.s[n][m], "global: {}", ctx());
        rescore(&r, query, refseq, sc, 0);

        let r = al.local_reference_end_alignment(query, refseq, None);
        let want = (0..=m).map(|j| dp.s[n][j]).max().unwrap();
        assert_eq!(r.score as i64, want, "local_reference_end: {}", ctx());
        rescore(&r, query, refseq, sc, 0);

        // The mirror modes are the same objective on reversed inputs, so the oracle runs reversed
        // too rather than growing a second recurrence.
        let (rq, rr) = (rev(query), rev(refseq));
        let dpr = Dp::new(&rq, &rr, sc);
        let r = al.local_reference_start_alignment(query, refseq, None);
        let want = (0..=m).map(|j| dpr.s[n][j]).max().unwrap();
        assert_eq!(r.score as i64, want, "local_reference_start: {}", ctx());
        rescore(&r, query, refseq, sc, 0);

        let r = al.local_end_alignment(query, refseq, None);
        let want = (0..=n)
            .flat_map(|i| (0..=m).map(move |j| (i, j)))
            .map(|(i, j)| dp.s[i][j] + if i == n { bonus } else { 0 })
            .max()
            .unwrap()
            .max(0);
        assert_eq!(r.score as i64, want, "local_end: {}", ctx());
        rescore(
            &r,
            query,
            refseq,
            sc,
            if r.query_end == n { bonus } else { 0 },
        );

        let r = al.local_start_alignment(query, refseq, None);
        let want = (0..=n)
            .flat_map(|i| (0..=m).map(move |j| (i, j)))
            .map(|(i, j)| dpr.s[i][j] + if i == n { bonus } else { 0 })
            .max()
            .unwrap()
            .max(0);
        assert_eq!(r.score as i64, want, "local_start: {}", ctx());
        rescore(
            &r,
            query,
            refseq,
            sc,
            if r.query_start == 0 { bonus } else { 0 },
        );

        // split_reference: the right arm runs forward against `refseq2`, the left arm reversed
        // against `refseq`, and the jump point maximises the sum.
        let right = Dp::new(query, refseq2, sc);
        let left = Dp::new(&rq, &rr, sc);
        let f = |k: usize| (0..=refseq2.len()).map(|j| right.s[k][j]).max().unwrap();
        let g = |k: usize| (0..=m).map(|j| left.s[n - k][j]).max().unwrap();
        let want = (0..=n).map(|k| f(k) + g(k)).max().unwrap();
        let r = al.split_reference_alignment(query, refseq, refseq2, None);
        assert_eq!(r.score as i64, want, "split_reference: {}", ctx());
        rescore(&r.right, query, refseq2, sc, 0);
        rescore(&r.left, query, refseq, sc, 0);
    }

    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn dna(&mut self, len: usize) -> Vec<u8> {
            const A: &[u8] = b"ACGTACGTACGTNU";
            (0..len)
                .map(|_| A[(self.next() >> 16) as usize % A.len()])
                .collect()
        }
    }

    /// Three schemes: the single-piece one, the two-piece default, and the awkward one from
    /// `backend/tests.rs`.
    fn schemes() -> [Scores; 3] {
        [
            Scores {
                match_: 2,
                mismatch: 8,
                gap_open: [12, 12],
                gap_extend: [1, 1],
                end_bonus: 10,
            },
            Scores::default(),
            Scores {
                match_: 1,
                mismatch: 60,
                gap_open: [4, 20],
                gap_extend: [3, 0],
                end_bonus: 0,
            },
        ]
    }

    /// A reference with a deleted block and a short run of its bases left in the query, which is
    /// the shape the long level exists for and the one random pairs never produce.
    fn deletion_with_an_island(rng: &mut Rng, del: usize, island: usize) -> (Vec<u8>, Vec<u8>) {
        let reference = rng.dna(60 + del + 60);
        let mut query = Vec::new();
        query.extend_from_slice(&reference[..60]);
        query.extend_from_slice(&reference[60 + del / 2..60 + del / 2 + island]);
        query.extend_from_slice(&reference[60 + del..]);
        (query, reference)
    }

    /// Every mode, every scheme, against the plain DP.
    ///
    /// The length sweep straddles a chunk boundary at both lane widths and includes the empty side,
    /// which skips the kernel entirely and has its own closed form.
    #[test]
    fn every_mode_matches_a_textbook_dp() {
        let mut rng = Rng(0x1234_5678_9ABC_DEF0);
        for sc in schemes() {
            for &qlen in &[0usize, 1, 2, 5, 16, 17, 33, 48, 70] {
                for &rlen in &[0usize, 1, 3, 16, 31, 32, 64, 90] {
                    let q = rng.dna(qlen);
                    let r = rng.dna(rlen);
                    let r2 = rng.dna(rlen);
                    check_all(&q, &r, &r2, sc);
                }
            }
            for &(del, island) in &[(40usize, 3usize), (60, 5), (100, 7), (120, 2)] {
                let (q, r) = deletion_with_an_island(&mut rng, del, island);
                let r2 = rng.dna(40);
                check_all(&q, &r, &r2, sc);
            }
        }
    }

    /// The band is a separate code path (its own seals, its own edge fixes, its own argmax scan),
    /// so it is checked against the same DP confined to the same strip.
    ///
    /// Scores only, not CIGARs: many band-optimal paths tie and the oracle has no tie-break rule.
    /// The CIGAR is still re-scored, so a traceback that walked out of the strip would show up as a
    /// score that does not match its own alignment.
    #[test]
    fn a_band_matches_a_banded_textbook_dp() {
        let mut rng = Rng(0xA5A5_5A5A_C3C3_3C3C);
        for sc in schemes() {
            for &qlen in &[1usize, 5, 17, 40, 70] {
                for &rlen in &[1usize, 7, 20, 50, 90] {
                    let q = rng.dna(qlen);
                    let r = rng.dna(rlen);
                    let mut al = SimdAligner::new(sc);
                    for w in [1usize, 4, 13, 40] {
                        // The modes that are pinned at a far corner widen `w` rather than return
                        // nothing, so the oracle has to be told the width the kernel actually used.
                        let wg = w.max(1).max(qlen.abs_diff(rlen));
                        let dp = Dp::banded(&q, &r, sc, wg);
                        let g = al.global_alignment(&q, &r, Some(w));
                        assert_eq!(g.score as i64, dp.s[qlen][rlen], "global w={w} sc={sc:?}");
                        rescore(&g, &q, &r, sc, 0);

                        let wr = w.max(1).max(qlen.saturating_sub(rlen));
                        let dp = Dp::banded(&q, &r, sc, wr);
                        let lr = al.local_reference_end_alignment(&q, &r, Some(w));
                        let want = (0..=rlen).map(|j| dp.s[qlen][j]).max().unwrap();
                        assert_eq!(lr.score as i64, want, "local_reference_end w={w} sc={sc:?}");
                        rescore(&lr, &q, &r, sc, 0);

                        // `local_end` is never widened: its anchor `(0, 0)` is in every band.
                        let dp = Dp::banded(&q, &r, sc, w.max(1));
                        let le = al.local_end_alignment(&q, &r, Some(w));
                        let bonus = sc.end_bonus as i64;
                        let want = (0..=qlen)
                            .flat_map(|i| (0..=rlen).map(move |j| (i, j)))
                            .map(|(i, j)| {
                                if dp.s[i][j] <= NEG {
                                    NEG
                                } else {
                                    dp.s[i][j] + if i == qlen { bonus } else { 0 }
                                }
                            })
                            .max()
                            .unwrap()
                            .max(0);
                        assert_eq!(le.score as i64, want, "local_end w={w} sc={sc:?}");
                        rescore(
                            &le,
                            &q,
                            &r,
                            sc,
                            if le.query_end == qlen { bonus } else { 0 },
                        );
                    }
                }
            }
        }
    }
}
