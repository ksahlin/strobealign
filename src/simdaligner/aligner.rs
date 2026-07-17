use super::Scores;

use super::avx2_u8::{U8Probe, fits_u8_lanes};
use super::{AlignmentResult, SplitReferenceAlignment};

/// A reusable aligner configured with a fixed [`Scores`] scheme.
///
/// Every alignment method takes `&mut self` because the scratch buffers are held internally and
/// reused; only the query and reference slices change between calls. Construct one per thread
/// and keep it - constructing one per alignment reallocates every buffer.
///
/// See the [module docs](super) for the guarantees, the alphabet, and what a band means.
pub struct SimdAligner {
    workspace: U8Probe,
    scores: Scores,
}

impl Default for SimdAligner {
    fn default() -> Self {
        SimdAligner::new(Scores::default())
    }
}

impl SimdAligner {
    /// Builds an aligner for the given scoring scheme.
    ///
    /// Every precondition is checked here, once per scheme, so that the alignment path carries
    /// no feature detection and no validation at all.
    ///
    /// # Panics
    ///
    /// - If the CPU does not support AVX2. See [`SimdAligner::is_supported`] to check first.
    /// - If `scores.gap_open < scores.gap_extend`. The kernel collapses the gap-open and
    ///   gap-extend cases into a single recurrence on that assumption.
    /// - If `match_ + 3 * gap_open - gap_extend > 255`, i.e. roughly `gap_open > 84` at
    ///   `match_ = 2`. Beyond that bound an intermediate value exceeds the kernel's 8-bit
    ///   lanes.
    pub fn new(scores: Scores) -> Self {
        assert!(
            Self::is_supported(),
            "this CPU doesn't support AVX2, which this aligner requires - check SimdAligner::is_supported() first"
        );
        // Hard asserts: once per scheme, never in the inner loop. The messages say why.
        assert!(
            scores.gap_open >= scores.gap_extend,
            "gap_open ({}) must be >= gap_extend ({}): the kernel collapses the open-vs-extend \
             recurrence on that assumption",
            scores.gap_open,
            scores.gap_extend
        );
        assert!(
            fits_u8_lanes(scores.match_, scores.gap_open, scores.gap_extend),
            "scores do not fit the u8 kernel: match ({}) + 3*gap_open ({}) - gap_extend ({}) \
             must be <= 255, i.e. gap_open <= ~84 at match=2. Past this an intermediate wraps \
             and the score comes back better than optimal, silently.",
            scores.match_,
            scores.gap_open,
            scores.gap_extend
        );
        SimdAligner {
            workspace: U8Probe::new(),
            scores,
        }
    }

    /// Whether this CPU supports the instruction set the aligner requires (AVX2).
    ///
    /// [`SimdAligner::new`] panics if this is false.
    pub fn is_supported() -> bool {
        is_x86_feature_detected!("avx2")
    }

    /// The scoring scheme this aligner was built with.
    pub fn scores(&self) -> Scores {
        self.scores
    }

    /// Global (end-to-end) alignment: both sequences are aligned in full, and terminal gaps are
    /// penalised.
    ///
    /// `query_start`/`query_end` and `ref_start`/`ref_end` always span the entire inputs.
    ///
    /// Empty inputs follow affine-gap semantics: aligning a length-`n` sequence against an empty
    /// one yields a single length-`n` gap scoring `-(gap_open + (n - 1) * gap_extend)`, and two
    /// empty inputs yield an empty alignment scoring 0.
    ///
    /// # `bandwidth`
    ///
    /// Global is the one mode pinned at **both** corners, and `(qlen, rlen)` sits `|qlen - rlen|`
    /// off the diagonal a band is centred on, so a narrower band contains no path at all and `w`
    /// is widened to reach it. Free when the two lengths are close - the shape global is normally
    /// handed - and on a lopsided pair it degrades to correct-but-slow rather than to wrong.
    pub fn global_alignment(
        &mut self,
        query: &[u8],
        reference: &[u8],
        bandwidth: Option<usize>,
    ) -> AlignmentResult {
        self.workspace.global_alignment(
            query,
            reference,
            self.scores.match_,
            self.scores.mismatch,
            self.scores.gap_open,
            self.scores.gap_extend,
            bandwidth,
        )
    }

    /// Alignment with a **local reference end**: the query is spanned in full, but the reference
    /// may end anywhere, and the unaligned reference tail costs nothing.
    ///
    /// Formally, the optimum over every end cell `(query.len(), j)` for `j` in
    /// `0..=reference.len()`. `query_start`/`query_end` therefore span the whole query,
    /// `ref_start` is 0, and `ref_end` is where the alignment chose to stop: the CIGAR covers
    /// exactly `reference[..ref_end]` and says nothing about the rest.
    ///
    /// `j = 0` is a legal end cell, and sometimes the winning one - spanning the query as one
    /// gap and aligning no reference at all can beat a run of mismatches, in which case the
    /// whole reference is the free end gap.
    ///
    /// # Tie-break
    ///
    /// Among equally-scoring end cells the *smallest* `ref_end` wins: don't consume reference
    /// you don't have to.
    ///
    /// # `bandwidth`
    ///
    /// Anchored at `(0, 0)`, like [`local_end_alignment`](SimdAligner::local_end_alignment), so
    /// `Some(w)` bounds `ref_end` to `[qlen - w, qlen + w]`. Size `w` from how much indel
    /// imbalance you expect, not from the reference's length.
    ///
    /// The query must still be spanned, which needs `w >= qlen - rlen`; `w` is widened to that
    /// where necessary, which only bites when the query is the longer side.
    pub fn local_reference_end_alignment(
        &mut self,
        query: &[u8],
        reference: &[u8],
        bandwidth: Option<usize>,
    ) -> AlignmentResult {
        self.workspace.local_reference_end_alignment(
            query,
            reference,
            self.scores.match_,
            self.scores.mismatch,
            self.scores.gap_open,
            self.scores.gap_extend,
            bandwidth,
        )
    }

    /// Alignment with a **local reference start**: the mirror of
    /// [`local_reference_end_alignment`](SimdAligner::local_reference_end_alignment). The query
    /// is spanned in full, but the reference may *begin* anywhere, and the unaligned reference
    /// prefix costs nothing.
    ///
    /// Formally, the optimum over every start cell `(0, s)` for `s` in `0..=reference.len()`,
    /// ending at `(query.len(), reference.len())`. `query_start`/`query_end` span the whole
    /// query, `ref_end` is always `reference.len()`, and `ref_start` is where the alignment
    /// chose to begin: the CIGAR covers exactly `reference[ref_start..]`.
    ///
    /// This is exact rather than an approximation: affine gap costs are symmetric under
    /// reversal, and so is literal byte identity, so the reversed problem has the same optimum.
    ///
    /// # Tie-break
    ///
    /// Among equally-scoring start cells the *largest* `ref_start` wins - the shortest reference
    /// span, mirroring the `_end` variant's rule.
    ///
    /// # `bandwidth`
    ///
    /// As [`local_reference_end_alignment`](SimdAligner::local_reference_end_alignment),
    /// **mirrored with the mode**: the band is anchored at the *ends* of the two inputs, where
    /// this alignment is pinned, and bounds `ref_start` to within `w` of `rlen - qlen`. `None`
    /// is exact.
    pub fn local_reference_start_alignment(
        &mut self,
        query: &[u8],
        reference: &[u8],
        bandwidth: Option<usize>,
    ) -> AlignmentResult {
        self.workspace.local_reference_start_alignment(
            query,
            reference,
            self.scores.match_,
            self.scores.mismatch,
            self.scores.gap_open,
            self.scores.gap_extend,
            bandwidth,
        )
    }

    /// **Local-end alignment**: both sequences start at 0, and the alignment ends wherever it
    /// scores best. Both trailing tails - query and reference - are free.
    ///
    /// Formally, the optimum over every end cell `(i, j)` of
    /// `H[i][j] + (if i == query.len() { scores.end_bonus } else { 0 })`. `query_end` and
    /// `ref_end` are wherever it stopped, and the CIGAR covers exactly `query[..query_end]`
    /// against `reference[..ref_end]`.
    ///
    /// # The score is never negative
    ///
    /// Cell `(0, 0)` is a legal end cell scoring 0, so a hopeless extension returns the empty
    /// alignment: spans `0..0` on both sequences, score 0, empty CIGAR. That is a normal result.
    ///
    /// # `end_bonus`
    ///
    /// Added to the score only if the alignment covers the whole query. It nudges the alignment
    /// towards spanning the query without forcing it, which is the difference between this and
    /// [`local_reference_end_alignment`](SimdAligner::local_reference_end_alignment). It never
    /// enters the recurrence: it is applied once, at the end, to choose between the best cell
    /// anywhere and the best cell that finishes the query.
    ///
    /// # Tie-break
    ///
    /// Longest extension: among equally-scoring end cells the largest `query_end` wins, then the
    /// largest `ref_end`.
    ///
    /// # How the CIGAR can end
    ///
    /// Both of these are properties of the objective, and worth knowing if you consume the
    /// coordinates:
    ///
    /// - A trailing insertion happens only if `end_bonus >= gap_open`. Ending with a `k`-base
    ///   insertion costs `gap_open + (k - 1) * gap_extend` to reach `(qlen, j)` from
    ///   `(qlen - k, j)`, and `(qlen - k, j)` is itself a legal end cell, so the bonus must
    ///   cover that or stopping short would score higher.
    /// - With `gap_open > 0` the CIGAR never ends with a deletion, by the same argument: a local
    ///   maximum cannot end in a gap.
    ///
    /// When a trailing insertion does occur, `ref_end` does not advance with it.
    ///
    /// # `bandwidth`
    ///
    /// This mode and [`local_start_alignment`](SimdAligner::local_start_alignment) are what the
    /// band was built for: they *finish* an alignment something else already anchored, so the
    /// start is known and the drift is bounded by biology rather than by the matrix. They are
    /// also the only two **never widened** - `(0, 0)` is in every band, so `Some(w)` is exactly
    /// `w`. Under a band the tie-break is still longest-extension and the score is still never
    /// negative.
    ///
    /// Below roughly `qlen ~ 1.2 * w` the band's bookkeeping costs more than it saves, so pass
    /// `None` for short extensions.
    pub fn local_end_alignment(
        &mut self,
        query: &[u8],
        reference: &[u8],
        bandwidth: Option<usize>,
    ) -> AlignmentResult {
        self.workspace.local_end_alignment(
            query,
            reference,
            self.scores.match_,
            self.scores.mismatch,
            self.scores.gap_open,
            self.scores.gap_extend,
            self.scores.end_bonus,
            bandwidth,
        )
    }

    /// **Local-start alignment**: the mirror of
    /// [`local_end_alignment`](SimdAligner::local_end_alignment). Both sequences end at their
    /// ends, and the alignment *begins* wherever it scores best; both leading prefixes are free.
    ///
    /// Formally, the optimum over every start cell `(a, b)` of the best alignment of
    /// `query[a..]` against `reference[b..]`, plus `end_bonus` if `a == 0`. `query_start` and
    /// `ref_start` are wherever it began; `query_end` and `ref_end` are always the ends of the
    /// inputs.
    ///
    /// # The score is never negative
    ///
    /// The start cell `(query.len(), reference.len())` is legal and scores 0, so a hopeless
    /// extension returns the empty alignment. Note it sits at the *end* of both sequences, the
    /// mirror of local-end's, which sits at the start.
    ///
    /// # Tie-break
    ///
    /// Longest extension, mirrored: among equally-scoring start cells the *smallest*
    /// `query_start` wins, then the smallest `ref_start`. Likewise a *leading* insertion
    /// requires `end_bonus >= gap_open`, and with `gap_open > 0` the CIGAR never begins with a
    /// deletion.
    ///
    /// # `bandwidth`
    ///
    /// As [`local_end_alignment`](SimdAligner::local_end_alignment), **mirrored with the mode**:
    /// the band is anchored at the *ends* of the two inputs, where this alignment is pinned, so
    /// it constrains start cells to `|(qlen - a) - (rlen - b)| <= w`. In both modes `w` bounds
    /// the drift from the end that something else already established.
    pub fn local_start_alignment(
        &mut self,
        query: &[u8],
        reference: &[u8],
        bandwidth: Option<usize>,
    ) -> AlignmentResult {
        self.workspace.local_start_alignment(
            query,
            reference,
            self.scores.match_,
            self.scores.mismatch,
            self.scores.gap_open,
            self.scores.gap_extend,
            self.scores.end_bonus,
            bandwidth,
        )
    }

    /// **Split-reference alignment**: one query aligned across two references, with a single
    /// jump between them. Built for detecting a tandem duplication that a linear aligner reports
    /// as an insertion.
    ///
    /// The alignment starts at `query[0]` against `right_reference[0]` and ends at the end of
    /// `left_reference`, jumping once from the right reference to the left. It returns two
    /// [`AlignmentResult`]s that **partition the query exactly**:
    ///
    /// - `right` covers `query[..k]` against a prefix of `right_reference`;
    /// - `left` covers `query[k..]` against a suffix of `left_reference`;
    ///
    /// where `k`, the jump point, is chosen to maximise the total score:
    ///
    /// ```text
    /// score = max over k in 0..=query.len() of  F[k] + G[k]
    ///   F[k] = the best alignment of query[..k] against any prefix of right_reference
    ///   G[k] = the best alignment of query[k..] against any suffix of left_reference
    /// ```
    ///
    /// Exact, like everything else here: both full matrices are computed and every jump point is
    /// considered. `end_bonus` is ignored - the query is spanned in full by construction, so a
    /// bonus for spanning it would be a constant added to every candidate.
    ///
    /// # The score can be negative, and an arm can be empty
    ///
    /// Unlike the local modes there is no "give up" option: the two arms must cover the query,
    /// so a query that neither reference explains still gets covered, at a cost. Two
    /// consequences worth knowing before consuming the coordinates:
    ///
    /// - `k = 0` and `k = query.len()` are legal. One arm then spans the whole query and the
    ///   other is the empty alignment.
    /// - An arm may consume no reference at all, reporting its query segment as one pure
    ///   insertion, with `ref_start == ref_end`. That is the right answer, not a failure: it is
    ///   how an unalignable segment gets covered.
    ///
    /// # Tie-break
    ///
    /// Total score decides. Among equal totals the **most balanced** split wins - the smallest
    /// `|left.score - right.score|` - because `5 + 5` describes a real duplication and
    /// `-20 + 30` does not, though both total 10. Among equally balanced splits, the smallest
    /// `k`. Each arm's own reference endpoint breaks ties the way its one-reference counterpart
    /// does.
    ///
    /// # `bandwidth`
    ///
    /// `Some(w)` is **one band, applied to both arms**, each measured from that arm's own anchor:
    /// the right arm from `(query[0], right_reference[0])`, the left arm from the ends of `query`
    /// and `left_reference`.
    ///
    /// **It does not constrain the jump point.** `k` is still chosen over every value the two
    /// bands admit, which is what this mode exists to find; the band says only that neither arm
    /// wanders far from its own diagonal on the way there.
    ///
    /// Both arms must together cover the query, and an arm reaches at most `reflen + w` query
    /// bases in band, so `w` is widened where necessary to
    /// `(qlen - left_reference.len() - right_reference.len()) / 2`.
    pub fn split_reference_alignment(
        &mut self,
        query: &[u8],
        left_reference: &[u8],
        right_reference: &[u8],
        bandwidth: Option<usize>,
    ) -> SplitReferenceAlignment {
        self.workspace.split_reference_alignment(
            query,
            left_reference,
            right_reference,
            self.scores.match_,
            self.scores.mismatch,
            self.scores.gap_open,
            self.scores.gap_extend,
            bandwidth,
        )
    }
}
