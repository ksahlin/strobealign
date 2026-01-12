//! A modified version of the `std::slice::partition_point` method that is
//! faster for our use case.

use std::cmp::Ordering::{Less, Greater, Equal};

/// This is a modified version of `std::slice::binary_search_by` that has
/// been modified to not be branchless, that is, the call to
/// `std::hint::select_unpredictable` has been replaced with a regular `if`.
///
/// This is quite a bit faster when searching buckets in the randstrobes vector.
///
/// The explanation is probably the following one that applies when the vector
/// to be searched is very large (quoting from
/// <https://en.algorithmica.org/hpc/data-structures/binary-search/>):
///
/// > The real question you need to ask is not why the branchless
/// > implementation is worse but why the branchy version is better.
/// > It happens because when you have branching, the CPU can speculate on one
/// > of the branches and start fetching either the left or the right key
/// > before it can even confirm that it is the right one — which effectively
/// > acts as implicit prefetching.
/// >
/// > For the branchless implementation, this doesn’t happen, as `cmov` is
/// > treated as every other instruction, and the branch predictor doesn’t try
/// > to peek into its operands to predict the future.
#[inline]
fn custom_binary_search_by<T, F>(s: &[T], mut f: F) -> Result<usize, usize>
where
    F: FnMut(&T) -> std::cmp::Ordering,
{
    let mut size = s.len();
    if size == 0 {
        return Err(0);
    }
    let mut base = 0usize;

    // This loop intentionally doesn't have an early exit if the comparison
    // returns Equal. We want the number of loop iterations to depend *only*
    // on the size of the input slice so that the CPU can reliably predict
    // the loop count.
    while size > 1 {
        let half = size / 2;
        let mid = base + half;

        // SAFETY: the call is made safe by the following invariants:
        // - `mid >= 0`: by definition
        // - `mid < size`: `mid = size / 2 + size / 4 + size / 8 ...`
        let cmp = f(unsafe { s.get_unchecked(mid) });

        // Changed compared to version in standard library:
        // Use a reguler if instead of
        //   std::hint::select_unpredictable(cmp == Greater, base, mid);
        // which was used to force a branchless comparison.
        base = if cmp == Greater { base } else { mid };

        // This is imprecise in the case where `size` is odd and the
        // comparison returns Greater: the mid element still gets included
        // by `size` even though it's known to be larger than the element
        // being searched for.
        //
        // This is fine though: we gain more performance by keeping the
        // loop iteration count invariant (and thus predictable) than we
        // lose from considering one additional element.
        size -= half;
    }

    // SAFETY: base is always in [0, size) because base <= mid.
    let cmp = f(unsafe { s.get_unchecked(base) });
    if cmp == Equal {
        // SAFETY: same as the `get_unchecked` above.
        unsafe { std::hint::assert_unchecked(base < s.len()) };
        Ok(base)
    } else {
        let result = base + (cmp == Less) as usize;
        // SAFETY: same as the `get_unchecked` above.
        // Note that this is `<=`, unlike the assume in the `Ok` path.
        unsafe { std::hint::assert_unchecked(result <= s.len()) };
        Err(result)
    }
}

#[inline]
pub fn custom_partition_point<T, P>(s: &[T], mut pred: P) -> usize
where
    P: FnMut(&T) -> bool,
{
    custom_binary_search_by(s, |x| if pred(x) { Less } else { Greater }).unwrap_or_else(|i| i)
}
