//! Parallel, (mostly) in-place MSD radix sort for large arrays whose sort key
//! has uniformly distributed high bits (such as randstrobe hashes).
//!
//! The algorithm follows the block-based approach of IPS²Ra/IPS⁴o:
//!
//! 1. Each thread classifies a stripe of the input into buckets by the top
//!    `L1_BITS` bits of the key, using a small buffer per bucket. Whenever a
//!    buffer is full, it is written back into the stripe as one block.
//! 2. Blocks are moved to their bucket's destination (a permutation of
//!    fixed-size blocks, done in parallel).
//! 3. Each bucket is completed with the elements that remained in the
//!    buffers and then sorted independently (in parallel) using one
//!    out-of-place counting sort pass on the next key bits, followed by
//!    comparison sorts of the (tiny) resulting sub-buckets.
//!
//! Extra memory is O(threads × buckets × block size) plus a scratch buffer of
//! the size of a bucket per thread. Unusually large buckets (from repetitive
//! keys) are instead sorted in place with a parallel comparison sort.

use std::sync::atomic::{AtomicBool, Ordering};

use rayon::prelude::*;

/// Number of key bits used for the first (in-place, parallel) level
const L1_BITS: u32 = 10;
const N_BUCKETS: usize = 1 << L1_BITS;
/// Number of elements per block
const BLOCK: usize = 256;
/// Inputs shorter than this are sorted with a comparison sort
const MIN_RADIX_LEN: usize = 1 << 16;
const NONE: u32 = u32::MAX;
/// Buckets larger than LARGE_BUCKET_FACTOR times the average bucket size (and
/// at least MIN_LARGE_BUCKET_LEN) are sorted with a parallel comparison sort
const LARGE_BUCKET_FACTOR: usize = 8;
const MIN_LARGE_BUCKET_LEN: usize = 1 << 20;

#[derive(Clone, Copy)]
struct SyncPtr<T>(*mut T);
unsafe impl<T> Send for SyncPtr<T> {}
unsafe impl<T> Sync for SyncPtr<T> {}

impl<T> SyncPtr<T> {
    /// Copy block `src` to block `dst`.
    ///
    /// # Safety
    /// Both blocks must be in bounds and no other thread may access them
    /// concurrently.
    unsafe fn copy_block(self, src: usize, dst: usize) {
        unsafe {
            std::ptr::copy_nonoverlapping(self.0.add(src * BLOCK), self.0.add(dst * BLOCK), BLOCK);
        }
    }
}

struct Stripe<T> {
    /// Index of the first block slot of this stripe
    first_slot: usize,
    /// Bucket of each block written back into the stripe, in order
    block_buckets: Vec<u32>,
    /// Per-bucket buffers (N_BUCKETS × BLOCK) with the leftover elements
    buffers: Vec<T>,
    /// Number of leftover elements per bucket
    fill: Vec<usize>,
}

impl<T> Stripe<T> {
    fn leftovers(&self, bucket: usize) -> &[T] {
        &self.buffers[bucket * BLOCK..bucket * BLOCK + self.fill[bucket]]
    }
}

/// Sort `v` in parallel (using the current rayon thread pool).
///
/// `key` must be monotone with respect to the ordering of `T`, that is,
/// `a <= b` must imply `key(a) <= key(b)`. The result is the same as that
/// of `v.sort_unstable()`. Performance is best if the high bits of the key
/// are uniformly distributed.
pub fn radix_sort_by_key<T, F>(v: &mut [T], key: F)
where
    T: Copy + Ord + Default + Send + Sync,
    F: Fn(&T) -> u64 + Sync,
{
    let n = v.len();
    if n < MIN_RADIX_LEN {
        v.par_sort_unstable();
        return;
    }
    let bucket_of = |x: &T| (key(x) >> (64 - L1_BITS)) as usize;

    // Phase 1: Classify each stripe. Stripes start at block boundaries.
    let n_full_slots = n / BLOCK;
    let n_stripes = rayon::current_num_threads().clamp(1, n_full_slots);
    let mut stripe_slices = Vec::with_capacity(n_stripes);
    let mut rest = &mut v[..];
    for i in 0..n_stripes {
        let first_slot = i * n_full_slots / n_stripes;
        let len = if i + 1 == n_stripes {
            rest.len()
        } else {
            ((i + 1) * n_full_slots / n_stripes - first_slot) * BLOCK
        };
        let (stripe, tail) = std::mem::take(&mut rest).split_at_mut(len);
        stripe_slices.push((first_slot, stripe));
        rest = tail;
    }
    let stripes: Vec<Stripe<T>> = stripe_slices
        .into_par_iter()
        .map(|(first_slot, stripe)| classify(stripe, first_slot, &bucket_of))
        .collect();

    // Phase 2: Compute bucket boundaries and block destinations
    let mut full_blocks = vec![0usize; N_BUCKETS];
    let mut counts = vec![0usize; N_BUCKETS];
    for stripe in &stripes {
        for &b in &stripe.block_buckets {
            full_blocks[b as usize] += 1;
        }
        for (count, fill) in counts.iter_mut().zip(&stripe.fill) {
            *count += fill;
        }
    }
    let mut bucket_starts = Vec::with_capacity(N_BUCKETS + 1);
    let mut sum = 0;
    for b in 0..N_BUCKETS {
        counts[b] += full_blocks[b] * BLOCK;
        bucket_starts.push(sum);
        sum += counts[b];
    }
    bucket_starts.push(sum);
    debug_assert_eq!(sum, n);

    // A bucket's blocks go to the block slots that start within the bucket.
    // There are always enough of them, but the last block may extend into the
    // next bucket (it "spills"). If it extends beyond the end of the array,
    // it goes to a separate overflow block instead.
    let first_dest_slot: Vec<usize> = bucket_starts[..N_BUCKETS]
        .iter()
        .map(|&start| start.div_ceil(BLOCK))
        .collect();
    let n_slots = n.div_ceil(BLOCK);
    let mut dest = vec![NONE; n_slots];
    let mut inv = vec![NONE; n_slots];
    let mut next = first_dest_slot.clone();
    let mut overflow: Option<(usize, Vec<T>)> = None;
    for stripe in &stripes {
        for (j, &b) in stripe.block_buckets.iter().enumerate() {
            let src = stripe.first_slot + j;
            let d = next[b as usize];
            next[b as usize] += 1;
            if (d + 1) * BLOCK > n {
                overflow = Some((b as usize, v[src * BLOCK..(src + 1) * BLOCK].to_vec()));
            } else {
                dest[src] = d as u32;
                inv[d] = src as u32;
            }
        }
    }

    // Phase 3: Permute blocks. Following `inv` backwards from a free slot
    // that is the destination of some block gives a chain of moves that
    // ends at a slot that no block moves into. Chains are disjoint and can
    // be processed in parallel. The remaining blocks form cycles.
    let moved: Vec<AtomicBool> = (0..n_slots).map(|_| AtomicBool::new(false)).collect();
    let ptr = SyncPtr(v.as_mut_ptr());
    (0..n_slots)
        .into_par_iter()
        .filter(|&slot| inv[slot] != NONE && dest[slot] == NONE)
        .for_each(|end| {
            let mut cur = end;
            while inv[cur] != NONE {
                let src = inv[cur] as usize;
                // SAFETY: Each slot is on at most one chain
                unsafe { ptr.copy_block(src, cur) };
                moved[src].store(true, Ordering::Relaxed);
                cur = src;
            }
        });
    let mut tmp = vec![T::default(); BLOCK];
    for start in 0..n_slots {
        if dest[start] == NONE
            || dest[start] as usize == start
            || moved[start].load(Ordering::Relaxed)
        {
            continue;
        }
        tmp.copy_from_slice(&v[start * BLOCK..(start + 1) * BLOCK]);
        let mut cur = start;
        loop {
            let src = inv[cur] as usize;
            moved[src].store(true, Ordering::Relaxed);
            if src == start {
                v[cur * BLOCK..(cur + 1) * BLOCK].copy_from_slice(&tmp);
                break;
            }
            // SAFETY: Only this thread is running
            unsafe { ptr.copy_block(src, cur) };
            cur = src;
        }
    }

    // Phase 4: Complete each bucket. First save the spilled parts of blocks
    // since they are overwritten when the next bucket is completed.
    let block_ranges: Vec<(usize, usize)> = (0..N_BUCKETS)
        .map(|b| {
            if full_blocks[b] == 0 {
                return (bucket_starts[b], bucket_starts[b]);
            }
            let start = first_dest_slot[b] * BLOCK;
            let mut end = start + full_blocks[b] * BLOCK;
            if matches!(overflow, Some((ob, _)) if ob == b) {
                end -= BLOCK;
            }
            (start, end)
        })
        .collect();
    let spills: Vec<Vec<T>> = (0..N_BUCKETS)
        .into_par_iter()
        .map(|b| {
            let (_, end) = block_ranges[b];
            let bucket_end = bucket_starts[b + 1];
            if end > bucket_end {
                v[bucket_end..end].to_vec()
            } else {
                vec![]
            }
        })
        .collect();

    let mut bucket_slices = Vec::with_capacity(N_BUCKETS);
    let mut rest = &mut v[..];
    for &count in &counts {
        let (bucket, tail) = std::mem::take(&mut rest).split_at_mut(count);
        bucket_slices.push(bucket);
        rest = tail;
    }
    let large_bucket_len = (LARGE_BUCKET_FACTOR * n / N_BUCKETS).max(MIN_LARGE_BUCKET_LEN);
    bucket_slices
        .into_par_iter()
        .enumerate()
        .for_each_init(Vec::new, |scratch, (b, bucket)| {
            let bucket_start = bucket_starts[b];
            let (start, end) = block_ranges[b];
            let (start, end) = (
                start - bucket_start,
                end.min(bucket_starts[b + 1]) - bucket_start,
            );
            let mut elements = spills[b]
                .iter()
                .chain(
                    overflow
                        .iter()
                        .filter(|(ob, _)| *ob == b)
                        .flat_map(|(_, block)| block),
                )
                .chain(stripes.iter().flat_map(|stripe| stripe.leftovers(b)));
            let (head, rest) = bucket.split_at_mut(start);
            let tail = &mut rest[end - start..];
            for slot in head.iter_mut().chain(tail.iter_mut()) {
                *slot = *elements.next().unwrap();
            }
            debug_assert!(elements.next().is_none());

            if bucket.len() > large_bucket_len {
                // Avoid a large scratch buffer and a single-threaded
                // bottleneck for unusually large buckets (repetitive keys)
                bucket.par_sort_unstable();
            } else {
                sort_bucket(bucket, scratch, &key);
            }
        });
}

fn classify<T: Copy + Default>(
    stripe: &mut [T],
    first_slot: usize,
    bucket_of: &impl Fn(&T) -> usize,
) -> Stripe<T> {
    let mut buffers = vec![T::default(); N_BUCKETS * BLOCK];
    let mut fill = vec![0usize; N_BUCKETS];
    let mut block_buckets = vec![];
    let mut written = 0;
    for i in 0..stripe.len() {
        let x = stripe[i];
        let b = bucket_of(&x);
        let f = fill[b];
        buffers[b * BLOCK + f] = x;
        if f + 1 == BLOCK {
            // Writing is safe because we have read at least as many elements
            // as we write.
            stripe[written..written + BLOCK].copy_from_slice(&buffers[b * BLOCK..(b + 1) * BLOCK]);
            written += BLOCK;
            fill[b] = 0;
            block_buckets.push(b as u32);
        } else {
            fill[b] = f + 1;
        }
    }

    Stripe {
        first_slot,
        block_buckets,
        buffers,
        fill,
    }
}

/// Sort a bucket in which all keys have the same top L1_BITS bits
fn sort_bucket<T, F>(bucket: &mut [T], scratch: &mut Vec<T>, key: &F)
where
    T: Copy + Ord + Default,
    F: Fn(&T) -> u64,
{
    let m = bucket.len();
    if m <= 256 {
        bucket.sort_unstable();
        return;
    }
    // Aim for ~16 elements per sub-bucket
    let bits = (m / 16).ilog2().clamp(1, 64 - L1_BITS);
    let shift = 64 - L1_BITS - bits;
    let mask = (1u64 << bits) - 1;
    let digit = |x: &T| ((key(x) >> shift) & mask) as usize;

    let mut offsets = vec![0usize; (1 << bits) + 1];
    for x in bucket.iter() {
        offsets[digit(x) + 1] += 1;
    }
    for i in 1..offsets.len() {
        offsets[i] += offsets[i - 1];
    }
    scratch.clear();
    scratch.resize(m, T::default());
    let mut pos = offsets.clone();
    for x in bucket.iter() {
        let d = digit(x);
        scratch[pos[d]] = *x;
        pos[d] += 1;
    }
    for w in offsets.windows(2) {
        scratch[w[0]..w[1]].sort_unstable();
    }
    bucket.copy_from_slice(scratch);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check(mut v: Vec<(u64, u64)>, n_threads: usize) {
        let mut expected = v.clone();
        expected.sort_unstable();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build()
            .unwrap();
        pool.install(|| radix_sort_by_key(&mut v, |x| x.0));
        assert!(v == expected);
    }

    fn random_vec(n: usize, seed: u64, key_mask: u64) -> Vec<(u64, u64)> {
        let mut rng = fastrand::Rng::with_seed(seed);
        (0..n)
            .map(|_| (rng.u64(..) & key_mask, rng.u64(..1000)))
            .collect()
    }

    #[test]
    fn sorts_random() {
        for &n in &[
            0,
            1,
            1000,
            MIN_RADIX_LEN,
            MIN_RADIX_LEN + 1,
            300_001,
            1_000_000,
        ] {
            for &threads in &[1, 3, 8] {
                check(random_vec(n, n as u64, u64::MAX), threads);
            }
        }
    }

    #[test]
    fn sorts_skewed() {
        // Only a few distinct top bits and all keys equal
        check(random_vec(500_000, 1, 0xC000_0000_0000_00FF), 4);
        check(random_vec(500_000, 2, 0), 4);
        check(random_vec(500_000, 3, 0x0000_0000_FFFF_FFFF), 4);
        // One bucket large enough to take the parallel comparison sort path
        let mut v = random_vec(3_000_000, 5, u64::MAX);
        for x in v.iter_mut().take(1_500_000) {
            x.0 &= 0x00FF_FFFF_FFFF_FFFF;
        }
        check(v, 4);
    }

    #[test]
    fn sorts_presorted() {
        let mut v = random_vec(400_000, 4, u64::MAX);
        v.sort_unstable();
        check(v.clone(), 5);
        v.reverse();
        check(v, 5);
    }
}
