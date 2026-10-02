use crate::index::{BucketIndex, RandstrobeHash, RefRandstrobe, StrobemerIndex};
use crate::refseq::RefSequence;
use crate::seeding::syncmers::SeqAccess;
use crate::seeding::{RandstrobeIterator, SeedingParameters, SyncmerIterator, SyncmerParameters};

use std::fmt::{Display, Formatter};
use std::mem::MaybeUninit;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

use log::{debug, trace};
use rayon;
use rayon::slice::ParallelSliceMut;

/// Create a StrobemerIndex
pub fn make_index(
    refseq: &RefSequence,
    parameters: SeedingParameters,
    bits: u8,
    filter_fraction: f64,
    n_threads: usize,
) -> (StrobemerIndex, IndexCreationStatistics) {
    let mut stats = IndexCreationStatistics::default();
    let estimated_number_of_randstrobes = parameters.syncmer.estimate_number_of_syncmers(refseq);
    let total_length: usize = refseq.total_length();
    let memory_bytes: usize = total_length / 4  // 2 bits per nucleotide
        + size_of::<RefRandstrobe>() * estimated_number_of_randstrobes
        + size_of::<BucketIndex>() * (1usize << bits);
    debug!(
        "  Estimated total memory usage: {:.1} GB",
        memory_bytes as f64 / 1E9
    );

    let timer = Instant::now();
    debug!("  Generating randstrobes ...");
    let mut randstrobes = make_randstrobes_parallel(
        refseq,
        &parameters,
        estimated_number_of_randstrobes,
        n_threads,
    );

    debug!("  Generating seeds: {:.2} s", timer.elapsed().as_secs_f64());
    // stats.elapsed_generating_seeds = randstrobes_timer.duration();

    let timer = Instant::now();
    debug!("  Sorting ...");
    // TODO
    // ensure comparison function is branchless
    // Comment from C++ code:
    // Compare both hash and position to ensure that the order of the
    // RefRandstrobes in the index is reproducible no matter which sorting
    // function is used. This branchless comparison is faster than the
    // equivalent one using std::tie.
    // __uint128_t lhs = (static_cast<__uint128_t>(m_hash_offset_flag) << 64) | ((static_cast<uint64_t>(m_position) << 32) | m_ref_index);
    // __uint128_t rhs = (static_cast<__uint128_t>(other.m_hash_offset_flag) << 64) | ((static_cast<uint64_t>(other.m_position) << 32) | m_ref_index);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .unwrap();
    pool.install(|| randstrobes.par_sort_unstable());

    debug!("    Took {:.2} s", timer.elapsed().as_secs_f64());

    // Remove the sentinels from the end
    while randstrobes
        .pop_if(|r| *r == RefRandstrobe::sentinel())
        .is_some()
    {}
    randstrobes.shrink_to_fit();
    let total_randstrobes = randstrobes.len();
    stats.tot_strobemer_count = total_randstrobes as u64;

    trace!(
        "  Estimated number of randstrobes vs actual: {:.6}",
        estimated_number_of_randstrobes as f64 / total_randstrobes as f64
    );

    // stats.elapsed_sorting_seeds = sorting_timer.duration();

    let timer = Instant::now();
    debug!("  Generating hash table index ...");

    let mut tot_high_ab = 0;
    let mut tot_mid_ab = 0;

    stats.tot_occur_once = 0;
    let mut randstrobe_start_indices = Vec::with_capacity((1usize << bits) + 1);
    let mut unique_mers = usize::from(!randstrobes.is_empty());

    let mut prev_hash: RandstrobeHash = if randstrobes.is_empty() {
        0
    } else {
        randstrobes[0].hash()
    };
    let mut count = 1;

    if !randstrobes.is_empty() {
        randstrobe_start_indices.push(0);
    }

    // strobemer_counts[i] is how many strobemers occur i times,
    // except that `strobemer_counts[1000]` is how many strobemers occur
    // 1000 times *or more*.
    let mut strobemer_counts = [0usize; 1001];
    #[allow(clippy::needless_range_loop)]
    for position in 1..randstrobes.len() {
        let cur_hash = randstrobes[position].hash();
        if cur_hash == prev_hash {
            count += 1;
            continue;
        }
        unique_mers += 1;

        if count == 1 {
            stats.tot_occur_once += 1;
        } else {
            if count > 100 {
                tot_high_ab += 1;
            } else {
                tot_mid_ab += 1;
            }
            strobemer_counts[count.min(strobemer_counts.len() - 1)] += 1;
        }
        count = 1;
        let cur_hash_n = cur_hash >> (64 - bits);
        while randstrobe_start_indices.len() <= cur_hash_n as usize {
            randstrobe_start_indices.push(position as BucketIndex);
        }
        prev_hash = cur_hash;
    }
    // wrap up last entry
    if count == 1 {
        stats.tot_occur_once += 1;
    } else {
        if count > 100 {
            tot_high_ab += 1;
        } else {
            tot_mid_ab += 1;
        }
        strobemer_counts[count.min(strobemer_counts.len() - 1)] += 1;
    }
    strobemer_counts[1] = unique_mers;
    while randstrobe_start_indices.len() < ((1usize << bits) + 1) {
        randstrobe_start_indices.push(randstrobes.len() as BucketIndex);
    }
    stats.tot_high_ab = tot_high_ab;
    stats.tot_mid_ab = tot_mid_ab;

    let index_cutoff = (unique_mers as f64 * filter_fraction) as usize;
    stats.index_cutoff = index_cutoff;

    let mut total = 0;
    let mut filter_cutoff = 1;
    for i in (1..strobemer_counts.len()).rev() {
        total += strobemer_counts[i];
        if total >= index_cutoff {
            filter_cutoff = i;
            break;
        }
    }

    trace!(
        "Filter cutoff before clamping to [30, 100]: {}",
        filter_cutoff
    );
    let filter_cutoff = usize::clamp(
        filter_cutoff,
        30, // cutoff is around 30-50 on hg38. No reason to have a lower cutoff than this if aligning to a smaller genome or contigs.
        100, // limit upper cutoff for normal NAM finding - use rescue mode instead
    );
    //stats.elapsed_hash_index = hash_index_timer.duration();
    debug!("    Took {:.2} s", timer.elapsed().as_secs_f64());
    stats.distinct_strobemers = unique_mers as u64;

    (
        StrobemerIndex::new(
            parameters,
            bits,
            filter_cutoff,
            randstrobes,
            randstrobe_start_indices,
            refseq.starts.clone(),
        ),
        stats,
    )
}

fn make_randstrobes_parallel(
    refseq: &RefSequence,
    parameters: &SeedingParameters,
    estimated_number_of_randstrobes: usize,
    n_threads: usize,
) -> Vec<RefRandstrobe> {
    // Reserve slightly more memory than we estimate we need
    let n = (estimated_number_of_randstrobes as f64 * 1.01) as usize;
    let mut randstrobes = Vec::with_capacity(n);

    const SLICE_LENGTH: usize = 10000;
    let output_slice_index = AtomicUsize::new(0);

    let uninit = randstrobes.spare_capacity_mut();
    let mut slices = vec![];
    {
        let mut slice = uninit;
        while slice.len() > SLICE_LENGTH {
            let (left, right) = slice.split_at_mut(SLICE_LENGTH);
            slices.push(Arc::new(Mutex::new(left)));
            slice = right;
        }
        if !slice.is_empty() {
            slices.push(Arc::new(Mutex::new(slice)));
        }
    }

    let contig_index = AtomicUsize::new(0);
    // If we did not allocate a large enough randstrobes vector,
    // additional randstrobes are stored here.
    let overflow = Arc::new(Mutex::new(Vec::new()));

    thread::scope(|s| {
        for _ in 0..n_threads {
            s.spawn(|| {
                let index = output_slice_index.fetch_add(1, Ordering::SeqCst);
                if index >= slices.len() {
                    return;
                }
                let mut output_slice = slices[index].lock().unwrap();

                let mut i = 0;
                let mut is_overflowing = false;
                loop {
                    // Get index of a contig to work on
                    let j = contig_index.fetch_add(1, Ordering::SeqCst);
                    if j >= refseq.names.len() {
                        break;
                    }
                    let start = refseq.contig_start(j);
                    let seq = &refseq.contig(j);

                    let mut iter = make_randstrobe_iter(seq, parameters);
                    if !is_overflowing {
                        for randstrobe in iter.by_ref() {
                            let pos = randstrobe.strobe1_pos + start;
                            let offset = randstrobe.strobe2_pos - randstrobe.strobe1_pos;
                            let randstrobe = RefRandstrobe::new(randstrobe.hash, pos, offset as u8);
                            output_slice[i].write(randstrobe);
                            i += 1;
                            if i == output_slice.len() {
                                let index = output_slice_index.fetch_add(1, Ordering::SeqCst);
                                if index >= slices.len() {
                                    is_overflowing = true;
                                    break;
                                }

                                output_slice = slices[index].lock().unwrap();
                                i = 0;
                            }
                        }
                    }
                    // Since the estimated number of randstrobes is usually quite
                    // close to the actual number and because we overallocated
                    // a little bit, we should in practice very rarely end up
                    // in this path where randstrobes are pushed one by one
                    // onto a shared Vec, which is very slow.
                    if is_overflowing {
                        for randstrobe in iter {
                            let pos = randstrobe.strobe1_pos + start;
                            let offset = randstrobe.strobe2_pos - randstrobe.strobe1_pos;
                            let randstrobe = RefRandstrobe::new(randstrobe.hash, pos, offset as u8);

                            overflow.lock().unwrap().push(randstrobe);
                        }
                    }
                }

                // Since we work with uninitialized memory,
                // we need to fill the rest of the slice with something.
                // We use sentinels that will end up at the end of the
                // randstrobes Vec after sorting.
                output_slice[i..].fill_with(|| MaybeUninit::new(RefRandstrobe::sentinel()));
            });
        }
    });

    let index = output_slice_index.fetch_add(1, Ordering::SeqCst);
    let capacity = randstrobes.capacity();
    unsafe {
        randstrobes.set_len((index * SLICE_LENGTH).min(capacity));
    }
    trace!(
        "Pre-allocated randstrobes vector was too short by {} randstrobes",
        overflow.lock().unwrap().len()
    );
    randstrobes.extend_from_slice(&overflow.lock().unwrap());

    randstrobes
}

fn make_randstrobe_iter<S: SeqAccess>(
    seq: S,
    parameters: &SeedingParameters,
) -> RandstrobeIterator<SyncmerIterator<S>> {
    let syncmer_iter = SyncmerIterator::new(
        seq,
        parameters.syncmer.k,
        parameters.syncmer.s,
        parameters.syncmer.t,
    );

    RandstrobeIterator::new(syncmer_iter, parameters.randstrobe.clone())
}

impl SyncmerParameters {
    /// Pick a suitable number of bits for indexing randstrobe start indices
    pub fn pick_bits(&self, refseq: &RefSequence) -> u8 {
        // Two randstrobes per bucket on average
        // TOOD checked_ilog2 or ilog2
        ((self.estimate_number_of_syncmers(refseq) as f64).log2() as u32).clamp(9, 32) as u8 - 1
    }

    pub fn estimate_number_of_syncmers(&self, refseq: &RefSequence) -> usize {
        let total_length: usize = refseq.total_length();

        total_length / (self.k - self.s + 1) + 1
    }
}

#[derive(Default)]
pub struct IndexCreationStatistics {
    tot_strobemer_count: u64,
    tot_occur_once: u64,
    tot_high_ab: u64,
    tot_mid_ab: u64,
    index_cutoff: usize,
    distinct_strobemers: u64,
}

impl Display for IndexCreationStatistics {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Index statistics")?;
        writeln!(f, "  Total strobemers:    {:14}", self.tot_strobemer_count)?;
        writeln!(
            f,
            "  Distinct strobemers: {:14} (100.00%)",
            self.distinct_strobemers
        )?;
        writeln!(
            f,
            "    1 occurrence:      {:14} ({:6.2}%)",
            self.tot_occur_once,
            100.0 * self.tot_occur_once as f64 / self.distinct_strobemers as f64
        )?;
        writeln!(
            f,
            "    2..100 occurrences:{:14} ({:6.2}%)",
            self.tot_mid_ab,
            100.0 * self.tot_mid_ab as f64 / self.distinct_strobemers as f64
        )?;
        writeln!(
            f,
            "    >100 occurrences:  {:14} ({:6.2}%)",
            self.tot_high_ab,
            100.0 * self.tot_high_ab as f64 / self.distinct_strobemers as f64
        )?;
        if self.tot_high_ab >= 1 {
            writeln!(
                f,
                "Ratio distinct to highly abundant: {}",
                self.distinct_strobemers / self.tot_high_ab
            )?;
        }
        if self.tot_mid_ab >= 1 {
            writeln!(
                f,
                "Ratio distinct to non distinct: {}",
                self.distinct_strobemers / (self.tot_high_ab + self.tot_mid_ab)
            )?;
        }
        write!(f, "Filtered cutoff index: {}", self.index_cutoff)?;

        Ok(())
    }
}

#[cfg(test)]
mod test {
    use crate::io::fasta::read_ref;
    use crate::packed_seq::PackedSeq;

    use super::*;

    #[test]
    fn pick_bits() {
        let parameters = SyncmerParameters::try_new(20, 16).unwrap();
        let refseq = read_ref("tests/phix.fasta").unwrap();
        assert_eq!(parameters.pick_bits(&refseq), 9);
    }

    #[test]
    fn index_phix() {
        let refseq = read_ref("tests/phix.fasta").unwrap();
        let parameters = SeedingParameters::new(150);
        let bits = parameters.syncmer.pick_bits(&refseq);
        let (_index, stats) = make_index(&refseq, parameters, bits, 0.1, 1);
        assert!(stats.distinct_strobemers > 0);
        assert_eq!(stats.tot_strobemer_count, 1090);
    }

    #[test]
    fn index_empty_reference() {
        let refseq = RefSequence::new(PackedSeq::new(), vec![0], vec!["name".to_string()]).unwrap();
        let parameters = SeedingParameters::new(150);
        let bits = parameters.syncmer.pick_bits(&refseq);
        let (_index2, stats) = make_index(&refseq, parameters, bits, 0.1, 1);
        assert_eq!(stats.distinct_strobemers, 0);
    }
}
