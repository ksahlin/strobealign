use crate::index::{BucketIndex, RandstrobeHash, RefRandstrobe, StrobemerIndex};
use crate::io::fasta::RefSequence;
use crate::seeding::{RandstrobeIterator, SeedingParameters, SyncmerIterator, SyncmerParameters};

use std::cmp::{Reverse, min};
use std::fmt::{Display, Formatter};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

use log::{debug, trace};
use rayon;
use rayon::slice::ParallelSliceMut;

/// Create a StrobemerIndex
pub fn make_index<'a>(
    references: &'a [RefSequence],
    parameters: SeedingParameters,
    bits: Option<u8>,
    filter_fraction: f64,
    n_threads: usize,
) -> (StrobemerIndex<'a>, IndexCreationStatistics) {
    let total_reference_length = references.iter().map(|r| r.sequence.len()).sum();
    let bits = bits.unwrap_or_else(|| parameters.syncmer.pick_bits(total_reference_length));

    let timer = Instant::now();
    let randstrobe_counts = count_all_randstrobes(references, &parameters, n_threads);
    debug!("  Counting hashes: {:.2} s", timer.elapsed().as_secs_f64());
    // stats.elapsed_counting_hashes = count_hash.duration();
    let mut stats = IndexCreationStatistics::default();

    let total_randstrobes: usize = randstrobe_counts.iter().sum();
    stats.tot_strobemer_count = total_randstrobes as u64;

    debug!("  Total number of randstrobes: {}", total_randstrobes);
    let total_length: usize = references.iter().map(|refseq| refseq.sequence.len()).sum();
    let memory_bytes: usize = total_length
        + size_of::<RefRandstrobe>() * total_randstrobes
        + size_of::<BucketIndex>() * (1usize << bits);
    debug!(
        "  Estimated total memory usage: {:.1} GB",
        memory_bytes as f64 / 1E9
    );

    if total_randstrobes > BucketIndex::MAX as usize {
        panic!("Too many randstrobes");
    }
    let timer = Instant::now();
    debug!("  Generating randstrobes ...");
    let mut randstrobes =
        make_randstrobes_parallel(references, &parameters, &randstrobe_counts, n_threads);
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
    // stats.elapsed_sorting_seeds = sorting_timer.duration();

    let timer = Instant::now();
    debug!("  Generating hash table index ...");

    let mut tot_high_ab = 0;
    let mut tot_mid_ab = 0;
    let mut strobemer_counts = Vec::<usize>::new();

    stats.tot_occur_once = 0;
    let mut randstrobe_start_indices = Vec::with_capacity((1usize << bits) + 1);
    let mut unique_mers = u64::from(!randstrobes.is_empty());

    let mut prev_hash: RandstrobeHash = if randstrobes.is_empty() {
        0
    } else {
        randstrobes[0].hash()
    };
    let mut count = 1;

    if !randstrobes.is_empty() {
        randstrobe_start_indices.push(0);
    }
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
            strobemer_counts.push(count);
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
        strobemer_counts.push(count);
    }
    while randstrobe_start_indices.len() < ((1usize << bits) + 1) {
        randstrobe_start_indices.push(randstrobes.len() as BucketIndex);
    }
    stats.tot_high_ab = tot_high_ab;
    stats.tot_mid_ab = tot_mid_ab;

    strobemer_counts.sort_unstable_by_key(|k| Reverse(*k));

    let index_cutoff = (unique_mers as f64 * filter_fraction) as usize;
    stats.index_cutoff = index_cutoff;
    let filter_cutoff = if index_cutoff < strobemer_counts.len() {
        strobemer_counts[index_cutoff]
    } else {
        *strobemer_counts.last().unwrap_or(&30)
    };
    trace!(
        "Filter cutoff before clamping to [30, 100]: {}",
        filter_cutoff
    );
    let filter_cutoff = usize::clamp(
        filter_cutoff,
        30, // cutoff is around 30-50 on hg38. No reason to have a lower cutoff than this if aligning to a smaller genome or contigs.
        100, // limit upper cutoff for normal NAM finding - use rescue mode instead
    );
    let rescue_cutoff = min(filter_cutoff * 2, 1000);
    //stats.elapsed_hash_index = hash_index_timer.duration();
    debug!("    Took {:.2} s", timer.elapsed().as_secs_f64());
    stats.distinct_strobemers = unique_mers;

    (
        StrobemerIndex {
            references,
            parameters,
            bits,
            filter_cutoff,
            partial_filter_cutoff: filter_cutoff,
            rescue_cutoff,
            randstrobes,
            randstrobe_start_indices,
        },
        stats,
    )
}

/*
    fn make_randstrobes(&self, randstrobe_counts: &[usize]) -> Vec<RefRandstrobe> {
        let mut randstrobes = vec![RefRandstrobe::default(); randstrobe_counts.iter().sum()];

        // Fill randstrobes vector
        let mut offset = 0;
        for ref_index in 0..self.references.len() {
            self.assign_randstrobes(ref_index, &mut randstrobes[offset..offset+randstrobe_counts[ref_index]]);
            offset += randstrobe_counts[ref_index];
        }
        randstrobes
    }
*/

fn make_randstrobes_parallel(
    references: &[RefSequence],
    parameters: &SeedingParameters,
    randstrobe_counts: &[usize],
    n_threads: usize,
) -> Vec<RefRandstrobe> {
    let mut randstrobes = vec![RefRandstrobe::default(); randstrobe_counts.iter().sum()];
    let mut slices = vec![];
    {
        let mut slice = &mut randstrobes[..];
        for &mid in randstrobe_counts.iter().take(references.len() - 1) {
            let (left, right) = slice.split_at_mut(mid);
            slices.push(Arc::new(Mutex::new(left)));
            slice = right;
        }
        slices.push(Arc::new(Mutex::new(slice)));
    }
    let ref_index = AtomicUsize::new(0);
    thread::scope(|s| {
        for _ in 0..n_threads {
            s.spawn(|| {
                loop {
                    let j = ref_index.fetch_add(1, Ordering::SeqCst);
                    if j >= references.len() {
                        break;
                    }
                    assign_randstrobes(references, parameters, j, *slices[j].lock().unwrap());
                }
            });
        }
    });

    randstrobes
}

/// Compute randstrobes of one reference contig and assign them to the provided slice
fn assign_randstrobes(
    references: &[RefSequence],
    parameters: &SeedingParameters,
    ref_index: usize,
    randstrobes: &mut [RefRandstrobe],
) {
    let seq = &references[ref_index].sequence;
    if seq.len() < parameters.randstrobe.w_max {
        return;
    }
    let syncmer_iter = SyncmerIterator::new(
        seq,
        parameters.syncmer.k,
        parameters.syncmer.s,
        parameters.syncmer.t,
    );

    let randstrobe_iter = RandstrobeIterator::new(syncmer_iter, parameters.randstrobe.clone());

    let mut n = 0;
    for (i, randstrobe) in randstrobe_iter.enumerate() {
        n += 1;
        let offset = randstrobe.strobe2_pos - randstrobe.strobe1_pos;
        randstrobes[i] = RefRandstrobe::new(
            randstrobe.hash,
            ref_index as u32,
            randstrobe.strobe1_pos as u32,
            offset as u8,
        );
    }
    debug_assert_eq!(n, randstrobes.len());
}

fn count_all_randstrobes(
    references: &[RefSequence],
    parameters: &SeedingParameters,
    n_threads: usize,
) -> Vec<usize> {
    let counts = vec![0; references.len()];
    let mutex = Mutex::new(counts);
    let ref_index = AtomicUsize::new(0);
    thread::scope(|s| {
        for _ in 0..n_threads {
            s.spawn(|| {
                loop {
                    let j = ref_index.fetch_add(1, Ordering::SeqCst);
                    if j >= references.len() {
                        break;
                    }
                    let count = count_randstrobes(&references[j].sequence, parameters);
                    mutex.lock().unwrap()[j] = count;
                }
            });
        }
    });

    mutex.into_inner().unwrap()
}

/// Count randstrobes by counting syncmers
fn count_randstrobes(seq: &[u8], parameters: &SeedingParameters) -> usize {
    let syncmer_iterator = SyncmerIterator::new(
        seq,
        parameters.syncmer.k,
        parameters.syncmer.s,
        parameters.syncmer.t,
    );

    syncmer_iterator.count()
}

impl SyncmerParameters {
    /// Pick a suitable number of bits for indexing randstrobe start indices
    pub fn pick_bits(&self, size: usize) -> u8 {
        let estimated_number_of_randstrobes = size / (self.k - self.s + 1) + 1;
        // Two randstrobes per bucket on average
        // TOOD checked_ilog2 or ilog2
        ((estimated_number_of_randstrobes as f64).log2() as u32).clamp(9, 32) as u8 - 1
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
    use crate::{
        indexer::make_index,
        io::fasta::{RefSequence, read_ref},
        revcomp::reverse_complement,
        seeding::{SeedingParameters, SyncmerParameters},
    };

    #[test]
    fn test_pick_bits() {
        let parameters = SyncmerParameters::try_new(20, 16).unwrap();
        assert_eq!(parameters.pick_bits(0), 8);
    }

    #[test]
    fn test_index_phix() {
        let references = read_ref("tests/phix.fasta").unwrap();
        let parameters = SeedingParameters::new(150);
        let (_index, stats) = make_index(&references, parameters, None, 0.1, 1);
        assert!(stats.distinct_strobemers > 0);
    }

    #[test]
    fn test_index_empty_reference() {
        let references = vec![RefSequence {
            name: "name".to_string(),
            sequence: vec![],
        }];
        let parameters = SeedingParameters::new(150);
        let (_index2, stats) = make_index(&references, parameters, None, 0.1, 1);
        assert_eq!(stats.distinct_strobemers, 0);
    }

    #[test]
    fn test_partial_orientation() {
        let references = read_ref("tests/phix.fasta").unwrap();
        let seq = &references[0].sequence;
        let rc_seq = reverse_complement(seq);
        let rc_references = vec![RefSequence {
            name: "phix_rc".to_string(),
            sequence: rc_seq,
        }];

        let parameters = SeedingParameters::new(300);

        let (fwd_index, _stats) = make_index(&references, parameters.clone(), None, 0.0000001, 1);

        let (rc_index, _stats) = make_index(&rc_references, parameters.clone(), None, 0.0000001, 1);

        assert_eq!(fwd_index.randstrobes.len(), rc_index.randstrobes.len());

        // Iterate over fwd index entries and look up each hash in the rc index
        for i in 0..fwd_index.randstrobes.len() {
            let fwd_hash = fwd_index.randstrobes[i].hash();

            let rev_pos_result = rc_index.get_partial(fwd_hash);
            assert!(rev_pos_result.is_some());
            let rev_pos = rev_pos_result.unwrap();

            let rc_partial_query = fwd_index.randstrobes[i].hash()
                ^ ((1 as u64) << rc_index.parameters.randstrobe.partial_orientation_pos);
            let rev_pos_forward = rc_index.get_partial_forward_from(rc_partial_query, rev_pos);
            assert!(rev_pos_forward.is_some());
        }
    }
}
