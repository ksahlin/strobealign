use std::cmp::{max, min, Reverse};
use rayon::prelude::*;
use crate::strobes::{RandstrobeIterator, RandstrobeParameters};
use crate::syncmers::{SyncmerParameters, SyncmerIterator};
use crate::fasta::RefSequence;

/* Settings that influence index creation */
pub struct IndexParameters {
    canonical_read_length: usize,
    pub syncmer: SyncmerParameters,
    pub randstrobe: RandstrobeParameters,
}

impl IndexParameters {
    pub fn new(canonical_read_length: usize, k: usize, s: usize, l: usize, u: usize, q: u64, max_dist: u8) -> Self {
        let w_min = max(1, k / (k - s + 1) + l);
        let w_max = k / (k - s + 1) + u;
        IndexParameters {
            canonical_read_length,
            syncmer: SyncmerParameters::new(k, s),
            randstrobe: RandstrobeParameters  { w_min, w_max, q, max_dist }
        }
    }
}

impl SyncmerParameters {
    /// Pick a suitable number of bits for indexing randstrobe start indices
    fn pick_bits(&self, size: usize) -> u8 {
        let estimated_number_of_randstrobes = size / (self.k - self.s + 1) as usize;
        // Two randstrobes per bucket on average
        // TOOD checked_ilog2 or ilog2
        ((estimated_number_of_randstrobes as f64).log2() as u32 - 1).clamp(8, 31) as u8
    }
}

type Packed = u32;
type RandstrobeHash = u64;
type bucket_index_t = u64;

#[derive(Default)]
struct IndexCreationStatistics {
    tot_strobemer_count: u64,
    tot_occur_once: u64,
    tot_high_ab: u64,
    tot_mid_ab: u64,
    index_cutoff: usize,
    distinct_strobemers: u64,
}

#[derive(PartialEq,Eq,PartialOrd,Ord)]
struct RefRandstrobe {
    pub hash: RandstrobeHash,
    pub position: u32,
    /// packed representation of ref_index and strobe offset
    packed: u32,
/*
    bool operator< (const RefRandstrobe& other) const {
        return hash < other.hash;
    }

    int reference_index() const {
        return m_packed >> bit_alloc;
    }

    int strobe2_offset() const {
        return m_packed & mask;
    }

private:
    static constexpr int bit_alloc = 8;
    static constexpr int mask = (1 << bit_alloc) - 1;
 */
}

impl RefRandstrobe {
    fn new(hash: RandstrobeHash, ref_index: u32, position: u32, offset: u8) -> Self {
        let packed: u32 = ((ref_index as u32) << 8) + offset as u32;
        RefRandstrobe { hash, position, packed }
    }
}

struct QueryRandstrobe {
    hash: RandstrobeHash,
    start: u32,
    end: u32,
    is_reverse: bool,
}

fn count_randstrobes(seq: &[u8], parameters: &IndexParameters) -> usize {
    let syncmer_iterator = SyncmerIterator::new(seq, parameters.syncmer.k, parameters.syncmer.s, parameters.syncmer.t);
    let n_syncmers = syncmer_iterator.count();

    // The last w_min syncmers do not result in a randstrobe
    if n_syncmers < parameters.randstrobe.w_min {
        0
    } else {
        n_syncmers - parameters.randstrobe.w_min
    }
}

fn count_all_randstrobes(references: &Vec<RefSequence>, parameters: &IndexParameters) -> Vec<usize> {
    references.par_iter().map(|refseq| count_randstrobes(&refseq.sequence, parameters)).collect()
}

// TODO UnpopulatedStrobemerIndex
pub struct StrobemerIndex<'a> {
    references: &'a Vec<RefSequence>,
    pub parameters: IndexParameters,
    stats: IndexCreationStatistics,

    /// no. of bits of the hash to use when indexing a randstrobe bucket
    bits: u8,

    filter_cutoff: usize,

    ///The randstrobes vector contains all randstrobes sorted by hash.
    /// The randstrobe_start_indices vector points to entries in the
    /// randstrobes vector. randstrobe_start_indices[x] is the index of the
    /// first entry in randstrobes whose top *bits* bits of its hash value are
    /// greater than or equal to x.
    ///
    /// randstrobe_start_indices has one extra guard entry at the end that
    /// is always randstrobes.size().
    randstrobes: Vec<RefRandstrobe>,
    randstrobe_start_indices: Vec<bucket_index_t>,
}

impl<'a> StrobemerIndex<'a> {
    pub fn new(references: &'a Vec<RefSequence>, parameters: IndexParameters, bits: Option<u8>) -> Self {
        let total_reference_length = references.iter().map(|r| r.sequence.len()).sum();
        let bits = bits.unwrap_or_else(|| parameters.syncmer.pick_bits(total_reference_length));
        let randstrobes = vec![];
        let randstrobe_start_indices = vec![];
        let stats = IndexCreationStatistics::default();
        let filter_cutoff = 0;
        StrobemerIndex { references, parameters, stats, bits, filter_cutoff, randstrobes, randstrobe_start_indices }
    }

    pub fn populate(&mut self, filter_fraction: f64) {
        // Timer count_hash;
        let randstrobe_counts = count_all_randstrobes(self.references, &self.parameters);
        // stats.elapsed_counting_hashes = count_hash.duration();

        let total_randstrobes = randstrobe_counts.iter().sum();
        // stats.tot_strobemer_count = total_randstrobes;

        // logger.debug() << "  Total number of randstrobes: " << total_randstrobes << '\n';
        // uint64_t memory_bytes = references.total_length() + sizeof(RefRandstrobe) * total_randstrobes + sizeof(bucket_index_t) * (1u << bits);
        // logger.debug() << "  Estimated total memory usage: " << memory_bytes / 1E9 << " GB\n";

        // if (total_randstrobes > std::numeric_limits<bucket_index_t>::max()) {
        //     throw std::range_error("Too many randstrobes");
        // }
        // Timer randstrobes_timer;
        // logger.debug() << "  Generating randstrobes ...\n";
        self.randstrobes.reserve_exact(total_randstrobes);
        // TODO  self.assign_all_randstrobes(&randstrobe_counts);

        // Fill randstrobes vector
        for (ref_index, refseq) in self.references.iter().enumerate() {
            let seq = &refseq.sequence;
            if seq.len() < self.parameters.randstrobe.w_max {
                continue;
            }
            let mut syncmer_iter = SyncmerIterator::new(seq, self.parameters.syncmer.k, self.parameters.syncmer.s, self.parameters.syncmer.t);
            let randstrobe_iter = RandstrobeIterator::new(&mut syncmer_iter, &self.parameters.randstrobe);

            for randstrobe in randstrobe_iter {
                let offset = randstrobe.strobe2_pos - randstrobe.strobe1_pos;
                self.randstrobes.push(
                    RefRandstrobe::new(
                        randstrobe.hash,
                        ref_index as u32,
                        randstrobe.strobe1_pos as u32,
                        offset as u8
                    )
                );
            }
        }
        // stats.elapsed_generating_seeds = randstrobes_timer.duration();

        // Timer sorting_timer;
        // logger.debug() << "  Sorting ...\n";
        // sort by hash values
        // TODO sort_unstable (stdlib) is also pdqsort
        pdqsort::sort(&mut self.randstrobes);
        // stats.elapsed_sorting_seeds = sorting_timer.duration();

        // Timer hash_index_timer;
        // logger.debug() << "  Indexing ...\n";

        let mut tot_high_ab = 0;
        let mut tot_mid_ab = 0;
        let mut strobemer_counts = Vec::<usize>::new();

        self.stats.tot_occur_once = 0;
        self.randstrobe_start_indices.reserve((1usize << self.bits) + 1);

        let mut unique_mers = if self.randstrobes.is_empty() { 0 } else { 1 };

        let mut prev_hash: RandstrobeHash = if self.randstrobes.is_empty() { 0 } else { self.randstrobes[0].hash };
        let mut count = 0;

        for position in 0..self.randstrobes.len() {
            let cur_hash = self.randstrobes[position].hash;
            if cur_hash == prev_hash {
                count += 1;
                continue;
            }
            unique_mers += 1;

            if count == 1 {
                self.stats.tot_occur_once += 1;
            } else {
                if count > 100 {
                    tot_high_ab += 1;
                } else {
                    tot_mid_ab += 1;
                }
                strobemer_counts.push(count);
            }
            count = 1;
            let cur_hash_n = cur_hash >> (64 - self.bits);
            while self.randstrobe_start_indices.len() <= cur_hash_n as usize {
                self.randstrobe_start_indices.push(position as bucket_index_t);
            }
            prev_hash = cur_hash;
        }
        // wrap up last entry
        if count == 1 {
            self.stats.tot_occur_once += 1;
        } else {
            if count > 100 {
                tot_high_ab += 1;
            } else {
                tot_mid_ab += 1;
            }
            strobemer_counts.push(count);
        }
        while self.randstrobe_start_indices.len() < ((1usize << self.bits) + 1) {
            self.randstrobe_start_indices.push(self.randstrobes.len() as bucket_index_t);
        }
        self.stats.tot_high_ab = tot_high_ab;
        self.stats.tot_mid_ab = tot_mid_ab;

        strobemer_counts.sort_unstable_by_key(|k| Reverse(*k));

        let index_cutoff = (unique_mers as f64 * filter_fraction) as usize;
        self.stats.index_cutoff = index_cutoff;
        self.filter_cutoff =
            if !strobemer_counts.is_empty() {
                usize::clamp(
                    if index_cutoff < strobemer_counts.len() { strobemer_counts[index_cutoff] } else { *strobemer_counts.last().unwrap() },
                    30,// cutoff is around 30-50 on hg38. No reason to have a lower cutoff than this if aligning to a smaller genome or contigs.
                    100 // limit upper cutoff for normal NAM finding - use rescue mode instead
                )
            } else {
                30usize
            };
        //stats.elapsed_hash_index = hash_index_timer.duration();
        self.stats.distinct_strobemers = unique_mers;
    }

    fn assign_all_randstrobes(&mut self, randstrobe_counts: &Vec<usize>) {
        let mut offset = 0;
        for ref_index in 0..self.references.len() {
            self.assign_randstrobes(ref_index, offset);
            offset += randstrobe_counts[ref_index];
        }
/*
        std::vector<std::thread> workers;
        std::atomic_size_t ref_index{0};
        for (size_t i = 0; i < n_threads; ++i) {
            workers.push_back(
                std::thread(
                    [&]() {
                        while (true) {
                            size_t j = ref_index.fetch_add(1);
                            if (j >= references.size()) {
                                break;
                            }
                            assign_randstrobes(j, offsets[j]);
                        }
                    })
            );
        }
        for (auto& worker : workers) {
            worker.join();
        }
 */
    }

    /// Compute randstrobes of one reference and assign them to the randstrobes
    /// vector starting from the given offset
    fn assign_randstrobes(&mut self, ref_index: usize, offset: usize) {
        let seq = &self.references[ref_index].sequence;
        if seq.len() < self.parameters.randstrobe.w_max {
            return;
        }
        let mut syncmer_iter = SyncmerIterator::new(seq, self.parameters.syncmer.k, self.parameters.syncmer.s, self.parameters.syncmer.t);
        let randstrobe_iter = RandstrobeIterator::new(&mut syncmer_iter, &self.parameters.randstrobe);

        /*for randstrobe in randstrobe_iter {

            RefRandstrobe::packed_t packed = ref_index << 8;
            packed = packed + (randstrobe.strobe2_pos - randstrobe.strobe1_pos);
            randstrobes[offset++] = RefRandstrobe{randstrobe.hash, randstrobe.strobe1_pos, packed};
        }*/
    }
}
