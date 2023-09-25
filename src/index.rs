use std::cmp::{max, min, Reverse};
use std::mem;
use std::time::Instant;
use log::debug;
use crate::strobes::{RandstrobeIterator, RandstrobeParameters};
use crate::syncmers::{SyncmerParameters, SyncmerIterator};
use crate::fasta::RefSequence;

/// Pre-defined index parameters that work well for a certain
/// "canonical" read length (and similar read lengths)
struct Profile {
    canonical_read_length: usize,
    r_threshold: usize,
    k: usize,
    s_offset: isize,
    l: isize,
    u: isize,
}

static PROFILES: [Profile; 7] = [
    Profile { canonical_read_length:  50, r_threshold:  90, k: 20, s_offset: -4, l: -3, u:  2},
    Profile { canonical_read_length: 100, r_threshold: 110, k: 20, s_offset: -4, l: -2, u:  2},
    Profile { canonical_read_length: 125, r_threshold: 135, k: 20, s_offset: -4, l: -1, u:  4},
    Profile { canonical_read_length: 150, r_threshold: 175, k: 20, s_offset: -4, l:  1, u:  7},
    Profile { canonical_read_length: 250, r_threshold: 275, k: 20, s_offset: -4, l:  4, u: 13},
    Profile { canonical_read_length: 300, r_threshold: 375, k: 22, s_offset: -4, l:  2, u: 12},
    Profile { canonical_read_length: 400, r_threshold: usize::MAX, k: 23, s_offset: -6,  l: 2, u: 12},
];

/* Settings that influence index creation */
#[derive(Debug)]
pub struct IndexParameters {
    canonical_read_length: usize,
    pub syncmer: SyncmerParameters,
    pub randstrobe: RandstrobeParameters,
}

impl IndexParameters {
    pub fn new(canonical_read_length: usize, k: usize, s: usize, l: isize, u: isize, q: u64, max_dist: u8) -> Self {
        let w_min = max(0, (k / (k - s + 1)) as isize + l) as usize;
        let w_max = ((k / (k - s + 1)) as isize + u) as usize;
        IndexParameters {
            canonical_read_length,
            syncmer: SyncmerParameters::new(k, s),
            randstrobe: RandstrobeParameters  { w_min, w_max, q, max_dist }
        }
    }
    /// Create an IndexParameters instance based on a given read length.
    /// k, s, l, u, c and max_seed_len can be used to override determined parameters
    /// by setting them to a value other than IndexParameters::DEFAULT.
    pub fn from_read_length(
        read_length: usize,
        mut k: Option<usize>,
        mut s: Option<usize>,
        mut l: Option<isize>,
        mut u: Option<isize>,
        c: Option<u32>,
        max_seed_len: Option<usize>,
    ) -> IndexParameters {
        let default_c = 8;
        let mut canonical_read_length = 50;
        for profile in &PROFILES {
            if read_length <= profile.r_threshold {
                if k.is_none() {
                    k = Some(profile.k);
                }
                if s.is_none() {
                    s = Some((k.unwrap() as isize + profile.s_offset) as usize);
                }
                if l.is_none() {
                    l = Some(profile.l);
                }
                if u.is_none() {
                    u = Some(profile.u);
                }
                canonical_read_length = profile.canonical_read_length;
                break;
            }
        }

        let k = k.unwrap();
        let s = s.unwrap();
        let l = l.unwrap();
        let u = u.unwrap();

        let max_dist = match max_seed_len {
            Some(max_seed_len) => { (max_seed_len - k) as u8 },
            None => { usize::clamp(canonical_read_length.saturating_sub(70), k, 255) as u8 }
        };
        let q = 2u64.pow(c.unwrap_or(default_c)) - 1;

        IndexParameters::new(canonical_read_length, k, s, l, u, q, max_dist)
    }
}

impl SyncmerParameters {
    /// Pick a suitable number of bits for indexing randstrobe start indices
    fn pick_bits(&self, size: usize) -> u8 {
        let estimated_number_of_randstrobes = size / (self.k - self.s + 1);
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
pub struct RefRandstrobe {
    hash: RandstrobeHash,
    position: u32,
    /// packed representation of reference index and offset of second strobe
    packed: u32,
}

const REF_RANDSTROBE_OFFSET_BITS: u32 = 8;
const REF_RANDSTROBE_MASK: u32 = (1 << REF_RANDSTROBE_OFFSET_BITS) - 1;

impl RefRandstrobe {
    fn new(hash: RandstrobeHash, ref_index: u32, position: u32, offset: u8) -> Self {
        let packed = (ref_index << 8) + offset as u32;
        RefRandstrobe { hash, position, packed }
    }

    pub fn hash(&self) -> RandstrobeHash { self.hash }

    pub fn position(&self) -> usize { self.position as usize }

    pub fn reference_index(&self) -> usize {
        (self.packed >> REF_RANDSTROBE_OFFSET_BITS) as usize
    }

    pub fn strobe2_offset(&self) -> usize {
        (self.packed & REF_RANDSTROBE_MASK) as usize
    }
}

struct QueryRandstrobe {
    hash: RandstrobeHash,
    start: usize,
    end: usize,
    is_reverse: bool,
}

fn count_randstrobes(seq: &[u8], parameters: &IndexParameters) -> usize {
    let syncmer_iterator = SyncmerIterator::new(seq, parameters.syncmer.k, parameters.syncmer.s, parameters.syncmer.t);
    let n_syncmers = syncmer_iterator.count();

    // The last w_min syncmers do not result in a randstrobe
    n_syncmers.saturating_sub(parameters.randstrobe.w_min)
}

fn count_all_randstrobes(references: &[RefSequence], parameters: &IndexParameters) -> Vec<usize> {
    references.iter().map(|refseq| count_randstrobes(&refseq.sequence, parameters)).collect()
}

// TODO UnpopulatedStrobemerIndex
pub struct StrobemerIndex<'a> {
    references: &'a [RefSequence],
    pub parameters: IndexParameters,
    stats: IndexCreationStatistics,

    /// no. of bits of the hash to use when indexing a randstrobe bucket
    bits: u8,

    /// Regular (non-rescue) NAM finding ignores randstrobes that occur more often than
    /// this (see StrobemerIndex::is_filtered())
    pub filter_cutoff: usize,

    pub rescue_cutoff: usize,

    /// The randstrobes vector contains all randstrobes sorted by hash.
    /// The randstrobe_start_indices vector points to entries in the
    /// randstrobes vector. `randstrobe_start_indices[x]` is the index of the
    /// first entry in randstrobes whose top *bits* bits of its hash value are
    /// greater than or equal to x.
    ///
    /// randstrobe_start_indices has one extra guard entry at the end that
    /// is always randstrobes.len().
    pub randstrobes: Vec<RefRandstrobe>,
    randstrobe_start_indices: Vec<bucket_index_t>,
}

impl<'a> StrobemerIndex<'a> {
    pub fn new(references: &'a [RefSequence], parameters: IndexParameters, bits: Option<u8>) -> Self {
        let total_reference_length = references.iter().map(|r| r.sequence.len()).sum();
        let bits = bits.unwrap_or_else(|| parameters.syncmer.pick_bits(total_reference_length));
        let randstrobes = vec![];
        let randstrobe_start_indices = vec![];
        let stats = IndexCreationStatistics::default();
        let filter_cutoff = 0;
        let rescue_cutoff = 0;
        StrobemerIndex { references, parameters, stats, bits, filter_cutoff, rescue_cutoff, randstrobes, randstrobe_start_indices }
    }

    pub fn populate(&mut self, filter_fraction: f64, rescue_level: usize) {
        let timer = Instant::now();
        let randstrobe_counts = count_all_randstrobes(self.references, &self.parameters);
        debug!("  Counting hashes: {:.2} s", timer.elapsed().as_secs_f64());
        // stats.elapsed_counting_hashes = count_hash.duration();

        let total_randstrobes = randstrobe_counts.iter().sum();
        // stats.tot_strobemer_count = total_randstrobes;

        debug!("  Total number of randstrobes: {}", total_randstrobes);
        let total_length: usize = self.references.iter().map(|refseq| refseq.sequence.len()).sum();
        let memory_bytes: usize = total_length + mem::size_of::<RefRandstrobe>() * total_randstrobes + mem::size_of::<bucket_index_t>() * (1usize << self.bits);
        debug!("  Estimated total memory usage: {:.1} GB", memory_bytes as f64 / 1E9);

        if total_randstrobes > bucket_index_t::MAX as usize {
            panic!("Too many randstrobes");
        }
        let timer = Instant::now();
        debug!("  Generating randstrobes ...");
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
        debug!("  Generating seeds: {:.2} s", timer.elapsed().as_secs_f64());
        // stats.elapsed_generating_seeds = randstrobes_timer.duration();

        let timer = Instant::now();
        debug!("  Sorting ...");
        self.randstrobes.sort_unstable_by_key(|r| r.hash);

        debug!("    Took {:.2} s", timer.elapsed().as_secs_f64());
        // stats.elapsed_sorting_seeds = sorting_timer.duration();

        let timer = Instant::now();
        debug!("  Indexing ...");

        let mut tot_high_ab = 0;
        let mut tot_mid_ab = 0;
        let mut strobemer_counts = Vec::<usize>::new();

        self.stats.tot_occur_once = 0;
        self.randstrobe_start_indices.reserve((1usize << self.bits) + 1);

        let mut unique_mers = u64::from(!self.randstrobes.is_empty());

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
            usize::clamp(
                if index_cutoff < strobemer_counts.len() { strobemer_counts[index_cutoff] } else { *strobemer_counts.last().unwrap_or(&30) },
                30,// cutoff is around 30-50 on hg38. No reason to have a lower cutoff than this if aligning to a smaller genome or contigs.
                100 // limit upper cutoff for normal NAM finding - use rescue mode instead
            );
        self.rescue_cutoff = min(self.filter_cutoff * 2, 1000);
        //stats.elapsed_hash_index = hash_index_timer.duration();
        debug!("    Took {:.2} s", timer.elapsed().as_secs_f64());
        self.stats.distinct_strobemers = unique_mers;
    }

    fn assign_all_randstrobes(&mut self, randstrobe_counts: &[usize]) {
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

    /// Find index of first entry in randstrobe table that has the given hash value
    pub fn get(&self, hash: RandstrobeHash) -> Option<usize> {
        const MAX_LINEAR_SEARCH: usize = 4;
        let top_n = (hash >> (64 - self.bits)) as usize;
        let position_start = self.randstrobe_start_indices[top_n];
        let position_end = self.randstrobe_start_indices[top_n + 1];
        let bucket = &self.randstrobes[position_start as usize .. position_end as usize];
        if bucket.is_empty() {
            return None;
        } else if bucket.len() < MAX_LINEAR_SEARCH {
            for (pos, randstrobe) in bucket.iter().enumerate() {
                if randstrobe.hash == hash { return Some(position_start as usize + pos); }
                if randstrobe.hash > hash { return None; }
            }
            return None;
        }

        let pos = bucket.partition_point(|h| h.hash < hash);
        if pos < bucket.len() && bucket[pos].hash == hash {
            Some(position_start as usize + pos)
        } else {
            None
        }
    }

    pub fn get_count(&self, position: usize) -> usize {
        const MAX_LINEAR_SEARCH: usize = 8;
        let key = self.randstrobes[position].hash;
        let top_n = (key >> (64 - self.bits)) as usize;
        let position_end = self.randstrobe_start_indices[top_n + 1] as usize;

        if position_end - position < MAX_LINEAR_SEARCH {
            let mut count = 1;
            for position_start in position+1..position_end {
                if self.randstrobes[position_start].hash == key {
                    count += 1;
                } else {
                    break;
                }
            }
            count
        } else {
            let bucket = &self.randstrobes[position..position_end];
            bucket.partition_point(|h| h.hash <= key)
        }
    }

    /// Return whether the randstrobe at given position occurs more often than cutoff
    pub fn is_too_frequent(&self, position: usize, cutoff: usize) -> bool {
        if position + self.filter_cutoff < self.randstrobes.len() {
            self.randstrobes[position].hash == self.randstrobes[position + cutoff].hash
        } else {
            false
        }
    }

}
