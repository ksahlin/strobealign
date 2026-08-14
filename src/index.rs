use std::fs::File;
use std::io::{BufRead, BufReader, Error, Read, Write};
use std::path::Path;

use thiserror::Error;

use crate::partition::custom_partition_point;
use crate::seeding::{
    InvalidSeedingParameter, RandstrobeParameters, SeedingParameters, SyncmerParameters,
};

pub type RandstrobeHash = u64;
pub type BucketIndex = u64;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Default, Clone)]
#[repr(C)]
pub struct RefRandstrobe {
    /// Packed representation of the hash and the strobe offset.  
    /// Has the following layout
    /// `<strobe1 hash><strobe1 orientation bit><strobe2 hash><strobe2 orientation bit><strobe2 offset from strobe1>`
    ///
    /// [`RandstrobeParameters::partial_orientation_pos`] specifies the position of strobe1 orientation bit in the layout  
    /// [`RandstrobeParameters::forward_main_hash_mask`] is the mask for strobe1 hash (without the orientation bit)  
    /// [`RandstrobeParameters::main_hash_mask`] is the mask for strobe1 hash which includes the orientation bit  
    /// [`REF_RANDSTROBE_HASH_MASK`] is the mask for strobe1 and strobe2 hashes (including their orientations)  
    /// The [`STROBE2_OFFSET_BITS`] constant specifies the number of bits reserved for the offset  
    hash_offset: u64,
    position: u32,
    ref_index: u32,
}

/// Mask for the part of the randstrobe hash that includes individual strobe hashes and orientations
pub const REF_RANDSTROBE_HASH_MASK: u64 = 0xFFFFFFFFFFFFFF00;
/// Number of bits reserved for offset between first and second strobe
pub const STROBE2_OFFSET_BITS: u32 = 8;
/// Mask for the part of the randstrobe hash that includes the offset between first and second strobe
pub const STROBE2_OFFSET_MASK: u64 = (1u64 << STROBE2_OFFSET_BITS) - 1;
pub const REF_RANDSTROBE_MAX_NUMBER_OF_REFERENCES: usize = u32::MAX as usize;

impl RefRandstrobe {
    pub fn new(hash: RandstrobeHash, ref_index: u32, position: u32, offset: u8) -> Self {
        let hash_offset = (hash & REF_RANDSTROBE_HASH_MASK) | (offset as u64);
        RefRandstrobe {
            hash_offset,
            position,
            ref_index,
        }
    }

    pub fn hash(&self) -> RandstrobeHash {
        self.hash_offset & REF_RANDSTROBE_HASH_MASK
    }

    pub fn position(&self) -> usize {
        self.position as usize
    }

    pub fn reference_index(&self) -> usize {
        self.ref_index as usize
    }

    pub fn strobe2_offset(&self) -> usize {
        (self.hash_offset & STROBE2_OFFSET_MASK) as usize
    }
}

pub struct StrobemerIndex {
    pub parameters: SeedingParameters,

    /// no. of bits of the hash to use when indexing a randstrobe bucket
    bits: u8,

    /// Regular (non-rescue) NAM finding ignores randstrobes that occur more often than
    /// this (see StrobemerIndex::is_too_frequent())
    filter_cutoff: usize,

    /// The randstrobes vector contains all randstrobes sorted by hash.
    /// The bucket_starts vector points to entries in the
    /// randstrobes vector. `bucket_starts[x]` is the index of the
    /// first entry in randstrobes whose top *bits* bits of its hash value are
    /// greater than or equal to x.
    ///
    /// bucket_starts has one extra guard entry at the end that
    /// is always randstrobes.len().
    randstrobes: Vec<RefRandstrobe>,
    bucket_starts: Vec<BucketIndex>,
}

pub struct IndexEntry<'a> {
    strobemer_index: &'a StrobemerIndex,
    pub position: usize,
}

#[allow(clippy::len_without_is_empty)]
impl StrobemerIndex {
    pub fn new(
        parameters: SeedingParameters,
        bits: u8,
        filter_cutoff: usize,
        randstrobes: Vec<RefRandstrobe>,
        bucket_starts: Vec<BucketIndex>,
    ) -> Self {
        StrobemerIndex {
            parameters,
            bits,
            filter_cutoff,
            randstrobes,
            bucket_starts,
        }
    }

    pub fn len(&self) -> usize {
        self.randstrobes.len()
    }
    pub fn filter_cutoff(&self) -> usize {
        self.filter_cutoff
    }

    pub fn entry(&'_ self, position: usize) -> IndexEntry<'_> {
        IndexEntry {
            strobemer_index: self,
            position,
        }
    }

    // Find the first entry that matches the forwald full hash (including orientation bits)
    pub fn get_full_forward(&'_ self, hash: RandstrobeHash) -> Option<IndexEntry<'_>> {
        self.get_masked(hash, REF_RANDSTROBE_HASH_MASK)
    }

    /// Find the first entry that matches the undirected main hash (without orientation bit)
    pub fn get_partial(&'_ self, hash: RandstrobeHash) -> Option<IndexEntry<'_>> {
        self.get_masked(hash, self.parameters.randstrobe.main_hash_mask)
    }

    /// Find the first entry matching the forward main hash
    pub fn get_partial_forward(&'_ self, hash: RandstrobeHash) -> Option<IndexEntry<'_>> {
        self.get_masked(hash, self.parameters.randstrobe.forward_main_hash_mask)
    }

    /// Find the first entry matching the forward main hash, starting from
    /// the undirected main position
    pub fn get_partial_forward_from(
        &'_ self,
        hash: RandstrobeHash,
        undirected_position: usize,
    ) -> Option<IndexEntry<'_>> {
        self.get_masked_from(
            hash,
            self.parameters.randstrobe.forward_main_hash_mask,
            Some(undirected_position),
        )
    }

    /// Find index of first entry in randstrobe table that has the given
    /// hash value masked by the `hash_mask`.
    /// If `start_position` is provided, search starts from there instead of
    /// the bucket start.
    pub fn get_masked_from(
        &'_ self,
        hash: RandstrobeHash,
        hash_mask: RandstrobeHash,
        start_position: Option<usize>,
    ) -> Option<IndexEntry<'_>> {
        let masked_hash = hash & hash_mask;
        const MAX_LINEAR_SEARCH: usize = 4;
        let top_n = (hash >> (64 - self.bits)) as usize;
        let position_start = start_position.unwrap_or(self.bucket_starts[top_n] as usize);
        let position_end = self.bucket_starts[top_n + 1] as usize;
        let bucket = &self.randstrobes[position_start..position_end];
        if bucket.is_empty() {
            return None;
        } else if bucket.len() < MAX_LINEAR_SEARCH {
            for (pos, randstrobe) in bucket.iter().enumerate() {
                if randstrobe.hash() & hash_mask == masked_hash {
                    return Some(IndexEntry {
                        strobemer_index: self,
                        position: position_start + pos,
                    });
                }
                if randstrobe.hash() & hash_mask > masked_hash {
                    return None;
                }
            }
            return None;
        }

        let pos = custom_partition_point(bucket, |h| h.hash() & hash_mask < masked_hash);
        if pos < bucket.len() && bucket[pos].hash() & hash_mask == masked_hash {
            Some(IndexEntry {
                strobemer_index: self,
                position: position_start + pos,
            })
        } else {
            None
        }
    }

    pub fn get_masked(
        &'_ self,
        hash: RandstrobeHash,
        hash_mask: RandstrobeHash,
    ) -> Option<IndexEntry<'_>> {
        self.get_masked_from(hash, hash_mask, None)
    }

    /// Collects the entries of start..end whose reference position falls in
    /// one of intervals, which have to be sorted and disjoint.
    pub fn entries_in_intervals(
        &self,
        start: usize,
        end: usize,
        intervals: &[(usize, usize)],
        positions: &mut Vec<usize>,
    ) {
        let entries = &self.randstrobes[start..end];
        if entries.is_empty() || intervals.is_empty() {
            return;
        }

        const MIN_SEARCHABLE: usize = 64;
        if entries.len() < MIN_SEARCHABLE {
            for (i, entry) in entries.iter().enumerate() {
                if intervals_contain(intervals, entry.position()) {
                    positions.push(start + i);
                }
            }
            return;
        }

        let mut budget = entries.len() / 2;
        let mut block = 0usize;
        let cost = entries.len().ilog2() as usize + 1;
        while block < entries.len() {
            if budget < cost + intervals.len() {
                break;
            }
            let key = entries[block].hash_offset;
            let rest = &entries[block..];
            let block_end = block + custom_partition_point(rest, |h| h.hash_offset <= key);
            budget = budget.saturating_sub(cost);

            let mut i = block;
            for &(interval_start, interval_end) in intervals {
                if i >= block_end {
                    break;
                }
                let mut probes = 0;
                i = advance_to_position(&entries[..block_end], i, interval_start, &mut probes);
                budget = budget.saturating_sub(probes);
                while i < block_end && entries[i].position() < interval_end {
                    positions.push(start + i);
                    i += 1;
                }
            }
            block = block_end;
        }
        for (i, entry) in entries[block..].iter().enumerate() {
            if intervals_contain(intervals, entry.position()) {
                positions.push(start + block + i);
            }
        }
    }

    pub fn k(&self) -> usize {
        self.parameters.syncmer.k
    }
}

impl<'a> IndexEntry<'a> {
    fn randstrobe(&self) -> &RefRandstrobe {
        &self.strobemer_index.randstrobes[self.position]
    }

    pub fn hash(&self) -> RandstrobeHash {
        self.randstrobe().hash()
    }

    pub fn position(&self) -> usize {
        self.randstrobe().position()
    }

    pub fn strobe2_offset(&self) -> usize {
        self.randstrobe().strobe2_offset()
    }

    pub fn reference_index(&self) -> usize {
        self.randstrobe().reference_index()
    }

    pub fn get_hash_partial(&self) -> RandstrobeHash {
        self.strobemer_index.randstrobes[self.position].hash()
            & self.strobemer_index.parameters.randstrobe.main_hash_mask
    }

    pub fn get_hash_partial_forward(&self) -> RandstrobeHash {
        self.strobemer_index.randstrobes[self.position].hash_offset
            & self
                .strobemer_index
                .parameters
                .randstrobe
                .forward_main_hash_mask
    }

    pub fn strobe_extent_partial(&self) -> (usize, usize) {
        let p = self.strobemer_index.randstrobes[self.position].position;

        (p as usize, p as usize + self.k())
    }

    /// Count number of hits for the randstrobe *and* its "reverse complement"
    pub fn get_count_full(&self, hash_revcomp: u64) -> usize {
        let reverse_count;
        if let Some(entry_revcomp) = self.strobemer_index.get_full_forward(hash_revcomp) {
            reverse_count = entry_revcomp.get_count_full_forward();
        } else {
            reverse_count = 0;
        }
        reverse_count + self.get_count_full_forward()
    }

    pub fn get_count_full_forward(&self) -> usize {
        self.get_count(REF_RANDSTROBE_HASH_MASK)
    }

    pub fn get_count_partial(&self) -> usize {
        self.get_count(self.strobemer_index.parameters.randstrobe.main_hash_mask)
    }

    pub fn get_count(&self, hash_mask: u64) -> usize {
        let position = self.position;
        const MAX_LINEAR_SEARCH: usize = 8;
        let key = self.strobemer_index.randstrobes[position].hash();
        let masked_key = key & hash_mask;
        let top_n = (key >> (64 - self.strobemer_index.bits)) as usize;
        let position_end = self.strobemer_index.bucket_starts[top_n + 1] as usize;

        if position_end - position < MAX_LINEAR_SEARCH {
            let mut count = 1;
            for position_start in position + 1..position_end {
                if self.strobemer_index.randstrobes[position_start].hash() & hash_mask == masked_key
                {
                    count += 1;
                } else {
                    break;
                }
            }
            count
        } else {
            let bucket = &self.strobemer_index.randstrobes[position..position_end];
            custom_partition_point(bucket, |h| h.hash() & hash_mask <= masked_key)
        }
    }

    /// Return whether the randstrobe at the given position occurs more often than cutoff
    pub fn is_too_frequent_forward(&self, cutoff: usize) -> bool {
        if self.position + self.strobemer_index.filter_cutoff
            < self.strobemer_index.randstrobes.len()
        {
            self.strobemer_index.randstrobes[self.position].hash()
                == self.strobemer_index.randstrobes[self.position + cutoff].hash()
        } else {
            false
        }
    }

    pub fn is_too_frequent(&self, cutoff: usize, hash_revcomp: u64) -> bool {
        if self.is_too_frequent_forward(cutoff) {
            return true;
        }
        if let Some(entry_revcomp) = self.strobemer_index.get_full_forward(hash_revcomp) {
            if entry_revcomp.is_too_frequent_forward(cutoff) {
                return true;
            }
            let count = self.get_count_full_forward() + entry_revcomp.get_count_full_forward();

            return count > cutoff;
        }

        false
    }

    pub fn is_too_frequent_forward_partial(&self, cutoff: usize) -> bool {
        if self.position + cutoff < self.strobemer_index.randstrobes.len() {
            self.strobemer_index.randstrobes[self.position].hash()
                & self.strobemer_index.parameters.randstrobe.main_hash_mask
                == self.strobemer_index.randstrobes[self.position + cutoff].hash()
                    & self.strobemer_index.parameters.randstrobe.main_hash_mask
        } else {
            false
        }
    }

    pub fn is_too_frequent_partial(&self, cutoff: usize) -> bool {
        self.is_too_frequent_forward_partial(cutoff)
    }

    fn k(&self) -> usize {
        self.strobemer_index.k()
    }
}

const STI_FILE_FORMAT_VERSION: u32 = 7;

#[derive(Error, Debug)]
pub enum IndexReadingError {
    #[error("When reading the .sti (index) file: {0}")]
    Io(#[from] std::io::Error),

    #[error("Signature at the beginning of the .sti file is incorrect.")]
    WrongMagic,

    #[error("The .sti file format version is {0}, but this version of strobealign can only read files that have version {ver}", ver = STI_FILE_FORMAT_VERSION)]
    WrongFileFormatVersion(u32),

    #[error("Index parameters in .sti file and those specified on command line differ")]
    ParameterMismatch,

    #[error(
        "The randstrobe starts vector has an unexpected size in the .sti file. Did you use the correct -b option?"
    )]
    RandstrobeStartIndicesWrongSize,

    #[error("The .sti (index) file uses an unknown profile")]
    WrongProfile,

    #[error("The .sti (index) file uses an invalid indexing parameter: {0}")]
    InvalidIndexParameter(#[from] InvalidSeedingParameter),
}

impl StrobemerIndex {
    pub fn write<P: AsRef<Path>>(&self, path: P) -> Result<(), Error> {
        let mut file = File::create(path)?;

        file.write_all(b"STI\x01")?; // Magic number
        file.write_all(&STI_FILE_FORMAT_VERSION.to_ne_bytes())?;

        // Variable-length chunk reserved for future use
        file.write_all(&8u64.to_ne_bytes())?; // length in bytes
        file.write_all(&0u64.to_ne_bytes())?; // contents

        file.write_all(&(self.filter_cutoff as u32).to_ne_bytes())?;
        file.write_all(&(self.bits as u32).to_ne_bytes())?;

        file.write_all(&(self.parameters.profile as u32).to_ne_bytes())?;
        let sp = &self.parameters.syncmer;
        for val in [sp.k, sp.s].iter() {
            file.write_all(&(*val as u32).to_ne_bytes())?
        }
        let rp = &self.parameters.randstrobe;
        for val in [
            rp.w_min as u32,
            rp.w_max as u32,
            rp.q as u32,
            rp.max_dist as u32,
        ]
        .iter()
        {
            file.write_all(&val.to_ne_bytes())?;
        }
        file.write_all(&rp.main_hash_mask.to_ne_bytes())?;

        write_vec(&mut file, &self.randstrobes)?;
        write_vec(&mut file, &self.bucket_starts)?;

        Ok(())
    }
}

/// Returns the first entry at or after from whose reference position reaches
/// the given position.
fn advance_to_position(
    entries: &[RefRandstrobe],
    from: usize,
    position: usize,
    probes: &mut usize,
) -> usize {
    if from >= entries.len() || entries[from].position() >= position {
        return from;
    }
    let mut base = from;
    let mut step = 1usize;
    while base + step < entries.len() && entries[base + step].position() < position {
        base += step;
        step *= 2;
        *probes += 1;
    }
    // `entries[base]` is below `position` and `entries[base + step]` is not (or
    // is past the end), so the answer lies in `(base, base + step]`.
    let end = (base + step + 1).min(entries.len());
    *probes += (end - base).ilog2() as usize + 1;

    base + custom_partition_point(&entries[base..end], |h| h.position() < position)
}

/// Returns whether a reference position falls in one of the sorted disjoint
/// intervals.
fn intervals_contain(intervals: &[(usize, usize)], position: usize) -> bool {
    if intervals.len() < 8 {
        return intervals
            .iter()
            .any(|&(start, end)| position >= start && position < end);
    }
    let i = intervals.partition_point(|&(start, _)| start <= position);

    i > 0 && position < intervals[i - 1].1
}

pub fn read_index<P: AsRef<Path>>(
    path: P,
    parameters: SeedingParameters,
    bits: u8,
) -> Result<StrobemerIndex, IndexReadingError> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);

    let mut magic = [0; 4];
    reader.read_exact(&mut magic)?;
    if &magic != b"STI\x01" {
        // Magic number
        return Err(IndexReadingError::WrongMagic);
    }

    let file_format_version = read_u32(&mut reader)?;
    if file_format_version != STI_FILE_FORMAT_VERSION {
        return Err(IndexReadingError::WrongFileFormatVersion(
            file_format_version,
        ));
    }

    // Skip over variable-length chunk reserved for future use
    let mut buf = [0; 8];
    reader.read_exact(&mut buf)?;
    let reserved_chunk_size = u64::from_ne_bytes(buf);
    let mut buf = vec![0; reserved_chunk_size as usize];
    reader.read_exact(&mut buf)?;

    let filter_cutoff = read_u32(&mut reader)? as usize;
    let sti_bits = read_u32(&mut reader)? as u8;
    if sti_bits != bits {
        return Err(IndexReadingError::ParameterMismatch);
    }

    let profile = read_u32(&mut reader)?
        .try_into()
        .map_err(|_| IndexReadingError::WrongProfile)?;

    let k = read_u32(&mut reader)? as usize;
    let s = read_u32(&mut reader)? as usize;

    let w_min = read_u32(&mut reader)? as usize;
    let w_max = read_u32(&mut reader)? as usize;
    let q = read_u32(&mut reader)? as u64;
    let max_dist = read_u32(&mut reader)? as u8;
    let main_hash_mask = read_u64(&mut reader)?;

    let syncmer_parameters = SyncmerParameters::try_new(k, s)?;
    let partial_orientation_pos = main_hash_mask.trailing_zeros() - 1;
    let randstrobe_parameters = RandstrobeParameters {
        w_min,
        w_max,
        q,
        max_dist,
        main_hash_mask,
        forward_main_hash_mask: main_hash_mask | (1u64 << partial_orientation_pos),
        partial_orientation_pos,
    };
    let sti_parameters = SeedingParameters {
        profile,
        syncmer: syncmer_parameters,
        randstrobe: randstrobe_parameters,
    };

    if parameters != sti_parameters {
        return Err(IndexReadingError::ParameterMismatch);
    }

    let randstrobes = read_vec(&mut reader)?;
    let bucket_starts = read_vec(&mut reader)?;

    if bucket_starts.len() != (1 << bits) + 1 {
        return Err(IndexReadingError::RandstrobeStartIndicesWrongSize);
    }

    Ok(StrobemerIndex {
        parameters,
        bits,
        filter_cutoff,
        randstrobes,
        bucket_starts,
    })
}

fn read_u32<T: BufRead>(file: &mut T) -> Result<u32, Error> {
    let mut buf = [0; 4];
    file.read_exact(&mut buf)?;

    Ok(u32::from_ne_bytes(buf))
}

fn read_u64<T: BufRead>(file: &mut T) -> Result<u64, Error> {
    let mut buf = [0; 8];
    file.read_exact(&mut buf)?;

    Ok(u64::from_ne_bytes(buf))
}

fn read_vec<T, R: Read>(file: &mut BufReader<R>) -> Result<Vec<T>, Error> {
    let length = read_u64(file)? as usize;
    let bytes = vec![0; length * size_of::<T>()];
    let mut bytes = std::mem::ManuallyDrop::new(bytes);
    file.read_exact(&mut bytes)?;

    let data = unsafe {
        let (pointer, length, capacity) = (bytes.as_mut_ptr(), bytes.len(), bytes.capacity());
        let pointer = pointer.cast::<T>();
        let length = length / size_of::<T>();
        let capacity = capacity / size_of::<T>();

        Vec::from_raw_parts(pointer, length, capacity)
    };

    Ok(data)
}

fn write_vec<T>(file: &mut File, data: &[T]) -> Result<(), Error> {
    // write length
    file.write_all(&(data.len().to_ne_bytes()))?;
    // write data
    let data: &[u8] = unsafe {
        std::slice::from_raw_parts(data.as_ptr() as *const u8, std::mem::size_of_val(data))
    };

    file.write_all(data)
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::indexer::make_index;
    use crate::io::fasta::read_ref;
    use crate::packed_seq::PackedSeq;
    use crate::refseq::RefSequence;
    use crate::revcomp::reverse_complement;

    #[test]
    fn ref_randstrobe() {
        let hash: u64 = 0x1234567890ABCDEFu64 & REF_RANDSTROBE_HASH_MASK;
        let ref_index: u32 = (REF_RANDSTROBE_MAX_NUMBER_OF_REFERENCES - 1) as u32;
        let offset = 255;
        let position = !0;
        let rr = RefRandstrobe::new(hash, ref_index, position, offset);

        assert_eq!(rr.hash(), hash);
        assert_eq!(rr.position(), position as usize);
        assert_eq!(rr.reference_index(), ref_index as usize);
        assert_eq!(rr.strobe2_offset(), offset as usize);
    }

    #[test]
    fn sti_parameters_mismatch() {
        use temp_dir::TempDir;

        let dir = TempDir::new().unwrap();
        let fasta_path = dir.path().join("phix.fasta");
        std::fs::copy("tests/phix.fasta", &fasta_path).unwrap();
        let references = read_ref(fasta_path).unwrap();

        let parameters = SeedingParameters::new(300);
        let bits = parameters.syncmer.pick_bits(&references);
        let index = make_index(&references, parameters, bits, 0.0002, 1).0;
        let sti_path = dir.path().join("index.sti");
        index.write(&sti_path).unwrap();

        let read_index_result = read_index(&sti_path, SeedingParameters::new(50), bits);

        match read_index_result {
            Err(IndexReadingError::ParameterMismatch) => {}
            _ => {
                panic!("Parameters are expected not to match");
            }
        }
    }

    #[test]
    fn partial_orientation() {
        let refseq = read_ref("tests/phix.fasta").unwrap();
        let seq_decoded = refseq.contig(0).decode_all();
        let rc_seq = reverse_complement(&seq_decoded);
        let mut rc_refseq = RefSequence::new();
        rc_refseq.push("phix_rc".to_string(), PackedSeq::from_slice(&rc_seq));

        let parameters = SeedingParameters::new(300);
        let bits = parameters.syncmer.pick_bits(&refseq);

        let (fwd_index, _stats) = make_index(&refseq, parameters.clone(), bits, 0.0000001, 1);

        let (rc_index, _stats) = make_index(&rc_refseq, parameters.clone(), bits, 0.0000001, 1);

        assert_eq!(fwd_index.randstrobes.len(), rc_index.randstrobes.len());

        // Iterate over fwd index entries and look up each hash in the rc index
        for i in 0..fwd_index.randstrobes.len() {
            let fwd_hash = fwd_index.randstrobes[i].hash();

            let rev_pos_result = rc_index.get_partial(fwd_hash);
            assert!(rev_pos_result.is_some());
            let rev_entry = rev_pos_result.unwrap();

            let rc_partial_query = fwd_index.randstrobes[i].hash()
                ^ (1u64 << rc_index.parameters.randstrobe.partial_orientation_pos);
            let rev_pos_forward =
                rc_index.get_partial_forward_from(rc_partial_query, rev_entry.position);
            assert!(rev_pos_forward.is_some());
        }
    }

    #[test]
    fn an_entry_inside_an_interval_is_found() {
        let hash = 0x0123_4567_89AB_CD00u64 & REF_RANDSTROBE_HASH_MASK;
        let randstrobes = [10u32, 20, 30, 40]
            .map(|position| RefRandstrobe::new(hash, 0, position, 0))
            .to_vec();
        let index = StrobemerIndex::new(
            SeedingParameters::new(300),
            1,
            100,
            randstrobes,
            vec![0, 4, 4],
        );
        let mut found = vec![];

        index.entries_in_intervals(0, 4, &[(15, 35)], &mut found);

        assert_eq!(found, vec![1, 2]);
    }

    #[test]
    fn an_entry_on_the_end_of_an_interval_is_left_out() {
        let hash = 0x0123_4567_89AB_CD00u64 & REF_RANDSTROBE_HASH_MASK;
        let randstrobes = [10u32, 20, 30]
            .map(|position| RefRandstrobe::new(hash, 0, position, 0))
            .to_vec();
        let index = StrobemerIndex::new(
            SeedingParameters::new(300),
            1,
            100,
            randstrobes,
            vec![0, 3, 3],
        );
        let mut found = vec![];

        index.entries_in_intervals(0, 3, &[(10, 30)], &mut found);

        assert_eq!(found, vec![0, 1]);
    }

    #[test]
    fn only_the_entries_of_the_given_run_are_searched() {
        let hash = 0x0123_4567_89AB_CD00u64 & REF_RANDSTROBE_HASH_MASK;
        let randstrobes = [10u32, 20, 30, 40]
            .map(|position| RefRandstrobe::new(hash, 0, position, 0))
            .to_vec();
        let index = StrobemerIndex::new(
            SeedingParameters::new(300),
            1,
            100,
            randstrobes,
            vec![0, 4, 4],
        );
        let mut found = vec![];

        index.entries_in_intervals(2, 4, &[(0, 100)], &mut found);

        assert_eq!(found, vec![2, 3]);
    }

    #[test]
    fn a_run_past_the_scan_threshold_is_searched() {
        let hash = 0x0123_4567_89AB_CD00u64 & REF_RANDSTROBE_HASH_MASK;
        let randstrobes: Vec<RefRandstrobe> = (0..80)
            .map(|i| RefRandstrobe::new(hash, 0, i * 3, 0))
            .collect();
        let index = StrobemerIndex::new(
            SeedingParameters::new(300),
            1,
            100,
            randstrobes,
            vec![0, 80, 80],
        );
        let mut found = vec![];

        index.entries_in_intervals(0, 80, &[(30, 45)], &mut found);

        assert_eq!(found, vec![10, 11, 12, 13, 14]);
    }

    #[test]
    fn entries_are_found_in_every_offset_block_of_a_run() {
        let hash = 0x0123_4567_89AB_CD00u64 & REF_RANDSTROBE_HASH_MASK;
        let randstrobes: Vec<RefRandstrobe> = [0u8, 1, 2, 3]
            .into_iter()
            .flat_map(|offset| (0..20).map(move |i| RefRandstrobe::new(hash, 0, i * 2, offset)))
            .collect();
        let index = StrobemerIndex::new(
            SeedingParameters::new(300),
            1,
            100,
            randstrobes,
            vec![0, 80, 80],
        );
        let mut found = vec![];

        index.entries_in_intervals(0, 80, &[(10, 20)], &mut found);

        #[rustfmt::skip]
        assert_eq!(found, vec![
             5,  6,  7,  8,  9,
            25, 26, 27, 28, 29,
            45, 46, 47, 48, 49,
            65, 66, 67, 68, 69,
        ]);
    }

    #[test]
    fn a_run_with_more_intervals_than_budget_loses_no_entry() {
        let hash = 0x0123_4567_89AB_CD00u64 & REF_RANDSTROBE_HASH_MASK;
        let randstrobes: Vec<RefRandstrobe> = (0..64)
            .map(|i| RefRandstrobe::new(hash, 0, i * 10, 0))
            .collect();
        let index = StrobemerIndex::new(
            SeedingParameters::new(300),
            1,
            100,
            randstrobes,
            vec![0, 64, 64],
        );
        #[rustfmt::skip]
        let intervals = [
            (5, 15), (25, 35), (45, 55), (65, 75), (85, 95),
            (105, 115), (125, 135), (145, 155), (165, 175), (185, 195),
        ];
        let mut found = vec![];

        index.entries_in_intervals(0, 64, &intervals, &mut found);

        assert_eq!(found, vec![1, 3, 5, 7, 9, 11, 13, 15, 17, 19]);
    }

    #[test]
    fn advancing_lands_on_the_first_entry_that_reaches_a_position() {
        let hash = 0x0123_4567_89AB_CD00u64 & REF_RANDSTROBE_HASH_MASK;
        let entries =
            [0u32, 5, 5, 12, 40, 41, 90].map(|position| RefRandstrobe::new(hash, 0, position, 0));
        let mut probes = 0;

        assert_eq!(advance_to_position(&entries, 0, 0, &mut probes), 0);
        assert_eq!(advance_to_position(&entries, 0, 5, &mut probes), 1);
        assert_eq!(advance_to_position(&entries, 0, 6, &mut probes), 3);
        assert_eq!(advance_to_position(&entries, 0, 41, &mut probes), 5);
        assert_eq!(advance_to_position(&entries, 0, 91, &mut probes), 7);
        assert_eq!(advance_to_position(&entries, 4, 5, &mut probes), 4);
        assert_eq!(advance_to_position(&entries, 7, 0, &mut probes), 7);
    }

    #[test]
    fn a_position_belongs_to_an_interval_up_to_its_end() {
        let few = [(10usize, 20), (30, 40), (50, 60)];

        assert!(intervals_contain(&few, 10));
        assert!(!intervals_contain(&few, 20));
        assert!(!intervals_contain(&few, 45));
        assert!(intervals_contain(&few, 59));

        #[rustfmt::skip]
        let many = [
            (0usize, 2), (4, 6), (8, 10), (12, 14),
            (16, 18), (20, 22), (24, 26), (28, 30),
        ];

        assert!(intervals_contain(&many, 0));
        assert!(!intervals_contain(&many, 3));
        assert!(intervals_contain(&many, 28));
        assert!(!intervals_contain(&many, 30));
    }
}
