use std::fs::File;
use std::io::{BufRead, BufReader, Error, Read, Write};
use std::path::Path;

use thiserror::Error;

use crate::io::fasta::RefSequence;
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

pub struct StrobemerIndex<'a> {
    pub references: &'a [RefSequence],
    pub parameters: SeedingParameters,

    /// no. of bits of the hash to use when indexing a randstrobe bucket
    pub bits: u8,

    /// Regular (non-rescue) NAM finding ignores randstrobes that occur more often than
    /// this (see StrobemerIndex::is_too_frequent())
    pub filter_cutoff: usize,

    /// Filter partial seeds that occur more often than this
    pub partial_filter_cutoff: usize,

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
    pub randstrobe_start_indices: Vec<BucketIndex>,
}

impl<'a> StrobemerIndex<'a> {
    // Find the first entry that matches the forwald full hash (including orientation bits)
    pub fn get_full_forward(&self, hash: RandstrobeHash) -> Option<usize> {
        self.get_masked(hash, REF_RANDSTROBE_HASH_MASK)
    }

    /// Find the first entry that matches the undirected main hash (without orientation bit)
    pub fn get_partial(&self, hash: RandstrobeHash) -> Option<usize> {
        self.get_masked(hash, self.parameters.randstrobe.main_hash_mask)
    }

    /// Find the first entry matching the forward main hash
    pub fn get_partial_forward(&self, hash: RandstrobeHash) -> Option<usize> {
        self.get_masked(hash, self.parameters.randstrobe.forward_main_hash_mask)
    }

    /// Find the first entry matching the forward main hash, starting from
    /// the undirected main position
    pub fn get_partial_forward_from(
        &self,
        hash: RandstrobeHash,
        undirected_position: usize,
    ) -> Option<usize> {
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
        &self,
        hash: RandstrobeHash,
        hash_mask: RandstrobeHash,
        start_position: Option<usize>,
    ) -> Option<usize> {
        let masked_hash = hash & hash_mask;
        const MAX_LINEAR_SEARCH: usize = 4;
        let top_n = (hash >> (64 - self.bits)) as usize;
        let position_start =
            start_position.unwrap_or(self.randstrobe_start_indices[top_n] as usize);
        let position_end = self.randstrobe_start_indices[top_n + 1];
        let bucket = &self.randstrobes[position_start as usize..position_end as usize];
        if bucket.is_empty() {
            return None;
        } else if bucket.len() < MAX_LINEAR_SEARCH {
            for (pos, randstrobe) in bucket.iter().enumerate() {
                if randstrobe.hash() & hash_mask == masked_hash {
                    return Some(position_start as usize + pos);
                }
                if randstrobe.hash() & hash_mask > masked_hash {
                    return None;
                }
            }
            return None;
        }

        let pos = custom_partition_point(bucket, |h| h.hash() & hash_mask < masked_hash);
        if pos < bucket.len() && bucket[pos].hash() & hash_mask == masked_hash {
            Some(position_start as usize + pos)
        } else {
            None
        }
    }

    pub fn get_masked(&self, hash: RandstrobeHash, hash_mask: RandstrobeHash) -> Option<usize> {
        self.get_masked_from(hash, hash_mask, None)
    }

    pub fn k(&self) -> usize {
        self.parameters.syncmer.k
    }

    pub fn get_hash_partial(&self, position: usize) -> RandstrobeHash {
        self.randstrobes[position].hash() & self.parameters.randstrobe.main_hash_mask
    }

    pub fn get_hash_partial_forward(&self, position: usize) -> RandstrobeHash {
        self.randstrobes[position].hash_offset & self.parameters.randstrobe.forward_main_hash_mask
    }

    pub fn strobe_extent_partial(&self, position: usize) -> (usize, usize) {
        let p = self.randstrobes[position].position;

        (p as usize, p as usize + self.k())
    }

    /// Count number of hits for the randstrobe *and* its "reverse complement"
    pub fn get_count_full(&self, position: usize, hash_revcomp: u64) -> usize {
        let reverse_count;
        if let Some(position_revcomp) = self.get_full_forward(hash_revcomp) {
            reverse_count = self.get_count_full_forward(position_revcomp);
        } else {
            reverse_count = 0;
        }
        reverse_count + self.get_count_full_forward(position)
    }

    pub fn get_count_full_forward(&self, position: usize) -> usize {
        self.get_count(position, REF_RANDSTROBE_HASH_MASK)
    }

    pub fn get_count_partial(&self, position: usize) -> usize {
        self.get_count(position, self.parameters.randstrobe.main_hash_mask)
    }

    pub fn get_count(&self, position: usize, hash_mask: u64) -> usize {
        const MAX_LINEAR_SEARCH: usize = 8;
        let key = self.randstrobes[position].hash();
        let masked_key = key & hash_mask;
        let top_n = (key >> (64 - self.bits)) as usize;
        let position_end = self.randstrobe_start_indices[top_n + 1] as usize;

        if position_end - position < MAX_LINEAR_SEARCH {
            let mut count = 1;
            for position_start in position + 1..position_end {
                if self.randstrobes[position_start].hash() & hash_mask == masked_key {
                    count += 1;
                } else {
                    break;
                }
            }
            count
        } else {
            let bucket = &self.randstrobes[position..position_end];
            custom_partition_point(bucket, |h| h.hash() & hash_mask <= masked_key)
        }
    }

    /// Return whether the randstrobe at the given position occurs more often than cutoff
    pub fn is_too_frequent_forward(&self, position: usize, cutoff: usize) -> bool {
        if position + self.filter_cutoff < self.randstrobes.len() {
            self.randstrobes[position].hash() == self.randstrobes[position + cutoff].hash()
        } else {
            false
        }
    }

    pub fn is_too_frequent(&self, position: usize, cutoff: usize, hash_revcomp: u64) -> bool {
        if self.is_too_frequent_forward(position, cutoff) {
            return true;
        }
        if let Some(position_revcomp) = self.get_full_forward(hash_revcomp) {
            if self.is_too_frequent_forward(position_revcomp, cutoff) {
                return true;
            }
            let count = self.get_count_full_forward(position)
                + self.get_count_full_forward(position_revcomp);

            return count > cutoff;
        }

        false
    }

    pub fn is_too_frequent_forward_partial(&self, position: usize, cutoff: usize) -> bool {
        if position + cutoff < self.randstrobes.len() {
            self.randstrobes[position].hash() & self.parameters.randstrobe.main_hash_mask
                == self.randstrobes[position + cutoff].hash()
                    & self.parameters.randstrobe.main_hash_mask
        } else {
            false
        }
    }

    pub fn is_too_frequent_partial(&self, position: usize, cutoff: usize) -> bool {
        self.is_too_frequent_forward_partial(position, cutoff)
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

impl<'a> StrobemerIndex<'a> {
    pub fn write<P: AsRef<Path>>(&self, path: P) -> Result<(), Error> {
        let mut file = File::create(path)?;

        file.write_all(b"STI\x01")?; // Magic number
        file.write_all(&STI_FILE_FORMAT_VERSION.to_ne_bytes())?;

        // Variable-length chunk reserved for future use
        file.write_all(&8u64.to_ne_bytes())?; // length in bytes
        file.write_all(&0u64.to_ne_bytes())?; // contents

        file.write_all(&(self.filter_cutoff as u32).to_ne_bytes())?;
        file.write_all(&(self.bits as u32).to_ne_bytes())?;

        file.write_all(&(self.parameters.profile.clone() as u32).to_ne_bytes())?;
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
        file.write_all(&(rp.main_hash_mask as u64).to_ne_bytes())?;

        write_vec(&mut file, &self.randstrobes)?;
        write_vec(&mut file, &self.randstrobe_start_indices)?;

        Ok(())
    }
}

pub fn read_index<'a, P: AsRef<Path>>(
    path: P,
    references: &'a [RefSequence],
    parameters: SeedingParameters,
    bits: Option<u8>,
) -> Result<StrobemerIndex<'a>, IndexReadingError> {
    // TODO these two lines are duplicated in make_index
    let total_reference_length = references.iter().map(|r| r.sequence.len()).sum();
    let bits = bits.unwrap_or_else(|| parameters.syncmer.pick_bits(total_reference_length));

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
    let rescue_cutoff = filter_cutoff;
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
    let randstrobe_start_indices = read_vec(&mut reader)?;

    if randstrobe_start_indices.len() != (1 << bits) + 1 {
        return Err(IndexReadingError::RandstrobeStartIndicesWrongSize);
    }

    Ok(StrobemerIndex {
        references,
        parameters,
        bits,
        filter_cutoff,
        partial_filter_cutoff: filter_cutoff,
        rescue_cutoff,
        randstrobes,
        randstrobe_start_indices,
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
        std::slice::from_raw_parts(data.as_ptr() as *const u8, data.len() * size_of::<T>())
    };

    file.write_all(data)
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::indexer::make_index;
    use crate::io::fasta::read_fasta;

    use std::fs::File;
    use std::io::BufReader;

    #[test]
    fn test_ref_randstrobe() {
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
    fn test_sti_parameters_mismatch() {
        use temp_dir::TempDir;

        let dir = TempDir::new().unwrap();
        let fasta_path = dir.path().join("phix.fasta");
        std::fs::copy("tests/phix.fasta", &fasta_path).unwrap();
        let f = File::open(fasta_path).unwrap();
        let references = read_fasta(&mut BufReader::new(f)).unwrap();

        let parameters = SeedingParameters::new(300);
        let index = make_index(&references, parameters, None, 0.0002, 1).0;
        let sti_path = dir.path().join("index.sti");
        index.write(&sti_path).unwrap();

        let read_index_result =
            read_index(&sti_path, &references, SeedingParameters::new(50), None);

        match read_index_result {
            Err(IndexReadingError::ParameterMismatch) => {}
            _ => {
                panic!("Parameters are expected not to match");
            }
        }
    }
}
