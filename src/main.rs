use std::cmp::min;
use std::collections::VecDeque;
use std::{env, io};
use std::collections::hash_map::RandomState;
use std::fs::File;
use std::io::{BufReader, BufWriter, BufRead, Error, Write};
use std::path::Path;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(long_about = None)]
struct Args {

   /// Print syncmers instead of randstrobes
   #[arg(long, default_value_t = false)]
   syncmers: bool,

   /// Path to input FASTA
   path: String,
}

#[derive(Debug)]
struct FastaRecord {
    name: String,
    sequence: Vec<u8>,
}

fn read_fasta<R: BufRead>(reader: &mut R) -> Result<Vec<FastaRecord>, Error> {
    let mut records = Vec::<FastaRecord>::new();
    let mut name = String::new();
    let mut sequence = Vec::new();
    let mut has_record = false;
    for line in reader.lines() {
        let line = line.unwrap();
        let line = line.as_bytes();
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            if has_record {
                records.push(FastaRecord{name, sequence});
            }
            name = String::from_utf8(line[1..].to_vec()).unwrap();
            if let Some(i) = name.find(|c: char| c.is_ascii_whitespace()) {
                name = name[..i].to_string();
            }
            sequence = Vec::new();
            has_record = true;
        } else {
            sequence.extend(line);
        }
    }
    if has_record {
        records.push(FastaRecord{name, sequence});
    }

    Ok(records)
}


#[derive(Debug, Clone, Copy)]
struct Syncmer {
    hash: u64,
    position: usize,
}

struct SyncmerIterator<'a> {
    seq: &'a [u8],
    k: usize,
    s: usize,
    t: usize,
    kmask: u64,
    smask: u64,
    kshift: usize,
    sshift: usize,
    qs: VecDeque<u64>, // s-mer hashes
    qs_min_val: u64,
    qs_min_pos: usize,
    l: usize,
    xk: [u64; 2],
    xs: [u64; 2],
    i: usize,
}

impl<'a> SyncmerIterator<'a> {
    fn new(seq: &'a [u8], k: usize, s: usize, t: usize) -> SyncmerIterator {
        SyncmerIterator {
            seq,
            k,
            s,
            t,
            kmask: (1 << 2*k) - 1,
            smask: (1 << 2*s) - 1,
            kshift: (k - 1) * 2,
            sshift: (s - 1) * 2,
            qs: VecDeque::new(),
            qs_min_val: u64::MAX,
            qs_min_pos: 0,
            l: 0,
            xk: [0, 0],
            xs: [0, 0],
            i: 0,
        }
    }
}

// a, A -> 0
// c, C -> 1
// g, G -> 2
// t, T, u, U -> 3
static NUCLEOTIDES: [u8; 256] = [
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
];

fn syncmer_kmer_hash(hash: u64) -> u64 {
    hash
}

impl<'a> Iterator for SyncmerIterator<'a> {
    type Item = Syncmer;

    fn next(&mut self) -> Option<Self::Item> {



        while self.i < self.seq.len() {
            let ch = self.seq[self.i];
            let c = NUCLEOTIDES[ch as usize];
            if c < 4 { // not an "N" base
                self.xk[0] = (self.xk[0] << 2 | (c as u64)) & self.kmask;                  // forward strand
                self.xk[1] = self.xk[1] >> 2 | ((3 - c) as u64) << self.kshift;  // reverse strand
                self.xs[0] = (self.xs[0] << 2 | (c as u64)) & self.smask;                  // forward strand
                self.xs[1] = self.xs[1] >> 2 | ((3 - c) as u64) << self.sshift;  // reverse strand
                self.l += 1;
                if self.l < self.s {
                    self.i += 1;
                    continue;
                }
                // we find an s-mer
                let ys = min(self.xs[0], self.xs[1]);
                let hash_s = ys;
                self.qs.push_back(hash_s);
                // not enough hashes in the queue, yet
                if self.qs.len() < self.k - self.s + 1 {
                    self.i += 1;
                    continue;
                }
                if self.qs.len() == (self.k - self.s + 1) { // We are at the last s-mer within the first k-mer, need to decide if we add it
                    // TODO use argmin
                    for j in 0..self.qs.len() {
                        if self.qs[j] < self.qs_min_val {
                            self.qs_min_val = self.qs[j];
                            self.qs_min_pos = self.i + j + 1 - self.k;
                        }
                    }
                } else {
                    // update queue and current minimum and position
                    self.qs.pop_front();
                    if self.qs_min_pos == self.i - self.k { // we popped the previous minimizer, find new brute force
                        self.qs_min_val = u64::MAX;
                        self.qs_min_pos = self.i - self.s + 1;
                        for j in (0..self.qs.len()).rev() { //Iterate in reverse to choose the rightmost minimizer in a window
                            if self.qs[j] < self.qs_min_val {
                                self.qs_min_val = self.qs[j];
                                self.qs_min_pos = self.i + j + 1 - self.k;
                            }
                        }
                    } else if hash_s < self.qs_min_val { // the new value added to queue is the new minimum
                        self.qs_min_val = hash_s;
                        self.qs_min_pos = self.i - self.s + 1;
                    }
                }
                if self.qs_min_pos == self.i + self.t - self.k { // occurs at t:th position in k-mer
                    let yk = min(self.xk[0], self.xk[1]);
                    let syncmer = Syncmer { hash: syncmer_kmer_hash(yk), position: self.i - self.k + 1 };
                    self.i += 1;
                    return Some(syncmer);
                }
            } else {
                // if there is an "N", restart
                self.qs_min_val = u64::MAX;
                self.qs_min_pos = 0;
                self.l = 0;
                self.xs = [0, 0];
                self.xk = [0, 0];
                self.qs.clear();
            }
            self.i += 1;
        }
        None
    }
}

#[derive(Debug, PartialEq, Eq)]
struct Randstrobe {
    hash: u64,
    strobe1_pos: usize,
    strobe2_pos: usize,
}

struct RandstrobeParameters {
    w_min: usize,
    w_max: usize,
    q: u64,
    max_dist: usize,
}

struct RandstrobeIterator<'a> {
    parameters: &'a RandstrobeParameters,
    syncmers: VecDeque<Syncmer>,
    syncmer_iterator: &'a mut SyncmerIterator<'a>,
    strobe1_index: usize,
}

impl<'a> RandstrobeIterator<'a> {
    pub fn new(syncmer_iterator: &'a mut SyncmerIterator<'a>, parameters: &'a RandstrobeParameters) -> RandstrobeIterator<'a> {
        RandstrobeIterator {
            parameters,
            syncmers: VecDeque::<Syncmer>::new(),
            syncmer_iterator,
            strobe1_index: 0,
        }
    }
}

impl<'a> Iterator for RandstrobeIterator<'a> {
    type Item = Randstrobe;
    fn next(&mut self) -> Option<Self::Item> {
        while self.syncmers.len() <= self.parameters.w_max {
            if let Some(syncmer) = self.syncmer_iterator.next() {
                self.syncmers.push_back(syncmer);
            } else {
                break;
            }
        }
        if self.syncmers.len() <= self.parameters.w_min {
            return None;
        }
        let strobe1 = self.syncmers[0];
        let max_position = strobe1.position + self.parameters.max_dist;
        let mut min_val = u64::MAX;
        let mut strobe2 = self.syncmers[0]; // Defaults if no nearby syncmer

        for i in self.parameters.w_min .. self.syncmers.len() {
            debug_assert!(i <= self.parameters.w_max);
            if self.syncmers[i].position > max_position {
                break;
            }
            let b = (strobe1.hash ^ self.syncmers[i].hash) & self.parameters.q;
            let ones = b.count_ones() as u64;
            if ones < min_val {
                min_val = ones;
                strobe2 = self.syncmers[i];
            }
        }
        self.syncmers.pop_front();
        return Some(Randstrobe {
            hash: strobe1.hash + strobe2.hash,
            strobe1_pos: strobe1.position,
            strobe2_pos: strobe2.position,
        })
    }
}

fn main() -> Result<(), Error> {
    let args = Args::parse();
    let path = Path::new(&args.path);
    let f = File::open(path)?;
    let mut reader = BufReader::new(f);
    let records = read_fasta(&mut reader).unwrap();
    let k = 20usize;
    let s = 16usize;

    let mut writer = BufWriter::new(io::stdout());

    // IndexParameters(r=150, k=20, s=16, l=1, u=7, q=255, max_dist=80, t_syncmer=3, w_min=5, w_max=11)
    let parameters = RandstrobeParameters {
        w_min: 5,
        w_max: 11,
        q: 255,
        max_dist: 80,
    };
    if args.syncmers {
        for record in &records {
            let name = &record.name;
            for syncmer in SyncmerIterator::new(&record.sequence, k, s, 3) {
                writeln!(writer, "{}\t{}\t{}", name, syncmer.position, syncmer.position + k)?;
            }
        }
    } else {
        for record in &records {
            let name = &record.name;
            let mut syncmers = SyncmerIterator::new(&record.sequence, k, s, 3);
            for randstrobe in RandstrobeIterator::new(&mut syncmers, &parameters) {
                writeln!(writer, "{}\t{}\t{}", name, randstrobe.strobe1_pos, randstrobe.strobe2_pos + k)?;
            }
        }
    }
    Ok(())
}


#[test]
fn verify_cli() {
    use clap::CommandFactory;
    Args::command().debug_assert()
}
