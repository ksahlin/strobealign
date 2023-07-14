use std::cmp::min;
use std::collections::VecDeque;
use std::{env, io};
use std::collections::hash_map::RandomState;
use std::fs::File;
use std::hash::Hasher;
use std::io::{BufReader, BufWriter, BufRead, Error, Write};
use std::path::Path;
use clap::Parser;
use fxhash;
use rstrobes::fastq::FastqReader;
use rstrobes::strobes::{SyncmerIterator,RandstrobeIterator,RandstrobeParameters};

#[derive(Parser, Debug)]
#[command(long_about = None)]
struct Args {

    /// Print syncmers instead of randstrobes
    //#[arg(long, default_value_t = false)]
    //syncmers: bool,

    /// Path to input FASTA
    ref_path: String,

    fastq_path: String,
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

fn main() -> Result<(), Error> {
    let args = Args::parse();
    let path = Path::new(&args.ref_path);
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
    for record in &records {
        let name = &record.name;
        let mut syncmers = SyncmerIterator::new(&record.sequence, k, s, 3);
        for randstrobe in RandstrobeIterator::new(&mut syncmers, &parameters) {
            writeln!(writer, "{}\t{}\t{}", name, randstrobe.strobe1_pos, randstrobe.strobe2_pos + k)?;
        }
    }
    Ok(())
}


#[test]
fn verify_cli() {
    use clap::CommandFactory;
    Args::command().debug_assert()
}
