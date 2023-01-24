use std::env;
use std::fs::File;
use std::io::{BufReader, BufRead, Error};
use std::path::Path;

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
    let args = &env::args().collect::<Vec<_>>();
    let path = Path::new(&args[1]);
    let f = File::open(path)?;
    let mut reader = BufReader::new(f);
    let records = read_fasta(&mut reader).unwrap();
    for record in &records {
        println!(">{}\n{}", record.name, String::from_utf8(record.sequence.to_vec()).unwrap());
    }
    Ok(())
}
