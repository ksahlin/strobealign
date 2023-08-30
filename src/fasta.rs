use std::io::{BufRead, Error};

#[derive(Debug)]
pub struct RefSequence {
    pub name: String,
    pub sequence: Vec<u8>,
}

pub fn read_fasta<R: BufRead>(reader: &mut R) -> Result<Vec<RefSequence>, Error> {
    let mut records = Vec::<RefSequence>::new();
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
                records.push(RefSequence {name, sequence});
            }
            name = String::from_utf8(line[1..].to_vec()).unwrap();
            if let Some(i) = name.find(|c: char| c.is_ascii_whitespace()) {
                name = name[..i].to_string();
            }
            sequence = Vec::new();
            has_record = true;
        } else {
            sequence.extend(line.iter().map(|&c| c.to_ascii_uppercase()));
        }
    }
    if has_record {
        records.push(RefSequence {name, sequence});
    }

    Ok(records)
}
