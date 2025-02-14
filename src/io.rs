use std::fs::File;
use std::io;
use std::io::{Error, Read};
use std::path::Path;
use flate2::read::MultiGzDecoder;

/// Open an uncompressed or a gzip-compressed file depending on the file name extension
pub fn xopen<P: AsRef<Path>>(path: P) -> Result<Box<dyn Read + Send>, Error> {
    let path = path.as_ref();
    if path == Path::new("-") {
        Ok(Box::new(io::stdin()))
    } else {
        let f = File::open(path)?;
        match path.extension() {
            Some(x) if x == "gz" => Ok(Box::new(MultiGzDecoder::new(f))),
            _ => Ok(Box::new(f)),
        }
    }
}
