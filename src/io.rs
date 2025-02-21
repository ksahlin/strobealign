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

#[cfg(test)]
mod test {
    use std::error::Error;
    use std::io::{Read, Write};
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use temp_file::TempFileBuilder;
    use super::xopen;
    
    #[test]
    fn test_open_multi_block_gzip() -> Result<(), Box<dyn Error>>{
        // Create a multi-block gzip
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(b"abc")?;
        let first = encoder.finish()?;
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(b"def")?;
        let second = encoder.finish()?;
        let both = vec![first, second].concat();
        let tmp = TempFileBuilder::new().suffix(".gz").build()?.with_contents(&both)?;
        
        let mut buf = vec![];
        xopen(tmp.path()).unwrap().read_to_end(&mut buf)?;
        
        assert_eq!(buf, b"abcdef");
        
        Ok(())
    }
}
