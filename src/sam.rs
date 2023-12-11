use std::fmt::{Display, Formatter};
use crate::cigar::Cigar;
use crate::details::Details;
use crate::fasta::RefSequence;

pub const PAIRED: u16 = 1;
pub const PROPER_PAIR: u16 = 2;
pub const UNMAP: u16 = 4;
pub const MUNMAP: u16 = 8;
pub const REVERSE: u16 = 0x10;
pub const MREVERSE: u16 = 0x20;
pub const READ1: u16 = 0x40;
pub const READ2: u16 = 0x80;
pub const SECONDARY: u16 = 0x100;
pub const QCFAIL: u16 = 0x200;
pub const DUP: u16 = 0x400;
pub const SUPPLEMENTARY: u16 = 0x800;

// TODO String and Vec<u8> fields should be references
#[derive(Default)]
pub struct SamRecord {
    pub query_name: String,
    pub flags: u16,
    pub reference_name: Option<String>,
    /// 0-based position, converted to 1-based for output
    pub pos: Option<u32>,
    pub mapq: u8,
    pub cigar: Option<Cigar>,
    pub mate_reference_name: Option<String>,
    pub mate_pos: Option<u32>,
    pub template_len: Option<i32>,
    pub query_sequence: Option<Vec<u8>>,
    pub query_qualities: Option<Vec<u8>>,
    pub edit_distance: Option<u32>,
    pub alignment_score: Option<u32>,
    pub details: Option<Details>,
    pub rg_id: Option<String>,
}

impl SamRecord {
    pub fn is_secondary(&self) -> bool {
        self.flags & SECONDARY != 0
    }

    pub fn is_mapped(&self) -> bool {
        self.flags & UNMAP == 0
    }
}

impl Display for SamRecord {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        // SAM columns
        //
        // 1 QNAME query template name
        // 2 FLAG  bitwise flag
        // 3 RNAME reference sequence name
        // 4 POS   1-based leftmost mapping position (set to 0 for unmapped read w/o pos)
        // 5 MAPQ  mapping quality (255 if unavailable)
        // 6 CIGAR
        // 7 RNEXT reference name of mate/next read
        // 8 PNEXT position of mate/next read
        // 9 TLEN template length
        // 10 SEQ
        // 11 QUAL
        // 12 optional fields

        let mate_pos = match self.mate_pos {
            Some(pos) => pos + 1,
            None => 0,
        };
        let query_sequence = match &self.query_sequence {
            Some(seq) if !self.is_secondary() && !seq.is_empty() => std::str::from_utf8(seq).unwrap(),
            _ => "*",
        };
        let query_qualities = match &self.query_qualities {
            Some(query_qualities) if !self.is_secondary() && !query_qualities.is_empty() => std::str::from_utf8(query_qualities).unwrap(),
            _ => "*",
        };
        let pos = match self.pos {
            Some(pos) => pos + 1,
            None => 0,
        };
        let cigar = match &self.cigar {
            Some(cigar) => cigar.to_string(),
            None => "*".to_string(),
        };
        write!(f, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_name,
            self.flags,
            self.reference_name.as_ref().unwrap_or(&"*".to_string()),
            pos,
            self.mapq,
            cigar,
            self.mate_reference_name.as_ref().unwrap_or(&"*".to_string()),
            mate_pos,
            self.template_len.unwrap_or(0),
            query_sequence,
            query_qualities,
        )?;
        if let Some(edit_distance) = self.edit_distance {
            write!(f, "\tNM:i:{}", edit_distance)?;
        }
        if let Some(alignment_score) = self.alignment_score {
            write!(f, "\tAS:i:{}", alignment_score)?;
        }
        if let Some(details) = &self.details {
            write!(f, "\tna:i:{}\tnr:i:{}\tal:i:{}\tga:i:{}",
                details.nams,
                details.nam_rescue as u8,
                details.tried_alignment,
                details.gapped
            )?;
        }
        if let Some(rg_id) = &self.rg_id {
            write!(f, "\tRG:Z:{}", rg_id)?;
        }
        Ok(())
    }
}

pub struct ReadGroup {
    id: String,
    fields: Vec<String>,
}

impl ReadGroup {
    pub fn new(id: &str, fields: Vec<String>) -> Self {
        ReadGroup {
            id: id.to_string(),
            fields,
        }
    }
}

pub struct SamHeader<'a> {
    references: &'a [RefSequence],
    cmd_line: Option<&'a str>,
    version: &'a str,
    read_group: Option<ReadGroup>,
}

impl<'a> SamHeader<'a> {
    pub fn new(
        references: &'a [RefSequence],
        cmd_line: Option<&'a str>,
        version: &'a str,
        read_group: Option<ReadGroup>,
    ) -> Self {
        SamHeader {
            references,
            cmd_line,
            version,
            read_group,
        }
    }
}

impl<'a> Display for SamHeader<'a> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "@HD\tVN:1.6\tSO:unsorted")?;

        for refseq in self.references {
            writeln!(f, "@SQ\tSN:{}\tLN:{}", refseq.name, refseq.sequence.len())?;
        }
        if let Some(read_group) = &self.read_group {
            write!(f, "@RG\tID:{}", read_group.id)?;
            for field in &read_group.fields {
                write!(f, "\t{}", &field)?;
            }
            f.write_str("\n")?;
        }
        if let Some(cmd_line) = self.cmd_line {
            writeln!(f, "@PG\tID:strobealign\tPN:strobealign\tVN:{}\tCL:{}", &self.version, &cmd_line)?;
        };

        Ok(())
    }
}
