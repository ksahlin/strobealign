use std::fmt::{Debug, Display, Formatter};
use std::str::FromStr;

#[derive(Debug,PartialEq, Eq, Clone, Copy)]
pub enum CigarOperation {
    Match = 0,
    Insertion = 1,
    Deletion = 2,
    Skip = 3,  // N
    Softclip = 4,
    Hardclip = 5,
    Pad = 6,
    Eq = 7,
    X = 8,
}

impl TryFrom<u8> for CigarOperation {
    type Error = ();
    fn try_from(value: u8) -> Result<CigarOperation, ()> {
        let value = value & 0xf;
        if value > 8 {
            return Err(());
        }
        unsafe { std::mem::transmute(value & 0xf) }
    }
}

impl Display for CigarOperation {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let ch = match self {
            CigarOperation::Match => 'M',
            CigarOperation::Insertion => 'I',
            CigarOperation::Deletion => 'D',
            CigarOperation::Skip => 'N',
            CigarOperation::Softclip => 'S',
            CigarOperation::Hardclip => 'H',
            CigarOperation::Pad => 'P',
            CigarOperation::Eq => '=',
            CigarOperation::X => 'X',
        };
        write!(f, "{}", ch)
    }
}

impl CigarOperation {
    fn from_char(ch: char) -> Result<Self, ()> {
        Ok(match ch {
            'M' => CigarOperation::Match,
            'I' => CigarOperation::Insertion,
            'D' => CigarOperation::Deletion,
            'N' => CigarOperation::Skip,
            'S' => CigarOperation::Softclip,
            'H' => CigarOperation::Hardclip,
            'P' => CigarOperation::Pad,
            '=' => CigarOperation::Eq,
            'X' => CigarOperation::X,
            _ => { return Err(()) },
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct OpLen {
    op: CigarOperation,
    len: usize,
}

impl TryFrom<u32> for OpLen {
    type Error = ();
    fn try_from(value: u32) -> Result<Self, ()> {
        let len = (value >> 4) as usize;
        let op = CigarOperation::try_from((value & 0xf) as u8)?;
        Ok(OpLen { len, op })
    }
}

#[derive(Default, Clone, PartialEq, Eq)]
pub struct Cigar {
    ops: Vec<OpLen>,
}

impl Cigar {
    pub fn new() -> Self {
        Cigar { ops: vec![] }
    }

    pub fn is_empty(&self) -> bool {
        self.ops.is_empty()
    }

    pub fn push(&mut self, op: CigarOperation, len: usize) {
        if len == 0 { return; }
        if self.ops.is_empty() || self.ops.last().unwrap().op != op {
            self.ops.push(OpLen { op, len });
        } else {
            self.ops.last_mut().unwrap().len += len;
        }
    }

    pub fn push_unnormalized(&mut self, op: CigarOperation, len: usize) {
        self.ops.push(OpLen { op, len });
    }

    /// Return a new Cigar that uses M operations instead of =/X
    pub fn with_m(&self) -> Self {
        let mut cigar = Cigar::new();
        for oplen in &self.ops {
            cigar.push(
                match oplen.op {
                    CigarOperation::Eq | CigarOperation::X => CigarOperation::Match,
                    op => op,
                }, oplen.len);
        }
        cigar
    }

    pub fn reversed(&self) -> Self {
        Cigar {
            ops: self.ops.iter().copied().rev().collect()
        }
    }

    pub fn extend(&mut self, other: &Cigar) {
        // TODO use push only for first, then extend
        for oplen in &other.ops {
            self.push(oplen.op, oplen.len)
        }
    }

    /// Return edit distance. Assumes =/X is used - M operations are not counted!
    pub fn edit_distance(&self) -> usize {
        let mut dist = 0;
        for op_len in &self.ops {
            match op_len.op {
                CigarOperation::X|CigarOperation::Insertion|CigarOperation::Deletion => { dist += op_len.len; },
                _ => {},
            }
        }
        dist
    }

    /// Return a new Cigar that uses = and X operations instead of M
    pub fn with_eqx(&self, query: &[u8], refseq: &[u8]) -> Cigar {
        let mut cigar = Cigar::new();
        let mut q_i = 0;
        let mut r_i = 0;
        for oplen in &self.ops {
            let (op, len) = (oplen.op, oplen.len);
            match op {
                CigarOperation::Match => {
                    for _ in 0..len {
                        if query[q_i] == refseq[r_i] {
                            cigar.push(CigarOperation::Eq, 1);
                        } else {
                            cigar.push(CigarOperation::X, 1);
                        }
                        q_i += 1;
                        r_i += 1;
                    }
                },
                CigarOperation::Insertion|CigarOperation::Softclip => {
                    cigar.push(op, len);
                    q_i += len;
                },
                CigarOperation::Deletion|CigarOperation::Skip => {
                    cigar.push(op, len);
                    r_i += len;
                },
                CigarOperation::Hardclip|CigarOperation::Pad => {
                    cigar.push(op, len);
                },
                CigarOperation::Eq|CigarOperation::X => {
                    cigar.push(op, len);
                    r_i += len;
                    q_i += len;
                },
            }
        }

        cigar
    }
}

impl Display for Cigar {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for op_len in &self.ops {
            write!(f, "{}{}", op_len.len, op_len.op)?;
        }
        Ok(())
    }
}

impl Debug for Cigar {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        Display::fmt(self, f)
    }
}

impl FromStr for Cigar {
    type Err = ();

    /// parses ([0-9]*[MIDNSHP=X])*
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut cigar = Cigar::new();
        let mut number = None;
        for ch in s.chars() {
            match ch {
                ' ' => {},
                ch if ch.is_ascii_digit() => {
                    let val = ch as usize - '0' as usize;
                    if let Some(n) = number {
                        number = Some(n * 10 + val)
                    } else {
                        number = Some(val)
                    }
                }
                ch => {
                    if let Ok(op) = CigarOperation::from_char(ch) {
                        cigar.push(op, number.unwrap_or(1));
                        number = None;
                    } else {
                        return Err(());
                    }
                },
            }
        }
        Ok(cigar)
    }
}

impl TryFrom<&[u32]> for Cigar {
    type Error = ();
    fn try_from(encoded: &[u32]) -> Result<Self, ()> {
        let ops: Result<Vec<_>, _> = encoded.iter().map(|ol| OpLen::try_from(*ol)).collect();
        Ok(Cigar { ops: ops? })
    }
}

#[cfg(test)]
mod test {
    use std::str::FromStr;
    use super::{Cigar, CigarOperation};

    #[test]
    fn is_empty() {
        assert!(Cigar::default().is_empty());
        assert!(Cigar::new().is_empty());
    }

    #[test]
    fn not_empty() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOperation::Eq, 1);
        assert!(!cigar.is_empty());
    }

    #[test]
    fn test_construct() {
        let cigar = Cigar::new();
        assert_eq!(cigar.to_string(), "");

        let mut cigar = Cigar::new();
        cigar.push(CigarOperation::Match, 1);
        assert_eq!(cigar.to_string(), "1M");

        cigar.push(CigarOperation::Match, 1);
        assert_eq!(cigar.to_string(), "2M");

        cigar.push(CigarOperation::Insertion, 3);
        assert_eq!(cigar.to_string(), "2M3I");
        /*
        uint32_t ops[4] = {
            3 << 4 | CIGAR_MATCH,
            5 << 4 | CIGAR_X,
            7 << 4 | CIGAR_INS,
            13 << 4 | CIGAR_DEL
        };
        Cigar c2{ops, 4};
        CHECK(c2.m_ops.size() == 4);
        CHECK(c2.to_string() == "3M5X7I13D");
        */
    }

    #[test]
    fn parse_cigar() {
        assert!("".parse::<Cigar>().unwrap().is_empty());
        assert_eq!("1M".parse::<Cigar>().unwrap().to_string(), "1M");
        assert_eq!("2M".parse::<Cigar>().unwrap().to_string(), "2M");
        assert_eq!("11M".parse::<Cigar>().unwrap().to_string(), "11M");
        assert_eq!("10M2I1D99=1X4P5S".parse::<Cigar>().unwrap().to_string(), "10M2I1D99=1X4P5S");
        // Not standard, only for convenience
        assert_eq!("M".parse::<Cigar>().unwrap().to_string(), "1M");
        assert_eq!("M M".parse::<Cigar>().unwrap().to_string(), "2M");
        assert_eq!("MMII".parse::<Cigar>().unwrap().to_string(), "2M2I");
        assert_eq!("M 2M X".parse::<Cigar>().unwrap().to_string(), "3M1X");
        assert_eq!("1M 2D 1M".parse::<Cigar>().unwrap().to_string(), "1M2D1M");
    }

    /*
    #[test]
    fn convert_m_to_eqx() {
        let cigar = Cigar::from_str("1M").unwrap();
        assert_eq!(cigar.to_eqx("A", "A").to_string(), "1=");
        assert_eq!(cigar.to_eqx("A", "G").to_string(), "1X");

        let cigar = Cigar::from_str("3M").unwrap();
        assert_eq!(cigar.to_eqx("AAA", "AAA").to_string(), "3=");
        assert_eq!(cigar.to_eqx("AAA", "ATA").to_string(), "1=1X1=");
        assert_eq!(cigar.to_eqx("AAA", "ATT").to_string(), "1=2X");

        let cigar = Cigar::from_str("2M 1D 4M 1I 3M");
        assert_eq!(cigar.to_eqx("ACTTTGCATT", "ACGTATGAAA").to_string(), "2=1D1=1X2=1I1=2X");
    }*/

    #[test]
    fn convert_eqx_to_m() {
        assert!(Cigar::from_str("").unwrap().with_m().is_empty());
        assert_eq!(Cigar::from_str("5=").unwrap().with_m().to_string(), "5M");
        assert_eq!(Cigar::from_str("5S3=1X2=4S").unwrap().with_m().to_string(), "5S6M4S");
    }
    
    #[test]
    fn test_reverse() {
        let c = Cigar::from_str("3=1X4D5I7=").unwrap();
        assert_eq!(c.reversed(), Cigar::from_str("7=5I4D1X3=").unwrap());
        
    }
/*
TEST_CASE("concatenate Cigar") {
    Cigar c{"3M"};
    c += Cigar{"2M1X"};
    CHECK(c.to_string() == "5M1X");
}

TEST_CASE("edit distance") {
    CHECK(Cigar("3=1X4D5I7=").edit_distance() == 10);
}
*/

}
