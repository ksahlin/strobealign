use std::fmt::{Display, Formatter};
use std::str::FromStr;

#[derive(Debug,PartialEq,Eq,Clone,Copy)]
enum CigarOperation {
    Match = 0,
    Insertion = 1,
    Deletion = 2,
    Skip = 3,
    Softclip = 4,
    Hardclip = 5,
    Pad = 6,
    Eq = 7,
    X = 8,
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

#[derive(Debug)]
struct OpLen {
    op: CigarOperation,
    len: u32,
}

#[derive(Default,Debug)]
struct Cigar {
    ops: Vec<OpLen>,
}

impl Cigar {
    pub fn new() -> Self {
        Cigar { ops: vec![] }
    }

    pub fn is_empty(&self) -> bool {
        self.ops.is_empty()
    }

    pub fn push(&mut self, op: CigarOperation, len: u32) {
        if self.ops.is_empty() || self.ops.last().unwrap().op != op {
            self.ops.push(OpLen { op, len });
        } else {
            self.ops.last_mut().unwrap().len += len;
        }
    }

    /// Return a new Cigar that uses M operations instead of =/X
    fn with_m(&self) -> Self {
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
    // void operator+=(const Cigar& other) {
    //     for (auto op_len : other.m_ops) {
    //         push(op_len & 0xf, op_len >> 4);
    //     }
    // }

    /* This works only if I, D, X, = are the only operations used */
    /*pub fn edit_distance(&self) -> usize {
        let dist = 0;
        for op_len in self.ops {
            let len = op_len >> 4;
            match op_len & 0xf {

                _ => {},
            }
            let op = op_len & 0xf;
            if op == CigarOperation::Insertion || op == CIGAR_DEL || op == CIGAR_X) {
                dist += len;
            }
        }
        dist
    }*/

}

    /*Cigar Cigar::to_eqx(const std::string& query, const std::string& ref) const {
    size_t i = 0, j = 0;
    Cigar cigar;
    for (auto op_len : m_ops) {
        auto op = op_len & 0xf;
        auto len = op_len >> 4;
        if (op == CIGAR_MATCH) {
            for (size_t u = 0; u < len; ++u) {
                if (query[i] == ref[j]) {
                    cigar.push(CIGAR_EQ, 1);
                } else {
                    cigar.push(CIGAR_X, 1);
                }
                i++;
                j++;
            }
        } else if (op == CIGAR_INS) {
            cigar.push(op, len);
            i += len;
        } else if (op == CIGAR_DEL) {
            cigar.push(op, len);
            j += len;
        }
    }
    return cigar;
}*/

impl Display for Cigar {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for op_len in &self.ops {
            write!(f, "{}{}", op_len.len, op_len.op)?;
        }
        Ok(())
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
                    let val = ch as u32 - '0' as u32;
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
/*
fn compress_cigar(ops: &String) -> String {
    char prev = 0;
    int count = 0;
    std::stringstream cigar;
    bool first = true;
    for (auto op : ops) {
        if (!first && op != prev) {
            cigar << count << prev;
            count = 0;
        }
        count++;
        prev = op;
        first = false;
    }
    if !first {
        cigar << count << prev;
    }
    return cigar.str();
}*/

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
    fn compress_cigar() {
        // CHECK(compress_cigar("") == "");
        // CHECK(compress_cigar("M") == "1M");
        // CHECK(compress_cigar("MM") == "2M");
        // CHECK(compress_cigar("MMI") == "2M1I");
        // CHECK(compress_cigar("MMII") == "2M2I");
        // CHECK(compress_cigar("MI") == "1M1I");
        // CHECK(compress_cigar("MII") == "1M2I");
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
/*
TEST_CASE("concatenate Cigar") {
    Cigar c{"3M"};
    c += Cigar{"2M1X"};
    CHECK(c.to_string() == "5M1X");
}

TEST_CASE("edit distance") {
    CHECK(Cigar("3=1X4D5I7=").edit_distance() == 10);
}

TEST_CASE("reverse") {
    Cigar c{"3=1X4D5I7="};
    c.reverse();
    CHECK(c.to_string() == "7=5I4D1X3=");
}*/

}
