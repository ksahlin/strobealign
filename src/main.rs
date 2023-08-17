use std::io;
use std::fs::File;
use std::io::{BufReader, BufWriter, Error, Write};
use std::path::Path;
use clap::Parser;
use rstrobes::aligner::{Aligner, Scores};
use rstrobes::fastq::FastqReader;
use rstrobes::strobes::RandstrobeIterator;
use rstrobes::syncmers::SyncmerIterator;
use rstrobes::fasta;
use rstrobes::fasta::RefSequence;
use rstrobes::index::{IndexParameters, StrobemerIndex};
use rstrobes::mapper::{map_single_end_read, MappingParameters};

#[derive(Parser, Debug)]
#[command(long_about = None)]
struct Args {

    /// Print syncmers instead of randstrobes
    //#[arg(long, default_value_t = false)]
    //syncmers: bool,

    /// Path to input FASTA
    ref_path: String,

    fastq_path: String,

    #[arg(short, default_value_t = 1)]
    threads: usize,

    // SeedingArguments(args::ArgumentParser& parser)
    // : parser(parser)
    // //n{parser, "INT", "Number of strobes [2]", {'n'}}
    // , r{parser, "INT",
    //     "Mean read length. This parameter is estimated from the first 500 "
    //     "records in each read file. No need to set this explicitly unless you have a reason.", {'r'}}
    // , m{parser, "INT",
    //     "Maximum seed length. Defaults to r - 50. For reasonable values on -l and -u, "
    //     "the seed length distribution is usually determined by parameters l and u. "
    //     "Then, this parameter is only active in regions where syncmers are very sparse.", {'m'}}
    // , k{parser, "INT", "Strobe length, has to be below 32. [20]", {'k'}}
    // , l{parser, "INT", "Lower syncmer offset from k/(k-s+1). Start sample second syncmer k/(k-s+1) + l syncmers downstream [0]", {'l'}}
    // , u{parser, "INT", "Upper syncmer offset from k/(k-s+1). End sample second syncmer k/(k-s+1) + u syncmers downstream [7]", {'u'}}
    // , c{parser, "INT", "Bitcount length between 2 and 63. [8]", {'c'}}
    // , s{parser, "INT",
    //     "Submer size used for creating syncmers [k-4]. Only even numbers on k-s allowed. "
    //     "A value of s=k-4 roughly represents w=10 as minimizer window [k-4]. "
    //     "It is recommended not to change this parameter unless you have a good "
    //     "understanding of syncmers as it will drastically change the memory usage and "
    //     "results with non default values.", {'s'}}
    // , bits{parser, "INT", ""
    //     "[determined from reference size]", {'b'}}

    /// No. of top bits of hash to use as bucket indices (8-31)
    #[arg(short)]
    bits: Option<u8>,

    /// Top fraction of repetitive strobemers to filter out from sampling
    #[arg(short, default_value_t = 0.0002)]
    filter_fraction: f64,
}

fn main() -> Result<(), Error> {
    let args = Args::parse();
    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
    let path = Path::new(&args.ref_path);
    let f = File::open(path)?;
    let mut reader = BufReader::new(f);
    let references = fasta::read_fasta(&mut reader).unwrap();

    // IndexParameters(r=150, k=20, s=16, l=1, u=7, q=255, max_dist=80, t_syncmer=3, w_min=5, w_max=11)
    let parameters = IndexParameters::new(150, 20, 16, 1, 7, 255, 80);
    debug_assert_eq!(parameters.randstrobe.w_min, 5);
    debug_assert_eq!(parameters.randstrobe.w_max, 11);

    let mut index = StrobemerIndex::new(&references, parameters, args.bits);
    index.populate(args.filter_fraction);
    let index = index;

    let mapping_parameters = MappingParameters::default();

    let aligner = Aligner::new(Scores::default());
    // let mapper = Mapper::new(index);
    let f = File::open(&args.fastq_path)?;
    for record in FastqReader::new(f) {
        let record = record?;
        println!("Processing record {}", record.name);
        map_single_end_read(&record, &index, &references, &mapping_parameters, &aligner);
    }

    Ok(())
}

fn dumpstrobes(references: &Vec<RefSequence>, parameters: IndexParameters) -> Result<(), Error> {
    let mut writer = BufWriter::new(io::stdout());

    let k = parameters.syncmer.k;
    let s = parameters.syncmer.s;
    for record in references {
        let name = &record.name;
        let mut syncmers = SyncmerIterator::new(&record.sequence, k, s, 3);
        for randstrobe in RandstrobeIterator::new(&mut syncmers, &parameters.randstrobe) {
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
