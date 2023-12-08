use std::env;
use std::fs::File;
use std::io::{Error, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::Instant;
use log::{debug, info};
use clap::Parser;
use rstrobes::aligner::{Aligner, Scores};
use rstrobes::fastq::FastqReader;
use rstrobes::fasta;
use rstrobes::index::{IndexParameters, StrobemerIndex};
use rstrobes::mapper::{map_single_end_read, MappingParameters, SamOutput};
use rstrobes::sam::SamHeader;

mod logger;

const VERSION: &str = env!("CARGO_PKG_VERSION");
const NAME: &str = env!("CARGO_PKG_NAME");


#[derive(Parser, Debug)]
#[command(version, long_about = None)]
struct Args {
    #[arg(short, default_value_t = 1)]
    threads: usize,

    //args::ValueFlag<int> chunk_size(parser, "INT", "Number of reads processed by a worker thread at once [10000]", {"chunk-size"}, args::Options::Hidden);
    //args::Flag v(parser, "v", "Verbose output", {'v'});
    //args::Flag no_progress(parser, "no-progress", "Disable progress report (enabled by default if output is a terminal)", {"no-progress"});

    /// Emit =/X instead of M CIGAR operations
    #[arg(long)]
    eqx: bool,

    //args::Flag x(parser, "x", "Only map reads, no base level alignment (produces PAF file)", {'x'});

    /// Suppress output of unmapped reads
    #[arg(short = 'U')]
    only_mapped: bool,

    // args::Flag interleaved(parser, "interleaved", "Interleaved input", {"interleaved"});
    // args::ValueFlag<std::string> rgid(parser, "ID", "Read group ID", {"rg-id"});
    // args::ValueFlagList<std::string> rg(parser, "TAG:VALUE", "Add read group metadata to SAM header (can be specified multiple times). Example: SM:samplename", {"rg"});

    /// Add debugging details to SAM records
    #[arg(long)]
    details: bool,

    // args::ValueFlag<int> N(parser, "INT", "Retain at most INT secondary alignments (is upper bounded by -M and depends on -S) [0]", {'N'});
    // args::ValueFlag<std::string> index_statistics(parser, "PATH", "Print statistics of indexing to PATH", {"index-statistics"});
    // args::Flag i(parser, "index", "Do not map reads; only generate the strobemer index and write it to disk. If read files are provided, they are used to estimate read length", {"create-index", 'i'});
    // args::Flag use_index(parser, "use_index", "Use a pre-generated index previously written with --create-index.", { "use-index" });

    // Seeding arguments

    /// Mean read length. Default: estimated from the first 500 records in the input file
    #[arg(short)]
    read_length: Option<usize>,

    /// Strobe length (must be less than 32). Default: chosen based on read length
    #[arg(short)]
    k: Option<usize>,

    /// Submer size for creating syncmers. k-s must be even. Default: k-4
    #[arg(short)]
    s: Option<usize>,

    /// Lower syncmer offset from k/(k-s+1). Start sample second syncmer k/(k-s+1) + l syncmers downstream [0]
    #[arg(short)]
    l: Option<isize>,

    /// Upper syncmer offset from k/(k-s+1). End sample second syncmer k/(k-s+1) + u syncmers downstream [7]
    #[arg(short)]
    u: Option<isize>,

    /// Bitcount length between 2 and 63. [8]
    #[arg(short)]
    c: Option<u32>,

    /// Maximum seed length. For reasonable values on -l and -u, the seed length distribution is
    /// usually determined by parameters l and u. Then this parameter is only active in regions
    /// where syncmers are very sparse. Default: read_length - 50
    #[arg(short)]
    max_seed_length: Option<usize>,

    /// No. of top bits of hash to use as bucket indices (8-31). Default is to determine automatically.
    #[arg(short)]
    bits: Option<u8>,

    // args::Group alignment(parser, "Alignment:");
    // args::ValueFlag<int> A(parser, "INT", "Matching score [2]", {'A'});
    // args::ValueFlag<int> B(parser, "INT", "Mismatch penalty [8]", {'B'});
    // args::ValueFlag<int> O(parser, "INT", "Gap open penalty [12]", {'O'});
    // args::ValueFlag<int> E(parser, "INT", "Gap extension penalty [1]", {'E'});
    // args::ValueFlag<int> end_bonus(parser, "INT", "Soft clipping penalty [10]", {'L'});
    //
    // args::Group search(parser, "Search parameters:");

    /// Top fraction of repetitive strobemers to filter out from sampling
    #[arg(short, default_value_t = 0.0002)]
    filter_fraction: f64,

    // args::ValueFlag<float> S(parser, "FLOAT", "Try candidate sites with mapping score at least S of maximum mapping score [0.5]", {'S'});
    // args::ValueFlag<int> M(parser, "INT", "Maximum number of mapping sites to try [20]", {'M'});

    /// Rescue level. Perform additional search for reads with many repetitive seeds filtered out.
    /// This search includes seeds of R*repetitive_seed_size_filter (default: R=2). Higher R than
    /// default makes strobealign significantly slower but more accurate.
    /// R <= 1 deactivates rescue and is the fastest
    #[arg(short = 'R', default_value_t = 2)]
    rescue_level: usize,

    /// Path to input reference (in FASTA format)
    ref_path: String,

    fastq_path: String,
}

fn main() -> Result<(), Error> {
    sigpipe::reset();
    logger::init().unwrap();
    let args = Args::parse();
    info!("This is {} {}", NAME, VERSION);
    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
    let path = Path::new(&args.ref_path);
    let f = File::open(path)?;
    let mut reader = BufReader::new(f);
    let timer = Instant::now();
    let references = fasta::read_fasta(&mut reader).unwrap();
    info!("Time reading references: {:.2} s", timer.elapsed().as_secs_f64());
    let total_ref_size = references.iter().map(|r| r.sequence.len()).sum::<usize>();
    let max_contig_size = references.iter().map(|r| r.sequence.len()).max().expect("No reference found");
    info!("Reference size: {:.2} Mbp ({} contig{}; largest: {:.2} Mbp)",
        total_ref_size as f64 / 1E6,
        references.len(),
        if references.len() != 1 { "s" } else { "" },
        max_contig_size as f64 / 1E6
    );

    let read_length = match args.read_length {
        Some(r) => r,
        None => estimate_read_length(&args.fastq_path)?,
    };
    let parameters = IndexParameters::from_read_length(read_length, args.k, args.s, args.l, args.u, args.c, args.max_seed_length);
    info!("Indexing ...");
    debug!("{:?}", parameters);

    let timer = Instant::now();
    let mut index = StrobemerIndex::new(&references, parameters, args.bits);
    index.populate(args.filter_fraction, args.rescue_level);
    let index = index;
    info!("Total time indexing: {:.2} s", timer.elapsed().as_secs_f64());

    let mapping_parameters = MappingParameters::default();

    let aligner = Aligner::new(Scores::default());

    let sam_output = SamOutput::new(args.details, args.eqx);
    let cmd_line = env::args().skip(1).collect::<Vec<_>>().join(" ");
    let read_group_fields = vec![];
    let header = SamHeader::new(
        &references,
        &cmd_line,
        VERSION,
        None,
        &read_group_fields,
    );
    print!("{}", header);
    let out = std::io::stdout().lock();
    let mut out = BufWriter::new(out);

    let f = File::open(&args.fastq_path)?;
    for record in FastqReader::new(f) {
        let record = record?;
        let sam_records = map_single_end_read(&record, &index, &references, &mapping_parameters, &sam_output, &aligner);
        for sam_record in sam_records {
            if sam_record.is_mapped() || !args.only_mapped {
                writeln!(out, "{}", sam_record)?;
            }
        }
    }

    Ok(())
}

fn estimate_read_length<P: AsRef<Path>>(path: P) -> Result<usize, Error> {
    let f = File::open(&path)?;
    let mut s = 0;
    let mut n = 0;
    for record in FastqReader::new(f).take(500) {
        let record = record?;
        s += record.sequence.len();
        n += 1;
    }
    Ok(if n == 0 { 0 } else { s / n})
}

#[test]
fn verify_cli() {
    use clap::CommandFactory;
    Args::command().debug_assert()
}
