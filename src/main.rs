use std::env;
use std::fs::File;
use std::io::{BufReader, BufWriter, Error, Write};
use std::path::Path;
use clap::Parser;
use rstrobes::aligner::{Aligner, Scores};
use rstrobes::fastq::FastqReader;
use rstrobes::fasta;
use rstrobes::index::{IndexParameters, StrobemerIndex};
use rstrobes::mapper::{map_single_end_read, MappingParameters};
use rstrobes::sam::SamHeader;

const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Parser, Debug)]
#[command(version, long_about = None)]
struct Args {
    #[arg(short, default_value_t = 1)]
    threads: usize,

    //args::ValueFlag<int> chunk_size(parser, "INT", "Number of reads processed by a worker thread at once [10000]", {"chunk-size"}, args::Options::Hidden);
    //args::Flag v(parser, "v", "Verbose output", {'v'});
    //args::Flag no_progress(parser, "no-progress", "Disable progress report (enabled by default if output is a terminal)", {"no-progress"});
    //args::Flag eqx(parser, "eqx", "Emit =/X instead of M CIGAR operations", {"eqx"});
    //args::Flag x(parser, "x", "Only map reads, no base level alignment (produces PAF file)", {'x'});
    // args::Flag U(parser, "U", "Suppress output of unmapped reads", {'U'});
    // args::Flag interleaved(parser, "interleaved", "Interleaved input", {"interleaved"});
    // args::ValueFlag<std::string> rgid(parser, "ID", "Read group ID", {"rg-id"});
    // args::ValueFlagList<std::string> rg(parser, "TAG:VALUE", "Add read group metadata to SAM header (can be specified multiple times). Example: SM:samplename", {"rg"});
    // args::Flag details(parser, "details", "Add debugging details to SAM records", {"details"});
    //
    // args::ValueFlag<int> N(parser, "INT", "Retain at most INT secondary alignments (is upper bounded by -M and depends on -S) [0]", {'N'});
    // args::ValueFlag<std::string> index_statistics(parser, "PATH", "Print statistics of indexing to PATH", {"index-statistics"});
    // args::Flag i(parser, "index", "Do not map reads; only generate the strobemer index and write it to disk. If read files are provided, they are used to estimate read length", {"create-index", 'i'});
    // args::Flag use_index(parser, "use_index", "Use a pre-generated index previously written with --create-index.", { "use-index" });
    //
    // args::Group seeding_group(parser, "Seeding:");
    // auto seeding = SeedingArguments{parser};

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

    /// No. of top bits of hash to use as bucket indices (8-31). Default is determine automatically.
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
    index.populate(args.filter_fraction, args.rescue_level);
    let index = index;

    let mapping_parameters = MappingParameters::default();

    let aligner = Aligner::new(Scores::default());

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
        let sam_records = map_single_end_read(&record, &index, &references, &mapping_parameters, &aligner);
        for sam_record in sam_records {
            writeln!(out, "{}", sam_record)?;
        }
    }

    Ok(())
}

#[test]
fn verify_cli() {
    use clap::CommandFactory;
    Args::command().debug_assert()
}
