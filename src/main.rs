use std::{env, io};
use std::fs::File;
use std::io::{Error, BufReader, BufWriter, Write, Read};
use std::path::Path;
use std::time::Instant;
use log::{debug, info};
use clap::Parser;
use fastrand::Rng;
use flate2::read::MultiGzDecoder;
use rstrobes::aligner::{Aligner, Scores};
use rstrobes::fastq::{FastqReader, PeekableFastqReader, SequenceRecord};
use rstrobes::fasta;
use rstrobes::index::{IndexParameters, StrobemerIndex};
use rstrobes::insertsize::InsertSizeDistribution;
use rstrobes::maponly::{map_paired_end_read, map_single_end_read};
use rstrobes::mapper::{align_paired_end_read, align_single_end_read, MappingParameters, SamOutput};
use rstrobes::sam::{ReadGroup, SamHeader};
use rstrobes::io::xopen;

mod logger;

const VERSION: &str = env!("CARGO_PKG_VERSION");
const NAME: &str = env!("CARGO_PKG_NAME");


#[derive(Parser, Debug)]
#[command(version, long_about = None)]
struct Args {
    /// Number of threads
    #[arg(short, default_value_t = 1, value_name = "N")]
    threads: usize,

    //args::ValueFlag<int> chunk_size(parser, "INT", "Number of reads processed by a worker thread at once [10000]", {"chunk-size"}, args::Options::Hidden);

    /// Write output to PATH instead of stdout
    #[arg(short = 'o', value_name = "PATH", help_heading = "Input/output")]
    output: Option<String>,

    //args::Flag v(parser, "v", "Verbose output", {'v'});
    //args::Flag no_progress(parser, "no-progress", "Disable progress report (enabled by default if output is a terminal)", {"no-progress"});

    /// Only map reads, no base level alignment (produces PAF file)
    #[arg(short = 'x', help_heading = "Input/output")]
    map_only: bool,

    // SAM output

    /// Emit =/X instead of M CIGAR operations
    #[arg(long, help_heading = "SAM output")]
    eqx: bool,

    /// Do not output the PG header line
    #[arg(long="no-PG", help_heading = "SAM output")]
    no_pg: bool,

    /// Suppress output of unmapped reads
    #[arg(short = 'U', help_heading = "SAM output")]
    only_mapped: bool,

    // args::Flag interleaved(parser, "interleaved", "Interleaved input", {"interleaved"});

    /// Read group ID
    #[arg(long, help_heading = "SAM output")]
    rg_id: Option<String>,

    /// Add read group metadata to SAM header (can be specified multiple times). Example: SM:samplename
    #[arg(long, help_heading = "SAM output")]
    rg: Vec<String>,

    /// Add extra details to SAM records (helpful for debugging)
    #[arg(long, help_heading = "SAM output")]
    details: bool,

    /// Append FASTQ comment to SAM record
    #[arg(short = 'C', help_heading = "SAM output")]
    fastq_comments: bool,

    /// Retain at most N secondary alignments (is upper bounded by -M and depends on -S)
    #[arg(short = 'N', default_value_t = 0, value_name = "N", help_heading = "SAM output")]
    max_secondary: usize,

    // args::ValueFlag<std::string> index_statistics(parser, "PATH", "Print statistics of indexing to PATH", {"index-statistics"});
    // args::Flag i(parser, "index", "Do not map reads; only generate the strobemer index and write it to disk. If read files are provided, they are used to estimate read length", {"create-index", 'i'});
    // args::Flag use_index(parser, "use_index", "Use a pre-generated index previously written with --create-index.", { "use-index" });

    // Seeding arguments

    /// Mean read length. Default: estimated from the first 500 records in the input file
    #[arg(short, help_heading = "Seeding")]
    read_length: Option<usize>,

    /// Maximum seed length. For reasonable values on -l and -u, the seed length distribution is
    /// usually determined by parameters l and u. Then this parameter is only active in regions
    /// where syncmers are very sparse. Default: read_length - 50
    #[arg(short, help_heading = "Seeding")]
    max_seed_length: Option<usize>,

    /// Strobe length (must be less than 32). Default: chosen based on read length
    #[arg(short, help_heading = "Seeding")]
    k: Option<usize>,

    /// Submer size for creating syncmers. k-s must be even. Default: k-4
    #[arg(short, help_heading = "Seeding")]
    s: Option<usize>,

    /// Lower syncmer offset from k/(k-s+1). Start sample second syncmer k/(k-s+1) + l syncmers downstream [0]
    #[arg(short, help_heading = "Seeding")]
    l: Option<isize>,

    /// Upper syncmer offset from k/(k-s+1). End sample second syncmer k/(k-s+1) + u syncmers downstream [7]
    #[arg(short, help_heading = "Seeding")]
    u: Option<isize>,

    /// Bitcount length between 2 and 63. [8]
    #[arg(short, help_heading = "Seeding")]
    c: Option<u32>,

    /// No. of top bits of hash to use as bucket indices (8-31). Default is to determine automatically.
    #[arg(short, help_heading = "Seeding")]
    bits: Option<u8>,

    // Alignment scores

    /// Match score
    #[arg(short = 'A', default_value_t = Scores::default().match_, value_name = "N", help_heading = "Alignment")]
    match_score: u8,

    /// Mismatch penalty
    #[arg(short = 'B', default_value_t = Scores::default().mismatch, value_name = "N", help_heading = "Alignment")]
    mismatch_score: u8,

    /// Gap open penalty
    #[arg(short = 'O', default_value_t = Scores::default().gap_open, value_name = "N", help_heading = "Alignment")]
    gap_open_penalty: u8,

    /// Gap extension penalty
    #[arg(short = 'E', default_value_t = Scores::default().gap_extend, value_name = "N", help_heading = "Alignment")]
    gap_extension_penalty: u8,

    /// Soft-clipping penalty
    #[arg(short = 'L', default_value_t = Scores::default().end_bonus, value_name = "N", help_heading = "Alignment")]
    end_bonus: u32,

    /// Top fraction of repetitive strobemers to filter out from sampling
    #[arg(short, default_value_t = 0.0002, help_heading = "Search parameters")]
    filter_fraction: f64,

    /// Try candidate sites with mapping score at least S of maximum mapping score
    #[arg(short = 'S', default_value_t = 0.5, help_heading = "Search parameters")]
    dropoff_threshold: f32,

    /// Maximum number of mapping sites to try
    #[arg(short = 'M', default_value_t = MappingParameters::default().max_tries, help_heading = "Search parameters")]
    max_tries: usize,

    /// Rescue level. Perform additional search for reads with many repetitive seeds filtered out.
    /// This search includes seeds of R*repetitive_seed_size_filter (default: R=2). Higher R than
    /// default makes strobealign significantly slower but more accurate.
    /// R <= 1 deactivates rescue and is the fastest
    #[arg(short = 'R', default_value_t = 2, help_heading = "Search parameters")]
    rescue_level: usize,

    /// Path to input reference (in FASTA format)
    ref_path: String,

    /// Path to input FASTQ
    fastq_path: String,

    /// Path to input FASTQ with R2 reads (if paired end)
    fastq_path2: Option<String>,
}

fn main() -> Result<(), Error> {
    sigpipe::reset();
    logger::init().unwrap();
    let args = Args::parse();
    info!("This is {} {}", NAME, VERSION);
    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
    let path = Path::new(&args.ref_path);
    let f = xopen(path)?;
    let mut reader = BufReader::new(f);
    let timer = Instant::now();
    let references = fasta::read_fasta(&mut reader)?;
    info!("Time reading references: {:.2} s", timer.elapsed().as_secs_f64());
    let total_ref_size = references.iter().map(|r| r.sequence.len()).sum::<usize>();
    let max_contig_size = references.iter().map(|r| r.sequence.len()).max().expect("No reference found");
    info!("Reference size: {:.2} Mbp ({} contig{}; largest: {:.2} Mbp)",
        total_ref_size as f64 / 1E6,
        references.len(),
        if references.len() != 1 { "s" } else { "" },
        max_contig_size as f64 / 1E6
    );

    let f1 = xopen(&args.fastq_path)?;
    let mut fastq_reader1 = PeekableFastqReader::new(f1);

    let read_length = match args.read_length {
        Some(r) => r,
        None => estimate_read_length(&fastq_reader1.peek(500)?),
    };
    let parameters = IndexParameters::from_read_length(read_length, args.k, args.s, args.l, args.u, args.c, args.max_seed_length);
    info!("Indexing ...");
    debug!("{:?}", parameters);

    let timer = Instant::now();
    let mut index = StrobemerIndex::new(&references, parameters.clone(), args.bits);
    index.populate(args.filter_fraction, args.rescue_level, args.threads);
    let index = index;
    info!("Total time indexing: {:.2} s", timer.elapsed().as_secs_f64());

    let mapping_parameters = MappingParameters {
        max_secondary: args.max_secondary,
        max_tries: args.max_tries,
        dropoff_threshold: args.dropoff_threshold,
        .. MappingParameters::default()
    };

    let scores = Scores {
        match_: args.match_score,
        mismatch: args.mismatch_score,
        gap_open: args.gap_open_penalty,
        gap_extend: args.gap_extension_penalty,
        end_bonus: args.end_bonus,
    };
    let aligner = Aligner::new(scores);

    let cmd_line = env::args().skip(1).collect::<Vec<_>>().join(" ");
    let rg_id = match args.rg_id {
        Some(rg_id) => Some(rg_id),
        None if !args.rg.is_empty() => Some("1".to_string()),
        None => None,
    };
    let sam_output = SamOutput::new(args.details, args.eqx, rg_id.clone(), args.fastq_comments);
    let read_group = rg_id.map(|s| ReadGroup::new(&s, args.rg));

    let header = SamHeader::new(
        &references,
        if args.no_pg { None } else { Some(&cmd_line) },
        VERSION,
        read_group,
    );

    let out: Box<dyn Write> = match args.output {
        Some(output) => Box::new(File::create(output)?),
        None => Box::new(io::stdout().lock()),
    };
    let mut out = BufWriter::new(out);

    if !args.map_only {
        write!(out, "{}", header)?;
    }

    let mut rng = Rng::with_seed(0);
    if let Some(r2_path) = args.fastq_path2 {
        // paired-end reads
        let f2 = xopen(r2_path)?;
        let mut isizedist = InsertSizeDistribution::new();

        for (r1, r2) in fastq_reader1.zip(FastqReader::new(f2)) {
            let r1 = r1?;
            let r2 = r2?;
            if !args.map_only {
                let sam_records = align_paired_end_read(
                    &r1, &r2, &index, &references, &mapping_parameters, &sam_output, &parameters, &mut isizedist, &aligner, &mut rng
                );
                for sam_record in sam_records {
                    if sam_record.is_mapped() || !args.only_mapped {
                        writeln!(out, "{}", sam_record)?;
                    }
                }
            } else {
                let paf_records = map_paired_end_read(
                    &r1, &r2, &index, &references, mapping_parameters.rescue_level, &mut isizedist, &mut rng
                );
                for paf_record in paf_records {
                    writeln!(out, "{}", paf_record)?;
                }

            }
        }

    } else {
        // single-end reads
        for record in fastq_reader1 {
            let record = record?;
            if !args.map_only {
                let sam_records = align_single_end_read(&record, &index, &references, &mapping_parameters, &sam_output, &aligner, &mut rng);
                for sam_record in sam_records {
                    if sam_record.is_mapped() || !args.only_mapped {
                        writeln!(out, "{}", sam_record)?;
                    }
                }
            } else {
                let paf_records = map_single_end_read(&record, &index, &references, mapping_parameters.rescue_level, &mut rng);
                for paf_record in paf_records {
                    writeln!(out, "{}", paf_record)?;
                }
            }
        }
    }

    out.flush()?;

    Ok(())
}

fn estimate_read_length(records: &[SequenceRecord]) -> usize {
    let mut s = 0;
    let mut n = 0;
    for record in records {
        s += record.sequence.len();
        n += 1;
    }

    if n == 0 { 0 } else { s / n }
}

#[test]
fn verify_cli() {
    use clap::CommandFactory;
    Args::command().debug_assert()
}
