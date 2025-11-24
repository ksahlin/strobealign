use rstrobes::strobes::DEFAULT_AUX_LEN;
use std::{env, io, thread, time};
use std::cmp::min;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, IsTerminal, Write};
use std::process::exit;
use std::sync::{Arc, Mutex};
use std::sync::mpsc::{channel, sync_channel, Receiver, Sender};
use std::time::Instant;
use clap::builder::Styles;
use clap::builder::styling::AnsiColor;
use log::{debug, error, info, trace};
use clap::Parser;
use fastrand::Rng;
use thiserror::Error;
use rstrobes::aligner::{Aligner, Scores};
use rstrobes::chainer::{Chainer, ChainingParameters};
use rstrobes::details::Details;
use rstrobes::fastq::{interleaved_record_iterator, record_iterator, PeekableSequenceReader, SequenceRecord};
use rstrobes::fasta;
use rstrobes::fasta::{FastaError, RefSequence};
use rstrobes::index::{IndexParameters, StrobemerIndex, REF_RANDSTROBE_MAX_NUMBER_OF_REFERENCES};
use rstrobes::insertsize::InsertSizeDistribution;
use rstrobes::maponly::{abundances_paired_end_read, abundances_single_end_read, map_paired_end_read, map_single_end_read};
use rstrobes::mapper::{align_paired_end_read, align_single_end_read, MappingParameters, SamOutput};
use rstrobes::sam::{ReadGroup, SamHeader};
use rstrobes::io::xopen;
use rstrobes::mcsstrategy::McsStrategy;

mod logger;

const VERSION: &str = env!("CARGO_PKG_VERSION");
const NAME: &str = env!("CARGO_PKG_NAME");
const STYLES: Styles = Styles::plain()
    .header(AnsiColor::Blue.on_default())
    .usage(AnsiColor::Blue.on_default())
    .placeholder(AnsiColor::Green.on_default());

#[derive(Parser, Debug)]
#[command(version, long_about = None, styles = STYLES)]
struct Args {
    /// Number of threads
    #[arg(short, default_value_t = 1, value_name = "N")]
    threads: usize,

    /// Number of threads for indexing (default: same as -t)
    #[arg(long = "ithreads", hide = true)]
    indexing_threads: Option<usize>,

    /// Number of nucleotides processed by a worker thread at once
    #[arg(long, default_value_t = 1_000_000, hide = true)]
    chunk_size: usize,

    /// Write output to PATH instead of stdout
    #[arg(short = 'o', value_name = "PATH", help_heading = "Input/output")]
    output: Option<String>,

    /// Verbose output
    #[arg(short = 'v', help_heading = "Input/output", conflicts_with = "trace")]
    verbose: bool,

    /// Highly verbose output
    #[arg(long = "trace", hide = true)]
    trace: bool,

    /// Disable progress report (enabled by default if output is a terminal)
    #[arg(long = "no-progress", default_value_t = false, help_heading = "Input/output")]
    no_progress: bool,

    /// Only map reads, no base level alignment (produces PAF file)
    #[arg(short = 'x', help_heading = "Input/output")]
    map_only: bool,

    /// Output the estimated abundance per contig (tabular output)
    #[arg(long, help_heading = "Input/output", conflicts_with = "map_only")]
    aemb: bool,

    /// Interleaved (or mixed single-end/paired-end) input
    #[arg(short = 'p', long, help_heading = "Input/output", conflicts_with = "fastq_path2")]
    interleaved: bool,

    // args::ValueFlag<std::string> index_statistics(parser, "PATH", "Print statistics of indexing to PATH", {"index-statistics"});

    /// Do not map reads; only generate the strobemer index and write it to disk.
    /// If read files are provided, they are used to estimate read length
    #[arg(short = 'i',  long = "create-index", help_heading = "Input/output")]
    create_index: bool,

    // args::Flag use_index(parser, "use_index", "Use a pre-generated index previously written with --create-index.", { "use-index" });

    // SAM output

    /// Emit =/X instead of M CIGAR operations
    #[arg(long, help_heading = "SAM output")]
    eqx: bool,

    /// Do not output the PG header line
    #[arg(long="no-PG", help_heading = "SAM output")]
    no_pg: bool,

    /// Do not output unmapped single-end reads. Do not output pairs where both reads are unmapped
    #[arg(short = 'U', help_heading = "SAM output")]
    only_mapped: bool,

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

    // Seeding arguments

    /// Mean read length. Default: estimated from the first 500 records in the input file
    #[arg(short, help_heading = "Seeding")]
    read_length: Option<usize>,

    /// Maximum seed length.
    /// For reasonable values on -l and -u, the seed length distribution is
    /// usually determined by parameters l and u. Then this parameter is only active in regions
    /// where syncmers are very sparse. Default: read_length - 50
    #[arg(short, help_heading = "Seeding")]
    max_seed_length: Option<usize>,

    /// Syncmer (strobe) length (must be less than 32). Default: chosen based on read length
    #[arg(short, help_heading = "Seeding")]
    k: Option<usize>,

    /// Submer size for creating syncmers. k-s must be even. Default: k-4
    #[arg(short, help_heading = "Seeding")]
    s: Option<usize>,

    /// Start of sampling window for second syncmer (i.e., second syncmer must be at least l syncmers downstream). Default: 5
    #[arg(short, help_heading = "Seeding")]
    l: Option<usize>,

    /// End of sampling window for second syncmer (i.e., second syncmer must be at most u syncmers downstream). Default: 11
    #[arg(short, help_heading = "Seeding")]
    u: Option<usize>,

    /// Bitcount length between 2 and 63. Default: 8
    #[arg(short, help_heading = "Seeding")]
    c: Option<u32>,

    /// No. of top bits of hash to use as bucket indices (8-31). Default is to determine automatically.
    #[arg(short, help_heading = "Seeding")]
    bits: Option<u8>,

    /// No. of bits to use from secondary strobe hash
    #[arg(long, default_value_t = DEFAULT_AUX_LEN, value_name = "N", help_heading = "Seeding")]
    aux_len: u8,

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

    //args::Flag nams(parser, "nams", "Use NAMs instead of collinear chaining for alignments", {"nams"});

    /// Collinear chaining look back heuristic
    #[arg(short = 'H', default_value_t = 50, value_name = "N", help_heading = "Collinear chaining")]
    max_lookback: usize,

    /// Collinear chaining diagonal gap cost
    #[arg(long = "gd", default_value_t = 0.1, help_heading = "Collinear chaining")]
    diag_diff_penalty: f32,

    /// Collinear chaining gap length cost
    #[arg(long = "gl", default_value_t = 0.05, help_heading = "Collinear chaining")]
    gap_length_penalty: f32,

    /// Collinear chaining best chain score threshold
    #[arg(long = "vp", default_value_t = 0.7, help_heading = "Collinear chaining")]
    valid_score_threshold: f32,

    /// Collinear chaining skip distance, how far on the reference do we allow anchors to chain
    #[arg(long = "sg", default_value_t = 10000, help_heading = "Collinear chaining")]
    max_ref_gap: usize,

    /// Weight given to the number of anchors for the final score of chains
    #[arg(long = "mw", default_value_t = 0.01, help_heading = "Collinear chaining")]
    matches_weight: f32,

    /// Multi-context seed strategy for finding hits
    #[arg(long = "mcs", value_enum, default_value_t = McsStrategy::default(), help_heading = "Search parameters")]
    mcs_strategy: McsStrategy,

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
    fastq_path: Option<String>,

    /// Path to input FASTQ with R2 reads (if paired end)
    fastq_path2: Option<String>,
}

#[derive(Debug, Error)]
enum CliError {
    #[error("{0}")]
    Io(#[from] io::Error),

    #[error("{0}")]
    FastaError(#[from] FastaError),
}

fn main() -> Result<(), CliError> {
    sigpipe::reset();
    let args = Args::parse();
    let level = if args.trace { log::Level::Trace } else if args.verbose { log::Level::Debug } else { log::Level::Info };
    logger::init(level).unwrap();
    info!("This is {} {}", NAME, VERSION);

    // Open R1 FASTQ file and estimate read length if necessary
    let read_length;
    let fastq_reader1;

    if let Some(fastq_path) = args.fastq_path {
        let f1 = xopen(&fastq_path)?;
        let mut fastq_reader = PeekableSequenceReader::new(f1);
        read_length = match args.read_length {
            Some(r) => r,
            None => {
                let r = estimate_read_length(&fastq_reader.peek(500)?);
                info!("Estimated read length: {} bp", r);
                r
            },
        };
        fastq_reader1 = Some(fastq_reader);
    } else {
        if !args.create_index {
            error!("FASTQ path is required");
            exit(1);
        }
        if let Some(rl) = args.read_length {
            read_length = rl;
        } else {
            error!("With --create-index, either provide a FASTQ path or specify the read length with -r");
            exit(1);
        }
        fastq_reader1 = None;
    }

    let parameters = IndexParameters::from_read_length(read_length, args.k, args.s, args.l, args.u, args.c, args.max_seed_length, args.aux_len);

    info!("Canonical read length: {} bp", parameters.canonical_read_length);
    debug!("  {:?}", parameters.syncmer);
    debug!("  {:?}", parameters.randstrobe);
    debug!("  Maximum seed length: {}", parameters.randstrobe.max_dist as usize + parameters.syncmer.k);
    {
        let d = parameters.syncmer.k - parameters.syncmer.s + 1;
        debug!("  Syncmers are on average sampled every k - s + 1 = {} nucleotides", d);
        debug!("  Sampling window for second syncmer (in syncmers): [{}, {}]", parameters.randstrobe.w_min, parameters.randstrobe.w_max);
        debug!("  Sampling window for second syncmer (in nucleotides): [{}, {}]", parameters.randstrobe.w_min * d, parameters.randstrobe.w_max * d);
    }

    // Read reference FASTA
    let timer = Instant::now();
    let mut reader = BufReader::new(xopen(&args.ref_path)?);
    let references = fasta::read_fasta(&mut reader)?;
    drop(reader);
    info!("Time reading reference: {:.2} s", timer.elapsed().as_secs_f64());
    let total_ref_size = references.iter().map(|r| r.sequence.len()).sum::<usize>();
    let max_contig_size = references.iter().map(|r| r.sequence.len()).max().expect("No reference found");
    info!("Reference size: {:.2} Mbp ({} contig{}; largest: {:.2} Mbp)",
        total_ref_size as f64 / 1E6,
        references.len(),
        if references.len() != 1 { "s" } else { "" },
        max_contig_size as f64 / 1E6
    );
    if references.len() > REF_RANDSTROBE_MAX_NUMBER_OF_REFERENCES {
        error!("Too many reference sequences. Current maximum is {}.", REF_RANDSTROBE_MAX_NUMBER_OF_REFERENCES);
        exit(1);
    }

    // Create the index
    let timer = Instant::now();
    let mut index = StrobemerIndex::new(&references, parameters.clone(), args.bits);
    debug!("Auxiliary hash length: {}", args.aux_len);
    info!("Multi-context seed strategy: {}", args.mcs_strategy);
    info!("Bits used to index buckets: {}", index.bits);
    let indexing_threads = args.indexing_threads.unwrap_or(args.threads);
    info!("Indexing the reference using {} thread{} ...", indexing_threads, if indexing_threads == 1 { "" } else { "s" });
    index.populate(args.filter_fraction, indexing_threads);
    let index = index;
    info!("Total time indexing: {:.2} s", timer.elapsed().as_secs_f64());
    debug!("{}", &index.stats);
    debug!("Filtered cutoff count: {}", index.filter_cutoff);
    debug!("Using rescue cutoff: {}", index.rescue_cutoff);

    if args.create_index {
        let timer = Instant::now();
        let sti_path = args.ref_path + &parameters.filename_extension();
        info!("Writing index to {}", sti_path);
        index.write(sti_path)?;
        info!("Total time writing index: {:.2} s", timer.elapsed().as_secs_f32());

        exit(0);
    }
    let fastq_reader1 = fastq_reader1.unwrap();

    let timer = Instant::now();
    let mapping_parameters = MappingParameters {
        max_secondary: args.max_secondary,
        max_tries: args.max_tries,
        dropoff_threshold: args.dropoff_threshold,
        rescue_level: args.rescue_level,
        output_unmapped: !args.only_mapped,
        .. MappingParameters::default()
    };

    let chaining_parameters = ChainingParameters {
        max_lookback: args.max_lookback,
        diag_diff_penalty: args.diag_diff_penalty,
        gap_length_penalty: args.gap_length_penalty,
        valid_score_threshold: args.valid_score_threshold,
        max_ref_gap: args.max_ref_gap,
        matches_weight: args.matches_weight,
    };
    let scores = Scores {
        match_: args.match_score,
        mismatch: args.mismatch_score,
        gap_open: args.gap_open_penalty,
        gap_extend: args.gap_extension_penalty,
        end_bonus: args.end_bonus,
    };
    debug!("{:?}", &mapping_parameters);
    debug!("{:?}", &chaining_parameters);
    debug!("{:?}", &scores);

    let chainer = Chainer::new(index.k(), chaining_parameters);
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

    let out: Box<dyn Write + Send> = match args.output {
        Some(output) => Box::new(File::create(output)?),
        None => Box::new(io::stdout()),
    };
    let mut out = BufWriter::new(out);

    let mode = if args.map_only { Mode::Paf } else if args.aemb { Mode::Abundances } else { Mode::Sam };
    if mode == Mode::Sam {
        write!(out, "{}", header)?;
    }

    let mut record_iter = if args.interleaved {
        interleaved_record_iterator(fastq_reader1)
    } else {
        record_iterator(fastq_reader1, args.fastq_path2.as_deref())?
    };

    let chunks_iter = std::iter::from_fn(move || {
        let mut nucleotides = 0;
        let mut chunk = vec![];
        for record in record_iter.by_ref() {
            nucleotides += record.as_ref().map_or(0, |r| r.0.len() + r.1.as_ref().map_or(0, |r2| r2.len()));
            chunk.push(record);
            if nucleotides > args.chunk_size {
                break;
            }
        }
        if chunk.is_empty() { None } else { Some(chunk) }
    });

    let mapper = Mapper {
        index: &index,
        references: &references,
        mapping_parameters: &mapping_parameters,
        index_parameters: &parameters,
        chainer: &chainer,
        sam_output: &sam_output,
        aligner,
        include_unmapped: !args.only_mapped,
        mode,
        abundances: vec![0f64; references.len()],
    };

    let mode_message = match mode {
        Mode::Sam => "",
        Mode::Paf => " in mapping-only mode",
        Mode::Abundances => " in abundance estimation mode"
    };
    info!("Processing reads{} using {} thread{}", mode_message, args.threads, if args.threads != 1 { "s" } else {""});

    let show_progress = !args.no_progress && io::stderr().is_terminal();
    // TODO channel size?
    let (chunks_tx, chunks_rx) = sync_channel(args.threads);
    let reader_thread = thread::spawn(move || {
        for chunk in chunks_iter.enumerate() {
            chunks_tx.send(chunk).unwrap();
        }
        drop(chunks_tx);
    });
    let chunks_rx = Arc::new(Mutex::new(chunks_rx));

    let (out_tx, out_rx): (Sender<(usize, (Vec<u8>, Details))>, Receiver<(usize, (Vec<u8>, Details))>) = channel();
    let (progress_tx, progress_rx): (Sender<usize>, Receiver<_>) = channel();
    let (details_tx, details_rx): (Sender<Details>, Receiver<Details>) = channel();
    let out_arc = Arc::new(Mutex::new(out));
    let writer_thread = thread::spawn(move || {
        let mut map = HashMap::new();
        let mut next_index = 0;
        let mut total_details = Details::default();
        for (index, msg) in out_rx {
            map.insert(index, msg);
            while let Some((data, details)) = map.remove(&next_index) {
                out_arc.lock().unwrap().write_all(data.as_ref()).unwrap();
                total_details += details;
                if show_progress {
                    progress_tx.send(total_details.nam.n_reads).unwrap();
                }
                next_index += 1;
            }
        }
        details_tx.send(total_details).unwrap();
        drop(progress_tx);

        out_arc
    });

    let (abundances_tx, abundances_rx) = channel();

    // TODO
    // this thread should be joined after the workers are done
    if show_progress {
        thread::spawn(move || {
            let timer = Instant::now();
            let mut time_to_wait = time::Duration::from_millis(10);
            let mut reported_time = Instant::now();
            let mut total_reads = 0;
            for n_reads in progress_rx {
                total_reads = n_reads;
                // Start with a small waiting time between updates so that there’s no delay if
                // there are very few reads to align.
                time_to_wait = min(time_to_wait * 2, time::Duration::from_millis(1000));
                let elapsed = timer.elapsed();
                if Instant::now() >= reported_time + time_to_wait {
                    eprint!(" Mapped {:12.3} M reads @ {:8.2} us/read                   \r", n_reads as f64 / 1E6, elapsed.as_secs_f64() * 1E6 / n_reads as f64);
                    reported_time = Instant::now();
                }
            }
            eprintln!(" Mapped {:12.3} M reads @ {:8.2} us/read                   \r", total_reads as f64 / 1E6, timer.elapsed().as_secs_f64() * 1E6 / total_reads as f64);
            eprintln!("Done!");
        });
    }

    // worker threads
    thread::scope(|s| {
        for _ in 0..args.threads {
            let mut mapper = mapper.clone();
            let out_tx = out_tx.clone();
            let chunks_rx = chunks_rx.clone();
            let abundances_tx = abundances_tx.clone();
            s.spawn(move || {
                loop {
                    let msg = chunks_rx.lock().unwrap().recv();
                    if let Ok((index, chunk)) = msg {
                        let mapped_chunk = mapper.map_chunk(chunk).unwrap();
                        out_tx.send((index, mapped_chunk)).unwrap();
                    } else {
                        break;
                    }
                }
                drop(chunks_rx);
                drop(out_tx);
                if mode == Mode::Abundances {
                    abundances_tx.send(mapper.abundances).unwrap();
                }
            });
        }
    });
    drop(out_tx);
    reader_thread.join().unwrap();
    let details = details_rx.recv().unwrap();
    let out = writer_thread.join().unwrap();

    let mut out = Arc::into_inner(out).unwrap().into_inner().unwrap();

    if mode == Mode::Abundances {
        let mut abundances = vec![0f64; references.len()];
        for _ in 0..args.threads {
            let worker_abundances = abundances_rx.recv().unwrap();
            for i in 0..abundances.len() {
                abundances[i] += worker_abundances[i];
            }
        }
        output_abundances(&abundances, &references, &mut out)?;
    }

    debug!("");
    debug!("# Statistics");
    debug!("");
    debug!("Reads:                                    {:12}", details.nam.n_reads);
    debug!("");
    debug!("## Randstrobe lookup (without rescue)");
    debug!("");
    debug!("Randstrobes                               {:12}     100.0 %   Per read: {:7.1}", details.nam.n_randstrobes, details.nam.n_randstrobes as f64 / details.nam.n_reads as f64);
    debug!("  Full randstrobe found                   {:12} {:9.1} %", details.nam.hits.full_found, 100f64 * details.nam.hits.full_found as f64 / details.nam.n_randstrobes as f64);
    debug!("  Full randstrobe found but filtered      {:12} {:9.1} %", details.nam.hits.full_filtered, 100f64 * details.nam.hits.full_filtered as f64 / details.nam.n_randstrobes as f64);
    debug!("  Full randstrobe not found               {:12} {:9.1} %", details.nam.hits.full_not_found, 100f64 * details.nam.hits.full_not_found as f64 / details.nam.n_randstrobes as f64);
    debug!("    Partial randstrobe found              {:12} {:9.1} %", details.nam.hits.partial_found, 100f64 * details.nam.hits.partial_found as f64 / details.nam.n_randstrobes as f64);
    debug!("    Partial randstrobe found but filtered {:12} {:9.1} %", details.nam.hits.partial_filtered, 100f64 * details.nam.hits.partial_filtered as f64 / details.nam.n_randstrobes as f64);
    debug!("    Partial randstrobe not found          {:12} {:9.1} %", details.nam.hits.partial_not_found, 100f64 * details.nam.hits.partial_not_found as f64 / details.nam.n_randstrobes as f64);
    debug!("");
    debug!("## Chaining");
    debug!("");
    debug!("Found anchors:                            {:12}               Per read: {:7.1}", details.nam.n_anchors, details.nam.n_anchors as f64 / details.nam.n_reads as f64);
    debug!("Found chains:                             {:12}               Per read: {:7.1}", details.nam.n_nams, details.nam.n_nams as f64 / details.nam.n_reads as f64);
    debug!("");
    debug!("## Rescue (-R)");
    debug!("");
    debug!("Rescue attempts: {:12}", details.nam.nam_rescue);
    debug!("Rescue hits:     {:12}", details.nam.n_rescue_hits);
    debug!("Rescued chains:  {:12}", details.nam.n_rescue_nams);
    debug!("");
    debug!("## Other");
    debug!("");
    debug!("Total mapping sites tried: {}", details.tried_alignment);
    debug!("Inconsistent NAM ends: {}", details.inconsistent_nams);
    debug!("Mates rescued by alignment: {}", details.mate_rescue);
    // Single-threaded:
    // for chunk in chunks_iter {
    //     mapper.map_chunk(&mut out, &mut rng, chunk)?;
    // }
    /*
    TODO
    let mut out = out.clone().lock().unwrap();;
    if mode == Mode::Abundances {
        mapper.output_abundances(&mut out)?;
    }
    */
    // TODO out.lock().unwrap().flush()?;

    info!("Total time mapping: {:.2} s", timer.elapsed().as_secs_f64());
    //info!("Total time reading read-file(s): {:.2} s", );
    info!("Total time creating strobemers: {:.2} s", details.nam.time_randstrobes);
    info!("Total time finding hits (non-rescue mode): {:.2} s", details.nam.time_find_hits);
    info!("Total time finding hits (rescue mode): {:.2} s", details.nam.time_rescue);
    info!("Total time chaining (non-rescue mode): {:.2} s", details.nam.time_chaining);
    info!("Total time sorting NAMs/chains by score: {:.2} s", details.nam.time_sort_nams);
    info!("Total time extending and pairing seeds: {:.2} s", details.time_extend);

    Ok(())
}

#[derive(PartialEq, Debug, Clone, Copy)]
enum Mode {
    Sam, Paf, Abundances
}

#[derive(Clone)]
struct Mapper<'a> {
    index: &'a StrobemerIndex<'a>,
    references: &'a [RefSequence],
    mapping_parameters: &'a MappingParameters,
    index_parameters: &'a IndexParameters,
    sam_output: &'a SamOutput,
    aligner: Aligner,
    chainer: &'a Chainer,
    include_unmapped: bool,
    mode: Mode,
    abundances: Vec<f64>,
}

impl Mapper<'_> {
    fn map_chunk(
        &mut self, // TODO only because of abundances
        chunk: Vec<io::Result<(SequenceRecord, Option<SequenceRecord>)>>,
    ) -> io::Result<(Vec<u8>, Details)> {
        let mut out = vec![];
        let mut rng = Rng::with_seed(0);
        let mut isizedist = InsertSizeDistribution::new();
        let mut cumulative_details = Details::new();
        for record in chunk {
            let (r1, r2) = record?;
            trace!("\nQuery: {}", r1.name);
            match self.mode {
                Mode::Sam => {
                    let (sam_records, details) =
                        if let Some(r2) = r2 {
                            let (mut records, details) = align_paired_end_read(
                                &r1, &r2, self.index, self.references, self.mapping_parameters, self.sam_output, self.index_parameters, &mut isizedist, &self.chainer, &self.aligner, &mut rng
                            );
                            if !self.include_unmapped && !records[0].is_mapped() && !records[1].is_mapped() {
                                records = vec![];
                            }

                            (records, details)
                        } else {
                            let (mut records, details) =
                                align_single_end_read(&r1, self.index, self.references, self.mapping_parameters, self.sam_output, &self.chainer, &self.aligner, &mut rng);
                            if !self.include_unmapped && !records[0].is_mapped() {
                                records = vec![];
                            }

                            (records, details)
                        };

                    for sam_record in sam_records {
                        writeln!(out, "{}", sam_record)?;
                    }
                    cumulative_details += details;
                }
                Mode::Paf => {
                    let (paf_records, details) =
                        if let Some(r2) = r2 {
                            map_paired_end_read(
                                &r1, &r2, self.index, self.references, self.mapping_parameters.rescue_level, &mut isizedist, self.mapping_parameters.mcs_strategy, &self.chainer, &mut rng
                            )
                        } else {
                            map_single_end_read(&r1, self.index, self.references, self.mapping_parameters.rescue_level, self.mapping_parameters.mcs_strategy, &self.chainer, &mut rng)
                        };
                    for paf_record in paf_records {
                        writeln!(out, "{}", paf_record)?;
                    }
                    cumulative_details += details;
                }
                Mode::Abundances => {
                    if let Some(r2) = r2 {
                        abundances_paired_end_read(&r1, &r2, self.index, &mut self.abundances, self.mapping_parameters.rescue_level, &mut isizedist, self.mapping_parameters.mcs_strategy, &self.chainer, &mut rng);
                    } else {
                        abundances_single_end_read(&r1, self.index, &mut self.abundances, self.mapping_parameters.rescue_level, self.mapping_parameters.mcs_strategy, &self.chainer, &mut rng);
                    }
                }
            }
        }
        Ok((out, cumulative_details))
    }
}

pub fn output_abundances<T: Write>(abundances: &[f64], references: &[RefSequence], out: &mut T) -> io::Result<()> {
    for i in 0..references.len() {
        let normalized = abundances[i] / references[i].sequence.len() as f64;
        writeln!(out, "{}\t{:.6}", references[i].name, normalized)?;
    }
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
