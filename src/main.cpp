#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <thread>
#include <cassert>
#include <iomanip>
#include <chrono>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include "refs.hpp"
#include "exceptions.hpp"
#include "cmdline.hpp"
#include "index.hpp"
#include "pc.hpp"
#include "aln.hpp"
#include "logger.hpp"
#include "timer.hpp"
#include "readlen.hpp"
#include "version.hpp"
#include "randstrobes.hpp"
#include "buildconfig.hpp"


static Logger& logger = Logger::get();

/*
 * Return formatted SAM header as a string
 */
std::string sam_header(
    const References& references,
    const std::string& read_group_id,
    const std::vector<std::string>& read_group_fields
) {
    std::stringstream out;
    out << "@HD\tVN:1.6\tSO:unsorted\n";
    for (size_t i = 0; i < references.size(); ++i) {
        out << "@SQ\tSN:" << references.names[i] << "\tLN:" << references.lengths[i] << "\n";
    }
    if (!read_group_id.empty()) {
        out << "@RG\tID:" << read_group_id;
        for (const auto& field : read_group_fields) {
           out << '\t' << field;
        }
        out << '\n';
    }
    return out.str();
}

std::string pg_header(const std::string& cmd_line) {
    std::stringstream out;
    out << "@PG\tID:strobealign\tPN:strobealign\tVN:" << version_string() << "\tCL:" << cmd_line << std::endl;
    return out.str();
}

void warn_if_no_optimizations() {
    if (std::string(CMAKE_BUILD_TYPE) == "Debug") {
        logger.info() << "\n    ***** Binary was compiled without optimizations - this will be very slow *****\n\n";
    }
}

bool avx2_enabled() {
#ifdef __AVX2__
    return true;
#else
    return false;
#endif
}

InputBuffer get_input_buffer(const CommandLineOptions& opt) {
    if (opt.is_SE) {
        return InputBuffer(opt.reads_filename1, "", opt.chunk_size, false);
    } else if (opt.is_interleaved) {
        if (opt.reads_filename2 != "") {
            throw BadParameter("Cannot specify both --interleaved and specify two read files");
        }
        return InputBuffer(opt.reads_filename1, "", opt.chunk_size, true);
    } else {
        return InputBuffer(opt.reads_filename1, opt.reads_filename2, opt.chunk_size, false);
    }
}

void output_abundance(std::ostream& out, const std::vector<double>& abundances, const References& references){
        for (size_t i = 0; i < references.size(); ++i) {
            out << references.names[i] << '\t' << std::fixed << std::setprecision(6) << abundances[i] / double(references.sequences[i].size()) << std::endl;
        }
}

void show_progress_until_done(std::vector<int>& worker_done, std::vector<AlignmentStatistics>& stats) {
    Timer timer;
    bool reported = false;
    bool done = false;
    // Waiting time between progress updates
    // Start with a small value so that there’s no delay if there are very few
    // reads to align.
    auto time_to_wait = std::chrono::milliseconds(1);
    while (!done) {
        std::this_thread::sleep_for(time_to_wait);
        // Ramp up waiting time
        time_to_wait = std::min(time_to_wait * 2, std::chrono::milliseconds(1000));
        done = true;
        for (auto is_done : worker_done) {
            if (!is_done) {
                done = false;
                continue;
            }
        }
        auto n_reads = 0ull;
        for (auto& stat : stats) {
            n_reads += stat.n_reads;
        }
        auto elapsed = timer.elapsed();
        if (elapsed >= 1.0) {
            std::cerr
                << " Mapped "
                << std::setw(12) << (n_reads / 1E6) << " M reads @ "
                << std::setw(8) << (timer.elapsed() * 1E6 / n_reads) << " us/read                   \r";
            reported = true;
        }
    }
    if (reported) {
        std::cerr << '\n';
    }
}

int run_strobealign(int argc, char **argv) {
    auto opt = parse_command_line_arguments(argc, argv);

    logger.set_level(opt.verbose ? LOG_DEBUG : LOG_INFO);
    logger.info() << std::setprecision(2) << std::fixed;
    logger.info() << "This is strobealign " << version_string() << '\n';
    logger.debug() << "Build type: " << CMAKE_BUILD_TYPE << '\n';
    warn_if_no_optimizations();
    logger.debug() << "AVX2 enabled: " << (avx2_enabled() ? "yes" : "no") << '\n';

    if (opt.c >= 64 || opt.c <= 0) {
        throw BadParameter("c must be greater than 0 and less than 64");
    }

    if (!opt.is_sam_out && opt.is_abundance_out){
        throw BadParameter("Can not use -x and --aemb at the same time");
    }

    InputBuffer input_buffer = get_input_buffer(opt);
    if (!opt.r_set && !opt.reads_filename1.empty()) {
        opt.r = estimate_read_length(input_buffer);
        logger.info() << "Estimated read length: " << opt.r << " bp\n";
        input_buffer.rewind_reset();
    }
    IndexParameters index_parameters = IndexParameters::from_read_length(
        opt.r,
        opt.k_set ? opt.k : IndexParameters::DEFAULT,
        opt.s_set ? opt.s : IndexParameters::DEFAULT,
        opt.l_set ? opt.l : IndexParameters::DEFAULT,
        opt.u_set ? opt.u : IndexParameters::DEFAULT,
        opt.c_set ? opt.c : IndexParameters::DEFAULT,
        opt.max_seed_len_set ? opt.max_seed_len : IndexParameters::DEFAULT,
        opt.aux_len ? opt.aux_len : IndexParameters::DEFAULT
    );
    AlignmentParameters aln_params;
    aln_params.match = opt.A;
    aln_params.mismatch = opt.B;
    aln_params.gap_open = opt.O;
    aln_params.gap_extend = opt.E;
    aln_params.end_bonus = opt.end_bonus;

    MappingParameters map_param;
    map_param.r = opt.r;
    map_param.max_secondary = opt.max_secondary;
    map_param.dropoff_threshold = opt.dropoff_threshold;
    map_param.rescue_level = opt.rescue_level;
    map_param.max_tries = opt.max_tries;
    map_param.use_mcs = opt.mcs;
    map_param.output_format = (
            opt.is_abundance_out ? OutputFormat::Abundance :
            opt.is_sam_out ? OutputFormat::SAM :
                OutputFormat::PAF);
    map_param.cigar_ops = opt.cigar_eqx ? CigarOps::EQX : CigarOps::M;
    map_param.output_unmapped = opt.output_unmapped;
    map_param.details = opt.details;
    map_param.fastq_comments = opt.fastq_comments;
    map_param.verify();

    logger.debug() << index_parameters << '\n';
    logger.debug()
        << "  Maximum seed length: " << index_parameters.randstrobe.max_dist + index_parameters.syncmer.k << '\n'
        << "  Expected [w_min, w_max] in #syncmers: [" << index_parameters.randstrobe.w_min << ", " << index_parameters.randstrobe.w_max << "]\n"
        << "  Expected [w_min, w_max] in #nucleotides: [" << (index_parameters.syncmer.k - index_parameters.syncmer.s + 1) * index_parameters.randstrobe.w_min << ", " << (index_parameters.syncmer.k - index_parameters.syncmer.s + 1) * index_parameters.randstrobe.w_max << "]\n";
    logger.debug() << aln_params << '\n';
    logger.debug() << "Rescue level (R): " << map_param.rescue_level << '\n';
    logger.debug() << "Threads: " << opt.n_threads << std::endl;

//    assert(k <= (w/2)*w_min && "k should be smaller than (w/2)*w_min to avoid creating short strobemers");

    // Create index
    References references;
    Timer read_refs_timer;
    references = References::from_fasta(opt.ref_filename);
    logger.info() << "Time reading reference: " << read_refs_timer.elapsed() << " s\n";

    logger.info() << "Reference size: " << references.total_length() / 1E6 << " Mbp ("
        << references.size() << " contig" << (references.size() == 1 ? "" : "s")
        << "; largest: "
        << (*std::max_element(references.lengths.begin(), references.lengths.end()) / 1E6) << " Mbp)\n";
    if (references.total_length() == 0) {
        throw InvalidFasta("No reference sequences found");
    }

    if (references.size() > RefRandstrobe::max_number_of_references) {
        throw InvalidFasta("Too many reference sequences. Current maximum is " + std::to_string(RefRandstrobe::max_number_of_references));
    }

    logger.debug() << "Auxiliary hash length: " << opt.aux_len << "\n";
    logger.info() << "Using multi-context seeds: " << (map_param.use_mcs ? "yes" : "no") << '\n';
    StrobemerIndex index(references, index_parameters, opt.bits);
    if (opt.use_index) {
        // Read the index from a file
        assert(!opt.only_gen_index);
        Timer read_index_timer;
        std::string sti_path = opt.ref_filename + index_parameters.filename_extension();
        logger.info() << "Reading index from " << sti_path << '\n';
        index.read(sti_path);
        logger.debug() << "Bits used to index buckets: " << index.get_bits() << "\n";
        logger.info() << "Total time reading index: " << read_index_timer.elapsed() << " s\n";
    } else {
        logger.debug() << "Bits used to index buckets: " << index.get_bits() << "\n";
        logger.info() << "Indexing ...\n";
        Timer index_timer;
        index.populate(opt.f, opt.n_threads);
        
        logger.info() << "  Time counting seeds: " << index.stats.elapsed_counting_hashes.count() << " s" <<  std::endl;
        logger.info() << "  Time generating seeds: " << index.stats.elapsed_generating_seeds.count() << " s" <<  std::endl;
        logger.info() << "  Time sorting seeds: " << index.stats.elapsed_sorting_seeds.count() << " s" <<  std::endl;
        logger.info() << "  Time generating hash table index: " << index.stats.elapsed_hash_index.count() << " s" <<  std::endl;
        logger.info() << "Total time indexing: " << index_timer.elapsed() << " s\n";

        logger.debug()
            << "Index statistics\n"
            << "  Total strobemers:    " << std::setw(14) << index.stats.tot_strobemer_count << '\n'
            << "  Distinct strobemers: " << std::setw(14) << index.stats.distinct_strobemers << " (100.00%)\n"
            << "    1 occurrence:      " << std::setw(14) << index.stats.tot_occur_once
                << " (" << std::setw(6) << (100.0 * index.stats.tot_occur_once / index.stats.distinct_strobemers) << "%)\n"
            << "    2..100 occurrences:" << std::setw(14) << index.stats.tot_mid_ab
                << " (" << std::setw(6) << (100.0 * index.stats.tot_mid_ab / index.stats.distinct_strobemers) << "%)\n"
            << "    >100 occurrences:  " << std::setw(14) << index.stats.tot_high_ab
                << " (" << std::setw(6) << (100.0 * index.stats.tot_high_ab / index.stats.distinct_strobemers) << "%)\n"
            ;
        if (index.stats.tot_high_ab >= 1) {
            logger.debug() << "Ratio distinct to highly abundant: " << index.stats.distinct_strobemers / index.stats.tot_high_ab << std::endl;
        }
        if (index.stats.tot_mid_ab >= 1) {
            logger.debug() << "Ratio distinct to non distinct: " << index.stats.distinct_strobemers / (index.stats.tot_high_ab + index.stats.tot_mid_ab) << std::endl;
        }
        logger.debug() << "Filtered cutoff index: " << index.stats.index_cutoff << std::endl;
        logger.debug() << "Filtered cutoff count: " << index.stats.filter_cutoff << std::endl;
        
        if (!opt.logfile_name.empty()) {
            index.print_diagnostics(opt.logfile_name, index_parameters.syncmer.k);
            logger.debug() << "Finished printing log stats" << std::endl;
        }
        if (opt.only_gen_index) {
            Timer index_writing_timer;
            std::string sti_path = opt.ref_filename + index_parameters.filename_extension();
            logger.info() << "Writing index to " << sti_path << '\n';
            index.write(opt.ref_filename + index_parameters.filename_extension());
            logger.info() << "Total time writing index: " << index_writing_timer.elapsed() << " s\n";
            return EXIT_SUCCESS;
        }
    }

    // Map/align reads
        
    Timer map_align_timer;
    map_param.rescue_cutoff = map_param.rescue_level < 100 ? map_param.rescue_level * index.filter_cutoff : 1000;
    logger.debug() << "Using rescue cutoff: " << map_param.rescue_cutoff << std::endl;

    std::streambuf* buf;
    std::ofstream of;

    if (!opt.write_to_stdout) {
        of.open(opt.output_file_name);
        buf = of.rdbuf();
    }
    else {
        buf = std::cout.rdbuf();
    }

    std::ostream out(buf);
    
    if (map_param.output_format == OutputFormat::SAM) {
            std::stringstream cmd_line;
            for(int i = 0; i < argc; ++i) {
                cmd_line << argv[i] << " ";
            }
            out << sam_header(references, opt.read_group_id, opt.read_group_fields);
            if (opt.pg_header) {
                out << pg_header(cmd_line.str());
            }
    }

    std::vector<AlignmentStatistics> worker_statistics(opt.n_threads);
    
    logger.info() << "Processing " << (opt.is_SE ? "single-end" : "paired-end") << " reads ";
    if (map_param.output_format == OutputFormat::PAF) {
        logger.info() << "in mapping-only mode ";
    }
    if (map_param.output_format == OutputFormat::Abundance) {
        logger.info() << "in abundance estimation mode ";
    }
    logger.info() << "using " << opt.n_threads << " thread" << (opt.n_threads != 1 ? "s" : "") << std::endl;

    OutputBuffer output_buffer(out);
    std::vector<std::thread> workers;
    std::vector<int> worker_done(opt.n_threads);  // each thread sets its entry to 1 when it’s done
    std::vector<std::vector<double>> worker_abundances(opt.n_threads, std::vector<double>(references.size(), 0));
    for (int i = 0; i < opt.n_threads; ++i) {
        std::thread consumer(perform_task, std::ref(input_buffer), std::ref(output_buffer),
            std::ref(worker_statistics[i]), std::ref(worker_done[i]), std::ref(aln_params),
            std::ref(map_param), std::ref(index_parameters), std::ref(references),
            std::ref(index), std::ref(opt.read_group_id), std::ref(worker_abundances[i]));
        workers.push_back(std::move(consumer));
    }
    if (opt.show_progress && isatty(2)) {
        show_progress_until_done(worker_done, worker_statistics);
    }
    for (auto& worker : workers) {
        worker.join();
    }
    logger.info() << "Done!\n";

    AlignmentStatistics statistics;
    for (auto& it : worker_statistics) {
        statistics += it;
    }

    if (map_param.output_format == OutputFormat::Abundance) {
        std::vector<double> abundances(references.size(), 0);
        for (size_t i = 0; i < worker_abundances.size(); ++i) {
            for (size_t j = 0; j < worker_abundances[i].size(); ++j) {
                abundances[j] += worker_abundances[i][j];
            }
        }
        output_abundance(out, abundances, references);
    }

    logger.debug()
        << "Number of reads:               " << std::setw(12) << statistics.n_reads << std::endl
        << "Number of randstrobes:         " << std::setw(12) << statistics.n_randstrobes
        << "  Per read: " << std::setw(7) << static_cast<float>(statistics.n_randstrobes) / statistics.n_reads << std::endl
        << "Number of partial hits:        " << std::setw(12) << statistics.n_partial_hits << '\n'
        << "Number of non-rescue hits:     " << std::setw(12) << statistics.n_hits
        << "  Per read: " << std::setw(7) << static_cast<float>(statistics.n_hits) / statistics.n_reads << std::endl
        << "Number of non-rescue NAMs:     " << std::setw(12) << statistics.n_nams
        << "  Per read: " << std::setw(7) << static_cast<float>(statistics.n_nams) / statistics.n_reads << std::endl
        << "Number of NAM rescue attempts: " << std::setw(12) << statistics.nam_rescue << std::endl
        << "Number of rescue hits:         " << std::setw(12) << statistics.n_rescue_hits
        << "  Per rescue attempt: " << std::setw(7) << static_cast<float>(statistics.n_rescue_hits) / statistics.nam_rescue << std::endl
        << "Number of rescue NAMs:         " << std::setw(12) << statistics.n_rescue_nams
        << "  Per rescue attempt: " << std::setw(7) << static_cast<float>(statistics.n_rescue_nams) / statistics.nam_rescue << std::endl;
    logger.info()
        << "Total mapping sites tried: " << statistics.tried_alignment << std::endl
        << "Total calls to ssw: " << statistics.tot_aligner_calls << std::endl
        << "Inconsistent NAM ends: " << statistics.inconsistent_nams << std::endl
        << "Mates rescued by alignment: " << statistics.tot_rescued << std::endl
        << "Total time mapping: " << map_align_timer.elapsed() << " s." << std::endl
        << "Total time reading read-file(s): " << statistics.tot_read_file.count() / opt.n_threads << " s." << std::endl
        << "Total time creating strobemers: " << statistics.tot_construct_strobemers.count() / opt.n_threads << " s." << std::endl
        << "Total time finding NAMs (non-rescue mode): " << statistics.tot_find_nams.count() / opt.n_threads << " s." << std::endl
        << "Total time finding NAMs (rescue mode): " << statistics.tot_time_rescue.count() / opt.n_threads << " s." << std::endl
        << "Total time sorting NAMs (candidate sites): " << statistics.tot_sort_nams.count() / opt.n_threads << " s." << std::endl
        << "Total time extending and pairing seeds: " << statistics.tot_extend.count() / opt.n_threads << " s." << std::endl;
    return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
    try {
        return run_strobealign(argc, argv);
    } catch (BadParameter& e) {
        logger.error() << "A parameter is invalid: " << e.what() << std::endl;
    } catch (const std::runtime_error& e) {
        logger.error() << "strobealign: " << e.what() << std::endl;
    }
    return EXIT_FAILURE;
}
