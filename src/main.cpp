#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <thread>
#include <assert.h>
#include <math.h>
#include <inttypes.h>
#include <iomanip>

#include "args.hxx"
#include "robin_hood.h"
#include "ssw_cpp.h"
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


static Logger& logger = Logger::get();


static void print_diagnostics(const StrobemerIndex &index, const std::string& logfile_name, int k) {
    // Prins to csv file the statistics on the number of seeds of a particular length and what fraction of them them are unique in the index:
    // format:
    // seed_length, count, percentage_unique

    size_t max_size = 100000;
    std::vector<int> log_count(max_size, 0);  // stores count and each index represents the length
    std::vector<int> log_unique(max_size, 0);  // stores count unique and each index represents the length
    std::vector<int> log_repetitive(max_size, 0);  // stores count unique and each index represents the length


    std::vector<uint64_t> log_count_squared(max_size,0);
    uint64_t tot_seed_count = 0;
    uint64_t tot_seed_count_sq = 0;

    std::vector<uint64_t> log_count_1000_limit(max_size, 0);  // stores count and each index represents the length
    uint64_t tot_seed_count_1000_limit = 0;

    size_t seed_length;
    for (auto &it : index.mers_index) {
        auto ref_mer = it.second;
        auto offset = ref_mer.offset;
        auto count = ref_mer.count;

        for (size_t j = offset; j < offset + count; ++j) {
            auto r = index.flat_vector[j];
            auto p = r.packed;
            const uint8_t bit_alloc = 8;
            const int mask = (1 << bit_alloc) - 1;
            const int offset = p & mask;
            seed_length =  offset + k;
            if (seed_length < max_size){
                log_count[seed_length] ++;
                log_count_squared[seed_length] += count;
                tot_seed_count ++;
                tot_seed_count_sq += count;
                if (count <= 1000){
                    log_count_1000_limit[seed_length] ++;
                    tot_seed_count_1000_limit ++;
                }
            } else {
               logger.info() << "Detected seed size over " << max_size << " bp (can happen, e.g., over centromere): " << seed_length << std::endl;
            }
        }

        if (count == 1 && seed_length < max_size) {
            log_unique[seed_length]++;
        }
        if (count >= 10 && seed_length < max_size) {
            log_repetitive[seed_length]++;
        }
    }

    // printing
    std::ofstream log_file;
    log_file.open(logfile_name);

    for (size_t i = 0; i < log_count.size(); ++i) {
        if (log_count[i] > 0) {
            double e_count = log_count_squared[i] / log_count[i];
            log_file << i << ',' << log_count[i] << ',' << e_count << std::endl;
        }
    }

    // Get median
    size_t n = 0;
    int median = 0;
    for (size_t i = 0; i < log_count.size(); ++i) {
        n += log_count[i];
        if (n >= tot_seed_count/2) {
            break;
        }
    }
    // Get median 1000 limit
    size_t n_lim = 0;
    for (size_t i = 0; i < log_count_1000_limit.size(); ++i) {
        n_lim += log_count_1000_limit[i];
        if (n_lim >= tot_seed_count_1000_limit/2) {
            break;
        }
    }

    log_file << "E_size for total seeding wih max seed size m below (m, tot_seeds, E_hits)" << std::endl;
    double e_hits = (double) tot_seed_count_sq/ (double) tot_seed_count;
    double fraction_masked = 1.0 - (double) tot_seed_count_1000_limit/ (double) tot_seed_count;
    log_file << median << ',' << tot_seed_count << ',' << e_hits << ',' << 100*fraction_masked << std::endl;
}

/*
 * Return formatted SAM header as a string
 */
std::string sam_header(const References& references, const std::string& read_group_id, const std::vector<std::string>& read_group_fields) {
    std::stringstream out;
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
    out << "@PG\tID:strobealign\tPN:strobealign\tVN:" << version_string() << "\tCL:strobealign\n";
    return out.str();
}

int run_strobealign(int argc, char **argv) {
    CommandLineOptions opt;
    mapping_params map_param;
    std::tie(opt, map_param) = parse_command_line_arguments(argc, argv);

    logger.set_level(opt.verbose ? LOG_DEBUG : LOG_INFO);
    logger.info() << std::setprecision(2) << std::fixed;
    logger.info() << "This is strobealign " << version_string() << '\n';

    if (!opt.r_set && !opt.reads_filename1.empty()) {
        map_param.r = estimate_read_length(opt.reads_filename1, opt.reads_filename2);
        logger.info() << "Estimated read length: " << map_param.r << " bp\n";
    }
    if (opt.c >= 64 || opt.c <= 0) {
        throw BadParameter("c must be greater than 0 and less than 64");
    }
    IndexParameters index_parameters = IndexParameters::from_read_length(
        map_param.r, opt.c, opt.k_set ? opt.k : -1, opt.s_set ? opt.s : -1, opt.max_seed_len_set ? opt.max_seed_len : -1);

    alignment_params aln_params;
    aln_params.match = opt.A;
    aln_params.mismatch = opt.B;
    aln_params.gap_open = opt.O;
    aln_params.gap_extend = opt.E;

    logger.debug() << "Using" << std::endl
        << "k: " << index_parameters.k << std::endl
        << "s: " << index_parameters.s << std::endl
        << "w_min: " << index_parameters.w_min << std::endl
        << "w_max: " << index_parameters.w_max << std::endl
        << "Read length (r): " << map_param.r << std::endl
        << "Maximum seed length: " << index_parameters.max_dist + index_parameters.k << std::endl
        << "Threads: " << opt.n_threads << std::endl
        << "R: " << map_param.R << std::endl
        << "Expected [w_min, w_max] in #syncmers: [" << index_parameters.w_min << ", " << index_parameters.w_max << "]" << std::endl
        << "Expected [w_min, w_max] in #nucleotides: [" << (index_parameters.k - index_parameters.s + 1) * index_parameters.w_min << ", " << (index_parameters.k - index_parameters.s + 1) * index_parameters.w_max << "]" << std::endl
        << "A: " << opt.A << std::endl
        << "B: " << opt.B << std::endl
        << "O: " << opt.O << std::endl
        << "E: " << opt.E << std::endl;

    map_param.verify();
    index_parameters.verify();

//    assert(k <= (w/2)*w_min && "k should be smaller than (w/2)*w_min to avoid creating short strobemers");

    // Create index
    References references;
    Timer read_refs_timer;
    references = References::from_fasta(opt.ref_filename);
    logger.info() << "Time reading references: " << read_refs_timer.elapsed() << " s\n";

    if (references.total_length() == 0) {
        throw InvalidFasta("No reference sequences found");
    }

    StrobemerIndex index(references, index_parameters);
    if (opt.use_index) {
        // Read the index from a file
        assert(!opt.only_gen_index);
        Timer read_index_timer;
        index.read(opt.ref_filename + ".sti");
        logger.info() << "Total time reading index: " << read_index_timer.elapsed() << " s\n";
    } else {
        logger.info() << "Indexing ...\n";
        Timer index_timer;
        index.populate(opt.f);
        
        logger.info() << "  Time generating seeds: " << index.stats.elapsed_generating_seeds.count() << " s" <<  std::endl;
        logger.info() << "  Time sorting seeds: " << index.stats.elapsed_sorting_seeds.count() << " s" <<  std::endl;
        logger.info() << "  Time generating flat vector: " << index.stats.elapsed_flat_vector.count() << " s" <<  std::endl;
        logger.info() << "  Time generating hash table index: " << index.stats.elapsed_hash_index.count() << " s" <<  std::endl;
        logger.info() << "Total time indexing: " << index_timer.elapsed() << " s\n";

        logger.debug()
        << "Unique strobemers: " << index.stats.unique_mers << std::endl
        << "Total strobemers count: " << index.stats.tot_strobemer_count << std::endl
        << "Total strobemers occur once: " << index.stats.tot_occur_once << std::endl
        << "Fraction Unique: " << index.stats.frac_unique << std::endl
        << "Total strobemers highly abundant > 100: " << index.stats.tot_high_ab << std::endl
        << "Total strobemers mid abundance (between 2-100): " << index.stats.tot_mid_ab << std::endl
        << "Total distinct strobemers stored: " << index.stats.tot_distinct_strobemer_count << std::endl;
        if (index.stats.tot_high_ab >= 1) {
            logger.debug() << "Ratio distinct to highly abundant: " << index.stats.tot_distinct_strobemer_count / index.stats.tot_high_ab << std::endl;
        }
        if (index.stats.tot_mid_ab >= 1) {
            logger.debug() << "Ratio distinct to non distinct: " << index.stats.tot_distinct_strobemer_count / (index.stats.tot_high_ab + index.stats.tot_mid_ab) << std::endl;
        }
        logger.debug() << "Filtered cutoff index: " << index.stats.index_cutoff << std::endl;
        logger.debug() << "Filtered cutoff count: " << index.stats.filter_cutoff << std::endl;
        
        if (!opt.logfile_name.empty()) {
            print_diagnostics(index, opt.logfile_name, index_parameters.k);
            logger.debug() << "Finished printing log stats" << std::endl;
        }
        if (opt.only_gen_index) {
            Timer index_writing_timer;
            index.write(opt.ref_filename + ".sti");
            logger.info() << "Total time writing index: " << index_writing_timer.elapsed() << " s\n";
            return EXIT_SUCCESS;
        }
    }

    // Map/align reads
        
    Timer map_align_timer;
    map_param.rescue_cutoff = map_param.R < 100 ? map_param.R * index.filter_cutoff : 1000;
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

    if (map_param.is_sam_out) {
        out << sam_header(references, opt.read_group_id, opt.read_group_fields);
    }

    std::vector<AlignmentStatistics> log_stats_vec(opt.n_threads);

    logger.info() << "Running in " << (opt.is_SE ? "single-end" : "paired-end") << " mode" << std::endl;

    if (opt.is_SE) {
        auto ks = open_fastq(opt.reads_filename1);

        InputBuffer input_buffer(ks, ks, opt.chunk_size);
        OutputBuffer output_buffer(out);

        std::vector<std::thread> workers;
        for (int i = 0; i < opt.n_threads; ++i) {
            std::thread consumer(perform_task_SE, std::ref(input_buffer), std::ref(output_buffer),
                std::ref(log_stats_vec[i]), std::ref(aln_params),
                std::ref(map_param), std::ref(index_parameters), std::ref(references),
                std::ref(index), std::ref(opt.read_group_id));
            workers.push_back(std::move(consumer));
        }

        for (size_t i = 0; i < workers.size(); ++i) {
            workers[i].join();
        }
    }
    else {
        auto ks1 = open_fastq(opt.reads_filename1);
        auto ks2 = open_fastq(opt.reads_filename2);

        InputBuffer input_buffer(ks1, ks2, opt.chunk_size);
        OutputBuffer output_buffer(out);

        std::vector<std::thread> workers;
        for (int i = 0; i < opt.n_threads; ++i) {
            std::thread consumer(perform_task_PE, std::ref(input_buffer), std::ref(output_buffer),
                std::ref(log_stats_vec[i]), std::ref(aln_params),
                std::ref(map_param), std::ref(index_parameters), std::ref(references),
                std::ref(index), std::ref(opt.read_group_id));
            workers.push_back(std::move(consumer));
        }

        for (size_t i = 0; i < workers.size(); ++i) {
            workers[i].join();
        }
    }
    logger.info() << "Done!\n";

    AlignmentStatistics tot_statistics;
    for (auto& it : log_stats_vec) {
        tot_statistics += it;
    }

    logger.info() << "Total mapping sites tried: " << tot_statistics.tot_all_tried << std::endl
        << "Total calls to ssw: " << tot_statistics.tot_ksw_aligned << std::endl
        << "Calls to ksw (rescue mode): " << tot_statistics.tot_rescued << std::endl
        << "Did not fit strobe start site: " << tot_statistics.did_not_fit << std::endl
        << "Tried rescue: " << tot_statistics.tried_rescue << std::endl
        << "Total time mapping: " << map_align_timer.elapsed() << " s." << std::endl
        << "Total time reading read-file(s): " << tot_statistics.tot_read_file.count() / opt.n_threads << " s." << std::endl
        << "Total time creating strobemers: " << tot_statistics.tot_construct_strobemers.count() / opt.n_threads << " s." << std::endl
        << "Total time finding NAMs (non-rescue mode): " << tot_statistics.tot_find_nams.count() / opt.n_threads << " s." << std::endl
        << "Total time finding NAMs (rescue mode): " << tot_statistics.tot_time_rescue.count() / opt.n_threads << " s." << std::endl;
    //<< "Total time finding NAMs ALTERNATIVE (candidate sites): " << tot_find_nams_alt.count()/opt.n_threads  << " s." <<  std::endl;
    logger.info() << "Total time sorting NAMs (candidate sites): " << tot_statistics.tot_sort_nams.count() / opt.n_threads << " s." << std::endl
        << "Total time base level alignment (ssw): " << tot_statistics.tot_extend.count() / opt.n_threads << " s." << std::endl
        << "Total time writing alignment to files: " << tot_statistics.tot_write_file.count() << " s." << std::endl;
    return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
    try {
        return run_strobealign(argc, argv);
    } catch (BadParameter& e) {
        logger.error() << "A mapping or seeding parameter is invalid: " << e.what() << std::endl;
    } catch (const std::runtime_error& e) {
        logger.error() << "strobealign: " << e.what() << std::endl;
    }
    return EXIT_FAILURE;
}
