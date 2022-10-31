#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <chrono>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <thread>
#include <assert.h>
#include <math.h>
#include <inttypes.h>

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
#include "readlen.hpp"
#include "version.hpp"

using std::chrono::high_resolution_clock;


static Logger& logger = Logger::get();


static void print_diagnostics(const StrobemerIndex &index, const std::string& logfile_name, int k) {
    // Prins to csv file the statistics on the number of seeds of a particular length and what fraction of them them are unique in the index:
    // format:
    // seed_length, count, percentage_unique
    //

    int max_size = 100000;
    std::vector<int> log_count(max_size,0);  // stores count and each index represents the length
    std::vector<int> log_unique(max_size,0); // stores count unique and each index represents the length
    std::vector<int> log_repetitive(max_size,0); // stores count unique and each index represents the length


    std::vector<uint64_t> log_count_squared(max_size,0);
    uint64_t tot_seed_count = 0;
    uint64_t tot_seed_count_sq = 0;

    std::vector<int> log_count_1000_limit(max_size,0);  // stores count and each index represents the length
    uint64_t tot_seed_count_1000_limit = 0;
    uint64_t tot_seed_count_sq_1000_limit = 0;

    int seed_length;
    for (auto &it : index.mers_index) {
        auto ref_mer = it.second;
        auto offset = ref_mer.offset;
        auto count = ref_mer.count;

        for (size_t j = offset; j < offset + count; ++j) {
            auto r = index.flat_vector[j];
            auto p = r.packed;
            int bit_alloc = 8;
            int mask=(1<<bit_alloc) - 1;
            int offset = (p & mask);
            seed_length =  offset + k;
            if (seed_length < max_size){
                log_count[seed_length] ++;
                log_count_squared[seed_length] += count;
                tot_seed_count ++;
                tot_seed_count_sq += count;
                if (count <= 1000){
                    log_count_1000_limit[seed_length] ++;
                    tot_seed_count_1000_limit ++;
                    tot_seed_count_sq_1000_limit += count;
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
    int n = 0;
    int median = 0;
    for (size_t i = 0; i < log_count.size(); ++i) {
        n += log_count[i];
        if (n >= tot_seed_count/2) {
            median = i;
            break;
        }
    }
    // Get median 1000 limit
    int n_lim = 0;
    int median_lim = 0;
    for (size_t i = 0; i < log_count_1000_limit.size(); ++i) {
        n_lim += log_count_1000_limit[i];
        if ( n_lim >= tot_seed_count_1000_limit/2){
            median_lim = i;
            break;
        }
    }

    log_file << "E_size for total seeding wih max seed size m below (m, tot_seeds, E_hits)" << std::endl;
    double e_hits = (double) tot_seed_count_sq/ (double) tot_seed_count;
    double fraction_masked = 1.0 - (double) tot_seed_count_1000_limit/ (double) tot_seed_count;
    log_file << median << ',' << tot_seed_count << ',' << e_hits << ',' << 100*fraction_masked << std::endl;

//    for (int i=0 ; i < log_count.size(); ++i) {
//        if (log_count[i] > 0) {
//            log_file << i << ',' << log_count[i] << ',' << (float) log_unique[i] / (float) log_count[i] << ',' << (float) log_repetitive[i] / (float) log_count[i] << std::endl;
//        }
//    }

}

/*
 * Return formatted SAM header as a string
 */
std::string sam_header(const References& references) {
    std::stringstream out;
    for (size_t i = 0; i < references.size(); ++i) {
        out << "@SQ\tSN:" << references.names[i] << "\tLN:" << references.lengths[i] << "\n";
    }
    out << "@PG\tID:strobealign\tPN:strobealign\tVN:" VERSION_STRING "\tCL:strobealign\n";
    return out.str();
}

int run_strobealign(int argc, char **argv) {
    CommandLineOptions opt;
    mapping_params map_param;
    std::tie(opt, map_param) = parse_command_line_arguments(argc, argv);

    logger.set_level(opt.verbose ? LOG_DEBUG : LOG_INFO);
    if (!opt.r_set && !opt.reads_filename1.empty()) {
        map_param.r = estimate_read_length(opt.reads_filename1, opt.reads_filename2);
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
    auto start_read_refs = high_resolution_clock::now();
    references = References::from_fasta(opt.ref_filename);
    std::chrono::duration<double> elapsed_read_refs = high_resolution_clock::now() - start_read_refs;
    logger.info() << "Time reading references: " << elapsed_read_refs.count() << " s" << std::endl;

    if (references.total_length() == 0) {
        throw InvalidFasta("No reference sequences found");
    }

    StrobemerIndex index(references, index_parameters);
    if (opt.use_index) {
        // Read the index from a file
        assert(!opt.only_gen_index);
        auto start_read_index = high_resolution_clock::now();
        index.read(opt.index_filename);
        std::chrono::duration<double> elapsed_read_index = high_resolution_clock::now() - start_read_index;
        logger.info() << "Total time reading index: " << elapsed_read_index.count() << " s\n" << std::endl;
    } else {
        // Generate the index
        auto start = high_resolution_clock::now();
        IndexCreationStatistics index_creation_stats = index.populate(opt.f);
        
        logger.info() << "Time copying flat vector: " << index_creation_stats.elapsed_copy_flat_vector.count() << " s" << std::endl;
        logger.debug() << "Unique strobemers: " << index_creation_stats.unique_mers << std::endl;
        logger.info() << "Total time generating flat vector: " << index_creation_stats.elapsed_flat_vector.count() << " s" <<  std::endl;

        logger.debug()
        << "Total strobemers count: " << index_creation_stats.tot_strobemer_count << std::endl
        << "Total strobemers occur once: " << index_creation_stats.tot_occur_once << std::endl
        << "Fraction Unique: " << index_creation_stats.frac_unique << std::endl
        << "Total strobemers highly abundant > 100: " << index_creation_stats.tot_high_ab << std::endl
        << "Total strobemers mid abundance (between 2-100): " << index_creation_stats.tot_mid_ab << std::endl
        << "Total distinct strobemers stored: " << index_creation_stats.tot_distinct_strobemer_count << std::endl;
        if (index_creation_stats.tot_high_ab >= 1) {
            logger.debug() << "Ratio distinct to highly abundant: " << index_creation_stats.tot_distinct_strobemer_count / index_creation_stats.tot_high_ab << std::endl;
        }
        if (index_creation_stats.tot_mid_ab >= 1) {
            logger.debug() << "Ratio distinct to non distinct: " << index_creation_stats.tot_distinct_strobemer_count / (index_creation_stats.tot_high_ab + index_creation_stats.tot_mid_ab) << std::endl;
        }
        logger.debug() << "Filtered cutoff index: " << index_creation_stats.index_cutoff << std::endl;
        logger.debug() << "Filtered cutoff count: " << index_creation_stats.filter_cutoff << std::endl << std::endl;
        
        logger.info() << "Total time generating hash table index: " << index_creation_stats.elapsed_hash_index.count() << " s" <<  std::endl;

        std::chrono::duration<double> elapsed = high_resolution_clock::now() - start;
        logger.info() << "Total time indexing: " << elapsed.count() << " s" << std::endl;

        if (!opt.logfile_name.empty()) {
            print_diagnostics(index, opt.logfile_name, index_parameters.k);
            logger.debug() << "Finished printing log stats" << std::endl;
        }
        if (opt.only_gen_index) {
            auto start_write_index = high_resolution_clock::now();
            index.write(opt.index_filename);
            std::chrono::duration<double> elapsed_write_index = high_resolution_clock::now() - start_write_index;
            logger.info() << "Total time writing index: " << elapsed_write_index.count() << " s\n" << std::endl;
            return EXIT_SUCCESS;
        }
    }

    // Map/align reads
        
    auto start_aln_part = high_resolution_clock::now();
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
        out << sam_header(references);
    }

    std::unordered_map<std::thread::id, AlignmentStatistics> log_stats_vec(opt.n_threads);

    logger.info() << "Running in " << (opt.is_SE ? "single-end" : "paired-end") << " mode" << std::endl;

    if (opt.is_SE) {
        auto ks = open_fastq(opt.reads_filename1);

        int input_chunk_size = 100000;
        InputBuffer input_buffer(ks, ks, input_chunk_size);
        OutputBuffer output_buffer(out);

        std::vector<std::thread> workers;
        for (int i = 0; i < opt.n_threads; ++i) {
            std::thread consumer(perform_task_SE, std::ref(input_buffer), std::ref(output_buffer),
                std::ref(log_stats_vec), std::ref(aln_params),
                std::ref(map_param), std::ref(index_parameters), std::ref(references),
                std::ref(index));
            workers.push_back(std::move(consumer));
        }

        for (size_t i = 0; i < workers.size(); ++i) {
            workers[i].join();
        }
    }
    else {
        auto ks1 = open_fastq(opt.reads_filename1);
        auto ks2 = open_fastq(opt.reads_filename2);
        std::unordered_map<std::thread::id, i_dist_est> isize_est_vec(opt.n_threads);

        int input_chunk_size = 100000;
        InputBuffer input_buffer(ks1, ks2, input_chunk_size);
        OutputBuffer output_buffer(out);

        std::vector<std::thread> workers;
        for (int i = 0; i < opt.n_threads; ++i) {
            std::thread consumer(perform_task_PE, std::ref(input_buffer), std::ref(output_buffer),
                std::ref(log_stats_vec), std::ref(isize_est_vec), std::ref(aln_params),
                std::ref(map_param), std::ref(index_parameters), std::ref(references),
                std::ref(index));
            workers.push_back(std::move(consumer));
        }

        for (size_t i = 0; i < workers.size(); ++i) {
            workers[i].join();
        }
    }
    logger.info() << "Done!\n";

    AlignmentStatistics tot_statistics;
    for (auto& it : log_stats_vec) {
        tot_statistics += it.second;
    }
    // Record mapping end time
    std::chrono::duration<double> tot_aln_part = high_resolution_clock::now() - start_aln_part;

    logger.info() << "Total mapping sites tried: " << tot_statistics.tot_all_tried << std::endl
        << "Total calls to ssw: " << tot_statistics.tot_ksw_aligned << std::endl
        << "Calls to ksw (rescue mode): " << tot_statistics.tot_rescued << std::endl
        << "Did not fit strobe start site: " << tot_statistics.did_not_fit << std::endl
        << "Tried rescue: " << tot_statistics.tried_rescue << std::endl
        << "Total time mapping: " << tot_aln_part.count() << " s." << std::endl
        << "Total time reading read-file(s): " << tot_statistics.tot_read_file.count() / opt.n_threads << " s." << std::endl
        << "Total time creating strobemers: " << tot_statistics.tot_construct_strobemers.count() / opt.n_threads << " s." << std::endl
        << "Total time finding NAMs (non-rescue mode): " << tot_statistics.tot_find_nams.count() / opt.n_threads << " s." << std::endl
        << "Total time finding NAMs (rescue mode): " << tot_statistics.tot_time_rescue.count() / opt.n_threads << " s." << std::endl;
    //<< "Total time finding NAMs ALTERNATIVE (candidate sites): " << tot_find_nams_alt.count()/opt.n_threads  << " s." <<  std::endl;
    logger.info() << "Total time sorting NAMs (candidate sites): " << tot_statistics.tot_sort_nams.count() / opt.n_threads << " s." << std::endl
        << "Total time reverse compl seq: " << tot_statistics.tot_rc.count() / opt.n_threads << " s." << std::endl
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
    } catch (const InvalidFile& e) {
        logger.error() << "strobealign: " << e.what() << std::endl;
    }
    return EXIT_FAILURE;
}
