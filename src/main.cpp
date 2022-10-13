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

#include <zlib.h>
#include "args.hxx"
#include "kseq++.hpp"
#include "robin_hood.h"
#include "ssw_cpp.h"
#include "refs.hpp"
#include "exceptions.hpp"
#include "cmdline.hpp"
#include "index.hpp"
#include "pc.hpp"
#include "aln.hpp"
#include "logger.hpp"
#include "version.hpp"

using namespace klibpp;
using std::chrono::high_resolution_clock;


static Logger& logger = Logger::get();


static void print_diagnostics(mers_vector &ref_mers, kmer_lookup &mers_index, std::string logfile_name, int k, int m) {
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
    for (auto &it : mers_index) {
        auto ref_mer = it.second;

        auto offset = std::get<0>(ref_mer);
        auto count = std::get<1>(ref_mer);


        for (size_t j = offset; j < offset + count; ++j) {
            auto r = ref_mers[j];
            auto p = std::get<1>(r);
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

        if ( (count == 1) & (seed_length < max_size) ) {
            log_unique[seed_length] ++;
        }
        if ( (count >= 10) & (seed_length < max_size) ) {
            log_repetitive[seed_length] ++;
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

int est_read_length( klibpp::KStream<gzFile_s*, int (*)(gzFile_s*, void*, unsigned int), klibpp::mode::In_> &ks, int n_reads){
    auto records = ks.read(n_reads);
    int tot_read_len = 0;
    for (size_t i = 0; i < records.size(); ++i) {
        auto record1 = records[i];
        tot_read_len += record1.seq.length();
    }
    int avg_read_len = tot_read_len/n_reads;
    return avg_read_len;
}

/*
 * Return average read length of single-end or paired-end reads.
 * Set filename2 to the empty string if data is single end.
 */
int estimate_read_length(std::string filename1, std::string filename2) {
    bool is_paired = filename2 != "";

    gzFile fp1_tmp = gzopen(filename1.c_str(), "r");
    auto ks1_tmp = make_ikstream(fp1_tmp, gzread);
    auto r1_tmp = est_read_length(ks1_tmp, 500);
    gzclose(fp1_tmp);

    if (is_paired) {
        gzFile fp2_tmp = gzopen(filename2.c_str(), "r");
        auto ks2_tmp = make_ikstream(fp2_tmp, gzread);
        auto r2_tmp = est_read_length(ks2_tmp, 500);
        gzclose(fp2_tmp);
        return (r1_tmp + r2_tmp) / 2;
    } else {
        return r1_tmp;
    }
}

void adjust_mapping_params_depending_on_read_length(mapping_params &map_param, const CommandLineOptions &opt) {
    struct settings {
        int r_threshold;
        int k;
        int s_offset;
        int l;
        int u;
    };
    std::vector<settings> d = {
        settings {75, 20, -4, -4, 2},
        settings {125, 20, -4, -2, 2},
        settings {175, 20, -4, 1, 7},
        settings {275, 20, -4, 4, 13},
        settings {375, 22, -4, 2, 12},
        settings {std::numeric_limits<int>::max(), 23, -6, 2, 12},
    };
    for (const auto& v : d) {
        if (map_param.r <= v.r_threshold) {
            if (!opt.k_set) {
                map_param.k = v.k;
            }
            if (!opt.s_set) {
                map_param.s = map_param.k + v.s_offset;
            }
            map_param.l = v.l;
            map_param.u = v.u;
            break;
        }
    }
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

int main (int argc, char **argv)
{
    CommandLineOptions opt;
    mapping_params map_param;
    std::tie(opt, map_param) = parse_command_line_arguments(argc, argv);

    logger.set_level(opt.verbose ? LOG_DEBUG : LOG_INFO);
    if (!opt.r_set) {
        map_param.r = estimate_read_length(opt.reads_filename1, opt.reads_filename2);
    }
    adjust_mapping_params_depending_on_read_length(map_param, opt);

    if (!opt.max_seed_len_set){
        map_param.max_dist = map_param.r - 70 > map_param.k ? map_param.r - 70 : map_param.k;
        if (map_param.max_dist > 255){
            map_param.max_dist = 255;
        }
    } else {
        map_param.max_dist = opt.max_seed_len - map_param.k; //convert to distance in start positions
    }


    if ( (map_param.c < 64) && (map_param.c > 0)){
        map_param.q = pow (2, map_param.c) - 1;
    } else{
        logger.warning() << "Warning wrong value for parameter c, setting c=8" << std::endl;
        map_param.q = pow (2, 8) - 1;
    }

    map_param.w_min = std::max(1, map_param.k/(map_param.k-map_param.s+1) + map_param.l);
    map_param.w_max = map_param.k/(map_param.k-map_param.s+1) + map_param.u;
    map_param.t_syncmer = (map_param.k-map_param.s)/2 + 1;

    alignment_params aln_params;
    aln_params.match = opt.A;
    aln_params.mismatch = opt.B;
    aln_params.gap_open = opt.O;
    aln_params.gap_extend = opt.E;

    logger.debug() << "Using" << std::endl
        << "k: " << map_param.k << std::endl
        << "s: " << map_param.s << std::endl
        << "w_min: " << map_param.w_min << std::endl
        << "w_max: " << map_param.w_max << std::endl
        << "Read length (r): " << map_param.r << std::endl
        << "Maximum seed length: " << map_param.max_dist + map_param.k << std::endl
        << "Threads: " << opt.n_threads << std::endl
        << "R: " << map_param.R << std::endl
        << "Expected [w_min, w_max] in #syncmers: [" << map_param.w_min << ", " << map_param.w_max << "]" << std::endl
        << "Expected [w_min, w_max] in #nucleotides: [" << (map_param.k-map_param.s+1)*map_param.w_min << ", " << (map_param.k-map_param.s+1)*map_param.w_max << "]" << std::endl
        << "A: " << opt.A << std::endl
        << "B: " << opt.B << std::endl
        << "O: " << opt.O << std::endl
        << "E: " << opt.E << std::endl;

    try {
        map_param.verify();
    }
    catch (BadMappingParameter& e) {
        logger.error() << "A mapping parameter is invalid: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
//    assert(k <= (w/2)*w_min && "k should be smaller than (w/2)*w_min to avoid creating short strobemers");

    //////////// CREATE INDEX OF REF SEQUENCES /////////////////

    StrobemerIndex index;
    References references;
    if (opt.ref_filename.substr(opt.ref_filename.length() - 4) != ".sti") { //assume it is a fasta file if not named ".sti"
        //Generate index from FASTA
    
        // Record index creation start time
        auto start = high_resolution_clock::now();
        auto start_read_refs = start;
        try {
            references = References::from_fasta(opt.ref_filename);
        } catch (const InvalidFasta& e) {
            logger.error() << "strobealign: " << e.what() << std::endl;
            return EXIT_FAILURE;
        }

        std::chrono::duration<double> elapsed_read_refs = high_resolution_clock::now() - start_read_refs;
        logger.info() << "Time reading references: " << elapsed_read_refs.count() << " s" << std::endl;

        if (references.total_length() == 0) {
            logger.error() << "No reference sequences found, aborting" << std::endl;
            return EXIT_FAILURE;
        }

        index.populate(references, map_param);
        index.filter_cutoff = map_param.filter_cutoff;

        // Record index creation end time
        std::chrono::duration<double> elapsed = high_resolution_clock::now() - start;
        logger.info() << "Total time indexing: " << elapsed.count() << " s" << std::endl;

        if (opt.index_log) {
            print_diagnostics(index.flat_vector, index.mers_index, opt.logfile_name, map_param.k, map_param.max_dist + map_param.k);
            logger.debug() << "Finished printing log stats" << std::endl;
        }
        
        // If the program was called with the -i flag, write the index to file
        if (opt.only_gen_index) { // If the program was called with the -i flag, we do not do the alignment
            auto start_write_index = high_resolution_clock::now();
            index.write(references, opt.index_out_filename);
            std::chrono::duration<double> elapsed_write_index = high_resolution_clock::now() - start_write_index;
            logger.info() << "Total time writing index: " << elapsed_write_index.count() << " s\n" << std::endl;
        }
    }
    else {
        //load index from file
        auto start_read_index = high_resolution_clock::now();
        index.read(references, opt.ref_filename);
        std::chrono::duration<double> elapsed_read_index = high_resolution_clock::now() - start_read_index;
        logger.info() << "Total time reading index: " << elapsed_read_index.count() << " s\n" << std::endl;

    }

    ///////////////////////////// MAP ///////////////////////////////////////
    
    if (!(opt.only_gen_index && opt.reads_filename1.empty())) { // If the program was called with the -i flag and the fastqs are not specified, we don't run any alignment
        
        map_param.filter_cutoff = index.filter_cutoff; //This is calculated when building the filter and needs to be filled in

        // Record matching time
        auto start_aln_part = high_resolution_clock::now();

        //    std::ifstream query_file(reads_filename);

        map_param.rescue_cutoff = map_param.R < 100 ? map_param.R * map_param.filter_cutoff : 1000;
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
        //    std::ofstream out;
        //    out.open(opt.output_file_name);

        //    std::stringstream sam_output;
        //    std::stringstream paf_output;

        if (map_param.is_sam_out) {
            out << sam_header(references);
        }

        std::unordered_map<std::thread::id, AlignmentStatistics> log_stats_vec(opt.n_threads);


        logger.info() << "Running in " << (opt.is_SE ? "single-end" : "paired-end") << " mode" << std::endl;

        if (opt.is_SE) {
            gzFile fp = gzopen(opt.reads_filename1.c_str(), "r");
            auto ks = make_ikstream(fp, gzread);

            ////////// ALIGNMENT START //////////

            int input_chunk_size = 100000;
            // Create Buffers
            InputBuffer input_buffer(ks, ks, input_chunk_size);
            OutputBuffer output_buffer(out);

            std::vector<std::thread> workers;
            for (int i = 0; i < opt.n_threads; ++i) {
                std::thread consumer(perform_task_SE, std::ref(input_buffer), std::ref(output_buffer),
                    std::ref(log_stats_vec), std::ref(aln_params),
                    std::ref(map_param), std::ref(references),
                    std::ref(index.mers_index), std::ref(index.flat_vector));
                workers.push_back(std::move(consumer));
            }

            for (size_t i = 0; i < workers.size(); ++i) {
                workers[i].join();
            }

            gzclose(fp);
        }
        else {
            gzFile fp1 = gzopen(opt.reads_filename1.c_str(), "r");
            auto ks1 = make_ikstream(fp1, gzread);
            gzFile fp2 = gzopen(opt.reads_filename2.c_str(), "r");
            auto ks2 = make_ikstream(fp2, gzread);
            std::unordered_map<std::thread::id, i_dist_est> isize_est_vec(opt.n_threads);

            ////////// ALIGNMENT START //////////
            /////////////////////////////////////

            int input_chunk_size = 100000;
            // Create Buffers
            InputBuffer input_buffer(ks1, ks2, input_chunk_size);
            OutputBuffer output_buffer(out);

            std::vector<std::thread> workers;
            for (int i = 0; i < opt.n_threads; ++i) {
                std::thread consumer(perform_task_PE, std::ref(input_buffer), std::ref(output_buffer),
                    std::ref(log_stats_vec), std::ref(isize_est_vec), std::ref(aln_params),
                    std::ref(map_param), std::ref(references),
                    std::ref(index.mers_index), std::ref(index.flat_vector));
                workers.push_back(std::move(consumer));
            }

            for (size_t i = 0; i < workers.size(); ++i) {
                workers[i].join();
            }

            /////////////////////////////////////
            /////////////////////////////////////

            gzclose(fp1);
            gzclose(fp2);
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

        /////////////////////// FIND AND OUTPUT NAMs ///////////////////////////////
    }
}
