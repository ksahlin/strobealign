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
#include "fasta.hpp"
#include "exceptions.hpp"
#include "index.hpp"
#include "pc.hpp"
#include "aln.hpp"
#include "version.hpp"

using namespace klibpp;
using std::chrono::high_resolution_clock;

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
        auto hash_refmer = it.first;
        auto ref_mer = it.second;

        auto offset = std::get<0>(ref_mer);
        auto count = std::get<1>(ref_mer);


        for (size_t j = offset; j < offset + count; ++j) {
            auto r = ref_mers[j];
            auto p = std::get<1>(r);
            int bit_alloc = 8;
            int r_id = (p >> bit_alloc);
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
               std::cerr << "Detected seed size over " << max_size << " bp (can happen, e.g., over centromere): " << seed_length << std::endl;
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

    for (int i=0 ; i < log_count.size(); ++i) {
        if (log_count[i] > 0) {
            double e_count = log_count_squared[i] / log_count[i];
            log_file << i << ',' << log_count[i] << ',' << e_count << std::endl;
        }
    }

    // Get median
    int n = 0;
    int median = 0;
    for (int i=0 ; i < log_count.size(); ++i) {
        n += log_count[i];
        if ( n >= tot_seed_count/2){
            median = i;
            break;
        }
    }
    // Get median 1000 limit
    int n_lim = 0;
    int median_lim = 0;
    for (int i=0 ; i < log_count_1000_limit.size(); ++i) {
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


struct CommandLineOptions {
    int A { 2 };
    int B { 8 };
    int O { 12 };
    int E { 1 };
    int max_seed_len;
    std::string output_file_name;
    std::string logfile_name { "log.csv" };
    int n_threads { 3 };
    std::string ref_filename; //This is either a fasta file or an index file - if fasta, indexing will be run
    std::string reads_filename1;
    std::string reads_filename2;
    bool is_SE { true };
    bool k_set { false };
    bool write_to_stdout { true };
    bool index_log { false };
    bool r_set { false };
    bool max_seed_len_set { false };
    bool s_set{ false };
    bool only_gen_index{ false };
    std::string index_out_filename;
};


std::pair<CommandLineOptions, mapping_params> parse_command_line_arguments(int argc, char **argv) {

    args::ArgumentParser parser("StrobeAlign " VERSION_STRING);
    parser.helpParams.showTerminator = false;
    parser.helpParams.helpindent = 20;
    parser.helpParams.width = 90;
    parser.helpParams.programName = "strobealign";
    parser.helpParams.shortSeparator = " ";

    args::HelpFlag help(parser, "help", "Print help and exit", {'h', "help"});
    args::ValueFlag<int> threads(parser, "threads", "Number of threads [3]", {'t', "threads"});

    args::Group io(parser, "Input/output:");
    args::ValueFlag<std::string> o(parser, "STR", "redirect output to file [stdout]", {'o'});
    args::Flag x(parser, "x", "Only map reads, no base level alignment (produces paf file)", {'x'});
    args::ValueFlag<int> N(parser, "INT", "retain at most INT secondary alignments (is upper bounded by -M, and depends on -S) [0]", {'N'});
    args::ValueFlag<std::string> L(parser, "STR", "Print statistics of indexing to logfile [log.csv]", {'L'});
    args::ValueFlag<std::string> i(parser, "index", "Generates an index (.sti) file", { 'i' });

    args::Group seeding(parser, "Seeding:");
    //args::ValueFlag<int> n(parser, "INT", "number of strobes [2]", {'n'});
    args::ValueFlag<int> r(parser, "INT", "Mean read length. This parameter is estimated from first 500 records in each read file. No need to set this explicitly unless you have a reason. [disabled]", {'r'});
    args::ValueFlag<int> m(parser, "INT", "Maximum seed length. Defaults to r - 50. For reasonable values on -l and -u, the seed length distribution is usually determined by parameters l and u. Then, this parameter is only active in regions where syncmers are very sparse.", {'m'});

    args::ValueFlag<int> k(parser, "INT", "strobe length, has to be below 32. [20]", {'k'});
    args::ValueFlag<int> l(parser, "INT", "Lower syncmer offset from k/(k-s+1). Start sample second syncmer k/(k-s+1) + l syncmers downstream [0]", {'l'});

    args::ValueFlag<int> u(parser, "INT", "Upper syncmer offset from k/(k-s+1). End sample second syncmer k/(k-s+1) + u syncmers downstream [7]", {'u'});
    args::ValueFlag<int> c(parser, "INT", "Bitcount length between 2 and 63. [8]", {'c'});
    args::ValueFlag<int> s(parser, "INT", "Submer size used for creating syncmers [k-4]. Only even numbers on k-s allowed. A value of s=k-4 roughly represents w=10 as minimizer window [k-4]. It is recommended not to change this parameter unless you have a good understanding of syncmers as it will drastically change the memory usage and results with non default values.", {'s'});

    args::Group alignment(parser, "Alignment:");
    args::ValueFlag<int> A(parser, "INT", "matching score [2]", {'A'});
    args::ValueFlag<int> B(parser, "INT", "mismatch penalty [8]", {'B'});
    args::ValueFlag<int> O(parser, "INT", "gap open penalty [12]", {'O'});
    args::ValueFlag<int> E(parser, "INT", "gap extension penalty [1]", {'E'});

    args::Group search(parser, "Search parameters:");
    args::ValueFlag<float> f(parser, "FLOAT", "top fraction of repetitive strobemers to filter out from sampling [0.0002]", {'f'});
    args::ValueFlag<float> S(parser, "FLOAT", "Try candidate sites with mapping score at least S of maximum mapping score [0.5]", {'S'});
    args::ValueFlag<int> M(parser, "INT", "Maximum number of mapping sites to try [20]", {'M'});
    args::ValueFlag<int> R(parser, "INT", "Rescue level. Perform additional search for reads with many repetitive seeds filtered out. This search includes seeds of R*repetitive_seed_size_filter (default: R=2). Higher R than default makes StrobeAlign significantly slower but more accurate. R <= 1 deactivates rescue and is the fastest.", {'R'});

    // <ref.fa> <reads1.fast[a/q.gz]> [reads2.fast[a/q.gz]]
    args::Positional<std::string> ref_filename(parser, "reference", "A pregenerated strobemers index file (.sti) or a reference in FASTA format", args::Options::Required);
    args::Positional<std::string> reads1_filename(parser, "reads1", "Reads 1 in FASTA or FASTQ format, optionally gzip compressed");
    args::Positional<std::string> reads2_filename(parser, "reads2", "Reads 2 in FASTA or FASTQ format, optionally gzip compressed");


    try {
        parser.ParseCLI(argc, argv);
    }
    catch (const args::Completion& e) {
        std::cout << e.what();
        exit(0);
    }
    catch (const args::Help&)
    {
        std::cout << parser;
        exit(0);
    }
    catch (const args::Error& e)
    {
        std::cerr << parser;
        std::cerr << "Error: " << e.what() << std::endl;
        exit(1);
    }

    CommandLineOptions opt;

    mapping_params map_param;
    map_param.max_secondary = 0;
    map_param.n = 2;
    map_param.k = 20;
    map_param.s = map_param.k - 4;
    map_param.f = 0.0002;
    map_param.R = 2;
    map_param.dropoff_threshold = 0.5;
    map_param.maxTries = 20;
    map_param.l = 0;
    map_param.u = 7;
    map_param.c = 8;
    map_param.r = 150;
    map_param.max_dist = map_param.r - 50 < 255 ? map_param.r-50 : 255;
    map_param.is_sam_out = true;  // true: align, false: map

    if (threads) { opt.n_threads = args::get(threads); }

    // Input/output
    if (o) { opt.output_file_name = args::get(o); opt.write_to_stdout = false; }
    if (x) { map_param.is_sam_out = false; }
    if (N) { map_param.max_secondary = args::get(N); }
    if (L) { opt.logfile_name = args::get(L); opt.index_log = true; }
    if (i) { opt.only_gen_index = true; opt.index_out_filename = args::get(i); }

    // Seeding
    if (r) { map_param.r = args::get(r); opt.r_set = true; }
    if (m) { opt.max_seed_len = args::get(m); opt.max_seed_len_set = true; }
    if (k) { map_param.k = args::get(k); opt.k_set = true; }
    if (l) { map_param.l = args::get(l); }
    if (u) { map_param.u = args::get(u); }
    if (c) { map_param.c = args::get(c); }
    if (s) { map_param.s = args::get(s); opt.s_set = true; }

    // Alignment
    // if (n) { n = args::get(n); }
    if (A) { opt.A = args::get(A); }
    if (B) { opt.B = args::get(B); }
    if (O) { opt.O = args::get(O); }
    if (E) { opt.E = args::get(E); }

    // Search parameters
    if (f) { map_param.f = args::get(f); }
    if (S) { map_param.dropoff_threshold = args::get(S); }
    if (M) { map_param.maxTries = args::get(M); }
    if (R) { map_param.R = args::get(R); }

    // Reference and read files
    opt.ref_filename = args::get(ref_filename);
    opt.reads_filename1 = args::get(reads1_filename);

    if (reads2_filename) {
        opt.reads_filename2 = args::get(reads2_filename);
        opt.is_SE = false;
    } else {
        opt.reads_filename2 = std::string();
        opt.is_SE = true;
    }

    //If not generating index, fastq1 is mandatory:
    if (opt.reads_filename1.empty() && !opt.only_gen_index) {
        std::cerr << "Reads file (fastq) must be specified." << std::endl;
//		exit(1);
    }

    return std::make_pair(opt, map_param);
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

std::pair<mers_vector, kmer_lookup>  create_index(mapping_params& map_param, std::vector<std::string>& ref_seqs, uint64_t total_ref_seq_size) {
    auto start_flat_vector = high_resolution_clock::now();

    mers_vector flat_vector;

    ind_mers_vector ind_flat_vector; //includes hash - for sorting, will be discarded later
    int approx_vec_size = total_ref_seq_size / (map_param.k-map_param.s+1);
    std::cerr << "ref vector approximate size: " << approx_vec_size << std::endl;
    ind_flat_vector.reserve(approx_vec_size);
    for(size_t i = 0; i < ref_seqs.size(); ++i)
    {
//        std::cerr << i << " " << i_mod << std::endl;
        seq_to_randstrobes2(ind_flat_vector, map_param.n, map_param.k, map_param.w_min, map_param.w_max, ref_seqs[i], i, map_param.s, map_param.t_syncmer, map_param.q, map_param.max_dist);
    }
    std::cerr << "Ref vector actual size: " << ind_flat_vector.size() << std::endl;
    //ind_flat_vector.shrink_to_fit(); //I think this costs performance and is no longer needed, it will soon be deallocated anyway

    std::chrono::duration<double> elapsed_generating_seeds = high_resolution_clock::now() - start_flat_vector;
    std::cerr << "Time generating seeds: " << elapsed_generating_seeds.count() << " s\n" <<  std::endl;

//    create vector of vectors here nr_threads
//    std::vector<std::vector<std::tuple<uint64_t, unsigned int, unsigned int, unsigned int>>> vector_per_ref_chr(opt.n_threads);
//    for(size_t i = 0; i < ref_seqs.size(); ++i)
//    {
//        mers_vector randstrobes2; // pos, chr_id, kmer hash value
//        std::cerr << "Started thread: " << omp_get_thread_num() << " chr size: " << ref_lengths[i] << " acc map:" << acc_map[i] << std::endl;
//        randstrobes2 = seq_to_randstrobes2(n, k, w_min, w_max, ref_seqs[i], i, s, t);
//        for (auto &t : randstrobes2)
//        {
//            vector_per_ref_chr[omp_get_thread_num()].push_back(t);
//        }
//        std::cerr << "Completed thread: " << omp_get_thread_num() << " chr size: " << ref_lengths[i] << " acc map:" << acc_map[i] << std::endl;
//    }

    auto start_sorting = high_resolution_clock::now();
//    std::cerr << "Reserving flat vector size: " << approx_vec_size << std::endl;
//    all_mers_vector_tmp.reserve(approx_vec_size); // reserve size corresponding to sum of lengths of all sequences divided by expected sampling
    std::sort(ind_flat_vector.begin(), ind_flat_vector.end());
    std::chrono::duration<double> elapsed_sorting_seeds = high_resolution_clock::now() - start_sorting;
    std::cerr << "Time sorting seeds: " << elapsed_sorting_seeds.count() << " s\n" <<  std::endl;

    //Split up the sorted vector into a vector with the hash codes and the flat vector to keep in the index.
    //The hash codes are only needed when generating the index and can be discarded afterwards.
    //We want to do this split-up before creating the hash table to avoid a memory peak - the flat_vector is 
    //smaller - doubling that size temporarily will not cause us to go above peak memory.
    auto start_copy_flat_vector = high_resolution_clock::now();
    hash_vector h_vector;
    flat_vector.reserve(ind_flat_vector.size());
    h_vector.reserve(ind_flat_vector.size());
    for (std::size_t i = 0; i < ind_flat_vector.size(); ++i) {
        flat_vector.push_back(std::make_tuple(std::get<1>(ind_flat_vector[i]), std::get<2>(ind_flat_vector[i])));
        h_vector.push_back(std::get<0>(ind_flat_vector[i]));
    }
    uint64_t unique_mers = count_unique_elements(h_vector);
    std::cerr << "Unique strobemers: " << unique_mers << std::endl;
    ind_mers_vector().swap(ind_flat_vector);   //deallocate the vector used for sorting - ind_flat_vector.clear() will not deallocate, 
                                               //swap with an empty vector will. I tested (in debug mode), the buffer gets deleted

    std::chrono::duration<double> elapsed_copy_flat_vector = high_resolution_clock::now() - start_copy_flat_vector;
    std::cerr << "Time copying flat vector: " << elapsed_copy_flat_vector.count() << " s\n" << std::endl;


    std::chrono::duration<double> elapsed_flat_vector = high_resolution_clock::now() - start_flat_vector;
    std::cerr << "Total time generating flat vector: " << elapsed_flat_vector.count() << " s\n" <<  std::endl;

    auto start_hash_index = high_resolution_clock::now();
    kmer_lookup mers_index; // k-mer -> (offset in flat_vector, occurence count )
    mers_index.reserve(unique_mers);
    // construct index over flat array
    map_param.filter_cutoff = index_vector(h_vector, mers_index, map_param.f);
    std::chrono::duration<double> elapsed_hash_index = high_resolution_clock::now() - start_hash_index;
    std::cerr << "Total time generating hash table index: " << elapsed_hash_index.count() << " s\n" <<  std::endl;

//    mers_vector_reduced all_mers_vector;
//    all_mers_vector = remove_kmer_hash_from_flat_vector(flat_vector);
    /* destroy vector */
//    flat_vector.clear();

    return make_pair(flat_vector, mers_index);
}

/*
 * Return formatted SAM header as a string
 */
std::string sam_header(idx_to_acc& acc_map, const std::vector<unsigned int>& ref_lengths) {
    std::stringstream out;
    int nr_ref = acc_map.size();
    for (int i = 0; i < nr_ref; ++i) {
//        for (auto &it : acc_map) {
//            out << "@SQ\tSN:" << it.second << "\tLN:" << ref_lengths[it.first] << "\n";
        out << "@SQ\tSN:" << acc_map[i] << "\tLN:" << ref_lengths[i] << "\n";
    }
    out << "@PG\tID:strobealign\tPN:strobealign\tVN:" VERSION_STRING "\tCL:strobealign\n";
    return out.str();
}

int main (int argc, char **argv)
{
    CommandLineOptions opt;
    mapping_params map_param;
    std::tie(opt, map_param) = parse_command_line_arguments(argc, argv);

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
        std::cerr << "Warning wrong value for parameter c, setting c=8" << std::endl;
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

    std::cerr << "Using" << std::endl
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
	std::cerr << "A mapping parameter is invalid: " << e.what() << std::endl;
	return EXIT_FAILURE;
    }
//    assert(k <= (w/2)*w_min && "k should be smaller than (w/2)*w_min to avoid creating short strobemers");

    //////////// CREATE INDEX OF REF SEQUENCES /////////////////

    st_index index;

    if (opt.ref_filename.substr(opt.ref_filename.length() - 4) != ".sti") { //assume it is a fasta file if not named ".sti"
        //Generate index from FASTA
    
        // Record index creation start time
        auto start = high_resolution_clock::now();
        auto start_read_refs = start;
        uint64_t total_ref_seq_size;
        try {
            total_ref_seq_size = read_references(index.ref_seqs, index.ref_lengths, index.acc_map, opt.ref_filename);
        } catch (const InvalidFasta& e) {
            std::cerr << "strobealign: " << e.what() << std::endl;
            return EXIT_FAILURE;
        }
        std::chrono::duration<double> elapsed_read_refs = high_resolution_clock::now() - start_read_refs;
        std::cerr << "Time reading references: " << elapsed_read_refs.count() << " s\n" << std::endl;

        if (total_ref_seq_size == 0) {
            std::cerr << "No reference sequences found, aborting.." << std::endl;
            return EXIT_FAILURE;
        }

        std::tie(index.flat_vector, index.mers_index) = create_index(map_param, index.ref_seqs, total_ref_seq_size);
        index.filter_cutoff = map_param.filter_cutoff;

        // Record index creation end time
        std::chrono::duration<double> elapsed = high_resolution_clock::now() - start;
        std::cerr << "Total time indexing: " << elapsed.count() << " s\n" << std::endl;

        if (opt.index_log) {
            std::cerr << "Printing log stats" << std::endl;
            print_diagnostics(index.flat_vector, index.mers_index, opt.logfile_name, map_param.k, map_param.max_dist + map_param.k);
            std::cerr << "Finished printing log stats" << std::endl;
        }
        
        // If the program was called with the -i flag, write the index to file
        if (opt.only_gen_index) { // If the program was called with the -i flag, we do not do the alignment
            auto start_write_index = high_resolution_clock::now();
            write_index(index, opt.index_out_filename);
            std::chrono::duration<double> elapsed_write_index = high_resolution_clock::now() - start_write_index;
            std::cerr << "Total time writing index: " << elapsed_write_index.count() << " s\n" << std::endl;
        }
    }
    else {
        //load index from file
        auto start_read_index = high_resolution_clock::now();
        read_index(index, opt.ref_filename);
        std::chrono::duration<double> elapsed_read_index = high_resolution_clock::now() - start_read_index;
        std::cerr << "Total time reading index: " << elapsed_read_index.count() << " s\n" << std::endl;

    }

    ///////////////////////////// MAP ///////////////////////////////////////
    
    if (!(opt.only_gen_index && opt.reads_filename1.empty())) { // If the program was called with the -i flag and the fastqs are not specified, we don't run any alignment
        
        map_param.filter_cutoff = index.filter_cutoff; //This is calculated when building the filter and needs to be filled in

        // Record matching time
        auto start_aln_part = high_resolution_clock::now();

        //    std::ifstream query_file(reads_filename);

        map_param.rescue_cutoff = map_param.R < 100 ? map_param.R * map_param.filter_cutoff : 1000;
        std::cerr << "Using rescue cutoff: " << map_param.rescue_cutoff << std::endl;

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
            out << sam_header(index.acc_map, index.ref_lengths);
        }

        std::unordered_map<std::thread::id, logging_variables> log_stats_vec(opt.n_threads);


        std::cerr << "Running in " << (opt.is_SE ? "single-end" : "paired-end") << " mode" << std::endl;

        if (opt.is_SE) {
            gzFile fp = gzopen(opt.reads_filename1.c_str(), "r");
            auto ks = make_ikstream(fp, gzread);

            ////////// ALIGNMENT START //////////

            int input_chunk_size = 100000;
            // Create Buffers
            InputBuffer input_buffer = { {}, {}, {}, {}, {}, ks, ks, false, 0, input_chunk_size };
            OutputBuffer output_buffer = { {}, {}, {}, 0, out };

            std::vector<std::thread> workers;
            for (int i = 0; i < opt.n_threads; ++i) {
                std::thread consumer(perform_task_SE, std::ref(input_buffer), std::ref(output_buffer),
                    std::ref(log_stats_vec), std::ref(aln_params),
                    std::ref(map_param), std::ref(index.ref_lengths), std::ref(index.ref_seqs),
                    std::ref(index.mers_index), std::ref(index.flat_vector), std::ref(index.acc_map));
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
            InputBuffer input_buffer = { {}, {}, {}, {}, {}, ks1, ks2, false, 0, input_chunk_size };
            OutputBuffer output_buffer = { {}, {}, {}, 0, out };

            std::vector<std::thread> workers;
            for (int i = 0; i < opt.n_threads; ++i) {
                std::thread consumer(perform_task_PE, std::ref(input_buffer), std::ref(output_buffer),
                    std::ref(log_stats_vec), std::ref(isize_est_vec), std::ref(aln_params),
                    std::ref(map_param), std::ref(index.ref_lengths), std::ref(index.ref_seqs),
                    std::ref(index.mers_index), std::ref(index.flat_vector), std::ref(index.acc_map));
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
        std::cerr << "Done!\n";

        logging_variables tot_log_vars;
        for (auto& it : log_stats_vec) {
            tot_log_vars += it.second;
        }
        // Record mapping end time
        std::chrono::duration<double> tot_aln_part = high_resolution_clock::now() - start_aln_part;

        std::cerr << "Total mapping sites tried: " << tot_log_vars.tot_all_tried << std::endl
            << "Total calls to ssw: " << tot_log_vars.tot_ksw_aligned << std::endl
            << "Calls to ksw (rescue mode): " << tot_log_vars.tot_rescued << std::endl
            << "Did not fit strobe start site: " << tot_log_vars.did_not_fit << std::endl
            << "Tried rescue: " << tot_log_vars.tried_rescue << std::endl
            << "Total time mapping: " << tot_aln_part.count() << " s." << std::endl
            << "Total time reading read-file(s): " << tot_log_vars.tot_read_file.count() / opt.n_threads << " s." << std::endl
            << "Total time creating strobemers: " << tot_log_vars.tot_construct_strobemers.count() / opt.n_threads << " s." << std::endl
            << "Total time finding NAMs (non-rescue mode): " << tot_log_vars.tot_find_nams.count() / opt.n_threads << " s." << std::endl
            << "Total time finding NAMs (rescue mode): " << tot_log_vars.tot_time_rescue.count() / opt.n_threads << " s." << std::endl;
        //<< "Total time finding NAMs ALTERNATIVE (candidate sites): " << tot_find_nams_alt.count()/opt.n_threads  << " s." <<  std::endl;
        std::cerr << "Total time sorting NAMs (candidate sites): " << tot_log_vars.tot_sort_nams.count() / opt.n_threads << " s." << std::endl
            << "Total time reverse compl seq: " << tot_log_vars.tot_rc.count() / opt.n_threads << " s." << std::endl
            << "Total time base level alignment (ssw): " << tot_log_vars.tot_extend.count() / opt.n_threads << " s." << std::endl
            << "Total time writing alignment to files: " << tot_log_vars.tot_write_file.count() << " s." << std::endl;

        /////////////////////// FIND AND OUTPUT NAMs ///////////////////////////////
    }
}
