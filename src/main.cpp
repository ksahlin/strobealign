#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>
#include <chrono>  // for high_resolution_clock
//#include <omp.h>
#include <zlib.h>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <inttypes.h>

#include "kseq++.hpp"
using namespace klibpp;
#include "robin_hood.h"
#include "index.hpp"
//#include "ksw2.h"
#include "ssw_cpp.h"
#include "pc.hpp"
#include "aln.hpp"

//develop
#include <chrono>
#include <thread>
#include <sstream>







static uint64_t read_references(std::vector<std::string> &seqs, std::vector<unsigned int> &lengths, idx_to_acc &acc_map, std::string fn)
{
    uint64_t total_ref_seq_size = 0;
    std::ifstream file(fn);
    std::string line, seq;
    int ref_index = 0;

    if (!file.good()) {
        std::cerr << "Unable to read from file " << fn << std::endl;
        return total_ref_seq_size;
    }

    if (file.peek() != '>') {
        std::cerr << fn << " is not a valid FASTA file." << std::endl;
        return total_ref_seq_size;
    }

    while (getline(file, line)) {
        if (line[0] == '>') {
//            std::cerr << ref_index << " " << line << std::endl;
            if (seq.length() > 0){
//                seqs[ref_index -1] = seq;
                seqs.push_back(seq);
                lengths.push_back(seq.length());
                total_ref_seq_size += seq.length();
//                std::cerr << ref_index - 1 << " here " << seq << " " << seq.length() << " " << seq.size() << std::endl;
//                generate_kmers(h, k, seq, ref_index);
            }
//            acc_map[ref_index] = line.substr(1, line.length() -1); //line;
            acc_map[ref_index] = line.substr(1, line.find(' ') -1); // cutting at first space;
            ref_index++;
            seq = "";
        }
        else {
            seq += line;
        }
    }
    if (seq.length() > 0){
//        seqs[ref_index -1] = seq;
        seqs.push_back(seq);
        lengths.push_back(seq.length());
        total_ref_seq_size += seq.length();
//        std::cerr << ref_index -1 << " here2 " << seq << std::endl;
//        generate_kmers(h, k, seq, ref_index);
    }

    return total_ref_seq_size;
}



static inline void print_diagnostics(mers_vector &ref_mers, kmer_lookup &mers_index, std::string logfile_name, int k, int m) {
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
            auto p = std::get<2>(r);
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


void print_usage() {
    std::cerr << "\n";
    std::cerr << "StrobeAlign VERSION 0.7.1 \n";
    std::cerr << "\n";
    std::cerr << "StrobeAlign [options] <ref.fa> <reads1.fast[a/q.gz]> [reads2.fast[a/q.gz]]\n";
    std::cerr << "options:\n";
    std::cerr << "\n";
    std::cerr << "\t-h Print help and exit\n";
    std::cerr << "\n";
    std::cerr << "Resources:\n";
    std::cerr << "\t-t INT number of threads [3]\n";

    std::cerr << "\n";
    std::cerr << "Input/output:\n";
    std::cerr << "\t-o STR redirect output to file [stdout]\n";
    std::cerr << "\t-x Only map reads, no base level alignment (produces paf file)\n";
    std::cerr << "\t-N INT retain at most INT secondary alignments (is upper bounded by -M, and depends on -S) [0]  \n";
    std::cerr << "\t-L STR Print statistics of indexing to logfie [log.csv] \n";

    std::cerr << "\n";
    std::cerr << "Seeding:\n";
//    std::cerr << "\t-n INT number of strobes [2]\n";
    std::cerr << "\t-r INT Mean read length. This parameter is estimated from first 500 records in each read file. No need to set this explicitly  \n\t   unless you have a reason. [disabled] \n";
    std::cerr << "\t-m INT Maximum seed length. Defaults to r - 50. For reasonable values on -l and -u, the seed length distribution is usually determined by\n\t   parameters l and u. Then, this parameter is only active in regions where syncmers are very sparse.\n";
    std::cerr << "\t-k INT strobe length, has to be below 32. [20]\n";
    std::cerr << "\t-l INT Lower syncmer offset from k/(k-s+1). Start sample second syncmer k/(k-s+1) + l syncmers downstream [0]\n";
    std::cerr << "\t-u INT Upper syncmer offset from k/(k-s+1). End sample second syncmer k/(k-s+1) + u syncmers downstream [7]\n";
    std::cerr << "\t-c INT Bitcount length between 2 and 63. [8]\n";
    std::cerr << "\t-s INT Submer size used for creating syncmers [k-4]. Only even numbers on k-s allowed.\n\t   A value of s=k-4 roughly represents w=10 as minimizer window [k-4]. It is recommended not to change this parameter\n\t   unless you have a good understanding of syncmenrs as it will drastically change the memory usage and results with non default values. \n";


    std::cerr << "\n";
    std::cerr << "Alignment:\n";
//    std::cerr << "\t-n INT number of strobes [2]\n";
    std::cerr << "\t-A INT matching score [2]\n";
    std::cerr << "\t-B INT mismatch penalty [8]\n";
    std::cerr << "\t-O INT gap open penalty [12]\n";
    std::cerr << "\t-E INT gap extension penalty [1]\n";


    std::cerr << "\n";
    std::cerr << "Search parameters:\n";
    std::cerr << "\t-f FLOAT top fraction of repetitive strobemers to filter out from sampling [0.0002]\n";
    std::cerr << "\t-S FLOAT Try candidate sites with mapping score at least S of maximum mapping score [0.5]\n";
    std::cerr << "\t-M INT Maximum number of mapping sites to try [20]\n";
    std::cerr << "\t-R INT Rescue level. Perform additional search for reads with many repetitive seeds filtered out.\n\t   This search includes seeds of R*repetitive_seed_size_filter (default: R=2). Higher R than default makes StrobeAlign\n\t   significantly slower but more accurate. R <= 1 deactivates rescue and is the fastest. \n";
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
    std::string ref_filename;
    const char *reads_filename1;
    const char *reads_filename2;
    bool is_SE { true };
    bool k_set { false };
    bool write_to_stdout { true };
    bool index_log { false };
    bool r_set { false };
    bool max_seed_len_set { false };
    bool s_set { false };
};


std::pair<CommandLineOptions, mapping_params> parse_command_line_arguments(int argc, char **argv) {

    // Default parameters
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

    int opn = 1;
    while (opn < argc) {
        bool flag = false;
        if (argv[opn][0] == '-') {
            char ch = argv[opn][1];
            flag = true;
            switch (ch) {
//                case 'n':
//                    n = std::stoi(argv[opn + 1]);
//                    opn += 2;
//                    break;
                case 'h':
                    print_usage();
                    exit(0);
                    break;
                case 't':
                    opt.n_threads = std::stoi(argv[opn + 1]);
                    opn += 2;
                    break;
                case 'k':
                    map_param.k = std::stoi(argv[opn + 1]);
                    opn += 2;
                    opt.k_set = true;
                    break;
                case 'o':
                    opt.output_file_name = argv[opn + 1];
                    opn += 2;
                    opt.write_to_stdout = false;
                    break;
                case 'L':
                    opt.logfile_name = argv[opn + 1];
                    opt.index_log = true;
                    opn += 2;
                    break;
                case 's':
                    map_param.s = std::stoi(argv[opn + 1]);
                    opn += 2;
                    opt.s_set = true;
                    break;
                case 'f':
                    map_param.f = std::stof(argv[opn + 1]);
                    opn += 2;
                    break;
                case 'x':
                    map_param.is_sam_out = false;
                    opn += 1;
                    break;
                case 'R':
                    map_param.R = std::stoi(argv[opn + 1]);
                    opn += 2;
                    break;
                case 'l':
                    map_param.l = std::stoi(argv[opn + 1]);
                    opn += 2;
                    break;
                case 'u':
                    map_param.u = std::stoi(argv[opn + 1]);
                    opn += 2;
                    break;
                case 'c':
                    map_param.c = std::stoi(argv[opn + 1]);
                    opn += 2;
                    break;
                case 'r':
                    map_param.r = std::stoi(argv[opn + 1]);
                    opn += 2;
                    opt.r_set = true;
                    break;
                case 'm':
                    opt.max_seed_len = std::stoi(argv[opn + 1]);
                    opn += 2;
                    opt.max_seed_len_set = true;
                    break;
                case 'S':
                    map_param.dropoff_threshold = std::stof(argv[opn + 1]);
                    opn += 2;
                    break;
                case 'M':
                    map_param.maxTries = std::stoi(argv[opn + 1]);
                    opn += 2;
                    break;
                case 'A':
                    opt.A = std::stoi(argv[opn + 1]);
                    opn += 2;
                    break;
                case 'B':
                    opt.B = std::stoi(argv[opn + 1]);
                    opn += 2;
                    break;
                case 'O':
                    opt.O = std::stoi(argv[opn + 1]);
                    opn += 2;
                    break;
                case 'E':
                    opt.E = std::stoi(argv[opn + 1]);
                    opn += 2;
                    break;
                case 'N':
                    map_param.max_secondary = std::stoi(argv[opn + 1]);
                    opn += 2;
                    break;
                default:
                    print_usage();
                    exit(1);
                    break;
            }
        }
        if (!flag)
            break;
    }

    // File name to reference
    opt.ref_filename = argv[opn];
    opn++;

    if (opn == argc - 1) {
        opt.reads_filename1 = argv[opn];
        opt.reads_filename2 = nullptr;
        opt.is_SE = true;
    } else if (opn == argc - 2) {
        opt.reads_filename1 = argv[opn];
        opn++;
        opt.reads_filename2 = argv[opn];
        opt.is_SE = false;
    } else {
        print_usage();
        exit(1);
    }

    return std::make_pair(opt, map_param);
}

int main (int argc, char **argv)
{
    std::string choice = "randstrobes";
    CommandLineOptions opt;
    mapping_params map_param;
    std::tie(opt, map_param) = parse_command_line_arguments(argc, argv);

    if (opt.reads_filename2 == nullptr) {
        if (!opt.r_set){
            gzFile fp1_tmp = gzopen(opt.reads_filename1, "r");
            auto ks1_tmp = make_ikstream(fp1_tmp, gzread);
            auto r1_tmp = est_read_length(ks1_tmp, 500);
            gzclose(fp1_tmp);
            map_param.r = r1_tmp;
        }
    } else {
        if (!opt.r_set) {
            gzFile fp1_tmp = gzopen(opt.reads_filename1, "r");
            auto ks1_tmp = make_ikstream(fp1_tmp, gzread);
            gzFile fp2_tmp = gzopen(opt.reads_filename2, "r");
            auto ks2_tmp = make_ikstream(fp2_tmp, gzread);
            auto r1_tmp = est_read_length(ks1_tmp, 500);
            auto r2_tmp = est_read_length(ks2_tmp, 500);
            gzclose(fp1_tmp);
            gzclose(fp2_tmp);
            map_param.r = (r1_tmp + r2_tmp) / 2;
        }
    }

    if (map_param.r <= 75) { // based on params for 100
        if (!opt.k_set){
            map_param.k = 20;
        }
        if ( (!opt.s_set ) ){
            map_param.s = map_param.k - 4; // Update default s to k - 4 if user has not set s parameter
        }
        map_param.l = -4;
        map_param.u = 2;
    } else if (map_param.r <= 125) { // based on params for 100
        if (!opt.k_set){
            map_param.k = 20;
        }
        if ( (!opt.s_set ) ){
            map_param.s = map_param.k - 4; // Update default s to k - 4 if user has not set s parameter
        }
        map_param.l = -2;
        map_param.u = 2;
    } else if ((map_param.r > 125) && (map_param.r <= 175)) { // based on params for 150
        if (!opt.k_set) {
            map_param.k = 20;
        }
        if ( (!opt.s_set ) ){
            map_param.s = map_param.k - 4; // Update default s to k - 4 if user has not set s parameter
        }
        map_param.l = 1;
        map_param.u = 7;
    } else if ((map_param.r > 175) && (map_param.r <= 275)) { // based on params for 200 and 250
        if (!opt.k_set) {
            map_param.k = 20;
        }
        if ( (!opt.s_set ) ){
            map_param.s = map_param.k - 4; // Update default s to k - 4 if user has not set s parameter
        }
        map_param.l = 4;
        map_param.u = 13;
    } else if ((map_param.r > 275) && (map_param.r <= 375)) { // based on params for 300
        if (!opt.k_set) {
            map_param.k = 22;
        }
        if ( (!opt.s_set ) ){
            map_param.s = map_param.k - 4; // Update default s to k - 4 if user has not set s parameter
        }
        map_param.l = 2;
        map_param.u = 12;
    } else{
        if (!opt.k_set) {
            map_param.k = 23;
        }
        if ( (!opt.s_set ) ){
            map_param.s = map_param.k - 6; // Update default s to k - 4 if user has not set s parameter
        }
        map_param.l = 2;
        map_param.u = 12;
    }


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

    map_param.w_min = map_param.k/(map_param.k-map_param.s+1) + map_param.l > 1 ? map_param.k/(map_param.k-map_param.s+1) + map_param.l : 1;
    map_param.w_max = map_param.k/(map_param.k-map_param.s+1) + map_param.u;
    map_param.t_syncmer = (map_param.k-map_param.s)/2 + 1;

    alignment_params aln_params;
    aln_params.match = opt.A;
    aln_params.mismatch = opt.B;
    aln_params.gap_open = opt.O;
    aln_params.gap_extend = opt.E;

    std::cerr << "Using" << std::endl;
    std::cerr << "k: " << map_param.k << std::endl;
    std::cerr << "s: " << map_param.s << std::endl;
    std::cerr << "w_min: " << map_param.w_min << std::endl;
    std::cerr << "w_max: " << map_param.w_max << std::endl;
    std::cerr << "Read length (r): " << map_param.r << std::endl;
    std::cerr << "Maximum seed length: " << map_param.max_dist + map_param.k << std::endl;
    std::cerr << "Threads: " << opt.n_threads << std::endl;
    std::cerr << "R: " << map_param.R << std::endl;
    std::cerr << "Expected [w_min, w_max] in #syncmers: [" << map_param.w_min << ", " << map_param.w_max << "]" << std::endl;
    std::cerr << "Expected [w_min, w_max] in #nucleotides: [" << (map_param.k-map_param.s+1)*map_param.w_min << ", " << (map_param.k-map_param.s+1)*map_param.w_max << "]" << std::endl;
    std::cerr << "A: " << opt.A << std::endl;
    std::cerr << "B: " << opt.B << std::endl;
    std::cerr << "O: " << opt.O << std::endl;
    std::cerr << "E: " << opt.E << std::endl;

//    assert(k <= (w/2)*w_min && "k should be smaller than (w/2)*w_min to avoid creating short strobemers");
    assert(map_param.k > 7 && "You should really not use too small strobe size!");
    assert(map_param.k <= 32 && "k have to be smaller than 32!");
    assert( ( map_param.s <= map_param.k ) && " s have to be smaller or equal to k!");
    assert( ( (map_param.k-map_param.s) % 2 == 0) && " k - s have to be an even number to create canonical syncmers. Set s to e.g., k-2, k-4, k-6, k-8.");
    assert(map_param.max_dist <= 255 && " -m (maximum seed length have to be smaller than 255 + k.");



    //////////// CREATE INDEX OF REF SEQUENCES /////////////////

    // Record index creation start time

    auto start = std::chrono::high_resolution_clock::now();
    auto start_read_refs = std::chrono::high_resolution_clock::now();
    std::vector<std::string> ref_seqs;
    std::vector<unsigned int> ref_lengths;
    uint64_t total_ref_seq_size;
    idx_to_acc acc_map;
    total_ref_seq_size = read_references(ref_seqs, ref_lengths, acc_map, opt.ref_filename);
    auto finish_read_refs = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_read_refs = finish_read_refs - start_read_refs;
    std::cerr << "Time reading references: " << elapsed_read_refs.count() << " s\n" <<  std::endl;

    if (total_ref_seq_size == 0) {
        std::cerr << "No reference sequences found, aborting.." << std::endl;
        return 1;
    }

    auto start_flat_vector = std::chrono::high_resolution_clock::now();

    mers_vector flat_vector;
    int approx_vec_size = total_ref_seq_size / (map_param.k-map_param.s+1);
    std::cerr << "ref vector approximate size: " << approx_vec_size << std::endl;
    flat_vector.reserve(approx_vec_size);
    unsigned int mer_cnt = 0;
    for(size_t i = 0; i < ref_seqs.size(); ++i)
    {
        mers_vector randstrobes2; // pos, chr_id, kmer hash value
//        std::cerr << i << " " << i_mod << std::endl;
        randstrobes2 = seq_to_randstrobes2(map_param.n, map_param.k, map_param.w_min, map_param.w_max, ref_seqs[i], i, map_param.s, map_param.t_syncmer, map_param.q, map_param.max_dist);
        for (auto &t : randstrobes2)
        {
            flat_vector.push_back(t);
        }
    }
    std::cerr << "Ref vector actual size: " << flat_vector.size() << std::endl;
    flat_vector.shrink_to_fit();

    auto finish_generating_seeds = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_generating_seeds = finish_generating_seeds - start_flat_vector;
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

    uint64_t unique_mers = 0;
    auto start_sorting = std::chrono::high_resolution_clock::now();
//    std::cerr << "Reserving flat vector size: " << approx_vec_size << std::endl;
//    all_mers_vector_tmp.reserve(approx_vec_size); // reserve size corresponding to sum of lengths of all sequences divided by expected sampling
    process_flat_vector(flat_vector, unique_mers);
    auto finish_sorting = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_sorting_seeds = finish_sorting - start_sorting;
    std::cerr << "Time sorting seeds: " << elapsed_sorting_seeds.count() << " s\n" <<  std::endl;
    std::cerr << "Unique strobemers: " << unique_mers  <<  std::endl;

    auto finish_flat_vector = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_flat_vector = finish_flat_vector - start_flat_vector;
    std::cerr << "Total time generating flat vector: " << elapsed_flat_vector.count() << " s\n" <<  std::endl;

    auto start_hash_index = std::chrono::high_resolution_clock::now();
    kmer_lookup mers_index; // k-mer -> (offset in flat_vector, occurence count )
    mers_index.reserve(unique_mers);
    map_param.filter_cutoff = index_vector(flat_vector, mers_index, map_param.f); // construct index over flat array
    auto finish_hash_index = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_hash_index = finish_hash_index - start_hash_index;
    std::cerr << "Total time generating hash table index: " << elapsed_hash_index.count() << " s\n" <<  std::endl;

//    mers_vector_reduced all_mers_vector;
//    all_mers_vector = remove_kmer_hash_from_flat_vector(flat_vector);
    /* destroy vector */
//    flat_vector.clear();


////////////////////////////////////////////////////////////////////////


//    std::cerr << "Wrote index to disc" << std::endl;

//    std::chrono::milliseconds timespan(10000); // or whatever
//    std::this_thread::sleep_for(timespan);
    //////////////////////////////////////////////////////////////////////////

    // Record index creation end time
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cerr << "Total time indexing: " << elapsed.count() << " s\n" <<  std::endl;

    if (opt.index_log){
        std::cerr << "Printing log stats" << std::endl;
        print_diagnostics(flat_vector, mers_index, opt.logfile_name, map_param.k, map_param.max_dist + map_param.k);
        std::cerr << "Finished printing log stats" << std::endl;

    }

//    std::chrono::milliseconds timespan2(1000000); // or whatever
//    std::this_thread::sleep_for(timespan2);

    ///////////////////////////// MAP ///////////////////////////////////////

    // Record matching time
    auto start_aln_part = std::chrono::high_resolution_clock::now();

    logging_variables tot_log_vars;
//    std::ifstream query_file(reads_filename);

    map_param.rescue_cutoff = map_param.R < 100 ? map_param.R*map_param.filter_cutoff : 1000;
    std::cerr << "Using rescue cutoff: " << map_param.rescue_cutoff <<  std::endl;

    std::streambuf *buf;
    std::ofstream of;

    if(!opt.write_to_stdout) {
        of.open(opt.output_file_name);
        buf = of.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }

    std::ostream out(buf);
//    std::ofstream out;
//    out.open(opt.output_file_name);

//    std::stringstream sam_output;
//    std::stringstream paf_output;

    if (map_param.is_sam_out) {
        int nr_ref = acc_map.size();
        for (int i = 0; i < nr_ref; ++i) {
//        for (auto &it : acc_map) {
//            out << "@SQ\tSN:" << it.second << "\tLN:" << ref_lengths[it.first] << "\n";
            out << "@SQ\tSN:" << acc_map[i] << "\tLN:" << ref_lengths[i] << "\n";
        }
        out << "@PG\tID:strobealign\tPN:strobealign\tVN:0.7.1\tCL:strobealign\n";
    }

    std::unordered_map<std::thread::id, logging_variables> log_stats_vec(opt.n_threads);


    if(opt.is_SE) {
        std::cerr << "Running SE mode" <<  std::endl;
        //    KSeq record;
        gzFile fp = gzopen(opt.reads_filename1, "r");
        auto ks = make_ikstream(fp, gzread);

        ////////// ALIGNMENT START //////////
        /////////////////////////////////////

        int input_chunk_size = 100000;
        // Create Buffers
        InputBuffer input_buffer = { {}, {}, {}, {}, {}, ks, ks, false, 0, input_chunk_size};
        OutputBuffer output_buffer = { {}, {}, {}, 0, out};

        std::vector<std::thread> workers;
        for (int i = 0; i < opt.n_threads; ++i) {
            std::thread consumer(perform_task_SE, std::ref(input_buffer), std::ref(output_buffer),
                                 std::ref(log_stats_vec), std::ref(aln_params),
                                 std::ref(map_param), std::ref(ref_lengths), std::ref(ref_seqs),
                                 std::ref(mers_index), std::ref(flat_vector), std::ref(acc_map) );
            workers.push_back(std::move(consumer));
        }

        for (size_t i = 0; i < workers.size(); ++i) {
            workers[i].join();
        }

        std::cerr << "Done!\n";
        /////////////////////////////////////
        /////////////////////////////////////

        gzclose(fp);
    }
    else{
        std::cerr << "Running PE mode" <<  std::endl;
        gzFile fp1 = gzopen(opt.reads_filename1, "r");
        auto ks1 = make_ikstream(fp1, gzread);
        gzFile fp2 = gzopen(opt.reads_filename2, "r");
        auto ks2 = make_ikstream(fp2, gzread);
        std::unordered_map<std::thread::id, i_dist_est> isize_est_vec(opt.n_threads);

        ////////// ALIGNMENT START //////////
        /////////////////////////////////////

        int input_chunk_size = 100000;
        // Create Buffers
        InputBuffer input_buffer = { {}, {}, {}, {}, {}, ks1, ks2, false, 0, input_chunk_size};
        OutputBuffer output_buffer = { {}, {}, {}, 0, out};

        std::vector<std::thread> workers;
        for (int i = 0; i < opt.n_threads; ++i) {
            std::thread consumer(perform_task_PE, std::ref(input_buffer), std::ref(output_buffer),
                                 std::ref(log_stats_vec), std::ref(isize_est_vec), std::ref(aln_params),
                                 std::ref(map_param), std::ref(ref_lengths), std::ref(ref_seqs),
                                 std::ref(mers_index), std::ref(flat_vector), std::ref(acc_map) );
            workers.push_back(std::move(consumer));
        }

        for (size_t i = 0; i < workers.size(); ++i) {
            workers[i].join();
        }

        std::cerr << "Done!\n";
        /////////////////////////////////////
        /////////////////////////////////////

        gzclose(fp1);
        gzclose(fp2);
    }

    for (auto &it : log_stats_vec) {
        auto thread_id = it.first;
        auto log_vars = it.second;
        tot_log_vars.tot_all_tried += log_vars.tot_all_tried;
        tot_log_vars.tot_ksw_aligned += log_vars.tot_ksw_aligned;
        tot_log_vars.tot_rescued += log_vars.tot_rescued;
        tot_log_vars.did_not_fit += log_vars.did_not_fit;
        tot_log_vars.tried_rescue += log_vars.tried_rescue;

        tot_log_vars.tot_read_file += log_vars.tot_read_file;
        tot_log_vars.tot_construct_strobemers += log_vars.tot_construct_strobemers;
        tot_log_vars.tot_find_nams += log_vars.tot_find_nams;
        tot_log_vars.tot_time_rescue += log_vars.tot_time_rescue;
        tot_log_vars.tot_sort_nams += log_vars.tot_sort_nams;
        tot_log_vars.tot_rc += log_vars.tot_rc;
        tot_log_vars.tot_extend += log_vars.tot_extend;
        tot_log_vars.tot_write_file += log_vars.tot_write_file;
    }

    std::cerr << "Total mapping sites tried: " << tot_log_vars.tot_all_tried << std::endl;
    std::cerr << "Total calls to ssw: " << tot_log_vars.tot_ksw_aligned << std::endl;
    std::cerr << "Calls to ksw (rescue mode): " << tot_log_vars.tot_rescued << std::endl;
    std::cerr << "Did not fit strobe start site: " << tot_log_vars.did_not_fit  << std::endl;
    std::cerr << "Tried rescue: " << tot_log_vars.tried_rescue  << std::endl;
    // Record mapping end time
    auto finish_aln_part = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tot_aln_part = finish_aln_part - start_aln_part;
    std::cerr << "Total time mapping: " << tot_aln_part.count() << " s." <<  std::endl;
    std::cerr << "Total time reading read-file(s): " << tot_log_vars.tot_read_file.count()/opt.n_threads << " s." <<  std::endl;
    std::cerr << "Total time creating strobemers: " << tot_log_vars.tot_construct_strobemers.count()/opt.n_threads << " s." <<  std::endl;
    std::cerr << "Total time finding NAMs (non-rescue mode): " << tot_log_vars.tot_find_nams.count()/opt.n_threads  << " s." <<  std::endl;
    std::cerr << "Total time finding NAMs (rescue mode): " << tot_log_vars.tot_time_rescue.count()/opt.n_threads  << " s." <<  std::endl;
//    std::cerr << "Total time finding NAMs ALTERNATIVE (candidate sites): " << tot_find_nams_alt.count()/opt.n_threads  << " s." <<  std::endl;
    std::cerr << "Total time sorting NAMs (candidate sites): " << tot_log_vars.tot_sort_nams.count()/opt.n_threads  << " s." <<  std::endl;
    std::cerr << "Total time reverse compl seq: " << tot_log_vars.tot_rc.count()/opt.n_threads  << " s." <<  std::endl;
    std::cerr << "Total time base level alignment (ssw): " << tot_log_vars.tot_extend.count()/opt.n_threads  << " s." <<  std::endl;
    std::cerr << "Total time writing alignment to files: " << tot_log_vars.tot_write_file.count() << " s." <<  std::endl;

    //////////////////////////////////////////////////////////////////////////


    /////////////////////// FIND AND OUTPUT NAMs ///////////////////////////////





}

