#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>
#include <chrono>  // for high_resolution_clock
#include <omp.h>
#include <zlib.h>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <inttypes.h>

#include "source/kseq++.hpp"
using namespace klibpp;
#include "source/robin_hood.h"
#include "source/index.hpp"
#include "source/ksw2.h"
#include "source/ssw_cpp.h"
#include "source/pc.hpp"
#include "source/aln.hpp"

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
    file.close();
    return total_ref_seq_size;
}



static inline void print_diagnostics(mers_vector &ref_mers, kmer_lookup &mers_index, std::string logfile_name, int k) {
    // Prins to csv file the statistics on the number of seeds of a particular length and what fraction of them them are unique in the index:
    // format:
    // seed_length, count, percentage_unique
    //

    int max_size = 100000;
    std::vector<int> log_count(max_size,0);  // stores count and each index represents the length
    std::vector<int> log_unique(max_size,0); // stores count unique and each index represents the length
    std::vector<int> log_repetitive(max_size,0); // stores count unique and each index represents the length


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

//    std::cerr << "Here" << std::endl;

    // printing
    std::ofstream log_file;
    log_file.open(logfile_name);
    for (int i=0 ; i < log_count.size(); ++i) {
        if (log_count[i] > 0) {
            log_file << i << ',' << log_count[i] << ',' << (float) log_unique[i] / (float) log_count[i] << ',' << (float) log_repetitive[i] / (float) log_count[i] << std::endl;
        }
    }
    log_file.close();

    }




void print_usage() {
    std::cerr << "\n";
    std::cerr << "StrobeAlign VERSION 0.6.2 \n";
    std::cerr << "\n";
    std::cerr << "StrobeAlign [options] <ref.fa> <reads1.fast[a/q.gz]> [reads2.fast[a/q.gz]]\n";
    std::cerr << "options:\n";
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
    std::cerr << "\t-r INT Approximate read length. Sets suitable parameters for -k, -l, -u and -q. [150] \n";
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



int main (int argc, char **argv)
{

    if (argc < 3) {
        print_usage();
        return 0;
    }

    // Default parameters
    std::string choice = "randstrobes";
    int A = 2;
    int B = 8;
    int O = 12;
    int E = 1;

    mapping_params map_param;
    map_param.max_secondary = 0;
    bool is_sam_out = true; // true = align, false=map, default is_sam_out is align
    int n_threads = 3;
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
    int max_seed_len;
    map_param.max_dist = map_param.r - 50 < 255 ? map_param.r-50 : 255;
    bool r_set = false;
    bool max_seed_len_set = false;
    bool write_to_stdout = true;
    std::string output_file_name;
    std::string logfile_name = "log.csv";
    bool s_set = false;
    bool index_log = false;
    int opn = 1;
    while (opn < argc) {
        bool flag = false;
        if (argv[opn][0] == '-') {
//            if (argv[opn][1] == 'n') {
//                n = std::stoi(argv[opn + 1]);
//                opn += 2;
//                flag = true;
//            } else
            if (argv[opn][1] == 't') {
                n_threads = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'k') {
                map_param.k = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'o') {
                output_file_name = argv[opn + 1];
                opn += 2;
                flag = true;
                write_to_stdout = false;
            } else if (argv[opn][1] == 'L') {
                logfile_name = argv[opn + 1];
                index_log = true;
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 's') {
                map_param.s = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
                s_set = true;
            } else if (argv[opn][1] == 'f') {
                map_param.f = std::stof(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'x') {
                is_sam_out = false;
                opn += 1;
                flag = true;
            } else if (argv[opn][1] == 'R') {
                map_param.R = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'l') {
                map_param.l = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'u') {
                map_param.u = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'c') {
                map_param.c = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'r') {
                map_param.r = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
                r_set = true;
            } else if (argv[opn][1] == 'm') {
                max_seed_len = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
                max_seed_len_set = true;
            } else if (argv[opn][1] == 'S') {
                map_param.dropoff_threshold = std::stof(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'M') {
                map_param.maxTries = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'A') {
                A = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'B') {
                B = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'O') {
                O = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'E') {
                E = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'N') {
                map_param.max_secondary = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            }

            else {
                print_usage();
            }
        }
        if (!flag)
            break;
    }

    if (r_set) {
        if (map_param.r <= 125) { // based on params for 100
            map_param.k = 20;
            map_param.l = -2;
            map_param.u = 2;
        } else if ((map_param.r > 125) && (map_param.r <= 175)) { // based on params for 150
            map_param.k = 20;
            map_param.l = 1;
            map_param.u = 7;
        } else if ((map_param.r > 175) && (map_param.r <= 275)) { // based on params for 200 and 250
            map_param.k = 20;
            map_param.l = 4;
            map_param.u = 13;
        } else { // based on params for 300
            map_param.k = 22;
            map_param.l = 2;
            map_param.u = 12;
        }
    }

    if (!max_seed_len_set){
        map_param.max_dist = map_param.r - 70 > map_param.k ? map_param.r - 70 : map_param.k;
    } else {
        map_param.max_dist = max_seed_len - map_param.k; //convert to distance in start positions
    }

    if ( (!s_set ) ){
        map_param.s = map_param.k - 4; // Update default s to k - 4 if user has not set s parameter
    }

    if ( (map_param.c < 64) && (map_param.c > 0)){
        map_param.q = pow (2, map_param.c) - 1;
    } else{
        std::cerr << "Warning wrong value for parameter c, setting c=8" << std::endl;
        map_param.q = pow (2, 8) - 1;
    }
    omp_set_num_threads(n_threads); // set number of threads in "parallel" blocks
    map_param.w_min = map_param.k/(map_param.k-map_param.s+1) + map_param.l > 1 ? map_param.k/(map_param.k-map_param.s+1) + map_param.l : 1;
    map_param.w_max = map_param.k/(map_param.k-map_param.s+1) + map_param.u;
    map_param.t_syncmer = (map_param.k-map_param.s)/2 + 1;
    map_param.is_sam_out = is_sam_out;

    alignment_params aln_params;
    aln_params.match = A;
    aln_params.mismatch = B;
    aln_params.gap_open = O;
    aln_params.gap_extend = E;

    std::cerr << "Using" << std::endl;
    std::cerr << "k: " << map_param.k << std::endl;
    std::cerr << "s: " << map_param.s << std::endl;
    std::cerr << "w_min: " << map_param.w_min << std::endl;
    std::cerr << "w_max: " << map_param.w_max << std::endl;
    std::cerr << "maximum seed length: " << map_param.max_dist + map_param.k << std::endl;
    std::cerr << "threads: " << n_threads << std::endl;
    std::cerr << "R: " << map_param.R << std::endl;
    std::cerr << "[w_min, w_max] under thinning w roughly corresponds to sampling from downstream read coordinates (expected values): [" << (map_param.k-map_param.s+1)*map_param.w_min << ", " << (map_param.k-map_param.s+1)*map_param.w_max << "]" << std::endl;
    std::cerr << "A: " << A << std::endl;
    std::cerr << "B: " << B << std::endl;
    std::cerr << "O: " << O << std::endl;
    std::cerr << "E: " << E << std::endl;

//    assert(k <= (w/2)*w_min && "k should be smaller than (w/2)*w_min to avoid creating short strobemers");
    assert(map_param.k > 7 && "You should really not use too small strobe size!");
    assert(map_param.k <= 32 && "k have to be smaller than 32!");
    assert( ( map_param.s <= map_param.k ) && " s have to be smaller or equal to k!");
    assert( ( (map_param.k-map_param.s) % 2 == 0) && " k - s have to be an even number to create canonical syncmers. Set s to e.g., k-2, k-4, k-6, k-8.");
    assert(map_param.max_dist <= 255 && " -m (maximum seed length have to be smaller than 255 + k in v0.4 and up. If you need longer seeds, use v0.3");

//    assert(n == 2 && "Currently only n=2 is implemented");
    // File name to reference
    std::string ref_filename = argv[opn];
//    opn++;
//    const char *reads_filename = argv[opn];
//    opn++;
//    const char *reads_filename_PE2 = argv[opn];
    opn++;
    const char *reads_filename1;
    const char *reads_filename2;
    bool is_SE = false;
    if (opn == argc - 1) {
        reads_filename1 = argv[opn];
        is_SE = true;
    } else if (opn == argc - 2) {
        reads_filename1 = argv[opn];
        opn++;
        reads_filename2 = argv[opn];
    } else {
        print_usage();
        return 0;
    }


    //////////// CREATE INDEX OF REF SEQUENCES /////////////////

    // Record index creation start time

    auto start = std::chrono::high_resolution_clock::now();
    auto start_read_refs = std::chrono::high_resolution_clock::now();
    std::vector<std::string> ref_seqs;
    std::vector<unsigned int> ref_lengths;
    uint64_t total_ref_seq_size;
    idx_to_acc acc_map;
    total_ref_seq_size = read_references(ref_seqs, ref_lengths, acc_map, ref_filename);
    auto finish_read_refs = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_read_refs = finish_read_refs - start_read_refs;
    std::cerr << "Time reading references: " << elapsed_read_refs.count() << " s\n" <<  std::endl;

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
//    std::vector<std::vector<std::tuple<uint64_t, unsigned int, unsigned int, unsigned int>>> vector_per_ref_chr(n_threads);
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

    if (index_log){
        std::cerr << "Printing log stats" << std::endl;
        print_diagnostics(flat_vector, mers_index, logfile_name, map_param.k);
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

    if(!write_to_stdout) {
        of.open(output_file_name);
        buf = of.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }

    std::ostream out(buf);
//    std::ofstream out;
//    out.open(output_file_name);

//    std::stringstream sam_output;
//    std::stringstream paf_output;

    if (is_sam_out) {
        for (auto &it : acc_map) {
            out << "@SQ\tSN:" << it.second << "\tLN:" << ref_lengths[it.first] << "\n";
        }
        out << "@PG\tID:strobealign\tPN:strobealign\tVN:0.6.2\tCL:strobealign\n";
    }

    if(is_SE) {
        ;
//        std::cerr << "Running SE mode" <<  std::endl;
//        //    KSeq record;
//        gzFile fp = gzopen(reads_filename1, "r");
//        auto ks = make_ikstream(fp, gzread);
//        int n_q_chunk_size = 1000000;
//        KSeq record;
//        std::string seq, seq_rc;
//        unsigned int q_id = 0;
//        std::pair<float, int> info;
//        mers_vector_read query_mers; // pos, chr_id, kmer hash value
//        std::vector<nam> nams; // (r_id, r_pos_start, r_pos_end, q_pos_start, q_pos_end)
////        std::vector<alignment> alignments;
////        alignments.reserve(maxTries);
//        robin_hood::unordered_map< unsigned int, std::vector<hit>> hits_per_ref;
//        std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_fw;
//        std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_rc;
//        hits_per_ref.reserve(100);
//        hits_fw.reserve(5000);
//        hits_rc.reserve(5000);
////    mers_vector_read query_mers_rc; // pos, chr_id, kmer hash value
//
////    std::vector<std::stringstream> output_streams(n_threads);
//        std::vector<std::string> output_streams(n_threads);
//        for (int i = 0; i < n_threads; ++i) {
//            output_streams[i].reserve((n_q_chunk_size / n_threads + 1) *
//                                      450); // Reserve sufficient space for appending multiple SAM records (400 is an upper setimate on the number of characters for each sam record of a 200-300bp read)
//        }
//
//        while (ks) {
//
//            auto read_start = std::chrono::high_resolution_clock::now();
//            auto records = ks.read(n_q_chunk_size);  // read a chunk of 1000000 records
//            auto read_finish = std::chrono::high_resolution_clock::now();
//            tot_log_vars.tot_read_file += read_finish - read_start;
////        std::chrono::duration<double> elapsed_read = read_finish - read_start;
////        std::cerr << "Total time reading from file: " << elapsed_read.count() << " s\n" <<  std::endl;
//
//            int n_it = records.size();
//            std::cerr << "Mapping chunk of " << n_it << " query sequences... " << std::endl;
//            #pragma omp parallel for num_threads(n_threads) shared(aln_params, output_streams, tot_log_vars, out, q_id) private(record, seq_rc, query_mers, nams, hits_per_ref, info)
//            for (int i = 0; i < n_it; ++i) {
//                auto record = records[i];
//                // generate mers here
//                auto strobe_start = std::chrono::high_resolution_clock::now();
//                query_mers = seq_to_randstrobes2_read(map_param.n, map_param.k, map_param.w_min, map_param.w_max, record.seq, q_id, map_param.s, map_param.t_syncmer, map_param.q, map_param.max_dist);
//                auto strobe_finish = std::chrono::high_resolution_clock::now();
//                tot_log_vars.tot_construct_strobemers += strobe_finish - strobe_start;
//
////            // Find NAMs alternative function
////            auto nam_alt_start = std::chrono::high_resolution_clock::now();
////            std::vector<nam> nams_alt; // (r_id, r_pos_start, r_pos_end, q_pos_start, q_pos_end)
////            nams_alt = find_nams_alt(query_mers, all_mers_vector, mers_index, k, ref_seqs, record.seq, hit_upper_window_lim,
////                                       filter_cutoff);
////            auto nam_alt_finish = std::chrono::high_resolution_clock::now();
////            tot_find_nams_alt += nam_alt_finish - nam_alt_start;
//
//                // Find NAMs
////                std::cerr << "mapping " << record.name << std::endl;
//                auto nam_start = std::chrono::high_resolution_clock::now();
//                info = find_nams(nams, hits_per_ref, query_mers, flat_vector, mers_index, map_param.k, ref_seqs, record.seq, map_param.filter_cutoff);
//                hits_per_ref.clear();
//                auto nam_finish = std::chrono::high_resolution_clock::now();
//                tot_log_vars.tot_find_nams += nam_finish - nam_start;
//
//                if (map_param.R > 1) {
//                    auto rescue_start = std::chrono::high_resolution_clock::now();
//                    if ((nams.size() == 0) || (info.first < 0.7)) {
//                        tot_log_vars.tried_rescue += 1;
//                        nams.clear();
////                    std::cerr << "Rescue is_sam_out: " << record.name <<  std::endl;
//                        info = find_nams_rescue(hits_fw, hits_rc, nams, hits_per_ref, query_mers, flat_vector, mers_index, map_param.k, ref_seqs,
//                                                record.seq, map_param.rescue_cutoff);
//                        hits_per_ref.clear();
//                        hits_fw.clear();
//                        hits_rc.clear();
////                    std::cerr << "Found: " << nams.size() <<  std::endl;
//                    }
//                    auto rescue_finish = std::chrono::high_resolution_clock::now();
//                    tot_log_vars.tot_time_rescue += rescue_finish - rescue_start;
//                }
//
//
//                //Sort hits on score
//                auto nam_sort_start = std::chrono::high_resolution_clock::now();
//                std::sort(nams.begin(), nams.end(), score);
//                auto nam_sort_finish = std::chrono::high_resolution_clock::now();
//                tot_log_vars.tot_sort_nams += nam_sort_finish - nam_sort_start;
//
//                auto extend_start = std::chrono::high_resolution_clock::now();
//                if (!is_sam_out) {
//                    output_hits_paf(output_streams[omp_get_thread_num()], nams, record.name, acc_map, map_param.k,
//                                    record.seq.length(), ref_lengths);
////                output_streams[omp_get_thread_num()].append(paf_output.str()); // << paf_output.str();
//                } else {
////                auto rc_start = std::chrono::high_resolution_clock::now();
////                auto rc_finish = std::chrono::high_resolution_clock::now();
////                tot_rc += rc_finish - rc_start;
////                    align_SE_multimap_deprecated(aln_params, output_streams[omp_get_thread_num()], nams, record.name, acc_map, k, record.seq.length(),
////                             ref_lengths, ref_seqs, record.seq,
////                             tot_ksw_aligned, tot_all_tried, dropoff_threshold, did_not_fit, maxTries, max_secondary);
//                    if (map_param.max_secondary > 0){
//                        // I created an entire new function here, dupliocating a lot of the code as outputting secondary hits is has some overhead to the
//                        // original align_SE function (storing a vector of hits and sorting them)
//                        // Such overhead is not present in align_PE - which implements both options in the same function.
//                        align_SE_secondary_hits(aln_params,output_streams[omp_get_thread_num()], nams, record.name, acc_map, map_param.k, record.seq.length(),
//                                 ref_lengths, ref_seqs, record.seq, record.qual, tot_log_vars, map_param.dropoff_threshold, map_param.maxTries, map_param.max_secondary);
//                    } else {
//                        align_SE(aln_params,output_streams[omp_get_thread_num()], nams, record.name, acc_map, map_param.k, record.seq.length(),
//                                 ref_lengths, ref_seqs, record.seq, record.qual, tot_log_vars,  map_param.dropoff_threshold, map_param.maxTries);
//                    }
//
////                output_streams[omp_get_thread_num()] << sam_output.str();
//                }
//                auto extend_finish = std::chrono::high_resolution_clock::now();
//                tot_log_vars.tot_extend += extend_finish - extend_start;
//                q_id++;
//                nams.clear();
//            }
//            // Output results
//            auto write_start = std::chrono::high_resolution_clock::now();
//            for (int i = 0; i < n_threads; ++i) {
//                out << output_streams[i];
//                output_streams[i].clear();
//            }
//            auto write_finish = std::chrono::high_resolution_clock::now();
//            tot_log_vars.tot_write_file += write_finish - write_start;
//        }
//        gzclose(fp);
    }
    else{
        std::cerr << "Running PE mode" <<  std::endl;
        //    KSeq record;
//        int max_nam_n_hits1,max_nam_n_hits2;
//        std::vector<int> isizes;
        gzFile fp1 = gzopen(reads_filename1, "r");
        auto ks1 = make_ikstream(fp1, gzread);
        gzFile fp2 = gzopen(reads_filename2, "r");
        auto ks2 = make_ikstream(fp2, gzread);
        int n_q_chunk_size = 1000000;
        i_dist_est isize_est;
        logging_variables log_vars;
//        std::vector<std::string> output_streams(n_threads);
        robin_hood::unordered_map<std::thread::id, logging_variables> log_stats_vec(n_threads);
        robin_hood::unordered_map<std::thread::id, i_dist_est> isize_est_vec(n_threads);
//        for (int i = 0; i < n_threads; ++i) {
////            output_streams[i].reserve((n_q_chunk_size / n_threads + 1) *
////                                      600); // Reserve sufficient space for appending multiple SAM records (600 is an estimate on the number of characters for each sam record of a 200-300bp read)
//            i_dist_est isize_est;
//            logging_variables log_vars;
//            log_stats_vec[i] = log_vars;
//            isize_est_vec[i] = isize_est;
//
//        }

        ////// REAL PRODUCER-CONSUMER START /////

//        int n_threads = 4;
        int input_chunk_size = 100000;
//        int read_len = 100;
        std::cerr << "Nr threads: " << n_threads << "\n";

        // Create Buffers
        InputBuffer input_buffer = { {}, {}, {}, {}, {}, ks1, ks2, false, 0, input_chunk_size};

//        int64_t reserve_size = (OUTPUT_BUFFER_CAPACITY)*input_chunk_size * 2 * 4 * map_param.r;
//        std::cerr << "OUTPUT_BUFFER_CAPACITY: " << OUTPUT_BUFFER_CAPACITY << " input_chunk_size: " << input_chunk_size <<   " reserve_size: " << reserve_size << "\n";
        OutputBuffer output_buffer;

        std::vector<std::thread> workers;
        for (int i = 0; i < n_threads; ++i) {
            std::thread consumer(perform_task, std::ref(input_buffer), std::ref(output_buffer),
                                 std::ref(log_stats_vec), std::ref(isize_est_vec), std::ref(aln_params),
                                 std::ref(map_param), std::ref(ref_lengths), std::ref(ref_seqs),
                                 std::ref(mers_index), std::ref(flat_vector), std::ref(acc_map) );
            workers.push_back(std::move(consumer));
        }

        for (size_t i = 0; i < workers.size(); ++i) {
            workers[i].join();
        }

//        std::cerr << "LAST BUFFER SIZE: " << output_buffer.buffer_size <<  std::endl;
//        if (output_buffer.buffer_size > 0 ){ // write last set of records
//            std::cerr << "Final writing to output " <<output_buffer.buffer_size << " " << input_buffer.buffer_size  << " threadID: " << std::this_thread::get_id() << " " << input_buffer.finished_reading << std::endl;
//            output_buffer.output_records(std::this_thread::get_id()); // Implement write here
//        }

        std::cerr << "Done!\n";
        /////////////////////////////////////


//        auto read_start = std::chrono::high_resolution_clock::now();
//        auto records1 = ks1.read(n_q_chunk_size);  // read a chunk of 1000000 records
//        auto records2 = ks2.read(n_q_chunk_size);  // read a chunk of 1000000 records
//        auto read_finish = std::chrono::high_resolution_clock::now();
//        log_vars.tot_read_file += read_finish - read_start;
//        std::vector<KSeq> new_records1;
//        std::vector<KSeq> new_records2;
//        while (!records1.empty()) {
//            int n_it = records1.size();
//            std::cerr << "Mapping chunk of " << n_it << " query sequences... " << std::endl;
//
//            #pragma omp parallel sections shared(output_streams, log_stats_vec, isize_est_vec) //private(aln_params, map_param)
//            {
//                // Producer (read-write IO)
//#pragma omp section
//                {
//                    auto read_start = std::chrono::high_resolution_clock::now();
//                    new_records1 = ks1.read(n_q_chunk_size);
//                    new_records2 = ks2.read(n_q_chunk_size);
//                    auto read_finish = std::chrono::high_resolution_clock::now();
//                    log_stats_vec[omp_get_thread_num()].tot_read_file += read_finish - read_start;
//
////                    // Output results
////                    auto write_start = std::chrono::high_resolution_clock::now();
////                    for (int i = 0; i < n_threads; ++i) {
////                        out << old_output_streams[i];
////                        old_output_streams[i].clear();
////                    }
////                    auto write_finish = std::chrono::high_resolution_clock::now();
////                    tot_write_file += write_finish - write_start;
//                }
//
//                // Consumers
//#pragma omp section
//                {
//
////                    #pragma omp parallel for num_threads(n_threads) shared(aln_params, output_streams, out, q_id, tot_all_tried, did_not_fit, tot_ksw_aligned, tried_rescue, sample_size, mu, sigma, V, SSE) private(record1, record2, seq_rc1, seq_rc2, query_mers1, query_mers2, nams1, nams2, hits_per_ref, joint_NAM_scores, info1, info2)
////                    #pragma omp for schedule(dynamic) private(record1, record2, seq_rc1, seq_rc2, query_mers1, query_mers2, nams1, nams2, hits_per_ref, joint_NAM_scores, info1, info2)
//                    for (int i = 0; i < n_it; ++i) {
//#pragma omp task firstprivate(i)  //private(aln_params, map_param)
//                        {
////#pragma omp critical
////                            {
////                                std::cerr << aln_params.match << aln_params.mismatch << aln_params.gap_open
////                                          << aln_params.gap_extend << std::endl;
////                                std::cerr << map_param.q << map_param.n << map_param.k << map_param.w_min
////                                          << map_param.w_max << map_param.s << map_param.t_syncmer << map_param.q
////                                          << map_param.max_dist << map_param.max_secondary
////                                          << map_param.dropoff_threshold << map_param.r << map_param.m << map_param.l
////                                          << map_param.u << map_param.c << map_param.f << map_param.S << map_param.M
////                                          << map_param.R << map_param.max_dist << map_param.maxTries
////                                          << map_param.max_seed_len << map_param.rescue_cutoff
////                                          << map_param.filter_cutoff << std::endl;
////                            }
//
//                            // Declare variables
//                            mers_vector_read query_mers1, query_mers2; // pos, chr_id, kmer hash value
//                            auto log_vars = log_stats_vec[omp_get_thread_num()];
//                            auto isize_est = isize_est_vec[omp_get_thread_num()];
//                            auto record1 = records1[i];
//                            auto record2 = records2[i];
//
//                            // generate mers here
//                            auto strobe_start = std::chrono::high_resolution_clock::now();
////                std::cerr << "Going in! " << std::endl;
//                            query_mers1 = seq_to_randstrobes2_read(map_param.n, map_param.k, map_param.w_min, map_param.w_max, record1.seq, 0, map_param.s, map_param.t_syncmer,
//                                                                   map_param.q,
//                                                                   map_param.max_dist);
////                std::cerr << "Lolz1 " << std::endl;
//                            query_mers2 = seq_to_randstrobes2_read(map_param.n, map_param.k, map_param.w_min, map_param.w_max, record2.seq, 0, map_param.s, map_param.t_syncmer,
//                                                                   map_param.q,
//                                                                   map_param.max_dist);
////                std::cerr << "Lolz2 " << std::endl;
//                            auto strobe_finish = std::chrono::high_resolution_clock::now();
//                            log_vars.tot_construct_strobemers += strobe_finish - strobe_start;
////                std::cerr << record1.name << " " << query_mers1.size() << std::endl;
////                std::cerr << record2.name << " " << query_mers2.size() << std::endl;
//
//                            // Find NAMs
//                            auto nam_start = std::chrono::high_resolution_clock::now();
//                            robin_hood::unordered_map< unsigned int, std::vector<hit>> hits_per_ref;
//                            hits_per_ref.reserve(100);
//                            std::vector<nam> nams1;
//                            std::vector<nam> nams2;
//                            std::pair<float, int> info1, info2;
//                            info1 = find_nams(nams1, hits_per_ref, query_mers1, flat_vector, mers_index, map_param.k, ref_seqs,
//                                              record1.seq, map_param.filter_cutoff);
//                            hits_per_ref.clear();
//                            info2 = find_nams(nams2, hits_per_ref, query_mers2, flat_vector, mers_index, map_param.k, ref_seqs,
//                                              record2.seq, map_param.filter_cutoff);
//                            hits_per_ref.clear();
//                            auto nam_finish = std::chrono::high_resolution_clock::now();
//                            log_vars.tot_find_nams += nam_finish - nam_start;
//
//
//                            if (map_param.R > 1) {
//                                auto rescue_start = std::chrono::high_resolution_clock::now();
//                                std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_fw;
//                                std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> hits_rc;
//                                if ((nams1.size() == 0) || (info1.first < 0.7)) {
//                                    hits_fw.reserve(5000);
//                                    hits_rc.reserve(5000);
//                                    log_vars.tried_rescue += 1;
//                                    nams1.clear();
////                        std::cerr << "Rescue is_sam_out read 1: " << record1.name << info1.first <<  std::endl;
//                                    info1 = find_nams_rescue(hits_fw, hits_rc, nams1, hits_per_ref, query_mers1,
//                                                             flat_vector,
//                                                             mers_index, map_param.k, ref_seqs,
//                                                             record1.seq, map_param.rescue_cutoff);
//                                    hits_per_ref.clear();
//                                    hits_fw.clear();
//                                    hits_rc.clear();
////                    std::cerr << "Found: " << nams.size() <<  std::endl;
//                                }
//
//                                if ((nams2.size() == 0) || (info2.first < 0.7)) {
//                                    hits_fw.reserve(5000);
//                                    hits_rc.reserve(5000);
//                                    log_vars.tried_rescue += 1;
//                                    nams2.clear();
////                        std::cerr << "Rescue is_sam_out read 2: " << record2.name << info2.first <<  std::endl;
//                                    info2 = find_nams_rescue(hits_fw, hits_rc, nams2, hits_per_ref, query_mers2,
//                                                             flat_vector,
//                                                             mers_index, map_param.k, ref_seqs,
//                                                             record2.seq, map_param.rescue_cutoff);
//                                    hits_per_ref.clear();
//                                    hits_fw.clear();
//                                    hits_rc.clear();
////                    std::cerr << "Found: " << nams.size() <<  std::endl;
//                                }
//                                auto rescue_finish = std::chrono::high_resolution_clock::now();
//                                log_vars.tot_time_rescue += rescue_finish - rescue_start;
//                            }
//
//
//
//                            //Sort hits based on start choordinate on query sequence
//                            auto nam_sort_start = std::chrono::high_resolution_clock::now();
//                            std::sort(nams1.begin(), nams1.end(), score);
//                            std::sort(nams2.begin(), nams2.end(), score);
//                            auto nam_sort_finish = std::chrono::high_resolution_clock::now();
//                            log_vars.tot_sort_nams += nam_sort_finish - nam_sort_start;
//
////                std::cerr << record1.name << std::endl;
////                for (auto &n : nams1){
////                    std::cerr << "NAM ORG: " << n.ref_id << ": (" << n.score << ", " << n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ")" << std::endl;
////                }
////                std::cerr << record2.name << std::endl;
////                for (auto &n : nams2){
////                    std::cerr << "NAM ORG: " << n.ref_id << ": (" << n.score << ", " << n.n_hits << ", " << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e  << ")" << std::endl;
////                }
//
//                            auto extend_start = std::chrono::high_resolution_clock::now();
//                            if (!is_sam_out) {
////                    output_hits_paf(output_streams[omp_get_thread_num()], nams1, record1.name, acc_map, k,
////                                    record1.seq.length(), ref_lengths);
//                                nam nam_read1;
//                                nam nam_read2;
//                                std::vector<std::tuple<int, nam, nam>> joint_NAM_scores;
//                                get_best_map_location(joint_NAM_scores, nams1, nams2, isize_est,
//                                                      nam_read1,
//                                                      nam_read2);
//                                output_hits_paf_PE(output_streams[omp_get_thread_num()], nam_read1, record1.name,
//                                                   acc_map,
//                                                   map_param.k,
//                                                   record1.seq.length(), ref_lengths);
//                                output_hits_paf_PE(output_streams[omp_get_thread_num()], nam_read2, record2.name,
//                                                   acc_map,
//                                                   map_param.k,
//                                                   record2.seq.length(), ref_lengths);
//                                joint_NAM_scores.clear();
//                            } else {
//                                align_PE(aln_params, output_streams[omp_get_thread_num()], nams1, nams2, record1,
//                                         record2,
//                                         acc_map, map_param.k,
//                                         ref_lengths, ref_seqs, log_vars,
//                                         map_param.dropoff_threshold, isize_est, map_param.maxTries, map_param.max_secondary);
//                            }
//                            auto extend_finish = std::chrono::high_resolution_clock::now();
//                            log_vars.tot_extend += extend_finish - extend_start;
//                            nams1.clear();
//                            nams2.clear();
//
//                            log_stats_vec[omp_get_thread_num()] = log_vars;
//                            isize_est_vec[omp_get_thread_num()] = isize_est;
//                        }
//                    }
//                }
//            }
//
//            //// DEBUG INFO, PRINT THE WORK PERFORMED BY EACH TREAD /////
//
//            for (int i = 0; i < n_threads; ++i) {
//                i_dist_est isize_est;
//                isize_est = isize_est_vec[i];
//                std::cerr << "THREAD " << i << " estimated diff in start coordinates b/t mates, (mean: " << isize_est.mu << ", stddev: " << isize_est.sigma << ") " << std::endl;
//            }
//
//            for (int i = 0; i < n_threads; ++i) {
//                logging_variables log_s;
//                log_s = log_stats_vec[i];
//                std::cerr << "THREAD " << i << ": " << log_s.tot_all_tried << " " << log_s.tot_ksw_aligned << " " << log_s.tot_rescued << " " <<  log_s.tot_construct_strobemers.count() << " " <<  log_s.tot_find_nams.count() << " " <<  log_s.tot_extend.count() << std::endl;
//            }
//
//            auto isize_est = isize_est_vec[omp_get_thread_num()];
//            std::cerr << "Estimated diff in start coordinates b/t mates, (mean: " << isize_est.mu << ", stddev: " << isize_est.sigma << ") " << std::endl;
//            records1 = new_records1;
//            records2 = new_records2;
//            // Reset isize estiamtions
//            for (int i = 0; i < n_threads; ++i) {
//                i_dist_est isize_est;
//                isize_est_vec[i] = isize_est;
//            }
//
//            auto write_start = std::chrono::high_resolution_clock::now();
//            for (int i = 0; i < n_threads; ++i) {
//                out << output_streams[i];
//                output_streams[i].clear();
//            }
//            auto write_finish = std::chrono::high_resolution_clock::now();
//            log_stats_vec[omp_get_thread_num()].tot_write_file += write_finish - write_start;
//        }

        gzclose(fp1);
        gzclose(fp2);

        for (int i = 0; i < n_threads; ++i) {
            //auto isize_est = isize_est_vec[i];
            auto log_vars = log_stats_vec[std::this_thread::get_id()];
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

    }


//    out.close();

    std::cerr << "Total mapping sites tried: " << tot_log_vars.tot_all_tried << std::endl;
    std::cerr << "Total calls to ssw: " << tot_log_vars.tot_ksw_aligned << std::endl;
    std::cerr << "Calls to ksw (rescue mode): " << tot_log_vars.tot_rescued << std::endl;
    std::cerr << "Did not fit strobe start site: " << tot_log_vars.did_not_fit  << std::endl;
    std::cerr << "Tried rescue: " << tot_log_vars.tried_rescue  << std::endl;
    // Record mapping end time
    auto finish_aln_part = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tot_aln_part = finish_aln_part - start_aln_part;
    std::cerr << "Total time mapping: " << tot_aln_part.count() << " s." <<  std::endl;
    std::cerr << "Total time reading read-file(s): " << tot_log_vars.tot_read_file.count() << " s." <<  std::endl;
    std::cerr << "Total time creating strobemers: " << tot_log_vars.tot_construct_strobemers.count()/n_threads << " s." <<  std::endl;
    std::cerr << "Total time finding NAMs (non-rescue mode): " << tot_log_vars.tot_find_nams.count()/n_threads  << " s." <<  std::endl;
    std::cerr << "Total time finding NAMs (rescue mode): " << tot_log_vars.tot_time_rescue.count()/n_threads  << " s." <<  std::endl;
//    std::cerr << "Total time finding NAMs ALTERNATIVE (candidate sites): " << tot_find_nams_alt.count()/n_threads  << " s." <<  std::endl;
    std::cerr << "Total time sorting NAMs (candidate sites): " << tot_log_vars.tot_sort_nams.count()/n_threads  << " s." <<  std::endl;
    std::cerr << "Total time reverse compl seq: " << tot_log_vars.tot_rc.count()/n_threads  << " s." <<  std::endl;
    std::cerr << "Total time base level alignment (ssw): " << tot_log_vars.tot_extend.count()/n_threads  << " s." <<  std::endl;
    std::cerr << "Total time writing alignment to files: " << tot_log_vars.tot_write_file.count() << " s." <<  std::endl;

    //////////////////////////////////////////////////////////////////////////


    /////////////////////// FIND AND OUTPUT NAMs ///////////////////////////////





}

