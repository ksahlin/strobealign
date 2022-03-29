#ifndef pc_hpp
#define pc_hpp

// The declarations in the pc.hpp file
#include <thread>
#include <condition_variable>
#include <mutex>
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <queue>
#include <vector>
#include <sstream>

#include "robin_hood.h"
#include <zlib.h>
#include "kseq++.hpp"
using namespace klibpp;
#include "ssw_cpp.h"
#include "index.hpp"
#include "aln.hpp"

#define INPUT_BUFFER_CAPACITY 100

#define OUTPUT_BUFFER_CAPACITY 200

//int dummy_align_read(int q_size);


class InputBuffer {
    // InputBuffer fields


public:
    // Fields for concurrency input
    std::mutex mtx;
    std::condition_variable not_empty;
    std::condition_variable not_full;

    std::queue<std::vector<KSeq>> q1;
    std::queue<std::vector<KSeq>> q2;
    klibpp::KStream<gzFile_s*, int (*)(gzFile_s*, void*, unsigned int), klibpp::mode::In_> &ks1;
    klibpp::KStream<gzFile_s*, int (*)(gzFile_s*, void*, unsigned int), klibpp::mode::In_> &ks2;
    bool finished_reading = false;
    int buffer_size = 0;
    int X = 5000; // input read chunk size

    void add_records(std::thread::id thread_id);
//    void init_file(std::thread::id thread_id, const char * fname);
    // bool align_read(std::thread::id  thread_id);

};


class OutputBuffer {
    // OutputBuffer fields


public:

    // Fields for concurrency input
    std::mutex mtx;
    std::condition_variable not_empty;
    std::condition_variable not_full;

    int buffer_size = 0;
    std::string out;
//    int n_threads;
//    std::vector<std::string> output_strings;
    OutputBuffer( int64_t reserve_size){ // Constructor
        out.reserve(reserve_size);
    }
    //OutputBuffer( int n) {     // Constructor
//        int n_threads = n;
//        output_strings.reserve(n_threads);
//        for (int i = 0; i < n_threads; ++i) {
//            output_strings[i].reserve(OUTPUT_BUFFER_CAPACITY *
//                                      600); // Reserve sufficient space for appending multiple SAM records (400 is an upper setimate on the number of characters for each sam record of a 200-300bp read)
//        }
//    }


    void add_aligned_reads(std::thread::id  thread_id,  std::string &sam_alignments);

    void output_records(std::thread::id thread_id);

};



void perform_task(InputBuffer &input_buffer, OutputBuffer &output_buffer,
                  robin_hood::unordered_map<std::thread::id, logging_variables> &log_stats_vec, robin_hood::unordered_map<std::thread::id, i_dist_est> &isize_est_vec, alignment_params &aln_params,
                  mapping_params &map_param, std::vector<unsigned int> &ref_lengths, std::vector<std::string> &ref_seqs, kmer_lookup &mers_index, mers_vector &flat_vector, idx_to_acc &acc_map );

#endif // pc_hpp_