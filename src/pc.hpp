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
#include <unordered_map>

#include "robin_hood.h"
#include <zlib.h>
#include "kseq++.hpp"
using namespace klibpp;
#include "ssw_cpp.h"
#include "index.hpp"
#include "aln.hpp"


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
    int X = 100000; // input read chunk size

    void read_records_PE(std::vector<KSeq> &records1, std::vector<KSeq> &records2, logging_variables &log_vars);
    void read_records_SE(std::vector<KSeq> &records1, logging_variables &log_vars);

};


class OutputBuffer {
    // OutputBuffer fields


public:

    // Fields for concurrency input
    std::mutex mtx;
    std::condition_variable not_empty;
    std::condition_variable not_full;

    int buffer_size = 0;
    std::ostream &out;

    void output_records(std::string &sam_alignments);

};



void perform_task_PE(InputBuffer &input_buffer, OutputBuffer &output_buffer,
                  std::unordered_map<std::thread::id, logging_variables> &log_stats_vec, std::unordered_map<std::thread::id, i_dist_est> &isize_est_vec, alignment_params &aln_params,
                  mapping_params &map_param, std::vector<unsigned int> &ref_lengths, std::vector<std::string> &ref_seqs, kmer_lookup &mers_index, mers_vector &flat_vector, idx_to_acc &acc_map );

void perform_task_SE(InputBuffer &input_buffer, OutputBuffer &output_buffer,
                     std::unordered_map<std::thread::id, logging_variables> &log_stats_vec, alignment_params &aln_params,
                     mapping_params &map_param, std::vector<unsigned int> &ref_lengths, std::vector<std::string> &ref_seqs, kmer_lookup &mers_index, mers_vector &flat_vector, idx_to_acc &acc_map );
#endif // pc_hpp_
