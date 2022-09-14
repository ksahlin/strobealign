#ifndef pc_hpp
#define pc_hpp

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
#include <zlib.h>

#include "robin_hood.h"
#include "kseq++.hpp"
#include "ssw_cpp.h"
#include "index.hpp"
#include "aln.hpp"
#include "refs.hpp"

class InputBuffer {

public:
    typedef klibpp::KStream<gzFile_s*, int (*)(gzFile_s*, void*, unsigned int), klibpp::mode::In_> input_stream_t;

    InputBuffer(input_stream_t& ks1, input_stream_t& ks2, int chunk_size)
    : ks1(ks1), ks2(ks2), chunk_size(chunk_size) { }

    // Fields for concurrency input
    std::mutex mtx;
    std::condition_variable not_empty;
    std::condition_variable not_full;

    std::queue<std::vector<klibpp::KSeq>> q1;
    std::queue<std::vector<klibpp::KSeq>> q2;
    input_stream_t &ks1;
    input_stream_t &ks2;
    bool finished_reading = false;
    int buffer_size = 0;
    int chunk_size = 100000;

    void read_records_PE(std::vector<klibpp::KSeq> &records1, std::vector<klibpp::KSeq> &records2, logging_variables &log_vars);
    void read_records_SE(std::vector<klibpp::KSeq> &records1, logging_variables &log_vars);
};


class OutputBuffer {

public:
    OutputBuffer(std::ostream& out) : out(out) { }

    std::mutex mtx;
    std::condition_variable not_empty;
    std::condition_variable not_full;

    int buffer_size = 0;
    std::ostream &out;

    void output_records(std::string &sam_alignments);
};



void perform_task_PE(InputBuffer &input_buffer, OutputBuffer &output_buffer,
                  std::unordered_map<std::thread::id, logging_variables> &log_stats_vec, std::unordered_map<std::thread::id, i_dist_est> &isize_est_vec, alignment_params &aln_params,
                  mapping_params &map_param, std::vector<unsigned int> &ref_lengths, std::vector<std::string> &ref_seqs, kmer_lookup &mers_index, mers_vector &flat_vector, ref_names &acc_map );

void perform_task_SE(InputBuffer &input_buffer, OutputBuffer &output_buffer,
                     std::unordered_map<std::thread::id, logging_variables> &log_stats_vec, alignment_params &aln_params,
                     mapping_params &map_param, std::vector<unsigned int> &ref_lengths, std::vector<std::string> &ref_seqs, kmer_lookup &mers_index, mers_vector &flat_vector, ref_names &acc_map );
#endif // pc_hpp_
