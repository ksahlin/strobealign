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

#include "robin_hood.h"
#include "ssw_cpp.h"
#include "index.hpp"
#include "aln.hpp"
#include "refs.hpp"
#include "fastq.hpp"

class InputBuffer {

public:

    InputBuffer(input_stream_t& ks1, input_stream_t& ks2, int chunk_size)
    : ks1(ks1), ks2(ks2), chunk_size(chunk_size) { }

    std::mutex mtx;

    input_stream_t &ks1;
    input_stream_t &ks2;
    bool finished_reading{false};
    int chunk_size;
    size_t chunk_index{0};

    size_t read_records_PE(std::vector<klibpp::KSeq> &records1, std::vector<klibpp::KSeq> &records2, AlignmentStatistics &statistics);
    size_t read_records_SE(std::vector<klibpp::KSeq> &records1, AlignmentStatistics &statistics);
};


class OutputBuffer {

public:
    OutputBuffer(std::ostream& out) : out(out) { }

    std::mutex mtx;
    std::ostream &out;
    std::unordered_map<size_t, std::string> chunks;
    size_t next_chunk_index{0};

    void output_records(std::string chunk, size_t chunk_index);
};


void perform_task_PE(InputBuffer &input_buffer, OutputBuffer &output_buffer,
                  AlignmentStatistics& statistics, const alignment_params &aln_params,
                  const mapping_params &map_param, const IndexParameters& index_parameters, const References& references, const StrobemerIndex& index, const std::string& read_group_id);

void perform_task_SE(InputBuffer &input_buffer, OutputBuffer &output_buffer,
                     AlignmentStatistics& statistics, const alignment_params &aln_params,
                     const mapping_params &map_param, const IndexParameters& index_parameters, const References& references, const StrobemerIndex& index, const std::string& read_group_id);

#endif
