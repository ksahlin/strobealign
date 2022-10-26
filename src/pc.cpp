//
// Created by Kristoffer Sahlin on 3/22/22.
//

// Using initial base format of Buffer classed from: https://andrew128.github.io/ProducerConsumer/

#include "pc.hpp"
#include <thread>
#include <iostream>
#include <chrono>
#include <queue>

#include "robin_hood.h"
#include "index.hpp"
#include "kseq++.hpp"

using namespace klibpp;


void InputBuffer::read_records_PE(std::vector<KSeq> &records1, std::vector<KSeq> &records2, AlignmentStatistics &statistics) {
     auto read_start = std::chrono::high_resolution_clock::now();
    // Acquire a unique lock on the mutex
    std::unique_lock<std::mutex> unique_lock(mtx);
    records1 = ks1.read(chunk_size);
    records2 = ks2.read(chunk_size);

    if (records1.empty()){
        finished_reading = true;
    }

    // Unlock unique lock
    unique_lock.unlock();
    // Notify a single thread that buffer isn't empty
    not_empty.notify_one();
    auto read_finish = std::chrono::high_resolution_clock::now();
    statistics.tot_read_file += read_finish - read_start;

}

void InputBuffer::read_records_SE(std::vector<KSeq> &records1, AlignmentStatistics &statistics) {
    auto read_start = std::chrono::high_resolution_clock::now();
    // Acquire a unique lock on the mutex
    std::unique_lock<std::mutex> unique_lock(mtx);
    records1 = ks1.read(chunk_size);

    if (records1.empty()){
        finished_reading = true;
    }

    // Unlock unique lock
    unique_lock.unlock();
    // Notify a single thread that buffer isn't empty
    not_empty.notify_one();
    auto read_finish = std::chrono::high_resolution_clock::now();
    statistics.tot_read_file += read_finish - read_start;

}

void OutputBuffer::output_records(std::string &sam_alignments) {
    // Acquire a unique lock on the mutex
    std::unique_lock<std::mutex> unique_lock(mtx);
    out << sam_alignments;
    // Update appropriate fields
    buffer_size = 0;
    // Unlock unique lock
    unique_lock.unlock();
    not_full.notify_one();
}


inline bool align_reads_PE(
    InputBuffer &input_buffer,
    OutputBuffer &output_buffer,
    std::vector<KSeq> &records1,
    std::vector<KSeq> &records2,
    AlignmentStatistics &statistics,
    i_dist_est &isize_est,
    alignment_params &aln_params,
    mapping_params &map_param,
    const IndexParameters& index_parameters,
    const References& references,
    const StrobemerIndex& index
) {

    // If no more reads to align
    if (records1.empty() && input_buffer.finished_reading){
        return true;
    }

    std::string sam_out;
    sam_out.reserve(7*map_param.r *records1.size());
    for (size_t i = 0; i < records1.size(); ++i) {
        auto record1 = records1[i];
        auto record2 = records2[i];

        align_PE_read(record1, record2, sam_out, statistics, isize_est, aln_params,
                      map_param, index_parameters, references, index);
    }
//    std::cerr << isize_est_vec[thread_id].mu << " " << isize_est_vec[thread_id].sigma << "\n";
//    std::cerr << log_stats_vec[thread_id].tot_all_tried << " " << log_stats_vec[thread_id].tot_ksw_aligned << "\n";
//    IMMEDIATELY PRINT TO STDOUT/FILE HERE
    output_buffer.output_records(sam_out);
    return false;
}


void perform_task_PE(
    InputBuffer &input_buffer,
    OutputBuffer &output_buffer,
    std::unordered_map<std::thread::id, AlignmentStatistics> &log_stats_vec,
    std::unordered_map<std::thread::id, i_dist_est> &isize_est_vec,
    alignment_params &aln_params,
    mapping_params &map_param,
    const IndexParameters& index_parameters,
    const References& references,
    const StrobemerIndex& index
) {
    bool eof = false;
    while (true){
        std::vector<KSeq> records1;
        std::vector<KSeq> records2;
        auto thread_id = std::this_thread::get_id();
        if (log_stats_vec.find(thread_id) == log_stats_vec.end()) { //  Not initialized
            AlignmentStatistics statistics;
            log_stats_vec[thread_id] = statistics;
        }
        input_buffer.read_records_PE(records1, records2, log_stats_vec[thread_id]);
        eof = align_reads_PE(input_buffer, output_buffer, records1, records2,
                           log_stats_vec[thread_id], isize_est_vec[thread_id],
                          aln_params, map_param, index_parameters, references, index);

        if (eof) {
            break;
        }
    }
}


inline bool align_reads_SE(
    InputBuffer &input_buffer,
    OutputBuffer &output_buffer,
    std::vector<KSeq> &records,
    AlignmentStatistics &statistics,
    alignment_params &aln_params,
    mapping_params &map_param,
    const IndexParameters& index_parameters,
    const References& references,
    const StrobemerIndex& index
) {
    // If no more reads to align
    if (records.empty() && input_buffer.finished_reading){
        return true;
    }

    std::string sam_out;
    sam_out.reserve(7*map_param.r *records.size());
    for (size_t i = 0; i < records.size(); ++i) {
        auto record1 = records[i];

        align_SE_read(record1, sam_out,  statistics, aln_params,
                      map_param, index_parameters, references, index);
    }
//    std::cerr << isize_est_vec[thread_id].mu << " " << isize_est_vec[thread_id].sigma << "\n";
//    std::cerr << log_stats_vec[thread_id].tot_all_tried << " " << log_stats_vec[thread_id].tot_ksw_aligned << "\n";
//    IMMEDIATELY PRINT TO STDOUT/FILE HERE
    output_buffer.output_records(sam_out);
    return false;
}


void perform_task_SE(
    InputBuffer &input_buffer,
    OutputBuffer &output_buffer,
    std::unordered_map<std::thread::id,
    AlignmentStatistics> &log_stats_vec,
    alignment_params &aln_params,
    mapping_params &map_param,
    const IndexParameters& index_parameters,
    const References& references,
    const StrobemerIndex& index
) {
    bool eof = false;
    while (true){
        std::vector<KSeq> records1;
        auto thread_id = std::this_thread::get_id();
        if (log_stats_vec.find(thread_id) == log_stats_vec.end()) { //  Not initialized
            AlignmentStatistics statistics;
            log_stats_vec[thread_id] = statistics;
        }
        input_buffer.read_records_SE(records1, log_stats_vec[thread_id]);
        eof = align_reads_SE(input_buffer, output_buffer, records1,
                             log_stats_vec[thread_id],
                             aln_params, map_param, index_parameters, references, index);

        if (eof){
            break;
        }
    }

}
