//
// Created by Kristoffer Sahlin on 3/22/22.
//

// Using initial base format of Buffer classed from: https://andrew128.github.io/ProducerConsumer/

#include "pc.hpp"
#include <thread>
#include <iostream>
#include <chrono>
#include <queue>

#include "timer.hpp"
#include "robin_hood.h"
#include "index.hpp"
#include "kseq++.hpp"
#include "sam.hpp"


size_t InputBuffer::read_records_PE(std::vector<klibpp::KSeq> &records1, std::vector<klibpp::KSeq> &records2, AlignmentStatistics &statistics) {
    Timer timer;
    // Acquire a unique lock on the mutex
    std::unique_lock<std::mutex> unique_lock(mtx);
    records1 = ks1.read(chunk_size);
    records2 = ks2.read(chunk_size);
    size_t current_chunk_index = chunk_index;
    chunk_index++;

    if (records1.empty()) {
        finished_reading = true;
    }

    unique_lock.unlock();
    statistics.tot_read_file += timer.duration();

    return current_chunk_index;
}

size_t InputBuffer::read_records_SE(std::vector<klibpp::KSeq> &records1, AlignmentStatistics &statistics) {
    Timer timer;
    // Acquire a unique lock on the mutex
    std::unique_lock<std::mutex> unique_lock(mtx);
    records1 = ks1.read(chunk_size);
    size_t current_chunk_index = chunk_index;
    chunk_index++;

    if (records1.empty()){
        finished_reading = true;
    }

    unique_lock.unlock();
    statistics.tot_read_file += timer.duration();

    return current_chunk_index;
}

void OutputBuffer::output_records(std::string &sam_alignments) {
    std::unique_lock<std::mutex> unique_lock(mtx);
    out << sam_alignments;
    unique_lock.unlock();
}

void perform_task_PE(
    InputBuffer &input_buffer,
    OutputBuffer &output_buffer,
    std::unordered_map<std::thread::id, AlignmentStatistics> &log_stats_vec,
    std::unordered_map<std::thread::id, i_dist_est> &isize_est_vec,
    const alignment_params &aln_params,
    const mapping_params &map_param,
    const IndexParameters& index_parameters,
    const References& references,
    const StrobemerIndex& index,
    const std::string& read_group_id
) {
    bool eof = false;
    while (!eof) {
        std::vector<klibpp::KSeq> records1;
        std::vector<klibpp::KSeq> records2;
        auto thread_id = std::this_thread::get_id();
        if (log_stats_vec.find(thread_id) == log_stats_vec.end()) { //  Not initialized
            AlignmentStatistics statistics;
            log_stats_vec[thread_id] = statistics;
        }
        input_buffer.read_records_PE(records1, records2, log_stats_vec[thread_id]);
        if (records1.empty() && input_buffer.finished_reading){
            break;
        }

        std::string sam_out;
        sam_out.reserve(7*map_param.r *records1.size());
        Sam sam{sam_out, references, read_group_id};
        auto& statistics{log_stats_vec[thread_id]};
        auto& isize_est{isize_est_vec[thread_id]};
        for (size_t i = 0; i < records1.size(); ++i) {
            auto record1 = records1[i];
            auto record2 = records2[i];

            align_PE_read(record1, record2, sam, sam_out, statistics, isize_est, aln_params,
                        map_param, index_parameters, references, index);
        }
        // std::cerr << isize_est_vec[thread_id].mu << " " << isize_est_vec[thread_id].sigma << "\n";
        // std::cerr << log_stats_vec[thread_id].tot_all_tried << " " << log_stats_vec[thread_id].tot_ksw_aligned << "\n";
        // IMMEDIATELY PRINT TO STDOUT/FILE HERE
        output_buffer.output_records(sam_out);
    }
}


void perform_task_SE(
    InputBuffer &input_buffer,
    OutputBuffer &output_buffer,
    std::unordered_map<std::thread::id, AlignmentStatistics> &log_stats_vec,
    const alignment_params &aln_params,
    const mapping_params &map_param,
    const IndexParameters& index_parameters,
    const References& references,
    const StrobemerIndex& index,
    const std::string& read_group_id
) {
    bool eof = false;
    while (!eof){
        std::vector<klibpp::KSeq> records;
        auto thread_id = std::this_thread::get_id();
        if (log_stats_vec.find(thread_id) == log_stats_vec.end()) { //  Not initialized
            AlignmentStatistics statistics;
            log_stats_vec[thread_id] = statistics;
        }
        input_buffer.read_records_SE(records, log_stats_vec[thread_id]);

        if (records.empty() && input_buffer.finished_reading){
            break;
        }

        std::string sam_out;
        sam_out.reserve(7*map_param.r *records.size());
        Sam sam{sam_out, references, read_group_id};
        auto& statistics{log_stats_vec[thread_id]};
        for (size_t i = 0; i < records.size(); ++i) {
            auto record = records[i];

            align_SE_read(record, sam, sam_out, statistics, aln_params, map_param, index_parameters, references, index);
        }
    //    std::cerr << isize_est_vec[thread_id].mu << " " << isize_est_vec[thread_id].sigma << "\n";
    //    std::cerr << log_stats_vec[thread_id].tot_all_tried << " " << log_stats_vec[thread_id].tot_ksw_aligned << "\n";
    //    IMMEDIATELY PRINT TO STDOUT/FILE HERE
        output_buffer.output_records(sam_out);
    }
}
