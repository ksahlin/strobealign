//
// Created by Kristoffer Sahlin on 3/22/22.
//

// Using initial base format of Buffer classed from: https://andrew128.github.io/ProducerConsumer/

#include "pc.hpp"
#include <mutex>
#include <iostream>
#include <chrono>
#include <queue>

#include "timer.hpp"
#include "robin_hood.h"
#include "index.hpp"
#include "kseq++.hpp"
#include "sam.hpp"

/* Strip the /1 or /2 suffix from a read name */
void strip_suffix(std::string& name) {
    auto len = name.length();
    if (
        len >= 2
        && name[len - 2] == '/'
        && (name[len - 1] == '1' || name[len - 1] == '2')
    ) {
        name.pop_back();
        name.pop_back();
        // C++20 would allow this:
        // name.resize(len - 2);
    }
}

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

void OutputBuffer::output_records(std::string chunk, size_t chunk_index) {
    std::unique_lock<std::mutex> unique_lock(mtx);

    // Ensure we print the chunks in the order in which they were read
    assert(chunks.count(chunk_index) == 0);
    chunks.emplace(std::make_pair(chunk_index, chunk));
    while (true) {
        const auto& item = chunks.find(next_chunk_index);
        if (item == chunks.end()) {
            break;
        }
        out << item->second;
        chunks.erase(item);
        next_chunk_index++;
    }
    unique_lock.unlock();
}

void perform_task_PE(
    InputBuffer &input_buffer,
    OutputBuffer &output_buffer,
    AlignmentStatistics& statistics,
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
        auto chunk_index = input_buffer.read_records_PE(records1, records2, statistics);
        i_dist_est isize_est;
        if (records1.empty() && input_buffer.finished_reading){
            break;
        }

        std::string sam_out;
        sam_out.reserve(7*map_param.r *records1.size());
        Sam sam{sam_out, references, read_group_id, map_param.output_unmapped};
        for (size_t i = 0; i < records1.size(); ++i) {
            auto record1 = records1[i];
            auto record2 = records2[i];
            to_uppercase(record1.seq);
            to_uppercase(record2.seq);
            strip_suffix(record1.name);
            strip_suffix(record2.name);
            align_PE_read(record1, record2, sam, sam_out, statistics, isize_est, aln_params,
                        map_param, index_parameters, references, index);
        }
        output_buffer.output_records(std::move(sam_out), chunk_index);
        assert(sam_out == "");
    }
}


void perform_task_SE(
    InputBuffer &input_buffer,
    OutputBuffer &output_buffer,
    AlignmentStatistics& statistics,
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
        auto chunk_index = input_buffer.read_records_SE(records, statistics);

        if (records.empty() && input_buffer.finished_reading){
            break;
        }

        std::string sam_out;
        sam_out.reserve(7*map_param.r *records.size());
        Sam sam{sam_out, references, read_group_id, map_param.output_unmapped};
        for (size_t i = 0; i < records.size(); ++i) {
            auto record = records[i];
            to_uppercase(record.seq);
            strip_suffix(record.name);
            align_SE_read(record, sam, sam_out, statistics, aln_params, map_param, index_parameters, references, index);
        }
        output_buffer.output_records(std::move(sam_out), chunk_index);
        assert(sam_out == "");
    }
}
