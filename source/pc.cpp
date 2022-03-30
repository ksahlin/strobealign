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


void InputBuffer::read_records(std::thread::id  thread_id, std::vector<KSeq> &records1, std::vector<KSeq> &records2) {
    // Acquire a unique lock on the mutex
    std::unique_lock<std::mutex> unique_lock(mtx);
    records1 = ks1.read(X);
    records2 = ks2.read(X);

    if (records1.empty()){
        finished_reading = true;
    }

    // Unlock unique lock
    unique_lock.unlock();
    // Notify a single thread that buffer isn't empty
    not_empty.notify_one();

}

void InputBuffer::add_records(std::thread::id  thread_id) {
    // Acquire a unique lock on the mutex
    std::unique_lock<std::mutex> unique_lock(mtx);

    // Wait if the buffer is full
    not_full.wait(unique_lock, [this]() {
        return buffer_size != INPUT_BUFFER_CAPACITY;
    });

    // Add input to buffer

    auto records1 = ks1.read(X);
    q1.push(records1);
    auto records2 = ks2.read(X);
    q2.push(records2);

    // Update appropriate fields
    buffer_size ++;
//    buffer_size += records.size(); //++;

    if (records1.empty()){
        finished_reading = true;
    }

    // Unlock unique lock
    unique_lock.unlock();
    // Notify a single thread that buffer isn't empty
    not_empty.notify_one();

}


//void OutputBuffer::add_aligned_reads(std::thread::id  thread_id, std::string &sam_alignments) {
////    output_strings[thread_id].append(sam_alignment);
////    output_strings[thread_id].append("\n");
//
//    // Acquire a unique lock on the mutex
//    std::unique_lock<std::mutex> unique_lock(mtx);
//
//    // Wait if the buffer is full
//    not_full.wait(unique_lock, [this]() {
//        return buffer_size < OUTPUT_BUFFER_CAPACITY;
//    });
//
//    // Add input to buffer
//    out.append(sam_alignments);
////    out.append("\n");
//    // Update appropriate fields
//    buffer_size++;
//    // Unlock unique lock
//    unique_lock.unlock();
//    // Notify a single thread that buffer isn't empty
//    not_empty.notify_one();
//
//}

void OutputBuffer::output_records(std::thread::id  thread_id, std::string &sam_alignments) {
    // Acquire a unique lock on the mutex
    std::unique_lock<std::mutex> unique_lock(mtx);

//    // Wait if buffer is empty
//    not_empty.wait(unique_lock, [this]() {
//        return buffer_size > 0;
//    });

//    TODO: WRITE TO EITHER FILE OR STDOUT HERE - NEED TO PASS out AS ARG
    std::cout << sam_alignments;
//    sam_alignments.clear();

//    int q_size = q.size();
    // Update appropriate fields
    buffer_size = 0;
    // Unlock unique lock
    unique_lock.unlock();
    not_full.notify_one();

}



//inline bool align_reads(std::thread::id thread_id, InputBuffer &input_buffer, OutputBuffer &output_buffer,
//                               robin_hood::unordered_map<std::thread::id, logging_variables> &log_stats_vec,
//                               robin_hood::unordered_map<std::thread::id, i_dist_est> &isize_est_vec, alignment_params &aln_params,
//                 mapping_params &map_param, std::vector<unsigned int> &ref_lengths, std::vector<std::string> &ref_seqs,
//                 kmer_lookup &mers_index, mers_vector &flat_vector, idx_to_acc &acc_map ) {
//
//    // Acquire a unique lock on the mutex
//    std::unique_lock<std::mutex> unique_input_lock(input_buffer.mtx);
//
//    // If no more reads to align
//    if (input_buffer.q1.empty() && input_buffer.finished_reading){
//        unique_input_lock.unlock();
//        return true;
//    }
//
//    // More reads to align but buffer is empty (wait)
//    input_buffer.not_empty.wait(unique_input_lock, [&input_buffer]() {
//        return input_buffer.buffer_size != 0;
//    });
//
//    // Get value from position to remove in buffer
//    auto records1 = input_buffer.q1.front();
//    input_buffer.q1.pop();
//    auto records2 = input_buffer.q2.front();
//    input_buffer.q2.pop();
//
//    // Update appropriate fields
//    input_buffer.buffer_size--;
//    // Unlock unique lock
//    unique_input_lock.unlock();
//
//    // Notify a single thread that the buffer isn't full
//    input_buffer.not_full.notify_one();
//
//
////    std::string ref = "TGCCGAGTGATATCGCTGACGTCATCCTTGAGGGTGAAGTTCAGGTCGTCGAGCAACTCGGCAACGAAACTCAAATCCATATCCAGATCCCTTCCATTCGTGACGGATCCAGACGGTACCGGAGACATTCGGTTGCTGGATAGCTGGTTGTTGTGTTGCTGAAATAGACGTATTTGCAGCCGGTGCTGGAGTCTGAATATACCTGGTTTTTGATTCTGGCAGGACTGGGGTTACGTAGCATTATTGCTAACCCGGAAGTGCTGCATGCGCTGAATCCAATGTGGGCGGTGCATTTCTTCCTCCAATTAGTGAATGAACAGGGAACCGTGCTTGTTCAGGATCTGGCGGGAGTATTTGCTGCCTCGGAAGCGACAATCCGTGCCGATTTGCGCTTTCTCGAATTTGCCGCTGCATTGATCGTCAGTGGCCTGCTCGTTGGCTGTAATCAACTCACCCAATACACCATCACCGAACAAGAAATTAACCAGTCGCTTGCGAAA";
////    std::string alignments;
////    alignments.reserve(400*records1.size());
////    for (size_t i = 0; i < records1.size(); ++i) {
////        auto record1 = records1[i];
////        auto record2 = records2[i];
////        auto aln = ssw_align(ref, record1.seq, record1.seq.size(), 2, 8, 12, 1);
////        auto aln2 = ssw_align(ref, record2.seq, record2.seq.size(), 2, 8, 12, 1);
////        alignments.append(record1.name);
////        alignments.append("\n");
////        alignments.append(record2.name);
////        alignments.append("\n");
////    }
////    output_buffer.add_aligned_reads(thread_id, alignments);
//
//    std::string sam_out;
//    sam_out.reserve(7*map_param.r *records1.size());
//    for (size_t i = 0; i < records1.size(); ++i) {
//        auto record1 = records1[i];
//        auto record2 = records2[i];
//        auto log_vars = log_stats_vec[thread_id];
//        auto isize_est = isize_est_vec[thread_id];
//
//        align_PE_read(thread_id, record1, record2, sam_out,  log_vars, isize_est, aln_params,
//                map_param, ref_lengths, ref_seqs,
//                mers_index, flat_vector, acc_map );
////        align_PE_read(thread_id, record1, record2, sam_out, log_vars, isize_est, aln_params,
////                map_param, ref_lengths, ref_seqs,
////                mers_index, flat_vector, acc_map );
//        log_stats_vec[thread_id] = log_vars;
//        isize_est_vec[thread_id] = isize_est;
//    }
////    IMMEDIATELY PRINT TO STDOUT/FILE HERE
//    output_buffer.output_records(thread_id, sam_out);
////    output_buffer.add_aligned_reads(thread_id, sam_out);
//
//    return false;
//}


inline bool align_reads2(std::thread::id thread_id, InputBuffer &input_buffer, OutputBuffer &output_buffer,  std::vector<KSeq> &records1,  std::vector<KSeq> &records2,
                        robin_hood::unordered_map<std::thread::id, logging_variables> &log_stats_vec,
                        robin_hood::unordered_map<std::thread::id, i_dist_est> &isize_est_vec, alignment_params &aln_params,
                        mapping_params &map_param, std::vector<unsigned int> &ref_lengths, std::vector<std::string> &ref_seqs,
                        kmer_lookup &mers_index, mers_vector &flat_vector, idx_to_acc &acc_map ) {

    // If no more reads to align
    if (records1.empty() && input_buffer.finished_reading){
        return true;
    }

    std::string sam_out;
    sam_out.reserve(7*map_param.r *records1.size());
    for (size_t i = 0; i < records1.size(); ++i) {
        auto record1 = records1[i];
        auto record2 = records2[i];
        auto log_vars = log_stats_vec[thread_id];
        auto isize_est = isize_est_vec[thread_id];

        align_PE_read(thread_id, record1, record2, sam_out,  log_vars, isize_est, aln_params,
                      map_param, ref_lengths, ref_seqs,
                      mers_index, flat_vector, acc_map );
//        align_PE_read(thread_id, record1, record2, sam_out, log_vars, isize_est, aln_params,
//                map_param, ref_lengths, ref_seqs,
//                mers_index, flat_vector, acc_map );
        log_stats_vec[thread_id] = log_vars;
        isize_est_vec[thread_id] = isize_est;
    }
//    IMMEDIATELY PRINT TO STDOUT/FILE HERE
    output_buffer.output_records(thread_id, sam_out);
//    output_buffer.add_aligned_reads(thread_id, sam_out);

    return false;
}



void perform_task(InputBuffer &input_buffer, OutputBuffer &output_buffer,
                  robin_hood::unordered_map<std::thread::id, logging_variables> &log_stats_vec, robin_hood::unordered_map<std::thread::id, i_dist_est> &isize_est_vec, alignment_params &aln_params,
                  mapping_params &map_param, std::vector<unsigned int> &ref_lengths, std::vector<std::string> &ref_seqs,
                  kmer_lookup &mers_index, mers_vector &flat_vector, idx_to_acc &acc_map ){
    bool eof = false;
    while (true){
//        if (output_buffer.buffer_size >= OUTPUT_BUFFER_CAPACITY ){ // first try write if buffer is full
//            output_buffer.output_records(std::this_thread::get_id()); // Implement write here
//        } else
//       if ( (!input_buffer.finished_reading) && (input_buffer.buffer_size < INPUT_BUFFER_CAPACITY/2) ){
//            input_buffer.add_records(std::this_thread::get_id());
//        }
//       else { // otherwise align
//            eof = align_reads(std::this_thread::get_id(), input_buffer, output_buffer, log_stats_vec, isize_est_vec,
//                              aln_params, map_param, ref_lengths, ref_seqs, mers_index, flat_vector,  acc_map);
//        }

        std::vector<KSeq> records1;
        std::vector<KSeq> records2;
        input_buffer.read_records(std::this_thread::get_id(), records1, records2);
        eof = align_reads2(std::this_thread::get_id(), input_buffer, output_buffer, records1, records2,
                           log_stats_vec, isize_est_vec,
                          aln_params, map_param, ref_lengths, ref_seqs, mers_index, flat_vector,  acc_map);

        if (eof){
            break;
        }
//        thread_states[std::this_thread::get_id()] ++;
    }

}



//
//int main (int argc, char **argv) {
////    std::string ref = "ACAGTGTTCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
////    std::string ref_segm = ref.substr(0, 0);
////    std::cout << ref_segm << "\n";
//
//    int n_threads = 2;
//    int input_chunk_size = 10;
//    std::cerr << "Nr threads: " << n_threads << "\n";
//    // Create Buffers
//    const char *fname = "/Users/kxs624/tmp/STROBEALIGN/test_MT.fq";
//    gzFile fp = gzopen(fname, "r");
//    auto ks = make_ikstream(fp, gzread);
//    InputBuffer input_buffer = { {}, {}, {}, {}, ks, false, 0, input_chunk_size};
////    input_buffer.init_file(fname)
//    OutputBuffer output_buffer;
//    output_buffer.out.reserve((OUTPUT_BUFFER_CAPACITY) * 600);
//
//    std::unordered_map<std::thread::id, int> thread_states(n_threads);
//
//    std::vector<std::thread> workers;
//    for (int i = 0; i < n_threads; ++i) {
//        std::thread consumer(perform_task, std::ref(input_buffer), std::ref(output_buffer), std::ref(thread_states));
//        workers.push_back(std::move(consumer));
//    }
//
//    for (size_t i = 0; i < workers.size(); ++i) {
//        workers[i].join();
//    }
//
//    std::cerr << "LAST BUFFER SIZE: " << output_buffer.buffer_size <<  std::endl;
//    if (output_buffer.buffer_size > 0 ){ // write last set of records
//        std::cerr << "Final writing to output " <<output_buffer.buffer_size << " " << input_buffer.buffer_size  << " threadID: " << std::this_thread::get_id() << " " << input_buffer.finished_reading << std::endl;
//        output_buffer.output_records(std::this_thread::get_id()); // Implement write here
//    }
//
//    std::cerr << "Done!\n";
//    return 0;
//}
//
