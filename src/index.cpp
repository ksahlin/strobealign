//
//  index.cpp
//  cpp_course
//
//  Created by Kristoffer Sahlin on 4/21/21.
//
#include "index.hpp"

#include <math.h>   /* pow */
#include <fstream>
#include <cassert>
#include <algorithm>
#include "pdqsort/pdqsort.h"
#include <iostream>
#include <thread>
#include <atomic>
#include <hyperloglog/hyperloglog.hpp>
#include "io.hpp"
#include "timer.hpp"
#include "logger.hpp"

static Logger& logger = Logger::get();
static const uint32_t STI_FILE_FORMAT_VERSION = 2;

const int MAX_LINEAR_SEARCH = 4;
unsigned int StrobemerIndex::find(uint64_t key) const{
    const unsigned int top_N = key >> (64 - N);
    int position_start = hash_positions[top_N];
    int position_end = hash_positions[top_N + 1];

    if(position_start == position_end){
        return -1;
    }

    if (position_end - position_start < MAX_LINEAR_SEARCH) {
          for ( ; position_start < position_end; ++position_start) {
              if ((randstrobes_vector[position_start].hash & hash_mask) == (key & hash_mask)) return position_start;
          }
          return -1;
      }

    auto cmp = [this](const RefRandstrobeWithHash lhs, const RefRandstrobeWithHash rhs) {return (lhs.hash & hash_mask) < (rhs.hash & hash_mask); };

    auto pos = std::lower_bound(randstrobes_vector.begin() + position_start,
                                               randstrobes_vector.begin() + position_end,
                                               RefRandstrobeWithHash{key & hash_mask, 0, 0},
                                               cmp);
    if ((pos->hash & hash_mask) == (key & hash_mask)) return pos - randstrobes_vector.begin();
    return -1;
}

uint64_t count_unique_hashes(const std::vector<RefRandstrobeWithHash>& mers){
    if (mers.empty()) {
        return 0;
    }
    uint64_t prev_k = mers.at(0).hash;
    uint64_t unique_elements = 1;
    for (auto &curr_k : mers) {
        if (curr_k.hash != prev_k) {
            unique_elements ++;
        }
        prev_k = curr_k.hash;
    }
    return unique_elements;
}

void StrobemerIndex::write(const std::string& filename) const {
    std::ofstream ofs(filename, std::ios::binary);

    ofs.write("STI\1", 4); // Magic number
    write_int_to_ostream(ofs, STI_FILE_FORMAT_VERSION);

    // Variable-length chunk reserved for future use
    std::vector<char> reserved_chunk{0, 0, 0, 0, 0, 0, 0, 0};
    write_vector(ofs, reserved_chunk);

    write_int_to_ostream(ofs, filter_cutoff);
    parameters.write(ofs);

    write_vector(ofs, randstrobes_vector);
    write_vector(ofs, hash_positions);
}

void StrobemerIndex::read(const std::string& filename) {
    errno = 0;
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs.is_open()) {
        throw InvalidIndexFile(filename + ": " + strerror(errno));
    }

    union {
        char s[4];
        uint32_t v;
    } magic;
    ifs.read(magic.s, 4);
    if (magic.v != 0x01495453) { // "STI\1"
        throw InvalidIndexFile("Index file has incorrect format (magic number mismatch)");
    }

    uint32_t file_format_version = read_int_from_istream(ifs);
    if (file_format_version != STI_FILE_FORMAT_VERSION) {
        std::stringstream s;
        s << "Can only read index file format version " << STI_FILE_FORMAT_VERSION
            << ", but found version " << file_format_version;
        throw InvalidIndexFile(s.str());
    }

    // Skip over variable-length chunk reserved for future use
    uint64_t reserved_chunk_size;
    ifs.read(reinterpret_cast<char*>(&reserved_chunk_size), sizeof(reserved_chunk_size));
    ifs.seekg(reserved_chunk_size, std::ios_base::cur);

    filter_cutoff = read_int_from_istream(ifs);
    const IndexParameters sti_parameters = IndexParameters::read(ifs);
    if (parameters != sti_parameters) {
        throw InvalidIndexFile("Index parameters in .sti file and those specified on command line differ");
    }

    read_vector(ifs, randstrobes_vector);
    read_vector(ifs, hash_positions);
}

int estimate_randstrobe_hashes(const std::string& seq, const IndexParameters& parameters) {
    int num = 0;

    auto randstrobe_iter = RandstrobeIterator2(seq, parameters.k, parameters.s, parameters.t_syncmer, parameters.w_min, parameters.w_max, parameters.q, parameters.max_dist);
    Randstrobe randstrobe;
    while ((randstrobe = randstrobe_iter.next()) != randstrobe_iter.end()) {
        num++;
    }
    return num;
}

size_t estimate_randstrobe_hashes_parallel(const References& references, const IndexParameters& parameters, size_t n_threads) {
    std::vector<std::thread> workers;
    int total = 0;
    std::atomic_size_t ref_index = 0;

    std::vector<int> estimators;
    for (size_t i = 0; i < n_threads; ++i) {
        estimators.push_back(0);
    }

    for (size_t i = 0; i < n_threads; ++i) {
        workers.push_back(
            std::thread(
                [&ref_index](const References& references, const IndexParameters& parameters, int& estimator) {
                    while (true) {
                        size_t j = ref_index.fetch_add(1);
                        if (j >= references.size()) {
                            break;
                        }
                        estimator += estimate_randstrobe_hashes(references.sequences[j], parameters);
                    }
                }, std::ref(references), std::ref(parameters), std::ref(estimators[i]))
        );
    }
    for (auto& worker : workers) {
        worker.join();
    }

    for (auto& estimator : estimators) {
        total += estimator;
    }
    return total;
}

void StrobemerIndex::populate(float f, size_t n_threads) {
    stats.tot_strobemer_count = 0;

    Timer estimate_unique;
    auto randstrobe_hashes = estimate_randstrobe_hashes_parallel(references, parameters, n_threads);
    stats.elapsed_unique_hashes = estimate_unique.duration();
    logger.debug() << "Estimated number of randstrobe hashes: " << randstrobe_hashes << '\n';

    Timer randstrobes_timer;
    add_randstrobes_to_vector(randstrobe_hashes);
    stats.elapsed_generating_seeds = randstrobes_timer.duration();

    Timer sorting_timer;
    // sort by hash valuesles
    pdqsort_branchless(randstrobes_vector.begin(), randstrobes_vector.end());
    stats.elapsed_sorting_seeds = sorting_timer.duration();

    Timer hash_index_timer;

    unsigned int offset = 0;
    unsigned int tot_high_ab = 0;
    unsigned int tot_mid_ab = 0;
    std::vector<unsigned int> strobemer_counts;

    unsigned int randstrobe_hash_size = 0;      
    uint64_t prev_hash = -1;
    
    stats.tot_occur_once = 0;
    // calculate stats.flat_vector_size, randstrobe_hash_size
    // add the top N bits of hash to the hash_positions
    // calculate the count of hash that exists more than one time

    hash_positions.reserve(1 << N);

    unsigned int count = 0;
    unsigned int position = 0;
    unsigned int prev_hash_N = 0;
    bool occur_once = false;
    for (auto &mer: randstrobes_vector){
        uint64_t hash_value = mer.hash;
        if (position == 0){
            randstrobe_hash_size += 1;
            occur_once = true;
            unsigned int existing_hash_N = hash_value >> (64 - N);
            for (unsigned int index_temp = 0; index_temp <= existing_hash_N; index_temp++){
                hash_positions.push_back(position);
            }
            prev_hash_N = existing_hash_N;
            position += 1;
            count += 1;
            prev_hash = hash_value;
            continue;
        }

        if (hash_value == prev_hash){
            occur_once = false;
            position += 1;
            count += 1;
            continue;
        }

        if (hash_value != prev_hash){
            randstrobe_hash_size += 1;
            if (occur_once == true){
                stats.tot_occur_once += 1;
            }else{
                occur_once = true;
            }

            if (count > 100){
                tot_high_ab++;
                strobemer_counts.push_back(count);
            }else if (count > 1){
                tot_mid_ab++;
                strobemer_counts.push_back(count);
            }
            count = 1;
            unsigned int existing_hash_N = hash_value >> (64 - N);
            if (existing_hash_N != prev_hash_N){
                for (unsigned int index_temp = prev_hash_N + 1; index_temp <= existing_hash_N; index_temp++){
                    hash_positions.push_back(position);
            }
            prev_hash_N = existing_hash_N;
            }

            position += 1;
            prev_hash = hash_value;
        }

        if (&mer == &randstrobes_vector.back()){
            if (occur_once == true){
                stats.tot_occur_once += 1;
            }
            if (count > 100){
                tot_high_ab++;
                strobemer_counts.push_back(count);
            }else if (count > 1){
                tot_mid_ab++;
                strobemer_counts.push_back(count);
            }
            for (unsigned int index_temp = prev_hash_N + 1; index_temp < (1 << N); index_temp++){
                    hash_positions.push_back(randstrobes_vector.size());
            }
        }
    }

    stats.frac_unique = 1.0 * stats.tot_occur_once / randstrobe_hash_size;
    stats.tot_high_ab = tot_high_ab;
    stats.tot_mid_ab = tot_mid_ab;
    stats.tot_distinct_strobemer_count = randstrobe_hash_size;

    std::sort(strobemer_counts.begin(), strobemer_counts.end(), std::greater<int>());

    // unsigned int index_cutoff = randstrobe_map.size()*f;

    unsigned int index_cutoff= randstrobe_hash_size *f;
    stats.index_cutoff = index_cutoff;
    if (!strobemer_counts.empty()){
        filter_cutoff = index_cutoff < strobemer_counts.size() ?  strobemer_counts[index_cutoff] : strobemer_counts.back();
        filter_cutoff = std::max(30U, filter_cutoff); // cutoff is around 30-50 on hg38. No reason to have a lower cutoff than this if aligning to a smaller genome or contigs.
        filter_cutoff = std::min(100U, filter_cutoff); // limit upper cutoff for normal NAM finding - use rescue mode instead
    } else {
        filter_cutoff = 30;
    }
    stats.filter_cutoff = filter_cutoff;
    stats.elapsed_hash_index = hash_index_timer.duration();
    stats.unique_mers = randstrobe_hash_size;
}

void StrobemerIndex::add_randstrobes_to_vector(int randstrobe_hashes){
    // size_t tot_occur_once = 0;
    randstrobes_vector.reserve(randstrobe_hashes);
    for (size_t ref_index = 0; ref_index < references.size(); ++ref_index) {
        auto seq = references.sequences[ref_index];
        if (seq.length() < parameters.w_max) {
            continue;
        }
        auto randstrobe_iter = RandstrobeIterator2(seq, parameters.k, parameters.s, parameters.t_syncmer, parameters.w_min, parameters.w_max, parameters.q, parameters.max_dist);
        std::vector<Randstrobe> chunk;
        // TODO
        // Chunking makes this function faster, but the speedup is achieved even
        // with a chunk size of 1.
        const size_t chunk_size = 4;
        chunk.reserve(chunk_size);
        bool end = false;
        while (!end) {
            // fill chunk
            Randstrobe randstrobe;
            while (chunk.size() < chunk_size) {
                randstrobe = randstrobe_iter.next();
                if (randstrobe == randstrobe_iter.end()) {
                    end = true;
                    break;
                }
                chunk.push_back(randstrobe);
            }
            stats.tot_strobemer_count += chunk.size();
            for (auto randstrobe : chunk) {
                RefRandstrobeWithHash::packed_t packed = ref_index << 8;
                packed = packed + (randstrobe.strobe2_pos - randstrobe.strobe1_pos);
                // try to insert as direct entry
                randstrobes_vector.push_back(RefRandstrobeWithHash{randstrobe.hash, randstrobe.strobe1_pos, packed});
                }
            chunk.clear();
            }
        }
    // stats.tot_occur_once = tot_occur_once;
}

void StrobemerIndex::print_diagnostics(const std::string& logfile_name, int k) const {
    // Prins to csv file the statistics on the number of seeds of a particular length and what fraction of them them are unique in the index:
    // format:
    // seed_length, count, percentage_unique

    size_t max_size = 100000;
    std::vector<int> log_count(max_size, 0);  // stores count and each index represents the length
    std::vector<int> log_unique(max_size, 0);  // stores count unique and each index represents the length
    std::vector<int> log_repetitive(max_size, 0);  // stores count unique and each index represents the length


    std::vector<uint64_t> log_count_squared(max_size,0);
    uint64_t tot_seed_count = 0;
    uint64_t tot_seed_count_sq = 0;

    std::vector<uint64_t> log_count_1000_limit(max_size, 0);  // stores count and each index represents the length
    uint64_t tot_seed_count_1000_limit = 0;

    size_t seed_length = 0;

    for (int it = 0; it < randstrobes_vector.size(); it++) {
        seed_length = strobe2_offset(it) + k;
        auto count = get_count(it);

        if (seed_length < max_size){
            log_count[seed_length] ++;
            log_count_squared[seed_length] += count;
            tot_seed_count ++;
            tot_seed_count_sq += count;
            if (count <= 1000){
                log_count_1000_limit[seed_length] ++;
                tot_seed_count_1000_limit ++;
            }
        } else {
            // TODO This function should not log anything
            logger.info() << "Detected seed size over " << max_size << " bp (can happen, e.g., over centromere): " << seed_length << std::endl;
        }

        if (count == 1 && seed_length < max_size) {
            log_unique[seed_length]++;
        }
        if (count >= 10 && seed_length < max_size) {
            log_repetitive[seed_length]++;
        }
    }

    // printing
    std::ofstream log_file;
    log_file.open(logfile_name);

    for (size_t i = 0; i < log_count.size(); ++i) {
        if (log_count[i] > 0) {
            double e_count = log_count_squared[i] / log_count[i];
            log_file << i << ',' << log_count[i] << ',' << e_count << std::endl;
        }
    }

    // Get median
    size_t n = 0;
    int median = 0;
    for (size_t i = 0; i < log_count.size(); ++i) {
        n += log_count[i];
        if (n >= tot_seed_count/2) {
            break;
        }
    }
    // Get median 1000 limit
    size_t n_lim = 0;
    for (size_t i = 0; i < log_count_1000_limit.size(); ++i) {
        n_lim += log_count_1000_limit[i];
        if (n_lim >= tot_seed_count_1000_limit/2) {
            break;
        }
    }

    log_file << "E_size for total seeding wih max seed size m below (m, tot_seeds, E_hits)" << std::endl;
    double e_hits = (double) tot_seed_count_sq/ (double) tot_seed_count;
    double fraction_masked = 1.0 - (double) tot_seed_count_1000_limit/ (double) tot_seed_count;
    log_file << median << ',' << tot_seed_count << ',' << e_hits << ',' << 100*fraction_masked << std::endl;
}

unsigned int StrobemerIndex::get_count_line_search(const unsigned int position) const {
    const auto hash = get_hash(position);

    unsigned int count = 0;
    // step can be any number that is a power of 2, but a large number works
    // very well
    unsigned int step = 512;
    while (get_hash(position + count + step) == hash) {
        count += step;
        step *= 2;
    }
    while (step > 1) {
        while (get_hash(position + count + step) == hash) {
            count += step;
        }
        step /= 2;
    }
    while (get_hash(position + count) == hash) {
        ++count;
    }
    return count;
}

