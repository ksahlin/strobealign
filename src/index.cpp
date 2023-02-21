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

static const uint32_t STI_FILE_FORMAT_VERSION = 1;


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

    write_vector(ofs, flat_vector);

    // write mers_index
    auto size = uint64_t(randstrobe_map.size());
    ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
    for (auto& p : randstrobe_map) {
        ofs.write(reinterpret_cast<const char*>(&p.first), sizeof(p.first));
        ofs.write(reinterpret_cast<const char*>(&p.second), sizeof(p.second));
    }
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

    read_vector(ifs, flat_vector);

    uint64_t sz;
    // read mers_index:
    randstrobe_map.clear();
    ifs.read(reinterpret_cast<char*>(&sz), sizeof(sz));
    randstrobe_map.reserve(sz);
    // read in big chunks
    const uint64_t chunk_size = pow(2,20);//4 M => chunks of ~10 MB - The chunk size seem not to be that important
    auto buf_size = std::min(sz, chunk_size) * (sizeof(RandstrobeMap::key_type) + sizeof(RandstrobeMap::mapped_type));
    std::unique_ptr<char> buf_ptr(new char[buf_size]);
    char* buf2 = buf_ptr.get();
    auto left_to_read = sz;
    while (left_to_read > 0) {
        auto to_read = std::min(left_to_read, chunk_size);
        ifs.read(buf2, to_read * (sizeof(RandstrobeMap::key_type) + sizeof(RandstrobeMap::mapped_type)));
        //Add the elements directly from the buffer
        for (size_t i = 0; i < to_read; ++i) {
            auto start = buf2 + i * (sizeof(RandstrobeMap::key_type) + sizeof(RandstrobeMap::mapped_type));
            randstrobe_map[*reinterpret_cast<RandstrobeMap::key_type*>(start)] = *reinterpret_cast<RandstrobeMap::mapped_type*>(start + sizeof(RandstrobeMap::key_type));
        }
        left_to_read -= to_read;
    }
}

hll::HyperLogLog estimate_unique_randstrobe_hashes(const std::string& seq, const IndexParameters& parameters) {
    hll::HyperLogLog hll(10);

    auto randstrobe_iter = RandstrobeIterator2(seq, parameters.k, parameters.s, parameters.t_syncmer, parameters.w_min, parameters.w_max, parameters.q, parameters.max_dist);
    Randstrobe randstrobe;
    while ((randstrobe = randstrobe_iter.next()) != randstrobe_iter.end()) {
        hll.add(reinterpret_cast<char*>(&randstrobe.hash), sizeof(randstrobe.hash));
    }
    return hll;
}

size_t estimate_unique_randstrobe_hashes_parallel(const References& references, const IndexParameters& parameters, size_t n_threads) {
    std::vector<std::thread> workers;
    std::vector<hll::HyperLogLog> estimators;
    for (size_t i = 0; i < n_threads; ++i) {
        estimators.push_back(hll::HyperLogLog(10));
    }
    std::atomic_size_t ref_index = 0;
    for (size_t i = 0; i < n_threads; ++i) {
        workers.push_back(
            std::thread(
                [&ref_index](const References& references, const IndexParameters& parameters, hll::HyperLogLog& estimator) {
                    while (true) {
                        size_t j = ref_index.fetch_add(1);
                        if (j >= references.size()) {
                            break;
                        }
                        estimator.merge(estimate_unique_randstrobe_hashes(references.sequences[j], parameters));
                    }
                }, std::ref(references), std::ref(parameters), std::ref(estimators[i]))
        );
    }
    for (auto& worker : workers) {
        worker.join();
    }

    hll::HyperLogLog hll(10);
    for (auto& estimator : estimators) {
        hll.merge(estimator);
    }
    return hll.estimate();
}

void StrobemerIndex::populate(float f, size_t n_threads) {
    stats.tot_strobemer_count = 0;

    Timer estimate_unique;
    auto randstrobe_hashes = estimate_unique_randstrobe_hashes_parallel(references, parameters, n_threads);
    stats.elapsed_unique_hashes = estimate_unique.duration();
    logger.debug() << "Estimated number of unique randstrobe hashes: " << randstrobe_hashes << '\n';
    randstrobe_map.reserve(randstrobe_hashes);

    Timer randstrobes_timer;
    auto ind_flat_vector = add_randstrobes_to_hash_table();
    stats.elapsed_generating_seeds = randstrobes_timer.duration();

    Timer sorting_timer;
    // sort by hash values
    pdqsort_branchless(ind_flat_vector.begin(), ind_flat_vector.end());
    stats.elapsed_sorting_seeds = sorting_timer.duration();

    Timer hash_index_timer;
    stats.flat_vector_size = ind_flat_vector.size();

    unsigned int offset = 0;
    unsigned int tot_high_ab = 0;
    unsigned int tot_mid_ab = 0;
    std::vector<unsigned int> strobemer_counts;

    uint64_t prev_hash = -1;
    flat_vector.reserve(ind_flat_vector.size());
    for (auto &mer : ind_flat_vector) {
        flat_vector.push_back(RefRandstrobe{mer.position, mer.packed});
        if (mer.hash != prev_hash) {
            auto mer_index_entry = randstrobe_map.find(mer.hash);
            assert(mer_index_entry != randstrobe_map.end());
            auto count = mer_index_entry->second.count();
            assert(count > 1);

            mer_index_entry->second.set_offset(offset);
            if (count > 100){
                tot_high_ab++;
                strobemer_counts.push_back(count);
            } else {
                tot_mid_ab++;
                strobemer_counts.push_back(count);
            }
            count = 1;
        }
        prev_hash = mer.hash;
        offset++;
    }
    stats.frac_unique = 1.0 * stats.tot_occur_once / randstrobe_map.size();
    stats.tot_high_ab = tot_high_ab;
    stats.tot_mid_ab = tot_mid_ab;
    stats.tot_distinct_strobemer_count = randstrobe_map.size();

    std::sort(strobemer_counts.begin(), strobemer_counts.end(), std::greater<int>());

    unsigned int index_cutoff = randstrobe_map.size()*f;
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
    stats.unique_mers = randstrobe_map.size();
}

/*
 * Generate randstrobes for all reference sequences and add them to the hash
 * table. Only those randstrobes which occur only once have correct, fully
 * filled-in entries in the hash table. For the others (with multiple
 * occurrences), only their count is correct. The offset needs to be filled in
 * later.
 *
 * Fills in
 * - stats.tot_occur_once
 * - stats.tot_strobemer_count
 */
std::vector<RefRandstrobeWithHash> StrobemerIndex::add_randstrobes_to_hash_table() {
    std::vector<RefRandstrobeWithHash> randstrobes_with_hash;
    size_t tot_occur_once = 0;
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
                RandstrobeMapEntry entry{randstrobe.strobe1_pos, packed | 0x8000'0000};
                auto result = randstrobe_map.insert({randstrobe.hash, entry});
                if (result.second) {
                    tot_occur_once++;
                } else {
                    // already exists in hash table
                    auto existing = result.first;
                    auto existing_count = existing->second.count();
                    if (existing_count == 1) {
                        // current entry is a direct one, convert to an indirect one
                        auto existing_randstrobe = existing->second.as_ref_randstrobe();
                        randstrobes_with_hash.push_back(RefRandstrobeWithHash{randstrobe.hash, existing_randstrobe.position, existing_randstrobe.packed()});
                        tot_occur_once--;
                    }
                    // offset is adjusted later after sorting
                    existing->second.set_count(existing_count + 1);

                    randstrobes_with_hash.push_back(RefRandstrobeWithHash{randstrobe.hash, randstrobe.strobe1_pos, packed});
                }
            }
            chunk.clear();
        }
    }
    stats.tot_occur_once = tot_occur_once;
    return randstrobes_with_hash;
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
    for (auto &it : randstrobe_map) {
        auto ref_mer = it.second;
        auto offset = ref_mer.offset();
        auto count = ref_mer.count();

        for (size_t j = offset; j < offset + count; ++j) {
            auto r = flat_vector[j];
            seed_length = r.strobe2_offset() + k;
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
               // logger.info() << "Detected seed size over " << max_size << " bp (can happen, e.g., over centromere): " << seed_length << std::endl;
            }
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
