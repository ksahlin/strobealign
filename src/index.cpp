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
#include "poolstl/poolstl.hpp"
#include <iostream>
#include <thread>
#include <atomic>
#include "io.hpp"
#include "timer.hpp"
#include "logger.hpp"
#include <sstream>

static Logger& logger = Logger::get();
static const uint32_t STI_FILE_FORMAT_VERSION = 3;


namespace {

uint64_t count_randstrobes(const std::string& seq, const IndexParameters& parameters) {
    uint64_t n_syncmers = 0;
    SyncmerIterator syncmer_iterator(seq, parameters.syncmer);
    Syncmer syncmer;
    while (!(syncmer = syncmer_iterator.next()).is_end()) {
        n_syncmers++;
    }
    // The last w_min syncmers do not result in a randstrobe
    if (n_syncmers < parameters.randstrobe.w_min) {
        return 0;
    } else {
        return n_syncmers - parameters.randstrobe.w_min;
    }
}

std::vector<uint64_t> count_all_randstrobes(const References& references, const IndexParameters& parameters, size_t n_threads) {
    std::vector<std::thread> workers;
    std::atomic_size_t ref_index{0};

    std::vector<uint64_t> counts;
    counts.assign(references.size(), 0);

    for (size_t i = 0; i < n_threads; ++i) {
        workers.push_back(
            std::thread(
                [&ref_index](const References& references, const IndexParameters& parameters, std::vector<uint64_t>& counts) {
                    while (true) {
                        size_t j = ref_index.fetch_add(1);
                        if (j >= references.size()) {
                            break;
                        }
                        counts[j] = count_randstrobes(references.sequences[j], parameters);
                    }
                }, std::ref(references), std::ref(parameters), std::ref(counts))
        );
    }
    for (auto& worker : workers) {
        worker.join();
    }

    return counts;
}

}

void StrobemerIndex::write(const std::string& filename) const {
    std::ofstream ofs(filename, std::ios::binary);

    ofs.write("STI\1", 4); // Magic number
    write_int_to_ostream(ofs, STI_FILE_FORMAT_VERSION);

    // Variable-length chunk reserved for future use
    std::vector<char> reserved_chunk{0, 0, 0, 0, 0, 0, 0, 0};
    write_vector(ofs, reserved_chunk);

    write_int_to_ostream(ofs, filter_cutoff);
    write_int_to_ostream(ofs, bits);
    parameters.write(ofs);

    write_vector(ofs, randstrobes);
    write_vector(ofs, randstrobe_start_indices);
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
    randstrobe_hash_t reserved_chunk_size;
    ifs.read(reinterpret_cast<char*>(&reserved_chunk_size), sizeof(reserved_chunk_size));
    ifs.seekg(reserved_chunk_size, std::ios_base::cur);

    filter_cutoff = read_int_from_istream(ifs);
    bits = read_int_from_istream(ifs);
    const IndexParameters sti_parameters = IndexParameters::read(ifs);
    if (parameters != sti_parameters) {
        throw InvalidIndexFile("Index parameters in .sti file and those specified on command line differ");
    }

    read_vector(ifs, randstrobes);
    read_vector(ifs, randstrobe_start_indices);
    if (randstrobe_start_indices.size() != (1u << bits) + 1) {
        throw InvalidIndexFile("randstrobe_start_indices vector is of the wrong size");
    }
}

/* Pick a suitable number of bits for indexing randstrobe start indices */
int StrobemerIndex::pick_bits(size_t size) const {
    size_t estimated_number_of_randstrobes = size / (parameters.syncmer.k - parameters.syncmer.s + 1);
    // Two randstrobes per bucket on average
    return std::clamp(static_cast<int>(log2(estimated_number_of_randstrobes)) - 1, 8, 31);
}

void StrobemerIndex::populate(float f, unsigned n_threads) {
    Timer count_hash;
    auto randstrobe_counts = count_all_randstrobes(references, parameters, n_threads);
    stats.elapsed_counting_hashes = count_hash.duration();

    uint64_t total_randstrobes = 0;
    for (auto& count : randstrobe_counts) {
        total_randstrobes += count;
    }
    stats.tot_strobemer_count = total_randstrobes;

    logger.debug() << "  Total number of randstrobes: " << total_randstrobes << '\n';
    uint64_t memory_bytes = references.total_length() + sizeof(RefRandstrobe) * total_randstrobes + sizeof(bucket_index_t) * (1u << bits);
    logger.debug() << "  Estimated total memory usage: " << memory_bytes / 1E9 << " GB\n";

    if (total_randstrobes > std::numeric_limits<bucket_index_t>::max()) {
        throw std::range_error("Too many randstrobes");
    }
    Timer randstrobes_timer;
    logger.debug() << "  Generating randstrobes ...\n";
    randstrobes.assign(total_randstrobes, RefRandstrobe{0, 0, 0});
    assign_all_randstrobes(randstrobe_counts, n_threads);
    stats.elapsed_generating_seeds = randstrobes_timer.duration();

    Timer sorting_timer;
    logger.debug() << "  Sorting ...\n";
    task_thread_pool::task_thread_pool pool{n_threads};
    poolstl::pluggable_sort(poolstl::par.on(pool), randstrobes.begin(), randstrobes.end(), pdqsort_branchless);
    stats.elapsed_sorting_seeds = sorting_timer.duration();

    Timer hash_index_timer;
    logger.debug() << "  Indexing ...\n";

    uint64_t tot_high_ab = 0;
    uint64_t tot_mid_ab = 0;
    std::vector<uint64_t> strobemer_counts;

    stats.tot_occur_once = 0;
    randstrobe_start_indices.reserve((1u << bits) + 1);

    uint64_t unique_mers = randstrobes.empty() ? 0 : 1;
    randstrobe_hash_t prev_hash = randstrobes.empty() ? 0 : randstrobes[0].hash();
    unsigned int count = 1;
    // first randstrobe index will always be the 0
    if(!randstrobes.empty()) {
        randstrobe_start_indices.push_back(0);
    }
    for (bucket_index_t position = 1; position < randstrobes.size(); ++position) {
        const randstrobe_hash_t cur_hash = randstrobes[position].hash();
        if (cur_hash == prev_hash) {
            ++count;
            continue;
        }

        ++unique_mers;

        if (count == 1) {
            ++stats.tot_occur_once;
        } else {
            if (count > 100) {
                ++tot_high_ab;
            } else {
                ++tot_mid_ab;
            }
            strobemer_counts.push_back(count);
        }
        count = 1;
        const unsigned int cur_hash_N = cur_hash >> (64 - bits);
        while (randstrobe_start_indices.size() <= cur_hash_N) {
            randstrobe_start_indices.push_back(position);
        }
        prev_hash = cur_hash;
    }
    // wrap up last entry
    if (count == 1) {
        ++stats.tot_occur_once;
    } else {
        if (count > 100) {
            tot_high_ab++;
        } else {
            tot_mid_ab++;
        }
        strobemer_counts.push_back(count);
    }
    while (randstrobe_start_indices.size() < ((1u << bits) + 1)) {
        randstrobe_start_indices.push_back(randstrobes.size());
    }
    stats.tot_high_ab = tot_high_ab;
    stats.tot_mid_ab = tot_mid_ab;

    std::sort(strobemer_counts.begin(), strobemer_counts.end(), std::greater<int>());

    uint64_t index_cutoff = unique_mers * f;
    stats.index_cutoff = index_cutoff;
    if (!strobemer_counts.empty()){
        filter_cutoff = index_cutoff < strobemer_counts.size() ?  strobemer_counts[index_cutoff] : strobemer_counts.back();
        filter_cutoff = std::max(30U, filter_cutoff); // cutoff is around 30-50 on hg38. No reason to have a lower cutoff than this if aligning to a smaller genome or contigs.
        filter_cutoff = std::min(100U, filter_cutoff); // limit upper cutoff for normal NAM finding - use rescue mode instead
    } else {
        filter_cutoff = 30;
    }
    stats.filter_cutoff = filter_cutoff;
    partial_filter_cutoff = filter_cutoff;
    stats.elapsed_hash_index = hash_index_timer.duration();
    stats.distinct_strobemers = unique_mers;
}

void StrobemerIndex::assign_all_randstrobes(const std::vector<uint64_t>& randstrobe_counts, size_t n_threads) {
    // Compute offsets
    std::vector<size_t> offsets;
    size_t offset = 0;
    for (size_t ref_index = 0; ref_index < references.size(); ++ref_index) {
        offsets.push_back(offset);
        offset += randstrobe_counts[ref_index];
    }

    std::vector<std::thread> workers;
    std::atomic_size_t ref_index{0};
    for (size_t i = 0; i < n_threads; ++i) {
        workers.push_back(
            std::thread(
                [&]() {
                    while (true) {
                        size_t j = ref_index.fetch_add(1);
                        if (j >= references.size()) {
                            break;
                        }
                        assign_randstrobes(j, offsets[j]);
                    }
                })
        );
    }
    for (auto& worker : workers) {
        worker.join();
    }
}

/*
 * Compute randstrobes of one reference and assign them to the randstrobes
 * vector starting from the given offset
 */
void StrobemerIndex::assign_randstrobes(size_t ref_index, size_t offset) {
    auto& seq = references.sequences[ref_index];
    if (seq.length() < parameters.randstrobe.w_max) {
        return;
    }
    RandstrobeGenerator randstrobe_iter{seq, parameters.syncmer, parameters.randstrobe};
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
        for (auto randstrobe : chunk) {
            RefRandstrobe::packed_t packed = (ref_index << 9) | (randstrobe.first_strobe_is_main << 8);
            packed = packed + (randstrobe.strobe2_pos - randstrobe.strobe1_pos);
            randstrobes[offset++] = RefRandstrobe{randstrobe.hash, randstrobe.strobe1_pos, packed};
        }
        chunk.clear();
    }
}

void StrobemerIndex::print_diagnostics(const std::string& logfile_name, int k) const {
    // Prins to csv file the statistics on the number of seeds of a particular length and what fraction of them them are unique in the index:
    // format:
    // seed_length, count, percentage_unique

    size_t max_size = 100000;
    std::vector<int> log_count(max_size, 0);  // stores count and each index represents the length
    std::vector<int> log_unique(max_size, 0);  // stores count unique and each index represents the length
    std::vector<int> log_repetitive(max_size, 0);  // stores count unique and each index represents the length


    std::vector<randstrobe_hash_t> log_count_squared(max_size,0);
    randstrobe_hash_t tot_seed_count = 0;
    randstrobe_hash_t tot_seed_count_sq = 0;

    std::vector<randstrobe_hash_t> log_count_1000_limit(max_size, 0);  // stores count and each index represents the length
    randstrobe_hash_t tot_seed_count_1000_limit = 0;

    size_t seed_length = 0;

    for (size_t it = 0; it < randstrobes.size(); it++) {
        seed_length = strobe2_offset(it) + k;
        auto count = get_count_full(find_full(get_hash(it)));

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
    size_t median = 0;
    for (size_t i = 0; i < log_count.size(); ++i) {
        n += log_count[i];
        if (n >= tot_seed_count/2) {
            median = i;
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

    log_file << "E_size for total seeding with median seed size m below (m, tot_seeds, E_hits, fraction_masked)"
             << std::endl;
    double e_hits = (double) tot_seed_count_sq/ (double) tot_seed_count;
    double fraction_masked = 1.0 - (double) tot_seed_count_1000_limit/ (double) tot_seed_count;
    log_file << median << ',' << tot_seed_count << ',' << e_hits << ',' << 100*fraction_masked << std::endl;
}
