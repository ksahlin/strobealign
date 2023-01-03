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

#include "timer.hpp"

/* Create an IndexParameters instance based on a given read length.
 * k and/or s can be specified explicitly by setting them to a value other than
 * -1, but otherwise reasonable defaults are used for them as well.
 */
IndexParameters IndexParameters::from_read_length(int read_length, int c, int k, int s, int max_seed_len) {
    int l, u;
    struct settings {
        int r_threshold;
        int k;
        int s_offset;
        int l;
        int u;
    };
    std::vector<settings> d = {
        settings {75, 20, -4, -4, 2},
        settings {125, 20, -4, -2, 2},
        settings {175, 20, -4, 1, 7},
        settings {275, 20, -4, 4, 13},
        settings {375, 22, -4, 2, 12},
        settings {std::numeric_limits<int>::max(), 23, -6, 2, 12},
    };
    for (const auto& v : d) {
        if (read_length <= v.r_threshold) {
            if (k == -1) {
                k = v.k;
            }
            if (s == -1) {
                s = k + v.s_offset;
            }
            l = v.l;
            u = v.u;
            break;
        }
    }

    int max_dist;
    if (max_seed_len == -1) {
        max_dist = std::max(read_length - 70, k);
        max_dist = std::min(255, max_dist);
    } else {
        max_dist = max_seed_len - k; // convert to distance in start positions
    }
    int q = std::pow(2, c) - 1;
    return IndexParameters(k, s, l, u, q, max_dist);
}

void write_int_to_ostream(std::ostream& os, int value) {
    int val;
    val = value;
    os.write(reinterpret_cast<const char*>(&val), sizeof(val));
}

void IndexParameters::write(std::ostream& os) const {
    write_int_to_ostream(os, k);
    write_int_to_ostream(os, s);
    write_int_to_ostream(os, l);
    write_int_to_ostream(os, u);
    write_int_to_ostream(os, q);
    write_int_to_ostream(os, max_dist);
}

int read_int_from_istream(std::istream& is) {
    int val;
    is.read(reinterpret_cast<char*>(&val), sizeof(val));
    return val;
}

IndexParameters IndexParameters::read(std::istream& is) {
    int k = read_int_from_istream(is);
    int s = read_int_from_istream(is);
    int l = read_int_from_istream(is);
    int u = read_int_from_istream(is);
    int q = read_int_from_istream(is);
    int max_dist = read_int_from_istream(is);
    return IndexParameters(k, s, l, u, q, max_dist);
}

bool IndexParameters::operator==(const IndexParameters& other) const {
    return
        this->k == other.k
        && this->s == other.s
        && this->l == other.l
        && this->u == other.u
        && this->q == other.q
        && this->max_dist == other.max_dist
        && this->t_syncmer == other.t_syncmer
        && this->w_min == other.w_min
        && this->w_max == other.w_max;
}

uint64_t count_unique_hashes(const ind_mers_vector& mers){
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

void StrobemerIndex::index_vector(const ind_mers_vector &mers, float f) {
    stats.flat_vector_size = mers.size();

    unsigned int offset = 0;
    unsigned int prev_offset = 0;
    unsigned int count = 0;

    unsigned int tot_occur_once = 0;
    unsigned int tot_high_ab = 0;
    unsigned int tot_mid_ab = 0;
    std::vector<unsigned int> strobemer_counts;

    auto prev_mer = mers[0];
    uint64_t prev_hash = mers[0].hash;
    uint64_t curr_hash;

    for (auto &mer : mers) {
        flat_vector.push_back(ReferenceMer{mer.position, mer.packed});

        curr_hash = mer.hash;
        if (curr_hash == prev_hash){
            count++;
        }
        else {
            if (count == 1) {
                tot_occur_once++;
            }
            else if (count > 100){
                tot_high_ab++;
                strobemer_counts.push_back(count);
            }
            else{
                tot_mid_ab++;
                strobemer_counts.push_back(count);
            }

            if (count == 1) {
                add_entry(prev_hash, prev_mer.position, prev_mer.packed | 0x8000'0000);
                flat_vector[flat_vector.size() - 2] = flat_vector[flat_vector.size() - 1];
                flat_vector.pop_back();
                offset--;
            } else {
                add_entry(prev_hash, prev_offset, count);
            }
            count = 1;
            prev_hash = curr_hash;
            prev_offset = offset;
            prev_mer = mer;
        }
        offset++;
    }

    // last k-mer
    add_entry(curr_hash, prev_offset, count);

    float frac_unique = ((float) tot_occur_once )/ mers_index.size();
    stats.tot_strobemer_count = offset;
    stats.tot_occur_once = tot_occur_once;
    stats.frac_unique = frac_unique;
    stats.tot_high_ab = tot_high_ab;
    stats.tot_mid_ab = tot_mid_ab;
    stats.tot_distinct_strobemer_count = mers_index.size();

    std::sort(strobemer_counts.begin(), strobemer_counts.end(), std::greater<int>());

    unsigned int index_cutoff = mers_index.size()*f;
    stats.index_cutoff = index_cutoff;
    unsigned int filter_cutoff;
    if (!strobemer_counts.empty()){
        filter_cutoff = index_cutoff < strobemer_counts.size() ?  strobemer_counts[index_cutoff] : strobemer_counts.back();
        filter_cutoff = std::max(30U, filter_cutoff); // cutoff is around 30-50 on hg38. No reason to have a lower cutoff than this if aligning to a smaller genome or contigs.
        filter_cutoff = std::min(100U, filter_cutoff); // limit upper cutoff for normal NAM finding - use rescue mode instead
    } else {
        filter_cutoff = 30;
    }
    stats.filter_cutoff = filter_cutoff;
}

void StrobemerIndex::write(const std::string& filename) const {
    std::ofstream ofs(filename, std::ios::binary);

    ofs.write("STI\1", 4); // magic number

    write_int_to_ostream(ofs, filter_cutoff);
    parameters.write(ofs);

    write_vector(ofs, flat_vector);

    //write mers_index:
    auto size = uint64_t(mers_index.size());
    ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
    for (auto& p : mers_index) {
        ofs.write(reinterpret_cast<const char*>(&p.first), sizeof(p.first));
        ofs.write(reinterpret_cast<const char*>(&p.second), sizeof(p.second));
    }
}

void StrobemerIndex::read(const std::string& filename) {
    std::ifstream ifs(filename, std::ios::binary);

    union {
        char s[4];
        uint32_t v;
    } magic;
    ifs.read(magic.s, 4);
    if (magic.v != 0x01495453) { // "STI\1"
        throw InvalidIndexFile("Index file has incorrect format (magic number mismatch)");
    }

    filter_cutoff = read_int_from_istream(ifs);
    const IndexParameters sti_parameters = IndexParameters::read(ifs);
    if (parameters != sti_parameters) {
        throw InvalidIndexFile("Index parameters in .sti file and those specified on command line differ");
    }

    read_vector(ifs, flat_vector);

    uint64_t sz;
    // read mers_index:
    mers_index.clear();
    ifs.read(reinterpret_cast<char*>(&sz), sizeof(sz));
    mers_index.reserve(sz);
    // read in big chunks
    const uint64_t chunk_size = pow(2,20);//4 M => chunks of ~10 MB - The chunk size seem not to be that important
    auto buf_size = std::min(sz, chunk_size) * (sizeof(kmer_lookup::key_type) + sizeof(kmer_lookup::mapped_type));
    std::unique_ptr<char> buf_ptr(new char[buf_size]);
    char* buf2 = buf_ptr.get();
    auto left_to_read = sz;
    while (left_to_read > 0) {
        auto to_read = std::min(left_to_read, chunk_size);
        ifs.read(buf2, to_read * (sizeof(kmer_lookup::key_type) + sizeof(kmer_lookup::mapped_type)));
        //Add the elements directly from the buffer
        for (size_t i = 0; i < to_read; ++i) {
            auto start = buf2 + i * (sizeof(kmer_lookup::key_type) + sizeof(kmer_lookup::mapped_type));
            mers_index[*reinterpret_cast<kmer_lookup::key_type*>(start)] = *reinterpret_cast<kmer_lookup::mapped_type*>(start + sizeof(kmer_lookup::key_type));
        }
        left_to_read -= to_read;
    }
}

void StrobemerIndex::populate(float f) {
    auto ind_flat_vector = generate_and_sort_seeds();
    Timer flat_vector_timer;
    stats.elapsed_flat_vector = flat_vector_timer.duration();

    Timer hash_index_timer;
    mers_index.reserve(count_unique_hashes(ind_flat_vector));
    index_vector(ind_flat_vector, f);
    filter_cutoff = stats.filter_cutoff;
    stats.elapsed_hash_index = hash_index_timer.duration();
    stats.unique_mers = mers_index.size();
}

ind_mers_vector StrobemerIndex::generate_and_sort_seeds() const
{
    Timer randstrobes_timer;
    ind_mers_vector ind_flat_vector; //includes hash - for sorting, will be discarded later
    int expected_sampling = parameters.k - parameters.s + 1;
    int approx_vec_size = references.total_length() / expected_sampling;
    ind_flat_vector.reserve(approx_vec_size);
    for(size_t i = 0; i < references.size(); ++i) {
        randstrobes_reference(ind_flat_vector, parameters.k, parameters.w_min, parameters.w_max, references.sequences[i], i, parameters.s, parameters.t_syncmer, parameters.q, parameters.max_dist);
    }
    stats.elapsed_generating_seeds = randstrobes_timer.duration();

    Timer sorting_timer;
    std::sort(ind_flat_vector.begin(), ind_flat_vector.end());
    stats.elapsed_sorting_seeds = sorting_timer.duration();

    return ind_flat_vector;
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

    size_t seed_length;
    for (auto &it : mers_index) {
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
