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
#include "logger.hpp"

typedef std::vector<uint64_t> hash_vector; //only used during index generation


static Logger& logger = Logger::get();


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

uint64_t count_unique_elements(const hash_vector& h_vector){
    assert(h_vector.size() > 0);
    uint64_t prev_k = h_vector[0];
    uint64_t unique_elements = 1;
    for (auto &curr_k : h_vector) {
        if (curr_k != prev_k) {
            unique_elements ++;
        }
        prev_k = curr_k;
    }
    return unique_elements;
}

IndexCreationStatistics index_vector(const hash_vector &h_vector, kmer_lookup &mers_index, float f) {
    IndexCreationStatistics index_stats;
    index_stats.flat_vector_size = h_vector.size();

    unsigned int offset = 0;
    unsigned int prev_offset = 0;
    unsigned int count = 0;

    unsigned int tot_occur_once = 0;
    unsigned int tot_high_ab = 0;
    unsigned int tot_mid_ab = 0;
    std::vector<unsigned int> strobemer_counts;

    uint64_t prev_k = h_vector[0];
    uint64_t curr_k;

    for ( auto &t : h_vector) {
        curr_k = t;
        if (curr_k == prev_k){
            count ++;
        }
        else {
            if (count == 1){
                tot_occur_once ++;
            }
            else if (count > 100){
                tot_high_ab ++;
                strobemer_counts.push_back(count);
            }
            else{
                tot_mid_ab ++;
                strobemer_counts.push_back(count);
            }

            KmerLookupEntry s{prev_offset, count};
            mers_index[prev_k] = s;
            count = 1;
            prev_k = curr_k;
            prev_offset = offset;
        }
        offset ++;
    }

    // last k-mer
    KmerLookupEntry s{prev_offset, count};
    mers_index[curr_k] = s;
    float frac_unique = ((float) tot_occur_once )/ mers_index.size();

    index_stats.tot_strobemer_count = offset;
    index_stats.tot_occur_once = tot_occur_once;
    index_stats.frac_unique = frac_unique;
    index_stats.tot_high_ab = tot_high_ab;
    index_stats.tot_mid_ab = tot_mid_ab;
    index_stats.tot_distinct_strobemer_count = mers_index.size();

    std::sort(strobemer_counts.begin(), strobemer_counts.end(), std::greater<int>());

    unsigned int index_cutoff = mers_index.size()*f;
    index_stats.index_cutoff = index_cutoff;
    unsigned int filter_cutoff;
    if (!strobemer_counts.empty()){
        filter_cutoff = index_cutoff < strobemer_counts.size() ?  strobemer_counts[index_cutoff] : strobemer_counts.back();
        filter_cutoff = std::max(30U, filter_cutoff); // cutoff is around 30-50 on hg38. No reason to have a lower cutoff than this if aligning to a smaller genome or contigs.
        filter_cutoff = std::min(100U, filter_cutoff); // limit upper cutoff for normal NAM finding - use rescue mode instead
    } else {
        filter_cutoff = 30;
    }
    index_stats.filter_cutoff = filter_cutoff;
    return index_stats;
}


void StrobemerIndex::write(const std::string& filename) const {
    std::ofstream ofs(filename, std::ios::binary);

    ofs.write("STI\1", 4); // magic number

    write_int_to_ostream(ofs, filter_cutoff);
    parameters.write(ofs);

    //write flat_vector:
    auto s1 = uint64_t(flat_vector.size());
    ofs.write(reinterpret_cast<char*>(&s1), sizeof(s1));
    ofs.write(reinterpret_cast<const char*>(&flat_vector[0]), flat_vector.size() * sizeof(flat_vector[0]));

    //write mers_index:
    s1 = uint64_t(mers_index.size());
    ofs.write(reinterpret_cast<char*>(&s1), sizeof(s1));
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

    // read flat_vector:
    uint64_t sz;
    flat_vector.clear();
    ifs.read(reinterpret_cast<char*>(&sz), sizeof(sz));
    flat_vector.resize(sz); //annoyingly, this initializes the memory to zero (which is a waste of performance), but let's ignore that for now
    ifs.read(reinterpret_cast<char*>(&flat_vector[0]), sz*sizeof(flat_vector[0]));

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
        for (int i = 0; i < to_read; ++i) {
            auto start = buf2 + i * (sizeof(kmer_lookup::key_type) + sizeof(kmer_lookup::mapped_type));
            mers_index[*reinterpret_cast<kmer_lookup::key_type*>(start)] = *reinterpret_cast<kmer_lookup::mapped_type*>(start + sizeof(kmer_lookup::key_type));
        }
        left_to_read -= to_read;
    }
}

IndexCreationStatistics StrobemerIndex::populate(float f) {
    Timer flat_vector_timer;
    hash_vector h_vector;

    {
        auto ind_flat_vector = generate_seeds();

        //Split up the sorted vector into a vector with the hash codes and the flat vector to keep in the index.
        //The hash codes are only needed when generating the index and can be discarded afterwards.
        //We want to do this split-up before creating the hash table to avoid a memory peak - the flat_vector is
        //smaller - doubling that size temporarily will not cause us to go above peak memory.
        flat_vector.reserve(ind_flat_vector.size());
        h_vector.reserve(ind_flat_vector.size());
        for (std::size_t i = 0; i < ind_flat_vector.size(); ++i) {
            flat_vector.push_back(ReferenceMer{ind_flat_vector[i].position, ind_flat_vector[i].packed});
            h_vector.push_back(ind_flat_vector[i].hash);
        }
        // ind_flat_vector is freed here
    }
    uint64_t unique_mers = count_unique_elements(h_vector);
    std::chrono::duration<double> elapsed_flat_vector = flat_vector_timer.duration();

    Timer hash_index_timer;
    mers_index.reserve(unique_mers);
    // construct index over flat array
    IndexCreationStatistics index_stats = index_vector(h_vector, mers_index, f);
    filter_cutoff = index_stats.filter_cutoff;
    index_stats.elapsed_hash_index = hash_index_timer.duration();
    index_stats.unique_mers = unique_mers;
    index_stats.elapsed_flat_vector = elapsed_flat_vector;

    return index_stats;
}

ind_mers_vector StrobemerIndex::generate_seeds() const
{
    Timer flat_vector_timer;
    ind_mers_vector ind_flat_vector; //includes hash - for sorting, will be discarded later
    int expected_sampling = parameters.k - parameters.s + 1;
    int approx_vec_size = references.total_length() / expected_sampling;
    logger.debug() << "ref vector approximate size: " << approx_vec_size << std::endl;
    ind_flat_vector.reserve(approx_vec_size);
    for(size_t i = 0; i < references.size(); ++i) {
        randstrobes_reference(ind_flat_vector, parameters.k, parameters.w_min, parameters.w_max, references.sequences[i], i, parameters.s, parameters.t_syncmer, parameters.q, parameters.max_dist);
    }
    logger.debug() << "Ref vector actual size: " << ind_flat_vector.size() << std::endl;

    std::chrono::duration<double> elapsed_generating_seeds = flat_vector_timer.duration();
    logger.info() << "Time generating seeds: " << elapsed_generating_seeds.count() << " s" <<  std::endl;

    Timer sorting_timer;
    std::sort(ind_flat_vector.begin(), ind_flat_vector.end());

    auto elapsed_sorting_seeds = sorting_timer.duration();
    logger.info() << "Time sorting seeds: " << elapsed_sorting_seeds.count() << " s" <<  std::endl;


    return ind_flat_vector;
}
