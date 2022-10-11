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

#include "logger.hpp"

using std::chrono::high_resolution_clock;

typedef std::vector< uint64_t > hash_vector; //only used during index generation


static Logger& logger = Logger::get();

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

unsigned int index_vector(const hash_vector &h_vector, kmer_lookup &mers_index, float f) {
    logger.debug() << "Flat vector size: " << h_vector.size() << std::endl;
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

            std::tuple<unsigned int, unsigned int> s(prev_offset, count);
            mers_index[prev_k] = s;
            count = 1;
            prev_k = curr_k;
            prev_offset = offset;
        }
        offset ++;
    }

    // last k-mer
    std::tuple<unsigned int, unsigned int> s(prev_offset, count);
    mers_index[curr_k] = s;
    float frac_unique = ((float) tot_occur_once )/ mers_index.size();
    logger.debug()
        << "Total strobemers count: " << offset << std::endl
        << "Total strobemers occur once: " << tot_occur_once << std::endl
        << "Fraction Unique: " << frac_unique << std::endl
        << "Total strobemers highly abundant > 100: " << tot_high_ab << std::endl
        << "Total strobemers mid abundance (between 2-100): " << tot_mid_ab << std::endl
        << "Total distinct strobemers stored: " << mers_index.size() << std::endl;
    if (tot_high_ab >= 1) {
        logger.debug() << "Ratio distinct to highly abundant: " << mers_index.size() / tot_high_ab << std::endl;
    }
    if (tot_mid_ab >= 1) {
        logger.debug() << "Ratio distinct to non distinct: " << mers_index.size() / (tot_high_ab + tot_mid_ab) << std::endl;
    }
    // get count for top -f fraction of strobemer count to filter them out
    std::sort(strobemer_counts.begin(), strobemer_counts.end(), std::greater<int>());

    unsigned int index_cutoff = mers_index.size()*f;
    logger.debug() << "Filtered cutoff index: " << index_cutoff << std::endl;
    unsigned int filter_cutoff;
    if (!strobemer_counts.empty()){
        filter_cutoff =  index_cutoff < strobemer_counts.size() ?  strobemer_counts[index_cutoff] : strobemer_counts.back();
        filter_cutoff = filter_cutoff > 30 ? filter_cutoff : 30; // cutoff is around 30-50 on hg38. No reason to have a lower cutoff than this if aligning to a smaller genome or contigs.
        filter_cutoff = filter_cutoff > 100 ? 100 : filter_cutoff; // limit upper cutoff for normal NAM finding - use rescue mode instead
    } else {
        filter_cutoff = 30;
    }
    logger.debug() << "Filtered cutoff count: " << filter_cutoff << std::endl << std::endl;
    return filter_cutoff;
}


void StrobemerIndex::write(const References& references, const std::string& filename) const {
    std::ofstream ofs(filename, std::ios::binary);

    //write filter_cutoff
    ofs.write(reinterpret_cast<const char*>(&filter_cutoff), sizeof(filter_cutoff));

    //write ref_seqs:
    uint64_t s1 = uint64_t(references.sequences.size());
    ofs.write(reinterpret_cast<char*>(&s1), sizeof(s1));
    //For each string, write length and then the string
    uint32_t s2 = 0;
    for (std::size_t i = 0; i < references.sequences.size(); ++i) {
        s2 = uint32_t(references.sequences[i].length());
        ofs.write(reinterpret_cast<char*>(&s2), sizeof(s2));
        ofs.write(references.sequences[i].c_str(), references.sequences[i].length());
    }
    
    //write ref_lengths:
    //write everything in one large chunk
    s1 = uint64_t(references.lengths.size());
    ofs.write(reinterpret_cast<char*>(&s1), sizeof(s1));
    ofs.write(reinterpret_cast<const char*>(&references.lengths[0]), references.lengths.size()*sizeof(references.lengths[0]));

    //write acc_map:
    s1 = uint64_t(references.names.size());
    ofs.write(reinterpret_cast<char*>(&s1), sizeof(s1));
    //For each string, write length and then the string
    for (std::size_t i = 0; i < references.names.size(); ++i) {
        s2 = uint32_t(references.names[i].length());
        ofs.write(reinterpret_cast<char*>(&s2), sizeof(s2));
        ofs.write(references.names[i].c_str(), references.names[i].length());
    }

    //write flat_vector:
    s1 = uint64_t(flat_vector.size());
    ofs.write(reinterpret_cast<char*>(&s1), sizeof(s1));
    ofs.write(reinterpret_cast<const char*>(&flat_vector[0]), flat_vector.size() * sizeof(flat_vector[0]));

    //write mers_index:
    s1 = uint64_t(mers_index.size());
    ofs.write(reinterpret_cast<char*>(&s1), sizeof(s1));
    for (auto& p : mers_index) {
        ofs.write(reinterpret_cast<const char*>(&p.first), sizeof(p.first));
        ofs.write(reinterpret_cast<const char*>(&p.second), sizeof(p.second));
    }
};

void StrobemerIndex::read(References& references, const std::string& filename) {
    std::ifstream ifs(filename, std::ios::binary);
    //read filter_cutoff
    ifs.read(reinterpret_cast<char*>(&filter_cutoff), sizeof(filter_cutoff));

    //read ref_seqs:
    references.sequences.clear();
    uint64_t sz = 0;
    ifs.read(reinterpret_cast<char*>(&sz), sizeof(sz));
    references.sequences.reserve(sz);
    uint32_t sz2 = 0;
    auto& refseqs = references.sequences;
    for (uint64_t i = 0; i < sz; ++i) {
        ifs.read(reinterpret_cast<char*>(&sz2), sizeof(sz2));
        std::unique_ptr<char> buf_ptr(new char[sz2]);//The vector is short with large strings, so allocating this way should be ok.
        ifs.read(buf_ptr.get(), sz2);
        //we could potentially use std::move here to avoid reallocation, something like std::string(std::move(buf), sz2), but it has to be investigated more
        refseqs.push_back(std::string(buf_ptr.get(), sz2));
    }

    //read ref_lengths
    ////////////////
    references.lengths.clear();
    ifs.read(reinterpret_cast<char*>(&sz), sizeof(sz));
    references.lengths.resize(sz); //annoyingly, this initializes the memory to zero (which is a waste of performance), but ignore that for now
    ifs.read(reinterpret_cast<char*>(&references.lengths[0]), sz*sizeof(references.lengths[0]));

    //read acc_map:
    references.names.clear();
    ifs.read(reinterpret_cast<char*>(&sz), sizeof(sz));
    references.names.reserve(sz);
    auto& acc_map = references.names;

    for (int i = 0; i < sz; ++i) {
        ifs.read(reinterpret_cast<char*>(&sz2), sizeof(sz2));
        std::unique_ptr<char> buf_ptr(new char[sz2]);
        ifs.read(buf_ptr.get(), sz2);
        acc_map.push_back(std::string(buf_ptr.get(), sz2));
    }

    //read flat_vector:
    flat_vector.clear();
    ifs.read(reinterpret_cast<char*>(&sz), sizeof(sz));
    flat_vector.resize(sz); //annoyingly, this initializes the memory to zero (which is a waste of performance), but let's ignore that for now
    ifs.read(reinterpret_cast<char*>(&flat_vector[0]), sz*sizeof(flat_vector[0]));

    //read mers_index:
    mers_index.clear();
    ifs.read(reinterpret_cast<char*>(&sz), sizeof(sz));
    mers_index.reserve(sz);
    //read in big chunks
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
};


void StrobemerIndex::populate(const References& references, mapping_params& map_param) {
    auto start_flat_vector = high_resolution_clock::now();
    hash_vector h_vector;
    {
        auto ind_flat_vector = generate_seeds(references, map_param);

        //Split up the sorted vector into a vector with the hash codes and the flat vector to keep in the index.
        //The hash codes are only needed when generating the index and can be discarded afterwards.
        //We want to do this split-up before creating the hash table to avoid a memory peak - the flat_vector is
        //smaller - doubling that size temporarily will not cause us to go above peak memory.
        auto start_copy_flat_vector = high_resolution_clock::now();
        flat_vector.reserve(ind_flat_vector.size());
        h_vector.reserve(ind_flat_vector.size());
        for (std::size_t i = 0; i < ind_flat_vector.size(); ++i) {
            flat_vector.push_back(std::make_tuple(ind_flat_vector[i].position, ind_flat_vector[i].packed));
            h_vector.push_back(ind_flat_vector[i].hash);
        }
        std::chrono::duration<double> elapsed_copy_flat_vector = high_resolution_clock::now() - start_copy_flat_vector;
        logger.info() << "Time copying flat vector: " << elapsed_copy_flat_vector.count() << " s" << std::endl;

        // ind_flat_vector is freed here
    }
    uint64_t unique_mers = count_unique_elements(h_vector);
    logger.debug() << "Unique strobemers: " << unique_mers << std::endl;

    std::chrono::duration<double> elapsed_flat_vector = high_resolution_clock::now() - start_flat_vector;
    logger.info() << "Total time generating flat vector: " << elapsed_flat_vector.count() << " s" <<  std::endl;

    auto start_hash_index = high_resolution_clock::now();

    mers_index.reserve(unique_mers);
    // construct index over flat array
    map_param.filter_cutoff = index_vector(h_vector, mers_index, map_param.f);
    std::chrono::duration<double> elapsed_hash_index = high_resolution_clock::now() - start_hash_index;
    logger.info() << "Total time generating hash table index: " << elapsed_hash_index.count() << " s" <<  std::endl;
}

ind_mers_vector StrobemerIndex::generate_seeds(const References& references, const mapping_params& map_param) const
{
    auto start_flat_vector = high_resolution_clock::now();

    ind_mers_vector ind_flat_vector; //includes hash - for sorting, will be discarded later
    int expected_sampling = map_param.k - map_param.s + 1;
    int approx_vec_size = references.total_length() / expected_sampling;
    logger.debug() << "ref vector approximate size: " << approx_vec_size << std::endl;
    ind_flat_vector.reserve(approx_vec_size);
    for(size_t i = 0; i < references.size(); ++i) {
        seq_to_randstrobes2(ind_flat_vector, map_param.n, map_param.k, map_param.w_min, map_param.w_max, references.sequences[i], i, map_param.s, map_param.t_syncmer, map_param.q, map_param.max_dist);
    }
    logger.debug() << "Ref vector actual size: " << ind_flat_vector.size() << std::endl;

    std::chrono::duration<double> elapsed_generating_seeds = high_resolution_clock::now() - start_flat_vector;
    logger.info() << "Time generating seeds: " << elapsed_generating_seeds.count() << " s" <<  std::endl;

    auto start_sorting = high_resolution_clock::now();
    std::sort(ind_flat_vector.begin(), ind_flat_vector.end());
    std::chrono::duration<double> elapsed_sorting_seeds = high_resolution_clock::now() - start_sorting;
    logger.info() << "Time sorting seeds: " << elapsed_sorting_seeds.count() << " s" <<  std::endl;

    return ind_flat_vector;
}
