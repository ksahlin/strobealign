#include <iostream>
#include <fstream>
#include <unordered_map>
#include "source/robin_hood.h"
#include <vector>
#include <string>
#include <chrono>  // for high_resolution_clock

#include "source/index.hpp"




//typedef robin_hood::unordered_map< std::string , std::string > queries;
typedef robin_hood::unordered_map< unsigned int , std::string > references;
typedef robin_hood::unordered_map< unsigned int, std::string > idx_to_acc;

typedef robin_hood::unordered_map< uint64_t, std::tuple<uint64_t, unsigned int >> vector_index;

typedef std::vector< std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int>> mers_vector;



static void read_references(references &seqs, idx_to_acc &acc_map, std::string fn)
{
    std::ifstream file(fn);
    std::string line, seq;
    int ref_index = 0;
    while (getline(file, line)) {
        if (line[0] == '>') {
//            std::cout << ref_index << " " << line << std::endl;
            if (seq.length() > 0){
                seqs[ref_index -1] = seq;
//                std::cout << ref_index - 1 << " here " << seq << " " << seq.length() << " " << seq.size() << std::endl;
//                generate_kmers(h, k, seq, ref_index);
            }
            acc_map[ref_index] = line.substr(1, line.length() -1); //line;
            ref_index++;
            seq = "";
        }
        else {
            seq += line;
        }
    }
    if (seq.length() > 0){
        seqs[ref_index -1] = seq;
//        std::cout << ref_index -1 << " here2 " << seq << std::endl;
//        generate_kmers(h, k, seq, ref_index);
    }
    file.close();
}



static inline void print_diagnostics_new4(mers_vector &mers_vector, vector_index mers_index ) {
    uint64_t tot_flat_vector_size = 0;
    for (size_t i = 0; i < mers_vector.size(); ++i)
    {
        // access using []
        auto t = mers_vector[i];
//        std::cout << "(" << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << ", " << std::get<3>(t) << ", " << std::get<4>(t) << "), ";
        tot_flat_vector_size += sizeof(t);
    }
    std::cout << "Total size of flat mers-vector : " << tot_flat_vector_size/1000000  << " Mb." << std::endl;

//    uint64_t tot_hashtable_index_size = 0;
//    for (auto &it : mers_index)
//    {
////        std::cout << it.first << ": (" << std::get<0>(it.second) << ", " << std::get<1>(it.second) << "), " ;
//        tot_hashtable_index_size += sizeof(it.first);
//        tot_hashtable_index_size += sizeof(it.second);
//    }
//    std::cout << "Total size of hash table index : " << tot_hashtable_index_size/1000000  << " Mb." << std::endl;

    std::cout << "Total size of hash table index : " << (mers_index.size() * sizeof(vector_index::value_type))/1000000 << " Mb." << "\n";
}



static inline std::vector<nam> find_nams(mers_vector &query_mers, mers_vector &mers_vector, vector_index &mers_index, int k){
//    std::cout << "ENTER FIND NAMS " <<  std::endl;
    robin_hood::unordered_map< unsigned int, std::vector<hit>> hits_per_ref; // [ref_id] -> vector( struct hit)
    uint64_t hit_count_reduced = 0;
    uint64_t hit_count_all = 0;
    for (auto &q : query_mers)
//    for (size_t i = 0; i < query_mers.size(); ++i)
    {
        hit h;
        h.query_s = std::get<2>(q);
        h.query_e = std::get<4>(q) + k;
//        std::cout << h.query_s << " " << h.query_e <<  std::endl;

        uint64_t mer_hashv = std::get<0>(q);
        if (mers_index.find(mer_hashv) != mers_index.end()){ //  In  index
//            std::cout << "Found: " << h.query_s <<  std::endl;
            std::tuple<uint64_t, unsigned int> mer;
            mer = mers_index[mer_hashv];
            uint64_t offset = std::get<0>(mer);
            unsigned int count = std::get<1>(mer);
            unsigned int prev_ref_s = 0;
            unsigned int prev_ref_e = 0;
            unsigned int prev_ref_id = 0;
            for(size_t j = offset; j < offset+count; ++j)
            {

                auto r = mers_vector[j];
                unsigned int ref_s = std::get<2>(r);
                unsigned int ref_e = std::get<4>(r) + k;
                unsigned int ref_id = std::get<1>(r);

                h.ref_s = ref_s;
                h.ref_e = ref_e;
                hits_per_ref[ref_id].push_back(h);

//                if (j== offset){
////                    std::cout << "INITIALZING! " << std::endl;
//                    prev_ref_s = ref_s;
//                    prev_ref_e = ref_e;
//                    prev_ref_id = ref_id;
//                }
//                else if ( (prev_ref_s < ref_s) && (ref_s < prev_ref_e) && (ref_id == prev_ref_id)  ){
//                    if (ref_e > prev_ref_e){
//                        prev_ref_e = ref_e;
//                    }
//
//                }
//                else{
//                    h.ref_s = prev_ref_s;
//                    h.ref_e = prev_ref_e;
//                    hits_per_ref[prev_ref_id].push_back(h);
//                    hit_count_reduced ++;
//                    prev_ref_s = ref_s;
//                    prev_ref_e = ref_e;
//                    prev_ref_id = ref_id;
////                    std::cout << "REDUCED Hit! " << h.query_s << ", " << h.query_e << ", " << h.ref_s << ", " << h.ref_e << ", " << std::endl;
////                    if ( (h.query_e - h.query_s) < (h.ref_e - h.ref_s) ){
////                        ;
////                        std::cout << "REDUCED Hit! " << h.query_s << ", " << h.query_e << ", " << h.ref_s << ", " << h.ref_e << ", " << std::endl;
////                    }
//                }


                hit_count_all ++;
//                std::cout << "Hit! " << h.query_s << ", " << h.query_e << ", " << ref_s << ", " << ref_e << ", " << std::endl;

            }
//            h.ref_s = prev_ref_s;
//            h.ref_e = prev_ref_e;
//            hits_per_ref[prev_ref_id].push_back(h);
//            hit_count_reduced ++;
//            std::cout << "REDUCED Hit! " << h.query_s << ", " << h.query_e << ", " << h.ref_s << ", " << h.ref_e << ", " << std::endl;
//            if ( (h.query_e - h.query_s) < (h.ref_e - h.ref_s) ){
//                ;
//                std::cout << "REDUCED Hit! " << h.query_s << ", " << h.query_e << ", " << h.ref_s << ", " << h.ref_e << ", " << std::endl;
//            }

        }
    }

//    std::cout << "NUMBER OF HITS GENERATED: " << hit_count_all << std::endl;
//    std::cout << "NUMBER OF REDUCED HITS GENERATED: " << hit_count_reduced << std::endl;
    std::vector<nam> open_nams;
    std::vector<nam> final_nams; // [ref_id] -> vector(struct nam)

    for (auto &it : hits_per_ref)
    {
        unsigned int ref_id = it.first;
        std::vector<hit> hits = it.second;
        open_nams = std::vector<nam> (); // Initialize vector
        uint64_t prev_q_start = 0;
//        uint64_t prev_q_end = 0;
//        uint64_t hit_copy_id = -1;
        for (auto &h : hits){
            bool is_added = false;
//            assert(prev_q_start <= h.query_s && "Not larger!");
//            if (h.query_s % 500000 == 0 ){
//                std::cout << h.query_s << std::endl;
//            }
//            if (h.query_s > 1110 ){
//                return final_nams;
//            }
//            else {
//                hit_copy_id = 0;
//            }

//            std::cout << "HIT " << h.query_s <<  ", " << h.query_e << ", " << h.ref_s <<  ", " << h.ref_e << std::endl;
            for (auto & o : open_nams) {

                // Extend NAM
                if ( ( o.previous_query_start < h.query_s) && (h.query_s <= o.query_e ) && ( o.previous_ref_start <= h.ref_s) && (h.ref_s <= o.ref_e) ){ // && (o.previous_ref_start <= h.ref_s)  && (o.previous_query_start <= h.query_s) && (hit_copy_id <= o.copy_id)
//                    if (o.previous_ref_start <= h.ref_s) {
//                    std::cout << "In " << o.query_s <<  ", " << o.query_e <<  ", " << o.ref_s <<  ", " << o.ref_e << std::endl;
                    if (h.query_e > o.query_e) {
//                        std::cout << "Changing query " << o.query_e <<  " to " << h.query_e << std::endl;
                        o.query_e = h.query_e;
//                        std::cout << o.query_e << std::endl;
                    }
                    if (h.ref_e > o.ref_e) {
//                        std::cout << "Changing ref " << o.ref_e <<  " to " << h.ref_e << std::endl;
                        o.ref_e = h.ref_e;
                    }
//                    if (o.previous_ref_start <)
//                    }
                    o.previous_query_start = h.query_s;
                    o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
                    is_added = true;
//                    o.copy_id = hit_copy_id;
                    break;
                }
//                //Output local match
//                else if (( o.query_s <= h.query_s) && (h.query_s <= o.query_e ) && ( o.ref_s <= h.ref_s) && (h.ref_s <= o.ref_e)){
//                    std::cout << "Output local match " << o.query_s <<  " to " << o.query_e  << " " << o.ref_s <<  " to " << o.ref_e  << "     " << o.previous_ref_start << " " << h.ref_s << std::endl;
//                }

            }
//            prev_q_start = h.query_s;
//            prev_q_end = h.query_e;

//            std::cout << is_added << std::endl;
//            std::cout << " " << std::endl;
            // Add the hit to open matches
            if (not is_added){
                nam n;
                n.query_s = h.query_s;
                n.query_e = h.query_e;
                n.ref_s = h.ref_s;
                n.ref_e = h.ref_e;
                n.ref_id = ref_id;
                n.previous_query_start = h.query_s;
                n.previous_ref_start = h.ref_s;
//                n.copy_id = hit_copy_id;
                open_nams.push_back(n);
            }


//            // Output matches with identical query and reference end coordinates to a longer match. This means that the match is a submatch to another match
//            int before = open_nams.size();
//            // TODO: Also output the matches,, not just remove them as is currently done..
//            // This may happen because of interleaved repeats. Example of an interleaved repeat (with k<= 30) from ecoli with identical strings [0-126] and from [96-222]: CTGTTGCTGTTCCAGCTTGCGCGCTTTGGCACGGGCAATAGCGGCTTCGACGGCGGCTTTGCGCGGATCGACCTGTTCTTCTGGTTCCGCATTAGCCTGTTGCTGTTCCAGCTTGCGCGCTTTGGCACGGGCAATAGCGGCTTCGACGGCGGCTTTGCGCGGATCGACCTGTTCTTCTGGTTCCGCATTAGCCTGTTGCTGTTCCAGCTTGCGCGCTTTGGCGCGGGCGATAGCTGCTTCAACGGCAGTTTTACGTGGATCAGCAACGGTTGCTGCGTCGTTAGTTTGCTGCA
//            auto comp = [] ( const nam& n1, const nam& n2 ) {return (n1.query_e == n2.query_e) && (n1.ref_e == n2.ref_e);};
//            auto pred = []( const nam& n1, const nam& n2 ) {return (n1.ref_e < n2.ref_e) || ((n1.ref_e == n2.ref_e ) && (n1.query_s < n2.query_s)) ;};
//            std::sort(open_nams.begin(),open_nams.end(), pred);
//            auto last = std::unique(open_nams.begin(), open_nams.end(),comp);
//            for (auto &r : open_nams){
//                std::cout << "Deleting: " << r.ref_id << ": (" << r.copy_id << ", " << r.query_s << ", " << r.query_e << ", " << r.ref_s << ", " << r.ref_e << ")" << std::endl;
//            }
//            open_nams.erase(last, open_nams.end());
//            int after = open_nams.size();
//            if (before > after){
//                std::cout << "Removed " << before - after <<  " matches." <<  std::endl;
//            }

            // Only filter if we have advanced at least k nucleotides
            if (h.query_s > prev_q_start + k) {

                // Output all NAMs from open_matches to final_nams that the current hit have passed
                for (auto &n : open_nams) {
                    if (n.query_e < h.query_s) {
                        final_nams.push_back(n);
                    }
                }

                // Remove all NAMs from open_matches that the current hit have passed
                unsigned int c = h.query_s;
                auto predicate = [c](decltype(open_nams)::value_type const &nam) { return nam.query_e < c; };
                open_nams.erase(std::remove_if(open_nams.begin(), open_nams.end(), predicate), open_nams.end());
                prev_q_start = h.query_s;
            }


        }

        // Add all current open_matches to final NAMs
        for (auto &n : open_nams){
//        for (size_t i = 0; i < open_nams.size(); ++i){
            final_nams.push_back(n);
        }
    }

//    for (auto &n : final_nams){
//        std::cout << n.ref_id << ": (" << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e << ")" << std::endl;
//    }


    return final_nams;
}


static inline bool compareByQueryCoord(const nam &a, const nam &b)
{
    // first sort on ref ID, then on query, then on reference
    return (a.ref_id < b.ref_id) ||
           ( (a.ref_id == b.ref_id) && (a.query_s < b.query_s) ) ||
           ((a.ref_id == b.ref_id) && (a.query_s == b.query_s ) && (a.ref_s < b.ref_s)) ;
}

static inline bool compareByQueryLength(const nam &a, const nam &b)
{
    return (a.query_e - a.query_s) < ( b.query_e - b.query_s);
}

static inline void output_nams(std::vector<nam> &nams, std::ofstream &output_file, std::string query_acc, idx_to_acc &acc_map) {
    //Sort hits based on start choordinate on query sequence
    std::sort(nams.begin(), nams.end(), compareByQueryCoord);
    // Output results
    output_file << "> " << query_acc << "\n";
    for (auto &n : nams) {
        output_file << "  " << acc_map[n.ref_id]  << " " << n.ref_s << " " << n.query_s << " " << n.ref_e - n.ref_s << "\n";
//      python: outfile.write("  {0} {1} {2} {3}\n".format(ref_acc, ref_p, q_pos, k))
    }
}

int main (int argc, char *argv[])
{

    ///////////////////// INPUT /////////////////////////
//    std::string filename  = "test_ploy2.txt";
//    std::string reads_filename  = "test_ploy2.txt";

//    std::string filename  = "example_repeats.txt";
//    std::string reads_filename  = "example_repeats.txt";

//    std::string filename  = "example3.txt";
//    std::string reads_filename  = "example3.txt";

//    std::string filename  = "ecoli_repeats.txt";
//    std::string reads_filename  = "ecoli_repeats.txt";

//    std::string filename  = "ecoli_bug.txt";
//    std::string reads_filename  = "ecoli_bug.txt";

//    std::string filename  = "ecoli_randmer_bug.txt";
//    std::string reads_filename  = "ecoli_randmer_bug.txt";

    std::string filename  = "ecoli.fa";
//    std::string reads_filename  = "ecoli.fa";
    std::string reads_filename  = "SRR8187994_1.fasta";

//    std::string filename  = "hg38_chr21.fa";
//    std::string reads_filename  = "hg38_chr21.fa";

//    std::string filename  = "hg21_bug.txt";
//    std::string reads_filename  = "hg21_bug.txt";

//    std::string choice = "kmers";
//    std::string choice = "minstrobes";
//   std::string choice = "hybridstrobes";
    std::string choice = "randstrobes";
    int n = 2;
    int k = 12;
    int w_min = 13;
    int w_max = 40;
    int filter_nams = 0;
    assert(k <= w_min && "k have to be smaller than w_min");
    assert(k <= 32 && "k have to be smaller than 32!");
    references ref_seqs;
    idx_to_acc acc_map;
    read_references(ref_seqs, acc_map, filename);

    //////////////////////////////////////////////////////


    //////////// CREATE INDEX OF REF SEQUENCES /////////////////

    // Record index creation start time
    auto start = std::chrono::high_resolution_clock::now();

    three_pos_index tmp_index; // hash table holding all reference mers

    if (choice == "kmers" ){
        for (auto x : ref_seqs){
            mers_vector kmers; //  kmer hash value, pos, chr_id
            kmers = seq_to_kmers(k, x.second, x.first);
            tmp_index[x.first] = kmers;
        }
    }
    else if (choice == "hybridstrobes" ){
        if (n == 2 ){
            for (auto x : ref_seqs){
                mers_vector hybridstrobes2; // pos, chr_id, kmer hash value
                hybridstrobes2 = seq_to_hybridstrobes2(n, k, w_min, w_max, x.second, x.first);
                tmp_index[x.first] = hybridstrobes2;
            }
        }
//        else if (n == 3){ for (auto x : ref_seqs){generate_hybridstrobe3_index(h, k, x.second, x.first);}}
    }
    else if (choice == "minstrobes" ){
        if (n == 2 ){
            for (auto x : ref_seqs){
//                generate_minstrobe2_index(h, n, k, w_min, w_max, x.second, x.first);
                mers_vector minstrobes2; // pos, chr_id, kmer hash value
                minstrobes2 = seq_to_minstrobes2(n, k, w_min, w_max, x.second, x.first);
                tmp_index[x.first] = minstrobes2;
            }
        }
//        else if (n == 3){ for (auto x : ref_seqs){generate_minstrobe3_index(h, k, x.second, x.first);}}
    }
    else if (choice == "randstrobes" ){
        if (n == 2 ){
            for (auto x : ref_seqs){
                mers_vector randstrobes2; // pos, chr_id, kmer hash value
                randstrobes2 = seq_to_randstrobes2(n, k, w_min, w_max, x.second, x.first);
                tmp_index[x.first] = randstrobes2;
            }
        }
        else if (n == 3){
            for (auto x : ref_seqs){
                mers_vector randstrobes3; // pos, chr_id, kmer hash value
                randstrobes3 = seq_to_randstrobes3(n, k, w_min, w_max, x.second, x.first);
                tmp_index[x.first] = randstrobes3;
            }
        }
    }
    else {
        std::cout << choice << "not implemented : " << std::endl;
    }


    mers_vector all_mers_vector;
    all_mers_vector = construct_flat_vector_three_pos(tmp_index);
    robin_hood::unordered_map< uint64_t, std::tuple<uint64_t, unsigned int >> mers_index; // k-mer -> (offset in flat_vector, occurence count )
    mers_index = index_vector_three_pos(all_mers_vector); // construct index over flat array
    tmp_index.clear();
    print_diagnostics_new4(all_mers_vector, mers_index);

    //////////////////////////////////////////////////////////////////////////

    // Record index creation end time
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Total time generating index: " << elapsed.count() << " s\n" <<  std::endl;



    ///////////////////////////// MAP ///////////////////////////////////////

    // Record matching time
    auto start_map = std::chrono::high_resolution_clock::now();

    std::ifstream query_file(reads_filename);
    std::ofstream output_file;
    output_file.open ("output.tsv");

    std::string line, seq, prev_acc;
    unsigned int q_id = 0;
    unsigned int read_cnt = 0;
    mers_vector query_mers; // pos, chr_id, kmer hash value
    while (getline(query_file, line)) {
        if (line[0] == '>') {
            read_cnt ++;
            if (seq.length() > 0){
                // generate mers here
                if (choice == "kmers" ){
                    query_mers = seq_to_kmers(k, seq, q_id);
                }
                else if (choice == "randstrobes" ){
                    if (n == 2 ){
                        query_mers = seq_to_randstrobes2(n, k, w_min, w_max, seq, q_id);
                    }

                    else if (n == 3){
                        query_mers = seq_to_randstrobes3(n, k, w_min, w_max, seq, q_id);
                    }
                }
//                std::cout << "HERE " << line << std::endl;
                // Find NAMs
//                std::cout << "Processing read: " << prev_acc << " kmers generated: " << query_mers.size() << ", read length: " <<  seq.length() << std::endl;
                std::vector<nam> nams; // (r_id, r_pos_start, r_pos_end, q_pos_start, q_pos_end)
                nams = find_nams(query_mers, all_mers_vector, mers_index, k);
//                std::cout <<  "NAMs generated: " << nams.size() << std::endl;
                // Output results
                output_nams(nams, output_file, prev_acc, acc_map);

//              output_file << "> " <<  prev_acc << "\n";
//              output_file << "  " << ref_acc << " " << ref_p << " " << q_pos << " " << "\n";
//              outfile.write("  {0} {1} {2} {3}\n".format(ref_acc, ref_p, q_pos, k))
                if (read_cnt % 10000 == 0){
                    std::cout << "Processed " << read_cnt << "reads. " << std::endl;
                }
            }
            prev_acc = line.substr(1, line.length() -1);
            seq = "";
            q_id ++;
        }
        else {
            seq += line;
        }
    }
    if (seq.length() > 0){
        if (choice == "kmers" ){
            query_mers = seq_to_kmers(k, seq, q_id);
        }
        else if (choice == "randstrobes" ){
            if (n == 2 ){
                query_mers = seq_to_randstrobes2(n, k, w_min, w_max, seq, q_id);
            }

            else if (n == 3){
                query_mers = seq_to_randstrobes3(n, k, w_min, w_max, seq, q_id);
            }
        }
        // Find NAMs
//        std::cout << "Processing read: " << prev_acc << " kmers generated: " << query_mers.size() << ", read length: " <<  seq.length() << std::endl;
        std::vector<nam> nams; // (r_id, r_pos_start, r_pos_end, q_pos_start, q_pos_end)
        nams = find_nams(query_mers, all_mers_vector, mers_index, k);
//        std::cout <<  "NAMs generated: " << nams.size() << std::endl;
        output_nams(nams, output_file, prev_acc, acc_map);
    }

    query_file.close();
    output_file.close();


    // Record mapping end time
    auto finish_map = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_map = finish_map - start_map;
    std::cout << "Total time mapping: " << elapsed_map.count() << " s\n" <<  std::endl;


    //////////////////////////////////////////////////////////////////////////


    /////////////////////// FIND AND OUTPUT NAMs ///////////////////////////////





}

