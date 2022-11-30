#ifndef CMDLINE_HPP
#define CMDLINE_HPP

#include <string>
#include <utility>

#include "aln.hpp"

struct CommandLineOptions {
    int n_threads { 3 };

    // Input/output
    std::string output_file_name;
    bool write_to_stdout { true };
    bool verbose { false };
    std::string read_group_id { "" };
    std::string logfile_name { "" };
    bool only_gen_index { false };
    bool use_index { false };

    // Seeding
    bool r_set { false };
    bool max_seed_len_set { false };
    bool k_set { false };
    bool s_set { false };
    int max_seed_len;
    int k { 20 };
    int l { 0 };
    int u { 7 };
    int s { 16 };
    int c { 8 };

    // Alignment
    int A { 2 };
    int B { 8 };
    int O { 12 };
    int E { 1 };

    // Search parameters
    float f { 0.0002 };

    // Reference and read files
    std::string ref_filename; // This is either a fasta file or an index file - if fasta, indexing will be run
    std::string reads_filename1;
    std::string reads_filename2;
    std::string index_out_filename;
    bool is_SE { true };
};

std::pair<CommandLineOptions, mapping_params> parse_command_line_arguments(int argc, char **argv);

#endif
