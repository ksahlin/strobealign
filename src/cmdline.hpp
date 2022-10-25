#ifndef CMDLINE_HPP
#define CMDLINE_HPP

#include <string>
#include <utility>

#include "index.hpp"

struct CommandLineOptions {
    // Index parameters
    bool k_set { false };
    bool s_set{ false };
    int k { 20 };
    int s { 16 };
    int u { 7 };
    int l { 0 };

    int A { 2 };
    int B { 8 };
    int O { 12 };
    int E { 1 };
    int c { 8 };
    float f { 0.0002 };
    int max_seed_len;
    std::string output_file_name;
    std::string logfile_name { "log.csv" };
    bool verbose { false };
    int n_threads { 3 };
    std::string ref_filename; //This is either a fasta file or an index file - if fasta, indexing will be run
    std::string reads_filename1;
    std::string reads_filename2;
    bool is_SE { true };
    bool write_to_stdout { true };
    bool index_log { false };
    bool r_set { false };
    bool max_seed_len_set { false };
    bool only_gen_index{ false };
    std::string index_out_filename;
};

std::pair<CommandLineOptions, mapping_params> parse_command_line_arguments(int argc, char **argv);

#endif
