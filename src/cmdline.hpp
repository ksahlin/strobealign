#ifndef STROBEALIGN_CMDLINE_HPP
#define STROBEALIGN_CMDLINE_HPP

#include <vector>
#include <string>
#include <utility>

struct CommandLineOptions {
    int n_threads { 1 };
    int chunk_size {
#ifdef TRACE
        1
#else
        10000
#endif
    };

    // Input/output
    std::string output_file_name;
    bool write_to_stdout { true };
    bool verbose { false };
    bool show_progress { true };
    std::string logfile_name { "" };
    bool only_gen_index { false };
    bool use_index { false };
    bool is_sam_out { true };
    bool is_abundance_out {false};

    // SAM output
    bool cigar_eqx { false };
    bool pg_header { true };
    bool output_unmapped { true };
    std::string read_group_id { "" };
    std::vector<std::string> read_group_fields;
    bool details{false};
    bool fastq_comments{false};
    int max_secondary { 0 };

    // Seeding
    int r { 150 };
    int bits { -1 };
    bool r_set { false };
    bool max_seed_len_set { false };
    bool k_set { false };
    bool s_set { false };
    bool l_set { false };
    bool u_set { false };
    bool c_set { false };
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
    int end_bonus { 10 };

    // Search parameters
    float f { 0.0002 };
    float dropoff_threshold { 0.5 };
    int max_tries { 20 };
    int rescue_level { 2 };

    // Reference and read files
    std::string ref_filename; // This is either a fasta file or an index file - if fasta, indexing will be run
    std::string reads_filename1;
    std::string reads_filename2;
    bool is_SE { true };
    bool is_interleaved { false };
};

CommandLineOptions parse_command_line_arguments(int argc, char **argv);

#endif
