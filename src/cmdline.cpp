#include "cmdline.hpp"

#include <args.hxx>
#include "version.hpp"


std::tuple<CommandLineOptions, mapping_params, IndexParameters> parse_command_line_arguments(int argc, char **argv) {

    args::ArgumentParser parser("StrobeAlign " VERSION_STRING);
    parser.helpParams.showTerminator = false;
    parser.helpParams.helpindent = 20;
    parser.helpParams.width = 90;
    parser.helpParams.programName = "strobealign";
    parser.helpParams.shortSeparator = " ";

    args::HelpFlag help(parser, "help", "Print help and exit", {'h', "help"});
    args::ValueFlag<int> threads(parser, "threads", "Number of threads [3]", {'t', "threads"});

    args::Group io(parser, "Input/output:");
    args::ValueFlag<std::string> o(parser, "STR", "redirect output to file [stdout]", {'o'});
    args::Flag v(parser, "v", "Verbose output", {'v'});
    args::Flag x(parser, "x", "Only map reads, no base level alignment (produces paf file)", {'x'});
    args::ValueFlag<int> N(parser, "INT", "retain at most INT secondary alignments (is upper bounded by -M, and depends on -S) [0]", {'N'});
    args::ValueFlag<std::string> L(parser, "STR", "Print statistics of indexing to logfile [log.csv]", {'L'});
    args::ValueFlag<std::string> i(parser, "index", "Generates an index (.sti) file", { 'i' });

    args::Group seeding(parser, "Seeding:");
    //args::ValueFlag<int> n(parser, "INT", "number of strobes [2]", {'n'});
    args::ValueFlag<int> r(parser, "INT", "Mean read length. This parameter is estimated from first 500 records in each read file. No need to set this explicitly unless you have a reason. [disabled]", {'r'});
    args::ValueFlag<int> m(parser, "INT", "Maximum seed length. Defaults to r - 50. For reasonable values on -l and -u, the seed length distribution is usually determined by parameters l and u. Then, this parameter is only active in regions where syncmers are very sparse.", {'m'});

    args::ValueFlag<int> k(parser, "INT", "strobe length, has to be below 32. [20]", {'k'});
    args::ValueFlag<int> l(parser, "INT", "Lower syncmer offset from k/(k-s+1). Start sample second syncmer k/(k-s+1) + l syncmers downstream [0]", {'l'});

    args::ValueFlag<int> u(parser, "INT", "Upper syncmer offset from k/(k-s+1). End sample second syncmer k/(k-s+1) + u syncmers downstream [7]", {'u'});
    args::ValueFlag<int> c(parser, "INT", "Bitcount length between 2 and 63. [8]", {'c'});
    args::ValueFlag<int> s(parser, "INT", "Submer size used for creating syncmers [k-4]. Only even numbers on k-s allowed. A value of s=k-4 roughly represents w=10 as minimizer window [k-4]. It is recommended not to change this parameter unless you have a good understanding of syncmers as it will drastically change the memory usage and results with non default values.", {'s'});

    args::Group alignment(parser, "Alignment:");
    args::ValueFlag<int> A(parser, "INT", "matching score [2]", {'A'});
    args::ValueFlag<int> B(parser, "INT", "mismatch penalty [8]", {'B'});
    args::ValueFlag<int> O(parser, "INT", "gap open penalty [12]", {'O'});
    args::ValueFlag<int> E(parser, "INT", "gap extension penalty [1]", {'E'});

    args::Group search(parser, "Search parameters:");
    args::ValueFlag<float> f(parser, "FLOAT", "top fraction of repetitive strobemers to filter out from sampling [0.0002]", {'f'});
    args::ValueFlag<float> S(parser, "FLOAT", "Try candidate sites with mapping score at least S of maximum mapping score [0.5]", {'S'});
    args::ValueFlag<int> M(parser, "INT", "Maximum number of mapping sites to try [20]", {'M'});
    args::ValueFlag<int> R(parser, "INT", "Rescue level. Perform additional search for reads with many repetitive seeds filtered out. This search includes seeds of R*repetitive_seed_size_filter (default: R=2). Higher R than default makes StrobeAlign significantly slower but more accurate. R <= 1 deactivates rescue and is the fastest.", {'R'});

    // <ref.fa> <reads1.fast[a/q.gz]> [reads2.fast[a/q.gz]]
    args::Positional<std::string> ref_filename(parser, "reference", "A pregenerated strobemers index file (.sti) or a reference in FASTA format", args::Options::Required);
    args::Positional<std::string> reads1_filename(parser, "reads1", "Reads 1 in FASTA or FASTQ format, optionally gzip compressed");
    args::Positional<std::string> reads2_filename(parser, "reads2", "Reads 2 in FASTA or FASTQ format, optionally gzip compressed");


    try {
        parser.ParseCLI(argc, argv);
    }
    catch (const args::Completion& e) {
        std::cout << e.what();
        exit(0);
    }
    catch (const args::Help&)
    {
        std::cout << parser;
        exit(0);
    }
    catch (const args::Error& e)
    {
        std::cerr << parser;
        std::cerr << "Error: " << e.what() << std::endl;
        exit(1);
    }

    CommandLineOptions opt;

    mapping_params map_param;
    IndexParameters index_parameters;
    map_param.max_secondary = 0;
    index_parameters.k = 20;
    index_parameters.s = index_parameters.k - 4;
    map_param.f = 0.0002;
    map_param.R = 2;
    map_param.dropoff_threshold = 0.5;
    map_param.maxTries = 20;
    index_parameters.l = 0;
    index_parameters.u = 7;
    map_param.c = 8;
    map_param.r = 150;
    map_param.max_dist = std::min(map_param.r - 50, 255);
    map_param.is_sam_out = true;  // true: align, false: map

    if (threads) { opt.n_threads = args::get(threads); }

    // Input/output
    if (o) { opt.output_file_name = args::get(o); opt.write_to_stdout = false; }
    if (v) { opt.verbose = true; }
    if (x) { map_param.is_sam_out = false; }
    if (N) { map_param.max_secondary = args::get(N); }
    if (L) { opt.logfile_name = args::get(L); opt.index_log = true; }
    if (i) { opt.only_gen_index = true; opt.index_out_filename = args::get(i); }

    // Seeding
    if (r) { map_param.r = args::get(r); opt.r_set = true; }
    if (m) { opt.max_seed_len = args::get(m); opt.max_seed_len_set = true; }
    if (k) { index_parameters.k = args::get(k); opt.k_set = true; }
    if (l) { index_parameters.l = args::get(l); }
    if (u) { index_parameters.u = args::get(u); }
    if (s) { index_parameters.s = args::get(s); opt.s_set = true; }
    if (c) { map_param.c = args::get(c); }

    // Alignment
    // if (n) { n = args::get(n); }
    if (A) { opt.A = args::get(A); }
    if (B) { opt.B = args::get(B); }
    if (O) { opt.O = args::get(O); }
    if (E) { opt.E = args::get(E); }

    // Search parameters
    if (f) { map_param.f = args::get(f); }
    if (S) { map_param.dropoff_threshold = args::get(S); }
    if (M) { map_param.maxTries = args::get(M); }
    if (R) { map_param.R = args::get(R); }

    // Reference and read files
    opt.ref_filename = args::get(ref_filename);
    opt.reads_filename1 = args::get(reads1_filename);

    if (reads2_filename) {
        opt.reads_filename2 = args::get(reads2_filename);
        opt.is_SE = false;
    } else {
        opt.reads_filename2 = std::string();
        opt.is_SE = true;
    }

    //If not generating index, fastq1 is mandatory:
    if (opt.reads_filename1.empty() && !opt.only_gen_index) {
        std::cerr << "Reads file (fastq) must be specified." << std::endl;
        // exit(1);
    }

    return std::make_tuple(opt, map_param, index_parameters);
}
