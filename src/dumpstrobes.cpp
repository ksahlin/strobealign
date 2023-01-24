/*
 * Generate randstrobes for a reference and output them in BED format
 */
#include <iostream>
#include <args.hxx>
#include "refs.hpp"
#include "arguments.hpp"
#include "randstrobes.hpp"
#include "index.hpp"
#include "logger.hpp"

static Logger& logger = Logger::get();

void dump_randstrobes(const References& references, const IndexParameters& parameters) {
    for (size_t i = 0; i < references.size(); ++i) {
        auto& seq = references.sequences[i];

        std::vector<uint64_t> string_hashes;
        std::vector<unsigned int> pos_to_seq_coordinate;
        std::tie(string_hashes, pos_to_seq_coordinate) = make_string_to_hashvalues_open_syncmers_canonical(seq, parameters.k, parameters.s, parameters.t_syncmer);

        unsigned int nr_hashes = string_hashes.size();
        if (nr_hashes == 0) {
            continue;
        }

        RandstrobeIterator randstrobe_iter { string_hashes, pos_to_seq_coordinate, parameters.w_min, parameters.w_max, parameters.q, parameters.max_dist };
        while (randstrobe_iter.has_next()) {
            auto randstrobe = randstrobe_iter.next();
            std::cout
                << references.names[i] << '\t'
                << randstrobe.strobe1_pos << '\t'
                << randstrobe.strobe2_pos + parameters.k << '\n';
        }
    }
}

void dump_randstrobes2(const References& references, const IndexParameters& parameters) {
    for (size_t i = 0; i < references.size(); ++i) {
        auto& seq = references.sequences[i];
        auto randstrobe_iter = RandstrobeIterator2(seq, parameters.k, parameters.s, parameters.t_syncmer, parameters.w_min, parameters.w_max, parameters.q, parameters.max_dist);
        Randstrobe randstrobe;
        while ((randstrobe = randstrobe_iter.next()) != randstrobe_iter.end()) {
            std::cout
                << references.names[i] << '\t'
                << randstrobe.strobe1_pos << '\t'
                << randstrobe.strobe2_pos + parameters.k << '\n';
        }
    }
}

int run_dumprandstrobes(int argc, char **argv) {
    args::ArgumentParser parser("dumprandstrobes");
    parser.helpParams.showTerminator = false;
    parser.helpParams.helpindent = 20;
    parser.helpParams.width = 90;
    parser.helpParams.programName = "dumprandstrobes";
    parser.helpParams.shortSeparator = " ";

    args::HelpFlag help(parser, "help", "Print help and exit", {'h', "help"});
    auto seeding = SeedingArguments{parser};
    args::Positional<std::string> ref_filename(parser, "reference", "Reference in FASTA format", args::Options::Required);

    try {
        parser.ParseCLI(argc, argv);
    }
    catch (const args::Help&) {
        std::cout << parser;
        exit(EXIT_SUCCESS);
    }
    catch (const args::Error& e) {
        std::cerr << parser;
        std::cerr << "Error: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    // Seeding
    int r{150}, k{20}, s{16}, c{8};
    int max_seed_len{};

    bool k_set{false}, s_set{false}, c_set{false}, max_seed_len_set{false};
    if (seeding.r) { r = args::get(seeding.r); }
    if (seeding.m) { max_seed_len = args::get(seeding.m); max_seed_len_set = true; }
    if (seeding.k) { k = args::get(seeding.k); k_set = true; }
    if (seeding.s) { s = args::get(seeding.s); s_set = true; }
    if (seeding.c) { c = args::get(seeding.c); c_set = true; }

    // Reference
    auto ref_path = args::get(ref_filename);
    if (c >= 64 || c <= 0) {
        throw BadParameter("c must be greater than 0 and less than 64");
    }
    IndexParameters index_parameters = IndexParameters::from_read_length(
        r, c_set ? c : -1, k_set ? k : -1, s_set ? s : -1, max_seed_len_set ? max_seed_len : -1);
    index_parameters.verify();

    logger.info() << index_parameters << '\n';
    logger.info() << "Reading reference ...\n";
    References references;
    references = References::from_fasta(ref_path);

    logger.info() << "Reference size: " << references.total_length() / 1E6 << " Mbp ("
        << references.size() << " contig" << (references.size() == 1 ? "" : "s")
        << "; largest: "
        << (*std::max_element(references.lengths.begin(), references.lengths.end()) / 1E6) << " Mbp)\n";
    if (references.total_length() == 0) {
        throw InvalidFasta("No reference sequences found");
    }

    dump_randstrobes2(references, index_parameters);

    return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
    try {
        return run_dumprandstrobes(argc, argv);
    } catch (BadParameter& e) {
        logger.error() << "A seeding parameter is invalid: " << e.what() << std::endl;
    } catch (const std::runtime_error& e) {
        logger.error() << "dumpstrobes: " << e.what() << std::endl;
    }
    return EXIT_FAILURE;
}
