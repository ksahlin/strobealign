/*
 * Generate syncmers or randstrobes for a reference and output them in BED format
 */
#include <iostream>
#include <args.hxx>
#include "refs.hpp"
#include "arguments.hpp"
#include "randstrobes.hpp"
#include "index.hpp"
#include "logger.hpp"

static Logger& logger = Logger::get();

struct BedRecord {
    const std::string& contig;
    size_t start;
    size_t end;
};

std::ostream& operator<<(std::ostream& os, const BedRecord& record) {
    os << record.contig << '\t' << record.start << '\t' << record.end << '\n';
    return os;
}

void dump_randstrobes(std::ostream& os, const std::string& name, const std::string& sequence, const IndexParameters& parameters) {
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_coordinate;
    std::tie(string_hashes, pos_to_seq_coordinate) = make_string_to_hashvalues_open_syncmers_canonical(
        sequence, parameters.k, parameters.s, parameters.t_syncmer);

    unsigned int nr_hashes = string_hashes.size();
    if (nr_hashes == 0) {
        return;
    }

    RandstrobeIterator randstrobe_iter { string_hashes, pos_to_seq_coordinate, parameters.w_min, parameters.w_max, parameters.q, parameters.max_dist };
    while (randstrobe_iter.has_next()) {
        auto randstrobe = randstrobe_iter.next();
        os << BedRecord{name, randstrobe.strobe1_pos, randstrobe.strobe2_pos + parameters.k};
    }
}

void dump_randstrobes2(std::ostream& os, const std::string& name, const std::string& sequence, const IndexParameters& parameters) {
    auto randstrobe_iter = RandstrobeIterator2(
        sequence, parameters.k, parameters.s, parameters.t_syncmer, parameters.w_min, parameters.w_max, parameters.q, parameters.max_dist);
    Randstrobe randstrobe;
    while ((randstrobe = randstrobe_iter.next()) != randstrobe_iter.end()) {
        os << BedRecord{name, randstrobe.strobe1_pos, randstrobe.strobe2_pos + parameters.k};
    }
}

void dump_syncmers(std::ostream& os, const std::string& name, const std::string& sequence, const IndexParameters& parameters) {
    SyncmerIterator syncmer_iterator(sequence, parameters.k, parameters.s, parameters.t_syncmer);
    Syncmer syncmer;
    while (!(syncmer = syncmer_iterator.next()).is_end()) {
        os << BedRecord{name, syncmer.position, syncmer.position + parameters.k};
    }
}

int run_dumpstrobes(int argc, char **argv) {
    args::ArgumentParser parser("dumpstrobes");
    parser.helpParams.showTerminator = false;
    parser.helpParams.helpindent = 20;
    parser.helpParams.width = 90;
    parser.helpParams.programName = "dumpstrobes";
    parser.helpParams.shortSeparator = " ";

    args::HelpFlag help(parser, "help", "Print help and exit", {'h', "help"});
    args::Flag syncmers(parser, "syncmers", "Dump syncmers instead of randstrobes", {"syncmers"});
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
    int r{150}, k{20}, s{16}, c{8}, l{1}, u{7};
    int max_seed_len{};

    bool k_set{false}, s_set{false}, c_set{false}, max_seed_len_set{false}, l_set{false}, u_set{false};
    if (seeding.r) { r = args::get(seeding.r); }
    if (seeding.m) { max_seed_len = args::get(seeding.m); max_seed_len_set = true; }
    if (seeding.k) { k = args::get(seeding.k); k_set = true; }
    if (seeding.l) { l = args::get(seeding.l); l_set = true; }
    if (seeding.u) { u = args::get(seeding.u); u_set = true; }
    if (seeding.s) { s = args::get(seeding.s); s_set = true; }
    if (seeding.c) { c = args::get(seeding.c); c_set = true; }

    // Reference
    auto ref_path = args::get(ref_filename);
    if (c >= 64 || c <= 0) {
        throw BadParameter("c must be greater than 0 and less than 64");
    }
    IndexParameters index_parameters = IndexParameters::from_read_length(
        r,
        c_set ? c : IndexParameters::DEFAULT,
        k_set ? k : IndexParameters::DEFAULT,
        s_set ? s : IndexParameters::DEFAULT,
        l_set ? l : IndexParameters::DEFAULT,
        u_set ? u : IndexParameters::DEFAULT,
        max_seed_len_set ? max_seed_len : IndexParameters::DEFAULT
    );

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

    for (size_t i = 0; i < references.size(); ++i) {
        auto& seq = references.sequences[i];
        auto& name = references.names[i];
        if (syncmers) {
            dump_syncmers(std::cout, name, seq, index_parameters);
        } else {
            dump_randstrobes2(std::cout, name, seq, index_parameters);
        }
    }

    return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
    try {
        return run_dumpstrobes(argc, argv);
    } catch (BadParameter& e) {
        logger.error() << "A seeding parameter is invalid: " << e.what() << std::endl;
    } catch (const std::runtime_error& e) {
        logger.error() << "dumpstrobes: " << e.what() << std::endl;
    }
    return EXIT_FAILURE;
}
