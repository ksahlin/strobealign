#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <regex>
#include <set>
#include "refs.hpp"
#include "zstr.hpp"


/* Convert string to uppercase in-place */
void to_uppercase(std::string& s) {
        std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c) {
            return c & ~32;
        }
    );
}

namespace {

template <typename T>
References references_from_stream(T& stream) {
    std::vector<std::string> sequences;
    std::vector<std::string> names;
    std::set<std::string> seen_names;

    if (!stream.good()) {
        throw InvalidFasta("Cannot read from FASTA file");
    }

    auto c = stream.peek();
    if (c != '>') {
        std::ostringstream oss;
        oss << "FASTA file must begin with '>' character, not '"
            << static_cast<unsigned char>(c) << "'";
        throw InvalidFasta(oss.str().c_str());
    }

    const std::regex ws_re(R"(\s+)");
    // From SAM specification 1.2.1
    const std::regex rname_re(R"([0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*$)");

    std::string line, seq, name;
    size_t line_num = 0;
    bool eof = false;
    do {
        line_num += 1;
        eof = !bool{getline(stream, line)};
        if (eof || (!line.empty() && line[0] == '>')) {
            if (seq.length() > 0) {
                to_uppercase(seq);
                sequences.push_back(seq);
                names.push_back(name);
                seen_names.insert(name);
            }
            if (!eof) {
                // Cut at the first whitespace
                auto itr = std::sregex_token_iterator(line.begin() + 1, line.end(), ws_re, -1);
                if (itr == std::sregex_token_iterator()) {
                    std::ostringstream oss;
                    oss << "FASTA file has invalid reference sequence name on line "
                        << line_num;
                    throw InvalidFasta(oss.str().c_str());
                }
                name = *itr;

                // Check the name is valid for the SAM output
                std::smatch m;
                if (!std::regex_match(name, m, rname_re)) {
                    std::ostringstream oss;
                    oss << "FASTA file has SAM-incompatible reference sequence name '"
                        << name << "' on line " << line_num;
                    throw InvalidFasta(oss.str().c_str());
                }

                // Check for duplicate names
                if (seen_names.count(name)) {
                    std::ostringstream oss;
                    oss << "FASTA file has duplicate reference sequence name '"
                        << name << "' on line " << line_num;
                    throw InvalidFasta(oss.str().c_str());
                }
            }
            seq = "";
        } else {
            seq += line;
        }
    } while (!eof);

    return References(std::move(sequences), std::move(names));
}

}

/* Read compressed or uncompressed reference */
References References::from_fasta(const std::string& filename) {
    if (filename.length() > 3 && filename.substr(filename.length() - 3, 3) == ".gz") {
        zstr::ifstream ifs(filename);
        return references_from_stream(ifs);
    } else {
        std::ifstream ifs(filename);
        return references_from_stream(ifs);
    }
}

void References::add(std::string&& name, std::string&& sequence) {
    names.push_back(name);
    sequences.push_back(sequence);
    lengths.push_back(sequence.size());
}
