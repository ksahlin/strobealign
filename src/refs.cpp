#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
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

    std::string line, seq, name;
    bool eof = false;
    do {
        eof = !bool{getline(stream, line)};
        if (eof || (!line.empty() && line[0] == '>')) {
            if (seq.length() > 0) {
                to_uppercase(seq);
                sequences.push_back(seq);
                names.push_back(name);
            }
            if (!eof) {
                name = line.substr(1, line.find(' ') - 1); // cut at first space
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
