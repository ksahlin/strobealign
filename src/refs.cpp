#include "refs.hpp"
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "exceptions.hpp"

References References::from_fasta(const std::string& filename) {
    std::vector<std::string> sequences;
    ref_names names;

    std::ifstream file(filename);

    if (!file.good()) {
        throw InvalidFasta("Cannot read from FASTA file");
    }

    auto c = file.peek();
    if (c != '>') {
        std::ostringstream oss;
        oss << "FASTA file must begin with '>' character, not '"
            << static_cast<unsigned char>(c) << "'";
        throw InvalidFasta(oss.str().c_str());
    }

    std::string line, seq, name;
    bool eof = false;
    do {
        eof = !bool{getline(file, line)};
        if (eof || (!line.empty() && line[0] == '>')) {
            if (seq.length() > 0) {
                std::transform(seq.begin(), seq.end(), seq.begin(),
                    [](unsigned char c) {
                        return c & ~32;  // convert to uppercase
                    }
                );
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

void References::add(std::string&& name, std::string&& sequence) {
    names.push_back(name);
    sequences.push_back(sequence);
    lengths.push_back(sequence.size());
}
