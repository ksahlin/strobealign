#include "refs.hpp"
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "exceptions.hpp"

References References::from_fasta(const std::string& filename) {
    std::vector<std::string> sequences;
    ref_names names;
    ref_lengths lengths;

    std::ifstream file(filename);
    std::string line, seq;

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

    while (getline(file, line)) {
        if (line[0] == '>') {
            if (seq.length() > 0) {
                sequences.push_back(seq);
                lengths.push_back(seq.length());
            }
            names.push_back(line.substr(1, line.find(' ') - 1)); // cut at first space
            seq = "";
        }
        else {
            seq += line;
        }
    }
    if (seq.length() > 0){
        std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

        sequences.push_back(seq);
        lengths.push_back(seq.length());
    }

    return References(std::move(sequences), std::move(names), std::move(lengths));
}
