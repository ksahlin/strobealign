#include "refs.hpp"
#include <vector>
#include <fstream>
#include <sstream>
#include "exceptions.hpp"

References References::from_fasta(std::string filename) {
    std::vector<std::string> sequences;
    std::vector<std::string> names;
    std::vector<unsigned int> lengths;

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
        sequences.push_back(seq);
        lengths.push_back(seq.length());
    }

    return References(std::move(sequences), std::move(names), std::move(lengths));
}
