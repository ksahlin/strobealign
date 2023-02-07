#include "refs.hpp"
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <zlib.h>

/* Convert string to uppercase in-place */
void to_uppercase(std::string& s) {
        std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c) {
            return c & ~32;
        });
    }

namespace {

unsigned char to_uppercase1(unsigned char c) {
    return c & ~32;
}
class GZFile {
public:
    GZFile(const std::string& filename, const char* mode) {
        file_ = gzopen(filename.c_str(), mode);
        if (!file_) {
            throw InvalidFasta("Could not open file " + filename);
        }
    }
    ~GZFile() {
        gzclose(file_);
    }
    int getc() {
        return gzgetc(file_);
    }

private:
    gzFile file_;
};

enum struct State {
    Start,
    Header,
    HeaderRest,
    Sequence,
    EOL
};

}

References References::from_fasta(const std::string& filename) {
    std::vector<std::string> sequences;
    ref_names names;

    GZFile gf(filename, "r");

    std::string name, seq;
    auto state = State::Start;
    while (true) {
        const int c = gf.getc();
        if (c == EOF) {
            if (seq != "") {
                names.push_back(std::move(name));
                sequences.push_back(std::move(seq));
            }
            break;
        }
        // If we see a newline, we move to EOL state and ignore the input
        // character (except if it happens to be the very first character in
        // the file).
        if (c == '\n' && state != State::Start) {
            state = State::EOL;
            continue;
        }
        switch (state) {
            case State::Start:
                if (c == '>') {
                    state = State::Header;
                } else {
                    std::ostringstream oss;
                    oss << "FASTA file ('" << filename << "') must begin with '>' character, not '"
                        << static_cast<unsigned char>(c) << "'";
                    throw InvalidFasta(oss.str());
                }
                break;
            case State::Header:
                if (c == ' ') {
                    state = State::HeaderRest;
                } else {
                    name += static_cast<char>(c);
                }
                break;
            case State::HeaderRest:
                // Ignore everything until the end of the line
                break;
            case State::Sequence:
                seq += to_uppercase1(static_cast<char>(c));
                break;
            case State::EOL:
                if (c == '>') {
                    state = State::Header;
                    if (seq != "") {
                        names.push_back(std::move(name));
                        sequences.push_back(std::move(seq));
                    }

                    name.clear();
                    seq.clear();
                } else {
                    seq += to_uppercase1(static_cast<char>(c));
                    state = State::Sequence;
                }
                break;
            }
    }

    return References(std::move(sequences), std::move(names));
}

void References::add(std::string&& name, std::string&& sequence) {
    names.push_back(name);
    sequences.push_back(sequence);
    lengths.push_back(sequence.size());
}
