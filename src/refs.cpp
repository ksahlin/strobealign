#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <regex>
#include <string_view>
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

void check_no_duplicates(const std::vector<std::string>& names) {
    std::vector<std::string_view> names_view{names.begin(), names.end()};
    std::sort(names_view.begin(), names_view.end());
    auto it = std::adjacent_find(names_view.begin(), names_view.end());
    if (it != names_view.end()) {
        std::ostringstream oss;
        oss << "FASTA file has duplicate reference sequence name '"
            << *it << "'";
        throw InvalidFasta(oss.str().c_str());
    }
}


// SAM spec says that reference names must match the following regex:
// [0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*

bool _is_valid_first_char(char c) {
    return (c >= '0' && c <= '9') ||
            (c >= 'A' && c <= 'Z') ||
            (c >= 'a' && c <= 'z') ||
            c == '!' || c == '#' || c == '$' || c == '%' ||
            c == '&' || c == '+' || c == '.' || c == '/' ||
            c == ':' || c == ';' || c == '?' || c == '@' ||
            c == '^' || c == '_' || c == '|' || c == '~' ||
            c == '-';
}

bool _is_valid_char(char c) {
    return _is_valid_first_char(c) || c == '*' || c == '=';
}

bool is_valid_sam_name(const std::string& name) {
    // From SAM specification 1.2.1
    if (name.empty()) {
        return false;
    }
    if (!_is_valid_first_char(name[0])) {
        return false;
    }
    for (size_t i = 1; i < name.size(); ++i) {
        if (!_is_valid_char(name[i])) {
            return false;
        }
    }
    return true;
}

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
    size_t line_num = 0;
    bool eof = false;
    do {
        line_num += 1;
        eof = !bool{getline(stream, line)};
        if (eof || (!line.empty() && line[0] == '>')) {
            if (seq.length() > 0) {
                to_uppercase(seq);
                sequences.push_back(std::move(seq));
                seq.clear();
                names.push_back(std::move(name));
            }
            if (!eof) {
                // Cut at the first whitespace
                std::string::size_type space = line.find_first_of(" \t\f\v\n\r");
                name = space < line.length() ? line.substr(1, space - 1) : line.substr(1);

                // Check the name is valid for the SAM output
                if (!is_valid_sam_name(name)) {
                    std::ostringstream oss;
                    oss << "FASTA file has SAM-incompatible reference sequence name '"
                        << name << "' on line " << line_num;
                    throw InvalidFasta(oss.str().c_str());
                }
            }
        } else {
            seq += line;
        }
    } while (!eof);

    check_no_duplicates(names);
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
