#include <vector>
#include <set>
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

/* Following the SAM specification, section 1.2.1,
   return true if a given character is allowed in a reference sequence name. */
bool is_rname_char(char c) {
    // The specification uses:
    //      [0-9A-Za-z!#$%&*+./:;=?@^_|~-]

    if ('0' <= c && c <= '9') {
        return true;
    }
    if ('A' <= c && c <= 'Z') {
        return true;
    }
    if ('a' <= c && c <= 'z') {
        return true;
    }
    static constexpr auto other = "!#$%&*+./:;=?@^_|~-";
    for (auto p = other; *p; ++p) {
        if (c == *p) {
            return true;
        }
    }
    return false;
}

/* Following the SAM specification, section 1.2.1,
   return true if the given character is allowed at the start of a reference sequence name */
bool is_rname_start_char(char c) {
    return is_rname_char(c) && c != '*' && c != '=';
}

/* Return the longest prefix of the given string that satisfies the SAM specification
   rules for sequence names. */
std::string get_rname_prefix(std::string::const_iterator begin, std::string::const_iterator end) {
    if (begin == end) {
        return std::string{};
    }
    if (!is_rname_start_char(*begin)) {
        return std::string{};
    }
    for (auto itr = begin + 1; itr != end; ++itr) {
        if (!is_rname_char(*itr)) {
            std::string res;
            res.insert(res.end(), begin, itr);
            return res;
        }
    }
    std::string res;
    res.insert(res.end(), begin, end);
    return res;
}

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

    std::string line, seq, name;
    bool eof = false;
    size_t line_num = 0;
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
                name = get_rname_prefix(line.begin() + 1, line.end());
                if (name.size() == 0) {
                    std::ostringstream oss;
                    oss << "Cannot extract a valid reference sequence name at line "
                        << line_num;
                    throw InvalidFasta(oss.str().c_str());
                }
                if (seen_names.count(name)) {
                    std::ostringstream oss;
                    oss << "Duplicate reference sequence name '"
                        << name << "' at line " << line_num;
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
