#ifndef SAM_HPP
#define SAM_HPP

#include <string>

#include "index.hpp"

struct alignment {
    std::string cigar;
    int ref_start;
    int ed;
    int global_ed;
    int sw_score;
    int aln_score;
    int ref_id;
    int mapq;
    int aln_length;
    bool not_proper;
    bool is_rc;
    bool is_unaligned = false;
};

enum SamFlags {
    PAIRED = 1,
    PROPER_PAIR = 2,
    UNMAP = 4,
    MUNMAP = 8,
    REVERSE = 0x10,
    MREVERSE = 0x20,
    READ1 = 0x40,
    READ2 = 0x80,
    SECONDARY = 0x100,
    QCFAIL = 0x200,
    DUP = 0x400,
    SUPPLEMENTARY = 0x800,
};

class Sam {

public:
    Sam(std::string& sam_string, idx_to_acc& acc_map) : sam_string(sam_string), acc_map(acc_map) { }

    /* Add an alignment */
    void add(const alignment& sam_aln, const std::string& sequence, const std::string& sequence_rc, const std::string& query_acc, const std::string& qual);

    /* Add an unmapped read */
    void unmapped(std::string& name, std::string& sequence, std::string& qualities);

private:
    std::string& sam_string;
    idx_to_acc& acc_map;
};


#endif
