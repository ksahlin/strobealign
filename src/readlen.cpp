#include "readlen.hpp"

int total_read_len(const std::vector<klibpp::KSeq> &records) {

    int tot_read_len = 0;
    for (size_t i = 0; i < records.size(); ++i) {
        tot_read_len += records[i].seq.length();
    }
    return tot_read_len;
}

/*
 * Return average read length of single-end or paired-end reads.
 * Set filename2 to the empty string if data is single end.
 */
int estimate_read_length(InputBuffer& input_buffer) {
    std::vector<klibpp::KSeq> records1;
    std::vector<klibpp::KSeq> records2;
    std::vector<klibpp::KSeq> records3;
    AlignmentStatistics stats;
    input_buffer.read_records(records1, records2, records3, stats, 500);
    if (records1.empty() && records3.empty()) {
        return 150;
    }
    auto tot_read_len = total_read_len(records1)
                    + total_read_len(records2)
                    + total_read_len(records3);
    auto tot_read_num = records1.size() + records2.size() + records3.size();
    return tot_read_len / tot_read_num;
}

