#include "readlen.hpp"

int est_read_length(klibpp::KStream<gzFile_s*, int (*)(gzFile_s*, void*, unsigned int), klibpp::mode::In_> &ks, int n_reads) {
    auto records = ks.read(n_reads);
    int tot_read_len = 0;
    for (size_t i = 0; i < records.size(); ++i) {
        auto record1 = records[i];
        tot_read_len += record1.seq.length();
    }
    int avg_read_len = tot_read_len/n_reads;

    return avg_read_len;
}

/*
 * Return average read length of single-end or paired-end reads.
 * Set filename2 to the empty string if data is single end.
 */
int estimate_read_length(const std::string& filename1, const std::string& filename2) {
    bool is_paired = filename2 != "";

    gzFile fp1_tmp = gzopen(filename1.c_str(), "r");
    auto ks1_tmp = klibpp::make_ikstream(fp1_tmp, gzread);
    auto r1_tmp = est_read_length(ks1_tmp, 500);
    gzclose(fp1_tmp);

    if (is_paired) {
        gzFile fp2_tmp = gzopen(filename2.c_str(), "r");
        auto ks2_tmp = klibpp::make_ikstream(fp2_tmp, gzread);
        auto r2_tmp = est_read_length(ks2_tmp, 500);
        gzclose(fp2_tmp);
        return (r1_tmp + r2_tmp) / 2;
    } else {
        return r1_tmp;
    }
}

