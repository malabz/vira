#include "sequence_stats.hpp"

SequenceStats calculate_sequence_stats(const FastaRecord& rec, size_t ref_len) {
    SequenceStats stats;
    stats.seq_len = rec.seq.size();

    if (stats.seq_len == 0) {
        stats.n_ratio = 1.0;
        stats.len_diff_ratio = ref_len > 0 ? -1.0 : 0.0;
        return stats;
    }

    size_t n_count = 0;
    for (char c : rec.seq) {
        if (c == 'N' || c == 'n') ++n_count;
    }

    stats.n_ratio = static_cast<double>(n_count) / static_cast<double>(stats.seq_len);
    stats.len_diff_ratio = ref_len > 0
        ? (static_cast<double>(stats.seq_len) - static_cast<double>(ref_len)) / static_cast<double>(ref_len)
        : 0.0;
    return stats;
}
