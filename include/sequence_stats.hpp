#pragma once
#include "common.hpp"

#include <cstddef>

struct SequenceStats {
    size_t seq_len = 0;
    double n_ratio = 0.0;
    double len_diff_ratio = 0.0;
};

SequenceStats calculate_sequence_stats(const FastaRecord& rec, size_t ref_len);
