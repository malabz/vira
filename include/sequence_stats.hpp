/**
 * @file sequence_stats.hpp
 * @brief Lightweight sequence statistics used by input quality control.
 */
#pragma once
#include "common.hpp"

#include <cstddef>

/** @brief Length and ambiguity metrics computed for one query record. */
struct SequenceStats {
    size_t seq_len = 0;          ///< Number of characters in the concatenated sequence.
    double n_ratio = 0.0;        ///< Fraction of characters equal to N/n.
    double len_diff_ratio = 0.0; ///< Signed `(seq_len - ref_len) / ref_len` value.
};

/**
 * @brief Computes sequence length, N content, and relative length difference.
 * @param rec Query record to inspect.
 * @param ref_len Reference-sequence length.
 * @return Calculated metrics for @p rec.
 * @note For an empty query, N ratio is 1.0 and length difference is -1.0 when
 *       @p ref_len is nonzero. If @p ref_len is zero, length difference is 0.0.
 */
SequenceStats calculate_sequence_stats(const FastaRecord& rec, size_t ref_len);
