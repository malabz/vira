/**
 * @file input_qc.hpp
 * @brief Input normalization and pre-alignment quality-control interfaces.
 */
#pragma once
#include "common.hpp"

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * @brief Splits a comma-separated keyword list into normalized search terms.
 * @param csv Comma-separated text; surrounding whitespace is ignored.
 * @return Non-empty terms converted to lowercase, in input order.
 */
std::vector<std::string> split_keywords(const std::string& csv);

/**
 * @brief Returns a lowercase copy using byte-wise character conversion.
 * @param s Input text.
 * @return Lowercase copy of @p s.
 */
std::string lowercase_copy(const std::string& s);

/**
 * @brief Replaces tabs and line terminators with spaces for safe TSV fields.
 * @param s Unescaped field text.
 * @return Sanitized copy suitable for one TSV cell.
 */
std::string sanitize_tsv_field(const std::string& s);

/**
 * @brief Joins stable QC reason codes with semicolons.
 * @param reasons Ordered reason codes.
 * @return Semicolon-separated codes, or an empty string for an empty vector.
 */
std::string join_reason_codes(const std::vector<std::string>& reasons);

/** @brief Record, report row, and dispatch decision produced by input QC. */
struct InputQcDecision {
    FastaRecord record;   ///< Effective record, including any duplicate-ID rename.
    InputQcRow row;       ///< Audit row describing the decision and its reasons.
    bool dispatch = false; ///< Whether the record should enter the MAFFT queue.
};

/**
 * @brief Applies duplicate, content, length, and header-based input QC rules.
 *
 * Duplicate counts in @p seen_ids are updated as part of the assessment. The
 * returned record may be renamed when the duplicate policy is `rename`.
 *
 * @param rec Record to assess; moved into the returned decision.
 * @param args Validated thresholds and policy settings.
 * @param ref_len Reference-sequence length used for length-ratio rules.
 * @param bucket Dispatch bucket recorded in reports.
 * @param seq_len Precomputed length of @p rec.
 * @param n_ratio Precomputed N/n fraction, expected in `[0, 1]`.
 * @param len_diff_ratio Precomputed signed length difference relative to the reference.
 * @param seen_ids Mutable occurrence count keyed by original record identifier.
 * @param non_target_keywords Lowercase header substrings indicating a non-target sample.
 * @param engineered_keywords Lowercase header substrings indicating an engineered construct.
 * @return Effective record, complete QC row, and dispatch flag.
 * @note An empty sequence, illegal character, extreme length, or excluded
 *       duplicate produces `EXCLUDE`. Suspected non-target or engineered
 *       records are not dispatched unless review alignment is enabled.
 */
InputQcDecision assess_input_record(
    FastaRecord rec,
    const Args& args,
    size_t ref_len,
    const std::string& bucket,
    size_t seq_len,
    double n_ratio,
    double len_diff_ratio,
    std::unordered_map<std::string, size_t>& seen_ids,
    const std::vector<std::string>& non_target_keywords,
    const std::vector<std::string>& engineered_keywords
);
