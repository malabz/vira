/**
 * @file qc.hpp
 * @brief MAFFT mapout parsing and post-alignment quality classification.
 */
#pragma once
#include "common.hpp"

#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

/** @brief Reports malformed or structurally invalid MAFFT mapout content. */
class MapoutFormatError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/** @brief Reports failures while opening or reading a MAFFT mapout file. */
class MapoutIoError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/**
 * @brief Parses MAFFT mapout data into deletion blocks and per-sample totals.
 * @param map_path Path to a MAFFT `.map` file.
 * @param expected_seq_ids Exact primary-identifier set dispatched in the chunk.
 * @return Pair containing all contiguous lost intervals and statistics keyed by
 *         sequence identifier.
 * @throws MapoutIoError If the file cannot be opened or a read operation fails.
 * @throws MapoutFormatError If a section identifier, field count, or numeric
 *         field is malformed, expected identifiers are not unique, or the
 *         parsed primary-identifier set differs from @p expected_seq_ids.
 *         Blank and comment lines are accepted; section order is irrelevant.
 */
std::pair<std::vector<DeletedBlock>, std::unordered_map<std::string, QcStats>> parse_map_to_deleted_blocks(
    const std::string& map_path,
    const std::vector<std::string>& expected_seq_ids
);

/**
 * @brief Classifies post-alignment loss metrics using configured thresholds.
 * @param lost_ratio Fraction of original sequence positions lost in alignment.
 * @param max_lost_block_len Longest contiguous lost interval.
 * @param args Validated WARN and FAIL thresholds.
 * @return `FAIL` when either fail threshold is met, otherwise `WARN` when either
 *         warning threshold is met, otherwise `PASS`.
 * @note Threshold comparisons are inclusive and FAIL takes precedence over WARN.
 */
std::string assign_qc_level(double lost_ratio, long long max_lost_block_len, const Args& args);

/**
 * @brief Combines input metadata and mapout statistics into sample QC rows.
 * @param sample_meta Ordered metadata keyed by effective sequence identifier.
 * @param sample_qc_stats Parsed loss statistics keyed by sequence identifier.
 * @param args Validated post-alignment QC thresholds.
 * @return One row per metadata entry that also has validated mapout statistics,
 *         preserving metadata map-key order.
 * @note Missing statistics are omitted rather than converted to zero-loss PASS.
 *       Empty sequences with present statistics receive a lost ratio of zero
 *       to avoid division by zero.
 */
std::vector<SampleRow> build_sample_qc_rows(
    const std::map<std::string, SampleMeta>& sample_meta,
    const std::unordered_map<std::string, QcStats>& sample_qc_stats,
    const Args& args
);
