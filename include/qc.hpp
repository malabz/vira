/**
 * @file qc.hpp
 * @brief MAFFT mapout parsing and post-alignment quality classification.
 */
#pragma once
#include "common.hpp"

#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

/**
 * @brief Parses MAFFT mapout data into deletion blocks and per-sample totals.
 * @param map_path Path to a MAFFT `.map` file.
 * @return Pair containing all contiguous lost intervals and statistics keyed by
 *         sequence identifier.
 * @note An unreadable file returns two empty collections. Blank, comment, and
 *       malformed data lines are skipped rather than treated as fatal errors.
 */
std::pair<std::vector<DeletedBlock>, std::unordered_map<std::string, QcStats>> parse_map_to_deleted_blocks(
    const std::string& map_path
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
 * @return One row per entry in @p sample_meta, preserving its map-key order.
 * @note Missing mapout statistics are treated as zero loss. Empty sequences
 *       receive a lost ratio of zero to avoid division by zero.
 */
std::vector<SampleRow> build_sample_qc_rows(
    const std::map<std::string, SampleMeta>& sample_meta,
    const std::unordered_map<std::string, QcStats>& sample_qc_stats,
    const Args& args
);
