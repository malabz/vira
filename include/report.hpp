/**
 * @file report.hpp
 * @brief Writers for alignment, QC, diagnostic, and filtering artifacts.
 */
#pragma once
#include "common.hpp"

#include <string>
#include <unordered_map>
#include <vector>

/**
 * @brief Writes every parsed deletion interval as a TSV report.
 * @param rows Rows to serialize in their current order.
 * @param path Destination path, truncated if it exists.
 */
void write_deleted_blocks_tsv(const std::vector<DeletedBlock>& rows, const std::string& path);

/**
 * @brief Writes one post-alignment QC summary row per dispatched sample.
 * @param rows Rows to serialize in their current order.
 * @param path Destination path, truncated if it exists.
 */
void write_sample_qc_summary_tsv(const std::vector<SampleRow>& rows, const std::string& path);

/**
 * @brief Writes deletion intervals belonging to WARN or FAIL samples.
 * @param block_rows Parsed deletion intervals.
 * @param sample_row_map Post-alignment result keyed by sequence identifier.
 * @param path Destination path, truncated if it exists.
 * @note Blocks without a matching sample row and blocks from PASS samples are omitted.
 */
void write_abnormal_deleted_blocks_tsv(
    const std::vector<DeletedBlock>& block_rows,
    const std::unordered_map<std::string, SampleRow>& sample_row_map,
    const std::string& path
);

/**
 * @brief Creates an MSA that excludes samples classified as FAIL.
 * @param full_msa_path Complete MSA whose first record is the shared reference.
 * @param filtered_msa_path Destination FASTA path, truncated if it exists.
 * @param sample_row_map Post-alignment QC rows keyed by effective identifier.
 * @throws std::runtime_error If the input or output file cannot be opened.
 * @note The first FASTA record is always retained. Later records are retained
 *       only when they have a matching QC row whose level is not FAIL.
 */
void make_filtered_msa(
    const std::string& full_msa_path,
    const std::string& filtered_msa_path,
    const std::unordered_map<std::string, SampleRow>& sample_row_map
);

/**
 * @brief Writes the complete input-QC audit joined with optional post-QC data.
 * @param rows Input records in original processing order.
 * @param post_qc Post-alignment rows keyed by effective identifier; absent rows
 *        are emitted with empty metrics and `NA` QC level.
 * @param path Destination path, truncated if it exists.
 */
void write_qc_summary_tsv(
    const std::vector<InputQcRow>& rows,
    const std::unordered_map<std::string, SampleRow>& post_qc,
    const std::string& path
);

/**
 * @brief Writes records excluded or withheld from dispatch by input QC.
 * @param rows Input-QC rows in processing order.
 * @param path Destination path, truncated if it exists.
 * @note Both `EXCLUDE` and `REVIEW_NOT_ALIGNED` decisions are included.
 */
void write_excluded_samples_tsv(const std::vector<InputQcRow>& rows, const std::string& path);

/**
 * @brief Writes every record whose final alignment status is not `ALIGNED`.
 * @param rows Input-QC rows after chunk completion statuses are applied.
 * @param path Destination path, truncated if it exists.
 */
void write_unaligned_samples_tsv(const std::vector<InputQcRow>& rows, const std::string& path);

/**
 * @brief Writes records with an abnormal input, alignment, or post-QC status.
 * @param rows Input-QC rows after chunk completion statuses are applied.
 * @param post_qc Post-alignment rows keyed by effective identifier.
 * @param path Destination path, truncated if it exists.
 * @note A record is abnormal when its decision is not `ALIGN`, its alignment
 *       status is not `ALIGNED`, or its post-QC level is WARN/FAIL.
 */
void write_abnormal_samples_tsv(
    const std::vector<InputQcRow>& rows,
    const std::unordered_map<std::string, SampleRow>& post_qc,
    const std::string& path
);

/**
 * @brief Writes diagnostics for failed and timed-out chunks.
 * @param rows Failed chunk results in collection order.
 * @param path Destination path, truncated if it exists.
 */
void write_failed_chunks_tsv(const std::vector<ChunkResult>& rows, const std::string& path);

/**
 * @brief Writes timing and status information for chunks marked as slow.
 * @param rows Slow chunk results in collection order.
 * @param path Destination path, truncated if it exists.
 */
void write_slow_chunks_tsv(const std::vector<ChunkResult>& rows, const std::string& path);
