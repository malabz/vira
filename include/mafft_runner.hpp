/**
 * @file mafft_runner.hpp
 * @brief Execution boundary for one MAFFT query chunk.
 */
#pragma once
#include "common.hpp"

#include <string>
#include <vector>

/**
 * @brief Aligns one query chunk to the shared reference with MAFFT.
 *
 * The function captures stdout/stderr, enforces an optional timeout, optionally
 * parses MAFFT mapout QC data, and removes the per-chunk query and map files
 * before returning.
 *
 * @param chunk_id Zero-based identifier within @p bucket_name.
 * @param bucket_name Logical dispatch bucket used in labels and reports.
 * @param chunk_query_path Temporary FASTA containing the query records.
 * @param ref_path FASTA path containing the single reference record.
 * @param threads_per_mafft Positive thread count passed to MAFFT.
 * @param temp_dir Existing directory used for optional failure diagnostics.
 * @param save_failed_stderr Persist stderr when MAFFT fails or throws.
 * @param enable_qc Request `--mapout` and parse the resulting map file.
 * @param mafft_timeout_sec Timeout in seconds; zero disables timeout handling.
 * @param slow_chunk_sec Inclusive elapsed-time threshold; non-positive disables slow marking.
 * @param seq_ids Effective query identifiers in chunk order; moved into the result.
 * @param total_bases Sum of query-sequence lengths in the chunk.
 * @return A ChunkResult describing success, timing, output, diagnostics, and
 *         optional post-alignment QC data.
 * @note Process-launch exceptions are converted into a failed result rather
 *       than propagated. The caller must ensure `mafft` is available on PATH.
 */
ChunkResult run_mafft_chunk(
    int chunk_id,
    const std::string& bucket_name,
    const std::string& chunk_query_path,
    const std::string& ref_path,
    int threads_per_mafft,
    const std::string& temp_dir,
    bool save_failed_stderr,
    bool enable_qc,
    int mafft_timeout_sec,
    double slow_chunk_sec,
    std::vector<std::string> seq_ids,
    size_t total_bases
);
