/**
 * @file common.hpp
 * @brief Shared configuration and data-transfer types for the alignment pipeline.
 */
#pragma once

#include <chrono>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

/** @brief Monotonic clock used for chunk and whole-run timing. */
using Clock = std::chrono::steady_clock;

/**
 * @brief Validated runtime configuration assembled from command-line options.
 *
 * Default member values define the CLI defaults. Ratio fields are expressed as
 * fractions rather than percentages. Input-QC ratios are range-checked by
 * parse_args(); post-alignment warning and failure thresholds are used as supplied.
 */
struct Args {
    std::string ref;                    ///< Reference FASTA path; exactly one sequence is required.
    std::string query;                  ///< Query FASTA path.
    std::string out;                    ///< Destination path for the complete MSA.
    int max_instances = 16;             ///< Concurrency basis for adaptive chunking and the pending-job limit.
    int threads_per_mafft = 1;          ///< Thread count passed to each MAFFT process.
    int chunk_size = 1000;              ///< Upper bound on records placed in one chunk.
    std::string temp_dir;               ///< User-supplied base path for the temporary directory.
    bool temp_dir_provided = false;     ///< Whether @ref temp_dir came from the command line.
    bool cleanup = false;               ///< Remove the per-run temporary directory after normal completion.
    bool save_failed_stderr = false;    ///< Persist stderr for failed chunks in the temporary directory.
    bool disable_auto_chunk = false;    ///< Use @ref chunk_size without adaptive reduction.
    int max_pending_factor = 2;         ///< Pending-job limit multiplier relative to @ref max_instances.
    bool enable_qc = false;             ///< Request MAFFT mapout parsing and post-alignment QC reports.
    double warn_lost_ratio = 0.001;     ///< Inclusive lost-base ratio threshold for WARN.
    double fail_lost_ratio = 0.01;      ///< Inclusive lost-base ratio threshold for FAIL.
    int warn_max_block_len = 20;        ///< Inclusive maximum-deletion-block threshold for WARN.
    int fail_max_block_len = 100;       ///< Inclusive maximum-deletion-block threshold for FAIL.

    /** @name Layer 1: input quality control */
    ///@{
    std::string duplicate_id_policy = "exclude"; ///< Duplicate handling mode: `exclude` or `rename`.
    std::string allowed_seq_chars = "ACGTURYSWKMBDHVN-."; ///< Case-insensitive sequence alphabet.
    double review_n_ratio = 0.05;        ///< Inclusive N-ratio threshold for manual review.
    double exclude_n_ratio = 0.50;       ///< Inclusive N-ratio threshold for exclusion.
    double exclude_min_len_ratio = 0.10; ///< Length-ratio boundary below which records are excluded.
    double exclude_max_len_ratio = 5.00; ///< Length-ratio boundary above which records are excluded.
    ///@}

    /** @name Layer 2: virology suitability */
    ///@{
    double reasonable_min_len_ratio = 0.85; ///< Lower ratio below which an otherwise valid sample is reviewed.
    double reasonable_max_len_ratio = 1.15; ///< Upper ratio above which an otherwise valid sample is reviewed.
    std::string non_target_keywords = "bacillus,bacteria,bacterium,chromosome,escherichia,salmonella,staphylococcus"; ///< Comma-separated header terms suggesting a non-target sequence.
    std::string engineered_keywords = "mutant,recombinant,recomb,egfp,gfp,nluc,luciferase,reporter,construct,vector,plasmid,clone"; ///< Comma-separated header terms suggesting engineering.
    bool align_review_samples = false; ///< Dispatch suspected non-target or engineered samples when true.
    ///@}

    /** @name Layer 3: unified load-aware dispatch */
    ///@{
    size_t target_chunk_bases = 0; ///< Maximum bases per chunk; zero selects `effective_chunk_size * ref_len`.
    int mafft_timeout_sec = 0;     ///< Per-chunk timeout in seconds; zero disables the timeout.
    double slow_chunk_sec = 300.0; ///< Inclusive slow threshold in seconds; zero disables slow marking.
    ///@}
};

/**
 * @brief One FASTA record in parsed, unwrapped form.
 *
 * The header is stored without the leading `>` and the identifier is the first
 * whitespace-delimited token in that header. Sequence lines are concatenated.
 */
struct FastaRecord {
    std::string id;     ///< Primary identifier used to correlate results and reports.
    std::string header; ///< Complete FASTA header without the leading `>`.
    std::string seq;    ///< Concatenated sequence characters, without line terminators.
};

/** @brief Sequence attributes retained for post-alignment QC reporting. */
struct SampleMeta {
    std::string bucket;          ///< Dispatch bucket name; currently `all` for the unified queue.
    size_t seq_len = 0;          ///< Original query length in characters.
    double n_ratio = 0.0;        ///< Fraction of sequence characters equal to N/n.
    double len_diff_ratio = 0.0; ///< Signed `(seq_len - ref_len) / ref_len` value.
};

/**
 * @brief End-to-end status row for one input record.
 *
 * This structure is populated during input QC, updated when a chunk completes,
 * and later serialized into the input and abnormal-sample reports.
 */
struct InputQcRow {
    std::string original_id;       ///< Identifier before any duplicate-renaming policy is applied.
    std::string seq_id;            ///< Effective identifier used by alignment and reports.
    std::string header;            ///< Effective FASTA header after any duplicate rename.
    std::string bucket;            ///< Dispatch bucket assigned to the sample.
    size_t seq_len = 0;            ///< Original sequence length.
    double n_ratio = 0.0;          ///< Fraction of sequence characters equal to N/n.
    double len_diff_ratio = 0.0;   ///< Signed length difference relative to the reference.
    std::string decision;          ///< `ALIGN`, `ALIGN_REVIEW`, `EXCLUDE`, or `REVIEW_NOT_ALIGNED`.
    std::string reasons;           ///< Semicolon-separated input-QC reason codes; empty when none apply.
    bool needs_manual_review = false; ///< Whether the record requires human review.
    std::string alignment_status;  ///< `ALIGNED`, `PENDING`, `NOT_DISPATCHED`, `CHUNK_FAILED`, or `CHUNK_TIMEOUT`.
    std::string alignment_detail;  ///< Chunk label and, for failures, the failure reason.
};

/** @brief One contiguous query interval that MAFFT mapout marks as lost. */
struct DeletedBlock {
    std::string seq_id;              ///< Query identifier reported by mapout.
    std::string bucket;              ///< Dispatch bucket copied from the originating chunk.
    long long block_start_orig = 0;  ///< Original-sequence coordinate at the start of the interval.
    long long block_end_orig = 0;    ///< Original-sequence coordinate at the end of the interval.
    long long block_len = 0;         ///< Inclusive interval length.
    std::string left_ref_pos;        ///< Nearest mapped reference position to the left, or empty.
    std::string right_ref_pos;       ///< Nearest mapped reference position to the right, or empty.
};

/** @brief Aggregated mapout loss statistics for one sample. */
struct QcStats {
    long long lost_count = 0;         ///< Total number of original positions mapped to `-`.
    long long lost_block_count = 0;   ///< Number of contiguous lost intervals.
    long long max_lost_block_len = 0; ///< Length of the longest contiguous lost interval.
};

/** @brief Final post-alignment QC row for one dispatched sample. */
struct SampleRow {
    std::string seq_id;                 ///< Effective sample identifier.
    std::string bucket;                 ///< Dispatch bucket assigned to the sample.
    size_t seq_len = 0;                 ///< Original sequence length.
    double len_diff_ratio = 0.0;        ///< Signed length difference relative to the reference.
    double n_ratio = 0.0;               ///< Fraction of sequence characters equal to N/n.
    long long lost_count = 0;           ///< Number of original positions lost by alignment.
    double lost_ratio = 0.0;            ///< `lost_count / seq_len`, or zero for an empty sequence.
    long long lost_block_count = 0;     ///< Number of contiguous lost intervals.
    long long max_lost_block_len = 0;   ///< Longest contiguous lost interval.
    std::string qc_level;               ///< Post-alignment level: `PASS`, `WARN`, or `FAIL`.
};

/**
 * @brief Complete result of processing one independently dispatched chunk.
 *
 * A failed result carries status and diagnostics but no MSA payload. Successful
 * results may additionally contain parsed mapout details when QC was enabled.
 */
struct ChunkResult {
    std::string bucket_name; ///< Bucket from which the chunk was dispatched.
    int chunk_id = 0;        ///< Zero-based identifier unique within the bucket.
    std::string stdout_data; ///< Complete MAFFT stdout for a successful process.
    bool success = false;    ///< Whether MAFFT exited successfully before the timeout.
    bool timed_out = false;  ///< Whether timeout handling terminated the process group.
    bool slow = false;       ///< Whether elapsed time reached the slow-chunk threshold.
    double elapsed_sec = 0.0; ///< Wall-clock duration of the external process.
    std::string err_path;    ///< Saved stderr path, or empty when stderr was not persisted.
    std::string failure_reason; ///< Stable failure category used in reports.
    size_t record_count = 0; ///< Number of query records assigned to the chunk.
    size_t total_bases = 0;  ///< Sum of query-sequence lengths in the chunk.
    std::vector<std::string> seq_ids; ///< Effective identifiers in dispatch order.
    std::vector<DeletedBlock> deleted_blocks; ///< Parsed mapout deletion intervals.
    std::unordered_map<std::string, QcStats> qc_sample_stats; ///< Parsed per-sample mapout statistics.
};

