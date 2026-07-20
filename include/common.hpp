#pragma once

#include <chrono>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

using Clock = std::chrono::steady_clock;

struct Args {
    std::string ref;
    std::string query;
    std::string out;
    int max_instances = 16;
    int threads_per_mafft = 1;
    int chunk_size = 1000;
    std::string temp_dir;
    bool temp_dir_provided = false;
    bool cleanup = false;
    bool save_failed_stderr = false;
    bool disable_auto_chunk = false;
    int max_pending_factor = 2;
    bool enable_qc = false;
    double warn_lost_ratio = 0.001;
    double fail_lost_ratio = 0.01;
    int warn_max_block_len = 20;
    int fail_max_block_len = 100;

    // Layer 1: input QC
    std::string duplicate_id_policy = "exclude"; // exclude | rename
    std::string allowed_seq_chars = "ACGTURYSWKMBDHVN-.";
    double review_n_ratio = 0.05;
    double exclude_n_ratio = 0.50;
    double exclude_min_len_ratio = 0.10;
    double exclude_max_len_ratio = 5.00;

    // Layer 2: virology suitability
    double reasonable_min_len_ratio = 0.85;
    double reasonable_max_len_ratio = 1.15;
    std::string non_target_keywords = "bacillus,bacteria,bacterium,chromosome,escherichia,salmonella,staphylococcus";
    std::string engineered_keywords = "mutant,recombinant,recomb,egfp,gfp,nluc,luciferase,reporter,construct,vector,plasmid,clone";
    bool align_review_samples = false;

    // Layer 3: unified load-aware dispatch
    size_t target_chunk_bases = 0;      // 0 => auto = effective_chunk_size * ref_len
    int mafft_timeout_sec = 0;          // 0 => disabled
    double slow_chunk_sec = 300.0;
};

struct FastaRecord {
    std::string id;
    std::string header;
    std::string seq;
};

struct SampleMeta {
    std::string bucket;
    size_t seq_len = 0;
    double n_ratio = 0.0;
    double len_diff_ratio = 0.0;
};

struct InputQcRow {
    std::string original_id;
    std::string seq_id;
    std::string header;
    std::string bucket;
    size_t seq_len = 0;
    double n_ratio = 0.0;
    double len_diff_ratio = 0.0;
    std::string decision;          // ALIGN | ALIGN_REVIEW | EXCLUDE | REVIEW_NOT_ALIGNED
    std::string reasons;           // semicolon-separated reason codes
    bool needs_manual_review = false;
    std::string alignment_status;  // ALIGNED | PENDING | NOT_DISPATCHED | CHUNK_FAILED | CHUNK_TIMEOUT
    std::string alignment_detail;
};

struct DeletedBlock {
    std::string seq_id;
    std::string bucket;
    long long block_start_orig = 0;
    long long block_end_orig = 0;
    long long block_len = 0;
    std::string left_ref_pos;
    std::string right_ref_pos;
};

struct QcStats {
    long long lost_count = 0;
    long long lost_block_count = 0;
    long long max_lost_block_len = 0;
};

struct SampleRow {
    std::string seq_id;
    std::string bucket;
    size_t seq_len = 0;
    double len_diff_ratio = 0.0;
    double n_ratio = 0.0;
    long long lost_count = 0;
    double lost_ratio = 0.0;
    long long lost_block_count = 0;
    long long max_lost_block_len = 0;
    std::string qc_level;
};

struct ChunkResult {
    std::string bucket_name;
    int chunk_id = 0;
    std::string stdout_data;
    bool success = false;
    bool timed_out = false;
    bool slow = false;
    double elapsed_sec = 0.0;
    std::string err_path;
    std::string failure_reason;
    size_t record_count = 0;
    size_t total_bases = 0;
    std::vector<std::string> seq_ids;
    std::vector<DeletedBlock> deleted_blocks;
    std::unordered_map<std::string, QcStats> qc_sample_stats;
};

