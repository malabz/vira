#pragma once
#include "common.hpp"

#include <string>
#include <vector>

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
