#pragma once
#include "common.hpp"

#include <string>
#include <unordered_map>
#include <vector>

void write_deleted_blocks_tsv(const std::vector<DeletedBlock>& rows, const std::string& path);
void write_sample_qc_summary_tsv(const std::vector<SampleRow>& rows, const std::string& path);
void write_abnormal_deleted_blocks_tsv(
    const std::vector<DeletedBlock>& block_rows,
    const std::unordered_map<std::string, SampleRow>& sample_row_map,
    const std::string& path
);
void make_filtered_msa(
    const std::string& full_msa_path,
    const std::string& filtered_msa_path,
    const std::unordered_map<std::string, SampleRow>& sample_row_map
);

void write_qc_summary_tsv(
    const std::vector<InputQcRow>& rows,
    const std::unordered_map<std::string, SampleRow>& post_qc,
    const std::string& path
);
void write_excluded_samples_tsv(const std::vector<InputQcRow>& rows, const std::string& path);
void write_unaligned_samples_tsv(const std::vector<InputQcRow>& rows, const std::string& path);
void write_abnormal_samples_tsv(
    const std::vector<InputQcRow>& rows,
    const std::unordered_map<std::string, SampleRow>& post_qc,
    const std::string& path
);
void write_failed_chunks_tsv(const std::vector<ChunkResult>& rows, const std::string& path);
void write_slow_chunks_tsv(const std::vector<ChunkResult>& rows, const std::string& path);
