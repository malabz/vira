#pragma once
#include "common.hpp"

#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

std::pair<std::vector<DeletedBlock>, std::unordered_map<std::string, QcStats>> parse_map_to_deleted_blocks(
    const std::string& map_path
);

std::string assign_qc_level(double lost_ratio, long long max_lost_block_len, const Args& args);

std::vector<SampleRow> build_sample_qc_rows(
    const std::map<std::string, SampleMeta>& sample_meta,
    const std::unordered_map<std::string, QcStats>& sample_qc_stats,
    const Args& args
);
