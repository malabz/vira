#pragma once
#include "common.hpp"

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

std::vector<std::string> split_keywords(const std::string& csv);
std::string lowercase_copy(const std::string& s);
std::string sanitize_tsv_field(const std::string& s);
std::string join_reason_codes(const std::vector<std::string>& reasons);

struct InputQcDecision {
    FastaRecord record;
    InputQcRow row;
    bool dispatch = false;
};

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
