#include "qc.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {
std::vector<std::string> split_csv3(const std::string& line) {
    std::vector<std::string> parts;
    size_t start = 0;
    for (int i = 0; i < 2; ++i) {
        size_t p = line.find(',', start);
        if (p == std::string::npos) return {};
        parts.push_back(line.substr(start, p - start));
        start = p + 1;
    }
    parts.push_back(line.substr(start));

    for (auto& s : parts) {
        size_t b = s.find_first_not_of(" \t\r\n");
        size_t e = s.find_last_not_of(" \t\r\n");
        s = (b == std::string::npos) ? "" : s.substr(b, e - b + 1);
    }
    return parts;
}

std::pair<std::vector<DeletedBlock>, QcStats> flush_map_sequence(
    const std::string& seq_id,
    const std::vector<std::pair<long long, std::string>>& recs) {

    std::vector<DeletedBlock> blocks;
    QcStats stats;
    if (seq_id.empty() || recs.empty()) return {blocks, stats};

    const long long n = static_cast<long long>(recs.size());
    for (const auto& p : recs) {
        if (p.second == "-") ++stats.lost_count;
    }

    long long i = 0;
    while (i < n) {
        if (recs[i].second != "-") {
            ++i;
            continue;
        }

        const long long start_i = i;
        const long long start_orig = recs[i].first;

        while (i + 1 < n && recs[i + 1].second == "-" && recs[i + 1].first == recs[i].first + 1) {
            ++i;
        }

        const long long end_i = i;
        const long long end_orig = recs[i].first;

        DeletedBlock b;
        b.seq_id = seq_id;
        b.block_start_orig = start_orig;
        b.block_end_orig = end_orig;
        b.block_len = end_orig - start_orig + 1;

        for (long long j = start_i - 1; j >= 0; --j) {
            if (recs[j].second != "-") {
                b.left_ref_pos = recs[j].second;
                break;
            }
        }
        for (long long j = end_i + 1; j < n; ++j) {
            if (recs[j].second != "-") {
                b.right_ref_pos = recs[j].second;
                break;
            }
        }

        stats.max_lost_block_len = std::max(stats.max_lost_block_len, b.block_len);
        blocks.push_back(std::move(b));
        ++i;
    }

    stats.lost_block_count = static_cast<long long>(blocks.size());
    return {blocks, stats};
}
}

std::pair<std::vector<DeletedBlock>, std::unordered_map<std::string, QcStats>> parse_map_to_deleted_blocks(
    const std::string& map_path) {

    std::ifstream in(map_path);
    std::vector<DeletedBlock> all_blocks;
    std::unordered_map<std::string, QcStats> sample_stats;
    if (!in) return {all_blocks, sample_stats};

    std::string current_seq_id;
    std::vector<std::pair<long long, std::string>> records;

    auto flush = [&]() {
        auto [blocks, stats] = flush_map_sequence(current_seq_id, records);
        all_blocks.insert(all_blocks.end(), blocks.begin(), blocks.end());
        if (!current_seq_id.empty()) sample_stats[current_seq_id] = stats;
        records.clear();
    };

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '>') {
            flush();
            current_seq_id = line.substr(1);
            continue;
        }
        if (line[0] == '#') continue;

        auto parts = split_csv3(line);
        if (parts.size() != 3) continue;

        try {
            long long orig_pos = std::stoll(parts[1]);
            if (parts[2] == "-") {
                records.push_back({orig_pos, "-"});
            } else {
                (void)std::stoll(parts[2]);
                records.push_back({orig_pos, parts[2]});
            }
        } catch (...) {
            continue;
        }
    }
    flush();

    return {all_blocks, sample_stats};
}

std::string assign_qc_level(double lost_ratio, long long max_lost_block_len, const Args& args) {
    if (lost_ratio >= args.fail_lost_ratio || max_lost_block_len >= args.fail_max_block_len) return "FAIL";
    if (lost_ratio >= args.warn_lost_ratio || max_lost_block_len >= args.warn_max_block_len) return "WARN";
    return "PASS";
}

std::vector<SampleRow> build_sample_qc_rows(
    const std::map<std::string, SampleMeta>& sample_meta,
    const std::unordered_map<std::string, QcStats>& sample_qc_stats,
    const Args& args) {

    std::vector<SampleRow> rows;
    rows.reserve(sample_meta.size());

    for (const auto& kv : sample_meta) {
        const auto& seq_id = kv.first;
        const auto& meta = kv.second;

        QcStats qc;
        auto it = sample_qc_stats.find(seq_id);
        if (it != sample_qc_stats.end()) qc = it->second;

        SampleRow row;
        row.seq_id = seq_id;
        row.bucket = meta.bucket;
        row.seq_len = meta.seq_len;
        row.len_diff_ratio = meta.len_diff_ratio;
        row.n_ratio = meta.n_ratio;
        row.lost_count = qc.lost_count;
        row.lost_ratio = meta.seq_len > 0 ? static_cast<double>(qc.lost_count) / static_cast<double>(meta.seq_len) : 0.0;
        row.lost_block_count = qc.lost_block_count;
        row.max_lost_block_len = qc.max_lost_block_len;
        row.qc_level = assign_qc_level(row.lost_ratio, row.max_lost_block_len, args);
        rows.push_back(std::move(row));
    }

    return rows;
}
