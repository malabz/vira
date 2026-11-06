/**
 * @file qc.cpp
 * @brief Implements parsing and classification of MAFFT `--mapout` quality data.
 *
 * Map records are converted into contiguous deleted blocks per sequence, then
 * summarized into lost-base counts and threshold-based PASS/WARN/FAIL levels.
 * Blank and comment lines are accepted, while malformed structure or numeric
 * fields are rejected with path and line-number context.
 */
#include "qc.hpp"

#include <algorithm>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace {
/**
 * @brief Splits a map-data row into exactly three trimmed CSV fields.
 * @param line Input row. Commas after the second delimiter remain in the third field.
 * @return Three trimmed fields, or an empty vector when fewer than two commas exist.
 */
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

/**
 * @brief Parses one integer field without accepting trailing characters.
 * @param value Trimmed field value.
 * @param field_name Human-readable field name for diagnostics.
 * @param map_path Source mapout path.
 * @param line_number One-based source line number.
 * @return Parsed signed integer.
 * @throws MapoutFormatError If the value is invalid, non-positive, unsafe for
 *         coordinate arithmetic, out of range, or only partially consumed.
 */
long long parse_mapout_integer(
    const std::string& value,
    const std::string& field_name,
    const std::string& map_path,
    size_t line_number) {

    size_t consumed = 0;
    long long parsed = 0;
    try {
        parsed = std::stoll(value, &consumed);
    } catch (const std::invalid_argument&) {
        throw MapoutFormatError(
            "Malformed MAFFT mapout file '" + map_path + "' at line " + std::to_string(line_number) +
            ": expected an integer for " + field_name + ", got '" + value + "'."
        );
    } catch (const std::out_of_range&) {
        throw MapoutFormatError(
            "Malformed MAFFT mapout file '" + map_path + "' at line " + std::to_string(line_number) +
            ": value out of range for " + field_name + ": '" + value + "'."
        );
    }

    if (consumed != value.size()) {
        throw MapoutFormatError(
            "Malformed MAFFT mapout file '" + map_path + "' at line " + std::to_string(line_number) +
            ": trailing characters in " + field_name + ": '" + value + "'."
        );
    }
    if (parsed <= 0 || parsed == std::numeric_limits<long long>::max()) {
        throw MapoutFormatError(
            "Malformed MAFFT mapout file '" + map_path + "' at line " + std::to_string(line_number) +
            ": " + field_name + " must be between 1 and " +
            std::to_string(std::numeric_limits<long long>::max() - 1) + ", got '" + value + "'."
        );
    }
    return parsed;
}

/**
 * @brief Converts one sequence's ordered map records into deleted blocks and statistics.
 *
 * Adjacent `-` records form one block only when their original coordinates are
 * consecutive. The nearest non-deleted mapping on either side is recorded when
 * available; edge blocks leave the absent boundary string empty.
 *
 * @param seq_id Identifier attached to every emitted DeletedBlock.
 * @param recs Ordered pairs of original coordinate and mapped reference value.
 * @return Deleted blocks together with aggregate lost-base statistics. Empty
 *         inputs or an empty identifier produce empty results.
 */
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

/**
 * @brief Parses a MAFFT map file with strict structural and numeric validation.
 *
 * Sequence sections begin with `>`. Blank and comment lines are ignored. Data
 * rows outside a section, empty section identifiers, invalid field counts, and
 * invalid coordinates terminate parsing with a contextual exception.
 */
std::pair<std::vector<DeletedBlock>, std::unordered_map<std::string, QcStats>> parse_map_to_deleted_blocks(
    const std::string& map_path) {

    std::ifstream in(map_path);
    std::vector<DeletedBlock> all_blocks;
    std::unordered_map<std::string, QcStats> sample_stats;
    if (!in) throw MapoutIoError("Cannot open MAFFT mapout file: " + map_path + ".");

    std::string current_seq_id;
    std::vector<std::pair<long long, std::string>> records;

    auto flush = [&]() {
        auto [blocks, stats] = flush_map_sequence(current_seq_id, records);
        all_blocks.insert(all_blocks.end(), blocks.begin(), blocks.end());
        if (!current_seq_id.empty()) sample_stats[current_seq_id] = stats;
        records.clear();
    };

    std::string line;
    size_t line_number = 0;
    bool saw_section = false;
    std::unordered_set<std::string> section_ids;
    while (std::getline(in, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.find('\0') != std::string::npos) {
            throw MapoutFormatError(
                "Malformed MAFFT mapout file '" + map_path + "' at line " + std::to_string(line_number) +
                ": embedded NUL byte detected; the file may be binary or corrupted."
            );
        }
        const size_t first_content = line.find_first_not_of(" \t");
        if (first_content == std::string::npos) continue;
        if (line[0] == '>') {
            const std::string next_seq_id = line.substr(1);
            if (next_seq_id.find_first_not_of(" \t") == std::string::npos) {
                throw MapoutFormatError(
                    "Malformed MAFFT mapout file '" + map_path + "' at line " + std::to_string(line_number) +
                    ": sequence section identifier is empty."
                );
            }
            if (!section_ids.insert(next_seq_id).second) {
                throw MapoutFormatError(
                    "Malformed MAFFT mapout file '" + map_path + "' at line " + std::to_string(line_number) +
                    ": duplicate sequence section identifier '" + next_seq_id + "'."
                );
            }
            flush();
            current_seq_id = next_seq_id;
            saw_section = true;
            continue;
        }
        if (line[first_content] == '#') continue;

        if (current_seq_id.empty()) {
            throw MapoutFormatError(
                "Malformed MAFFT mapout file '" + map_path + "' at line " + std::to_string(line_number) +
                ": data row appears before the first sequence section."
            );
        }

        auto parts = split_csv3(line);
        if (parts.size() != 3) {
            throw MapoutFormatError(
                "Malformed MAFFT mapout file '" + map_path + "' at line " + std::to_string(line_number) +
                ": expected exactly three comma-separated fields."
            );
        }

        if (parts[0].empty()) {
            throw MapoutFormatError(
                "Malformed MAFFT mapout file '" + map_path + "' at line " + std::to_string(line_number) +
                ": sequence-symbol field is empty."
            );
        }
        const long long orig_pos = parse_mapout_integer(parts[1], "original position", map_path, line_number);
        if (parts[2] == "-") {
            records.push_back({orig_pos, "-"});
        } else {
            (void)parse_mapout_integer(parts[2], "reference position", map_path, line_number);
            records.push_back({orig_pos, parts[2]});
        }
    }
    if (in.bad()) throw MapoutIoError("I/O error while reading MAFFT mapout file: " + map_path + ".");
    if (!saw_section) {
        throw MapoutFormatError("Malformed MAFFT mapout file '" + map_path + "': no sequence sections were found.");
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
