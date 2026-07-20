#include "input_qc.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <unordered_set>

namespace {
bool contains_keyword(const std::string& lower_header, const std::vector<std::string>& keywords) {
    for (const auto& kw : keywords) {
        if (!kw.empty() && lower_header.find(kw) != std::string::npos) return true;
    }
    return false;
}

bool contains_invalid_sequence_char(const std::string& seq, const std::string& allowed_chars) {
    std::unordered_set<char> allowed;
    allowed.reserve(allowed_chars.size() * 2 + 1);
    for (char c : allowed_chars) allowed.insert(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));

    for (char c : seq) {
        const char upper = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (allowed.find(upper) == allowed.end()) return true;
    }
    return false;
}

void add_reason(std::vector<std::string>& reasons, const std::string& reason) {
    if (std::find(reasons.begin(), reasons.end(), reason) == reasons.end()) reasons.push_back(reason);
}

std::string renamed_duplicate_id(const std::string& id, size_t duplicate_index) {
    std::ostringstream oss;
    oss << id << "__dup" << duplicate_index;
    return oss.str();
}
}

std::vector<std::string> split_keywords(const std::string& csv) {
    std::vector<std::string> out;
    std::string token;
    std::istringstream iss(csv);
    while (std::getline(iss, token, ',')) {
        size_t start = 0;
        while (start < token.size() && std::isspace(static_cast<unsigned char>(token[start]))) ++start;
        size_t end = token.size();
        while (end > start && std::isspace(static_cast<unsigned char>(token[end - 1]))) --end;
        if (end > start) out.push_back(lowercase_copy(token.substr(start, end - start)));
    }
    return out;
}

std::string lowercase_copy(const std::string& s) {
    std::string out = s;
    std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return out;
}

std::string sanitize_tsv_field(const std::string& s) {
    std::string out = s;
    for (char& c : out) {
        if (c == '\t' || c == '\n' || c == '\r') c = ' ';
    }
    return out;
}

std::string join_reason_codes(const std::vector<std::string>& reasons) {
    std::ostringstream oss;
    for (size_t i = 0; i < reasons.size(); ++i) {
        if (i) oss << ';';
        oss << reasons[i];
    }
    return oss.str();
}

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
    const std::vector<std::string>& engineered_keywords) {

    InputQcDecision decision;
    decision.record = std::move(rec);

    InputQcRow row;
    row.original_id = decision.record.id;
    row.seq_id = decision.record.id;
    row.header = decision.record.header;
    row.bucket = bucket;
    row.seq_len = seq_len;
    row.n_ratio = n_ratio;
    row.len_diff_ratio = len_diff_ratio;

    std::vector<std::string> reasons;
    bool hard_exclude = false;
    bool review_not_aligned = false;
    bool soft_review = false;

    auto seen_it = seen_ids.find(decision.record.id);
    if (seen_it == seen_ids.end()) {
        seen_ids[decision.record.id] = 1;
    } else {
        seen_it->second += 1;
        if (args.duplicate_id_policy == "rename") {
            const std::string new_id = renamed_duplicate_id(decision.record.id, seen_it->second);
            decision.record.id = new_id;
            decision.record.header = new_id + " " + decision.record.header;
            row.seq_id = new_id;
            row.header = decision.record.header;
            add_reason(reasons, "DUPLICATE_ID_RENAMED");
            soft_review = true;
        } else {
            add_reason(reasons, "DUPLICATE_ID");
            hard_exclude = true;
        }
    }

    if (seq_len == 0) {
        add_reason(reasons, "EMPTY_SEQUENCE");
        hard_exclude = true;
    }

    if (contains_invalid_sequence_char(decision.record.seq, args.allowed_seq_chars)) {
        add_reason(reasons, "ILLEGAL_SEQUENCE_CHARACTER");
        hard_exclude = true;
    }

    if (seq_len > 0 && n_ratio >= args.exclude_n_ratio) {
        add_reason(reasons, "N_RATIO_EXCLUDE");
        hard_exclude = true;
    } else if (seq_len > 0 && n_ratio >= args.review_n_ratio) {
        add_reason(reasons, "N_RATIO_REVIEW");
        soft_review = true;
    }

    const double len_ratio = ref_len > 0 ? static_cast<double>(seq_len) / static_cast<double>(ref_len) : 0.0;
    if (ref_len > 0 && (len_ratio < args.exclude_min_len_ratio || len_ratio > args.exclude_max_len_ratio)) {
        add_reason(reasons, "EXTREME_LENGTH_EXCLUDE");
        hard_exclude = true;
    } else if (ref_len > 0 && (len_ratio < args.reasonable_min_len_ratio || len_ratio > args.reasonable_max_len_ratio)) {
        add_reason(reasons, "LENGTH_OUTSIDE_REASONABLE_RANGE");
        soft_review = true;
    }

    const std::string lower_header = lowercase_copy(decision.record.header);
    if (contains_keyword(lower_header, non_target_keywords)) {
        add_reason(reasons, "SUSPECTED_NON_TARGET_HEADER");
        review_not_aligned = true;
    }
    if (contains_keyword(lower_header, engineered_keywords)) {
        add_reason(reasons, "SUSPECTED_ENGINEERED_CONSTRUCT");
        review_not_aligned = true;
    }

    if (hard_exclude) {
        row.decision = "EXCLUDE";
        row.alignment_status = "NOT_DISPATCHED";
        decision.dispatch = false;
    } else if (review_not_aligned && !args.align_review_samples) {
        row.decision = "REVIEW_NOT_ALIGNED";
        row.alignment_status = "NOT_DISPATCHED";
        row.needs_manual_review = true;
        decision.dispatch = false;
    } else if (review_not_aligned || soft_review) {
        row.decision = "ALIGN_REVIEW";
        row.alignment_status = "PENDING";
        row.needs_manual_review = true;
        decision.dispatch = true;
    } else {
        row.decision = "ALIGN";
        row.alignment_status = "PENDING";
        decision.dispatch = true;
    }

    if (!row.needs_manual_review && (soft_review || review_not_aligned)) row.needs_manual_review = true;
    row.reasons = join_reason_codes(reasons);
    decision.row = std::move(row);
    return decision;
}
