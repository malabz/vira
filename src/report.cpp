/**
 * @file report.cpp
 * @brief Implements tabular QC reports and filtered multiple-alignment output.
 *
 * Report writers preserve the decision and status vocabulary defined by the
 * core data model. Input-QC free-text fields and failure diagnostic paths are
 * sanitized before TSV output, while filtered FASTA output always preserves
 * the first (reference) record.
 */
#include "report.hpp"

#include "fasta_io.hpp"
#include "input_qc.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace {
/**
 * @brief Writes the shared input-QC TSV header.
 * @param f Open output stream receiving the header.
 * @param include_post_qc Whether post-alignment metric columns are appended.
 */
void write_input_qc_header(std::ofstream& f, bool include_post_qc) {
    f << "original_id\tseq_id\theader\tbucket\tseq_len\tlen_diff_ratio\tn_ratio\tinput_decision\tinput_reasons\tneeds_manual_review\talignment_status\talignment_detail";
    if (include_post_qc) {
        f << "\tlost_count\tlost_ratio\tlost_block_count\tmax_lost_block_len\tpost_alignment_qc";
    }
    f << '\n';
}

/**
 * @brief Serializes one input-QC row with optional post-alignment metrics.
 * @param f Open output stream receiving the row.
 * @param r Input decision and alignment status to serialize.
 * @param post_qc Post-alignment metrics keyed by final sequence identifier.
 * @param include_post_qc Whether post-alignment columns are emitted. When no
 *        matching metrics exist, numeric fields remain empty and the QC level
 *        is written as `NA`.
 */
void write_input_qc_row(
    std::ofstream& f,
    const InputQcRow& r,
    const std::unordered_map<std::string, SampleRow>& post_qc,
    bool include_post_qc) {
    f << sanitize_tsv_field(r.original_id) << '\t'
      << sanitize_tsv_field(r.seq_id) << '\t'
      << sanitize_tsv_field(r.header) << '\t'
      << r.bucket << '\t'
      << r.seq_len << '\t'
      << r.len_diff_ratio << '\t'
      << r.n_ratio << '\t'
      << r.decision << '\t'
      << r.reasons << '\t'
      << (r.needs_manual_review ? "true" : "false") << '\t'
      << r.alignment_status << '\t'
      << sanitize_tsv_field(r.alignment_detail);

    if (include_post_qc) {
        auto it = post_qc.find(r.seq_id);
        if (it == post_qc.end()) {
            f << "\t\t\t\t\tNA";
        } else {
            const auto& q = it->second;
            f << '\t' << q.lost_count
              << '\t' << q.lost_ratio
              << '\t' << q.lost_block_count
              << '\t' << q.max_lost_block_len
              << '\t' << q.qc_level;
        }
    }
    f << '\n';
}

/**
 * @brief Determines whether a sample belongs in the abnormal-sample report.
 * @param r Input-QC and alignment state for the sample.
 * @param post_qc Optional post-alignment QC rows keyed by sequence identifier.
 * @return `true` for any non-ALIGN decision, non-ALIGNED status, or available
 *         post-alignment WARN/FAIL result; otherwise `false`.
 */
bool is_abnormal(
    const InputQcRow& r,
    const std::unordered_map<std::string, SampleRow>& post_qc) {
    if (r.decision != "ALIGN") return true;
    if (r.alignment_status != "ALIGNED") return true;
    auto it = post_qc.find(r.seq_id);
    return it != post_qc.end() && (it->second.qc_level == "WARN" || it->second.qc_level == "FAIL");
}
}

void write_deleted_blocks_tsv(const std::vector<DeletedBlock>& rows, const std::string& path) {
    std::ofstream f(path);
    f << "seq_id\tbucket\tblock_start_orig\tblock_end_orig\tblock_len\tleft_ref_pos\tright_ref_pos\n";
    for (const auto& r : rows) {
        f << r.seq_id << '\t' << r.bucket << '\t' << r.block_start_orig << '\t' << r.block_end_orig << '\t'
          << r.block_len << '\t' << r.left_ref_pos << '\t' << r.right_ref_pos << '\n';
    }
}

void write_sample_qc_summary_tsv(const std::vector<SampleRow>& rows, const std::string& path) {
    std::ofstream f(path);
    f << "seq_id\tbucket\tseq_len\tlen_diff_ratio\tn_ratio\tlost_count\tlost_ratio\tlost_block_count\tmax_lost_block_len\tqc_level\n";
    f << std::fixed << std::setprecision(8);
    for (const auto& r : rows) {
        f << r.seq_id << '\t' << r.bucket << '\t' << r.seq_len << '\t' << r.len_diff_ratio << '\t' << r.n_ratio << '\t'
          << r.lost_count << '\t' << r.lost_ratio << '\t' << r.lost_block_count << '\t' << r.max_lost_block_len << '\t' << r.qc_level << '\n';
    }
}

void write_abnormal_deleted_blocks_tsv(
    const std::vector<DeletedBlock>& block_rows,
    const std::unordered_map<std::string, SampleRow>& sample_row_map,
    const std::string& path) {

    std::ofstream f(path);
    f << "seq_id\tbucket\tqc_level\tblock_start_orig\tblock_end_orig\tblock_len\tleft_ref_pos\tright_ref_pos\n";
    for (const auto& r : block_rows) {
        auto it = sample_row_map.find(r.seq_id);
        if (it == sample_row_map.end()) continue;
        if (it->second.qc_level == "WARN" || it->second.qc_level == "FAIL") {
            f << r.seq_id << '\t' << r.bucket << '\t' << it->second.qc_level << '\t' << r.block_start_orig << '\t'
              << r.block_end_orig << '\t' << r.block_len << '\t' << r.left_ref_pos << '\t' << r.right_ref_pos << '\n';
        }
    }
}

/**
 * @brief Copies the reference and all known non-FAIL query records to a new MSA.
 *
 * The first FASTA record is treated as the reference and is always retained.
 * Subsequent records are retained only when their identifiers exist in
 * @p sample_row_map and their level is not `FAIL`; records without QC metadata
 * are therefore omitted.
 *
 * @throws std::runtime_error If the source MSA cannot be opened or the filtered
 *         output cannot be created.
 */
void make_filtered_msa(
    const std::string& full_msa_path,
    const std::string& filtered_msa_path,
    const std::unordered_map<std::string, SampleRow>& sample_row_map) {

    std::ifstream in(full_msa_path);
    if (!in) throw std::runtime_error("Cannot open full MSA for filtering: " + full_msa_path);

    std::ofstream out(filtered_msa_path);
    if (!out) throw std::runtime_error("Cannot write filtered MSA: " + filtered_msa_path);

    FastaRecord rec;
    std::string pending;
    bool first_record = true;

    while (read_next_fasta_record(in, rec, pending)) {
        if (first_record) {
            write_fasta_record(out, rec);
            first_record = false;
            continue;
        }

        auto it = sample_row_map.find(rec.id);
        if (it != sample_row_map.end() && it->second.qc_level != "FAIL") {
            write_fasta_record(out, rec);
        }
    }
}

void write_qc_summary_tsv(
    const std::vector<InputQcRow>& rows,
    const std::unordered_map<std::string, SampleRow>& post_qc,
    const std::string& path) {
    std::ofstream f(path);
    f << std::fixed << std::setprecision(8);
    write_input_qc_header(f, true);
    for (const auto& r : rows) write_input_qc_row(f, r, post_qc, true);
}

void write_excluded_samples_tsv(const std::vector<InputQcRow>& rows, const std::string& path) {
    std::ofstream f(path);
    f << std::fixed << std::setprecision(8);
    std::unordered_map<std::string, SampleRow> empty;
    write_input_qc_header(f, false);
    for (const auto& r : rows) {
        if (r.decision == "EXCLUDE" || r.decision == "REVIEW_NOT_ALIGNED") write_input_qc_row(f, r, empty, false);
    }
}

void write_unaligned_samples_tsv(const std::vector<InputQcRow>& rows, const std::string& path) {
    std::ofstream f(path);
    f << std::fixed << std::setprecision(8);
    std::unordered_map<std::string, SampleRow> empty;
    write_input_qc_header(f, false);
    for (const auto& r : rows) {
        if (r.alignment_status != "ALIGNED") write_input_qc_row(f, r, empty, false);
    }
}

void write_abnormal_samples_tsv(
    const std::vector<InputQcRow>& rows,
    const std::unordered_map<std::string, SampleRow>& post_qc,
    const std::string& path) {
    std::ofstream f(path);
    f << std::fixed << std::setprecision(8);
    write_input_qc_header(f, true);
    for (const auto& r : rows) {
        if (is_abnormal(r, post_qc)) write_input_qc_row(f, r, post_qc, true);
    }
}

void write_failed_chunks_tsv(const std::vector<ChunkResult>& rows, const std::string& path) {
    std::ofstream f(path);
    f << "bucket\tchunk_id\trecord_count\ttotal_bases\telapsed_sec\ttimed_out\tfailure_reason\tstderr_path\n";
    f << std::fixed << std::setprecision(3);
    for (const auto& r : rows) {
        f << r.bucket_name << '\t' << r.chunk_id << '\t' << r.record_count << '\t' << r.total_bases << '\t'
          << r.elapsed_sec << '\t' << (r.timed_out ? "true" : "false") << '\t' << r.failure_reason << '\t'
          << sanitize_tsv_field(r.err_path) << '\n';
    }
}

void write_slow_chunks_tsv(const std::vector<ChunkResult>& rows, const std::string& path) {
    std::ofstream f(path);
    f << "bucket\tchunk_id\trecord_count\ttotal_bases\telapsed_sec\tsuccess\ttimed_out\n";
    f << std::fixed << std::setprecision(3);
    for (const auto& r : rows) {
        f << r.bucket_name << '\t' << r.chunk_id << '\t' << r.record_count << '\t' << r.total_bases << '\t'
          << r.elapsed_sec << '\t' << (r.success ? "true" : "false") << '\t' << (r.timed_out ? "true" : "false") << '\n';
    }
}
