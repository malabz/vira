#include "args.hpp"
#include "sequence_stats.hpp"
#include "common.hpp"
#include "fasta_io.hpp"
#include "input_qc.hpp"
#include "mafft_runner.hpp"
#include "path_utils.hpp"
#include "qc.hpp"
#include "report.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
double seconds_since(const Clock::time_point& start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

constexpr const char* DEFAULT_BUCKET = "all";

std::string chunk_label(const ChunkResult& res) {
    std::ostringstream oss;
    oss << res.bucket_name << ".chunk_" << std::setw(8) << std::setfill('0') << res.chunk_id;
    return oss.str();
}

void print_configuration(
    size_t ref_len,
    size_t total_queries,
    double count_elapsed,
    const Args& args,
    int effective_chunk_size,
    size_t rough_chunk_count,
    size_t target_chunk_bases,
    const std::string& temp_dir) {

    std::cout << "=== Configuration ===\n";
    std::cout << "Reference length: " << ref_len << " bp\n";
    std::cout << "Query sequences: " << total_queries << "\n";
    std::cout << "Sequence counting time: " << std::fixed << std::setprecision(2) << count_elapsed << " seconds\n";
    std::cout << "Concurrent MAFFT processes: " << args.max_instances << "\n";
    std::cout << "Threads per MAFFT process: " << args.threads_per_mafft << "\n";
    std::cout << "Requested chunk size: " << args.chunk_size << "\n";
    std::cout << "Effective chunk size: " << effective_chunk_size << "\n";
    std::cout << "Unified dispatch bucket: " << DEFAULT_BUCKET << "\n";
    std::cout << "Target chunk bases: " << target_chunk_bases << "\n";
    std::cout << "Rough pre-QC chunk count estimate: " << rough_chunk_count << "\n";
    std::cout << "Review N-ratio threshold: " << args.review_n_ratio << "\n";
    std::cout << "Exclude N-ratio threshold: " << args.exclude_n_ratio << "\n";
    std::cout << "Reasonable length ratio range: [" << args.reasonable_min_len_ratio << ", " << args.reasonable_max_len_ratio << "]\n";
    std::cout << "Extreme exclusion length ratio range: [" << args.exclude_min_len_ratio << ", " << args.exclude_max_len_ratio << "]\n";
    std::cout << "Duplicate ID policy: " << args.duplicate_id_policy << "\n";
    std::cout << "Align review-only records: " << (args.align_review_samples ? "true" : "false") << "\n";
    std::cout << "MAFFT timeout per chunk: " << args.mafft_timeout_sec << " seconds\n";
    std::cout << "Slow chunk warning threshold: " << args.slow_chunk_sec << " seconds\n";
    std::cout << "Temporary directory: " << temp_dir << "\n";
    std::cout << "Remove temporary directory: " << (args.cleanup ? "true" : "false") << "\n";
    std::cout << "Post-alignment QC enabled: " << (args.enable_qc ? "true" : "false") << "\n";
    std::cout << "Output order is completion-based: true\n";
    std::cout << "=====================\n\n";
}
}

int main(int argc, char** argv) {
    const auto total_start = Clock::now();

    try {
        Args args = parse_args(argc, argv);
        const std::string ref_file = absolute_path(args.ref);
        const std::string query_file = absolute_path(args.query);
        const std::string final_out = absolute_path(args.out);

        if (!file_exists(ref_file)) throw std::runtime_error("Reference file not found: " + ref_file);
        if (!file_exists(query_file)) throw std::runtime_error("Query file not found: " + query_file);

        const std::string temp_dir = prepare_temp_dir(args, final_out);
        FastaRecord ref_record = read_single_fasta(ref_file);
        const size_t ref_len = ref_record.seq.size();
        if (ref_len == 0) throw std::runtime_error("Reference FASTA sequence is empty.");

        const auto count_start = Clock::now();
        const size_t total_queries = count_fasta_sequences(query_file);
        const double count_elapsed = seconds_since(count_start);

        if (total_queries == 0) {
            std::ofstream out(final_out);
            write_fasta_record(out, ref_record);
            std::cout << "The query FASTA contains no sequences; only the reference sequence was written.\n";
            std::cout << "Output written to: " << final_out << "\n";
            return 0;
        }

        int effective_chunk_size;
        if (args.disable_auto_chunk) {
            effective_chunk_size = args.chunk_size;
        } else {
            const int adaptive = std::max(1, static_cast<int>((total_queries + args.max_instances - 1) / args.max_instances));
            effective_chunk_size = std::min(args.chunk_size, adaptive);
        }

        const size_t target_chunk_bases = args.target_chunk_bases > 0
            ? args.target_chunk_bases
            : static_cast<size_t>(effective_chunk_size) * ref_len;
        const size_t rough_chunk_count = (total_queries + static_cast<size_t>(effective_chunk_size) - 1) / static_cast<size_t>(effective_chunk_size);

        print_configuration(
            ref_len,
            total_queries,
            count_elapsed,
            args,
            effective_chunk_size,
            rough_chunk_count,
            target_chunk_bases,
            temp_dir
        );

        const auto qc_and_dispatch_start = Clock::now();

        std::map<std::string, std::vector<FastaRecord>> buffers;
        std::map<std::string, size_t> buffer_bases;
        std::map<std::string, int> chunk_counters;
        std::map<std::string, size_t> bucket_counts;
        buffers[DEFAULT_BUCKET] = {};
        buffer_bases[DEFAULT_BUCKET] = 0;
        chunk_counters[DEFAULT_BUCKET] = 0;
        bucket_counts[DEFAULT_BUCKET] = 0;

        std::vector<InputQcRow> input_rows;
        input_rows.reserve(total_queries);
        std::unordered_map<std::string, size_t> row_index_by_seq_id;
        row_index_by_seq_id.reserve(total_queries * 2 + 1);
        std::unordered_map<std::string, size_t> seen_ids;
        seen_ids.reserve(total_queries * 2 + 1);

        std::map<std::string, SampleMeta> sample_meta;
        std::vector<DeletedBlock> all_deleted_blocks;
        std::unordered_map<std::string, QcStats> sample_qc_stats;
        std::vector<ChunkResult> failed_chunks;
        std::vector<ChunkResult> slow_chunks;

        const size_t max_pending_futures = static_cast<size_t>(args.max_instances * args.max_pending_factor);
        std::vector<std::future<ChunkResult>> futures;

        const auto non_target_keywords = split_keywords(args.non_target_keywords);
        const auto engineered_keywords = split_keywords(args.engineered_keywords);

        size_t dispatchable_count = 0;
        size_t excluded_count = 0;
        size_t review_not_aligned_count = 0;
        size_t align_review_count = 0;

        auto handle_result = [&](ChunkResult&& res, std::ostream& out_f) {
            const std::string label = chunk_label(res);
            if (res.success) {
                append_mafft_stdout_without_ref(res.stdout_data, out_f);
                for (const auto& seq_id : res.seq_ids) {
                    auto it = row_index_by_seq_id.find(seq_id);
                    if (it != row_index_by_seq_id.end()) {
                        input_rows[it->second].alignment_status = "ALIGNED";
                        input_rows[it->second].alignment_detail = label;
                    }
                }
                if (args.enable_qc) {
                    for (auto& row : res.deleted_blocks) {
                        row.bucket = res.bucket_name;
                        all_deleted_blocks.push_back(std::move(row));
                    }
                    for (auto& kv : res.qc_sample_stats) sample_qc_stats[kv.first] = kv.second;
                }
            } else {
                for (const auto& seq_id : res.seq_ids) {
                    auto it = row_index_by_seq_id.find(seq_id);
                    if (it != row_index_by_seq_id.end()) {
                        input_rows[it->second].alignment_status = res.timed_out ? "CHUNK_TIMEOUT" : "CHUNK_FAILED";
                        input_rows[it->second].alignment_detail = label + ":" + res.failure_reason;
                    }
                }
                failed_chunks.push_back(res);
            }
            if (res.slow) slow_chunks.push_back(res);
        };

        auto drain_one_finished_or_first = [&](std::ostream& out_f) {
            while (true) {
                for (auto it = futures.begin(); it != futures.end(); ++it) {
                    if (it->wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
                        handle_result(it->get(), out_f);
                        futures.erase(it);
                        return;
                    }
                }
                if (!futures.empty()) {
                    handle_result(futures.front().get(), out_f);
                    futures.erase(futures.begin());
                    return;
                }
            }
        };

        std::ofstream out_f(final_out);
        if (!out_f) throw std::runtime_error("Cannot write final output: " + final_out);
        write_fasta_record(out_f, ref_record);

        auto dispatch_buffer = [&](const std::string& bucket) {
            auto& records = buffers[bucket];
            if (records.empty()) return;

            const int chunk_id = chunk_counters[bucket]++;
            std::ostringstream name;
            name << bucket << ".chunk_" << std::setw(8) << std::setfill('0') << chunk_id << ".query.fasta";
            const std::string chunk_path = join_path(temp_dir, name.str());

            std::vector<std::string> seq_ids;
            seq_ids.reserve(records.size());
            for (const auto& r : records) seq_ids.push_back(r.id);
            const size_t total_bases = buffer_bases[bucket];

            write_chunk_file(records, chunk_path);
            records.clear();
            buffer_bases[bucket] = 0;

            futures.push_back(std::async(
                std::launch::async,
                run_mafft_chunk,
                chunk_id,
                bucket,
                chunk_path,
                ref_file,
                args.threads_per_mafft,
                temp_dir,
                args.save_failed_stderr,
                args.enable_qc,
                args.mafft_timeout_sec,
                args.slow_chunk_sec,
                std::move(seq_ids),
                total_bases
            ));

            if (futures.size() >= max_pending_futures) drain_one_finished_or_first(out_f);
        };

        std::cout << ">> Step 1/2: Streaming query records, applying input QC, and dispatching unified load-aware MAFFT jobs...\n";
        std::ifstream qin(query_file);
        if (!qin) throw std::runtime_error("Cannot open query FASTA: " + query_file);

        FastaRecord rec;
        std::string pending;
        size_t read_count = 0;

        while (read_next_fasta_record(qin, rec, pending)) {
            const SequenceStats stats = calculate_sequence_stats(rec, ref_len);
            const std::string bucket = DEFAULT_BUCKET;
            const size_t seq_len = stats.seq_len;
            const double n_ratio = stats.n_ratio;
            const double len_diff_ratio = stats.len_diff_ratio;

            InputQcDecision input_decision = assess_input_record(
                std::move(rec),
                args,
                ref_len,
                bucket,
                seq_len,
                n_ratio,
                len_diff_ratio,
                seen_ids,
                non_target_keywords,
                engineered_keywords
            );

            const size_t row_index = input_rows.size();
            input_rows.push_back(input_decision.row);

            if (input_decision.dispatch) {
                dispatchable_count++;
                if (input_rows[row_index].decision == "ALIGN_REVIEW") align_review_count++;
                row_index_by_seq_id[input_decision.record.id] = row_index;
                sample_meta[input_decision.record.id] = SampleMeta{bucket, seq_len, n_ratio, len_diff_ratio};

                buffers[bucket].push_back(std::move(input_decision.record));
                buffer_bases[bucket] += seq_len;
                bucket_counts[bucket]++;

                if (static_cast<int>(buffers[bucket].size()) >= effective_chunk_size || buffer_bases[bucket] >= target_chunk_bases) {
                    dispatch_buffer(bucket);
                }
            } else {
                if (input_rows[row_index].decision == "EXCLUDE") excluded_count++;
                if (input_rows[row_index].decision == "REVIEW_NOT_ALIGNED") review_not_aligned_count++;
            }

            read_count++;
            if (read_count % 100000 == 0) std::cout << "Processed " << read_count << " query sequences...\n";
        }

        dispatch_buffer(DEFAULT_BUCKET);

        const double qc_and_dispatch_elapsed = seconds_since(qc_and_dispatch_start);

        std::cout << "\n>> Step 2/2: Waiting for remaining MAFFT jobs and writing completed chunks...\n";
        const auto wait_start = Clock::now();
        const size_t remaining_total = futures.size();
        size_t completed = 0;

        while (!futures.empty()) {
            drain_one_finished_or_first(out_f);
            completed++;
            if (completed % 10 == 0 || completed == remaining_total) {
                std::cout << "Completed " << completed << "/" << remaining_total << " remaining chunks...\n";
            }
        }
        out_f.close();
        const double wait_elapsed = seconds_since(wait_start);

        const std::string prefix = output_prefix(final_out);
        const std::string filtered_msa_path = prefix + ".filtered.fasta";
        const std::string deleted_blocks_path = prefix + ".deleted_blocks.tsv";
        const std::string sample_qc_summary_path = prefix + ".sample_qc_summary.tsv";
        const std::string abnormal_deleted_blocks_path = prefix + ".abnormal_deleted_blocks.tsv";
        const std::string qc_summary_path = prefix + ".qc_summary.tsv";
        const std::string excluded_samples_path = prefix + ".excluded_samples.tsv";
        const std::string unaligned_samples_path = prefix + ".unaligned_samples.tsv";
        const std::string abnormal_samples_path = prefix + ".abnormal_samples.tsv";
        const std::string failed_chunks_path = prefix + ".failed_chunks.tsv";
        const std::string slow_chunks_path = prefix + ".slow_chunks.tsv";

        std::unordered_map<std::string, SampleRow> sample_row_map;
        if (args.enable_qc) {
            write_deleted_blocks_tsv(all_deleted_blocks, deleted_blocks_path);
            auto sample_rows = build_sample_qc_rows(sample_meta, sample_qc_stats, args);
            write_sample_qc_summary_tsv(sample_rows, sample_qc_summary_path);
            sample_row_map.reserve(sample_rows.size());
            for (const auto& row : sample_rows) sample_row_map[row.seq_id] = row;
            write_abnormal_deleted_blocks_tsv(all_deleted_blocks, sample_row_map, abnormal_deleted_blocks_path);
            make_filtered_msa(final_out, filtered_msa_path, sample_row_map);
        }

        write_qc_summary_tsv(input_rows, sample_row_map, qc_summary_path);
        write_excluded_samples_tsv(input_rows, excluded_samples_path);
        write_unaligned_samples_tsv(input_rows, unaligned_samples_path);
        write_abnormal_samples_tsv(input_rows, sample_row_map, abnormal_samples_path);
        write_failed_chunks_tsv(failed_chunks, failed_chunks_path);
        write_slow_chunks_tsv(slow_chunks, slow_chunks_path);

        const double total_elapsed = seconds_since(total_start);

        std::cout << "\nSequences dispatched:\n";
        std::cout << "  " << DEFAULT_BUCKET << ": " << bucket_counts[DEFAULT_BUCKET] << "\n";

        std::cout << "\nInput QC summary:\n";
        std::cout << "  Dispatchable for alignment: " << dispatchable_count << "\n";
        std::cout << "  ALIGN_REVIEW records: " << align_review_count << "\n";
        std::cout << "  Excluded records: " << excluded_count << "\n";
        std::cout << "  Review-not-aligned records: " << review_not_aligned_count << "\n";
        std::cout << "  Failed/timeout chunks: " << failed_chunks.size() << "\n";
        std::cout << "  Slow chunks: " << slow_chunks.size() << "\n";

        std::cout << "\nProcessing complete. Full MSA written to: " << final_out << "\n";
        std::cout << "Sequence counting time: " << std::fixed << std::setprecision(2) << count_elapsed << " seconds\n";
        std::cout << "Input QC/dispatch time: " << qc_and_dispatch_elapsed << " seconds\n";
        std::cout << "Final wait-and-write time: " << wait_elapsed << " seconds\n";
        std::cout << "Total runtime: " << total_elapsed << " seconds\n";
        std::cout << "QC summary: " << qc_summary_path << "\n";
        std::cout << "Excluded sample list: " << excluded_samples_path << "\n";
        std::cout << "Unaligned sample list: " << unaligned_samples_path << "\n";
        std::cout << "Abnormal sample list: " << abnormal_samples_path << "\n";
        std::cout << "Failed chunk list: " << failed_chunks_path << "\n";
        std::cout << "Slow chunk list: " << slow_chunks_path << "\n";

        if (args.enable_qc) {
            std::cout << "MAFFT mapout deleted blocks: " << deleted_blocks_path << "\n";
            std::cout << "Post-alignment QC summary: " << sample_qc_summary_path << "\n";
            std::cout << "Abnormal deletion-block summary: " << abnormal_deleted_blocks_path << "\n";
            std::cout << "Recommended filtered MSA: " << filtered_msa_path << "\n";
        }

        if (!failed_chunks.empty()) {
            std::cerr << "\nWarning: " << failed_chunks.size() << " MAFFT chunks failed or timed out. First 20 failures:\n";
            const size_t limit = std::min<size_t>(20, failed_chunks.size());
            for (size_t i = 0; i < limit; ++i) {
                const auto& res = failed_chunks[i];
                std::cerr << "  bucket=" << res.bucket_name
                          << ", chunk=" << res.chunk_id
                          << ", reason=" << res.failure_reason
                          << ", timed_out=" << (res.timed_out ? "true" : "false")
                          << ", stderr=" << res.err_path << "\n";
            }
            if (args.save_failed_stderr) std::cerr << "stderr for failed chunks has been saved in the temporary directory.\n";
        }

        if (!slow_chunks.empty()) {
            std::cerr << "\nWarning: " << slow_chunks.size() << " chunks exceeded the slow-task threshold. See: " << slow_chunks_path << "\n";
        }

        if (args.cleanup) {
            const int rc = std::system((std::string("rm -rf ") + shell_quote(temp_dir)).c_str());
            if (rc != 0) {
                std::cerr << "Warning: failed to remove temporary directory: " << temp_dir << "\n";
            } else {
                std::cout << "Temporary directory removed: " << temp_dir << "\n";
            }
        }

        return failed_chunks.empty() ? 0 : 2;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
