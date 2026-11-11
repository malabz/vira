/**
 * @file main.cpp
 * @brief Orchestrates ViRA input validation, chunked MAFFT execution, QC aggregation,
 *        and report generation.
 *
 * The entry point streams query records through input QC, bounds the number of
 * outstanding asynchronous MAFFT jobs, and writes successful chunks in completion
 * order. Algorithmic work is delegated to the module interfaces declared in
 * `include/`.
 */

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
#include <cerrno>
#include <cctype>
#include <chrono>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <sys/stat.h>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {
/**
 * @brief Measures elapsed monotonic time from a previously captured time point.
 * @param start Start time obtained from the shared steady `Clock`.
 * @return Elapsed seconds; non-negative when `start` is not later than the current time.
 */
double seconds_since(const Clock::time_point& start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

/// Logical bucket used by the current unified dispatch strategy.
constexpr const char* DEFAULT_BUCKET = "all";

/**
 * @brief Builds the stable human-readable identifier used in status reports.
 * @param res Chunk result containing the bucket name and numeric chunk identifier.
 * @return A label of the form `BUCKET.chunk_NNNNNNNN`.
 */
std::string chunk_label(const ChunkResult& res) {
    std::ostringstream oss;
    oss << res.bucket_name << ".chunk_" << std::setw(8) << std::setfill('0') << res.chunk_id;
    return oss.str();
}

/** @brief Returns whether two existing paths resolve to the same filesystem object. */
bool paths_refer_to_same_file(const std::string& a, const std::string& b) {
    struct stat a_stat {};
    struct stat b_stat {};
    return stat(a.c_str(), &a_stat) == 0 &&
           stat(b.c_str(), &b_stat) == 0 &&
           a_stat.st_dev == b_stat.st_dev &&
           a_stat.st_ino == b_stat.st_ino;
}

struct OutputTarget {
    std::string name;
    std::string path;
};

/**
 * @brief Rejects symbolic-link targets and aliases among inputs and active outputs.
 * @throws std::invalid_argument If an output could overwrite an input or another output.
 * @throws std::runtime_error If an existing output target cannot be inspected.
 */
void validate_distinct_output_targets(
    const std::vector<OutputTarget>& targets,
    const std::string& ref_file,
    const std::string& query_file) {

    for (size_t i = 0; i < targets.size(); ++i) {
        const auto& target = targets[i];
        struct stat link_info {};
        if (lstat(target.path.c_str(), &link_info) == 0) {
            if (S_ISLNK(link_info.st_mode)) {
                throw std::invalid_argument(
                    target.name + " path must not be a symbolic link: " + target.path + "."
                );
            }
            if (!S_ISREG(link_info.st_mode)) {
                throw std::invalid_argument(
                    target.name + " path must be a regular-file target: " + target.path + "."
                );
            }
            if (access(target.path.c_str(), W_OK) != 0) {
                throw std::runtime_error(
                    target.name + " file is not writable: " + target.path + " (" +
                    std::strerror(errno) + ")."
                );
            }
        } else if (errno != ENOENT) {
            throw std::runtime_error(
                "Cannot inspect " + target.name + " path '" + target.path + "': " +
                std::strerror(errno) + "."
            );
        }

        if (target.path == ref_file || target.path == query_file ||
            paths_refer_to_same_file(target.path, ref_file) ||
            paths_refer_to_same_file(target.path, query_file)) {
            throw std::invalid_argument(
                target.name + " path would overwrite an input file: " + target.path + "."
            );
        }

        for (size_t j = 0; j < i; ++j) {
            const auto& other = targets[j];
            if (target.path == other.path || paths_refer_to_same_file(target.path, other.path)) {
                throw std::invalid_argument(
                    target.name + " path aliases " + other.name + " path; output paths must be distinct: " +
                    target.path + "."
                );
            }
        }
    }
}

/**
 * @brief Rejects destructive output aliases and validates the output directory.
 * @throws std::invalid_argument If the output aliases an input or is not a regular-file target.
 * @throws std::runtime_error If the output location cannot be inspected or written.
 */
void validate_final_output_path(
    const std::string& final_out,
    const std::string& ref_file,
    const std::string& query_file) {

    if (final_out == ref_file || final_out == query_file ||
        paths_refer_to_same_file(final_out, ref_file) ||
        paths_refer_to_same_file(final_out, query_file)) {
        throw std::invalid_argument(
            "Output path must differ from reference and query input paths: " + final_out +
            ". Refusing to overwrite an input file."
        );
    }

    struct stat link_info {};
    if (lstat(final_out.c_str(), &link_info) == 0) {
        if (S_ISLNK(link_info.st_mode)) {
            throw std::invalid_argument("Final output path must not be a symbolic link: " + final_out + ".");
        }
    } else if (errno != ENOENT) {
        throw std::runtime_error(
            "Cannot inspect final output path: " + final_out + " (" + std::strerror(errno) + ")."
        );
    }

    struct stat output_stat {};
    if (stat(final_out.c_str(), &output_stat) == 0) {
        if (!S_ISREG(output_stat.st_mode)) {
            throw std::invalid_argument("Expected final output to be a regular file path: " + final_out + ".");
        }
        if (access(final_out.c_str(), W_OK) != 0) {
            throw std::runtime_error(
                "Final output file is not writable: " + final_out + " (" + std::strerror(errno) + ")."
            );
        }
    } else if (errno != ENOENT) {
        throw std::runtime_error(
            "Cannot inspect final output path: " + final_out + " (" + std::strerror(errno) + ")."
        );
    }

    const std::string parent = path_dirname(final_out);
    struct stat parent_stat {};
    if (stat(parent.c_str(), &parent_stat) != 0) {
        throw std::runtime_error(
            "Output directory does not exist or cannot be inspected: " + parent + " (" + std::strerror(errno) + ")."
        );
    }
    if (!S_ISDIR(parent_stat.st_mode)) {
        throw std::invalid_argument("Output parent path is not a directory: " + parent + ".");
    }
    if (access(parent.c_str(), W_OK) != 0) {
        throw std::runtime_error(
            "Output directory is not writable: " + parent + " (" + std::strerror(errno) + ")."
        );
    }
}

/** @brief Closes the final MSA stream and reports delayed write failures. */
void close_final_output_or_throw(std::ofstream& out, const std::string& path) {
    out.close();
    if (!out) {
        throw std::runtime_error(
            "Failed while writing final output file: " + path + ". The output may be incomplete."
        );
    }
}

/** @brief Validates reference symbols against the configured case-insensitive alphabet. */
void validate_reference_symbols(const FastaRecord& reference, const Args& args, const std::string& path) {
    for (size_t i = 0; i < reference.seq.size(); ++i) {
        const unsigned char raw = static_cast<unsigned char>(reference.seq[i]);
        const char normalized = static_cast<char>(std::toupper(raw));
        const bool allowed = std::any_of(
            args.allowed_seq_chars.begin(),
            args.allowed_seq_chars.end(),
            [&](unsigned char candidate) {
                return static_cast<char>(std::toupper(candidate)) == normalized;
            }
        );
        if (allowed) continue;

        std::ostringstream symbol;
        if (std::isprint(raw)) {
            symbol << "'" << reference.seq[i] << "'";
        } else {
            symbol << "byte 0x" << std::hex << std::uppercase << static_cast<int>(raw);
        }
        throw std::invalid_argument(
            "Reference FASTA contains unsupported sequence symbol " + symbol.str() +
            " at position " + std::to_string(i + 1) + " in " + path +
            ". Allowed symbols: " + args.allowed_seq_chars + "."
        );
    }
}

/**
 * @brief Prints the effective runtime configuration before records are dispatched.
 * @param ref_len Reference-sequence length in bases.
 * @param total_queries Number of query records counted before streaming begins.
 * @param count_elapsed Time spent counting query records, in seconds.
 * @param args Validated command-line configuration.
 * @param effective_chunk_size Maximum records permitted in a dispatched chunk.
 * @param rough_chunk_count Pre-QC estimate based on record count alone.
 * @param target_chunk_bases Maximum aggregate sequence bases permitted per chunk.
 * @param temp_dir Concrete temporary directory allocated for this run.
 */
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

/**
 * @brief Runs the complete ViRA command-line workflow.
 * @param argc Number of command-line arguments.
 * @param argv Argument vector; required paths and options are validated by `parse_args`.
 * @return `0` on success, `1` for configuration or I/O errors, or `2` when one or
 *         more MAFFT chunks fail or time out.
 * @note An empty query FASTA is a successful boundary case and produces an output
 *       containing only the reference record.
 */
int main(int argc, char** argv) {
    const auto total_start = Clock::now();

    try {
        Args args = parse_args(argc, argv);
        const std::string ref_file = absolute_path(args.ref);
        const std::string query_file = absolute_path(args.query);
        const std::string final_out = absolute_path(args.out);

        validate_final_output_path(final_out, ref_file, query_file);

        FastaRecord ref_record = read_single_fasta(ref_file);
        const size_t ref_len = ref_record.seq.size();
        validate_reference_symbols(ref_record, args, ref_file);

        const auto count_start = Clock::now();
        const size_t total_queries = count_fasta_sequences(query_file);
        const double count_elapsed = seconds_since(count_start);

        if (total_queries == 0) {
            std::ofstream out(final_out);
            if (!out) {
                throw std::runtime_error(
                    "Cannot create final output file: " + final_out + ". Check the parent directory and permissions."
                );
            }
            write_fasta_record(out, ref_record);
            close_final_output_or_throw(out, final_out);
            std::cout << "The query FASTA contains no sequences; only the reference sequence was written.\n";
            std::cout << "Output written to: " << final_out << "\n";
            return 0;
        }

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

        std::vector<OutputTarget> output_targets = {
            {"Final output", final_out},
            {"QC summary report", qc_summary_path},
            {"Excluded-samples report", excluded_samples_path},
            {"Unaligned-samples report", unaligned_samples_path},
            {"Abnormal-samples report", abnormal_samples_path},
            {"Failed-chunks report", failed_chunks_path},
            {"Slow-chunks report", slow_chunks_path}
        };
        if (args.enable_qc) {
            output_targets.push_back({"Filtered MSA", filtered_msa_path});
            output_targets.push_back({"Deleted-block report", deleted_blocks_path});
            output_targets.push_back({"Sample-QC report", sample_qc_summary_path});
            output_targets.push_back({"Abnormal deleted-block report", abnormal_deleted_blocks_path});
        }
        validate_distinct_output_targets(output_targets, ref_file, query_file);

        int effective_chunk_size;
        if (args.disable_auto_chunk) {
            effective_chunk_size = args.chunk_size;
        } else {
            const size_t instance_count = static_cast<size_t>(args.max_instances);
            const size_t adaptive_size = total_queries / instance_count + (total_queries % instance_count != 0 ? 1 : 0);
            if (adaptive_size > static_cast<size_t>(std::numeric_limits<int>::max())) {
                throw std::overflow_error("Automatic chunk size exceeds the supported integer range.");
            }
            const int adaptive = std::max(1, static_cast<int>(adaptive_size));
            effective_chunk_size = std::min(args.chunk_size, adaptive);
        }

        if (args.target_chunk_bases == 0 &&
            ref_len > std::numeric_limits<size_t>::max() / static_cast<size_t>(effective_chunk_size)) {
            throw std::overflow_error("Automatic target chunk base count exceeds the supported size range.");
        }
        const size_t target_chunk_bases = args.target_chunk_bases > 0
            ? args.target_chunk_bases
            : static_cast<size_t>(effective_chunk_size) * ref_len;
        const size_t effective_chunk_size_value = static_cast<size_t>(effective_chunk_size);
        const size_t rough_chunk_count = total_queries / effective_chunk_size_value +
            (total_queries % effective_chunk_size_value != 0 ? 1 : 0);
        const std::string temp_dir = prepare_temp_dir(args, final_out);

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
        if (total_queries > (std::numeric_limits<size_t>::max() - 1) / 2) {
            throw std::overflow_error("Query record count is too large for internal lookup-table sizing.");
        }
        input_rows.reserve(total_queries);
        std::unordered_map<std::string, size_t> row_index_by_seq_id;
        row_index_by_seq_id.reserve(total_queries * 2 + 1);
        std::unordered_map<std::string, size_t> seen_ids;
        seen_ids.reserve(total_queries * 2 + 1);
        std::unordered_set<std::string> used_effective_ids;
        used_effective_ids.reserve(total_queries * 2 + 1);
        used_effective_ids.insert(ref_record.id);

        std::map<std::string, SampleMeta> sample_meta;
        std::vector<DeletedBlock> all_deleted_blocks;
        std::unordered_map<std::string, QcStats> sample_qc_stats;
        std::vector<ChunkResult> failed_chunks;
        std::vector<ChunkResult> slow_chunks;

        const size_t max_pending_futures = static_cast<size_t>(args.max_instances) *
            static_cast<size_t>(args.max_pending_factor);
        std::vector<std::future<ChunkResult>> futures;

        const auto non_target_keywords = split_keywords(args.non_target_keywords);
        const auto engineered_keywords = split_keywords(args.engineered_keywords);

        size_t dispatchable_count = 0;
        size_t excluded_count = 0;
        size_t review_not_aligned_count = 0;
        size_t align_review_count = 0;

        /**
         * Integrates one completed chunk into the output stream and shared QC state.
         * Successful chunks append query alignments and mark their samples aligned;
         * failed chunks retain per-sample failure status for downstream reports.
         */
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

        /**
         * Consumes exactly one pending future. A ready future is preferred; otherwise
         * the oldest pending future is awaited to keep the queue within its bound.
         * The caller guarantees that `futures` is non-empty.
         */
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

        /**
         * Materializes and asynchronously dispatches the buffered records for a bucket.
         * Empty buffers are a no-op. Ownership of sequence identifiers is transferred
         * to the chunk task, while the reusable record buffer is cleared immediately.
         */
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
                ref_record.id,
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

        while (read_next_fasta_record(qin, rec, pending, query_file)) {
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
                used_effective_ids,
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

                if (seq_len > std::numeric_limits<size_t>::max() - buffer_bases[bucket]) {
                    throw std::overflow_error("Chunk base count exceeds the supported size range.");
                }
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
        close_final_output_or_throw(out_f, final_out);
        const double wait_elapsed = seconds_since(wait_start);

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
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    } catch (const std::ios_base::failure& e) {
        std::cerr << "Error: I/O stream failure: " << e.what() << "\n";
        return 1;
    } catch (const std::system_error& e) {
        std::cerr << "Error: System operation failed: " << e.what() << "\n";
        return 1;
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    } catch (const std::bad_alloc&) {
        std::cerr << "Error: Insufficient memory to complete the requested operation.\n";
        return 1;
    } catch (const std::length_error& e) {
        std::cerr << "Error: Requested input exceeds supported container limits: " << e.what() << "\n";
        return 1;
    }
}
