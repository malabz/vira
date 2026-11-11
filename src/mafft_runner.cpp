/**
 * @file mafft_runner.cpp
 * @brief Implements execution and result collection for one MAFFT query chunk.
 *
 * The runner translates a chunk request into a MAFFT command, captures the
 * child process output, optionally parses `--mapout` quality data, and records
 * process failures in a ChunkResult. Per-chunk query and map files are removed
 * on the normal return path regardless of the MAFFT exit status. Known runner
 * failures are logged with the chunk label and converted into explicit result
 * categories.
 */
#include "mafft_runner.hpp"

#include "fasta_io.hpp"
#include "path_utils.hpp"
#include "process_utils.hpp"
#include "qc.hpp"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <utility>
#include <vector>

namespace {
class MafftOutputFormatError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

std::string format_identifier_set(const std::set<std::string>& identifiers) {
    std::ostringstream output;
    output << '[';
    bool first = true;
    for (const auto& identifier : identifiers) {
        if (!first) output << ',';
        output << identifier;
        first = false;
    }
    output << ']';
    return output.str();
}

void validate_mafft_fasta_output(
    const std::string& output,
    const std::string& expected_ref_id,
    const std::vector<std::string>& expected_query_ids,
    const std::string& prefix) {

    if (expected_ref_id.empty()) {
        throw MafftOutputFormatError(
            "MAFFT chunk '" + prefix + "' has an empty expected reference identifier."
        );
    }

    std::set<std::string> expected_queries;
    std::set<std::string> duplicate_expected_queries;
    for (const auto& identifier : expected_query_ids) {
        if (identifier.empty()) {
            throw MafftOutputFormatError(
                "MAFFT chunk '" + prefix + "' has an empty expected query identifier."
            );
        }
        if (identifier == expected_ref_id) {
            throw MafftOutputFormatError(
                "MAFFT chunk '" + prefix + "' has a query identifier that conflicts with the reference identifier '" +
                expected_ref_id + "'."
            );
        }
        if (!expected_queries.insert(identifier).second) duplicate_expected_queries.insert(identifier);
    }
    if (!duplicate_expected_queries.empty()) {
        throw MafftOutputFormatError(
            "MAFFT chunk '" + prefix + "' has duplicate expected query identifiers: duplicates=" +
            format_identifier_set(duplicate_expected_queries) + "."
        );
    }

    std::istringstream in(output);
    std::string line;
    size_t header_count = 0;
    bool seen_header = false;
    bool current_record_has_sequence = false;
    std::vector<std::string> actual_record_ids;

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.find('\0') != std::string::npos) {
            throw MafftOutputFormatError(
                "MAFFT chunk '" + prefix + "' returned binary or corrupted FASTA output containing a NUL byte."
            );
        }
        const bool blank = std::all_of(line.begin(), line.end(), [](unsigned char c) {
            return std::isspace(c) != 0;
        });
        if (blank) continue;
        if (line[0] == '>') {
            if (seen_header && !current_record_has_sequence) {
                throw MafftOutputFormatError(
                    "MAFFT chunk '" + prefix + "' returned a FASTA record with an empty sequence."
                );
            }
            const std::string header = line.substr(1);
            const bool empty_header = std::all_of(header.begin(), header.end(), [](unsigned char c) {
                return std::isspace(c) != 0;
            });
            if (header.empty() || empty_header) {
                throw MafftOutputFormatError(
                    "MAFFT chunk '" + prefix + "' returned FASTA output with an empty header."
                );
            }
            const std::string identifier = fasta_primary_id(header);
            if (identifier.empty()) {
                throw MafftOutputFormatError(
                    "MAFFT chunk '" + prefix + "' returned FASTA output with an empty primary identifier."
                );
            }
            actual_record_ids.push_back(identifier);
            ++header_count;
            seen_header = true;
            current_record_has_sequence = false;
            continue;
        }
        if (!seen_header) {
            throw MafftOutputFormatError(
                "MAFFT chunk '" + prefix + "' returned non-FASTA output before the first header."
            );
        }
        current_record_has_sequence = true;
    }

    if (seen_header && !current_record_has_sequence) {
        throw MafftOutputFormatError(
            "MAFFT chunk '" + prefix + "' returned a FASTA record with an empty sequence."
        );
    }
    if (expected_query_ids.size() == std::numeric_limits<size_t>::max()) {
        throw MafftOutputFormatError("MAFFT output record count exceeds the supported size range.");
    }
    const size_t expected_headers = expected_query_ids.size() + 1;
    if (header_count != expected_headers) {
        throw MafftOutputFormatError(
            "MAFFT chunk '" + prefix + "' returned " + std::to_string(header_count) +
            " FASTA records; expected " + std::to_string(expected_headers) +
            " (one reference plus all query records)."
        );
    }

    if (actual_record_ids.empty() || actual_record_ids.front() != expected_ref_id) {
        const std::string actual = actual_record_ids.empty() ? "<missing>" : actual_record_ids.front();
        throw MafftOutputFormatError(
            "MAFFT chunk '" + prefix + "' returned reference identifier '" + actual +
            "' in the first record; expected '" + expected_ref_id + "'."
        );
    }

    std::set<std::string> actual_queries;
    std::set<std::string> duplicate_actual_queries;
    for (size_t i = 1; i < actual_record_ids.size(); ++i) {
        if (!actual_queries.insert(actual_record_ids[i]).second) {
            duplicate_actual_queries.insert(actual_record_ids[i]);
        }
    }

    std::set<std::string> missing;
    std::set<std::string> unexpected;
    std::set_difference(
        expected_queries.begin(), expected_queries.end(),
        actual_queries.begin(), actual_queries.end(),
        std::inserter(missing, missing.end())
    );
    std::set_difference(
        actual_queries.begin(), actual_queries.end(),
        expected_queries.begin(), expected_queries.end(),
        std::inserter(unexpected, unexpected.end())
    );
    if (!missing.empty() || !unexpected.empty() || !duplicate_actual_queries.empty()) {
        throw MafftOutputFormatError(
            "MAFFT chunk '" + prefix + "' returned a query identifier mismatch: missing=" +
            format_identifier_set(missing) + ", unexpected=" + format_identifier_set(unexpected) +
            ", duplicates=" + format_identifier_set(duplicate_actual_queries) + "."
        );
    }
}

/**
 * @brief Persists failure diagnostics and records the path only after a verified write.
 * @param result Chunk result receiving a valid diagnostic path.
 * @param temp_dir Directory for diagnostic output.
 * @param prefix Stable chunk label used in the filename and warnings.
 * @param content Diagnostic content to write.
 */
void save_failure_diagnostic(
    ChunkResult& result,
    const std::string& temp_dir,
    const std::string& prefix,
    const std::string& content) {

    const std::string err_path = join_path(temp_dir, prefix + ".stderr.txt");
    std::ofstream ferr(err_path, std::ios::out | std::ios::trunc);
    if (!ferr) {
        std::cerr << "Warning: Could not open stderr diagnostics for MAFFT chunk '" << prefix
                  << "' at: " << err_path << "\n";
        return;
    }

    ferr << content;
    ferr.close();
    if (!ferr) {
        std::cerr << "Warning: Could not write complete stderr diagnostics for MAFFT chunk '" << prefix
                  << "' at: " << err_path << "\n";
        return;
    }
    result.err_path = err_path;
}

void remove_temporary_file(const std::string& path, const std::string& description, const std::string& prefix) {
    if (std::remove(path.c_str()) == 0) return;
    const int remove_error = errno;
    if (remove_error == ENOENT) return;
    std::cerr << "Warning: Could not remove " << description << " for MAFFT chunk '" << prefix
              << "': " << path << " (" << std::strerror(remove_error) << ").\n";
}

/**
 * @brief Converts one known exception into a contextual failed-chunk result.
 */
void record_runner_failure(
    ChunkResult& result,
    const std::string& failure_reason,
    const std::string& message,
    const std::string& prefix,
    const std::string& temp_dir,
    bool save_failed_stderr) {

    result.success = false;
    result.stdout_data.clear();
    result.deleted_blocks.clear();
    result.qc_sample_stats.clear();
    result.failure_reason = failure_reason;
    std::cerr << "Error: MAFFT chunk '" << prefix << "' failed [" << failure_reason << "]: "
              << message << "\n";
    if (save_failed_stderr) save_failure_diagnostic(result, temp_dir, prefix, message + "\n");
}
}

/**
 * @brief Executes MAFFT for one prepared chunk and packages its observable result.
 *
 * Process exit failures and timeouts are represented in the returned
 * ChunkResult rather than propagated as exceptions. Known mapout, system, and
 * runtime failures are converted to distinct failure categories. When requested,
 * diagnostic stderr is written below @p temp_dir and its path is reported only
 * after a successful write. Map data is parsed only after a successful process
 * exit; a requested but missing map file is an explicit failure.
 *
 * @param chunk_id Stable identifier assigned by the dispatcher.
 * @param bucket_name Logical bucket recorded in reports and diagnostic names.
 * @param chunk_query_path Existing FASTA file containing the query chunk.
 * @param ref_path Existing reference FASTA passed to MAFFT.
 * @param expected_ref_id Primary reference identifier expected as the first
 *        MAFFT output record.
 * @param threads_per_mafft Thread count forwarded to MAFFT.
 * @param temp_dir Directory used for optional failure diagnostics.
 * @param save_failed_stderr Whether failed-process stderr should be persisted.
 * @param enable_qc Whether `--mapout` should be requested and parsed.
 * @param mafft_timeout_sec Timeout in seconds; zero disables the timeout.
 * @param slow_chunk_sec Elapsed-time threshold for marking the chunk slow;
 *        non-positive values disable slow marking.
 * @param seq_ids Query identifiers contained in the chunk; ownership is moved
 *        into the result.
 * @param total_bases Sum of query bases in the chunk.
 * @return A complete ChunkResult containing status, captured output, timing,
 *         identifiers, and optional post-alignment QC data.
 * @note The query chunk and optional map file are removed before normal return.
 */
ChunkResult run_mafft_chunk(
    int chunk_id,
    const std::string& bucket_name,
    const std::string& chunk_query_path,
    const std::string& ref_path,
    const std::string& expected_ref_id,
    int threads_per_mafft,
    const std::string& temp_dir,
    bool save_failed_stderr,
    bool enable_qc,
    int mafft_timeout_sec,
    double slow_chunk_sec,
    std::vector<std::string> seq_ids,
    size_t total_bases) {

    ChunkResult result;
    result.bucket_name = bucket_name;
    result.chunk_id = chunk_id;
    result.seq_ids = std::move(seq_ids);
    result.record_count = result.seq_ids.size();
    result.total_bases = total_bases;

    std::ostringstream prefix_oss;
    prefix_oss << bucket_name << ".chunk_" << std::setw(8) << std::setfill('0') << chunk_id;
    const std::string prefix = prefix_oss.str();

    std::vector<std::string> cmd = {"mafft", "--anysymbol"};
    if (enable_qc) cmd.push_back("--mapout");
    cmd.insert(cmd.end(), {
        "--6merpair",
        "--thread", std::to_string(threads_per_mafft),
        "--keeplength",
        "--addfragments", chunk_query_path,
        ref_path
    });

    const std::string map_path = chunk_query_path + ".map";

    try {
        CommandOutput co = run_process_capture(cmd, mafft_timeout_sec);
        result.elapsed_sec = co.elapsed_sec;
        result.timed_out = co.timed_out;
        result.slow = slow_chunk_sec > 0.0 && co.elapsed_sec >= slow_chunk_sec;

        if (co.exit_code == 0 && !co.timed_out) {
            validate_mafft_fasta_output(co.out, expected_ref_id, result.seq_ids, prefix);

            std::vector<DeletedBlock> validated_deleted_blocks;
            std::unordered_map<std::string, QcStats> validated_qc_sample_stats;

            if (enable_qc) {
                if (!file_exists(map_path)) {
                    throw MapoutIoError("MAFFT mapout file was not created: " + map_path + ".");
                }
                auto parsed = parse_map_to_deleted_blocks(map_path, result.seq_ids);
                validated_deleted_blocks = std::move(parsed.first);
                validated_qc_sample_stats = std::move(parsed.second);
            }

            result.stdout_data = std::move(co.out);
            result.deleted_blocks = std::move(validated_deleted_blocks);
            result.qc_sample_stats = std::move(validated_qc_sample_stats);
            result.success = true;
        } else {
            result.success = false;
            if (co.timed_out) {
                result.failure_reason = "TIMEOUT";
            } else if (co.exit_code == 127) {
                result.failure_reason = "MAFFT_EXEC_NOT_FOUND";
            } else if (co.exit_code == 126) {
                result.failure_reason = "MAFFT_EXEC_NOT_EXECUTABLE";
            } else {
                result.failure_reason = "MAFFT_EXIT_NONZERO";
            }
            std::cerr << "Error: MAFFT chunk '" << prefix << "' failed [" << result.failure_reason
                      << "] with exit code " << co.exit_code << ".\n";
            if (save_failed_stderr) save_failure_diagnostic(result, temp_dir, prefix, co.err);
        }
    } catch (const MafftOutputFormatError& e) {
        record_runner_failure(result, "MAFFT_OUTPUT_FORMAT_ERROR", e.what(), prefix, temp_dir, save_failed_stderr);
    } catch (const MapoutFormatError& e) {
        record_runner_failure(result, "MAPOUT_FORMAT_ERROR", e.what(), prefix, temp_dir, save_failed_stderr);
    } catch (const MapoutIoError& e) {
        record_runner_failure(result, "MAPOUT_IO_ERROR", e.what(), prefix, temp_dir, save_failed_stderr);
    } catch (const std::invalid_argument& e) {
        record_runner_failure(result, "RUNNER_INVALID_ARGUMENT", e.what(), prefix, temp_dir, save_failed_stderr);
    } catch (const std::system_error& e) {
        record_runner_failure(result, "RUNNER_SYSTEM_ERROR", e.what(), prefix, temp_dir, save_failed_stderr);
    } catch (const std::runtime_error& e) {
        record_runner_failure(result, "RUNNER_IO_ERROR", e.what(), prefix, temp_dir, save_failed_stderr);
    }

    remove_temporary_file(chunk_query_path, "temporary query file", prefix);
    if (enable_qc) remove_temporary_file(map_path, "temporary mapout file", prefix);
    return result;
}
