#include "mafft_runner.hpp"

#include "path_utils.hpp"
#include "process_utils.hpp"
#include "qc.hpp"

#include <cstdio>
#include <exception>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <utility>
#include <vector>

ChunkResult run_mafft_chunk(
    int chunk_id,
    const std::string& bucket_name,
    const std::string& chunk_query_path,
    const std::string& ref_path,
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
            result.stdout_data = std::move(co.out);
            result.success = true;

            if (enable_qc && file_exists(map_path)) {
                auto parsed = parse_map_to_deleted_blocks(map_path);
                result.deleted_blocks = std::move(parsed.first);
                result.qc_sample_stats = std::move(parsed.second);
            }
        } else {
            result.success = false;
            result.failure_reason = co.timed_out ? "TIMEOUT" : "MAFFT_EXIT_NONZERO";
            if (save_failed_stderr) {
                const std::string err_path = join_path(temp_dir, prefix + ".stderr.txt");
                std::ofstream ferr(err_path);
                ferr << co.err;
                result.err_path = err_path;
            }
        }
    } catch (const std::exception& e) {
        result.success = false;
        result.failure_reason = "RUNNER_EXCEPTION";
        if (save_failed_stderr) {
            const std::string err_path = join_path(temp_dir, prefix + ".stderr.txt");
            std::ofstream ferr(err_path);
            ferr << e.what() << '\n';
            result.err_path = err_path;
        }
    }

    std::remove(chunk_query_path.c_str());
    if (enable_qc) std::remove(map_path.c_str());
    return result;
}
