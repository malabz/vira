/**
 * @file process_utils.hpp
 * @brief POSIX child-process execution and captured-output interfaces.
 */
#pragma once

#include <string>
#include <vector>

/** @brief Captured termination state and output for one external command. */
struct CommandOutput {
    int exit_code = -1;       ///< Exit status, `128 + signal`, `124` on timeout, or `-1` if unavailable.
    std::string out;          ///< Bytes captured from standard output.
    std::string err;          ///< Bytes captured from standard error, including a timeout annotation.
    bool timed_out = false;   ///< Whether the timeout path terminated the child process group.
    double elapsed_sec = 0.0; ///< Wall-clock duration from launch monitoring through termination.
};

/**
 * @brief Executes a command without a shell and captures stdout and stderr.
 *
 * The child is placed in its own process group so timeout handling can signal
 * its descendants. Output pipes are drained concurrently to avoid deadlock.
 *
 * @param args Argument vector whose first element names the executable.
 * @param timeout_sec Timeout in seconds; zero or a negative value disables it.
 * @return Captured output, normalized exit status, timeout flag, and duration.
 * @throws std::invalid_argument If @p args has no non-empty executable name.
 * @throws std::runtime_error If process setup, waiting, or captured pipe I/O fails.
 * @note This interface requires POSIX process, pipe, signal, and wait APIs.
 */
CommandOutput run_process_capture(const std::vector<std::string>& args, int timeout_sec = 0);
