/**
 * @file process_utils.cpp
 * @brief Implements POSIX child-process execution with output capture and timeout control.
 *
 * The child is placed in a dedicated process group so a timeout can terminate
 * MAFFT and any descendants. Separate reader threads drain stdout and stderr to
 * prevent either pipe from blocking the child when output is large.
 */
#include "process_utils.hpp"

#include <cerrno>
#include <chrono>
#include <csignal>
#include <cstddef>
#include <cstring>
#include <new>
#include <stdexcept>
#include <string>
#include <system_error>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
#include <vector>

namespace {
/**
 * @brief Measures monotonic elapsed wall-clock time from @p start.
 * @param start Earlier steady-clock time point.
 * @return Elapsed seconds as a floating-point value.
 */
double elapsed_seconds(const std::chrono::steady_clock::time_point& start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}
}

/**
 * @brief Runs one executable without a shell and captures both output streams.
 *
 * @param args Argument vector whose first element names the executable. The
 *        vector must be non-empty; entries are passed directly to `execvp`.
 * @param timeout_sec Maximum runtime in seconds. Values less than or equal to
 *        zero disable timeout enforcement.
 * @return Exit status, captured stdout/stderr, elapsed time, and timeout state.
 *         A timeout is reported with exit code 124; signal termination is
 *         reported as `128 + signal`.
 * @throws std::invalid_argument If no executable name is supplied.
 * @throws std::runtime_error If process setup, waiting, or pipe I/O fails.
 * @note On timeout, SIGTERM is sent to the process group, followed by SIGKILL
 *       after a 500 ms grace period.
 */
CommandOutput run_process_capture(const std::vector<std::string>& args, int timeout_sec) {
    if (args.empty() || args.front().empty()) {
        throw std::invalid_argument("Cannot execute an empty command: an executable name is required.");
    }

    int out_pipe[2];
    int err_pipe[2];
    if (pipe(out_pipe) != 0) {
        throw std::runtime_error("Failed to create the child stdout pipe: " + std::string(std::strerror(errno)));
    }
    if (pipe(err_pipe) != 0) {
        const int pipe_error = errno;
        close(out_pipe[0]);
        close(out_pipe[1]);
        throw std::runtime_error("Failed to create the child stderr pipe: " + std::string(std::strerror(pipe_error)));
    }

    pid_t pid = fork();
    if (pid < 0) {
        const int fork_error = errno;
        close(out_pipe[0]);
        close(out_pipe[1]);
        close(err_pipe[0]);
        close(err_pipe[1]);
        throw std::runtime_error("Failed to fork the child process for '" + args.front() + "': " + std::string(std::strerror(fork_error)));
    }

    if (pid == 0) {
        // Put the child in its own process group so timeout cleanup can terminate descendants.
        setpgid(0, 0);
        dup2(out_pipe[1], STDOUT_FILENO);
        dup2(err_pipe[1], STDERR_FILENO);
        close(out_pipe[0]);
        close(out_pipe[1]);
        close(err_pipe[0]);
        close(err_pipe[1]);

        std::vector<char*> argv;
        argv.reserve(args.size() + 1);
        for (const auto& s : args) argv.push_back(const_cast<char*>(s.c_str()));
        argv.push_back(nullptr);

        execvp(argv[0], argv.data());
        const int exec_error = errno;
        if (exec_error == ENOENT || exec_error == ENOTDIR) {
            static constexpr char not_found[] =
                "[vira] Failed to execute child process: executable was not found.\n";
            (void)write(STDERR_FILENO, not_found, sizeof(not_found) - 1);
            _exit(127);
        }
        static constexpr char not_executable[] =
            "[vira] Failed to execute child process: executable could not be executed.\n";
        (void)write(STDERR_FILENO, not_executable, sizeof(not_executable) - 1);
        _exit(126);
    }

    close(out_pipe[1]);
    close(err_pipe[1]);

    std::string stdout_data;
    std::string stderr_data;
    int stdout_read_error = 0;
    int stderr_read_error = 0;
    int stdout_capture_failure = 0;
    int stderr_capture_failure = 0;

    auto read_stdout = [&]() {
        char buf[65536];
        try {
            while (true) {
                const ssize_t n = read(out_pipe[0], buf, sizeof(buf));
                if (n > 0) {
                    stdout_data.append(buf, static_cast<size_t>(n));
                    continue;
                }
                if (n == 0) break;
                if (errno == EINTR) continue;
                stdout_read_error = errno;
                break;
            }
        } catch (const std::bad_alloc&) {
            stdout_capture_failure = 1;
        } catch (const std::length_error&) {
            stdout_capture_failure = 2;
        }
        close(out_pipe[0]);
    };

    auto read_stderr = [&]() {
        char buf[65536];
        try {
            while (true) {
                const ssize_t n = read(err_pipe[0], buf, sizeof(buf));
                if (n > 0) {
                    stderr_data.append(buf, static_cast<size_t>(n));
                    continue;
                }
                if (n == 0) break;
                if (errno == EINTR) continue;
                stderr_read_error = errno;
                break;
            }
        } catch (const std::bad_alloc&) {
            stderr_capture_failure = 1;
        } catch (const std::length_error&) {
            stderr_capture_failure = 2;
        }
        close(err_pipe[0]);
    };

    std::thread t_out;
    std::thread t_err;
    try {
        t_out = std::thread(read_stdout);
        t_err = std::thread(read_stderr);
    } catch (const std::system_error& e) {
        kill(-pid, SIGKILL);
        kill(pid, SIGKILL);
        if (!t_out.joinable()) close(out_pipe[0]);
        if (!t_err.joinable()) close(err_pipe[0]);
        while (waitpid(pid, nullptr, 0) < 0 && errno == EINTR) {
        }
        if (t_out.joinable()) t_out.join();
        if (t_err.joinable()) t_err.join();
        throw std::runtime_error(
            "Failed to start output-capture threads for child process '" + args.front() + "': " + e.what()
        );
    }

    const auto start = std::chrono::steady_clock::now();
    int status = 0;
    bool timed_out = false;
    int wait_error = 0;

    while (true) {
        const pid_t rc = waitpid(pid, &status, WNOHANG);
        if (rc == pid) break;
        if (rc < 0) {
            if (errno == EINTR) continue;
            wait_error = errno;
            kill(-pid, SIGKILL);
            kill(pid, SIGKILL);
            while (waitpid(pid, &status, 0) < 0 && errno == EINTR) {
            }
            break;
        }

        if (timeout_sec > 0 && elapsed_seconds(start) >= static_cast<double>(timeout_sec)) {
            timed_out = true;
            kill(-pid, SIGTERM);
            kill(pid, SIGTERM);
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
            kill(-pid, SIGKILL);
            kill(pid, SIGKILL);
            while (waitpid(pid, &status, 0) < 0) {
                if (errno == EINTR) continue;
                wait_error = errno;
                break;
            }
            break;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    t_out.join();
    t_err.join();

    if (wait_error != 0) {
        throw std::runtime_error("Failed while waiting for child process '" + args.front() + "': " + std::string(std::strerror(wait_error)));
    }
    if (stdout_read_error != 0) {
        throw std::runtime_error("Failed while reading stdout from child process '" + args.front() + "': " + std::string(std::strerror(stdout_read_error)));
    }
    if (stderr_read_error != 0) {
        throw std::runtime_error("Failed while reading stderr from child process '" + args.front() + "': " + std::string(std::strerror(stderr_read_error)));
    }
    if (stdout_capture_failure != 0) {
        throw std::runtime_error(
            "Cannot capture stdout from child process '" + args.front() +
            (stdout_capture_failure == 1 ? "': insufficient memory." : "': output exceeds supported string limits.")
        );
    }
    if (stderr_capture_failure != 0) {
        throw std::runtime_error(
            "Cannot capture stderr from child process '" + args.front() +
            (stderr_capture_failure == 1 ? "': insufficient memory." : "': output exceeds supported string limits.")
        );
    }

    CommandOutput co;
    co.timed_out = timed_out;
    co.elapsed_sec = elapsed_seconds(start);
    if (timed_out) {
        co.exit_code = 124;
        stderr_data += "\n[vira] MAFFT chunk timed out.\n";
    } else if (WIFEXITED(status)) {
        co.exit_code = WEXITSTATUS(status);
    } else if (WIFSIGNALED(status)) {
        co.exit_code = 128 + WTERMSIG(status);
    }
    co.out = std::move(stdout_data);
    co.err = std::move(stderr_data);
    return co;
}
