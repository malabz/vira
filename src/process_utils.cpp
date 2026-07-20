#include "process_utils.hpp"

#include <chrono>
#include <csignal>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
#include <vector>

namespace {
double elapsed_seconds(const std::chrono::steady_clock::time_point& start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}
}

CommandOutput run_process_capture(const std::vector<std::string>& args, int timeout_sec) {
    int out_pipe[2];
    int err_pipe[2];
    if (pipe(out_pipe) != 0 || pipe(err_pipe) != 0) {
        throw std::runtime_error("Failed to create pipes for the child process.");
    }

    pid_t pid = fork();
    if (pid < 0) {
        close(out_pipe[0]);
        close(out_pipe[1]);
        close(err_pipe[0]);
        close(err_pipe[1]);
        throw std::runtime_error("Failed to fork the child process.");
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
        _exit(127);
    }

    close(out_pipe[1]);
    close(err_pipe[1]);

    std::string stdout_data;
    std::string stderr_data;

    std::thread t_out([&]() {
        char buf[65536];
        ssize_t n;
        while ((n = read(out_pipe[0], buf, sizeof(buf))) > 0) {
            stdout_data.append(buf, static_cast<size_t>(n));
        }
        close(out_pipe[0]);
    });

    std::thread t_err([&]() {
        char buf[65536];
        ssize_t n;
        while ((n = read(err_pipe[0], buf, sizeof(buf))) > 0) {
            stderr_data.append(buf, static_cast<size_t>(n));
        }
        close(err_pipe[0]);
    });

    const auto start = std::chrono::steady_clock::now();
    int status = 0;
    bool timed_out = false;

    while (true) {
        const pid_t rc = waitpid(pid, &status, WNOHANG);
        if (rc == pid) break;
        if (rc < 0) {
            status = 0;
            break;
        }

        if (timeout_sec > 0 && elapsed_seconds(start) >= static_cast<double>(timeout_sec)) {
            timed_out = true;
            kill(-pid, SIGTERM);
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
            kill(-pid, SIGKILL);
            waitpid(pid, &status, 0);
            break;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    t_out.join();
    t_err.join();

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
