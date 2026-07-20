#pragma once

#include <string>
#include <vector>

struct CommandOutput {
    int exit_code = -1;
    std::string out;
    std::string err;
    bool timed_out = false;
    double elapsed_sec = 0.0;
};

CommandOutput run_process_capture(const std::vector<std::string>& args, int timeout_sec = 0);
