/**
 * @file path_utils.cpp
 * @brief Implements POSIX path handling and temporary-directory naming helpers.
 *
 * These helpers deliberately avoid filesystem-library dependencies. Path
 * normalization is best effort: unresolved output paths are made relative to
 * the current working directory, while existing paths are canonicalized with
 * `realpath`.
 */
#include "path_utils.hpp"

#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

std::string path_basename(const std::string& path) {
    size_t pos = path.find_last_of('/');
    return (pos == std::string::npos) ? path : path.substr(pos + 1);
}

std::string path_dirname(const std::string& path) {
    size_t pos = path.find_last_of('/');
    if (pos == std::string::npos) return ".";
    if (pos == 0) return "/";
    return path.substr(0, pos);
}

std::string path_without_extension(const std::string& path) {
    std::string dir = path_dirname(path);
    std::string base = path_basename(path);
    size_t dot = base.find_last_of('.');
    std::string stem = (dot == std::string::npos) ? base : base.substr(0, dot);
    return (dir == ".") ? stem : (dir + "/" + stem);
}

std::string join_path(const std::string& a, const std::string& b) {
    if (a.empty() || a == ".") return b;
    if (a.back() == '/') return a + b;
    return a + "/" + b;
}

bool file_exists(const std::string& path) {
    return access(path.c_str(), F_OK) == 0;
}

std::string absolute_path(const std::string& path) {
    char* resolved = realpath(path.c_str(), nullptr);
    if (resolved) {
        std::string out(resolved);
        free(resolved);
        return out;
    }
    if (!path.empty() && path[0] == '/') return path;
    char cwd[4096];
    if (!getcwd(cwd, sizeof(cwd))) return path;
    return join_path(cwd, path);
}

namespace {
/**
 * @brief Creates one directory with POSIX mode 0755.
 * @param path Directory path to create; parent directories must already exist.
 * @throws std::runtime_error If `mkdir` fails, including when @p path already
 *         exists or its parent is unavailable.
 */
void make_dir_or_throw(const std::string& path) {
    if (mkdir(path.c_str(), 0755) != 0) {
        throw std::runtime_error("Failed to create temporary directory: " + path + " (" + std::strerror(errno) + ")");
    }
}

/**
 * @brief Generates a zero-padded, twelve-digit decimal suffix.
 * @return A string in the inclusive lexical range `000000000000` through
 *         `999999999999`.
 * @note Randomness reduces name collisions but does not provide a security or
 *       uniqueness guarantee; directory creation remains the final check.
 */
std::string numeric_suffix_12() {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<unsigned long long> dist(0, 999999999999ULL);
    std::ostringstream oss;
    oss << std::setw(12) << std::setfill('0') << dist(gen);
    return oss.str();
}
}

std::string prepare_temp_dir(const Args& args, const std::string& final_out) {
    const std::string suffix = numeric_suffix_12();
    std::string dir;
    if (args.temp_dir_provided) dir = args.temp_dir + "." + suffix;
    else dir = path_basename(final_out) + ".tmp_chunks." + suffix;
    make_dir_or_throw(dir);
    return absolute_path(dir);
}

std::string output_prefix(const std::string& final_out) {
    return path_without_extension(final_out);
}

std::string shell_quote(const std::string& s) {
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') out += "'\\''";
        else out += c;
    }
    out += "'";
    return out;
}
