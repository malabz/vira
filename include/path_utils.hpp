/**
 * @file path_utils.hpp
 * @brief POSIX-oriented path, temporary-directory, and shell-quoting helpers.
 */
#pragma once
#include "common.hpp"

#include <string>

/**
 * @brief Joins two path fragments with one `/` when needed.
 * @param a Directory-like prefix; empty and `.` return @p b unchanged.
 * @param b Suffix to append.
 * @return Lexically joined path without normalization or existence checks.
 */
std::string join_path(const std::string& a, const std::string& b);

/**
 * @brief Extracts text after the final `/`.
 * @param path Input path.
 * @return Basename text; a trailing slash yields an empty string.
 */
std::string path_basename(const std::string& path);

/**
 * @brief Extracts the directory portion before the final `/`.
 * @param path Input path.
 * @return `.`, `/`, or the lexical directory portion as appropriate.
 */
std::string path_dirname(const std::string& path);

/**
 * @brief Removes the final filename extension while preserving its directory.
 * @param path Input path.
 * @return Path with text from the basename's final `.` onward removed.
 * @note No filesystem access or special handling for dotfiles is performed.
 */
std::string path_without_extension(const std::string& path);

/**
 * @brief Tests whether a filesystem entry is visible to the process.
 * @param path Entry to test.
 * @return `true` when POSIX `access(F_OK)` succeeds; directories also qualify.
 */
bool file_exists(const std::string& path);

/**
 * @brief Resolves a path to an absolute spelling when possible.
 * @param path Existing or prospective path.
 * @return Canonical `realpath` output for an existing entry; otherwise the
 *         absolute input, a lexical current-directory join, or the original
 *         path if the current directory cannot be read.
 * @note A nonexistent result is not canonicalized and may still contain `.` or `..`.
 */
std::string absolute_path(const std::string& path);

/**
 * @brief Creates a randomly suffixed per-run temporary directory.
 * @param args Configuration controlling the optional base directory name.
 * @param final_out Final MSA path used to derive the default directory name.
 * @return Absolute path to the newly created directory.
 * @throws std::runtime_error If the directory cannot be created.
 * @note A random 12-digit suffix is always appended; existing directories are
 *       never reused.
 */
std::string prepare_temp_dir(const Args& args, const std::string& final_out);

/**
 * @brief Derives the prefix used by auxiliary output files.
 * @param final_out Final MSA path.
 * @return @p final_out without the basename's final extension.
 */
std::string output_prefix(const std::string& final_out);

/**
 * @brief Produces a single-quoted POSIX shell word.
 * @param s Arbitrary argument text, including embedded single quotes.
 * @return Quoted text suitable for interpolation as one shell argument.
 */
std::string shell_quote(const std::string& s);
