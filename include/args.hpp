/**
 * @file args.hpp
 * @brief Command-line parsing and usage-text interfaces.
 */
#pragma once
#include "common.hpp"

/**
 * @brief Parses and validates command-line arguments.
 *
 * Required input, output, parallelism, chunking, and QC options are converted
 * into an Args value. Numeric values must be complete, finite, in range, and
 * consistent with the CLI threshold relationships.
 *
 * @param argc Number of entries in @p argv, including the program name.
 * @param argv Null-terminated argument strings supplied by the process.
 * @return A fully populated and validated argument set.
 * @throws std::invalid_argument If a required option or value is missing, an
 *         option is unknown, a numeric value is malformed or out of range, or
 *         a threshold combination is invalid.
 * @note Passing `-h` or `--help` prints usage and terminates the process with
 *       a successful exit status rather than returning.
 */
Args parse_args(int argc, char** argv);

/**
 * @brief Writes the command-line usage reference to standard output.
 * @param prog Program name to display in usage examples; must point to a valid
 *        null-terminated string.
 */
void print_usage(const char* prog);
