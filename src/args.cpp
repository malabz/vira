/**
 * @file args.cpp
 * @brief Implements ViRA command-line help, option parsing, and configuration validation.
 *
 * Required paths are collected without touching the filesystem. Numeric and policy
 * constraints are normalized or rejected before the main workflow starts, while a
 * small set of deprecated options remains accepted for command-line compatibility.
 */

#include "args.hpp"

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>

void print_usage(const char* prog) {
    std::printf("Usage:\n");
    std::printf("  %s -r REF.fasta -q QUERY.fasta -o OUT.fasta [options]\n\n", prog);
    std::printf("Required arguments:\n");
    std::printf("  -r, --ref FILE                       Reference FASTA containing exactly one sequence\n");
    std::printf("  -q, --query FILE                     Query FASTA file\n");
    std::printf("  -o, --out FILE                       Final full MSA FASTA output\n\n");
    std::printf("Parallelism and chunking:\n");
    std::printf("  --max-instances INT                  Concurrent MAFFT processes [16]\n");
    std::printf("  --threads-per-mafft INT              Threads per MAFFT process [1]\n");
    std::printf("  --chunk-size INT                     Maximum records per chunk [1000]\n");
    std::printf("  --target-chunk-bases INT             Maximum total bases per chunk [auto]\n");
    std::printf("  --disable-auto-chunk                 Disable automatic chunk-size reduction\n");
    std::printf("  --max-pending-factor INT             Max queued tasks = max-instances * factor [2]\n\n");
    std::printf("MAFFT runtime controls:\n");
    std::printf("  --mafft-timeout-sec INT              Per-chunk timeout; 0 disables timeout [0]\n");
    std::printf("  --slow-chunk-sec FLOAT               Report chunks taking at least this many seconds [300]\n\n");
    std::printf("Temporary files:\n");
    std::printf("  --temp-dir DIR                       Base temp dir; random numeric suffix is appended\n");
    std::printf("  --cleanup                            Remove the temporary directory after completion\n\n");
    std::printf("Layer 1: input QC:\n");
    std::printf("  --duplicate-id-policy MODE           exclude|rename [exclude]\n");
    std::printf("  --allowed-seq-chars STRING           Allowed IUPAC-like sequence symbols [ACGTURYSWKMBDHVN-.]\n");
    std::printf("  --review-n-ratio FLOAT               Mark for review if N ratio >= value [0.05]\n");
    std::printf("  --exclude-n-ratio FLOAT              Exclude if N ratio >= value [0.50]\n");
    std::printf("  --exclude-min-len-ratio FLOAT        Exclude if seq/ref length ratio < value [0.10]\n");
    std::printf("  --exclude-max-len-ratio FLOAT        Exclude if seq/ref length ratio > value [5.00]\n\n");
    std::printf("Layer 2: virology suitability:\n");
    std::printf("  --reasonable-min-len-ratio FLOAT     Review if seq/ref length ratio < value [0.85]\n");
    std::printf("  --reasonable-max-len-ratio FLOAT     Review if seq/ref length ratio > value [1.15]\n");
    std::printf("  --non-target-keywords CSV            Header keywords suggesting non-target species\n");
    std::printf("  --engineered-keywords CSV            Header keywords suggesting engineered constructs\n");
    std::printf("  --align-review-samples               Also align suspected non-target/engineered records\n\n");
    std::printf("MAFFT --mapout post-alignment QC:\n");
    std::printf("  --enable-qc                          Enable QC based on MAFFT --mapout\n");
    std::printf("  --warn-lost-ratio FLOAT              WARN if lost_ratio >= value [0.001]\n");
    std::printf("  --fail-lost-ratio FLOAT              FAIL if lost_ratio >= value [0.01]\n");
    std::printf("  --warn-max-block-len INT             WARN if max_lost_block_len >= value [20]\n");
    std::printf("  --fail-max-block-len INT             FAIL if max_lost_block_len >= value [100]\n");
    std::printf("  --save-failed-stderr                 Save stderr for failed MAFFT chunks\n\n");
    std::printf("Output naming:\n");
    std::printf("  Auxiliary files use the --out prefix:\n");
    std::printf("    PREFIX.qc_summary.tsv\n");
    std::printf("    PREFIX.excluded_samples.tsv\n");
    std::printf("    PREFIX.unaligned_samples.tsv\n");
    std::printf("    PREFIX.abnormal_samples.tsv\n");
    std::printf("    PREFIX.failed_chunks.tsv\n");
    std::printf("    PREFIX.slow_chunks.tsv\n");
    std::printf("  With --enable-qc, additional files are written:\n");
    std::printf("    PREFIX.filtered.fasta\n");
    std::printf("    PREFIX.deleted_blocks.tsv\n");
    std::printf("    PREFIX.sample_qc_summary.tsv\n");
    std::printf("    PREFIX.abnormal_deleted_blocks.tsv\n\n");
    std::printf("Examples:\n");
    std::printf("  %s -r ref.fasta -q query.fasta -o result.fasta\n", prog);
    std::printf("  %s -r ref.fasta -q query.fasta -o result.fasta --enable-qc --mafft-timeout-sec 1800\n\n", prog);
    std::printf("Other:\n");
    std::printf("  -h, --help                           Show this help message and exit\n");
}

namespace {
/**
 * @brief Consumes and returns the value following the current command-line option.
 * @param i In/out index of the current argument; advanced to the consumed value.
 * @param argc Total number of arguments.
 * @param argv Argument vector.
 * @return The next argument as an owning string.
 * @throws std::runtime_error If the option is the final argument and has no value.
 */
std::string read_value(int& i, int argc, char** argv) {
    if (i + 1 >= argc) throw std::runtime_error("Missing option value.");
    ++i;
    return std::string(argv[i]);
}

/// Converts an integer option using the command-line parser's permissive C semantics.
int parse_int(const std::string& s) { return std::atoi(s.c_str()); }
/// Converts a size option using the command-line parser's unsigned C semantics.
size_t parse_size(const std::string& s) { return static_cast<size_t>(std::strtoull(s.c_str(), nullptr, 10)); }
/// Converts a floating-point option using the command-line parser's permissive C semantics.
double parse_double(const std::string& s) { return std::atof(s.c_str()); }
}

Args parse_args(int argc, char** argv) {
    Args args;

    if (argc == 1) {
        print_usage(argv[0]);
        throw std::runtime_error("No arguments were provided.");
    }

    for (int i = 1; i < argc; ++i) {
        std::string opt(argv[i]);
        if (opt == "-h" || opt == "--help") { print_usage(argv[0]); std::exit(0); }
        if (opt == "-r" || opt == "--ref") { args.ref = read_value(i, argc, argv); continue; }
        if (opt == "-q" || opt == "--query") { args.query = read_value(i, argc, argv); continue; }
        if (opt == "-o" || opt == "--out") { args.out = read_value(i, argc, argv); continue; }
        if (opt == "--max-instances") { args.max_instances = parse_int(read_value(i, argc, argv)); continue; }
        if (opt == "--threads-per-mafft") { args.threads_per_mafft = parse_int(read_value(i, argc, argv)); continue; }
        if (opt == "--chunk-size") { args.chunk_size = parse_int(read_value(i, argc, argv)); continue; }
        if (opt == "--target-chunk-bases") { args.target_chunk_bases = parse_size(read_value(i, argc, argv)); continue; }
        if (opt == "--long-chunk-size") { (void)read_value(i, argc, argv); continue; } // deprecated; unified dispatch ignores this option
        if (opt == "--long-target-chunk-bases") { (void)read_value(i, argc, argv); continue; } // deprecated; unified dispatch ignores this option
        if (opt == "--single-long-ratio") { (void)read_value(i, argc, argv); continue; } // deprecated; unified dispatch ignores this option
        if (opt == "--mafft-timeout-sec") { args.mafft_timeout_sec = parse_int(read_value(i, argc, argv)); continue; }
        if (opt == "--slow-chunk-sec") { args.slow_chunk_sec = parse_double(read_value(i, argc, argv)); continue; }
        if (opt == "--temp-dir") { args.temp_dir = read_value(i, argc, argv); args.temp_dir_provided = true; continue; }
        if (opt == "--cleanup") { args.cleanup = true; continue; }
        if (opt == "--len-tolerance") { (void)read_value(i, argc, argv); continue; } // deprecated; no near/short/long bucketing
        if (opt == "--n-threshold") { (void)read_value(i, argc, argv); continue; } // deprecated; no highN/lowN bucketing
        if (opt == "--save-failed-stderr") { args.save_failed_stderr = true; continue; }
        if (opt == "--disable-auto-chunk") { args.disable_auto_chunk = true; continue; }
        if (opt == "--max-pending-factor") { args.max_pending_factor = parse_int(read_value(i, argc, argv)); continue; }
        if (opt == "--enable-qc") { args.enable_qc = true; continue; }
        if (opt == "--warn-lost-ratio") { args.warn_lost_ratio = parse_double(read_value(i, argc, argv)); continue; }
        if (opt == "--fail-lost-ratio") { args.fail_lost_ratio = parse_double(read_value(i, argc, argv)); continue; }
        if (opt == "--warn-max-block-len") { args.warn_max_block_len = parse_int(read_value(i, argc, argv)); continue; }
        if (opt == "--fail-max-block-len") { args.fail_max_block_len = parse_int(read_value(i, argc, argv)); continue; }
        if (opt == "--duplicate-id-policy") { args.duplicate_id_policy = read_value(i, argc, argv); continue; }
        if (opt == "--allowed-seq-chars") { args.allowed_seq_chars = read_value(i, argc, argv); continue; }
        if (opt == "--review-n-ratio") { args.review_n_ratio = parse_double(read_value(i, argc, argv)); continue; }
        if (opt == "--exclude-n-ratio") { args.exclude_n_ratio = parse_double(read_value(i, argc, argv)); continue; }
        if (opt == "--exclude-min-len-ratio") { args.exclude_min_len_ratio = parse_double(read_value(i, argc, argv)); continue; }
        if (opt == "--exclude-max-len-ratio") { args.exclude_max_len_ratio = parse_double(read_value(i, argc, argv)); continue; }
        if (opt == "--reasonable-min-len-ratio") { args.reasonable_min_len_ratio = parse_double(read_value(i, argc, argv)); continue; }
        if (opt == "--reasonable-max-len-ratio") { args.reasonable_max_len_ratio = parse_double(read_value(i, argc, argv)); continue; }
        if (opt == "--non-target-keywords") { args.non_target_keywords = read_value(i, argc, argv); continue; }
        if (opt == "--engineered-keywords") { args.engineered_keywords = read_value(i, argc, argv); continue; }
        if (opt == "--align-review-samples") { args.align_review_samples = true; continue; }
        throw std::runtime_error("Unknown option: " + opt + ". Use --help to see available options.");
    }

    if (args.ref.empty() || args.query.empty() || args.out.empty()) {
        print_usage(argv[0]);
        throw std::runtime_error("Missing required arguments.");
    }
    if (args.max_instances < 1) args.max_instances = 1;
    if (args.threads_per_mafft < 1) args.threads_per_mafft = 1;
    if (args.chunk_size < 1) args.chunk_size = 1;
    if (args.max_pending_factor < 1) args.max_pending_factor = 1;
    if (args.mafft_timeout_sec < 0) args.mafft_timeout_sec = 0;
    if (args.slow_chunk_sec < 0.0) args.slow_chunk_sec = 0.0;
    if (args.duplicate_id_policy != "exclude" && args.duplicate_id_policy != "rename") {
        throw std::runtime_error("--duplicate-id-policy must be exclude or rename.");
    }
    if (args.review_n_ratio < 0.0 || args.exclude_n_ratio < 0.0 || args.review_n_ratio > 1.0 || args.exclude_n_ratio > 1.0) {
        throw std::runtime_error("N-ratio thresholds must be within [0, 1].");
    }
    if (args.review_n_ratio > args.exclude_n_ratio) {
        throw std::runtime_error("--review-n-ratio cannot exceed --exclude-n-ratio.");
    }
    if (args.exclude_min_len_ratio < 0.0 || args.exclude_max_len_ratio <= 0.0 || args.exclude_min_len_ratio > args.exclude_max_len_ratio) {
        throw std::runtime_error("Invalid exclude length-ratio thresholds.");
    }
    if (args.reasonable_min_len_ratio < 0.0 || args.reasonable_max_len_ratio <= 0.0 || args.reasonable_min_len_ratio > args.reasonable_max_len_ratio) {
        throw std::runtime_error("Invalid reasonable length-ratio thresholds.");
    }
    return args;
}
