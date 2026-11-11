/**
 * @file fasta_io.hpp
 * @brief Streaming FASTA input and alignment-output helpers.
 */
#pragma once
#include "common.hpp"

#include <fstream>
#include <iosfwd>
#include <string>
#include <vector>

/**
 * @brief Extracts the primary identifier from a FASTA header.
 * @param header Complete header text without the leading `>`.
 * @return The first whitespace-delimited token, or @p header when no token can
 *         be extracted.
 * @note This is the canonical identifier rule used by FASTA input, MAFFT
 *       output validation, and mapout section validation.
 */
std::string fasta_primary_id(const std::string& header);

/**
 * @brief Reads the next FASTA record from a stream.
 *
 * Sequence lines are concatenated, and a trailing carriage return from CRLF
 * input is removed. @p pending_header carries a previously encountered
 * next-record header between calls.
 *
 * @param in Open input stream positioned at or before the next record.
 * @param rec Destination record; reset before parsing.
 * @param pending_header Header text without `>` retained between calls.
 * @param source_path Optional source path included in validation and I/O errors.
 * @return `true` when a record was parsed; `false` when no further non-empty
 *         header is available.
 * @throws std::invalid_argument If non-FASTA content appears before the first
 *         header, a header has no identifier, or an embedded NUL byte is found.
 * @throws std::runtime_error If the stream reports an I/O failure.
 */
bool read_next_fasta_record(
    std::ifstream& in,
    FastaRecord& rec,
    std::string& pending_header,
    const std::string& source_path = ""
);

/**
 * @brief Reads a FASTA file that must contain exactly one record.
 * @param path Input FASTA path.
 * @return The single parsed record.
 * @throws std::invalid_argument If the path is not a regular uncompressed FASTA
 *         file, the content is malformed, the reference sequence is empty, or
 *         the record count is not exactly one.
 * @throws std::runtime_error If the file cannot be opened or read.
 */
FastaRecord read_single_fasta(const std::string& path);

/**
 * @brief Counts FASTA headers without loading sequence bodies into memory.
 * @param path Input FASTA path.
 * @return Number of validated FASTA records. A truly empty or whitespace-only
 *         query file returns zero.
 * @throws std::invalid_argument If the path is not a regular uncompressed FASTA
 *         file or the content is malformed.
 * @throws std::runtime_error If the file cannot be opened or read.
 * @note Leading whitespace before `>` is rejected as malformed content.
 */
size_t count_fasta_sequences(const std::string& path);

/**
 * @brief Writes one FASTA record with sequence lines wrapped at 80 characters.
 * @param out Destination stream.
 * @param rec Record whose stored header and sequence are written verbatim.
 */
void write_fasta_record(std::ostream& out, const FastaRecord& rec);

/**
 * @brief Creates a query-chunk FASTA file from an ordered record collection.
 * @param records Records to serialize in their current order; an empty vector
 *        produces an empty file.
 * @param path Destination path, truncated if it already exists.
 * @throws std::runtime_error If the destination cannot be opened or fully written.
 */
void write_chunk_file(const std::vector<FastaRecord>& records, const std::string& path);

/**
 * @brief Appends MAFFT FASTA output while omitting its first record.
 *
 * Each MAFFT chunk includes the shared reference as its first record. This
 * helper drops that record so the caller can keep a single reference in the
 * combined MSA.
 *
 * @param mafft_output Complete FASTA text captured from MAFFT stdout.
 * @param out Destination stream for records after the first header.
 * @note If the text contains fewer than two headers, nothing is appended.
 */
void append_mafft_stdout_without_ref(const std::string& mafft_output, std::ostream& out);
