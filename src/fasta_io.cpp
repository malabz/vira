/**
 * @file fasta_io.cpp
 * @brief Implements streaming FASTA parsing, normalized FASTA output, and chunk-result merging.
 *
 * Parsing is intentionally streaming: only the current record and one look-ahead
 * header are retained. Sequence lines are preserved apart from trailing CR removal
 * and line joining, and emitted records are wrapped at 80 characters.
 */

#include "fasta_io.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sys/stat.h>

namespace {
bool is_blank(const std::string& text) {
    return std::all_of(text.begin(), text.end(), [](unsigned char c) { return std::isspace(c) != 0; });
}

bool has_header_identifier(const std::string& header) {
    return !header.empty() && !is_blank(header);
}

std::string lowercase_extension(const std::string& path) {
    const size_t slash = path.find_last_of("/\\");
    const size_t dot = path.find_last_of('.');
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) return "";

    std::string extension = path.substr(dot);
    std::transform(extension.begin(), extension.end(), extension.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return extension;
}

std::string input_description(const char* role) {
    return std::string(role) + " input";
}

void validate_fasta_path(const std::string& path, const char* role) {
    const std::string extension = lowercase_extension(path);
    const std::array<const char*, 9> archive_extensions = {
        ".gz", ".zip", ".bz2", ".xz", ".zst", ".7z", ".rar", ".tgz", ".tar"
    };

    if (std::find(archive_extensions.begin(), archive_extensions.end(), extension) != archive_extensions.end()) {
        throw std::invalid_argument(
            "Expected an uncompressed FASTA file for " + input_description(role) +
            ", got compressed or archive extension '" + extension + "': " + path +
            ". Please decompress it first."
        );
    }
    if (extension == ".fastq" || extension == ".fq") {
        throw std::invalid_argument(
            "Expected a FASTA file for " + input_description(role) +
            ", got FASTQ extension '" + extension + "': " + path + "."
        );
    }
    if (extension != ".fa" && extension != ".fasta" && extension != ".fas" && extension != ".fna") {
        const std::string actual = extension.empty() ? "no extension" : "extension '" + extension + "'";
        throw std::invalid_argument(
            "Expected a FASTA file for " + input_description(role) +
            " with extension .fa, .fasta, .fas, or .fna; got " + actual + ": " + path + "."
        );
    }

    struct stat info {};
    if (stat(path.c_str(), &info) != 0) {
        throw std::runtime_error(
            "Cannot access " + input_description(role) + " FASTA '" + path + "': " + std::strerror(errno) + "."
        );
    }
    if (!S_ISREG(info.st_mode)) {
        throw std::invalid_argument(
            "Expected a regular FASTA file for " + input_description(role) +
            ", got a non-regular path: " + path + "."
        );
    }
}

std::string compressed_signature_name(const std::array<unsigned char, 8>& bytes, std::streamsize count) {
    const auto has = [&](std::initializer_list<unsigned char> signature) {
        if (count < static_cast<std::streamsize>(signature.size())) return false;
        return std::equal(signature.begin(), signature.end(), bytes.begin());
    };

    if (has({0x1f, 0x8b})) return "gzip-compressed";
    if (has({0x50, 0x4b, 0x03, 0x04}) || has({0x50, 0x4b, 0x05, 0x06}) || has({0x50, 0x4b, 0x07, 0x08})) {
        return "ZIP archive";
    }
    if (has({0x42, 0x5a, 0x68})) return "bzip2-compressed";
    if (has({0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00})) return "XZ-compressed";
    if (has({0x28, 0xb5, 0x2f, 0xfd})) return "Zstandard-compressed";
    if (has({0x37, 0x7a, 0xbc, 0xaf, 0x27, 0x1c})) return "7-Zip archive";
    if (has({0x52, 0x61, 0x72, 0x21, 0x1a, 0x07})) return "RAR archive";
    return "";
}

void reject_compressed_content(std::ifstream& in, const std::string& path, const char* role) {
    std::array<unsigned char, 8> bytes {};
    in.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    const std::streamsize count = in.gcount();
    if (in.bad()) {
        throw std::runtime_error("I/O error while inspecting " + input_description(role) + " FASTA: " + path + ".");
    }

    in.clear();
    in.seekg(0, std::ios::beg);
    if (!in) {
        throw std::runtime_error("Cannot seek " + input_description(role) + " FASTA after inspection: " + path + ".");
    }

    const std::string format = compressed_signature_name(bytes, count);
    if (!format.empty()) {
        throw std::invalid_argument(
            "Expected an uncompressed FASTA file for " + input_description(role) +
            ", detected " + format + " content in '" + path + "'. Please decompress it first."
        );
    }
}

std::string source_label(const std::string& source_path) {
    return source_path.empty() ? "the input stream" : "'" + source_path + "'";
}

[[noreturn]] void throw_invalid_fasta(const std::string& source_path, const std::string& detail) {
    throw std::invalid_argument("Invalid FASTA format in " + source_label(source_path) + ": " + detail);
}

void reject_nul(const std::string& line, const std::string& source_path) {
    if (line.find('\0') != std::string::npos) {
        throw_invalid_fasta(source_path, "embedded NUL byte detected; the file may be binary or corrupted.");
    }
}
}

std::string fasta_primary_id(const std::string& header) {
    std::istringstream input(header);
    std::string id;
    input >> id;
    return id.empty() ? header : id;
}

bool read_next_fasta_record(
    std::ifstream& in,
    FastaRecord& rec,
    std::string& pending_header,
    const std::string& source_path) {

    rec = FastaRecord{};
    std::string line;

    // `pending_header` carries the one-record look-ahead between successive calls.
    if (pending_header.empty()) {
        while (std::getline(in, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            reject_nul(line, source_path);
            if (is_blank(line)) continue;
            if (!line.empty() && line[0] == '>') {
                pending_header = line.substr(1);
                if (!has_header_identifier(pending_header)) {
                    throw_invalid_fasta(source_path, "header is missing a sequence identifier.");
                }
                break;
            }
            throw_invalid_fasta(source_path, "expected a header line starting with '>'.");
        }
        if (in.bad()) throw std::runtime_error("I/O error while reading FASTA input " + source_label(source_path) + ".");
    }

    if (pending_header.empty()) return false;
    if (!has_header_identifier(pending_header)) {
        throw_invalid_fasta(source_path, "header is missing a sequence identifier.");
    }

    rec.header = pending_header;
    rec.id = fasta_primary_id(rec.header);
    pending_header.clear();

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        reject_nul(line, source_path);
        if (!line.empty() && line[0] == '>') {
            pending_header = line.substr(1);
            if (!has_header_identifier(pending_header)) {
                throw_invalid_fasta(source_path, "header is missing a sequence identifier.");
            }
            break;
        }
        rec.seq += line;
    }
    if (in.bad()) throw std::runtime_error("I/O error while reading FASTA input " + source_label(source_path) + ".");

    return true;
}

FastaRecord read_single_fasta(const std::string& path) {
    validate_fasta_path(path, "reference");
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error(
            "Cannot open reference input FASTA for reading: " + path + ". Check the file permissions."
        );
    }
    reject_compressed_content(in, path, "reference");

    FastaRecord rec;
    std::string pending;
    if (!read_next_fasta_record(in, rec, pending, path)) {
        throw std::invalid_argument("Reference FASTA is empty or contains no sequence record: " + path + ".");
    }
    if (rec.seq.empty()) {
        throw std::invalid_argument("Reference FASTA sequence is empty: " + path + ".");
    }

    FastaRecord extra;
    if (read_next_fasta_record(in, extra, pending, path)) {
        throw std::invalid_argument("Reference FASTA must contain exactly one sequence record: " + path + ".");
    }
    return rec;
}

size_t count_fasta_sequences(const std::string& path) {
    validate_fasta_path(path, "query");
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error(
            "Cannot open query input FASTA for reading: " + path + ". Check the file permissions."
        );
    }
    reject_compressed_content(in, path, "query");

    size_t count = 0;
    size_t line_number = 0;
    bool seen_header = false;
    std::string line;
    while (std::getline(in, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.find('\0') != std::string::npos) {
            throw std::invalid_argument(
                "Invalid FASTA format in query file '" + path + "' at line " + std::to_string(line_number) +
                ": embedded NUL byte detected; the file may be binary or corrupted."
            );
        }
        if (is_blank(line)) continue;
        if (line[0] == '>') {
            const std::string header = line.substr(1);
            if (!has_header_identifier(header)) {
                throw std::invalid_argument(
                    "Invalid FASTA format in query file '" + path + "' at line " + std::to_string(line_number) +
                    ": header is missing a sequence identifier."
                );
            }
            ++count;
            seen_header = true;
            continue;
        }
        if (!seen_header) {
            throw std::invalid_argument(
                "Invalid FASTA format in query file '" + path + "' at line " + std::to_string(line_number) +
                ": expected a header line starting with '>'."
            );
        }
    }
    if (in.bad()) throw std::runtime_error("I/O error while reading query FASTA: " + path + ".");
    return count;
}

void write_fasta_record(std::ostream& out, const FastaRecord& rec) {
    out << '>' << rec.header << '\n';
    constexpr size_t width = 80;
    for (size_t i = 0; i < rec.seq.size(); i += width) {
        out << rec.seq.substr(i, std::min(width, rec.seq.size() - i)) << '\n';
    }
}

void write_chunk_file(const std::vector<FastaRecord>& records, const std::string& path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write temporary chunk file: " + path);
    for (const auto& r : records) write_fasta_record(out, r);
    out.close();
    if (!out) throw std::runtime_error("Failed while writing temporary chunk FASTA: " + path + ".");
}

void append_mafft_stdout_without_ref(const std::string& mafft_output, std::ostream& out) {
    bool first_header_seen = false;
    bool started_writing = false;
    size_t pos = 0;

    // MAFFT emits the reference as the first FASTA record for every chunk. Suppress
    // that record so the caller can keep a single reference at the start of the MSA.
    while (pos < mafft_output.size()) {
        size_t next = mafft_output.find('\n', pos);
        if (next == std::string::npos) next = mafft_output.size();
        std::string_view line(mafft_output.data() + pos, next - pos);
        bool has_newline = next < mafft_output.size();

        if (!line.empty() && line[0] == '>') {
            if (!first_header_seen) {
                first_header_seen = true;
            } else {
                started_writing = true;
                out.write(line.data(), static_cast<std::streamsize>(line.size()));
                if (has_newline) out << '\n';
            }
        } else if (started_writing) {
            out.write(line.data(), static_cast<std::streamsize>(line.size()));
            if (has_newline) out << '\n';
        }
        pos = has_newline ? next + 1 : next;
    }
}
