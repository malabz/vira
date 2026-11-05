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
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string_view>

bool read_next_fasta_record(std::ifstream& in, FastaRecord& rec, std::string& pending_header) {
    rec = FastaRecord{};
    std::string line;

    // `pending_header` carries the one-record look-ahead between successive calls.
    if (pending_header.empty()) {
        while (std::getline(in, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (!line.empty() && line[0] == '>') {
                pending_header = line.substr(1);
                break;
            }
        }
    }

    if (pending_header.empty()) return false;

    rec.header = pending_header;
    std::istringstream iss(rec.header);
    iss >> rec.id;
    if (rec.id.empty()) rec.id = rec.header;
    pending_header.clear();

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line[0] == '>') {
            pending_header = line.substr(1);
            break;
        }
        rec.seq += line;
    }

    return true;
}

FastaRecord read_single_fasta(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open reference FASTA: " + path);

    FastaRecord rec;
    std::string pending;
    if (!read_next_fasta_record(in, rec, pending)) throw std::runtime_error("Reference FASTA is empty: " + path);

    FastaRecord extra;
    if (read_next_fasta_record(in, extra, pending)) {
        throw std::runtime_error("Reference FASTA must contain exactly one sequence: " + path);
    }
    return rec;
}

size_t count_fasta_sequences(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open query FASTA: " + path);

    size_t count = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') ++count;
    }
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
