#pragma once
#include "common.hpp"

#include <fstream>
#include <iosfwd>
#include <string>
#include <vector>

bool read_next_fasta_record(std::ifstream& in, FastaRecord& rec, std::string& pending_header);
FastaRecord read_single_fasta(const std::string& path);
size_t count_fasta_sequences(const std::string& path);
void write_fasta_record(std::ostream& out, const FastaRecord& rec);
void write_chunk_file(const std::vector<FastaRecord>& records, const std::string& path);
void append_mafft_stdout_without_ref(const std::string& mafft_output, std::ostream& out);
