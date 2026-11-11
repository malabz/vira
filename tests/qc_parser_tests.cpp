/**
 * @file qc_parser_tests.cpp
 * @brief Standalone negative and boundary tests for MAFFT mapout parsing.
 *
 * This executable intentionally uses no third-party test framework. It links
 * directly against `src/qc.cpp`, creates isolated temporary fixtures, and
 * returns a non-zero status when any expectation fails.
 */

#include "qc.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <system_error>
#include <unistd.h>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {
struct TempDirectoryCleanup {
    fs::path path;

    ~TempDirectoryCleanup() {
        std::error_code cleanup_error;
        fs::remove_all(path, cleanup_error);
    }
};

class TestSuite {
public:
    explicit TestSuite(fs::path fixture_dir) : fixture_dir_(std::move(fixture_dir)) {}

    int run() {
        test_valid_mapout();
        test_mapout_identity_validation();
        test_missing_stats_are_omitted();

        expect_format_error(
            "invalid-original-integer.map",
            ">sample\nA,not-an-integer,1\n",
            {"line 2", "expected an integer for original position", "not-an-integer"}
        );
        expect_format_error(
            "invalid-reference-integer.map",
            ">sample\nA,1,2x\n",
            {"line 2", "trailing characters in reference position", "2x"}
        );
        expect_format_error(
            "out-of-range-integer.map",
            ">sample\nA,999999999999999999999999999999,1\n",
            {"line 2", "value out of range for original position"}
        );
        expect_format_error(
            "missing-field.map",
            ">sample\nA,1\n",
            {"line 2", "expected exactly three comma-separated fields"}
        );
        expect_format_error(
            "extra-field.map",
            ">sample\nA,1,1,unexpected\n",
            {"line 2", "reference position", "1,unexpected"}
        );
        expect_format_error(
            "empty-symbol.map",
            ">sample\n,1,1\n",
            {"line 2", "sequence-symbol field is empty"}
        );
        expect_format_error(
            "data-before-section.map",
            "A,1,1\n>sample\nC,2,2\n",
            {"line 1", "data row appears before the first sequence section"}
        );
        expect_format_error(
            "empty-section.map",
            "> \t\nA,1,1\n",
            {"line 1", "sequence section identifier is empty"}
        );
        expect_format_error(
            "duplicate-section.map",
            ">sample\nA,1,1\n>sample\nC,2,2\n",
            {"line 3", "duplicate sequence section identifier", "sample"}
        );

        std::string nul_content = ">sample\nA,1,1";
        nul_content.push_back('\0');
        nul_content += "corrupted\n";
        expect_format_error(
            "embedded-nul.map",
            nul_content,
            {"line 2", "embedded NUL byte", "binary or corrupted"}
        );

        expect_format_error(
            "zero-original-coordinate.map",
            ">sample\nA,0,1\n",
            {"line 2", "original position must be between 1 and", "got '0'"}
        );
        expect_format_error(
            "negative-reference-coordinate.map",
            ">sample\nA,1,-2\n",
            {"line 2", "reference position must be between 1 and", "got '-2'"}
        );
        expect_format_error(
            "max-original-coordinate.map",
            ">sample\nA," + std::to_string(std::numeric_limits<long long>::max()) + ",-\n",
            {"line 2", "original position must be between 1 and"}
        );
        expect_format_error(
            "no-sections.map",
            "# comment only\n   \n\t# indented comment\n",
            {"no sequence sections were found"}
        );

        test_missing_file();

        if (failures_ == 0) {
            std::cout << "PASS: " << checks_ << " mapout parser checks completed.\n";
            return 0;
        }
        std::cerr << "FAIL: " << failures_ << " of " << checks_ << " mapout parser checks failed.\n";
        return 1;
    }

private:
    fs::path fixture_dir_;
    int checks_ = 0;
    int failures_ = 0;

    void check(bool condition, const std::string& message) {
        ++checks_;
        if (condition) return;
        ++failures_;
        std::cerr << "[FAIL] " << message << '\n';
    }

    bool write_fixture(const fs::path& path, const std::string& content) {
        std::ofstream out(path, std::ios::binary | std::ios::trunc);
        if (!out) {
            check(false, "Cannot create test fixture: " + path.string());
            return false;
        }
        out.write(content.data(), static_cast<std::streamsize>(content.size()));
        out.close();
        if (!out) {
            check(false, "Cannot fully write test fixture: " + path.string());
            return false;
        }
        return true;
    }

    void expect_message_parts(
        const std::string& test_name,
        const std::string& message,
        const std::vector<std::string>& expected_parts) {

        for (const auto& part : expected_parts) {
            check(
                message.find(part) != std::string::npos,
                test_name + ": expected diagnostic to contain '" + part + "', got: " + message
            );
        }
    }

    void expect_format_error(
        const std::string& filename,
        const std::string& content,
        const std::vector<std::string>& expected_parts,
        std::vector<std::string> expected_ids = {"sample"}) {

        const fs::path path = fixture_dir_ / filename;
        if (!write_fixture(path, content)) return;

        bool caught_expected = false;
        try {
            (void)parse_map_to_deleted_blocks(path.string(), expected_ids);
        } catch (const MapoutFormatError& error) {
            caught_expected = true;
            const std::string message = error.what();
            check(message.find(path.string()) != std::string::npos, filename + ": diagnostic must include the path");
            expect_message_parts(filename, message, expected_parts);
        } catch (const MapoutIoError& error) {
            check(false, filename + ": expected MapoutFormatError, got MapoutIoError: " + error.what());
        }
        check(caught_expected, filename + ": expected MapoutFormatError");
    }

    void test_valid_mapout() {
        const fs::path path = fixture_dir_ / "valid.map";
        const std::string content =
            "# leading comment\n"
            "   \n"
            ">sample1\n"
            "A,1,1\n"
            "C,2,-\n"
            "G,3,-\n"
            "T,4,4\n"
            "\t# indented comment\n"
            ">sample2\n"
            "A,1,-\n"
            "C,2,2\n";
        if (!write_fixture(path, content)) return;

        auto parsed = parse_map_to_deleted_blocks(path.string(), {"sample1", "sample2"});
        const auto& blocks = parsed.first;
        const auto& stats = parsed.second;

        check(blocks.size() == 2, "valid mapout should produce two deleted blocks");
        check(stats.size() == 2, "valid mapout should produce statistics for two samples");

        const auto sample1 = stats.find("sample1");
        check(sample1 != stats.end(), "valid mapout should include sample1 statistics");
        if (sample1 != stats.end()) {
            check(sample1->second.lost_count == 2, "sample1 lost_count should be 2");
            check(sample1->second.lost_block_count == 1, "sample1 should have one lost block");
            check(sample1->second.max_lost_block_len == 2, "sample1 maximum block length should be 2");
        }

        const auto sample2 = stats.find("sample2");
        check(sample2 != stats.end(), "valid mapout should include sample2 statistics");
        if (sample2 != stats.end()) {
            check(sample2->second.lost_count == 1, "sample2 lost_count should be 1");
            check(sample2->second.lost_block_count == 1, "sample2 should have one lost block");
            check(sample2->second.max_lost_block_len == 1, "sample2 maximum block length should be 1");
        }

        if (blocks.size() == 2) {
            check(blocks[0].seq_id == "sample1", "first block should belong to sample1");
            check(blocks[0].block_start_orig == 2, "sample1 block should start at original position 2");
            check(blocks[0].block_end_orig == 3, "sample1 block should end at original position 3");
            check(blocks[0].block_len == 2, "sample1 block length should be 2");
            check(blocks[0].left_ref_pos == "1", "sample1 block left reference position should be 1");
            check(blocks[0].right_ref_pos == "4", "sample1 block right reference position should be 4");
            check(blocks[1].seq_id == "sample2", "second block should belong to sample2");
            check(blocks[1].left_ref_pos.empty(), "edge block should have an empty left reference position");
            check(blocks[1].right_ref_pos == "2", "sample2 block right reference position should be 2");
        }
    }

    void test_mapout_identity_validation() {
        const fs::path descriptive_path = fixture_dir_ / "descriptive-sections.map";
        const std::string descriptive_content =
            ">q2 second query description\n"
            "C,1,1\n"
            ">q1 first query description\n"
            "A,1,1\n";
        if (write_fixture(descriptive_path, descriptive_content)) {
            auto parsed = parse_map_to_deleted_blocks(descriptive_path.string(), {"q1", "q2"});
            check(parsed.second.size() == 2, "descriptive mapout sections should preserve both primary IDs");
            check(parsed.second.find("q1") != parsed.second.end(), "descriptive q1 section should use primary ID q1");
            check(parsed.second.find("q2") != parsed.second.end(), "reordered descriptive q2 section should be accepted");
        }

        expect_format_error(
            "missing-expected-section.map",
            ">q1\nA,1,1\n",
            {"missing=[q2]", "unexpected=[]"},
            {"q1", "q2"}
        );
        expect_format_error(
            "replaced-expected-section.map",
            ">q1\nA,1,1\n>rogue\nC,1,1\n",
            {"missing=[q2]", "unexpected=[rogue]"},
            {"q1", "q2"}
        );
        expect_format_error(
            "case-changed-section.map",
            ">Q1\nA,1,1\n>q2\nC,1,1\n",
            {"missing=[q1]", "unexpected=[Q1]"},
            {"q1", "q2"}
        );
        expect_format_error(
            "extra-section.map",
            ">q1\nA,1,1\n>q2\nC,1,1\n>rogue\nG,1,1\n",
            {"missing=[]", "unexpected=[rogue]"},
            {"q1", "q2"}
        );
        expect_format_error(
            "sorted-mismatch-diagnostic.map",
            ">zeta\nA,1,1\n>alpha\nC,1,1\n",
            {"missing=[q1,q2]", "unexpected=[alpha,zeta]"},
            {"q2", "q1"}
        );
        expect_format_error(
            "duplicate-primary-section.map",
            ">q1 first description\nA,1,1\n>q1 second description\nC,1,1\n",
            {"line 3", "duplicate sequence section identifier", "q1"},
            {"q1"}
        );
        expect_format_error(
            "duplicate-expected-ids.map",
            ">q1\nA,1,1\n",
            {"Invalid expected MAFFT mapout identifier set", "duplicates=[q1]"},
            {"q1", "q1"}
        );
    }

    void test_missing_stats_are_omitted() {
        std::map<std::string, SampleMeta> metadata;
        metadata["q1"] = SampleMeta{"all", 10, 0.0, 0.0};
        metadata["q2"] = SampleMeta{"all", 10, 0.0, 0.0};

        std::unordered_map<std::string, QcStats> statistics;
        statistics["q1"] = QcStats{};

        const Args args;
        const auto rows = build_sample_qc_rows(metadata, statistics, args);
        check(rows.size() == 1, "missing statistics must not create a zero-loss PASS SampleRow");
        if (rows.size() == 1) {
            check(rows[0].seq_id == "q1", "only the sample with present statistics should produce a SampleRow");
            check(rows[0].qc_level == "PASS", "present zero-loss statistics should retain PASS semantics");
        }
    }

    void test_missing_file() {
        const fs::path path = fixture_dir_ / "missing.map";
        bool caught_expected = false;
        try {
            (void)parse_map_to_deleted_blocks(path.string(), {"sample"});
        } catch (const MapoutIoError& error) {
            caught_expected = true;
            const std::string message = error.what();
            check(message.find("Cannot open MAFFT mapout file") != std::string::npos, "missing file should report an open failure");
            check(message.find(path.string()) != std::string::npos, "missing-file diagnostic should include the path");
        } catch (const MapoutFormatError& error) {
            check(false, "missing file should raise MapoutIoError, got MapoutFormatError: " + std::string(error.what()));
        }
        check(caught_expected, "missing file should raise MapoutIoError");
    }
};

bool create_unique_temp_directory(fs::path& output) {
    std::error_code temp_error;
    const fs::path base = fs::temp_directory_path(temp_error);
    if (temp_error) {
        std::cerr << "Cannot locate the system temporary directory: " << temp_error.message() << '\n';
        return false;
    }

    const auto timestamp = std::chrono::steady_clock::now().time_since_epoch().count();
    const std::string stem = "vira-qc-parser-tests-" + std::to_string(static_cast<long long>(getpid())) +
        "-" + std::to_string(timestamp);
    for (int attempt = 0; attempt < 100; ++attempt) {
        output = base / (stem + "-" + std::to_string(attempt));
        std::error_code create_error;
        if (fs::create_directory(output, create_error)) return true;
        if (create_error && create_error != std::errc::file_exists) {
            std::cerr << "Cannot create test temporary directory '" << output.string()
                      << "': " << create_error.message() << '\n';
            return false;
        }
    }
    std::cerr << "Cannot create a unique temporary directory for mapout parser tests.\n";
    return false;
}
}

int main() {
    fs::path fixture_dir;
    if (!create_unique_temp_directory(fixture_dir)) return 2;
    TempDirectoryCleanup cleanup{fixture_dir};
    TestSuite suite(fixture_dir);
    return suite.run();
}
