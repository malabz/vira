"""End-to-end tests for MAFFT process boundaries and output protection."""

from __future__ import annotations

import csv
import os
from pathlib import Path
import stat
import subprocess
import tempfile
import textwrap
import unittest


VIRA_BINARY = os.environ.get("VIRA_TEST_BINARY")


MOCK_MAFFT = r"""#!/usr/bin/env python3
from pathlib import Path
import os
import sys
import time


def read_fasta(path):
    records = []
    header = None
    sequence = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if raw_line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence)))
            header = raw_line[1:]
            sequence = []
        else:
            sequence.append(raw_line.strip())
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def primary_id(header):
    return header.split()[0]


def write_fasta(records):
    for header, sequence in records:
        sys.stdout.write(f">{header}\n{sequence}\n")

arguments = sys.argv[1:]
mode = os.environ.get("VIRA_MOCK_MAFFT_MODE", "normal")
if mode == "sleep":
    time.sleep(4)
if mode == "nonzero":
    sys.stderr.write("mock MAFFT rejected the input\n")
    raise SystemExit(9)

try:
    query_path = Path(arguments[arguments.index("--addfragments") + 1])
    reference_path = Path(arguments[-1])
except (ValueError, IndexError):
    sys.stderr.write("mock MAFFT received an invalid command\n")
    raise SystemExit(64)

if mode == "malformed-output":
    sys.stdout.write("not FASTA output\n")
    raise SystemExit(0)

reference_records = read_fasta(reference_path)
query_records = read_fasta(query_path)

if "--mapout" in arguments:
    map_path = Path(str(query_path) + ".map")
    if mode == "malformed-mapout":
        map_path.write_text(">query1\nA,not-an-integer,1\n", encoding="utf-8")
    else:
        map_records = [(primary_id(header), sequence) for header, sequence in query_records]
        if mode == "mapout-missing-section":
            map_records = map_records[:1]
        elif mode == "mapout-replaced-section":
            map_records = map_records[:1] + [("rogue", map_records[-1][1])]
        elif mode == "mapout-extra-section":
            map_records.append(("rogue", map_records[-1][1]))
        elif mode == "mapout-reverse-sections":
            map_records.reverse()
        elif mode == "mapout-duplicate-section":
            map_records.append(map_records[0])
        with map_path.open("w", encoding="utf-8") as stream:
            for record_id, record_sequence in map_records:
                stream.write(f">{record_id}\n")
                for position, symbol in enumerate(record_sequence, start=1):
                    stream.write(f"{symbol},{position},{position}\n")

if mode in {
    "normal",
    "malformed-mapout",
    "mapout-missing-section",
    "mapout-replaced-section",
    "mapout-extra-section",
    "mapout-reverse-sections",
    "mapout-duplicate-section",
}:
    sys.stdout.buffer.write(reference_path.read_bytes())
    sys.stdout.buffer.write(query_path.read_bytes())
else:
    output_records = list(reference_records) + list(query_records)
    if mode == "reference-replaced":
        output_records[0] = ("rogue_reference", output_records[0][1])
    elif mode == "query-replaced":
        output_records[1] = ("renamed_q1", output_records[1][1])
    elif mode == "case-changed-query":
        output_records[1] = ("Q1", output_records[1][1])
    elif mode == "unexpected-query":
        output_records[-1] = ("rogue", output_records[-1][1])
    elif mode == "multiple-unexpected-queries":
        output_records[1] = ("zeta", output_records[1][1])
        output_records[2] = ("alpha", output_records[2][1])
    elif mode == "duplicate-query":
        output_records[2] = (primary_id(output_records[1][0]), output_records[2][1])
    elif mode == "reference-not-first":
        output_records = [output_records[1], output_records[0], *output_records[2:]]
    elif mode == "reverse-query-order":
        output_records = [output_records[0], *reversed(output_records[1:])]
    elif mode == "descriptive-headers":
        output_records = [
            (primary_id(header) + " mock description", sequence)
            for header, sequence in output_records
        ]
    write_fasta(output_records)
"""


@unittest.skipUnless(VIRA_BINARY, "VIRA_TEST_BINARY is set by tests/run_tests.py")
class WorkflowTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="vira-workflow-")
        self.root = Path(self.temporary.name)
        self.reference = self.root / "reference.fasta"
        self.query = self.root / "query.fasta"
        self.reference.write_text(">reference\nACGTNACGT\n", encoding="utf-8")
        self.query.write_text(">query1\nACGTNACGT\n", encoding="utf-8")

        self.mock_bin = self.root / "mock-bin"
        self.mock_bin.mkdir()
        self.mock_mafft = self.mock_bin / "mafft"
        self.mock_mafft.write_text(textwrap.dedent(MOCK_MAFFT), encoding="utf-8")
        self.mock_mafft.chmod(
            stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
            stat.S_IRGRP | stat.S_IXGRP |
            stat.S_IROTH | stat.S_IXOTH
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def run_vira(
        self,
        output: Path,
        *extra_arguments: str,
        mode: str = "normal",
        path_override: str | None = None,
        timeout: float = 8.0,
        reference: Path | None = None,
        query: Path | None = None,
    ) -> subprocess.CompletedProcess[str]:
        environment = os.environ.copy()
        environment["VIRA_MOCK_MAFFT_MODE"] = mode
        environment["PATH"] = path_override or (
            str(self.mock_bin) + os.pathsep + environment.get("PATH", "")
        )
        command = [
            str(VIRA_BINARY),
            "-r", str(reference or self.reference),
            "-q", str(query or self.query),
            "-o", str(output),
            *extra_arguments,
        ]
        return subprocess.run(
            command,
            cwd=self.root,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )

    def assert_safe_failure(self, result: subprocess.CompletedProcess[str], message: str) -> None:
        self.assertGreater(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(message, result.stderr)
        self.assertNotIn("Segmentation fault", result.stderr)

    @staticmethod
    def fasta_ids(path: Path) -> list[str]:
        return [
            line[1:].split()[0]
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.startswith(">")
        ]

    @staticmethod
    def report_path(output: Path, suffix: str) -> Path:
        return Path(str(output.with_suffix("")) + suffix)

    @staticmethod
    def tsv_rows(path: Path) -> list[dict[str, str]]:
        with path.open(encoding="utf-8", newline="") as stream:
            return list(csv.DictReader(stream, delimiter="\t"))

    def make_two_query_fixture(self) -> Path:
        path = self.root / "two-query.fasta"
        path.write_text(">q1\nACGTNACGT\n>q2\nACGTNACGA\n", encoding="utf-8")
        return path

    def run_two_query_qc_case(
        self,
        output: Path,
        query: Path,
        *,
        mode: str = "normal",
    ) -> subprocess.CompletedProcess[str]:
        return self.run_vira(
            output,
            "--enable-qc",
            "--disable-auto-chunk",
            "--chunk-size", "1000",
            "--max-instances", "1",
            mode=mode,
            query=query,
        )

    def assert_chunk_failure_is_atomic(
        self,
        result: subprocess.CompletedProcess[str],
        output: Path,
        failure_reason: str,
        query_ids: tuple[str, ...] = ("q1", "q2"),
    ) -> None:
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn(failure_reason, result.stderr)

        failed_chunks = self.report_path(output, ".failed_chunks.tsv").read_text(encoding="utf-8")
        self.assertIn(failure_reason, failed_chunks)
        self.assertEqual(self.fasta_ids(output), ["reference"])
        self.assertEqual(
            self.fasta_ids(self.report_path(output, ".filtered.fasta")),
            ["reference"],
        )

        qc_rows = {row["seq_id"]: row for row in self.tsv_rows(self.report_path(output, ".qc_summary.tsv"))}
        unaligned_rows = {
            row["seq_id"]: row
            for row in self.tsv_rows(self.report_path(output, ".unaligned_samples.tsv"))
        }
        abnormal_rows = {
            row["seq_id"]: row
            for row in self.tsv_rows(self.report_path(output, ".abnormal_samples.tsv"))
        }
        sample_qc_ids = {
            row["seq_id"]
            for row in self.tsv_rows(self.report_path(output, ".sample_qc_summary.tsv"))
        }
        for query_id in query_ids:
            self.assertEqual(qc_rows[query_id]["alignment_status"], "CHUNK_FAILED")
            self.assertEqual(qc_rows[query_id]["post_alignment_qc"], "NA")
            self.assertEqual(unaligned_rows[query_id]["alignment_status"], "CHUNK_FAILED")
            self.assertEqual(abnormal_rows[query_id]["alignment_status"], "CHUNK_FAILED")
            self.assertEqual(abnormal_rows[query_id]["post_alignment_qc"], "NA")
            self.assertNotIn(query_id, sample_qc_ids)

    def test_normal_mock_alignment_succeeds(self) -> None:
        output = self.root / "normal.fasta"
        result = self.run_vira(output)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(
            output.read_text(encoding="utf-8"),
            self.reference.read_text(encoding="utf-8") + self.query.read_text(encoding="utf-8"),
        )
        self.assertTrue((self.root / "normal.qc_summary.tsv").is_file())

    def test_mock_alignment_with_mapout_qc_succeeds(self) -> None:
        output = self.root / "with-qc.fasta"
        result = self.run_vira(output, "--enable-qc")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        filtered = self.root / "with-qc.filtered.fasta"
        summary = self.root / "with-qc.sample_qc_summary.tsv"
        self.assertEqual(filtered.read_text(encoding="utf-8").count(">"), 2)
        self.assertIn("PASS", summary.read_text(encoding="utf-8"))
        self.assertEqual(
            self.report_path(output, ".qc_summary.tsv").read_text(encoding="utf-8").splitlines()[0],
            "original_id\tseq_id\theader\tbucket\tseq_len\tlen_diff_ratio\tn_ratio\tinput_decision\tinput_reasons\tneeds_manual_review\talignment_status\talignment_detail\tlost_count\tlost_ratio\tlost_block_count\tmax_lost_block_len\tpost_alignment_qc",
        )
        self.assertEqual(
            summary.read_text(encoding="utf-8").splitlines()[0],
            "seq_id\tbucket\tseq_len\tlen_diff_ratio\tn_ratio\tlost_count\tlost_ratio\tlost_block_count\tmax_lost_block_len\tqc_level",
        )

    def test_stdout_identity_exact_match_succeeds(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "identity-match.fasta"
        result = self.run_two_query_qc_case(output, query)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.fasta_ids(output), ["reference", "q1", "q2"])

    def test_stdout_reordered_queries_succeed(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "reordered-query-output.fasta"
        result = self.run_two_query_qc_case(output, query, mode="reverse-query-order")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.fasta_ids(output), ["reference", "q2", "q1"])

    def test_stdout_descriptive_headers_use_primary_ids(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "descriptive-output.fasta"
        result = self.run_two_query_qc_case(output, query, mode="descriptive-headers")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.fasta_ids(output), ["reference", "q1", "q2"])

    def test_stdout_replaced_reference_fails_closed(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "replaced-reference.fasta"
        result = self.run_two_query_qc_case(output, query, mode="reference-replaced")
        self.assert_chunk_failure_is_atomic(result, output, "MAFFT_OUTPUT_FORMAT_ERROR")
        self.assertIn("expected 'reference'", result.stderr)

    def test_stdout_replaced_query_fails_closed(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "replaced-query.fasta"
        result = self.run_two_query_qc_case(output, query, mode="query-replaced")
        self.assert_chunk_failure_is_atomic(result, output, "MAFFT_OUTPUT_FORMAT_ERROR")
        self.assertIn("missing=[q1]", result.stderr)
        self.assertIn("unexpected=[renamed_q1]", result.stderr)

    def test_stdout_identifiers_are_case_sensitive(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "case-changed-query.fasta"
        result = self.run_two_query_qc_case(output, query, mode="case-changed-query")
        self.assert_chunk_failure_is_atomic(result, output, "MAFFT_OUTPUT_FORMAT_ERROR")
        self.assertIn("missing=[q1]", result.stderr)
        self.assertIn("unexpected=[Q1]", result.stderr)

    def test_stdout_duplicate_query_fails_closed(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "duplicate-query-output.fasta"
        result = self.run_two_query_qc_case(output, query, mode="duplicate-query")
        self.assert_chunk_failure_is_atomic(result, output, "MAFFT_OUTPUT_FORMAT_ERROR")
        self.assertIn("missing=[q2]", result.stderr)
        self.assertIn("duplicates=[q1]", result.stderr)

    def test_stdout_unexpected_query_fails_closed(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "unexpected-query.fasta"
        result = self.run_two_query_qc_case(output, query, mode="unexpected-query")
        self.assert_chunk_failure_is_atomic(result, output, "MAFFT_OUTPUT_FORMAT_ERROR")
        self.assertIn("missing=[q2]", result.stderr)
        self.assertIn("unexpected=[rogue]", result.stderr)

    def test_stdout_mismatch_diagnostics_are_sorted(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "sorted-stdout-diagnostic.fasta"
        result = self.run_two_query_qc_case(output, query, mode="multiple-unexpected-queries")
        self.assert_chunk_failure_is_atomic(result, output, "MAFFT_OUTPUT_FORMAT_ERROR")
        self.assertIn("missing=[q1,q2]", result.stderr)
        self.assertIn("unexpected=[alpha,zeta]", result.stderr)

    def test_stdout_reference_must_be_first(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "reference-not-first.fasta"
        result = self.run_two_query_qc_case(output, query, mode="reference-not-first")
        self.assert_chunk_failure_is_atomic(result, output, "MAFFT_OUTPUT_FORMAT_ERROR")
        self.assertIn("first record", result.stderr)

    def test_mapout_reordered_sections_succeed(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "reordered-mapout.fasta"
        result = self.run_two_query_qc_case(output, query, mode="mapout-reverse-sections")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(
            {row["seq_id"] for row in self.tsv_rows(self.report_path(output, ".sample_qc_summary.tsv"))},
            {"q1", "q2"},
        )

    def test_mapout_missing_section_fails_closed(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "missing-mapout-section.fasta"
        result = self.run_two_query_qc_case(output, query, mode="mapout-missing-section")
        self.assert_chunk_failure_is_atomic(result, output, "MAPOUT_FORMAT_ERROR")
        self.assertIn("missing=[q2]", result.stderr)
        self.assertIn("unexpected=[]", result.stderr)

    def test_mapout_replaced_section_fails_closed(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "replaced-mapout-section.fasta"
        result = self.run_two_query_qc_case(output, query, mode="mapout-replaced-section")
        self.assert_chunk_failure_is_atomic(result, output, "MAPOUT_FORMAT_ERROR")
        self.assertIn("missing=[q2]", result.stderr)
        self.assertIn("unexpected=[rogue]", result.stderr)

    def test_mapout_extra_section_fails_closed(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "extra-mapout-section.fasta"
        result = self.run_two_query_qc_case(output, query, mode="mapout-extra-section")
        self.assert_chunk_failure_is_atomic(result, output, "MAPOUT_FORMAT_ERROR")
        self.assertIn("missing=[]", result.stderr)
        self.assertIn("unexpected=[rogue]", result.stderr)

    def test_mapout_duplicate_section_fails_closed(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "duplicate-mapout-section.fasta"
        result = self.run_two_query_qc_case(output, query, mode="mapout-duplicate-section")
        self.assert_chunk_failure_is_atomic(result, output, "MAPOUT_FORMAT_ERROR")
        self.assertIn("duplicate sequence section identifier 'q1'", result.stderr)

    def test_duplicate_rename_collision_produces_globally_unique_ids(self) -> None:
        query = self.root / "rename-collision.fasta"
        query.write_text(
            ">x\nACGTNACGT\n>x\nACGTNACGA\n>x__dup2\nACGTNACGG\n",
            encoding="utf-8",
        )
        output = self.root / "rename-collision-output.fasta"
        result = self.run_vira(
            output,
            "--duplicate-id-policy", "rename",
            "--disable-auto-chunk",
            "--chunk-size", "1000",
            "--max-instances", "1",
            query=query,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(
            self.fasta_ids(output),
            ["reference", "x", "x__dup2", "x__dup2__dup2"],
        )

    def test_query_identifier_conflicting_with_reference_is_excluded(self) -> None:
        query = self.root / "reference-conflict.fasta"
        query.write_text(">reference\nACGTNACGT\n", encoding="utf-8")
        output = self.root / "reference-conflict-output.fasta"
        result = self.run_vira(output, query=query)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.fasta_ids(output), ["reference"])
        row = self.tsv_rows(self.report_path(output, ".qc_summary.tsv"))[0]
        self.assertEqual(row["input_decision"], "EXCLUDE")
        self.assertEqual(row["alignment_status"], "NOT_DISPATCHED")
        self.assertIn("DUPLICATE_ID", row["input_reasons"])

    def test_empty_query_is_a_successful_boundary_case(self) -> None:
        empty_query = self.root / "empty.fasta"
        empty_query.write_text("", encoding="utf-8")
        output = self.root / "empty-result.fasta"
        result = self.run_vira(output, query=empty_query)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(output.read_text(encoding="utf-8").count(">"), 1)

    def test_missing_mafft_is_reported_without_a_crash(self) -> None:
        empty_path = self.root / "empty-path"
        empty_path.mkdir()
        result = self.run_vira(self.root / "missing.fasta", path_override=str(empty_path))
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("MAFFT_EXEC_NOT_FOUND", result.stderr)

    def test_non_executable_mafft_is_distinguished(self) -> None:
        non_executable_path = self.root / "non-executable-path"
        non_executable_path.mkdir()
        candidate = non_executable_path / "mafft"
        candidate.write_text("not executable\n", encoding="utf-8")
        candidate.chmod(stat.S_IRUSR | stat.S_IWUSR)
        result = self.run_vira(
            self.root / "not-executable.fasta",
            path_override=str(non_executable_path),
        )
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("MAFFT_EXEC_NOT_EXECUTABLE", result.stderr)

    def test_mafft_timeout_terminates_the_process_group(self) -> None:
        result = self.run_vira(
            self.root / "timeout.fasta",
            "--mafft-timeout-sec", "1",
            mode="sleep",
            timeout=6.0,
        )
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("TIMEOUT", result.stderr)

    def test_malformed_mafft_stdout_becomes_a_chunk_failure(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "bad-output.fasta"
        result = self.run_two_query_qc_case(output, query, mode="malformed-output")
        self.assert_chunk_failure_is_atomic(result, output, "MAFFT_OUTPUT_FORMAT_ERROR")

    def test_malformed_mapout_becomes_a_chunk_failure(self) -> None:
        query = self.make_two_query_fixture()
        output = self.root / "bad-mapout.fasta"
        result = self.run_two_query_qc_case(output, query, mode="malformed-mapout")
        self.assert_chunk_failure_is_atomic(result, output, "MAPOUT_FORMAT_ERROR")
        self.assertIn("line 2", result.stderr)

    def test_nonzero_mafft_exit_is_reported(self) -> None:
        result = self.run_vira(self.root / "nonzero.fasta", mode="nonzero")
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("MAFFT_EXIT_NONZERO", result.stderr)

    def test_output_cannot_overwrite_an_input(self) -> None:
        before = self.reference.read_bytes()
        result = self.run_vira(self.reference)
        self.assert_safe_failure(result, "Refusing to overwrite an input file")
        self.assertEqual(self.reference.read_bytes(), before)

    def test_auxiliary_hardlink_cannot_overwrite_final_output(self) -> None:
        output = self.root / "aliased.fasta"
        output.write_text(">sentinel\nACGT\n", encoding="utf-8")
        report = self.root / "aliased.qc_summary.tsv"
        os.link(output, report)
        before = output.read_bytes()
        result = self.run_vira(output)
        self.assert_safe_failure(result, "output paths must be distinct")
        self.assertEqual(output.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
