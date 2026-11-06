"""End-to-end tests for MAFFT process boundaries and output protection."""

from __future__ import annotations

import os
from pathlib import Path
import stat
import subprocess
import tempfile
import textwrap
import unittest


VIRA_BINARY = os.environ.get("VIRA_TEST_BINARY")


MOCK_MAFFT = r"""#!/usr/bin/env python3
import os
from pathlib import Path
import sys
import time

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

if "--mapout" in arguments:
    map_path = Path(str(query_path) + ".map")
    if mode == "malformed-mapout":
        map_path.write_text(">query1\nA,not-an-integer,1\n", encoding="utf-8")
    else:
        records = []
        header = None
        sequence = []
        for raw_line in query_path.read_text(encoding="utf-8").splitlines():
            if raw_line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence)))
                header = raw_line[1:].split()[0]
                sequence = []
            else:
                sequence.append(raw_line.strip())
        if header is not None:
            records.append((header, "".join(sequence)))
        with map_path.open("w", encoding="utf-8") as stream:
            for record_id, record_sequence in records:
                stream.write(f">{record_id}\n")
                for position, symbol in enumerate(record_sequence, start=1):
                    stream.write(f"{symbol},{position},{position}\n")

sys.stdout.buffer.write(reference_path.read_bytes())
sys.stdout.buffer.write(query_path.read_bytes())
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

    def test_normal_mock_alignment_succeeds(self) -> None:
        output = self.root / "normal.fasta"
        result = self.run_vira(output)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(output.read_text(encoding="utf-8").count(">"), 2)
        self.assertTrue((self.root / "normal.qc_summary.tsv").is_file())

    def test_mock_alignment_with_mapout_qc_succeeds(self) -> None:
        output = self.root / "with-qc.fasta"
        result = self.run_vira(output, "--enable-qc")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        filtered = self.root / "with-qc.filtered.fasta"
        summary = self.root / "with-qc.sample_qc_summary.tsv"
        self.assertEqual(filtered.read_text(encoding="utf-8").count(">"), 2)
        self.assertIn("PASS", summary.read_text(encoding="utf-8"))

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
        result = self.run_vira(
            self.root / "bad-output.fasta",
            mode="malformed-output",
        )
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("MAFFT_OUTPUT_FORMAT_ERROR", result.stderr)

    def test_malformed_mapout_becomes_a_chunk_failure(self) -> None:
        result = self.run_vira(
            self.root / "bad-mapout.fasta",
            "--enable-qc",
            mode="malformed-mapout",
        )
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("MAPOUT_FORMAT_ERROR", result.stderr)
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
