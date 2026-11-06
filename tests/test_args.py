"""Black-box edge-case tests for ViRA command-line argument validation."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import tempfile
import unittest


class ArgsValidationTests(unittest.TestCase):
    """Exercise strict numeric parsing and cross-option validation."""

    TIMEOUT_SECONDS = 5

    @classmethod
    def setUpClass(cls) -> None:
        binary = os.environ.get("VIRA_TEST_BINARY")
        if not binary:
            raise RuntimeError("VIRA_TEST_BINARY must point to the compiled ViRA executable")

        cls.binary = Path(binary).resolve()
        if not cls.binary.is_file() or not os.access(cls.binary, os.X_OK):
            raise RuntimeError(f"VIRA_TEST_BINARY is not executable: {cls.binary}")

    def setUp(self) -> None:
        self._temporary_directory = tempfile.TemporaryDirectory(prefix="vira-args-test-")
        self.workdir = Path(self._temporary_directory.name)
        self.reference = self.workdir / "reference.fasta"
        self.query = self.workdir / "query.fasta"
        self.output = self.workdir / "output.fasta"
        self.reference.write_text(">reference\nACGT\n", encoding="ascii")
        self.query.write_text("", encoding="ascii")

    def tearDown(self) -> None:
        self._temporary_directory.cleanup()

    def run_cli(self, *extra_args: str) -> subprocess.CompletedProcess[str]:
        command = [
            str(self.binary),
            "--ref",
            str(self.reference),
            "--query",
            str(self.query),
            "--out",
            str(self.output),
            *extra_args,
        ]
        try:
            return subprocess.run(
                command,
                cwd=self.workdir,
                capture_output=True,
                text=True,
                timeout=self.TIMEOUT_SECONDS,
                check=False,
            )
        except subprocess.TimeoutExpired as error:
            self.fail(f"ViRA hung while validating arguments: {error.cmd}")

    def assert_invalid(self, args: tuple[str, ...], expected_message: str) -> None:
        """Assert a parser failure is contextual, prompt, and not signal-driven."""
        result = self.run_cli(*args)
        diagnostics = result.stdout + result.stderr
        self.assertEqual(result.returncode, 1, diagnostics)
        self.assertGreater(result.returncode, 0, "process was terminated by a signal")
        self.assertIn("Error:", diagnostics)
        self.assertIn(expected_message, diagnostics)
        self.assertNotIn("Segmentation fault", diagnostics)

    def test_rejects_nonfinite_and_overflowing_floats(self) -> None:
        for value in ("NaN", "Inf", "-Inf", "1e309"):
            with self.subTest(value=value):
                self.assert_invalid(
                    ("--slow-chunk-sec", value),
                    "Invalid number for --slow-chunk-sec",
                )

    def test_rejects_partially_parsed_numeric_values(self) -> None:
        cases = (
            (("--slow-chunk-sec", "1.5seconds"), "Invalid number for --slow-chunk-sec"),
            (("--max-instances", "12x"), "Invalid integer for --max-instances"),
            (("--threads-per-mafft", "+1"), "Invalid integer for --threads-per-mafft"),
            (("--warn-lost-ratio", " 0.1"), "Invalid number for --warn-lost-ratio"),
        )
        for args, expected in cases:
            with self.subTest(args=args):
                self.assert_invalid(args, expected)

    def test_rejects_integer_and_size_overflow(self) -> None:
        cases = (
            (
                ("--max-instances", "2147483648"),
                "Value for --max-instances is outside the supported integer range",
            ),
            (
                ("--target-chunk-bases", "18446744073709551616"),
                "Value for --target-chunk-bases is outside the supported non-negative integer range",
            ),
            (
                ("--target-chunk-bases", "-1"),
                "Invalid non-negative integer for --target-chunk-bases",
            ),
            (
                ("--max-instances", "2147483647", "--max-pending-factor", "2"),
                "product of --max-instances and --max-pending-factor exceeds the supported range",
            ),
        )
        for args, expected in cases:
            with self.subTest(args=args):
                self.assert_invalid(args, expected)

    def test_rejects_zero_for_positive_integer_options(self) -> None:
        cases = (
            (("--max-instances", "0"), "--max-instances must be at least 1"),
            (("--threads-per-mafft", "0"), "--threads-per-mafft must be at least 1"),
            (("--chunk-size", "0"), "--chunk-size must be at least 1"),
            (("--max-pending-factor", "0"), "--max-pending-factor must be at least 1"),
        )
        for args, expected in cases:
            with self.subTest(args=args):
                self.assert_invalid(args, expected)

    def test_rejects_negative_runtime_thresholds(self) -> None:
        cases = (
            (("--mafft-timeout-sec", "-1"), "--mafft-timeout-sec must be 0 or greater"),
            (("--slow-chunk-sec", "-0.01"), "--slow-chunk-sec must be 0 or greater"),
        )
        for args, expected in cases:
            with self.subTest(args=args):
                self.assert_invalid(args, expected)

    def test_rejects_invalid_input_qc_ratio_ranges(self) -> None:
        cases = (
            (("--review-n-ratio", "-0.1"), "N-ratio thresholds must be within [0, 1]"),
            (("--exclude-n-ratio", "1.1"), "N-ratio thresholds must be within [0, 1]"),
            (
                ("--review-n-ratio", "0.9", "--exclude-n-ratio", "0.8"),
                "--review-n-ratio cannot exceed --exclude-n-ratio",
            ),
            (
                ("--exclude-min-len-ratio", "2", "--exclude-max-len-ratio", "1"),
                "Invalid exclude length-ratio thresholds",
            ),
            (
                ("--reasonable-min-len-ratio", "2", "--reasonable-max-len-ratio", "1"),
                "Invalid reasonable length-ratio thresholds",
            ),
            (
                ("--reasonable-min-len-ratio", "0.05"),
                "Reasonable length-ratio thresholds must lie within the exclusion range",
            ),
        )
        for args, expected in cases:
            with self.subTest(args=args):
                self.assert_invalid(args, expected)

    def test_rejects_invalid_post_alignment_qc_thresholds(self) -> None:
        cases = (
            (("--warn-lost-ratio", "-0.1"), "lost-ratio thresholds must be within [0, 1]"),
            (("--fail-lost-ratio", "1.1"), "lost-ratio thresholds must be within [0, 1]"),
            (
                ("--warn-lost-ratio", "0.5", "--fail-lost-ratio", "0.4"),
                "--warn-lost-ratio cannot exceed --fail-lost-ratio",
            ),
            (("--warn-max-block-len", "-1"), "block-length thresholds must be non-negative"),
            (
                ("--warn-max-block-len", "101", "--fail-max-block-len", "100"),
                "--warn-max-block-len cannot exceed --fail-max-block-len",
            ),
        )
        for args, expected in cases:
            with self.subTest(args=args):
                self.assert_invalid(args, expected)

    def test_rejects_missing_values_unknown_options_and_empty_alphabet(self) -> None:
        cases = (
            (("--chunk-size",), "Missing value for option --chunk-size"),
            (("--not-a-real-option",), "Unknown option: --not-a-real-option"),
            (("--allowed-seq-chars", ""), "--allowed-seq-chars must not be empty"),
            (("--duplicate-id-policy", "invalid"), "--duplicate-id-policy must be exclude or rename"),
        )
        for args, expected in cases:
            with self.subTest(args=args):
                self.assert_invalid(args, expected)

    def test_accepts_exact_numeric_boundaries(self) -> None:
        result = self.run_cli(
            "--target-chunk-bases",
            "0",
            "--mafft-timeout-sec",
            "0",
            "--slow-chunk-sec",
            "0",
            "--review-n-ratio",
            "0",
            "--exclude-n-ratio",
            "1",
            "--warn-lost-ratio",
            "0",
            "--fail-lost-ratio",
            "1",
            "--warn-max-block-len",
            "0",
            "--fail-max-block-len",
            "0",
        )
        diagnostics = result.stdout + result.stderr
        self.assertEqual(result.returncode, 0, diagnostics)
        self.assertEqual(self.output.read_text(encoding="ascii"), ">reference\nACGT\n")


if __name__ == "__main__":
    unittest.main()
