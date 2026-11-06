"""Black-box negative tests for sequence-file boundary validation."""

from __future__ import annotations

import gzip
import io
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import zipfile


class NegativeSequenceInputTests(unittest.TestCase):
    """Verify malformed inputs fail cleanly before alignment starts."""

    @classmethod
    def setUpClass(cls) -> None:
        binary = os.environ.get("VIRA_TEST_BINARY")
        if not binary:
            raise RuntimeError("VIRA_TEST_BINARY must point to the compiled ViRA executable")
        cls.binary = Path(binary).resolve()
        if not cls.binary.is_file():
            raise RuntimeError(f"ViRA test executable does not exist: {cls.binary}")

    def setUp(self) -> None:
        self._temporary_directory = tempfile.TemporaryDirectory(prefix="vira-negative-")
        self.workdir = Path(self._temporary_directory.name)
        self.valid_reference = self.workdir / "reference.fasta"
        self.valid_query = self.workdir / "query.fasta"
        self.valid_reference.write_text(">reference\nACGTNACGT\n", encoding="ascii")
        self.valid_query.write_text(">query1\nACGTNACGT\n", encoding="ascii")

    def tearDown(self) -> None:
        self._temporary_directory.cleanup()

    def assert_rejected(
        self,
        *,
        reference: Path | None = None,
        query: Path | None = None,
        expected_messages: tuple[str, ...],
    ) -> str:
        """Run ViRA and assert a bounded, non-signal failure with useful context."""

        output = self.workdir / "result.fasta"
        command = [
            str(self.binary),
            "--ref",
            str(reference or self.valid_reference),
            "--query",
            str(query or self.valid_query),
            "--out",
            str(output),
        ]
        try:
            completed = subprocess.run(
                command,
                cwd=self.workdir,
                capture_output=True,
                text=True,
                timeout=10,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            self.fail(f"ViRA hung for malformed input: {exc}")

        diagnostic = completed.stdout + completed.stderr
        self.assertGreater(
            completed.returncode,
            0,
            msg=(
                "Malformed input must produce a positive application exit code, "
                f"not success or signal termination. Output:\n{diagnostic}"
            ),
        )
        for message in expected_messages:
            self.assertIn(message, diagnostic)
        self.assertFalse(output.exists(), "A rejected input must not create the final output")
        return diagnostic

    def test_rejects_compressed_extensions(self) -> None:
        gzip_query = self.workdir / "query.fasta.gz"
        gzip_query.write_bytes(gzip.compress(b">query1\nACGT\n"))
        zip_query = self.workdir / "query.fasta.zip"
        with zipfile.ZipFile(zip_query, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr("query.fasta", ">query1\nACGT\n")

        for extension, query in ((".gz", gzip_query), (".zip", zip_query)):
            with self.subTest(extension=extension):
                self.assert_rejected(
                    query=query,
                    expected_messages=(
                        f"compressed or archive extension '{extension}'",
                        "Please decompress it first.",
                    ),
                )

    def test_rejects_compressed_content_with_fasta_suffix(self) -> None:
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr("query.fasta", ">query1\nACGT\n")
        signatures = {
            "gzip": (gzip.compress(b">query1\nACGT\n"), "gzip-compressed"),
            "zip": (zip_buffer.getvalue(), "ZIP archive"),
        }
        for label, (content, expected_format) in signatures.items():
            with self.subTest(format=label):
                query = self.workdir / f"disguised-{label}.fasta"
                query.write_bytes(content)
                self.assert_rejected(
                    query=query,
                    expected_messages=(
                        f"detected {expected_format} content",
                        "Please decompress it first.",
                    ),
                )

    def test_rejects_unexpected_suffix_even_when_content_is_fasta(self) -> None:
        query = self.workdir / "query.txt"
        query.write_text(">query1\nACGT\n", encoding="ascii")
        self.assert_rejected(
            query=query,
            expected_messages=(
                "Expected a FASTA file",
                "with extension .fa, .fasta, .fas, or .fna",
            ),
        )

    def test_rejects_sequence_data_before_first_header(self) -> None:
        query = self.workdir / "damaged.fasta"
        query.write_text("query1\nACGT\n", encoding="ascii")
        self.assert_rejected(
            query=query,
            expected_messages=(
                "Invalid FASTA format",
                "at line 1",
                "expected a header line starting with '>'",
            ),
        )

    def test_rejects_empty_reference_inputs(self) -> None:
        cases = {
            "empty-file": (b"", "Reference FASTA is empty or contains no sequence record"),
            "header-only": (b">reference\n", "Reference FASTA sequence is empty"),
        }
        for label, (content, expected_message) in cases.items():
            with self.subTest(case=label):
                reference = self.workdir / f"{label}.fasta"
                reference.write_bytes(content)
                self.assert_rejected(
                    reference=reference,
                    expected_messages=(expected_message,),
                )

    def test_rejects_illegal_reference_symbol(self) -> None:
        reference = self.workdir / "illegal-symbol.fasta"
        reference.write_text(">reference\nACGT?ACGT\n", encoding="ascii")
        self.assert_rejected(
            reference=reference,
            expected_messages=(
                "Reference FASTA contains unsupported sequence symbol '?'",
                "at position 5",
                "Allowed symbols:",
            ),
        )

    def test_rejects_embedded_nul_byte(self) -> None:
        query = self.workdir / "nul-byte.fasta"
        query.write_bytes(b">query1\nACG\x00T\n")
        self.assert_rejected(
            query=query,
            expected_messages=(
                "Invalid FASTA format",
                "embedded NUL byte detected",
                "binary or corrupted",
            ),
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
