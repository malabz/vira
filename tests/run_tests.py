#!/usr/bin/env python3
"""Build ViRA in an isolated directory and run the complete test suite."""

from __future__ import annotations

import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile


sys.dont_write_bytecode = True


def run_command(command: list[str], *, cwd: Path, env: dict[str, str] | None = None) -> int:
    """Run one visible test command and return its process status."""
    print(f"+ {shlex.join(command)}", flush=True)
    try:
        completed = subprocess.run(command, cwd=cwd, env=env, check=False, timeout=180)
    except subprocess.TimeoutExpired:
        print(f"Error: command timed out: {shlex.join(command)}", file=sys.stderr)
        return 124
    return completed.returncode


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    test_dir = repo_root / "tests"
    compiler = shlex.split(os.environ.get("CXX", "c++"))
    common_flags = [
        "-std=c++17",
        "-O0",
        "-g",
        "-pthread",
        "-Wall",
        "-Wextra",
        "-Wpedantic",
        "-I",
        str(repo_root / "include"),
    ]

    with tempfile.TemporaryDirectory(prefix="vira-tests-") as build_text:
        build_dir = Path(build_text)
        vira_binary = build_dir / "vira"
        qc_test_binary = build_dir / "qc_parser_tests"

        application_sources = sorted(str(path) for path in (repo_root / "src").glob("*.cpp"))
        build_application = compiler + common_flags + application_sources + ["-o", str(vira_binary)]
        if run_command(build_application, cwd=repo_root) != 0:
            return 1

        build_qc_tests = compiler + common_flags + [
            str(test_dir / "qc_parser_tests.cpp"),
            str(repo_root / "src" / "qc.cpp"),
            "-o",
            str(qc_test_binary),
        ]
        if run_command(build_qc_tests, cwd=repo_root) != 0:
            # GCC 8 exposes C++17 filesystem behind a separate compatibility library.
            build_qc_tests_with_legacy_filesystem = build_qc_tests + ["-lstdc++fs"]
            if run_command(build_qc_tests_with_legacy_filesystem, cwd=repo_root) != 0:
                return 1

        if run_command([str(qc_test_binary)], cwd=build_dir) != 0:
            return 1

        test_environment = os.environ.copy()
        test_environment["VIRA_TEST_BINARY"] = str(vira_binary)
        test_environment["PYTHONDONTWRITEBYTECODE"] = "1"
        python_tests = [
            sys.executable,
            "-B",
            "-m",
            "unittest",
            "discover",
            "-s",
            str(test_dir),
            "-p",
            "test_*.py",
            "-v",
        ]
        return 0 if run_command(python_tests, cwd=repo_root, env=test_environment) == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
