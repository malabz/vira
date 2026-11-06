# ViRA automated tests

The test suite is self-contained and uses only a C++17 compiler and Python 3.8+.
MAFFT is not required: process-boundary tests create an isolated mock executable.
All binaries, fixtures, temporary directories, and generated reports are created
under the operating system's temporary directory and removed after the run.

Run every test from the repository root:

```bash
python3 -B tests/run_tests.py
```

Use a specific compiler when needed:

```bash
CXX=clang++ python3 -B tests/run_tests.py
```

The runner performs the following steps:

1. Builds the complete `vira` executable with warnings enabled.
2. Builds and runs the standalone C++ Mapout parser tests.
3. Discovers and runs every `tests/test_*.py` integration test.

The suite covers strict CLI parsing, malformed and compressed sequence inputs,
Mapout parser failures, MAFFT launch/timeout/output failures, output alias
protection, and successful alignment/QC boundary cases.

To run only the Python integration tests against an existing binary:

```bash
VIRA_TEST_BINARY=/absolute/path/to/vira \
  python3 -B -m unittest discover -s tests -p 'test_*.py' -v
```
