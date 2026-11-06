# Known Issues

This document records known limitations of the current ViRA development snapshot.
It is a development candidate, not a stable release.

## Open P1 integrity issues

The following issues remain intentionally unfixed in this snapshot. They were
identified during an independent acceptance audit and are not resolved by the
successful regression runs described below.

### 1. MAFFT output record identities are not verified

**Severity:** P1 — output integrity

`validate_mafft_fasta_output()` verifies that MAFFT stdout:

- is structurally valid FASTA;
- contains no NUL bytes;
- has non-empty headers and sequences; and
- contains exactly one reference record plus the expected number of query records.

It does not verify that the returned reference and query identifiers match the
identifiers submitted for that chunk, nor does it reject duplicate or unexpected
output identifiers.

A crafted or corrupted MAFFT-compatible process can therefore return the correct
number of FASTA records under incorrect or duplicate headers and still pass the
current validation. `handle_result()` then marks the expected input identifiers as
`ALIGNED`, even when those identifiers are absent from the appended alignment.

**Potential impact:** sequences can be misidentified or duplicated in the final
MSA while the corresponding expected samples are reported as successfully aligned.

**Conceptual reproduction:** replace MAFFT with a compatible executable that exits
with status 0 and emits the expected number of valid FASTA records, but changes or
duplicates the query headers. The chunk can be accepted because its record count
and FASTA structure remain valid.

**Current partial safeguards:** non-zero MAFFT exits, timeouts, malformed FASTA,
empty records, NUL bytes, and incorrect record counts are rejected. These checks do
not establish record identity.

**Operational mitigation:** use a trusted MAFFT installation and independently
compare the final MSA identifiers with the input identifiers. Do not treat the
`ALIGNED` status alone as proof of output identity.

### 2. Mapout section identities are not matched to expected chunk IDs

**Severity:** P1 — post-alignment QC integrity

The mapout parser validates section syntax, field counts, numeric values, NUL bytes,
and duplicate section identifiers. After parsing, however, the set of mapout
section identifiers is not compared with the set of query identifiers dispatched
in the corresponding chunk.

If a syntactically valid mapout file is truncated or omits one or more expected
sections, the chunk can still be accepted. During sample-QC row construction, a
sample with no parsed mapout statistics receives default zero-valued statistics
and can consequently be reported as `PASS`. Unexpected section identifiers are
likewise not rejected by an exact expected-set comparison.

**Potential impact:** missing or mismatched mapout data can cause affected samples
to receive an incorrect zero-loss `PASS` result.

**Conceptual reproduction:** for a chunk containing multiple query identifiers,
provide a syntactically valid mapout file containing only one expected section.
Parsing can succeed while omitted samples receive default zero-loss statistics.

**Current partial safeguards:** missing mapout files, files with no sections,
duplicate sections, malformed rows, invalid coordinates, and I/O failures are
rejected. These checks do not prove that the section set is complete and belongs
to the submitted chunk.

**Operational mitigation:** independently verify that mapout data covers every
aligned query represented in `*.sample_qc_summary.tsv`. Treat zero-loss `PASS`
results cautiously when mapout completeness has not been independently established.

## Real-MAFFT acceptance coverage

An independent, read-only acceptance run compared the pre-improvement baseline
commit `4e6a03e0ba18ba49480a5e6c9d65c62d92828604` with this development working tree.

The runs used:

- Conda environment `msa`;
- MAFFT `v7.526 (2024/Apr/26)`;
- identical inputs, ViRA arguments, MAFFT version, and thread settings for both
  revisions; and
- separately compiled baseline and development binaries in an external temporary
  directory.

The exercised scenarios included:

- three queries in one chunk with one MAFFT process;
- eight queries split across four chunks with one MAFFT process;
- the same four-chunk input with three concurrent MAFFT processes;
- three repeated multi-process runs;
- inputs containing substitutions, insertions, and terminally missing sequence;
- input-QC classifications including duplicate identifiers, high-N content,
  illegal characters, and an engineered-sequence header;
- real MAFFT `--mapout` processing across multiple chunks; and
- non-trivial mapout deletion blocks of 4 bp and 6 bp, including QC failure and
  filtered-MSA generation.

For these fixtures:

- baseline and development runs exited successfully;
- no MAFFT chunk failed or timed out;
- compared final MSAs were byte-for-byte identical;
- single-process and multi-process results were byte-for-byte identical;
- input-QC and post-alignment QC reports were byte-for-byte identical; and
- the non-trivial mapout case produced identical deletion blocks, QC levels, and
  filtered alignments.

### Acceptance limitations

This evidence is intentionally limited:

- it covers MAFFT 7.526 only;
- it uses small synthetic fixtures rather than a representative production corpus;
- it does not establish behavior across other operating systems, compilers, or
  MAFFT versions;
- it does not cover large-scale performance, resource exhaustion, every timeout
  path, or every possible malformed child-process output;
- the real-MAFFT regression matrix was executed as an external acceptance exercise
  and is not part of the repository's automated test suite or CI; and
- successful regression output does not resolve either P1 integrity issue listed
  above.

In the audit environment, MAFFT's shell wrapper could not access `/dev/stderr`
inside the restricted sandbox. Valid acceptance runs were therefore executed in
an approved normal process context while keeping generated files in an external
temporary directory. The sandbox-specific failure was not classified as a ViRA or
MAFFT functional regression.
