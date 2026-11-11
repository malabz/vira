# Known Issues

This document records known limitations of the current ViRA development snapshot.
It is a development candidate, not a stable release.

## Resolved P1 integrity issues

The following issues were identified during an independent acceptance audit and
are resolved in the current development snapshot. Both fixes use fail-closed chunk
handling: no query alignment or post-alignment QC payload is committed until every
enabled identity check succeeds.

### 1. MAFFT output record identities are verified

**Former severity:** P1 — output integrity

MAFFT stdout retains the existing structural, NUL, empty-record, and record-count
checks. It now additionally enforces these invariants:

- the first record's primary identifier exactly matches the parsed reference ID;
- query primary identifiers form exactly the set dispatched in the chunk;
- identifier comparison is case-sensitive and independent of query output order;
- duplicate, missing, unexpected, renamed, and reference-conflicting IDs fail;
- expected query IDs must themselves be unique; and
- full header descriptions are ignored after the first whitespace-delimited token.

FASTA input, MAFFT stdout, and mapout sections share one primary-ID extraction
rule. Input QC also reserves effective IDs globally, including the reference ID,
so duplicate renaming cannot create collisions such as `x`, `x__dup2`, and a later
original `x__dup2`.

Identity failures use `MAFFT_OUTPUT_FORMAT_ERROR`. The chunk remains failed, its
stdout is not appended, and its expected samples are reported as `CHUNK_FAILED`
rather than `ALIGNED`.

Automated workflow coverage includes correct IDs, replaced reference/query IDs,
duplicate queries, missing and unexpected IDs, reference-first enforcement,
order-independent query sets, descriptive headers, reference/query conflicts, and
duplicate-rename collisions. Failure tests also verify that no query from the
failed chunk reaches the full or filtered MSA.

### 2. Mapout section identities are matched to expected chunk IDs

**Former severity:** P1 — post-alignment QC integrity

When `--enable-qc` is active, mapout section primary identifiers must form exactly
the query-ID set dispatched in that chunk. Section order is irrelevant, while
missing, unexpected, replaced, extra, or duplicate primary IDs raise
`MapoutFormatError` and retain the stable `MAPOUT_FORMAT_ERROR` failure category.
Diagnostics render sorted `missing=[...]` and `unexpected=[...]` sets.

The chunk success path is atomic: stdout validation, mapout parsing, and mapout
identity validation all finish before stdout, deleted blocks, statistics, or the
success flag are stored in `ChunkResult`. A failed chunk therefore contributes no
accepted alignment or QC payload.

`build_sample_qc_rows()` now omits samples without validated mapout statistics.
Consequently, failed samples have `NA` post-alignment QC in `qc_summary.tsv`, do
not appear as zero-loss `PASS` rows in `sample_qc_summary.tsv`, and cannot enter the
filtered MSA through missing metadata.

Automated parser and workflow coverage includes complete, missing, replaced,
extra, reordered, descriptive, and duplicate sections, duplicate expected IDs,
atomic failure reports, and the absence of a default PASS when statistics are
missing.

### Remaining mapout limitation

Exact section-identity validation proves that every expected section exists and no
unexpected section is present. It does not prove that every coordinate row inside
an otherwise valid section is complete. Internal-row truncation remains outside
the scope of this fix because its validation requires confirmed MAFFT coordinate
semantics for query sequences containing `-` and `.`.

## Real-MAFFT acceptance coverage

The original independent, read-only acceptance run compared the pre-improvement
baseline commit `4e6a03e0ba18ba49480a5e6c9d65c62d92828604` with the pre-fix development
snapshot. After the identity fixes above, the same scenario matrix was rerun by
comparing committed pre-fix snapshot
`e0b5a7ed1a5923856f36105a80ac02e393734fda` with the fixed working tree.

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

In the post-fix rerun, single-process outputs and all QC artifacts were
byte-for-byte identical to the pre-fix snapshot. Concurrent final MSAs were either
byte-for-byte identical or differed only in the existing completion-based record
order; strict normalization by primary ID confirmed identical full headers and
aligned sequence content. Real mapout sections passed the new exact-ID validation,
and the 4 bp and 6 bp deletion-block results remained unchanged.

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
- exact mapout section identity does not establish completeness of coordinate rows
  inside each section, as described above.

In the audit environment, MAFFT's shell wrapper could not access `/dev/stderr`
inside the restricted sandbox. Valid acceptance runs were therefore executed in
an approved normal process context while keeping generated files in an external
temporary directory. The sandbox-specific failure was not classified as a ViRA or
MAFFT functional regression.
