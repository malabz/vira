# ViRA: A scalable viral reference aligner for accurate large-scale genome alignment and variant analysis

vira is a tool written in C++17 for viral-genome multiple sequence alignment.

## 🔨Installation and Usage

<!--
### 1.1 Linux/WSL(Windows Subsystem for Linux) - install dependency

vira requires MAFFT. We recommend using Conda to install MAFFT.

```bash
# 1 Create and activate a conda environment
conda create -n vira_env
conda activate vira_env

# 2 Install MAFFT
conda install -c bioconda mafft
```
-->

### 1.1 Linux/WSL(Windows Subsystem for Linux) - from the source code

1.Download and compile the source code. Make sure your compiler supports C++17.

```bash
# 1 Download
git clone https://github.com/malabz/ViRA.git

# 2 Open the folder
cd vira

# 3 Compile
make

# 4 Test vira
./vira -h
```

2.You can also compile manually.

```bash
g++ -std=c++17 -O2 -pthread -Iinclude src/*.cpp -o vira
```

For an optimized build:

```bash
make release
```

### 2 Usage

```text
Usage:
  ./vira -r REF.fasta -q QUERY.fasta -o OUT.fasta [options]

Required arguments:
  -r, --ref FILE                       Uncompressed FASTA (.fa/.fasta/.fas/.fna), exactly one sequence
  -q, --query FILE                     Uncompressed FASTA (.fa/.fasta/.fas/.fna)
  -o, --out FILE                       Final full MSA FASTA output

Parallelism and chunking:
  --max-instances INT                  Concurrent MAFFT processes [16]
  --threads-per-mafft INT              Threads per MAFFT process [1]
  --chunk-size INT                     Maximum records per chunk [1000]
  --target-chunk-bases INT             Maximum total bases per chunk [auto]
  --disable-auto-chunk                 Disable automatic chunk-size reduction
  --max-pending-factor INT             Max queued tasks = max-instances * factor [2]

MAFFT runtime controls:
  --mafft-timeout-sec INT              Per-chunk timeout; 0 disables timeout [0]
  --slow-chunk-sec FLOAT               Report chunks taking at least this many seconds [300]

Temporary files:
  --temp-dir DIR                       Base temp dir; random numeric suffix is appended
  --cleanup                            Remove the temporary directory after completion

Layer 1: input QC:
  --duplicate-id-policy MODE           exclude|rename [exclude]
  --allowed-seq-chars STRING           Allowed IUPAC-like sequence symbols [ACGTURYSWKMBDHVN-.]
  --review-n-ratio FLOAT               Mark for review if N ratio >= value [0.05]
  --exclude-n-ratio FLOAT              Exclude if N ratio >= value [0.50]
  --exclude-min-len-ratio FLOAT        Exclude if seq/ref length ratio < value [0.10]
  --exclude-max-len-ratio FLOAT        Exclude if seq/ref length ratio > value [5.00]

Layer 2: virology suitability:
  --reasonable-min-len-ratio FLOAT     Review if seq/ref length ratio < value [0.85]
  --reasonable-max-len-ratio FLOAT     Review if seq/ref length ratio > value [1.15]
  --non-target-keywords CSV            Header keywords suggesting non-target species
  --engineered-keywords CSV            Header keywords suggesting engineered constructs
  --align-review-samples               Also align suspected non-target/engineered records

MAFFT --mapout post-alignment QC:
  --enable-qc                          Enable QC based on MAFFT --mapout
  --warn-lost-ratio FLOAT              WARN if lost_ratio >= value [0.001]
  --fail-lost-ratio FLOAT              FAIL if lost_ratio >= value [0.01]
  --warn-max-block-len INT             WARN if max_lost_block_len >= value [20]
  --fail-max-block-len INT             FAIL if max_lost_block_len >= value [100]
  --save-failed-stderr                 Save stderr for failed MAFFT chunks

Examples:
  ./vira -r ref.fasta -q query.fasta -o result.fasta
  ./vira -r ref.fasta -q query.fasta -o result.fasta --enable-qc --mafft-timeout-sec 1800
```

## 🔬Output files

Auxiliary files share the same prefix as `--out`.

Output|Description
:---:|:---
`PREFIX.qc_summary.tsv`|Summary of input QC and sample classification
`PREFIX.excluded_samples.tsv`|Samples excluded before alignment
`PREFIX.unaligned_samples.tsv`|Samples not dispatched for alignment
`PREFIX.abnormal_samples.tsv`|Samples marked as abnormal or requiring review
`PREFIX.failed_chunks.tsv`|Failed or timed-out MAFFT chunks
`PREFIX.slow_chunks.tsv`|Chunks slower than the configured threshold

When `--enable-qc` is used, the post-alignment MAFFT `--mapout` outputs are also written.

Output|Description
:---:|:---
`PREFIX.filtered.fasta`|Filtered final MSA
`PREFIX.deleted_blocks.tsv`|Deleted block information parsed from MAFFT mapout
`PREFIX.sample_qc_summary.tsv`|Per-sample post-alignment QC summary
`PREFIX.abnormal_deleted_blocks.tsv`|Abnormal deleted blocks detected by post-alignment QC

## 📍Reminder

1. MAFFT installation is required for the utilization of vira.
2. The reference FASTA file must contain exactly one sequence.
3. The query FASTA file should contain unique sequence IDs. Duplicate IDs are excluded by default and can be renamed with `--duplicate-id-policy rename`.
4. Suspected non-target species or engineered constructs are not aligned by default. Use `--align-review-samples` only when you explicitly want these records aligned.
5. The previous `near/short/long × highN/lowN` bucketing logic has been removed. All QC-passing records enter a unified dispatch queue.
6. Reference and query inputs must be uncompressed FASTA files using `.fa`, `.fasta`, `.fas`, or `.fna`. Compressed archives such as `.gz` or `.zip` must be decompressed first.

## 🖥️Environment

System|Compiler
:---:|:---
Linux|GCC with C++17 support
WSL|GCC with C++17 support
macOS|Clang or GCC with C++17 support

## 🙏References & Acknowledgements

We would like to acknowledge the following MSA tool that contributed to the development of vira:

- **[MAFFT](https://mafft.cbrc.jp/alignment/software/)**: This is a widely used multiple sequence alignment tool known for its high accuracy and scalability.

## 🔖Citation

If you use vira in your research, please cite the related publication or repository.

## 👋Contacts

If you find any bug, welcome to contact us on the issues page of this repository.
