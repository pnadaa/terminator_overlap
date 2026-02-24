# BacTermFinder Coverage Analysis

Quantifies strand-specific coverage of bacterial genomic regions (e.g. stem-loops,
insertion sequence targets) by above-threshold terminator windows predicted by
[BacTermFinder](https://github.com/BioinformaticsLabAtMUN/BacTermFinder), and
benchmarks observed coverage against a randomised genomic background.

---

## Overview

Given a multi-FASTA file of query regions (coordinates encoded in each sequence
ID) and a directory tree of per-genome BacTermFinder `*_mean.csv` outputs, this
script:

1. Parses query coordinates and infers strand from raw coordinate order.
2. Loads BacTermFinder windows that meet a user-defined `probability_mean`
   threshold and builds strand-specific interval trees.
3. Computes the overlap in bp and the percentage of each query region covered
   by above-threshold terminator windows (strand-specifically).
4. Optionally places *N* random same-length windows across the genome to
   estimate background coverage (mean ± SD).
5. Writes TSV results, a suite of debug files, and an optional coverage bar chart.

---

## Requirements

### Python ≥ 3.10

| Package | Purpose |
|---|---|
| `numpy` | Vectorised random-control draws, statistics |
| `pandas` | CSV/TSV I/O and tabular processing |
| `biopython` | FASTA parsing (`Bio.SeqIO`) |
| `intervaltree` | Strand-specific interval overlap queries |
| `matplotlib` | Optional coverage plot (`--plot`) |

Install with:

```bash
pip install numpy pandas biopython intervaltree matplotlib
```

## **External tool**

-   **BLAST+** `blastdbcmd` — required only when `--blast-db` is provided for\
    random-control genome-length lookups.

------------------------------------------------------------------------

## **Input formats**

## **Query FASTA (`--query-fasta`)**

Each record ID must follow the pattern:

```text
`{accession}_{raw_start}-{raw_end}`
```


Strand is inferred from coordinate order:

| **Condition**         | **Inferred strand** |
|:----------------------|:--------------------|
| `raw_start ≤ raw_end` | `+` (forward)       |
| `raw_start > raw_end` | `−` (reverse)       |

Examples:

```text
>CP043804_723996-724438    # forward strand
>CP043804_724236-724198    # reverse strand (start > end)`
```

> **Important: If all your FASTA coordinates are written low→high regardless\
> of strand, all records will be inferred as `+`. Reverse-strand entries must\
> genuinely encode `raw_start > raw_end`.**

## **BacTermFinder mean CSV files (`--mean-root`)**

Expected directory layout:

```text
`<mean_root>/
  <ACCESSION>/
    <ACCESSION>_mean.csv`
```
Each CSV must contain at minimum:

| **Column** | **Description** |
|:---|:---|
| `SampleName` | Encodes window coordinates: `window_{genome}_{low}_{high}` or `window_{genome}_{low}_{high}_{strand}` |
| `probability_mean` | Ensemble CNN terminator probability (0–1) |

------------------------------------------------------------------------

## **Usage**

```bash
`python btf_coverage.py \
  --query-fasta    regions.fa \
  --mean-root      /path/to/btf_output/ \
  --threshold      0.3 \
  --blast-db       /path/to/blastdb/nt \
  --random-n       200 \
  --seed           1 \
  --debug-dir      debug/ \
  --out-prefix     results/coverage \
  --plot`
```
## **All arguments**

| **Argument** | **Default** | **Description** |
|:---|:---|:---|
| `--query-fasta` | *(required)* | Multi-FASTA with coordinate-encoded IDs |
| `--mean-root` | *(required)* | Root directory containing per-genome `*_mean.csv` files |
| `--mean-suffix` | `_mean.csv` | Filename suffix for BacTermFinder mean files |
| `--threshold` | `0.3` | Minimum `probability_mean` to include a window |
| `--blast-db` | *None* | BLAST DB prefix for genome-length queries (enables random control) |
| `--blastdbcmd` | `blastdbcmd` | Path to `blastdbcmd` executable |
| `--random-n` | `200` | Random control iterations per query; `0` disables |
| `--seed` | `1` | NumPy RNG seed for reproducibility |
| `--debug-dir` | *(required)* | Directory for debug/troubleshooting TSVs |
| `--debug-max-overlap-windows` | `50` | Max raw overlapping windows listed per query in debug output; `0` disables |
| `--plot` | *off* | Write a PNG bar chart of observed vs random coverage |
| `--out-prefix` | *(required)* | Prefix for main output files |
| `--log-level` | `INFO` | Logging verbosity: `DEBUG`, `INFO`, `WARNING`, `ERROR` |

------------------------------------------------------------------------

## **Outputs**

## **Main results**

**`{out-prefix}.overlap.tsv`**

One row per query region with the following columns:

| **Column** | **Description** |
|:---|:---|
| `qseqid` | FASTA record ID |
| `accession_from_fasta` | Accession parsed from the ID |
| `region_start` / `region_end` | Normalised coordinates (always low→high) |
| `query_strand` | Inferred strand (`+` or `−`) |
| `region_len` | Length of query region (bp) |
| `genome_key_used` | Accession used to locate the mean CSV (may add `.1`) |
| `mean_csv` | Path to the mean CSV that was loaded |
| `windows_rows_kept` | Number of windows passing the threshold |
| `windows_total_bp_this_strand` | Total merged window bp on the query strand |
| `overlap_bp` | Overlap between query region and terminator windows (bp) |
| `percent_region_covered` | `overlap_bp / region_len × 100` |
| `random_mean_percent` | Mean random-control coverage (%) |
| `random_sd_percent` | SD of random-control coverage (%) |
| `blast_entry_used` | BLAST DB entry used for genome-length lookup |
| `blast_entry_len` | Genome length returned by `blastdbcmd` |
| `kept_window_min/max_coord` | Coordinate extent of all kept windows |
| `note` | Warnings or skip reasons |

## **Debug outputs (`{debug-dir}/`)**

| **File** | **Contents** |
|:---|:---|
| `debug.per_query.tsv` | Verbose per-query diagnostics including all paths tried |
| `debug.mean_missing.tsv` | Queries for which no mean CSV was found |
| `debug.length_failed.tsv` | Queries for which `blastdbcmd` length lookup failed |
| `debug.overlapping_windows.tsv` | Raw (pre-merge) windows overlapping each query region, with individual `probability_mean` values |

## **Optional plot**

**`{out-prefix}.overlap.png`** — horizontal bar chart of per-query observed\
coverage with random mean ± SD error bars (produced with `--plot`).

------------------------------------------------------------------------

## **How accession resolution works**

The script attempts to locate mean CSVs using a `.1`-preferring strategy:

-   Input `CP043804` → tries `CP043804.1/`, then `CP043804/`

-   Input `CP043804.1` → tries `CP043804.1/`, then `CP043804/`

This avoids silent mismatches between versioned and unversioned accessions in\
FASTA IDs and directory names.

------------------------------------------------------------------------

## **Random control**

When `--blast-db` and `--random-n > 0` are provided, for each query the script:

1.  Looks up the genome length via `blastdbcmd`.

2.  Draws `--random-n` uniformly random start positions across the genome using\
    a seeded NumPy generator.

3.  Computes terminator-window coverage for each random window of the same\
    length as the query, on the same strand.

4.  Reports mean ± SD as `random_mean_percent` / `random_sd_percent`.

The random control is skipped per query if `--blast-db` is absent, the genome\
length lookup fails, or the genome is shorter than the query region.

------------------------------------------------------------------------

## **Logging and diagnostics**

The script emits `WARNING`-level messages for common data quality issues:

-   Windows with no strand information in `SampleName` (inserted into both strand\
    trees; strand-specific queries become strand-naive).

-   Windows whose embedded genome ID does not match the expected accession after\
    version-stripping.

-   Queries where all FASTA records infer to the same strand (possible coordinate\
    encoding issue).

-   Genomes where no windows pass the threshold.

Run with `--log-level DEBUG` to see cache hits and per-genome window counts.

------------------------------------------------------------------------

## **Notes**

-   Overlapping windows are merged before coverage calculation to avoid\
    double-counting.

-   All interval arithmetic uses half-open intervals internally (`[begin, end)`);\
    FASTA coordinates are treated as fully inclusive (`[start, end]`).

-   The script caches loaded interval trees per genome, so genomes shared across\
    multiple query entries are only loaded once.