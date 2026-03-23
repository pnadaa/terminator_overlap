**BacTermFinder Coverage Analysis Script**

Quantifies the overlap of bacterial genomic regions of interest against
[BacTermFinder](https://github.com/BioinformaticsLabAtMUN/BacTermFinder)
terminator window predictions, and benchmarks observed overlap against a
random genomic background using an empirical permutation test.

---

## Overview

Given a multi-FASTA file of query regions (coordinates encoded in the
sequence IDs) and a directory of per-genome BacTermFinder `_mean.csv`
output files, this script:

1. Parses each FASTA record to extract genomic coordinates and strand.
2. Loads per-genome terminator windows from `_mean.csv` into strand-aware
   `IntervalTree` structures, filtered by a `probability_mean` threshold.
3. Computes strand-specific overlap (bp and %) between each query region
   and predicted terminator windows.
4. Extracts per-query probability metrics (max, mean, median, and
   best-qualifying-window `probability_mean`).
5. Runs a randomisation control (default: 1000 iterations) by placing
   the query region at random genomic positions and computing overlap,
   yielding empirical p-values and z-scores.
6. Writes TSV results and optional seaborn visualisation plots.

---

## Requirements

### Python >= 3.9

### Python packages
```
numpy
pandas
biopython
intervaltree
matplotlib
seaborn
```

Install via pip:
```bash
pip install numpy pandas biopython intervaltree matplotlib seaborn
```

### External tools (optional, required for random control)
- `blastdbcmd` (NCBI BLAST+) — used to retrieve genome lengths from a
  local BLAST database for the randomisation test. If not provided, the
  random control is disabled.

---

## Input Files

### `--query-fasta` (required)
A multi-FASTA file where each record ID encodes the genomic accession,
start, and end coordinates:

```
>{accession}_{start}-{end}
```

**Strand is inferred from coordinate order:**
- `start <= end` → forward strand (`+`)
- `start > end` → reverse strand (`−`), coordinates are then normalised

**Examples:**
```
>CP043804_723996-724438       # forward strand
>CP043804_724438-723996       # reverse strand
```

The FASTA sequence itself is not used; only the ID is parsed.

### `--mean-root` (required)
Root directory containing per-genome BacTermFinder output subdirectories:

```
<mean-root>/
  <GENOME_ACCESSION>/
    <GENOME_ACCESSION>_mean.csv
```

Each `_mean.csv` file must contain at least two columns:
- `SampleName` — window identifier with encoded coordinates and strand
  in the format `window_{genome}_{low}_{high}` or
  `window_{genome}_{low}_{high}_{strand}`
- `probability_mean` — BacTermFinder terminator probability score

---

## Usage

```bash
python terminator_overlap.py \
  --query-fasta  regions.fasta \
  --mean-root    /path/to/bactermfinder_output/ \
  --out-prefix   results/my_analysis \
  --debug-dir    results/debug/ \
  [OPTIONS]
```

---

## Arguments

### Required
| Argument         | Description |
|------------------|-------------|
| `--query-fasta`  | Multi-FASTA file with coordinate-encoded sequence IDs |
| `--mean-root`    | Root directory for per-genome `_mean.csv` files |
| `--out-prefix`   | Output file prefix (directory will be created if needed) |
| `--debug-dir`    | Directory for debug/troubleshooting TSV outputs |

### Thresholds
| Argument                         | Default | Description |
|----------------------------------|---------|-------------|
| `--coverage-threshold`           | `0.3`   | `probability_mean` cutoff for coverage/overlap_bp calculation |
| `--probability-threshold`        | `0.0`   | `probability_mean` cutoff for raw probability metric queries |
| `--plot-min-window-overlap-frac` | `0.99`  | Minimum fraction of the query region a window must overlap to qualify for best-window scoring |

### Random Control
| Argument        | Default       | Description |
|-----------------|---------------|-------------|
| `--random-n`    | `1000`        | Number of random placements per query (set `0` to disable) |
| `--seed`        | `1`           | NumPy RNG seed for reproducibility |
| `--blast-db`    | *(none)*      | Local BLAST DB prefix for genome length lookups |
| `--blastdbcmd`  | `blastdbcmd`  | Path to `blastdbcmd` executable |

### Performance
| Argument     | Default | Description |
|--------------|---------|-------------|
| `--workers`  | `1`     | Worker processes for genome-level parallelism; `0` = use all available CPUs (respects `PBS_NCPUS` on HPC) |

### Input Format
| Argument        | Default      | Description |
|-----------------|--------------|-------------|
| `--mean-suffix` | `_mean.csv`  | Filename suffix for BacTermFinder mean CSV files |

### Output / Plotting
| Argument                        | Description |
|---------------------------------|-------------|
| `--plot`                        | Write the new 3-panel probability-focused seaborn figure |
| `--legacy-plot`                 | Write the previous coverage/probability legacy figure |
| `--debug-max-overlap-windows`   | Max overlapping windows to list per query in debug TSV (default: `50`; `0` disables) |
| `--log-level`                   | Logging verbosity: `DEBUG`, `INFO` (default), `WARNING`, `ERROR` |

---

## Output Files

### Main results
| File | Description |
|------|-------------|
| `{out_prefix}.overlap.tsv` | Per-query overlap statistics and probability metrics |

### Plots (optional)
| File | Description |
|------|-------------|
| `{out_prefix}.overlap.png` / `.svg` | New 3-panel probability plot (`--plot`) |
| `{out_prefix}.legacy.overlap.png` / `.svg` | Legacy coverage plot (`--legacy-plot`) |

### Debug outputs (always written to `--debug-dir`)
| File | Description |
|------|-------------|
| `debug.per_query.tsv` | Detailed per-query diagnostic fields |
| `debug.mean_missing.tsv` | Queries where no `_mean.csv` was found |
| `debug.length_failed.tsv` | Queries where BLAST genome length lookup failed |
| `debug.overlapping_windows.tsv` | Individual overlapping windows per query (up to `--debug-max-overlap-windows`) |

### Output TSV columns (`*.overlap.tsv`)

| Column | Description |
|--------|-------------|
| `qseqid` | Query FASTA record ID |
| `accession_from_fasta` | Genome accession parsed from FASTA ID |
| `region_start` / `region_end` | Normalised genomic coordinates (1-based inclusive) |
| `query_strand` | Inferred strand (`+` or `−`) |
| `region_len` | Query region length in bp |
| `genome_key_used` | Resolved genome accession used to locate `_mean.csv` |
| `windows_rows_kept` | Number of windows above `--coverage-threshold` |
| `windows_total_bp_this_strand` | Total bp of merged windows on the query strand |
| `n_windows_above_threshold` | Number of windows overlapping the query region |
| `max_probability_mean` | Maximum `probability_mean` of overlapping windows |
| `mean_probability_mean` | Mean `probability_mean` of overlapping windows |
| `median_probability_mean` | Median `probability_mean` of overlapping windows |
| `best_overlap_probability_mean` | Highest `probability_mean` among qualifying windows (overlap > `--plot-min-window-overlap-frac`) |
| `best_overlap_window_*` | Metadata for the best-qualifying window (SampleName, coords, overlap bp/frac) |
| `overlap_bp` | Total bp of overlap between query region and merged windows |
| `percent_region_covered` | Percentage of the query region covered by windows |
| `random_mean_percent` | Mean coverage (%) across random placements |
| `random_sd_percent` | SD of coverage (%) across random placements |
| `random_best_probability_mean_mean` | Mean best-window probability across random placements |
| `random_best_probability_mean_median` | Median best-window probability across random placements |
| `random_best_probability_mean_sd` | SD of best-window probability across random placements |
| `z_score` | Z-score of observed best-window probability vs. random distribution |
| `empirical_pvalue_best_overlap` | Empirical p-value: fraction of random placements with best-window probability >= observed |
| `blast_entry_used` | BLAST DB entry used for genome length lookup |
| `blast_entry_len` | Genome length (bp) from BLAST DB |
| `note` | Warnings or skip reasons for this query |

---

## Visualisations

### New probability-focused plot (`--plot`)
A 3-panel seaborn figure:
1. **Panel A** — Observed vs. random best-qualifying window `probability_mean`
   (boxplot + stripplot with per-query connecting lines)
2. **Panel B** — Distribution of per-query overlapping-window probability
   summaries (violin or boxplot depending on n)
3. **Panel C** — Ranked empirical p-values (−log₁₀ scale) from the random
   best-window null distribution, with significance thresholds marked

### Legacy coverage plot (`--legacy-plot`)
A 3-panel figure showing:
1. Observed vs. random mean coverage (%)
2. Coverage (%) vs. max `probability_mean` scatter plot
3. Empirical p-value distribution histogram + KDE

---

## Examples

### Basic run (no random control)
```bash
python terminator_overlap.py \
  --query-fasta   IS_flanking_regions.fasta \
  --mean-root     /data/bactermfinder/ \
  --out-prefix    results/IS_terminators \
  --debug-dir     results/debug \
  --random-n      0
```

### Full run with random control and plots
```bash
python terminator_overlap.py \
  --query-fasta        IS_flanking_regions.fasta \
  --mean-root          /data/bactermfinder/ \
  --blast-db           /databases/nt_bacteria \
  --random-n           1000 \
  --coverage-threshold 0.3 \
  --workers            8 \
  --out-prefix         results/IS_terminators \
  --debug-dir          results/debug \
  --plot \
  --legacy-plot \
  --log-level          INFO
```

### HPC (PBS/Torque) — auto CPU detection
```bash
#PBS -l ncpus=16
python terminator_overlap.py \
  --query-fasta   regions.fasta \
  --mean-root     /scratch/bactermfinder/ \
  --blast-db      /scratch/blast/nt_bacteria \
  --workers       0 \
  --out-prefix    $PBS_O_WORKDIR/results/analysis \
  --debug-dir     $PBS_O_WORKDIR/results/debug
```

---

## Notes

- **Strand specificity**: Strand is inferred solely from the raw coordinate
  order in the FASTA ID. Reverse-strand queries must encode `raw_start > raw_end`.
- **Genome accession versioning**: The script automatically tries `.1`-versioned
  accession variants (e.g., `CP043804.1`) when resolving `_mean.csv` paths and
  BLAST entries.
- **Parallelism**: Work is parallelised at the genome level (one genome per
  worker process). Each worker independently loads interval trees for its
  assigned genome, avoiding memory duplication. On Linux, `forkserver` context
  is used for clean worker processes.
- **Randomisation**: The empirical p-value uses the formula
  `(count_random_>=_observed + 1) / (n_random + 1)` to avoid p = 0 artefacts.
  Minimum detectable p-value is `1 / (random_n + 1)`.
