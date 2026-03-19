# BacTermFinder Coverage Analysis

A strand-aware Python script that quantifies coverage of bacterial genomic regions by
[BacTermFinder](https://github.com/BioinformaticsLabAtMUN/BacTermFinder) terminator
windows, and compares observed regions against a random genomic background to assess
statistical enrichment.

---

## Overview

Given a set of query genomic regions (encoded as FASTA coordinate IDs) and per-genome
BacTermFinder output (`*_mean.csv`), this script:

1. Loads above-threshold BacTermFinder windows into strand-aware interval trees
2. Computes overlap and probability metrics between each query region and overlapping windows
3. Generates a random genomic background (Monte Carlo) per query to derive empirical p-values
4. Writes TSV results, debug outputs, and optional publication-quality plots

---

## Requirements

- Python ≥ 3.9
- [Biopython](https://biopython.org/) (`Bio.SeqIO`)
- [intervaltree](https://github.com/chaimleib/intervaltree)
- [NumPy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [matplotlib](https://matplotlib.org/)
- [seaborn](https://seaborn.pydata.org/)
- NCBI `blastdbcmd` (only required for `--blast-db` random control)

Install Python dependencies:

```bash
pip install biopython intervaltree numpy pandas matplotlib seaborn
```

---

## Input Formats

### Query FASTA (`--query-fasta`)

A multi-FASTA file where each record ID encodes the accession and genomic coordinates:

```
>{accession}_{start}-{end}
```

Strand is inferred from coordinate order:
- `raw_start <= raw_end` → forward strand (`+`)
- `raw_start > raw_end`  → reverse strand (`-`), coordinates are then normalised to `start < end`

**Examples:**
```
>CP043804_723996-724438     # forward strand
>CP043804_724236-724198     # reverse strand (raw_start > raw_end)
```

### Mean CSV root (`--mean-root`)

A directory tree where each genome's BacTermFinder mean output is located at:

```
<mean-root>/<GENOME>/<GENOME>_mean.csv
```

Each `*_mean.csv` must contain at minimum two columns:
- `SampleName` — encodes window coordinates and optionally strand, e.g.  
  `window_CP043804.1_10000_10200_+` or `window_CP043804.1_10000_10200`
- `probability_mean` — BacTermFinder mean terminator probability for that window

Both plain `.csv` and gzip-compressed `.csv.gz` files are supported.

---

## Usage

```bash
python bactermfinder_coverage_fixed.py \
    --query-fasta    queries.fasta \
    --mean-root      /path/to/mean_csvs/ \
    --out-prefix     results/my_run \
    --debug-dir      results/debug/ \
    --blast-db       /path/to/blastdb/nt \
    --threshold      0.3 \
    --random-n       1000 \
    --plot \
    --legacy-plot
```

### All arguments

| Argument | Default | Description |
|---|---|---|
| `--query-fasta` | *(required)* | Multi-FASTA with coordinate-encoded IDs |
| `--mean-root` | *(required)* | Root directory for per-genome `*_mean.csv` files |
| `--out-prefix` | *(required)* | Output path prefix for TSV and PNG/SVG files |
| `--debug-dir` | *(required)* | Directory for troubleshooting TSV outputs |
| `--mean-suffix` | `_mean.csv` | Filename suffix for mean CSV files |
| `--threshold` | `0.3` | `probability_mean` threshold for window inclusion |
| `--blast-db` | `None` | BLAST DB prefix for genome length lookup (enables random control) |
| `--blastdbcmd` | `blastdbcmd` | Path to `blastdbcmd` executable |
| `--random-n` | `1000` | Monte Carlo iterations per query (set `0` to disable) |
| `--seed` | `1` | NumPy RNG seed for reproducibility |
| `--plot` | flag | Write the new 3-panel probability-focused plot |
| `--legacy-plot` | flag | Write the legacy 3-panel coverage/probability plot |
| `--plot-min-window-overlap-frac` | `0.10` | Minimum overlap fraction for a window to qualify for best-window scoring |
| `--debug-max-overlap-windows` | `50` | Max overlapping windows listed per query in debug TSV (`0` disables) |
| `--log-level` | `INFO` | Logging verbosity: `DEBUG`, `INFO`, `WARNING`, `ERROR` |

---

## Outputs

### Main results TSV

**`<out-prefix>.overlap.tsv`** — one row per query region, with columns:

| Column | Description |
|---|---|
| `qseqid` | Query FASTA record ID |
| `accession_from_fasta` | Genome accession parsed from the FASTA ID |
| `region_start` / `region_end` | Normalised (start ≤ end) genomic coordinates |
| `query_strand` | Inferred strand (`+` or `-`) |
| `region_len` | Query region length in bp |
| `genome_key_used` | Genome key used to locate the `*_mean.csv` file |
| `n_windows_above_threshold` | Count of overlapping windows above `--threshold` |
| `max_probability_mean` | Maximum `probability_mean` of all overlapping windows |
| `mean_probability_mean` | Mean `probability_mean` across overlapping windows |
| `median_probability_mean` | Median `probability_mean` across overlapping windows |
| `best_overlap_probability_mean` | Highest `probability_mean` among windows passing `--plot-min-window-overlap-frac` |
| `best_overlap_window_sample` | `SampleName` of that best-qualifying window |
| `overlap_bp` | Total bp of merged window coverage over the query region |
| `percent_region_covered` | `overlap_bp / region_len × 100` |
| `random_mean_percent` | Mean coverage (%) across random control iterations |
| `random_sd_percent` | SD of random coverage (%) |
| `random_best_probability_mean_mean` | Mean best-window probability across random iterations |
| `random_best_probability_mean_sd` | SD of random best-window probability |
| `z_score` | Z-score of observed best-window probability vs random distribution |
| `empirical_pvalue_best_overlap` | Empirical p-value: fraction of random iterations ≥ observed best probability |
| `blast_entry_len` | Genome length retrieved from BLAST DB (for random control) |
| `note` | Any warnings or skip reasons for this query |

### Debug TSVs (`--debug-dir`)

| File | Description |
|---|---|
| `debug.per_query.tsv` | Full per-query summary including paths tried and all metrics |
| `debug.mean_missing.tsv` | Queries where no `*_mean.csv` could be found |
| `debug.length_failed.tsv` | Queries where `blastdbcmd` genome length lookup failed |
| `debug.overlapping_windows.tsv` | Per-window overlap details for each query (up to `--debug-max-overlap-windows` per query) |

### Plots

Both plot functions save **PNG (300 dpi)** and **SVG** (for vector editing):

#### New probability-focused plot (`--plot`)

`<out-prefix>.overlap.png` / `.overlap.svg` — 3-panel figure:

- **Panel A** — Paired boxplot+strip of observed vs random best-window `probability_mean`,
  with per-query connecting lines, mean/median annotations, and threshold line
- **Panel B** — Violin (or boxplot for small n) of four per-query probability summary metrics:
  mean, median, max overlap probability, and best qualifying window probability
- **Panel C** — Ranked -log₁₀(empirical p-value) scatter plot; significant hits (p ≤ 0.05)
  coloured red and annotated with query IDs; p = 0.05 and p = 0.10 threshold lines shown

#### Legacy coverage plot (`--legacy-plot`)

`<out-prefix>.legacy.overlap.png` / `.legacy.overlap.svg` — 3-panel figure:

- **Panel A** — Paired boxplot+strip of observed vs random mean coverage (%)
- **Panel B** — Scatter of max `probability_mean` vs coverage (%), coloured by a continuous
  viridis scale reflecting coverage magnitude
- **Panel C** — Histogram + KDE of empirical p-value distribution

---

## Strand Handling

Strand is inferred purely from the FASTA coordinate order — no strand field is required in
the FASTA header. Windows in `*_mean.csv` are placed into strand-specific interval trees
only if their `SampleName` encodes a strand suffix (`+` or `-`). Unstranded windows are
inserted into both trees. If all windows lack strand information, the script will log a
warning and strand queries will behave as strand-naive.

---

## Random Control Details

For each query, `--random-n` random regions of the same length are drawn uniformly across
the genome (bounded by the genome length from `blastdbcmd`). For each random placement:
- Coverage overlap is computed against the merged interval trees
- The best-qualifying window probability is computed identically to the observed case

The empirical p-value uses the Laplace-smoothed estimator:

```
p = (count(random_best_prob >= observed_best_prob) + 1) / (random_n + 1)
```

This avoids p = 0 and is conservative. Set `--random-n 0` or omit `--blast-db` to skip
the random control entirely (p-values will be `NaN`).

---

## Example

```bash
python bactermfinder_coverage_fixed.py \
    --query-fasta    IS_flanking_regions.fasta \
    --mean-root      /data/bactermfinder_output/ \
    --out-prefix     output/IS_terminator_overlap \
    --debug-dir      output/debug/ \
    --blast-db       /databases/blast/nt \
    --threshold      0.5 \
    --random-n       2000 \
    --seed           42 \
    --plot \
    --legacy-plot \
    --plot-min-window-overlap-frac 0.20 \
    --log-level      INFO
```

This will produce:
```
output/
  IS_terminator_overlap.overlap.tsv
  IS_terminator_overlap.overlap.png
  IS_terminator_overlap.overlap.svg
  IS_terminator_overlap.legacy.overlap.png
  IS_terminator_overlap.legacy.overlap.svg
  debug/
    debug.per_query.tsv
    debug.mean_missing.tsv
    debug.length_failed.tsv
    debug.overlapping_windows.tsv
```

---

## Notes

- Accession version suffixes (e.g. `.2`) are stripped when matching FASTA IDs to `*_mean.csv`
  filenames; the script always tries `ACCESSION.1` first, then the bare accession
- Multiple genomes are cached so each `*_mean.csv` is loaded only once per run
- For very large FASTA inputs, consider increasing `--random-n` only after profiling runtime,
  as it scales linearly with the number of queries
