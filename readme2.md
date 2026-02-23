```         
# BacTermFinder Coverage Analysis

Quantifies the coverage of bacterial stem-loop regions by
[BacTermFinder](https://github.com/BioinformaticsLabAtMUN/BacTermFinder)
terminator prediction windows and statistically compares observed coverage
against a random genomic background.

BacTermFinder is a CNN-based tool that predicts both intrinsic and
factor-dependent bacterial transcription terminators by sliding a fixed-length
window across a genome and returning a `probability_mean` score per window.
This script asks: *do your stem-loop regions of interest fall preferentially
within high-confidence terminator windows, relative to chance?*

---

## Features

- Parses a multi-FASTA file of stem-loop sequences whose IDs encode genomic
  coordinates (`accession_start_end`)
- Loads and caches BacTermFinder `_mean.csv` result files per accession,
  with automatic accession version normalisation (e.g. `.1` suffixes)
- Computes per-stem-loop **coverage percentage** by windows above a
  user-defined probability threshold, with strand-aware filtering
- Builds an `IntervalTree` index automatically for large window sets
  (&gt;500 k rows) for efficient overlap queries
- Performs a **random control** analysis by sampling *n* genomic positions of
  the same length as each stem-loop and computing mean coverage, standard
  deviation, Z-score, and empirical p-value
- Outputs a tidy **CSV** of all metrics and a three-panel **PNG figure**

---

## Installation

```bash
pip install pandas numpy matplotlib seaborn intervaltree biopython
```

Or with conda/micromamba:

```         
bash
```

`conda install -c conda-forge pandas numpy matplotlib seaborn intervaltree biopython`

Python ≥ 3.9 is recommended.

------------------------------------------------------------------------

## **Input Requirements**

## **1. Stem-loop FASTA (`-f`)**

A multi-FASTA file where each record ID follows the pattern:

```         
text
```

`<accession>_<start>-<end>`

For example:

```         
text
```

`>NC_002952_14231-14289
GCTAGCTAGCTAGCT...`

Coordinates are 1-based and inclusive. Accession version suffixes (e.g.\
`NC_002952.1`) are handled automatically.

## **2. BacTermFinder results directory (`-r`)**

A directory containing per-genome subdirectories of BacTermFinder output.\
The script searches for `_mean.csv` files in the following paths (tried in\
order):

```         
text
```

`<results_dir>/<acc>/<acc>_mean.csv
<results_dir>/<acc>.1/<acc>.1_mean.csv
<results_dir>/<norm_acc>/<norm_acc>_mean.csv`

Each `_mean.csv` must contain at least:

| **Column**         | **Description**                                 |
|:-------------------|:------------------------------------------------|
| `SampleName`       | `<label>, <acc>_<win_start>_<win_end>_<strand>` |
| `probability_mean` | Mean CNN probability score for the window       |

------------------------------------------------------------------------

## **Usage**

```         
bash
```

`python btf_coverage.py \
  -f stem_loops.fasta \
  -r /path/to/btf_results/ \
  -t 0.3 \
  -n 1000 \
  --strand both \
  --output-csv terminator_coverage.csv \
  --output-plot terminator_coverage.png \
  --seed 42`

## **Arguments**

| **Flag** | **Default** | **Description** |
|:---|:---|:---|
| `-f`, `--fasta` | *(required)* | Multi-FASTA of stem-loop sequences |
| `-r`, `--results-dir` | *(required)* | BacTermFinder results directory |
| `-t`, `--threshold` | `0.3` | Minimum `probability_mean` for a window to count |
| `-n`, `--n-controls` | `1000` | Random genomic positions sampled per stem-loop |
| `--strand` | `both` | Windows to include: `both`, `plus`, or `minus` |
| `--output-csv` | `terminator_coverage.csv` | Output CSV path |
| `--output-plot` | `terminator_coverage.png` | Output PNG path |
| `--seed` | `42` | NumPy random seed for reproducibility |
| `-v`, `--verbose` | `False` | Enable DEBUG-level logging |

------------------------------------------------------------------------

## **Output**

## **CSV (`terminator_coverage.csv`)**

One row per input stem-loop:

| **Column** | **Description** |
|:---|:---|
| `seq_id` | Original FASTA record ID |
| `accession` | Parsed genome accession |
| `start` / `end` | Genomic coordinates (1-based, inclusive) |
| `length` | Stem-loop length in bp |
| `n_windows_above_threshold` | Count of overlapping windows ≥ threshold |
| `coverage_pct` | \% of stem-loop bases covered by merged windows ≥ threshold |
| `max_probability_mean` | Highest `probability_mean` among overlapping windows |
| `mean_probability_mean` | Mean `probability_mean` among overlapping windows |
| `control_mean_coverage_pct` | Mean coverage across *n* random genomic positions |
| `control_std_coverage_pct` | Std dev of random control coverages |
| `z_score` | `(observed − control_mean) / control_std` |
| `empirical_pvalue` | Fraction of random controls with coverage ≥ observed |

`NaN` is reported for any stem-loop whose accession CSV was not found, or\
where the genome is shorter than the stem-loop.

## **Figure (`terminator_coverage.png`)**

A 300 dpi three-panel figure:

-   **Panel A — Observed vs. Random Coverage**: paired boxplot + strip plot\
    comparing each stem-loop's observed coverage to its random control mean,\
    with connecting lines and median annotations.

-   **Panel B — Coverage vs. Terminator Probability**: scatter plot of\
    `max_probability_mean` against `coverage_pct`, coloured by whether the\
    stem-loop is covered at all; the threshold is marked with a dashed line.

-   **Panel C — Empirical P-value Distribution**: histogram of empirical\
    p-values across all stem-loops; p = 0.05 is marked with a dashed line.

------------------------------------------------------------------------

## **How the Statistics Work**

**Coverage** is computed by:

1.  Collecting all BacTermFinder windows that overlap the query interval

2.  Optionally filtering by strand and by `probability_mean ≥ threshold`

3.  Clipping each surviving window to the query boundaries

4.  Merging overlapping clipped intervals

5.  Dividing total covered bases by stem-loop length × 100

**Random control** draws `--n-controls` random start positions uniformly from\
`[1, genome_length − stem_length + 1]` (genome length is estimated as the\
maximum window end coordinate in the CSV). The same coverage calculation is\
applied to each random position. The **empirical p-value** is the proportion\
of random positions whose coverage ≥ observed coverage — a value near 0\
indicates that the observed stem-loop is covered far more than expected by\
chance.

------------------------------------------------------------------------

## **Notes**

-   The script caches BacTermFinder windows per unique accession, so\
    multiple stem-loops from the same genome only trigger one file read.

-   For very large CSVs (\>500 k windows), an `IntervalTree` is built\
    automatically to speed up overlap queries.

-   If a BacTermFinder CSV cannot be resolved for an accession, all metrics\
    for the affected stem-loops are set to `NaN` and a warning is logged.

-   The random control genome length is inferred from `max(win_end)` in the CSV.\
    If this is a poor estimate for your genome, consider pre-filtering your\
    BacTermFinder CSVs to full-genome runs only.