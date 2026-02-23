"""
BacTermFinder Coverage Analysis Script
Quantifies coverage of bacterial stem-loop regions by BacTermFinder terminator windows
and performs statistical comparison against a random genomic background.
"""

import os
import re
import logging
import argparse
import warnings
from pathlib import Path
from typing import Tuple, List, Dict, Any, Optional

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from intervaltree import IntervalTree
from Bio import SeqIO


def normalise_accession(acc: str) -> str:
    """Strip a trailing version suffix of the form '.<integer>' from acc."""
    return re.sub(r'\.\d+$', '', acc)


def parse_fasta_id(record_id: str) -> Tuple[str, int, int]:
    """Parse FASTA ID string to extract accession, start, and end coordinates."""
    # Updated to match {accessionID}_{start}-{end}
    pattern = re.compile(r'^(.+)_(\d+)-(\d+)$')
    match = pattern.match(record_id)
    if not match:
        raise ValueError(f"Regex does not match for record ID: {record_id}")
    return match.group(1), int(match.group(2)), int(match.group(3))


def resolve_csv_path(results_dir: Path, acc: str) -> Optional[Path]:
    """Find the corresponding CSV for an accession among valid candidate paths."""
    norm_acc = normalise_accession(acc)
    paths = [
        results_dir / acc / f"{acc}_mean.csv",
        results_dir / f"{acc}.1" / f"{acc}.1_mean.csv",
        results_dir / norm_acc / f"{norm_acc}_mean.csv"
    ]
    
    # Use a set to handle cases where different candidates map to the exact same file string
    existing_paths = {p for p in paths if p.exists() and p.is_file()}
    
    if len(existing_paths) > 1:
        raise RuntimeError(f"Conflicting paths found for accession {acc}: {existing_paths}")
        
    if not existing_paths:
        return None
        
    return list(existing_paths)[0]


def load_bactermfinder_windows(csv_path: Path) -> Tuple[pd.DataFrame, int]:
    """
    Load BacTermFinder results from CSV, returning a validated DataFrame and estimated genome length.
    Automatically builds an IntervalTree cache for > 500k rows.
    """
    df = pd.read_csv(csv_path, sep=None, engine='python')
    if 'SampleName' not in df.columns or 'probability_mean' not in df.columns:
        logging.error(f"Missing required columns in {csv_path}")
        return pd.DataFrame(columns=['win_start', 'win_end', 'strand', 'probability_mean']), 0

    valid_rows = []
    # Updated regex to match format: windowID_accessionID_startcoordinate_endcoordinate_strand
    # [^_]+ matches the windowID up to the first underscore.
    # (.+) greedily matches the accessionID (which may contain underscores).
    pattern = re.compile(r'^[^_]+_(.+)_(\d+)_(\d+)_([+-])$')

    for idx, row in df.iterrows():
        sample_name = str(row['SampleName'])

        match = pattern.match(sample_name)
        if not match:
            logging.warning(f"Row {idx} in {csv_path.name}: Regex failed on SampleName '{sample_name}'")
            continue

        win_accession, win_start, win_end, strand = match.groups()
        valid_rows.append({
            'win_start': int(win_start),
            'win_end': int(win_end),
            'strand': strand,
            'probability_mean': float(row['probability_mean'])
        })

    out_df = pd.DataFrame(valid_rows)
    max_end = out_df['win_end'].max() if not out_df.empty else 0

    if len(out_df) > 500_000:
        logging.info(f"Window data for {csv_path.name} exceeds 500,000 rows. Building IntervalTree...")
        tree = IntervalTree()
        for row in out_df.itertuples(index=False):
            tree[row.win_start : row.win_end + 1] = {
                'win_start': row.win_start,
                'win_end': row.win_end,
                'strand': row.strand,
                'probability_mean': row.probability_mean
            }
        out_df.attrs['interval_tree'] = tree

    return out_df, max_end


def compute_coverage(
    query_start: int,
    query_end: int,
    windows_df: pd.DataFrame,
    threshold: float,
    strand_filter: str,
) -> Tuple[float, int, float, float]:
    """
    Calculates coverage of a query interval by high-probability terminator windows.
    Returns: (coverage_pct, n_windows_above_thresh, max_probability_mean, mean_probability_mean)
    """
    # a. Filter overlaps
    tree = windows_df.attrs.get('interval_tree')
    cols = ['win_start', 'win_end', 'strand', 'probability_mean']
    
    if tree is not None:
        overlapping = tree.overlap(query_start, query_end + 1)
        if not overlapping:
            overlap_df = pd.DataFrame(columns=cols)
        else:
            overlap_df = pd.DataFrame([iv.data for iv in overlapping])
    else:
        overlap_df = windows_df[
            (windows_df['win_start'] <= query_end) & 
            (windows_df['win_end'] >= query_start)
        ]

    # b. Strand filter
    if strand_filter == "plus":
        overlap_df = overlap_df[overlap_df['strand'] == "+"]
    elif strand_filter == "minus":
        overlap_df = overlap_df[overlap_df['strand'] == "-"]

    # c. Filter threshold
    survivors = overlap_df[overlap_df['probability_mean'] >= threshold]

    # d. If no windows remain
    if survivors.empty:
        return (0.0, 0, float('nan'), float('nan'))

    # e. Clip each surviving window
    clipped_intervals = []
    for _, row in survivors.iterrows():
        c_start = max(row['win_start'], query_start)
        c_end = min(row['win_end'], query_end)
        if c_start <= c_end:
            clipped_intervals.append([c_start, c_end])

    if not clipped_intervals:
        return (0.0, 0, float('nan'), float('nan'))

    # f. Merge overlapping clipped intervals
    clipped_intervals.sort(key=lambda x: x[0])
    merged = []
    for start, end in clipped_intervals:
        if not merged:
            merged.append([start, end])
        else:
            last_start, last_end = merged[-1]
            if start <= last_end:
                merged[-1][1] = max(last_end, end)
            else:
                merged.append([start, end])

    # g. Calculate coverage
    coverage_bp = sum(end - start + 1 for start, end in merged)
    coverage_pct = 100.0 * coverage_bp / (query_end - query_start + 1)

    # h. Calculate probability metrics
    max_prob = float(survivors['probability_mean'].max())
    mean_prob = float(survivors['probability_mean'].mean())
    n_windows = len(survivors)

    # i. Return
    return (coverage_pct, n_windows, max_prob, mean_prob)


def compute_control(
    stem_length: int,
    genome_length: int,
    windows_df: pd.DataFrame,
    threshold: float,
    strand_filter: str,
    n_controls: int,
    rng: np.random.Generator,
    observed_coverage_pct: float
) -> Tuple[float, float, float, float]:
    """
    Computes the random control coverage metrics against a null genomic background.
    Returns: (control_mean_coverage_pct, control_std_coverage_pct, z_score, empirical_pvalue)
    """
    # a. Check boundaries
    if genome_length < stem_length:
        logging.warning(f"Genome length ({genome_length}) < stem length ({stem_length}). Cannot compute control.")
        return (float('nan'), float('nan'), float('nan'), float('nan'))

    # b. Draw random starts
    max_start = genome_length - stem_length + 1
    random_starts = rng.integers(1, max_start + 1, size=n_controls)

    # c. Compute control coverages
    control_coverages = []
    for r in random_starts:
        r_end = r + stem_length - 1
        cov_pct, _, _, _ = compute_coverage(r, r_end, windows_df, threshold, strand_filter)
        control_coverages.append(cov_pct)

    # d. Calculate summary stats
    control_mean = float(np.mean(control_coverages))
    control_std = float(np.std(control_coverages, ddof=1)) if n_controls > 1 else 0.0

    # e. Calculate Z-score
    if control_std == 0.0:
        z_score = float('nan')
    else:
        z_score = float((observed_coverage_pct - control_mean) / control_std)

    # f. Calculate Empirical P-value
    count_greater_equal = sum(1 for c in control_coverages if c >= observed_coverage_pct)
    empirical_pvalue = float(count_greater_equal / n_controls)

    # g. Return
    return (control_mean, control_std, z_score, empirical_pvalue)


def main():
    """Main CLI execution logic."""
    parser = argparse.ArgumentParser(description="Quantify BacTermFinder coverage over stem loops.")
    parser.add_argument("-f", "--fasta", required=True, type=Path, help="Path to stem-loop multifasta file")
    parser.add_argument("-r", "--results-dir", required=True, type=Path, help="Path to BacTermFinder results directory")
    parser.add_argument("-t", "--threshold", type=float, default=0.3, help="probability_mean threshold (default: 0.3)")
    parser.add_argument("-n", "--n-controls", type=int, default=1000, help="Random control iterations per stem loop")
    parser.add_argument("--strand", choices=["both", "plus", "minus"], default="both", help="Windows to include (default: both)")
    parser.add_argument("--output-csv", type=Path, default=Path("terminator_coverage.csv"), help="Output CSV path")
    parser.add_argument("--output-plot", type=Path, default=Path("terminator_coverage.png"), help="Output PNG path")
    parser.add_argument("--seed", type=int, default=42, help="NumPy random seed")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable DEBUG logging")
    
    args = parser.parse_args()

    # Initialise logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )

    if not args.results_dir.exists() or not args.results_dir.is_dir():
        raise FileNotFoundError(f"Results directory does not exist or is not a directory: {args.results_dir}")

    rng = np.random.default_rng(args.seed)

    # Step 1: Parse FASTA
    logging.info(f"Parsing FASTA file: {args.fasta}")
    stem_loops = []
    for record in SeqIO.parse(args.fasta, "fasta"):
        try:
            acc, start, end = parse_fasta_id(record.id)
        except ValueError as e:
            logging.warning(f"Skipping record: {e}")
            continue
            
        length = end - start + 1
        if end < start or length < 1:
            logging.error(f"Invalid coordinates for {record.id}: start={start}, end={end}. Skipping.")
            continue
            
        stem_loops.append({
            "seq_id": record.id,
            "accession": acc,
            "start": start,
            "end": end,
            "length": length
        })
    
    df_stemloops = pd.DataFrame(stem_loops)
    if df_stemloops.empty:
        logging.error("No valid FASTA records parsed. Exiting.")
        return

    # Step 2: Resolve & Cache Windows per unique accession
    unique_accessions = df_stemloops['accession'].unique()
    windows_dict: Dict[str, pd.DataFrame] = {}
    genome_lengths: Dict[str, int] = {}

    for acc in unique_accessions:
        norm_acc = normalise_accession(acc)
        if norm_acc in windows_dict:
            continue

        try:
            csv_path = resolve_csv_path(args.results_dir, acc)
        except RuntimeError as e:
            logging.error(str(e))
            raise

        if csv_path is None:
            logging.warning(f"No BacTermFinder CSV found for accession {acc} (normalised: {norm_acc})")
            windows_dict[norm_acc] = pd.DataFrame() # Sentinel for missing
            genome_lengths[norm_acc] = 0
        else:
            logging.debug(f"Loading windows from {csv_path.name}")
            win_df, max_end = load_bactermfinder_windows(csv_path)
            windows_dict[norm_acc] = win_df
            genome_lengths[norm_acc] = max_end

    # Step 3 & 4: Calculate Coverage & Random Control iteratively
    logging.info("Computing coverage and random controls...")
    results = []

    for row in df_stemloops.itertuples(index=False):
        norm_acc = normalise_accession(row.accession)
        win_df = windows_dict.get(norm_acc)

        if win_df is None or win_df.empty:
            # Handle Edge Case 1: Missing CSV
            results.append({
                'seq_id': row.seq_id, 'accession': row.accession,
                'start': row.start, 'end': row.end, 'length': row.length,
                'n_windows_above_threshold': np.nan, 'coverage_pct': np.nan,
                'max_probability_mean': np.nan, 'mean_probability_mean': np.nan,
                'control_mean_coverage_pct': np.nan, 'control_std_coverage_pct': np.nan,
                'z_score': np.nan, 'empirical_pvalue': np.nan
            })
            continue

        cov_pct, n_win, max_p, mean_p = compute_coverage(
            row.start, row.end, win_df, args.threshold, args.strand
        )

        gen_len = genome_lengths[norm_acc]
        ctrl_mean, ctrl_std, z_score, emp_pval = compute_control(
            row.length, gen_len, win_df, args.threshold, args.strand,
            args.n_controls, rng, cov_pct
        )

        results.append({
            'seq_id': row.seq_id,
            'accession': row.accession,
            'start': row.start,
            'end': row.end,
            'length': row.length,
            'n_windows_above_threshold': n_win,
            'coverage_pct': cov_pct,
            'max_probability_mean': max_p,
            'mean_probability_mean': mean_p,
            'control_mean_coverage_pct': ctrl_mean,
            'control_std_coverage_pct': ctrl_std,
            'z_score': z_score,
            'empirical_pvalue': emp_pval
        })

    df_results = pd.DataFrame(results)

    # Output generation
    float_cols = [
        'coverage_pct', 'max_probability_mean', 'mean_probability_mean',
        'control_mean_coverage_pct', 'control_std_coverage_pct',
        'z_score', 'empirical_pvalue'
    ]
    df_results[float_cols] = df_results[float_cols].round(4)
    
    logging.info(f"Saving results to {args.output_csv}")
    df_results.to_csv(args.output_csv, index=False)

    # Visualisation
    logging.info(f"Generating plot: {args.output_plot}")
    sns.set_theme(style="ticks", rc={"axes.labelsize": 10, "xtick.labelsize": 10, "ytick.labelsize": 10})
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    valid_results = df_results.dropna(subset=['coverage_pct', 'control_mean_coverage_pct']).copy()
    
    if valid_results.empty:
        logging.warning("No valid data for plotting. Generating empty figure.")
        fig.text(0.5, 0.5, 'No valid data to plot', ha='center', va='center')
        plt.tight_layout()
        plt.savefig(args.output_plot, dpi=300)
        plt.close()
        return

    # Panel A: Observed vs. Random Coverage
    axA = axes[0]
    melted = valid_results[['seq_id', 'coverage_pct', 'control_mean_coverage_pct']].melt(
        id_vars='seq_id', value_vars=['coverage_pct', 'control_mean_coverage_pct'],
        var_name='Group', value_name='Coverage'
    )
    melted['Group'] = melted['Group'].replace({
        'coverage_pct': 'Stem Loops',
        'control_mean_coverage_pct': 'Random Control'
    })

    sns.boxplot(data=melted, x='Group', y='Coverage', color='white', ax=axA, showfliers=False, width=0.4)
    sns.stripplot(data=melted, x='Group', y='Coverage', color='black', alpha=0.5, ax=axA, jitter=False)

    # Connect paired points
    for _, row in valid_results.iterrows():
        axA.plot([0, 1], [row['coverage_pct'], row['control_mean_coverage_pct']], color='grey', alpha=0.3, linewidth=0.5)

    median_obs = valid_results['coverage_pct'].median()
    median_ctrl = valid_results['control_mean_coverage_pct'].median()
    axA.text(0.05, median_obs, f"Median: {median_obs:.1f}%", va='bottom', ha='left', color='red', fontsize=9)
    axA.text(1.05, median_ctrl, f"Median: {median_ctrl:.1f}%", va='bottom', ha='left', color='red', fontsize=9)

    axA.set_ylabel("BacTermFinder Coverage (%)")
    axA.set_xlabel("")
    axA.set_title("Observed vs. Random Coverage")

    # Panel B: Coverage vs. Terminator Probability
    axB = axes[1]
    scatter_data = valid_results.copy()
    scatter_data['Covered'] = scatter_data['coverage_pct'] > 0
    palette = {}
    if True in scatter_data['Covered'].values: palette[True] = 'blue'
    if False in scatter_data['Covered'].values: palette[False] = 'red'

    sns.scatterplot(
        data=scatter_data, x='max_probability_mean', y='coverage_pct',
        hue='Covered', ax=axB, palette=palette, alpha=0.7
    )
    axB.axvline(args.threshold, color='black', linestyle='--', label='Threshold')
    axB.set_xlabel("Max Terminator Probability")
    axB.set_ylabel("BacTermFinder Coverage (%)")
    axB.set_title("Coverage vs. Terminator Probability")
    
    # Safely handle the legend in Panel B
    handles, labels = axB.get_legend_handles_labels()
    if handles:
        axB.legend()

    # Panel C: Empirical P-value Distribution
    axC = axes[2]
    pvals = valid_results['empirical_pvalue'].dropna()
    sns.histplot(pvals, bins=20, binrange=(0, 1), ax=axC, color='skyblue')
    axC.axvline(0.05, color='black', linestyle='--', label='p=0.05')
    axC.set_xlabel("Empirical P-value (vs. random)")
    axC.set_ylabel("Count")
    axC.set_title("Empirical P-value Distribution")
    axC.legend()

    plt.tight_layout()
    plt.savefig(args.output_plot, dpi=300)
    plt.close()

    logging.info("Analysis completed successfully.")

if __name__ == "__main__":
    main()
