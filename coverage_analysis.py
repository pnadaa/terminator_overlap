"""
BacTermFinder Coverage Analysis Script (strand-aware)

Quantifies coverage of bacterial genomic regions by BacTermFinder terminator
windows and compares observed regions to a random genomic background.

Plotting:
- --plot writes the new probability-focused 3-panel seaborn figure
- --legacy-plot writes the previous coverage/probability legacy figure

New probability-focused plot:
1) Observed vs random best qualifying BacTermFinder window probability_mean
   using windows that overlap the input region by > --plot-min-window-overlap-frac
2) Distribution of per-query overlapping-window probability_mean summaries
3) Ranked empirical p-values based on the random best-window null distribution
"""

import argparse
import logging
import multiprocessing as mp
import os
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from intervaltree import IntervalTree


_WORKER_STATE = {}


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def run_cmd(cmd, *, text=True):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=text)
    if p.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit {p.returncode}): {' '.join(cmd)}\n"
            f"STDERR:\n{p.stderr}\nSTDOUT:\n{p.stdout}"
        )
    return p.stdout

def get_available_cpus() -> int:
    """
    Resolve available CPU count for the current process in priority order:
      1. PBS_NCPUS env var (set by PBS/Torque, most explicit)
      2. os.sched_getaffinity (Linux kernel affinity mask)
      3. os.cpu_count() (total node CPUs, last resort)
    """
    pbs_ncpus = os.environ.get("PBS_NCPUS") or os.environ.get("NCPUS")
    if pbs_ncpus:
        try:
            return int(pbs_ncpus)
        except ValueError:
            pass
    if hasattr(os, "sched_getaffinity"):
        try:
            n = len(os.sched_getaffinity(0))
            if n > 0:
                return n
        except OSError:
            pass
    return os.cpu_count() or 1


def normalise_accession(acc: str) -> str:
    return re.sub(r"\.\d+$", "", acc)


def prefer_dot1_candidates(genome_id: str):
    g = genome_id.strip()
    if not g:
        return []
    if "." in g:
        base = g.rsplit(".", 1)[0]
        return [f"{base}.1", base]
    return [f"{g}.1", g]


def resolve_mean_csv(acc_raw: str, mean_root: Path, mean_suffix: str):
    mean_tried = []
    for g in prefer_dot1_candidates(acc_raw):
        for candidate_suffix in [mean_suffix, mean_suffix + ".gz"]:
            p = mean_root / g / f"{g}{candidate_suffix}"
            mean_tried.append(str(p))
            if p.exists():
                return {
                    "genome_key_used": g,
                    "mean_csv": p,
                    "mean_tried": mean_tried,
                }
    return {
        "genome_key_used": None,
        "mean_csv": None,
        "mean_tried": mean_tried,
    }


def blastdb_length(blastdbcmd: str, blast_db: str, entry: str) -> int:
    cmd = [blastdbcmd, "-db", blast_db, "-entry", entry, "-outfmt", "%l"]
    out = run_cmd(cmd).strip()
    try:
        return int(out)
    except ValueError:
        raise RuntimeError(f"Could not parse blastdbcmd length for entry={entry!r}: {out!r}")


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_fasta_coord_id(record_id: str):
    """
    Parse FASTA ID: {accession}_{start}-{end}

    Strand is inferred from raw coordinate order:
      raw_start <= raw_end -> '+'
      raw_start > raw_end  -> '-'

    Returns: (acc, start, end, strand) with start <= end
    """
    if "_" not in record_id:
        raise ValueError(f"FASTA id missing '_' separator: {record_id}")
    acc, coord_part = record_id.split("_", 1)
    if "-" not in coord_part:
        raise ValueError(f"FASTA id missing '-' in coords: {record_id}")
    a, b = coord_part.split("-", 1)
    raw_start = int(a)
    raw_end = int(b)
    if raw_start > raw_end:
        strand = "-"
        raw_start, raw_end = raw_end, raw_start
    else:
        strand = "+"
    return acc, min(raw_start, raw_end), max(raw_start, raw_end), strand


def parse_window_samplename(sample_name: str):
    """
    Supports:
      window_genome_low_high
      window_genome_low_high_strand

    Returns: (window_id, genome, start, end, strand_or_None)
    """
    parts = sample_name.split("_")
    if len(parts) < 4:
        raise ValueError(f"SampleName too short to parse window coords: {sample_name}")

    maybe_strand = parts[-1]
    if maybe_strand in {"+", "-"}:
        high = int(parts[-2])
        low = int(parts[-3])
        genome = "_".join(parts[1:-3])
        window = parts[0]
        strand = maybe_strand
    else:
        high = int(parts[-1])
        low = int(parts[-2])
        genome = "_".join(parts[1:-2])
        window = parts[0]
        strand = None

    start = min(low, high)
    end = max(low, high)
    return window, genome, start, end, strand


# ---------------------------------------------------------------------------
# IntervalTree loading
# ---------------------------------------------------------------------------

def load_positive_windows(
    mean_csv: Path,
    threshold: float,
    genome_label: str,
    *,
    merge: bool = True,
    keep_debug_data: bool = False,
    log_diagnostics: bool = True,
):
    """
    Load mean_csv and build per-strand IntervalTrees for windows with
    probability_mean >= threshold.

    Returns:
      trees       : dict {'+': IntervalTree, '-': IntervalTree}
      total_bp    : dict {'+': int, '-': int}
      kept_rows   : int
      kept_minmax : (min_coord, max_coord)
    """
    df = pd.read_csv(mean_csv)
    if "SampleName" not in df.columns or "probability_mean" not in df.columns:
        raise ValueError(f"{mean_csv} missing required columns SampleName/probability_mean")

    df["probability_mean"] = df["probability_mean"].astype(float)
    kept = df[df["probability_mean"] >= threshold].copy()

    trees = {"+": IntervalTree(), "-": IntervalTree()}
    kept_min = None
    kept_max = None
    n_total = len(kept)
    n_stranded = 0
    n_genome_mismatch = 0
    norm_label = normalise_accession(genome_label)

    for sn, pm in zip(kept["SampleName"].astype(str), kept["probability_mean"].astype(float)):
        try:
            _, genome_in_sn, s, e, strand = parse_window_samplename(sn)
        except Exception as exc:
            logging.warning("%s: could not parse SampleName %r: %s", genome_label, sn, exc)
            continue

        if normalise_accession(genome_in_sn) != norm_label:
            n_genome_mismatch += 1

        kept_min = s if kept_min is None else min(kept_min, s)
        kept_max = e if kept_max is None else max(kept_max, e)

        data = None
        if keep_debug_data:
            data = {
                "SampleName": sn,
                "probability_mean": float(pm),
                "Genome": genome_in_sn,   # restore
                "Strand": strand,          # restore
            }

        if strand is not None:
            n_stranded += 1

        target_strands = [strand] if strand is not None else ["+", "-"]
        for ts in target_strands:
            trees[ts].addi(s, e + 1, data)

    if log_diagnostics:
        if n_total == 0:
            logging.warning("%s: No windows passed the threshold — trees are empty.", genome_label)
        elif n_stranded == 0:
            logging.warning(
                "%s: 0/%d windows carry strand information; strand-specific queries "
                "will behave like strand-naive queries because both trees contain the "
                "same windows. Verify SampleName strand encoding.",
                genome_label, n_total,
            )
        elif n_stranded < n_total:
            logging.warning(
                "%s: %d/%d windows carry strand information; %d unstranded window(s) "
                "inserted into both '+' and '-' trees.",
                genome_label, n_stranded, n_total, n_total - n_stranded,
            )
        else:
            logging.debug("%s: All %d windows carry strand information.", genome_label, n_total)

        if n_genome_mismatch > 0:
            logging.warning(
                "%s: %d/%d windows have a genome ID in SampleName that does not match "
                "the expected accession %r after version-stripping.",
                genome_label, n_genome_mismatch, n_total, genome_label,
            )

    if merge:
        for tree in trees.values():
            tree.merge_overlaps(strict=False, data_reducer=lambda a, b: None)

    total_bp = {s: int(sum(iv.end - iv.begin for iv in t)) for s, t in trees.items()}
    return trees, total_bp, int(n_total), (kept_min, kept_max)


# ---------------------------------------------------------------------------
# NumPy-backed vectorized overlap (for random control)
# ---------------------------------------------------------------------------

def build_merged_arrays(merged_trees: dict, strand: str):
    """Convert a merged IntervalTree strand to sorted NumPy arrays (begin, end exclusive)."""
    ivs = sorted(merged_trees[strand], key=lambda iv: iv.begin)
    if not ivs:
        return np.empty((0,), dtype=np.int64), np.empty((0,), dtype=np.int64)
    starts = np.array([iv.begin for iv in ivs], dtype=np.int64)
    ends = np.array([iv.end for iv in ivs], dtype=np.int64)
    return starts, ends


def vectorized_overlap_bp(
    iv_starts: np.ndarray,
    iv_ends: np.ndarray,
    r_starts: np.ndarray,
    region_len: int,
) -> np.ndarray:
    """
    Compute overlap_bp for all random starts simultaneously.

    iv_starts, iv_ends : merged interval arrays, end-exclusive
    r_starts           : (N,) int64 array of random region starts (1-based inclusive)
    region_len         : region length in bp

    Returns:
      (N,) int64 array of total overlap bp per random placement
    """
    r_starts = np.asarray(r_starts, dtype=np.int64)
    if region_len <= 0 or r_starts.size == 0 or iv_starts.size == 0:
        return np.zeros(r_starts.shape[0], dtype=np.int64)

    r_ends = r_starts + region_len
    a0 = np.maximum(iv_starts[None, :], r_starts[:, None])
    a1 = np.minimum(iv_ends[None, :], r_ends[:, None])
    overlaps = np.maximum(0, a1 - a0).sum(axis=1)
    return overlaps.astype(np.int64, copy=False)


# ---------------------------------------------------------------------------
# Overlap and probability metrics
# ---------------------------------------------------------------------------

def interval_overlap_bp(iv_begin: int, iv_end_exclusive: int, start_inclusive: int, end_inclusive: int) -> int:
    q0, q1 = int(start_inclusive), int(end_inclusive) + 1
    a0 = max(iv_begin, q0)
    a1 = min(iv_end_exclusive, q1)
    return int(max(0, a1 - a0))


def overlap_bp(trees_by_strand: dict, query_strand: str, start_inclusive: int, end_inclusive: int) -> int:
    tree = trees_by_strand.get(query_strand)
    if tree is None:
        return 0
    total = 0
    q0, q1 = int(start_inclusive), int(end_inclusive) + 1
    if q1 <= q0:
        return 0
    for iv in tree.overlap(q0, q1):
        total += interval_overlap_bp(iv.begin, iv.end, start_inclusive, end_inclusive)
    return int(total)


def get_raw_overlaps(raw_trees_by_strand: dict, query_strand: str, start_inclusive: int, end_inclusive: int):
    tree = raw_trees_by_strand.get(query_strand)
    if tree is None:
        return []
    q0, q1 = int(start_inclusive), int(end_inclusive) + 1
    return sorted(tree.overlap(q0, q1), key=lambda iv: (iv.begin, iv.end))


def compute_probability_metrics(raw_trees_by_strand: dict, query_strand: str, start_inclusive: int, end_inclusive: int):
    """
    Returns:
      (n_windows, max_probability_mean, mean_probability_mean, median_probability_mean)
    """
    overlaps = get_raw_overlaps(raw_trees_by_strand, query_strand, start_inclusive, end_inclusive)
    if not overlaps:
        return 0, float("nan"), float("nan"), float("nan")

    probs = [
        float(iv.data["probability_mean"])
        for iv in overlaps
        if iv.data is not None and "probability_mean" in iv.data
    ]
    if not probs:
        return len(overlaps), float("nan"), float("nan"), float("nan")

    arr = np.array(probs, dtype=float)
    return len(arr), float(np.max(arr)), float(np.mean(arr)), float(np.median(arr))


def best_overlap_probability_window(
    raw_trees_by_strand: dict,
    query_strand: str,
    start_inclusive: int,
    end_inclusive: int,
    min_overlap_frac: float = 0.10,
):
    """
    Identify the overlapping raw window with the highest probability_mean among
    windows overlapping the input region by > min_overlap_frac of the input region.

    Returns a dict with score and window metadata.
    If no qualifying overlap exists, probability_mean is 0.0 and metadata is None.
    """
    region_len = int(end_inclusive) - int(start_inclusive) + 1
    best = {
        "probability_mean": 0.0,
        "window_sample": None,
        "window_start": None,
        "window_end_inclusive": None,
        "window_overlap_bp": 0,
        "window_overlap_frac": 0.0,
    }
    if region_len <= 0:
        return best

    overlaps = get_raw_overlaps(raw_trees_by_strand, query_strand, start_inclusive, end_inclusive)
    for iv in overlaps:
        ov_bp = interval_overlap_bp(iv.begin, iv.end, start_inclusive, end_inclusive)
        ov_frac = ov_bp / region_len if region_len > 0 else 0.0
        if ov_bp <= 0 or ov_frac <= min_overlap_frac:
            continue
        if iv.data is None or "probability_mean" not in iv.data:
            continue

        prob = float(iv.data["probability_mean"])
        candidate = {
            "probability_mean": prob,
            "window_sample": iv.data.get("SampleName"),
            "window_start": int(iv.begin),
            "window_end_inclusive": int(iv.end - 1),
            "window_overlap_bp": int(ov_bp),
            "window_overlap_frac": float(ov_frac),
        }

        if (
            candidate["probability_mean"] > best["probability_mean"]
            or (
                candidate["probability_mean"] == best["probability_mean"]
                and candidate["window_overlap_bp"] > best["window_overlap_bp"]
            )
            or (
                candidate["probability_mean"] == best["probability_mean"]
                and candidate["window_overlap_bp"] == best["window_overlap_bp"]
                and best["window_start"] is not None
                and candidate["window_start"] < best["window_start"]
            )
            or (
                candidate["probability_mean"] == best["probability_mean"]
                and candidate["window_overlap_bp"] == best["window_overlap_bp"]
                and best["window_start"] is None
            )
        ):
            best = candidate

    return best


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def annotate_mean_and_median(ax, xloc, series, color_mean="darkgreen", color_median="firebrick"):
    vals = pd.Series(series).dropna().astype(float)
    if vals.empty:
        return
    mean_v = float(vals.mean())
    median_v = float(vals.median())
    # Use blended transform: x in axes fraction, y in data coordinates
    trans = ax.get_yaxis_transform()
    offset = 0.02  # small fixed axes-fraction nudge
    ax.text(
        xloc + offset, mean_v, f"Mean: {mean_v:.3f}",
        va="bottom", ha="left", color=color_mean, fontsize=8,
        clip_on=True, transform=trans,
    )
    va_median = "top" if abs(mean_v - median_v) > 0.02 else "bottom"
    ax.text(
        xloc + offset, median_v, f"Median: {median_v:.3f}",
        va=va_median, ha="left", color=color_median, fontsize=8,
        clip_on=True, transform=trans,
    )



_PALETTE_OBS_RAND = {
    "Observed best window": "#4C72B0",
    "Random best window mean": "#DD8452",
}

_PALETTE_METRICS = {
    "Mean overlap prob": "#4C72B0",
    "Median overlap prob": "#55A868",
    "Max overlap prob": "#C44E52",
    "Best qualifying prob": "#8172B2",
}

_PALETTE_COV = {
    "Observed": "#4C72B0",
    "Random mean": "#DD8452",
}


def write_new_plot(res_df: pd.DataFrame, out_prefix: Path, threshold: float, min_overlap_frac: float, random_n: int = 0):
    import matplotlib
    import matplotlib.pyplot as plt
    import seaborn as sns

    matplotlib.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"]})
    sns.set_theme(
        style="ticks",
        rc={"axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9},
    )

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), constrained_layout=True)
    fig.patch.set_facecolor("white")

    # ------------------------------------------------------------------
    # Panel A: Observed vs random best qualifying probability_mean
    # ------------------------------------------------------------------
    axA = axes[0]
    validA = res_df.dropna(
        subset=["best_overlap_probability_mean", "random_best_probability_mean_mean"]
    ).copy()

    if validA.empty:
        axA.text(
            0.5, 0.5, "No valid observed/random\nbest-window data",
            ha="center", va="center", transform=axA.transAxes
        )
        axA.set_axis_off()
    else:
        panelA = validA[
            ["qseqid", "best_overlap_probability_mean", "random_best_probability_mean_mean"]
        ].melt(
            id_vars="qseqid",
            value_vars=["best_overlap_probability_mean", "random_best_probability_mean_mean"],
            var_name="Group",
            value_name="ProbabilityMean",
        )
        panelA["Group"] = panelA["Group"].replace({
            "best_overlap_probability_mean": "Observed best window",
            "random_best_probability_mean_mean": "Random best window mean",
        })

        sns.boxplot(
            data=panelA, x="Group", y="ProbabilityMean", hue="Group",
            palette=_PALETTE_OBS_RAND, legend=False, width=0.45, showfliers=False, ax=axA,
            boxprops=dict(alpha=0.35),
        )
        sns.stripplot(
            data=panelA, x="Group", y="ProbabilityMean", hue="Group",
            palette=_PALETTE_OBS_RAND, legend=False, alpha=0.75, jitter=0.08, ax=axA,
        )

        for _, row in validA.iterrows():
            axA.plot(
                [0, 1],
                [row["best_overlap_probability_mean"], row["random_best_probability_mean_mean"]],
                color="grey", alpha=0.10, linewidth=0.3,
            )

        annotate_mean_and_median(axA, 0, validA["best_overlap_probability_mean"])
        annotate_mean_and_median(axA, 1, validA["random_best_probability_mean_mean"])

        axA.axhline(
            threshold, linestyle="--", color="black", linewidth=1,
            label=f"Threshold ({threshold})"
        )
        axA.legend(fontsize=8, loc="upper right")
        axA.set_xlabel("")
        axA.set_ylabel("Best qualifying window probability_mean")
        n_zero = (validA["best_overlap_probability_mean"] == 0.0).sum()
        axA.set_title(
            f"Observed vs random best window\n"
            f"(qualifying overlap > {min_overlap_frac:.0%}; "
            f"n_no_overlap={n_zero})"
        )

    # ------------------------------------------------------------------
    # Panel B: Overlapping-window probability_mean summaries
    # ------------------------------------------------------------------
    axB = axes[1]
    validB = res_df.dropna(
        subset=["mean_probability_mean", "median_probability_mean", "max_probability_mean"]
    ).copy()

    if validB.empty:
        axB.text(
            0.5, 0.5, "No overlapping-window\nprobability_mean data",
            ha="center", va="center", transform=axB.transAxes
        )
        axB.set_axis_off()
    else:
        panelB = validB[
            ["qseqid", "mean_probability_mean", "median_probability_mean",
             "max_probability_mean", "best_overlap_probability_mean"]
        ].melt(
            id_vars="qseqid",
            value_vars=[
                "mean_probability_mean",
                "median_probability_mean",
                "max_probability_mean",
                "best_overlap_probability_mean",
            ],
            var_name="Metric",
            value_name="ProbabilityMean",
        )
        panelB["Metric"] = panelB["Metric"].replace({
            "mean_probability_mean": "Mean overlap prob",
            "median_probability_mean": "Median overlap prob",
            "max_probability_mean": "Max overlap prob",
            "best_overlap_probability_mean": "Best qualifying prob",
        })

        n_per_metric = panelB.groupby("Metric")["ProbabilityMean"].count().min()
        if n_per_metric >= 5:
            sns.violinplot(
                data=panelB, x="Metric", y="ProbabilityMean",  hue="Metric",
                palette=_PALETTE_METRICS, legend=False, inner="box", cut=0, ax=axB,
            )
        else:
            sns.boxplot(
                data=panelB, x="Metric", y="ProbabilityMean", hue="Metric",
                palette=_PALETTE_METRICS, legend=False, width=0.55, showfliers=False, ax=axB,
                boxprops=dict(alpha=0.35),
            )

        sns.stripplot(
            data=panelB, x="Metric", y="ProbabilityMean",  hue="Metric",
            palette=_PALETTE_METRICS, legend=False, alpha=0.6, jitter=0.12, size=3, ax=axB,
        )

        axB.axhline(
            threshold, linestyle="--", color="black", linewidth=1,
            label=f"Threshold ({threshold})"
        )
        axB.legend(fontsize=8, loc="upper right")
        axB.set_xlabel("")
        axB.set_ylabel("probability_mean")
        axB.set_title("Per-query overlapping-window probability summaries")
        axB.tick_params(axis="x", rotation=30)

    # ------------------------------------------------------------------
    # Panel C: Ranked empirical p-values from random best-window null
    # ------------------------------------------------------------------
    axC = axes[2]
    validC = res_df.dropna(subset=["empirical_pvalue_best_overlap"]).copy()

    if validC.empty:
        axC.text(
            0.5, 0.5, "No empirical p-values\n(random control unavailable)",
            ha="center", va="center", transform=axC.transAxes
        )
        axC.set_axis_off()
    else:
        panelC = validC.sort_values("empirical_pvalue_best_overlap").copy()
        panelC["rank"] = np.arange(1, len(panelC) + 1)
        panelC["neglog10_p"] = -np.log10(panelC["empirical_pvalue_best_overlap"])
        panelC["Significant"] = panelC["empirical_pvalue_best_overlap"] <= 0.05

        sig_df    = panelC[panelC["Significant"]]
        nonsig_df = panelC[~panelC["Significant"]]

        # ── Non-significant: rasterized so 100k+ overlapping points don't
        #    accumulate alpha and wash out. Rendered as a bitmap layer inside
        #    the SVG — no transparency stacking, tiny file size. ──────────────
        axC.scatter(
            nonsig_df["rank"],
            nonsig_df["neglog10_p"],
            color="steelblue",
            s=4,                   # small markers reduce apparent clutter
            alpha=0.5,
            linewidths=0,
            label=f"Not significant (n={len(nonsig_df):,})",
            rasterized=True,       # ← key fix: bitmap layer, no alpha accumulation
            zorder=3,
        )

        # ── Significant: kept as crisp vectors on top ────────────────────────
        if not sig_df.empty:
            axC.scatter(
                sig_df["rank"],
                sig_df["neglog10_p"],
                color="firebrick",
                s=30,
                alpha=1.0,
                linewidths=0.4,
                edgecolors="darkred",
                label=f"Significant (n={len(sig_df):,})",
                rasterized=False,  # vector — stays sharp at any zoom
                zorder=4,
            )

        axC.axhline(-np.log10(0.05), linestyle="--", color="black", linewidth=1, label="p = 0.05")
        axC.axhline(-np.log10(0.10), linestyle=":", color="grey", linewidth=0.8, label="p = 0.10")
        if random_n > 0:
            min_detectable_p = 1.0 / (random_n + 1)
            axC.axhline(
                -np.log10(min_detectable_p), linestyle="-.",
                color="purple", linewidth=1.0,
                label=f"Min detectable p (n={random_n})",
            )

        axC.grid(True, axis="y", alpha=0.3, linestyle="--", zorder=0)
        axC.set_xlabel(f"Query rank (n={len(panelC):,})")
        axC.set_ylabel("-log₁₀(empirical p-value)")
        axC.set_title("Empirical p-values from best-window randomisation")
        handles, labels = axC.get_legend_handles_labels()
        if handles:
            axC.legend(fontsize=8)

    sns.despine(fig=fig)

    plot_path = out_prefix.with_suffix(".overlap.png")
    fig.savefig(plot_path, dpi=300, bbox_inches="tight", facecolor="white")
    svg_path = out_prefix.with_suffix(".overlap.svg")
    fig.savefig(svg_path, bbox_inches="tight", facecolor="white", dpi=300)
    plt.close(fig)
    logging.info("Wrote new plot: %s (+ .svg)", plot_path)


def write_legacy_plot(res_df: pd.DataFrame, out_prefix: Path, threshold: float):
    import matplotlib
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    import seaborn as sns

    matplotlib.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"]})
    sns.set_theme(
        style="ticks",
        rc={"axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9},
    )

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), constrained_layout=True)
    fig.patch.set_facecolor("white")

    valid = res_df.dropna(subset=["percent_region_covered", "random_mean_percent"]).copy()

    if valid.empty:
        fig.text(0.5, 0.5, "No valid data to plot", ha="center", va="center")
    else:
        # Panel A: coverage observed vs random
        axA = axes[0]
        melted = valid[
            ["qseqid", "percent_region_covered", "random_mean_percent"]
        ].melt(
            id_vars="qseqid",
            value_vars=["percent_region_covered", "random_mean_percent"],
            var_name="Group",
            value_name="Coverage",
        )
        melted["Group"] = melted["Group"].replace({
            "percent_region_covered": "Observed",
            "random_mean_percent": "Random mean",
        })

        sns.boxplot(
            data=melted, x="Group", y="Coverage", hue="Group",
            palette=_PALETTE_COV, legend=False, ax=axA, showfliers=False, width=0.4,
            boxprops=dict(alpha=0.35),
        )
        sns.stripplot(
            data=melted, x="Group", y="Coverage", hue="Group",
            palette=_PALETTE_COV, legend=False, alpha=0.65, ax=axA, jitter=0.08,
        )

        for _, row in valid.iterrows():
            axA.plot(
                [0, 1],
                [row["percent_region_covered"], row["random_mean_percent"]],
                color="grey", alpha=0.3, linewidth=0.5,
            )

        annotate_mean_and_median(axA, 0, valid["percent_region_covered"])
        annotate_mean_and_median(axA, 1, valid["random_mean_percent"])
        axA.set_xlabel("")
        axA.set_ylabel("Coverage (%)")
        axA.set_title("Legacy: observed vs random coverage")

        # Panel B: coverage vs max probability
        axB = axes[1]
        scatter_data = valid.dropna(subset=["max_probability_mean"]).copy()
        if scatter_data.empty:
            axB.text(0.5, 0.5, "No max probability data", ha="center", va="center", transform=axB.transAxes)
        else:
            norm = mcolors.Normalize(
                vmin=scatter_data["percent_region_covered"].min(),
                vmax=scatter_data["percent_region_covered"].max(),
            )
            cmap = cm.viridis
            sc = axB.scatter(
                scatter_data["max_probability_mean"],
                scatter_data["percent_region_covered"],
                c=scatter_data["percent_region_covered"],
                cmap=cmap,
                alpha=0.75,
                s=40,
                zorder=3,
            )
            cbar = fig.colorbar(sc, ax=axB, shrink=0.85)
            cbar.set_label("Coverage (%)", fontsize=8)
            axB.axvline(threshold, color="black", linestyle="--", label=f"Threshold ({threshold})")
            axB.legend(fontsize=8)
            axB.set_xlabel("Max probability_mean")
            axB.set_ylabel("Coverage (%)")
            axB.set_title("Legacy: coverage vs max probability")

        # Panel C: p-value histogram + KDE
        axC = axes[2]
        pvals = valid["empirical_pvalue_best_overlap"].dropna()
        if not pvals.empty:
            sns.histplot(
                pvals, bins=20, binrange=(0, 1), ax=axC, color="skyblue",
                stat="density", kde=True, line_kws={"linewidth": 1.5}
            )
            axC.axvline(0.05, color="black", linestyle="--", label="p = 0.05")
            axC.legend(fontsize=8)
        else:
            axC.text(0.5, 0.5, "No p-values", ha="center", va="center", transform=axC.transAxes)
        axC.set_xlabel("Empirical p-value")
        axC.set_ylabel("Density")
        axC.set_title("Legacy: p-value distribution")

    sns.despine(fig=fig)

    legacy_path = out_prefix.with_suffix(".legacy.overlap.png")
    fig.savefig(legacy_path, dpi=300, bbox_inches="tight", facecolor="white")
    svg_path = out_prefix.with_suffix(".legacy.overlap.svg")
    fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    logging.info("Wrote legacy plot: %s (+ .svg)", legacy_path)


# ---------------------------------------------------------------------------
# Parallel worker state + per-query processing
# ---------------------------------------------------------------------------

def set_worker_state(state: dict):
    global _WORKER_STATE
    _WORKER_STATE = state


def process_one_query(job):
    state = _WORKER_STATE

    index, qid, acc_raw, region_start, region_end, query_strand = job
    region_len = region_end - region_start + 1
    nanv = float("nan")

    mean_resolution = state["mean_resolution_cache"][acc_raw]
    genome_key_used = mean_resolution["genome_key_used"]
    mean_csv = mean_resolution["mean_csv"]
    mean_tried = mean_resolution["mean_tried"]

    if mean_csv is None:
        note = f"mean.csv not found for accession={acc_raw!r}"
        result = {
            "qseqid": qid,
            "accession_from_fasta": acc_raw,
            "region_start": region_start,
            "region_end": region_end,
            "query_strand": query_strand,
            "region_len": region_len,
            "genome_key_used": None,
            "mean_csv": None,
            "windows_rows_kept": None,
            "windows_total_bp_this_strand": None,
            "n_windows_above_threshold": None,
            "max_probability_mean": nanv,
            "mean_probability_mean": nanv,
            "median_probability_mean": nanv,
            "best_overlap_probability_mean": nanv,
            "best_overlap_window_sample": None,
            "best_overlap_window_start": None,
            "best_overlap_window_end_inclusive": None,
            "best_overlap_window_overlap_bp": None,
            "best_overlap_window_overlap_frac": nanv,
            "overlap_bp": 0,
            "percent_region_covered": nanv,
            "random_mean_percent": nanv,
            "random_sd_percent": nanv,
            "random_best_probability_mean_mean": nanv,
            "random_best_probability_mean_median": nanv,
            "random_best_probability_mean_sd": nanv,
            "z_score": nanv,
            "empirical_pvalue": nanv,
            "empirical_pvalue_best_overlap": nanv,
            "blast_entry_used": None,
            "blast_entry_len": None,
            "kept_window_min_coord": None,
            "kept_window_max_coord": None,
            "note": note,
        }
        debug = {
            "qseqid": qid,
            "accession": acc_raw,
            "mean_paths_tried": ";".join(mean_tried),
            "note": note,
        }
        mean_missing_row = {
            "qseqid": qid,
            "accession": acc_raw,
            "mean_paths_tried": ";".join(mean_tried),
        }
        return {
            "result": result,
            "debug": debug,
            "mean_missing_row": mean_missing_row,
            "length_failed_row": None,
            "overlap_window_rows": [],
        }

    merged_entry = state["merged_cache"][genome_key_used]
    raw_trees_for_genome = state["raw_cache"][genome_key_used]

    merged_trees = merged_entry["merged_trees"]
    total_bp_by_strand = merged_entry["total_bp_by_strand"]
    kept_rows = merged_entry["kept_rows"]
    mean_csv_str = merged_entry["mean_csv_str"]
    kept_min = merged_entry["kept_min"]
    kept_max = merged_entry["kept_max"]
    merged_arrays = merged_entry["merged_arrays"]

    windows_total_bp = total_bp_by_strand.get(query_strand, 0)

    ov = overlap_bp(merged_trees, query_strand, region_start, region_end)
    percent_region_covered = 100.0 * ov / region_len if region_len > 0 else float("nan")

    n_windows, max_prob, mean_prob, median_prob = compute_probability_metrics(
        raw_trees_for_genome, query_strand, region_start, region_end
    )

    best_obs = best_overlap_probability_window(
        raw_trees_for_genome,
        query_strand,
        region_start,
        region_end,
        min_overlap_frac=state["plot_min_window_overlap_frac"],
    )

    subj_len = None
    blast_entry_used = None
    note = ""
    rand_percents = []
    rand_best_probs = []
    length_failed_row = None

    random_n = state["random_n"]
    if random_n > 0:
        if state["blast_db"] is None:
            note = "Random control disabled (no --blast-db provided)."
        else:
            len_info = state["length_resolution_cache"].get(acc_raw, {})
            subj_len = len_info.get("subj_len")
            blast_entry_used = len_info.get("blast_entry_used")

            if subj_len is None:
                candidates = len_info.get("candidates", prefer_dot1_candidates(acc_raw))
                note = len_info.get("note") or f"blastdbcmd length lookup failed for candidates={candidates}"
                length_failed_row = {
                    "qseqid": qid,
                    "accession": acc_raw,
                    "candidates": ";".join(candidates),
                }
            elif subj_len >= region_len > 0:
                max_start = subj_len - region_len + 1
                rng = np.random.default_rng(state["seed"] + index)
                r_starts = rng.integers(1, max_start + 1, size=random_n, dtype=np.int64)

                iv_starts, iv_ends = merged_arrays[query_strand]
                rand_bp = vectorized_overlap_bp(iv_starts, iv_ends, r_starts, region_len)
                rand_percents = (100.0 * rand_bp / region_len).astype(float).tolist()

                for r_start in r_starts:
                    r_start = int(r_start)
                    r_end = r_start + region_len - 1
                    best_rand = best_overlap_probability_window(
                        raw_trees_for_genome,
                        query_strand,
                        r_start,
                        r_end,
                        min_overlap_frac=state["plot_min_window_overlap_frac"],
                    )
                    rand_best_probs.append(float(best_rand["probability_mean"]))
            else:
                note = (
                    f"Genome length {subj_len} < region_len {region_len}; "
                    "random control skipped."
                )

    if rand_percents:
        arr_cov = np.array(rand_percents, dtype=float)
        mean_r = float(np.mean(arr_cov))
        sd_r = float(np.std(arr_cov, ddof=1)) if len(arr_cov) > 1 else 0.0
    else:
        mean_r = float("nan")
        sd_r = float("nan")

    if rand_best_probs:
        arr_best = np.array(rand_best_probs, dtype=float)
        rand_best_mean = float(np.mean(arr_best))
        rand_best_median = float(np.median(arr_best))
        rand_best_sd = float(np.std(arr_best, ddof=1)) if len(arr_best) > 1 else 0.0
        obs_best_prob = float(best_obs["probability_mean"])
        z_score = float((obs_best_prob - rand_best_mean) / rand_best_sd) if rand_best_sd > 0.0 else float("nan")
        count_ge = int(np.sum(arr_best >= obs_best_prob))
        empirical_pvalue = float((count_ge + 1) / (len(arr_best) + 1))
    else:
        rand_best_mean = float("nan")
        rand_best_median = float("nan")
        rand_best_sd = float("nan")
        z_score = float("nan")
        empirical_pvalue = float("nan")

    result = {
        "qseqid": qid,
        "accession_from_fasta": acc_raw,
        "region_start": region_start,
        "region_end": region_end,
        "query_strand": query_strand,
        "region_len": region_len,
        "genome_key_used": genome_key_used,
        "mean_csv": mean_csv_str,
        "windows_rows_kept": kept_rows,
        "windows_total_bp_this_strand": windows_total_bp,
        "n_windows_above_threshold": n_windows,
        "max_probability_mean": max_prob,
        "mean_probability_mean": mean_prob,
        "median_probability_mean": median_prob,
        "best_overlap_probability_mean": float(best_obs["probability_mean"]),
        "best_overlap_window_sample": best_obs["window_sample"],
        "best_overlap_window_start": best_obs["window_start"],
        "best_overlap_window_end_inclusive": best_obs["window_end_inclusive"],
        "best_overlap_window_overlap_bp": best_obs["window_overlap_bp"],
        "best_overlap_window_overlap_frac": best_obs["window_overlap_frac"],
        "overlap_bp": ov,
        "percent_region_covered": percent_region_covered,
        "random_mean_percent": mean_r,
        "random_sd_percent": sd_r,
        "random_best_probability_mean_mean": rand_best_mean,
        "random_best_probability_mean_median": rand_best_median,
        "random_best_probability_mean_sd": rand_best_sd,
        "z_score": z_score,
        "empirical_pvalue": empirical_pvalue,
        "empirical_pvalue_best_overlap": empirical_pvalue,
        "blast_entry_used": blast_entry_used,
        "blast_entry_len": subj_len,
        "kept_window_min_coord": kept_min,
        "kept_window_max_coord": kept_max,
        "note": note,
    }

    debug = {
        "qseqid": qid,
        "accession": acc_raw,
        "region_start": region_start,
        "region_end": region_end,
        "query_strand": query_strand,
        "region_len": region_len,
        "genome_key_used": genome_key_used,
        "mean_csv_used": mean_csv_str,
        "mean_paths_tried": ";".join(mean_tried),
        "kept_rows": kept_rows,
        "kept_min_coord": kept_min,
        "kept_max_coord": kept_max,
        "windows_total_bp_this_strand": windows_total_bp,
        "n_windows_above_threshold": n_windows,
        "max_probability_mean": max_prob,
        "mean_probability_mean": mean_prob,
        "median_probability_mean": median_prob,
        "best_overlap_probability_mean": float(best_obs["probability_mean"]),
        "best_overlap_window_sample": best_obs["window_sample"],
        "best_overlap_window_overlap_bp": best_obs["window_overlap_bp"],
        "best_overlap_window_overlap_frac": best_obs["window_overlap_frac"],
        "overlap_bp": ov,
        "percent_region_covered": percent_region_covered,
        "random_mean_percent": mean_r,
        "random_best_probability_mean_mean": rand_best_mean,
        "random_best_probability_mean_median": rand_best_median,
        "z_score": z_score,
        "empirical_pvalue_best_overlap": empirical_pvalue,
        "blast_entry_used": blast_entry_used,
        "blast_entry_len": subj_len,
        "note": note,
    }

    overlap_window_rows = []
    debug_max_overlap_windows = state["debug_max_overlap_windows"]
    if debug_max_overlap_windows and debug_max_overlap_windows > 0:
        raw_overlaps = get_raw_overlaps(raw_trees_for_genome, query_strand, region_start, region_end)
        for iv in raw_overlaps[:debug_max_overlap_windows]:
            d = iv.data or {}
            ov_bp_single = interval_overlap_bp(iv.begin, iv.end, region_start, region_end)
            ov_frac_single = ov_bp_single / region_len if region_len > 0 else float("nan")
            qualifies = ov_frac_single > state["plot_min_window_overlap_frac"]
            overlap_window_rows.append({
                "qseqid": qid,
                "genome_key_used": genome_key_used,
                "region_start": region_start,
                "region_end": region_end,
                "query_strand": query_strand,
                "window_start": iv.begin,
                "window_end_inclusive": iv.end - 1,
                "window_sample": d.get("SampleName"),
                "window_probability_mean": d.get("probability_mean"),
                "window_overlap_bp": ov_bp_single,
                "window_overlap_frac": ov_frac_single,
                "window_passes_best_overlap_filter": qualifies,
                "window_genome_from_sample": d.get("Genome"),
                "window_strand": d.get("Strand"),
            })

    return {
        "result": result,
        "debug": debug,
        "mean_missing_row": None,
        "length_failed_row": length_failed_row,
        "overlap_window_rows": overlap_window_rows,
    }

def process_genome_jobs(args_tuple):
    """
    Top-level worker function: loads one genome's trees then processes
    all of its queries. Runs entirely inside the worker process.
    """
    genome_key, mean_csv_str, genome_jobs, shared_state = args_tuple
    mean_csv = Path(mean_csv_str)

    merged_trees, total_bp_by_strand, kept_rows, kept_minmax = load_positive_windows(
        mean_csv, shared_state["coverage_threshold"], genome_key,
        merge=True, keep_debug_data=False, log_diagnostics=False,
    )
    merged_arrays = {s: build_merged_arrays(merged_trees, s) for s in ("+", "-")}
    merged_entry = {
        "merged_trees": merged_trees,
        "total_bp_by_strand": total_bp_by_strand,
        "kept_rows": kept_rows,
        "mean_csv_str": mean_csv_str,
        "kept_min": kept_minmax[0],
        "kept_max": kept_minmax[1],
        "merged_arrays": merged_arrays,
    }

    raw_trees, _, _, _ = load_positive_windows(
        mean_csv, shared_state["probability_threshold"], genome_key,
        merge=False, keep_debug_data=True, log_diagnostics=False,
    )

    local_state = dict(shared_state)
    local_state["merged_cache"] = {genome_key: merged_entry}
    local_state["raw_cache"]    = {genome_key: raw_trees}
    set_worker_state(local_state)

    return [process_one_query(job) for job in genome_jobs]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Compute strand-specific overlap between FASTA-coordinate regions and "
            "above-threshold window intervals from per-genome *_mean.csv."
        )
    )
    ap.add_argument(
        "--query-fasta",
        required=True,
        help=(
            "Multi-FASTA with ids like CP043804_723996-724438 "
            "(reverse strand: CP043804_724236-724198, i.e. raw_start > raw_end)"
        ),
    )
    ap.add_argument(
        "--mean-root",
        required=True,
        help="Root directory containing per-genome subdirs: /GENOME/GENOME_mean.csv",
    )
    ap.add_argument(
        "--mean-suffix",
        default="_mean.csv",
        help="Filename suffix for mean CSV files (default: _mean.csv)",
    )
    ap.add_argument(
        "--coverage-threshold",
        type=float,
        default=0.3,
        help="probability_mean threshold for coverage/overlap_bp calculation (default: 0.3)",
    )
    ap.add_argument(
        "--probability-threshold",
        type=float,
        default=0.0,
        help="probability_mean threshold for raw probability metric queries (default: 0.0)",
    )
    ap.add_argument(
        "--blast-db",
        default=None,
        help="Local BLAST DB prefix used to retrieve exact genome lengths for random control",
    )
    ap.add_argument(
        "--blastdbcmd",
        default="blastdbcmd",
        help="Path to blastdbcmd executable (default: blastdbcmd in PATH)",
    )
    ap.add_argument(
        "--random-n",
        type=int,
        default=1000,
        help="Random control iterations per query (default: 1000; set 0 to disable)",
    )
    ap.add_argument(
        "--seed", 
        type=int,
        default=1,
        help="NumPy RNG seed (default: 1)"
    )
    ap.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of worker processes for per-query parallelism "
             "(default: 1; use 0 for all available CPUs; Linux/fork recommended)",
    )
    ap.add_argument(
        "--debug-dir",
        required=True,
        help="Directory for troubleshooting TSV outputs"
    )
    ap.add_argument(
        "--debug-max-overlap-windows",
        type=int,
        default=50,
        help="Max overlapping raw windows to list per query in debug TSV (default: 50; 0 disables)",
    )
    ap.add_argument(
        "--plot",
        action="store_true",
        help="Write the new 3-panel probability-focused seaborn plot",
    )
    ap.add_argument(
        "--legacy-plot",
        action="store_true",
        help="Write the legacy 3-panel coverage/probability seaborn plot",
    )
    ap.add_argument(
        "--plot-min-window-overlap-frac",
        type=float,
        default=0.99,
        help=(
            "Minimum fraction of the input sequence that a BacTermFinder window must overlap "
            "to qualify for best-window probability scoring (default: 0.99)"
        ),
    )
    ap.add_argument(
        "--out-prefix", 
        required=True, 
        help="Output file prefix (TSV + optional PNG)"
    )
    ap.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity",
    )
    args = ap.parse_args()

    if args.random_n < 0:
        raise ValueError("--random-n must be >= 0")
    if args.workers < 0:
        raise ValueError("--workers must be >= 0")
    if not (0.0 <= args.plot_min_window_overlap_frac < 1.0):
        raise ValueError("--plot-min-window-overlap-frac must be >= 0 and < 1")

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s: %(message)s",
    )

    mean_root = Path(args.mean_root)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    dbg = Path(args.debug_dir)
    dbg.mkdir(parents=True, exist_ok=True)

    raw_records = list(SeqIO.parse(args.query_fasta, "fasta"))
    if not raw_records:
        raise ValueError("No sequences found in --query-fasta")

    parsed_records = []
    for rec in raw_records:
        try:
            acc_raw, region_start, region_end, query_strand = parse_fasta_coord_id(rec.id)
            parsed_records.append((rec.id, acc_raw, region_start, region_end, query_strand))
        except ValueError as exc:
            logging.warning("Skipping unparseable FASTA record %r: %s", rec.id, exc)

    n_plus = sum(1 for *_, s in parsed_records if s == "+")
    n_minus = sum(1 for *_, s in parsed_records if s == "-")
    logging.info(
        "Parsed %d/%d FASTA records successfully: %d on '+' and %d on '-'.",
        len(parsed_records), len(raw_records), n_plus, n_minus,
    )

    if parsed_records and (n_plus == 0 or n_minus == 0):
        dominant = "+" if n_minus == 0 else "-"
        logging.warning(
            "All %d records inferred as strand %r. Verify that reverse-strand entries "
            "encode raw_start > raw_end if this is unexpected.",
            len(parsed_records), dominant,
        )

    # ------------------------------------------------------------------
    # Resolve mean CSVs per accession
    # ------------------------------------------------------------------
    unique_accessions = sorted({acc_raw for _, acc_raw, _, _, _ in parsed_records})
    mean_resolution_cache = {
        acc_raw: resolve_mean_csv(acc_raw, mean_root, args.mean_suffix)
        for acc_raw in unique_accessions
    }

    # ------------------------------------------------------------------
    # Build genome index (no trees loaded yet)
    # ------------------------------------------------------------------
    genomes_to_build = {}
    for acc_raw, resolution in mean_resolution_cache.items():
        genome_key = resolution["genome_key_used"]
        mean_csv = resolution["mean_csv"]
        if genome_key is not None and mean_csv is not None and genome_key not in genomes_to_build:
            genomes_to_build[genome_key] = mean_csv

    logging.info(
        "Will stream interval-tree caches for %d unique genome(s) (one at a time).",
        len(genomes_to_build),
    )

    # ------------------------------------------------------------------
    # Precompute genome lengths once per accession, if random control enabled
    # ------------------------------------------------------------------
    length_resolution_cache = {}
    if args.random_n > 0 and args.blast_db is not None:
        blast_length_cache = {}
        logging.info("Precomputing BLAST genome lengths for %d accession(s).", len(unique_accessions))
        for acc_raw in unique_accessions:
            candidates = prefer_dot1_candidates(acc_raw)
            subj_len = None
            blast_entry_used = None
            note = ""
            for cand in candidates:
                if cand in blast_length_cache:
                    subj_len = blast_length_cache[cand]
                    blast_entry_used = cand
                    break
            if subj_len is None:
                for cand in candidates:
                    try:
                        subj_len = blastdb_length(args.blastdbcmd, args.blast_db, cand)
                        blast_length_cache[cand] = subj_len
                        blast_entry_used = cand
                        break
                    except Exception as exc:
                        logging.debug("blastdbcmd failed for %r: %s", cand, exc)
            if subj_len is None:
                note = f"blastdbcmd length lookup failed for candidates={candidates}"
            length_resolution_cache[acc_raw] = {
                "subj_len": subj_len,
                "blast_entry_used": blast_entry_used,
                "note": note,
                "candidates": candidates,
            }

    # ------------------------------------------------------------------
    #     # ------------------------------------------------------------------
    # Group jobs by genome
    # ------------------------------------------------------------------
    from collections import defaultdict

    genome_to_jobs = defaultdict(list)
    missing_jobs = []

    all_jobs = [
        (i, qid, acc_raw, region_start, region_end, query_strand)
        for i, (qid, acc_raw, region_start, region_end, query_strand) in enumerate(parsed_records)
    ]

    for job in all_jobs:
        _, _, acc_raw, _, _, _ = job
        genome_key = mean_resolution_cache[acc_raw]["genome_key_used"]
        if genome_key is None:
            missing_jobs.append(job)
        else:
            genome_to_jobs[genome_key].append(job)

    n_workers = args.workers if args.workers > 0 else get_available_cpus()

    # Shared state passed to every worker — no trees (each worker loads its own)
    shared_state = {
        "merged_cache": {},
        "raw_cache": {},
        "mean_resolution_cache": mean_resolution_cache,
        "length_resolution_cache": length_resolution_cache,
        "seed": args.seed,
        "random_n": args.random_n,
        "blast_db": args.blast_db,
        "coverage_threshold": args.coverage_threshold,
        "probability_threshold": args.probability_threshold,
        "plot_min_window_overlap_frac": args.plot_min_window_overlap_frac,
        "debug_max_overlap_windows": args.debug_max_overlap_windows,
    }

    genome_batch_args = [
        (genome_key, str(genomes_to_build[genome_key]), genome_jobs, shared_state)
        for genome_key, genome_jobs in genome_to_jobs.items()
    ]

    # ------------------------------------------------------------------
    # Run: genome-level parallelism
    # Each worker loads one genome and processes all its queries
    # ------------------------------------------------------------------
    results = []
    debug_rows = []
    mean_missing_rows = []
    overlap_window_rows = []
    length_failed_rows = []

    outputs_by_index = {}

    # Handle no-CSV queries immediately in the main process
    set_worker_state(shared_state)
    for job in missing_jobs:
        out = process_one_query(job)
        outputs_by_index[job[0]] = out

    logging.info(
        "Processing %d genome batch(es) with %d worker(s).",
        len(genome_batch_args), n_workers,
    )

    if n_workers > 1 and genome_batch_args:
        # forkserver: workers fork from a clean server process,
        # avoiding COW page-fault amplification of parent memory
        try:
            ctx = mp.get_context("forkserver")
        except ValueError:
            ctx = mp.get_context("spawn")

        with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as pool:
            for (genome_key, _, genome_jobs, _), genome_outputs in zip(
                genome_batch_args,
                pool.map(process_genome_jobs, genome_batch_args, chunksize=1),
            ):
                for job, out in zip(genome_jobs, genome_outputs):
                    outputs_by_index[job[0]] = out
    else:
        for args_tuple in genome_batch_args:
            genome_key, _, genome_jobs, _ = args_tuple
            genome_outputs = process_genome_jobs(args_tuple)
            for job, out in zip(genome_jobs, genome_outputs):
                outputs_by_index[job[0]] = out

    # Restore original FASTA order
    for i in range(len(all_jobs)):
        out = outputs_by_index[i]
        results.append(out["result"])
        debug_rows.append(out["debug"])
        if out["mean_missing_row"] is not None:
            mean_missing_rows.append(out["mean_missing_row"])
        if out["length_failed_row"] is not None:
            length_failed_rows.append(out["length_failed_row"])
        overlap_window_rows.extend(out["overlap_window_rows"])



    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    res_df = pd.DataFrame(results)
    res_path = out_prefix.with_suffix(".overlap.tsv")
    res_df.to_csv(res_path, sep="\t", index=False)

    pd.DataFrame(debug_rows).to_csv(dbg / "debug.per_query.tsv", sep="\t", index=False)
    pd.DataFrame(mean_missing_rows).to_csv(dbg / "debug.mean_missing.tsv", sep="\t", index=False)
    pd.DataFrame(length_failed_rows).to_csv(dbg / "debug.length_failed.tsv", sep="\t", index=False)

    ow_cols = [
        "qseqid", "genome_key_used", "region_start", "region_end", "query_strand",
        "window_start", "window_end_inclusive", "window_sample", "window_probability_mean",
        "window_overlap_bp", "window_overlap_frac", "window_passes_best_overlap_filter",
        "window_genome_from_sample", "window_strand",
    ]
    ow_df = pd.DataFrame(overlap_window_rows) if overlap_window_rows else pd.DataFrame(columns=ow_cols)
    ow_df.to_csv(dbg / "debug.overlapping_windows.tsv", sep="\t", index=False)

    if args.plot:
        write_new_plot(
            res_df, out_prefix,
            threshold=args.coverage_threshold,
            min_overlap_frac=args.plot_min_window_overlap_frac,
            random_n=args.random_n,
        )

    if args.legacy_plot:
        write_legacy_plot(
            res_df,
            out_prefix,
            threshold=args.coverage_threshold,
        )

    logging.info("Wrote main results: %s", res_path)
    logging.info("Wrote debug dir: %s", dbg)


if __name__ == "__main__":
    main()
