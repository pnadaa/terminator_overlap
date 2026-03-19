"""
BacTermFinder Coverage Analysis Script (strand-aware)

Quantifies coverage of bacterial stem-loop regions by BacTermFinder terminator windows
and performs statistical comparison against a random genomic background.

Improvements over previous iteration:
  - Per-query strand inference from FASTA coordinate order (raw_start > raw_end => '-')
  - Per-strand IntervalTrees; unstranded windows inserted into both trees (no silent drops)
  - blastdbcmd genome-length lookup for accurate random-control placement
  - Extensive debug TSV outputs (per_query, mean_missing, length_failed, overlapping_windows)
  - Z-score and empirical p-value statistics
  - Per-query probability metrics (n_windows, max_probability_mean, mean_probability_mean)
  - 3-panel seaborn visualisation (observed vs random, coverage vs probability, p-value dist)
  - n_controls default 1000; NumPy Generator RNG for reproducibility
"""

import argparse
import logging
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from intervaltree import Interval, IntervalTree

ACC_RE = re.compile(r"^(?:[A-Z]{1,3}_)?[A-Z]{1,4}\d+(?:\.\d+)?$")


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def run_cmd(cmd, *, text=True):
    """Run a subprocess command and return stdout, raising on non-zero exit."""
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=text)
    if p.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit {p.returncode}): {' '.join(cmd)}\n"
            f"STDERR:\n{p.stderr}\nSTDOUT:\n{p.stdout}"
        )
    return p.stdout


def normalise_accession(acc: str) -> str:
    """Strip a trailing version suffix of the form '.<n>' from acc."""
    return re.sub(r"\.\d+$", "", acc)


def prefer_dot1_candidates(genome_id: str):
    """
    Return candidate accession strings to try, always preferring '.1' first.
      CP043804   -> [CP043804.1, CP043804]
      CP043804.1 -> [CP043804.1, CP043804]
    """
    g = genome_id.strip()
    if not g:
        return []
    if "." in g:
        base = g.rsplit(".", 1)[0]
        return [f"{base}.1", base]
    return [f"{g}.1", g]


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_fasta_coord_id(record_id: str):
    """
    Parse a FASTA header of the form {accession}_{start}-{end}.

    Strand is inferred from raw coordinate order:
      raw_start <= raw_end  ->  '+' (forward)
      raw_start >  raw_end  ->  '-' (reverse)

    Returns: (acc, normalised_start, normalised_end, strand)
    where normalised_start <= normalised_end always.
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
        raw_start, raw_end = int(b), int(a)
    else:
        strand = "+"
    return acc, min(raw_start, raw_end), max(raw_start, raw_end), strand


def parse_window_samplename(sample_name: str):
    """
    Decode window coordinates from a SampleName field.

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
        high   = int(parts[-2])
        low    = int(parts[-3])
        genome = "_".join(parts[1:-3])
        window = parts[0]
        strand = maybe_strand
    else:
        high   = int(parts[-1])
        low    = int(parts[-2])
        genome = "_".join(parts[1:-2])
        window = parts[0]
        strand = None
    start = min(low, high)
    end   = max(low, high)
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

    genome_label  is the expected accession (used for cross-checking genome IDs
                  embedded in each SampleName and for per-genome diagnostic logging).
    log_diagnostics=False suppresses warnings on secondary/debug loads of the same CSV.

    Returns:
      trees       : dict {"+": IntervalTree, "-": IntervalTree}
      total_bp    : dict {"+": int, "-": int}
      kept_rows   : int  rows passing threshold (strand-agnostic count)
      kept_minmax : (int | None, int | None)  global coordinate extent
    """
    df = pd.read_csv(mean_csv)
    if "SampleName" not in df.columns or "probability_mean" not in df.columns:
        raise ValueError(f"{mean_csv} missing required columns SampleName/probability_mean")

    df["probability_mean"] = df["probability_mean"].astype(float)
    kept = df[df["probability_mean"] >= threshold].copy()

    trees: dict = {"+": IntervalTree(), "-": IntervalTree()}
    kept_min = kept_max = None
    n_total = len(kept)
    n_stranded = 0
    n_genome_mismatch = 0
    norm_label = normalise_accession(genome_label)

    for sn, pm in zip(kept["SampleName"].astype(str), kept["probability_mean"].astype(float)):
        _, genome_in_sn, s, e, strand = parse_window_samplename(sn)

        # Fix 3: Cross-check genome ID embedded in SampleName vs expected accession
        if normalise_accession(genome_in_sn) != norm_label:
            n_genome_mismatch += 1

        kept_min = s if kept_min is None else min(kept_min, s)
        kept_max = e if kept_max is None else max(kept_max, e)

        data = None
        if keep_debug_data:
            data = {
                "SampleName": sn,
                "probability_mean": float(pm),
                "Genome": genome_in_sn,
                "Strand": strand,
            }

        if strand is not None:
            n_stranded += 1
        # Strand-unaware windows inserted into both trees so they are never silently dropped
        for ts in (["+", "-"] if strand is None else [strand]):
            trees[ts].addi(s, e + 1, data)  # half-open interval

    # Fix 2: Per-genome strand-information diagnostics
    if log_diagnostics:
        if n_total == 0:
            logging.warning("%s: No windows passed the threshold — trees are empty.", genome_label)
        elif n_stranded == 0:
            logging.warning(
                "%s: 0/%d windows carry strand information. "
                "Strand-specific overlap queries will be identical to strand-naive queries "
                "because both trees contain the same windows. "
                "Verify that BacTermFinder output for this genome encodes strand in SampleName.",
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

        # Fix 3: Genome ID cross-check diagnostic
        if n_genome_mismatch > 0:
            logging.warning(
                "%s: %d/%d windows have a genome ID in SampleName that does not match "
                "the expected accession '%s' (after version-stripping). "
                "These windows are still included but may indicate a CSV/accession mismatch.",
                genome_label, n_genome_mismatch, n_total, genome_label,
            )

    if merge:
        for tree in trees.values():
            tree.merge_overlaps(strict=False, data_reducer=lambda a, b: None)

    total_bp = {s: int(sum(iv.end - iv.begin for iv in t)) for s, t in trees.items()}
    return trees, total_bp, int(n_total), (kept_min, kept_max)


# ---------------------------------------------------------------------------
# Overlap and probability metrics
# ---------------------------------------------------------------------------

def overlap_bp(
    trees_by_strand: dict,
    query_strand: str,
    start_inclusive: int,
    end_inclusive: int,
) -> int:
    """Compute overlap in bp between [start_inclusive, end_inclusive] and the tree for query_strand."""
    tree = trees_by_strand.get(query_strand)
    if tree is None:
        return 0
    q0, q1 = int(start_inclusive), int(end_inclusive) + 1
    if q1 <= q0:
        return 0
    total = 0
    for iv in tree.overlap(q0, q1):
        a0 = max(iv.begin, q0)
        a1 = min(iv.end,   q1)
        if a1 > a0:
            total += a1 - a0
    return int(total)


def compute_probability_metrics(
    raw_trees_by_strand: dict,
    query_strand: str,
    start_inclusive: int,
    end_inclusive: int,
) -> tuple:
    """
    From the raw (unmerged, keep_debug_data=True) per-strand IntervalTree, extract
    probability_mean values for windows overlapping [start_inclusive, end_inclusive]
    on query_strand.  All intervals already passed the threshold filter.

    Returns: (n_windows, max_probability_mean, mean_probability_mean)
    """
    tree = raw_trees_by_strand.get(query_strand)
    if tree is None:
        return 0, float("nan"), float("nan")
    q0, q1 = int(start_inclusive), int(end_inclusive) + 1
    overlaps = list(tree.overlap(q0, q1))
    if not overlaps:
        return 0, float("nan"), float("nan")
    probs = [
        iv.data["probability_mean"]
        for iv in overlaps
        if iv.data and "probability_mean" in iv.data
    ]
    if not probs:
        return len(overlaps), float("nan"), float("nan")
    return len(probs), float(max(probs)), float(sum(probs) / len(probs))


# ---------------------------------------------------------------------------
# BLAST DB genome-length helper
# ---------------------------------------------------------------------------

def blastdb_length(blastdbcmd: str, blast_db: str, entry: str) -> int:
    """Return the sequence length for a BLAST DB entry via blastdbcmd."""
    cmd = [blastdbcmd, "-db", blast_db, "-entry", entry, "-outfmt", "%l"]
    out = run_cmd(cmd).strip()
    try:
        return int(out)
    except ValueError:
        raise RuntimeError(
            f"Could not parse blastdbcmd length for entry={entry!r}: {out!r}"
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    """Main CLI entry point."""
    ap = argparse.ArgumentParser(
        description=(
            "Compute strand-specific overlap between FASTA-coordinate regions and "
            "above-threshold window intervals from per-genome *_mean.csv."
        )
    )
    ap.add_argument(
        "--query-fasta", required=True,
        help=(
            "Multi-FASTA with ids like CP043804_723996-724438 "
            "(reverse strand: CP043804_724236-724198, i.e. raw_start > raw_end)"
        ),
    )
    ap.add_argument(
        "--mean-root", required=True,
        help="Root directory containing per-genome subdirs: <root>/GENOME/GENOME_mean.csv",
    )
    ap.add_argument(
        "--mean-suffix", default="_mean.csv",
        help="Filename suffix for mean CSV files (default: _mean.csv)",
    )
    ap.add_argument(
        "--threshold", type=float, default=0.3,
        help="probability_mean threshold (default: 0.3)",
    )
    ap.add_argument(
        "--blast-db", default=None,
        help="Local BLAST DB prefix used to retrieve exact genome lengths for random control",
    )
    ap.add_argument(
        "--blastdbcmd", default="blastdbcmd",
        help="Path to blastdbcmd executable (default: blastdbcmd in PATH)",
    )
    ap.add_argument(
        "--random-n", type=int, default=1000,
        help="Random control iterations per query (default: 1000; set 0 to disable)",
    )
    ap.add_argument("--seed", type=int, default=1, help="NumPy RNG seed (default: 1)")
    ap.add_argument("--debug-dir", required=True, help="Directory for troubleshooting TSV outputs")
    ap.add_argument(
        "--debug-max-overlap-windows", type=int, default=50,
        help="Max overlapping raw windows to list per query in debug TSV (default: 50; 0 disables)",
    )
    ap.add_argument(
        "--plot", action="store_true",
        help="Write a 3-panel PNG plot (requires matplotlib + seaborn)",
    )
    ap.add_argument("--out-prefix", required=True, help="Output file prefix (TSV + optional PNG)")
    ap.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    args = ap.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s: %(message)s",
    )

    # Fix 4: NumPy Generator for reproducibility and vectorised integer draws
    rng = np.random.default_rng(args.seed)

    mean_root  = Path(args.mean_root)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    dbg = Path(args.debug_dir)
    dbg.mkdir(parents=True, exist_ok=True)

    raw_records = list(SeqIO.parse(args.query_fasta, "fasta"))
    if not raw_records:
        raise ValueError("No sequences found in --query-fasta")

    # -----------------------------------------------------------------------
    # Fix 1: Pre-pass — parse all FASTA IDs, collect strand statistics early
    # -----------------------------------------------------------------------
    parsed_records = []
    for rec in raw_records:
        try:
            acc_raw, region_start, region_end, query_strand = parse_fasta_coord_id(rec.id)
            parsed_records.append((rec, acc_raw, region_start, region_end, query_strand))
        except ValueError as exc:
            logging.warning("Skipping unparseable FASTA record %r: %s", rec.id, exc)

    n_plus  = sum(1 for *_, s in parsed_records if s == "+")
    n_minus = sum(1 for *_, s in parsed_records if s == "-")
    logging.info(
        "Parsed %d/%d FASTA records successfully: %d on '+' (raw start <= end), "
        "%d on '-' (raw start > end).",
        len(parsed_records), len(raw_records), n_plus, n_minus,
    )
    if parsed_records and (n_plus == 0 or n_minus == 0):
        dominant = "+" if n_minus == 0 else "-"
        logging.warning(
            "All %d records inferred as strand '%s'. "
            "If your FASTA coordinates are always written low→high regardless of strand, "
            "this inference will always yield '+'. "
            "Verify that reversed-strand entries genuinely encode raw_start > raw_end.",
            len(parsed_records), dominant,
        )

    # Caches keyed by genome_key_used
    merged_cache:  dict = {}   # genome_key -> (merged_trees, total_bp, kept_rows, csv_str, minmax)
    raw_cache:     dict = {}   # genome_key -> raw_trees (unmerged, keep_debug_data=True)
    length_cache:  dict = {}   # accession  -> genome length (int)

    results             = []
    debug_rows          = []
    mean_missing_rows   = []
    overlap_window_rows = []
    length_failed_rows  = []

    for rec, acc_raw, region_start, region_end, query_strand in parsed_records:
        qid        = rec.id
        region_len = region_end - region_start + 1

        # -------------------------------------------------------------------
        # Locate mean.csv — prefer .1 accession variant
        # -------------------------------------------------------------------
        genome_key_used = None
        mean_csv        = None
        mean_tried      = []
        for g in prefer_dot1_candidates(acc_raw):
            # First try the plain, then the gzipped suffix mean_csv
            for candidate_suffix in [args.mean_suffix, args.mean_suffix + ".gz"]:
                p = mean_root / g / f"{g}{candidate_suffix}"
                mean_tried.append(str(p))
                if p.exists():
                    genome_key_used = g
                    mean_csv        = p
                    break
            if mean_csv is not None:
                break

        if mean_csv is None:
            note = f"mean.csv not found for accession={acc_raw!r}"
            _nan = float("nan")
            results.append({
                "qseqid":                    qid,
                "accession_from_fasta":      acc_raw,
                "region_start":              region_start,
                "region_end":                region_end,
                "query_strand":              query_strand,
                "region_len":                region_len,
                "genome_key_used":           None,
                "mean_csv":                  None,
                "windows_rows_kept":         None,
                "windows_total_bp_this_strand": None,
                "n_windows_above_threshold": None,
                "max_probability_mean":      _nan,
                "mean_probability_mean":     _nan,
                "overlap_bp":                0,
                "percent_region_covered":    _nan,
                "random_mean_percent":       _nan,
                "random_sd_percent":         _nan,
                "z_score":                   _nan,
                "empirical_pvalue":          _nan,
                "blast_entry_used":          None,
                "blast_entry_len":           None,
                "note":                      note,
            })
            mean_missing_rows.append({
                "qseqid":            qid,
                "accession":         acc_raw,
                "mean_paths_tried":  ";".join(mean_tried),
            })
            debug_rows.append({
                "qseqid":           qid,
                "accession":        acc_raw,
                "mean_paths_tried": ";".join(mean_tried),
                "note":             note,
            })
            continue

        # -------------------------------------------------------------------
        # Load / use cached per-strand IntervalTrees
        # -------------------------------------------------------------------
        if genome_key_used not in merged_cache:
            merged_trees, total_bp_by_strand, kept_rows, kept_minmax = load_positive_windows(
                mean_csv, args.threshold, genome_key_used,
                merge=True, keep_debug_data=False, log_diagnostics=True,
            )
            merged_cache[genome_key_used] = (
                merged_trees, total_bp_by_strand, kept_rows, str(mean_csv), kept_minmax
            )
            # Secondary load: unmerged + debug data for probability metrics and window listing
            raw_trees, _, _, _ = load_positive_windows(
                mean_csv, args.threshold, genome_key_used,
                merge=False, keep_debug_data=True, log_diagnostics=False,
            )
            raw_cache[genome_key_used] = raw_trees

        merged_trees, total_bp_by_strand, kept_rows, mean_csv_str, (kept_min, kept_max) = (
            merged_cache[genome_key_used]
        )
        windows_total_bp   = total_bp_by_strand.get(query_strand, 0)
        raw_trees_for_genome = raw_cache[genome_key_used]

        # -------------------------------------------------------------------
        # Strand-specific overlap (bp and %)
        # -------------------------------------------------------------------
        ov                   = overlap_bp(merged_trees, query_strand, region_start, region_end)
        percent_region_covered = (
            100.0 * ov / region_len if region_len > 0 else float("nan")
        )

        # -------------------------------------------------------------------
        # Per-query probability metrics from raw unmerged tree
        # -------------------------------------------------------------------
        n_windows, max_prob, mean_prob = compute_probability_metrics(
            raw_trees_for_genome, query_strand, region_start, region_end
        )

        # -------------------------------------------------------------------
        # Random control — Fix 4: vectorised NumPy integer draw
        # -------------------------------------------------------------------
        subj_len:       int | None   = None
        blast_entry_used: str | None = None
        note          = ""
        rand_percents: list = []

        if args.random_n > 0:
            if args.blast_db is None:
                note = "Random control disabled (no --blast-db provided)."
            else:
                candidates = prefer_dot1_candidates(acc_raw)
                # Check length cache first
                for cand in candidates:
                    if cand in length_cache:
                        subj_len         = length_cache[cand]
                        blast_entry_used = cand
                        break
                # Query blastdbcmd if not cached
                if subj_len is None:
                    for cand in candidates:
                        try:
                            subj_len = blastdb_length(args.blastdbcmd, args.blast_db, cand)
                            length_cache[cand] = subj_len
                            blast_entry_used = cand
                            break
                        except Exception as exc:
                            logging.debug("blastdbcmd failed for %r: %s", cand, exc)
                            continue
                if subj_len is None:
                    note = f"blastdbcmd length lookup failed for candidates={candidates}"
                    length_failed_rows.append({
                        "qseqid":     qid,
                        "accession":  acc_raw,
                        "candidates": ";".join(candidates),
                    })
                elif subj_len >= region_len > 0:
                    max_start = subj_len - region_len + 1
                    r_starts  = rng.integers(1, max_start + 1, size=args.random_n)
                    for r_start in r_starts:
                        r_end = int(r_start) + region_len - 1
                        rov   = overlap_bp(merged_trees, query_strand, int(r_start), r_end)
                        rand_percents.append(100.0 * rov / region_len)
                else:
                    note = (
                        f"Genome length {subj_len} < region_len {region_len}; "
                        "random control skipped."
                    )

        # -------------------------------------------------------------------
        # Summarise random control + z-score + empirical p-value
        # -------------------------------------------------------------------
        if rand_percents:
            arr    = np.array(rand_percents, dtype=float)
            mean_r = float(np.mean(arr))
            sd_r   = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
            z_score = (
                float((percent_region_covered - mean_r) / sd_r)
                if sd_r > 0.0 else float("nan")
            )
            count_ge        = sum(1 for c in rand_percents if c >= percent_region_covered)
            empirical_pvalue = float(count_ge / len(rand_percents))
        else:
            mean_r = sd_r = z_score = empirical_pvalue = float("nan")

        # -------------------------------------------------------------------
        # Collect results
        # -------------------------------------------------------------------
        results.append({
            "qseqid":                       qid,
            "accession_from_fasta":         acc_raw,
            "region_start":                 region_start,
            "region_end":                   region_end,
            "query_strand":                 query_strand,
            "region_len":                   region_len,
            "genome_key_used":              genome_key_used,
            "mean_csv":                     mean_csv_str,
            "windows_rows_kept":            kept_rows,
            "windows_total_bp_this_strand": windows_total_bp,
            "n_windows_above_threshold":    n_windows,
            "max_probability_mean":         max_prob,
            "mean_probability_mean":        mean_prob,
            "overlap_bp":                   ov,
            "percent_region_covered":       percent_region_covered,
            "random_mean_percent":          mean_r,
            "random_sd_percent":            sd_r,
            "z_score":                      z_score,
            "empirical_pvalue":             empirical_pvalue,
            "blast_entry_used":             blast_entry_used,
            "blast_entry_len":              subj_len,
            "kept_window_min_coord":        kept_min,
            "kept_window_max_coord":        kept_max,
            "note":                         note,
        })

        debug_rows.append({
            "qseqid":                       qid,
            "accession":                    acc_raw,
            "region_start":                 region_start,
            "region_end":                   region_end,
            "query_strand":                 query_strand,
            "region_len":                   region_len,
            "genome_key_used":              genome_key_used,
            "mean_csv_used":                mean_csv_str,
            "mean_paths_tried":             ";".join(mean_tried),
            "kept_rows":                    kept_rows,
            "kept_min_coord":               kept_min,
            "kept_max_coord":               kept_max,
            "windows_total_bp_this_strand": windows_total_bp,
            "n_windows_above_threshold":    n_windows,
            "max_probability_mean":         max_prob,
            "mean_probability_mean":        mean_prob,
            "overlap_bp":                   ov,
            "percent_region_covered":       percent_region_covered,
            "z_score":                      z_score,
            "empirical_pvalue":             empirical_pvalue,
            "blast_entry_used":             blast_entry_used,
            "blast_entry_len":              subj_len,
            "note":                         note,
        })

        # Optional: per-query raw overlapping window listing for troubleshooting
        if args.debug_max_overlap_windows and args.debug_max_overlap_windows > 0:
            raw_tree = raw_trees_for_genome.get(query_strand, IntervalTree())
            overlaps_sorted = sorted(
                raw_tree.overlap(region_start, region_end + 1),
                key=lambda iv: (iv.begin, iv.end),
            )
            for iv in overlaps_sorted[: args.debug_max_overlap_windows]:
                d = iv.data or {}
                ov_bp_single = overlap_bp(
                    {query_strand: IntervalTree([Interval(iv.begin, iv.end)])},
                    query_strand, region_start, region_end,
                )
                overlap_window_rows.append({
                    "qseqid":                  qid,
                    "genome_key_used":         genome_key_used,
                    "region_start":            region_start,
                    "region_end":              region_end,
                    "query_strand":            query_strand,
                    "window_start":            iv.begin,
                    "window_end_inclusive":    iv.end - 1,
                    "window_sample":           d.get("SampleName"),
                    "window_probability_mean": d.get("probability_mean"),
                    "window_overlap_bp":       ov_bp_single,
                    "window_genome_from_sample": d.get("Genome"),
                    "window_strand":           d.get("Strand"),
                })

    # -----------------------------------------------------------------------
    # Write main results TSV
    # -----------------------------------------------------------------------
    res_df   = pd.DataFrame(results)
    res_path = out_prefix.with_suffix(".overlap.tsv")
    res_df.to_csv(res_path, sep="\t", index=False)

    # -----------------------------------------------------------------------
    # Write debug TSVs
    # -----------------------------------------------------------------------
    pd.DataFrame(debug_rows).to_csv(
        dbg / "debug.per_query.tsv", sep="\t", index=False
    )
    pd.DataFrame(mean_missing_rows).to_csv(
        dbg / "debug.mean_missing.tsv", sep="\t", index=False
    )
    pd.DataFrame(length_failed_rows).to_csv(
        dbg / "debug.length_failed.tsv", sep="\t", index=False
    )
    _ow_cols = [
        "qseqid", "genome_key_used", "region_start", "region_end", "query_strand",
        "window_start", "window_end_inclusive", "window_sample", "window_probability_mean",
        "window_overlap_bp", "window_genome_from_sample", "window_strand",
    ]
    ow_df = (
        pd.DataFrame(overlap_window_rows)
        if overlap_window_rows
        else pd.DataFrame(columns=_ow_cols)
    )
    ow_df.to_csv(dbg / "debug.overlapping_windows.tsv", sep="\t", index=False)

    # -----------------------------------------------------------------------
    # 3-panel seaborn visualisation
    # -----------------------------------------------------------------------
    if args.plot:
        import matplotlib.pyplot as plt
        import seaborn as sns

        dfp   = res_df.dropna(subset=["percent_region_covered"]).copy()
        valid = dfp.dropna(subset=["percent_region_covered", "random_mean_percent"]).copy()

        sns.set_theme(
            style="ticks",
            rc={"axes.labelsize": 10, "xtick.labelsize": 10, "ytick.labelsize": 10},
        )
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))

        if valid.empty:
            logging.warning("No valid data for plotting — generating empty figure.")
            fig.text(0.5, 0.5, "No valid data to plot", ha="center", va="center")
        else:
            # ------------------------------------------------------------------
            # Panel A: Observed vs. Random Coverage (boxplot + stripplot + lines)
            # ------------------------------------------------------------------
            axA = axes[0]
            melted = valid[
                ["qseqid", "percent_region_covered", "random_mean_percent"]
            ].melt(
                id_vars="qseqid",
                value_vars=["percent_region_covered", "random_mean_percent"],
                var_name="Group", value_name="Coverage",
            )
            melted["Group"] = melted["Group"].replace({
                "percent_region_covered": "Stem Loops",
                "random_mean_percent":    "Random Control",
            })
            sns.boxplot(
                data=melted, x="Group", y="Coverage",
                color="white", ax=axA, showfliers=False, width=0.4,
            )
            sns.stripplot(
                data=melted, x="Group", y="Coverage",
                color="black", alpha=0.5, ax=axA, jitter=False,
            )
            for _, row in valid.iterrows():
                axA.plot(
                    [0, 1],
                    [row["percent_region_covered"], row["random_mean_percent"]],
                    color="grey", alpha=0.3, linewidth=0.5,
                )
            median_obs  = valid["percent_region_covered"].median()
            median_ctrl = valid["random_mean_percent"].median()
            axA.text(
                0.05, median_obs,
                f"Median: {median_obs:.1f}%", va="bottom", ha="left",
                color="red", fontsize=9,
            )
            axA.text(
                1.05, median_ctrl,
                f"Median: {median_ctrl:.1f}%", va="bottom", ha="left",
                color="red", fontsize=9,
            )
            axA.set_ylabel("BacTermFinder Coverage (%)")
            axA.set_xlabel("")
            axA.set_title("Observed vs. Random Coverage (strand-specific)")

            # ------------------------------------------------------------------
            # Panel B: Coverage vs. Max Terminator Probability
            # ------------------------------------------------------------------
            axB = axes[1]
            scatter_data          = valid.dropna(subset=["max_probability_mean"]).copy()
            scatter_data["Covered"] = scatter_data["percent_region_covered"] > 0
            palette = {}
            if True  in scatter_data["Covered"].values: palette[True]  = "steelblue"
            if False in scatter_data["Covered"].values: palette[False] = "tomato"
            sns.scatterplot(
                data=scatter_data,
                x="max_probability_mean", y="percent_region_covered",
                hue="Covered", ax=axB, palette=palette, alpha=0.7,
            )
            axB.axvline(
                args.threshold, color="black", linestyle="--",
                label=f"Threshold ({args.threshold})",
            )
            axB.set_xlabel("Max Terminator Probability (above threshold)")
            axB.set_ylabel("BacTermFinder Coverage (%)")
            axB.set_title("Coverage vs. Terminator Probability")
            handles, labels = axB.get_legend_handles_labels()
            if handles:
                axB.legend(fontsize=8)

            # ------------------------------------------------------------------
            # Panel C: Empirical P-value Distribution
            # ------------------------------------------------------------------
            axC  = axes[2]
            pvals = valid["empirical_pvalue"].dropna()
            if not pvals.empty:
                sns.histplot(pvals, bins=20, binrange=(0, 1), ax=axC, color="skyblue")
            axC.axvline(0.05, color="black", linestyle="--", label="p = 0.05")
            axC.set_xlabel("Empirical P-value (vs. random)")
            axC.set_ylabel("Count")
            axC.set_title("Empirical P-value Distribution")
            axC.legend()

        plt.tight_layout()
        plot_path = out_prefix.with_suffix(".overlap.png")
        fig.savefig(plot_path, dpi=300)
        plt.close(fig)
        logging.info("Wrote plot: %s", plot_path)

    logging.info("Wrote main results: %s", res_path)
    logging.info("Wrote debug dir:    %s", str(dbg))


if __name__ == "__main__":
    main()
