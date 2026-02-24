"""
BacTermFinder Coverage Analysis Script
Quantifies coverage of bacterial stem-loop regions by BacTermFinder terminator windows
and performs statistical comparison against a random genomic background.
"""

import argparse
import logging
import re
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from intervaltree import Interval, IntervalTree


ACC_RE = re.compile(r"^(?:[A-Z]{1,3}_)?[A-Z]{1,4}\d+(?:\.\d+)?$")


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
    """Strip a trailing version suffix of the form '.<integer>' from acc."""
    return re.sub(r'\.\d+$', '', acc)


def prefer_dot1_candidates(genome_id: str):
    """
    Return candidate accession strings to try, always preferring '.1' first.

    - CP043804   -> [CP043804.1, CP043804]
    - CP043804.1 -> [CP043804.1, CP043804]
    """
    g = genome_id.strip()
    if not g:
        return []
    if "." in g:
        base = g.rsplit(".", 1)[0]
        return [f"{base}.1", base]
    return [f"{g}.1", g]


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
    strand = "-" if raw_start > raw_end else "+"
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

    genome_label is the expected accession (used for cross-checking genome IDs
    embedded in each SampleName and for per-genome diagnostic logging).
    Set log_diagnostics=False on secondary/debug loads of the same CSV to
    avoid duplicate warnings.

    Returns:
      trees        : dict {"+": IntervalTree, "-": IntervalTree}
      total_bp     : dict {"+": int, "-": int}
      kept_rows    : int   rows passing threshold (strand-agnostic count)
      kept_minmax  : (int | None, int | None)  global coordinate extent
    """
    df = pd.read_csv(mean_csv)
    if "SampleName" not in df.columns or "probability_mean" not in df.columns:
        raise ValueError(f"{mean_csv} missing required columns SampleName/probability_mean")

    df["probability_mean"] = df["probability_mean"].astype(float)
    kept = df[df["probability_mean"] >= threshold].copy()

    trees: dict[str, IntervalTree] = {"+": IntervalTree(), "-": IntervalTree()}
    kept_min: int | None = None
    kept_max: int | None = None

    n_total = len(kept)
    n_stranded = 0
    n_genome_mismatch = 0
    norm_label = normalise_accession(genome_label)

    for sn, pm in zip(kept["SampleName"].astype(str), kept["probability_mean"].astype(float)):
        _, genome_in_sn, s, e, strand = parse_window_samplename(sn)

        # Fix 3: Cross-check the genome ID embedded in SampleName vs expected accession
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
        target_strands = ["+", "-"] if strand is None else [strand]
        for ts in target_strands:
            trees[ts].addi(s, e + 1, data)  # half-open interval

    # Fix 2: Per-genome strand-information diagnostics
    if log_diagnostics:
        if n_total == 0:
            logging.warning(
                "%s: No windows passed the threshold — trees are empty.", genome_label
            )
        elif n_stranded == 0:
            logging.warning(
                "%s: 0/%d windows carry strand information. "
                "Strand-specific overlap queries will be identical to strand-naive queries "
                "because both trees contain the same windows. "
                "Verify that BacTermFinder output for this genome encodes strand in SampleName.",
                genome_label, n_total,
            )
        elif n_stranded < n_total:
            n_unstranded = n_total - n_stranded
            logging.warning(
                "%s: %d/%d windows carry strand information; %d unstranded window(s) "
                "inserted into both '+' and '-' trees.",
                genome_label, n_stranded, n_total, n_unstranded,
            )
        else:
            logging.debug(
                "%s: All %d windows carry strand information.", genome_label, n_total
            )

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


def overlap_bp(
    trees_by_strand: dict,
    query_strand: str,
    start_inclusive: int,
    end_inclusive: int,
) -> int:
    """Compute the overlap in bp between [start_inclusive, end_inclusive] and the tree for query_strand."""
    tree = trees_by_strand.get(query_strand)
    if tree is None:
        return 0
    q0 = int(start_inclusive)
    q1 = int(end_inclusive) + 1
    if q1 <= q0:
        return 0
    total = 0
    for iv in tree.overlap(q0, q1):
        a0 = max(iv.begin, q0)
        a1 = min(iv.end, q1)
        if a1 > a0:
            total += a1 - a0
    return int(total)


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
        help="Multi-FASTA with ids like CP043804_723996-724438 (reverse: CP043804_724236-724198)",
    )
    ap.add_argument("--mean-root", required=True, help="Root dir: output/GENOME/GENOME_mean.csv")
    ap.add_argument("--mean-suffix", default="_mean.csv", help="Suffix for mean files (default: _mean.csv)")
    ap.add_argument("--threshold", type=float, default=0.3, help="probability_mean threshold (default: 0.3)")

    ap.add_argument(
        "--blast-db", default=None,
        help="Local BLAST DB prefix (only used to get genome lengths for random control)",
    )
    ap.add_argument("--blastdbcmd", default="blastdbcmd", help="Path to blastdbcmd (default: in PATH)")
    ap.add_argument(
        "--random-n", type=int, default=200,
        help="Random control iterations per query (default: 200; set 0 to disable)",
    )
    ap.add_argument("--seed", type=int, default=1, help="RNG seed (default: 1)")

    ap.add_argument("--debug-dir", required=True, help="Directory to write troubleshooting outputs")
    ap.add_argument(
        "--debug-max-overlap-windows", type=int, default=50,
        help="Max overlapping raw windows to list per query (default: 50; 0 disables)",
    )

    ap.add_argument("--plot", action="store_true", help="Write a PNG plot of observed vs random control")
    ap.add_argument("--out-prefix", required=True, help="Prefix for output files (TSV + optional PNG)")
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = ap.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s")

    # Fix 4: Replace Python random module with NumPy Generator for reproducibility and vectorisation
    rng = np.random.default_rng(args.seed)

    mean_root = Path(args.mean_root)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    dbg = Path(args.debug_dir)
    dbg.mkdir(parents=True, exist_ok=True)

    raw_records = list(SeqIO.parse(args.query_fasta, "fasta"))
    if not raw_records:
        raise ValueError("No sequences found in --query-fasta")

    # Fix 1: Pre-pass to parse all FASTA IDs, collect strand counts, and warn early
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
    if len(parsed_records) > 0 and (n_plus == 0 or n_minus == 0):
        dominant = "+" if n_minus == 0 else "-"
        logging.warning(
            "All %d records inferred as strand '%s'. "
            "If your FASTA coordinates are always written low→high regardless of strand, "
            "this inference will always yield '+'. "
            "Verify that reversed-strand entries genuinely encode raw_start > raw_end.",
            len(parsed_records), dominant,
        )

    # Caches keyed by genome_key_used
    merged_cache: dict[str, tuple] = {}
    raw_cache:    dict[str, dict]  = {}
    length_cache: dict[str, int]   = {}

    results            = []
    debug_rows         = []
    mean_missing_rows  = []
    overlap_window_rows = []
    length_failed_rows = []

    for rec, acc_raw, region_start, region_end, query_strand in parsed_records:
        qid = rec.id
        region_len = region_end - region_start + 1

        # Locate mean.csv: try .1 candidate first
        genome_key_used = None
        mean_csv        = None
        mean_tried      = []
        for g in prefer_dot1_candidates(acc_raw):
            p = mean_root / g / f"{g}{args.mean_suffix}"
            mean_tried.append(str(p))
            if p.exists():
                genome_key_used = g
                mean_csv = p
                break

        if mean_csv is None:
            note = f"mean.csv not found for accession={acc_raw!r}"
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
                "overlap_bp":                0,
                "percent_region_covered":    float("nan"),
                "random_mean_percent":       float("nan"),
                "random_sd_percent":         float("nan"),
                "blast_entry_used":          None,
                "blast_entry_len":           None,
                "note":                      note,
            })
            mean_missing_rows.append({
                "qseqid":           qid,
                "accession":        acc_raw,
                "mean_paths_tried": ";".join(mean_tried),
            })
            debug_rows.append({
                "qseqid":           qid,
                "accession":        acc_raw,
                "mean_paths_tried": ";".join(mean_tried),
                "note":             note,
            })
            continue

        # Load / use cached per-strand trees
        if genome_key_used not in merged_cache:
            merged_trees, total_bp_by_strand, kept_rows, kept_minmax = load_positive_windows(
                mean_csv, args.threshold, genome_key_used,
                merge=True, keep_debug_data=False, log_diagnostics=True,
            )
            merged_cache[genome_key_used] = (
                merged_trees, total_bp_by_strand, kept_rows, str(mean_csv), kept_minmax
            )

            # Secondary load for debug window listing — suppress duplicate warnings
            raw_trees, _, _, _ = load_positive_windows(
                mean_csv, args.threshold, genome_key_used,
                merge=False, keep_debug_data=True, log_diagnostics=False,
            )
            raw_cache[genome_key_used] = raw_trees

        merged_trees, total_bp_by_strand, kept_rows, mean_csv_str, (kept_min, kept_max) = (
            merged_cache[genome_key_used]
        )
        windows_total_bp = total_bp_by_strand.get(query_strand, 0)

        # Strand-specific overlap
        ov = overlap_bp(merged_trees, query_strand, region_start, region_end)
        percent_region_covered = 100.0 * (ov / region_len) if region_len > 0 else float("nan")

        # Random control — Fix 4: vectorised NumPy integer draw
        subj_len: int | None = None
        blast_entry_used: str | None = None
        note = ""
        rand_percents: list[float] = []

        if args.random_n > 0:
            if args.blast_db is None:
                note = "Random control disabled (no --blast-db provided)."
            else:
                candidates = prefer_dot1_candidates(acc_raw)
                for cand in candidates:
                    if cand in length_cache:
                        subj_len = length_cache[cand]
                        blast_entry_used = cand
                        break
                if subj_len is None:
                    for cand in candidates:
                        try:
                            subj_len = blastdb_length(args.blastdbcmd, args.blast_db, cand)
                            length_cache[cand] = subj_len
                            blast_entry_used = cand
                            break
                        except Exception:
                            continue

                if subj_len is None:
                    note = f"blastdbcmd length lookup failed for candidates={candidates}"
                    length_failed_rows.append({
                        "qseqid":     qid,
                        "accession":  acc_raw,
                        "candidates": ";".join(candidates),
                    })
                elif subj_len >= region_len and region_len > 0:
                    # Draw all random starts at once; rng.integers high is exclusive
                    max_start = subj_len - region_len + 1
                    r_starts = rng.integers(1, max_start + 1, size=args.random_n)
                    for r_start in r_starts:
                        r_end = int(r_start) + region_len - 1
                        rov = overlap_bp(merged_trees, query_strand, int(r_start), r_end)
                        rand_percents.append(100.0 * (rov / region_len))
                else:
                    note = (
                        f"Genome length {subj_len} < region_len {region_len}; "
                        "random control skipped."
                    )

        if rand_percents:
            arr = np.array(rand_percents, dtype=float)
            mean_r = float(np.mean(arr))
            sd_r   = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
        else:
            mean_r = float("nan")
            sd_r   = float("nan")

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
            "overlap_bp":                   ov,
            "percent_region_covered":       percent_region_covered,
            "random_mean_percent":          mean_r,
            "random_sd_percent":            sd_r,
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
            "overlap_bp":                   ov,
            "percent_region_covered":       percent_region_covered,
            "blast_entry_used":             blast_entry_used,
            "blast_entry_len":              subj_len,
            "note":                         note,
        })

        # Optional: list raw overlapping windows per query for troubleshooting
        if args.debug_max_overlap_windows and args.debug_max_overlap_windows > 0:
            raw_trees   = raw_cache[genome_key_used]
            raw_tree    = raw_trees.get(query_strand, IntervalTree())
            overlaps    = sorted(
                raw_tree.overlap(region_start, region_end + 1),
                key=lambda iv: (iv.begin, iv.end),
            )
            for iv in overlaps[: args.debug_max_overlap_windows]:
                d = iv.data or {}
                ov_bp_single = overlap_bp(
                    {query_strand: IntervalTree([Interval(iv.begin, iv.end)])},
                    query_strand,
                    region_start,
                    region_end,
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

    # Write outputs
    res_df   = pd.DataFrame(results)
    res_path = out_prefix.with_suffix(".overlap.tsv")
    res_df.to_csv(res_path, sep="\t", index=False)

    pd.DataFrame(debug_rows).to_csv(dbg / "debug.per_query.tsv",    sep="\t", index=False)
    pd.DataFrame(mean_missing_rows).to_csv(dbg / "debug.mean_missing.tsv", sep="\t", index=False)
    pd.DataFrame(length_failed_rows).to_csv(dbg / "debug.length_failed.tsv", sep="\t", index=False)

    _ow_cols = [
        "qseqid", "genome_key_used", "region_start", "region_end", "query_strand",
        "window_start", "window_end_inclusive", "window_sample", "window_probability_mean",
        "window_overlap_bp", "window_genome_from_sample", "window_strand",
    ]
    ow_df = pd.DataFrame(overlap_window_rows) if overlap_window_rows else pd.DataFrame(columns=_ow_cols)
    ow_df.to_csv(dbg / "debug.overlapping_windows.tsv", sep="\t", index=False)

    if args.plot:
        import matplotlib.pyplot as plt
        dfp   = res_df.dropna(subset=["percent_region_covered"]).copy()
        fig_h = max(3.0, 0.35 * len(dfp) + 1.5)
        fig, ax = plt.subplots(figsize=(11, fig_h))
        y = list(range(len(dfp)))
        ax.barh(y, dfp["percent_region_covered"].values, label="Observed", alpha=0.85)
        if dfp["random_mean_percent"].notna().any():
            ax.errorbar(
                dfp["random_mean_percent"].values,
                y,
                xerr=dfp["random_sd_percent"].values,
                fmt="o",
                color="black",
                label="Random mean ± SD",
                markersize=4,
                capsize=2,
                linewidth=1,
            )
        ax.set_yticks(y)
        ax.set_yticklabels(dfp["qseqid"].values)
        ax.set_xlabel("Percent of region covered by above-threshold windows")
        ax.set_title("Region coverage by predicted-terminator windows (strand-specific)")
        ax.grid(axis="x", linestyle=":", alpha=0.4)
        ax.legend(loc="lower right")
        plt.tight_layout()
        plot_path = out_prefix.with_suffix(".overlap.png")
        fig.savefig(plot_path, dpi=200)
        plt.close(fig)

    logging.info("Wrote main results: %s", res_path)
    logging.info("Wrote debug dir:    %s", str(dbg))


if __name__ == "__main__":
    main()
