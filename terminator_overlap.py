#!/usr/bin/env python3
import argparse
import logging
import random
import re
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from intervaltree import Interval, IntervalTree


# For debug: accessions like CP043804, CP043804.1, AB930127.1, NZ_CP012345.1 etc.
ACC_RE = re.compile(r"^(?:[A-Z]{1,3}_)?[A-Z]{1,4}\d+(?:\.\d+)?$")


def run_cmd(cmd, *, text=True):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=text)
    if p.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit {p.returncode}): {' '.join(cmd)}\n"
            f"STDERR:\n{p.stderr}\nSTDOUT:\n{p.stdout}"
        )
    return p.stdout


def prefer_dot1_candidates(genome_id: str):
    """
    Always try '.1' first, then unversioned.
    - CP043804      -> [CP043804.1, CP043804]
    - CP043804.1    -> [CP043804.1, CP043804]
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
    FASTA header format:
      {accession}_{start}-{end}
    Example:
      CP043804_723996-724438

    accession may contain dots (e.g. CP043804.1), but not underscores in your example.
    """
    if "_" not in record_id:
        raise ValueError(f"FASTA id missing '_' separator: {record_id}")
    acc, coord_part = record_id.split("_", 1)
    if "-" not in coord_part:
        raise ValueError(f"FASTA id missing '-' in coords: {record_id}")
    a, b = coord_part.split("-", 1)
    start = int(a)
    end = int(b)
    return acc, min(start, end), max(start, end)


def parse_window_samplename(sample_name: str):
    """
    Decode window coordinates encoded in SampleName.

    Supports either:
      window_genome_low_high
      window_genome_low_high_strand

    Examples:
      26_CP054550.1_125_226
      1_AB930127.1_0_101_+
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


def load_positive_windows(mean_csv: Path, threshold: float, *, merge=True, keep_debug_data=False):
    """
    Load mean_csv and build an IntervalTree of windows with probability_mean >= threshold.

    If merge=True, intervals are merged into a union coverage set (best for overlap_bp).
    If keep_debug_data=True, interval data stores SampleName and probability_mean (only meaningful when merge=False).
    """
    df = pd.read_csv(mean_csv)
    if "SampleName" not in df.columns or "probability_mean" not in df.columns:
        raise ValueError(f"{mean_csv} missing required columns SampleName/probability_mean")

    df["probability_mean"] = df["probability_mean"].astype(float)
    kept = df[df["probability_mean"] >= threshold].copy()

    tree = IntervalTree()
    kept_min = None
    kept_max = None

    for sn, pm in zip(kept["SampleName"].astype(str), kept["probability_mean"].astype(float)):
        _, genome_in_sn, s, e, strand = parse_window_samplename(sn)
        kept_min = s if kept_min is None else min(kept_min, s)
        kept_max = e if kept_max is None else max(kept_max, e)
        data = None
        if keep_debug_data:
            data = {"SampleName": sn, "probability_mean": float(pm), "Genome": genome_in_sn, "Strand": strand}
        tree.addi(s, e + 1, data)  # half-open

    if merge:
        tree.merge_overlaps(strict=False, data_reducer=lambda a, b: None)

    total_bp = int(sum(iv.end - iv.begin for iv in tree))
    return tree, total_bp, int(len(kept)), (kept_min, kept_max)


def overlap_bp(tree: IntervalTree, start_inclusive: int, end_inclusive: int):
    q0 = int(start_inclusive)
    q1 = int(end_inclusive) + 1
    if q1 <= q0:
        return 0
    total = 0
    for iv in tree.overlap(q0, q1):
        a0 = max(iv.begin, q0)
        a1 = min(iv.end, q1)
        if a1 > a0:
            total += (a1 - a0)
    return int(total)


def blastdb_length(blastdbcmd: str, blast_db: str, entry: str):
    # %l is sequence length in blastdbcmd custom outfmt. [web:23]
    cmd = [blastdbcmd, "-db", blast_db, "-entry", entry, "-outfmt", "%l"]
    out = run_cmd(cmd).strip()
    try:
        return int(out)
    except ValueError:
        raise RuntimeError(f"Could not parse blastdbcmd length for entry={entry!r}: {out!r}")


def main():
    ap = argparse.ArgumentParser(
        description="Compute overlap between FASTA-coordinate regions and above-threshold window intervals from per-genome *_mean.csv (no BLAST mapping)."
    )
    ap.add_argument("--query-fasta", required=True, help="Multi-FASTA with ids like CP043804_723996-724438")
    ap.add_argument("--mean-root", required=True, help="Root dir: output/GENOME/GENOME_mean.csv")
    ap.add_argument("--mean-suffix", default="_mean.csv", help="Suffix for mean files (default: _mean.csv)")
    ap.add_argument("--threshold", type=float, default=0.3, help="probability_mean threshold (default: 0.3)")

    # Random control (optional)
    ap.add_argument("--blast-db", default=None, help="Local BLAST DB prefix (only used to get genome lengths for random control)")
    ap.add_argument("--blastdbcmd", default="blastdbcmd", help="Path to blastdbcmd (default: blastdbcmd in PATH)")
    ap.add_argument("--random-n", type=int, default=200, help="Random control iterations per query (default: 200; set 0 to disable)")
    ap.add_argument("--seed", type=int, default=1, help="RNG seed (default: 1)")

    # Debugging outputs
    ap.add_argument("--debug-dir", required=True, help="Directory to write troubleshooting outputs")
    ap.add_argument("--debug-max-overlap-windows", type=int, default=50,
                    help="Max overlapping raw windows to list per query (default: 50; 0 disables)")

    ap.add_argument("--plot", action="store_true", help="Write a PNG plot of observed vs random control")
    ap.add_argument("--out-prefix", required=True, help="Prefix for output files (TSV + optional PNG)")
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = ap.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s")
    random.seed(args.seed)

    mean_root = Path(args.mean_root)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    dbg = Path(args.debug_dir)
    dbg.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(args.query_fasta, "fasta"))
    if not records:
        raise ValueError("No sequences found in --query-fasta")

    # Caches
    merged_cache = {}  # genome_key_used -> (merged_tree, total_bp, kept_rows, mean_csv_str, kept_minmax)
    raw_cache = {}     # genome_key_used -> raw_tree (unmerged, with data)
    length_cache = {}  # blastdbcmd entry -> length

    results = []
    debug_rows = []
    mean_missing_rows = []
    overlap_window_rows = []
    length_failed_rows = []

    for rec in records:
        qid = rec.id
        acc_raw, region_start, region_end = parse_fasta_coord_id(qid)
        region_len = region_end - region_start + 1

        # Locate mean.csv: try accession.1 first, then accession
        genome_key_used = None
        mean_csv = None
        mean_tried = []
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
                "qseqid": qid,
                "accession_from_fasta": acc_raw,
                "region_start": region_start,
                "region_end": region_end,
                "region_len": region_len,
                "genome_key_used": None,
                "mean_csv": None,
                "windows_rows_kept": None,
                "windows_total_bp": None,
                "overlap_bp": 0,
                "percent_region_covered": float("nan"),
                "random_mean_percent": float("nan"),
                "random_sd_percent": float("nan"),
                "blast_entry_used": None,
                "blast_entry_len": None,
                "note": note,
            })
            mean_missing_rows.append({"qseqid": qid, "accession": acc_raw, "mean_paths_tried": ";".join(mean_tried)})
            debug_rows.append({"qseqid": qid, "accession": acc_raw, "mean_paths_tried": ";".join(mean_tried), "note": note})
            continue

        # Load/cached windows
        if genome_key_used not in merged_cache:
            merged_tree, windows_total_bp, kept_rows, kept_minmax = load_positive_windows(
                mean_csv, args.threshold, merge=True, keep_debug_data=False
            )
            merged_cache[genome_key_used] = (merged_tree, windows_total_bp, kept_rows, str(mean_csv), kept_minmax)

            raw_tree, _, _, _ = load_positive_windows(
                mean_csv, args.threshold, merge=False, keep_debug_data=True
            )
            raw_cache[genome_key_used] = raw_tree

        merged_tree, windows_total_bp, kept_rows, mean_csv_str, (kept_min, kept_max) = merged_cache[genome_key_used]

        ov = overlap_bp(merged_tree, region_start, region_end)
        percent_region_covered = 100.0 * (ov / region_len) if region_len > 0 else float("nan")

        # Random control (optional)
        subj_len = None
        blast_entry_used = None
        note = ""

        if args.random_n > 0:
            if args.blast_db is None:
                note = "Random control disabled (no --blast-db provided)."
            else:
                # Try .1 first, then base, consistent with your preference
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
                        except Exception as e:
                            continue

                if subj_len is None:
                    note = f"blastdbcmd length lookup failed for candidates={candidates}"
                    length_failed_rows.append({"qseqid": qid, "accession": acc_raw, "candidates": ";".join(candidates)})
                    rand_percents = []
                else:
                    rand_percents = []
                    if subj_len >= region_len and region_len > 0:
                        for _ in range(args.random_n):
                            r_start = random.randint(1, subj_len - region_len + 1)
                            r_end = r_start + region_len - 1
                            rov = overlap_bp(merged_tree, r_start, r_end)
                            rand_percents.append(100.0 * (rov / region_len))
                    else:
                        note = f"Genome length {subj_len} < region_len {region_len}; random control skipped."
        else:
            rand_percents = []

        if rand_percents:
            mean_r = sum(rand_percents) / len(rand_percents)
            sd_r = (sum((x - mean_r) ** 2 for x in rand_percents) / max(1, (len(rand_percents) - 1))) ** 0.5
        else:
            mean_r = float("nan")
            sd_r = float("nan")

        results.append({
            "qseqid": qid,
            "accession_from_fasta": acc_raw,
            "region_start": region_start,
            "region_end": region_end,
            "region_len": region_len,
            "genome_key_used": genome_key_used,
            "mean_csv": mean_csv_str,
            "windows_rows_kept": kept_rows,
            "windows_total_bp": windows_total_bp,
            "overlap_bp": ov,
            "percent_region_covered": percent_region_covered,
            "random_mean_percent": mean_r,
            "random_sd_percent": sd_r,
            "blast_entry_used": blast_entry_used,
            "blast_entry_len": subj_len,
            "kept_window_min_coord": kept_min,
            "kept_window_max_coord": kept_max,
            "note": note,
        })

        debug_rows.append({
            "qseqid": qid,
            "accession": acc_raw,
            "region_start": region_start,
            "region_end": region_end,
            "region_len": region_len,
            "genome_key_used": genome_key_used,
            "mean_csv_used": mean_csv_str,
            "mean_paths_tried": ";".join(mean_tried),
            "kept_rows": kept_rows,
            "kept_min_coord": kept_min,
            "kept_max_coord": kept_max,
            "overlap_bp": ov,
            "percent_region_covered": percent_region_covered,
            "blast_entry_used": blast_entry_used,
            "blast_entry_len": subj_len,
            "note": note,
        })

        # Optional: list overlapping raw windows for troubleshooting
        if args.debug_max_overlap_windows and args.debug_max_overlap_windows > 0:
            raw_tree = raw_cache[genome_key_used]
            overlaps = list(raw_tree.overlap(region_start, region_end + 1))
            overlaps.sort(key=lambda iv: (iv.begin, iv.end))
            for iv in overlaps[: args.debug_max_overlap_windows]:
                d = iv.data or {}
                ov_bp = overlap_bp(IntervalTree([Interval(iv.begin, iv.end)]), region_start, region_end)
                overlap_window_rows.append({
                    "qseqid": qid,
                    "genome_key_used": genome_key_used,
                    "region_start": region_start,
                    "region_end": region_end,
                    "window_start": iv.begin,
                    "window_end_inclusive": iv.end - 1,
                    "window_sample": d.get("SampleName"),
                    "window_probability_mean": d.get("probability_mean"),
                    "window_overlap_bp": ov_bp,
                    "window_genome_from_sample": d.get("Genome"),
                    "window_strand": d.get("Strand"),
                })

    # Write outputs
    res_df = pd.DataFrame(results)
    res_path = out_prefix.with_suffix(".overlap.tsv")
    res_df.to_csv(res_path, sep="\t", index=False)

    pd.DataFrame(debug_rows).to_csv(dbg / "debug.per_query.tsv", sep="\t", index=False)
    pd.DataFrame(mean_missing_rows).to_csv(dbg / "debug.mean_missing.tsv", sep="\t", index=False)
    pd.DataFrame(length_failed_rows).to_csv(dbg / "debug.length_failed.tsv", sep="\t", index=False)

    if overlap_window_rows:
        pd.DataFrame(overlap_window_rows).to_csv(dbg / "debug.overlapping_windows.tsv", sep="\t", index=False)
    else:
        pd.DataFrame(columns=[
            "qseqid","genome_key_used","region_start","region_end","window_start","window_end_inclusive",
            "window_sample","window_probability_mean","window_overlap_bp","window_genome_from_sample","window_strand"
        ]).to_csv(dbg / "debug.overlapping_windows.tsv", sep="\t", index=False)

    if args.plot:
        import matplotlib.pyplot as plt
        dfp = res_df.dropna(subset=["percent_region_covered"]).copy()
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
        ax.set_title("Region coverage by predicted-terminator windows (FASTA coords)")
        ax.grid(axis="x", linestyle=":", alpha=0.4)
        ax.legend(loc="lower right")
        plt.tight_layout()
        plot_path = out_prefix.with_suffix(".overlap.png")
        fig.savefig(plot_path, dpi=200)
        plt.close(fig)

    logging.info("Wrote main: %s", res_path)
    logging.info("Wrote debug dir: %s", str(dbg))


if __name__ == "__main__":
    main()
