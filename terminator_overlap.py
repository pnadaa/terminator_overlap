#!/usr/bin/env python3
import argparse
import logging
import math
import random
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from intervaltree import IntervalTree


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
    - If input is AB930127      -> [AB930127.1, AB930127]
    - If input is AB930127.1    -> [AB930127.1, AB930127]
    """
    g = genome_id.strip()
    if not g:
        return []
    if "." in g:
        base = g.rsplit(".", 1)[0]
        return [f"{base}.1", base]
    else:
        return [f"{g}.1", g]


def parse_fasta_header(record_id: str):
    """
    Accept:
      CP047845_1576955-1577393
      CP047845_1576955_1577393  (fallback)
    genomename may contain underscores; we split on the last '_' and parse coords.
    """
    if "_" not in record_id:
        raise ValueError(f"FASTA id missing '_' before coordinates: {record_id}")
    genome, coord_part = record_id.rsplit("_", 1)

    if "-" in coord_part:
        a, b = coord_part.split("-", 1)
    elif "_" in coord_part:
        a, b = coord_part.split("_", 1)
    else:
        raise ValueError(f"FASTA id coordinate part must contain '-' or '_': {record_id}")

    start = int(a)
    end = int(b)
    return genome, min(start, end), max(start, end)


def parse_samplename(sample_name: str):
    """
    SampleName:
      index_genomename_upstreamend_downstreamend_strand
    where strand is + or -, and genomename may contain underscores.
    """
    parts = sample_name.split("_")
    if len(parts) < 5:
        raise ValueError(f"SampleName does not have >=5 underscore-separated fields: {sample_name}")
    strand = parts[-1]
    downstream = int(parts[-2])
    upstream = int(parts[-3])
    genome = "_".join(parts[1:-3])
    idx = parts[0]
    if strand not in {"+", "-"}:
        raise ValueError(f"Strand must be '+' or '-', got {strand} in {sample_name}")
    start = min(upstream, downstream)
    end = max(upstream, downstream)
    return idx, genome, start, end, strand


def load_predicted_intervals(mean_csv: Path, threshold: float):
    """
    Returns:
      tree: IntervalTree with merged predicted intervals (half-open)
      pred_total_bp: total bp covered by predicted intervals (after merging)
      n_rows_kept: number of rows kept after threshold filter
    """
    df = pd.read_csv(mean_csv)
    if "SampleName" not in df.columns or "probability_mean" not in df.columns:
        raise ValueError(f"{mean_csv} missing required columns SampleName/probability_mean")

    df["probability_mean"] = df["probability_mean"].astype(float)
    df = df[df["probability_mean"] >= threshold].copy()

    tree = IntervalTree()
    for sn in df["SampleName"].astype(str).tolist():
        _, _, s, e, _ = parse_samplename(sn)
        tree.addi(s, e + 1)  # store as half-open [s, e+1)

    tree.merge_overlaps(strict=False)

    pred_total_bp = sum((iv.end - iv.begin) for iv in tree)
    return tree, int(pred_total_bp), int(len(df))


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
    cmd = [blastdbcmd, "-db", blast_db, "-entry", entry, "-outfmt", "%l"]
    out = run_cmd(cmd).strip()
    try:
        return int(out)
    except ValueError:
        raise RuntimeError(f"Could not parse blastdbcmd length for entry={entry!r}: {out!r}")


def compute_percent(overlap: int, denom_bp: int):
    if denom_bp <= 0:
        return float("nan")
    return 100.0 * (overlap / denom_bp)


def main():
    ap = argparse.ArgumentParser(
        description="Compute overlap between per-genome predicted terminators (mean.csv) and FASTA regions; includes random-region control."
    )
    ap.add_argument("--query-fasta", required=True, help="Multi-FASTA with headers like >GENOME_start-end (sequence can be AUGC and include '-')")
    ap.add_argument("--mean-root", required=True, help="Root directory containing per-genome folders, e.g. output/GENOME/GENOME_mean.csv")
    ap.add_argument("--mean-suffix", default="_mean.csv", help="Suffix for mean CSV files (default: _mean.csv)")
    ap.add_argument("--threshold", type=float, default=0.3, help="Keep rows with probability_mean >= threshold (default: 0.3)")
    ap.add_argument("--blast-db", required=True, help="Local BLAST DB prefix used with blastdbcmd -db")
    ap.add_argument("--blastdbcmd", default="blastdbcmd", help="Path to blastdbcmd (default: blastdbcmd in PATH)")
    ap.add_argument("--overlap-denom", choices=["predicted", "region", "union"], default="predicted",
                    help="Denominator for percent overlap (default: predicted)")
    ap.add_argument("--random-n", type=int, default=200, help="Random control iterations per query (default: 200)")
    ap.add_argument("--seed", type=int, default=1, help="RNG seed for random control (default: 1)")
    ap.add_argument("--plot", action="store_true", help="Write a PNG plot of observed overlap vs random control")
    ap.add_argument("--out-prefix", required=True, help="Prefix for output files (TSV + optional PNG)")
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = ap.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s")
    random.seed(args.seed)

    mean_root = Path(args.mean_root)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    # Caches
    interval_cache = {}  # genome_key_used -> (tree, pred_total_bp, n_rows_kept, mean_csv_path)
    length_cache = {}    # blast_entry_used -> length

    results = []
    random_rows = []

    records = list(SeqIO.parse(args.query_fasta, "fasta"))
    if not records:
        raise ValueError("No sequences found in --query-fasta")

    for rec in records:
        qid = rec.id
        genome_raw, region_start, region_end = parse_fasta_header(qid)
        region_len = region_end - region_start + 1

        # Find mean.csv by trying .1 first, then base.
        mean_candidates = []
        genome_candidates = prefer_dot1_candidates(genome_raw)
        for g in genome_candidates:
            mean_candidates.append(mean_root / g / f"{g}{args.mean_suffix}")

        mean_csv = None
        genome_key_used = None
        for g, p in zip(genome_candidates, mean_candidates):
            if p.exists():
                mean_csv = p
                genome_key_used = g
                break

        if mean_csv is None:
            results.append({
                "qseqid": qid,
                "genome_from_fasta": genome_raw,
                "genome_key_used": None,
                "mean_csv": None,
                "region_start": region_start,
                "region_end": region_end,
                "region_len": region_len,
                "pred_rows_kept": None,
                "predicted_total_bp": None,
                "overlap_bp": 0,
                "percent_overlap": float("nan"),
                "random_mean_percent": float("nan"),
                "random_sd_percent": float("nan"),
                "note": f"mean.csv not found under {mean_root} for candidates: {', '.join(str(x) for x in mean_candidates)}"
            })
            continue

        # Load / cache predicted intervals for this genome_key_used
        if genome_key_used not in interval_cache:
            tree, pred_total_bp, n_rows_kept = load_predicted_intervals(mean_csv, args.threshold)
            interval_cache[genome_key_used] = (tree, pred_total_bp, n_rows_kept, str(mean_csv))
        tree, pred_total_bp, n_rows_kept, mean_csv_str = interval_cache[genome_key_used]

        ov = overlap_bp(tree, region_start, region_end)

        if args.overlap_denom == "predicted":
            denom = pred_total_bp
        elif args.overlap_denom == "region":
            denom = region_len
        else:  # union
            denom = pred_total_bp + region_len - ov

        percent = compute_percent(ov, denom)

        # Random control: sample intervals of same length from the genome sequence length in BLAST DB.
        # Try same .1-first logic for blastdbcmd entry.
        blast_entry_used = None
        subj_len = None
        blast_candidates = prefer_dot1_candidates(genome_raw)
        for cand in blast_candidates:
            if cand in length_cache:
                blast_entry_used = cand
                subj_len = length_cache[cand]
                break

        if subj_len is None:
            for cand in blast_candidates:
                try:
                    subj_len = blastdb_length(args.blastdbcmd, args.blast_db, cand)
                    length_cache[cand] = subj_len
                    blast_entry_used = cand
                    break
                except Exception:
                    continue

        rand_percents = []
        if subj_len is None:
            note = "blastdbcmd length lookup failed for both .1 and unversioned genome id"
        else:
            note = ""
            if args.random_n > 0 and subj_len >= region_len and region_len > 0:
                for i in range(args.random_n):
                    r_start = random.randint(1, subj_len - region_len + 1)
                    r_end = r_start + region_len - 1
                    rov = overlap_bp(tree, r_start, r_end)

                    if args.overlap_denom == "predicted":
                        rden = pred_total_bp
                    elif args.overlap_denom == "region":
                        rden = region_len
                    else:
                        rden = pred_total_bp + region_len - rov

                    rpercent = compute_percent(rov, rden)
                    rand_percents.append(rpercent)
                    random_rows.append({
                        "qseqid": qid,
                        "genome_key_used": genome_key_used,
                        "blast_entry_used": blast_entry_used,
                        "rand_start": r_start,
                        "rand_end": r_end,
                        "rand_overlap_bp": rov,
                        "rand_percent_overlap": rpercent
                    })

        if rand_percents:
            mean_r = sum(rand_percents) / len(rand_percents)
            sd_r = (sum((x - mean_r) ** 2 for x in rand_percents) / max(1, (len(rand_percents) - 1))) ** 0.5
        else:
            mean_r = float("nan")
            sd_r = float("nan")

        results.append({
            "qseqid": qid,
            "genome_from_fasta": genome_raw,
            "genome_key_used": genome_key_used,
            "mean_csv": mean_csv_str,
            "region_start": region_start,
            "region_end": region_end,
            "region_len": region_len,
            "pred_rows_kept": n_rows_kept,
            "predicted_total_bp": pred_total_bp,
            "overlap_bp": ov,
            "percent_overlap": percent,
            "random_mean_percent": mean_r,
            "random_sd_percent": sd_r,
            "blast_entry_used": blast_entry_used,
            "blast_entry_len": subj_len,
            "note": note
        })

    res_df = pd.DataFrame(results)
    res_path = out_prefix.with_suffix(".overlap.tsv")
    res_df.to_csv(res_path, sep="\t", index=False)

    rnd_path = out_prefix.with_suffix(".random.tsv")
    if random_rows:
        pd.DataFrame(random_rows).to_csv(rnd_path, sep="\t", index=False)
    else:
        pd.DataFrame(columns=[
            "qseqid","genome_key_used","blast_entry_used","rand_start","rand_end","rand_overlap_bp","rand_percent_overlap"
        ]).to_csv(rnd_path, sep="\t", index=False)

    if args.plot:
        import matplotlib.pyplot as plt

        dfp = res_df.dropna(subset=["percent_overlap"]).copy()
        fig_h = max(3.0, 0.35 * len(dfp) + 1.5)
        fig, ax = plt.subplots(figsize=(11, fig_h))

        y = list(range(len(dfp)))
        ax.barh(y, dfp["percent_overlap"].values, label="Observed", alpha=0.85)

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
        ax.set_xlabel(f"Percent overlap (denom={args.overlap_denom})")
        ax.set_title("Predicted terminators vs FASTA region overlap")
        ax.grid(axis="x", linestyle=":", alpha=0.4)
        ax.legend(loc="lower right")

        plt.tight_layout()
        plot_path = out_prefix.with_suffix(".overlap.png")
        fig.savefig(plot_path, dpi=200)
        plt.close(fig)

    logging.info("Wrote: %s", res_path)
    logging.info("Wrote: %s", rnd_path)
    if args.plot:
        logging.info("Wrote: %s", out_prefix.with_suffix(".overlap.png"))


if __name__ == "__main__":
    main()
