#!/usr/bin/env python3
import argparse
import csv
import logging
import math
import os
import random
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from intervaltree import Interval, IntervalTree


@dataclass(frozen=True)
class BlastHit:
    qseqid: str
    sseqid: str
    sstart: int
    send: int
    bitscore: float
    evalue: float

_ACCESSION_LIKE = re.compile(r"^[A-Za-z]{1,3}_?\d+(?:\.\d+)?$")

def subject_key_candidates(sseqid: str):
    """
    Produce lookup keys in the requested order:
    - If sseqid has no explicit version, try '.1' first then unversioned.
    - If sseqid has a version, try as-is first then strip version.
    Also strips common BLAST decorations like 'ref|...|' by taking the last
    pipe-delimited token that looks like an accession.
    """
    raw = sseqid.split()[0]  # drop any trailing description

    token = raw
    if "|" in raw:
        # choose a pipe field that looks accession-like, prefer with version if present
        fields = [f for f in raw.split("|") if f]
        acc_like = [f for f in fields if _ACCESSION_LIKE.match(f)]
        token = acc_like[-1] if acc_like else fields[-1]

    if "." in token:
        base = token.rsplit(".", 1)[0]
        return [token, base]
    else:
        return [f"{token}.1", token]


def run_cmd(cmd, *, text=True):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=text)
    if p.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit {p.returncode}): {' '.join(cmd)}\n"
            f"STDERR:\n{p.stderr}\nSTDOUT:\n{p.stdout}"
        )
    return p.stdout


def parse_samplename(sample_name: str):
    """
    SampleName format:
      index_genomename_upstreamend_downstreamend_strand
    where strand is +/-, and genomename may contain underscores.
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


def parse_fasta_header(record_id: str):
    """
    Accept:
      CP047845_1576955-1577393
      CP047845_1576955_1577393
    where genome may contain underscores (we split on the last '_' then parse coords).
    """
    if "_" not in record_id:
        raise ValueError(f"FASTA id missing '_' before coordinates: {record_id}")

    genome, coord_part = record_id.rsplit("_", 1)

    # coord_part could be "start-end" or (rarely) "start_end" if someone used underscores
    if "-" in coord_part:
        a, b = coord_part.split("-", 1)
    elif "_" in coord_part:
        a, b = coord_part.split("_", 1)
    else:
        raise ValueError(f"FASTA id coordinate part must contain '-' or '_': {record_id}")

    start = int(a)
    end = int(b)
    return genome, min(start, end), max(start, end)



def to_intervaltree(df: pd.DataFrame):
    """
    Build IntervalTrees keyed by genome/contig name.
    Uses half-open intervals [start, end+1) to represent inclusive coordinates.
    """
    trees = defaultdict(IntervalTree)
    for _, row in df.iterrows():
        genome = row["Genome"]
        start = int(row["Start"])
        end = int(row["End"])
        trees[genome].addi(start, end + 1, row["probability_mean"])
    for g in list(trees.keys()):
        trees[g].merge_overlaps(strict=False)
    return trees


def interval_len(iv: Interval):
    return int(iv.end - iv.begin)


def overlap_bp(tree: IntervalTree, start_inclusive: int, end_inclusive: int):
    """
    Compute total bp overlap between merged intervals in tree and a query interval [start,end] inclusive.
    tree intervals are half-open.
    """
    q0 = start_inclusive
    q1 = end_inclusive + 1
    if q1 <= q0:
        return 0
    total = 0
    for iv in tree.overlap(q0, q1):
        a0 = max(iv.begin, q0)
        a1 = min(iv.end, q1)
        if a1 > a0:
            total += (a1 - a0)
    return int(total)


def total_bp(tree: IntervalTree):
    return int(sum(interval_len(iv) for iv in tree))


def pick_best_hit(hits):
    if not hits:
        return None
    # Prefer higher bitscore, then lower evalue, then longer span.
    def key(h: BlastHit):
        span = abs(h.send - h.sstart) + 1
        return (h.bitscore, -math.log10(h.evalue + 1e-300), span)

    return max(hits, key=key)


def blast_all_queries(query_fasta: Path, blast_db: str, blastn: str, *, evalue: float,
                      threads: int, task: str, max_target_seqs: int):
    outfmt = "6 qseqid sseqid qstart qend sstart send evalue bitscore"
    cmd = [
        blastn,
        "-query", str(query_fasta),
        "-db", blast_db,
        "-outfmt", outfmt,
        "-evalue", str(evalue),
        "-num_threads", str(threads),
        "-task", task,
        "-max_target_seqs", str(max_target_seqs),
        "-max_hsps", "1",
    ]
    tsv = run_cmd(cmd)
    per_q = defaultdict(list)
    for line in tsv.splitlines():
        if not line.strip():
            continue
        qseqid, sseqid, qstart, qend, sstart, send, evalue_s, bitscore_s = line.split("\t")
        hit = BlastHit(
            qseqid=qseqid,
            sseqid=sseqid,
            sstart=int(sstart),
            send=int(send),
            bitscore=float(bitscore_s),
            evalue=float(evalue_s),
        )
        per_q[qseqid].append(hit)
    return per_q


def blastdb_seq_length(sseqid: str, blast_db: str, blastdbcmd: str):
    # blastdbcmd custom outfmt supports %l for sequence length
    cmd = [blastdbcmd, "-db", blast_db, "-entry", sseqid, "-outfmt", "%l"]
    out = run_cmd(cmd).strip()
    try:
        return int(out)
    except ValueError:
        raise RuntimeError(f"Could not parse blastdbcmd length for {sseqid}: {out!r}")


def compute_percent(overlap: int, denom_bp: int):
    if denom_bp <= 0:
        return float("nan")
    return 100.0 * (overlap / denom_bp)


def main():
    ap = argparse.ArgumentParser(
        description="Filter predicted terminators, BLAST query sequences, and compute overlap vs BLAST-hit genomic intervals."
    )
    ap.add_argument("--pred-csv", required=True, help="CSV with header SampleName,probability_ENAC,probability_PS2,probability_NCP,probability_binary,probability_mean")
    ap.add_argument("--query-fasta", required=True, help="Multi-FASTA: >genomename_start_end then sequence in AUGC (U will be converted to T for BLAST)")
    ap.add_argument("--blast-db", required=True, help="Local BLAST DB prefix (as used with -db)")
    ap.add_argument("--threshold", type=float, default=0.3, help="Keep rows with probability_mean >= threshold (default: 0.3)")
    ap.add_argument("--out-prefix", required=True, help="Prefix for output files (TSV + optional PNG)")
    ap.add_argument("--plot", action="store_true", help="Create a plot (PNG) of observed overlap vs random control")
    ap.add_argument("--random-n", type=int, default=200, help="Random control iterations per query (default: 200)")
    ap.add_argument("--seed", type=int, default=1, help="RNG seed for random control (default: 1)")
    ap.add_argument("--overlap-denom", choices=["predicted", "blast", "union"], default="predicted",
                    help="Denominator for reported percent overlap (default: predicted)")
    ap.add_argument("--blastn", default="blastn", help="Path to blastn (default: blastn in PATH)")
    ap.add_argument("--blastdbcmd", default="blastdbcmd", help="Path to blastdbcmd (default: blastdbcmd in PATH)")
    ap.add_argument("--evalue", type=float, default=1e-5, help="BLAST e-value cutoff (default: 1e-5)")
    ap.add_argument("--threads", type=int, default=4, help="BLAST threads (default: 4)")
    ap.add_argument("--task", default="blastn", help="BLAST task, e.g. blastn or blastn-short (default: blastn)")
    ap.add_argument("--max-target-seqs", type=int, default=50, help="BLAST -max_target_seqs (default: 50)")
    ap.add_argument("--subject-id-regex", default=None,
                    help="Optional regex with one capture group to map BLAST sseqid -> genome key used in CSV. Example: '^([^|]+)\\|'")
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = ap.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s")
    random.seed(args.seed)

    pred = pd.read_csv(args.pred_csv)
    required = ["SampleName", "probability_mean"]
    missing = [c for c in required if c not in pred.columns]
    if missing:
        raise ValueError(f"Missing required columns in pred-csv: {missing}")

    pred = pred[pred["probability_mean"].astype(float) >= args.threshold].copy()

    parsed_rows = []
    for sn, pm in zip(pred["SampleName"].astype(str), pred["probability_mean"].astype(float)):
        idx, genome, start, end, strand = parse_samplename(sn)
        parsed_rows.append((sn, idx, genome, start, end, strand, float(pm)))

    pred2 = pd.DataFrame(parsed_rows, columns=["SampleName", "Index", "Genome", "Start", "End", "Strand", "probability_mean"])
    trees = to_intervaltree(pred2)

    # Write all queries to a temp fasta with U->T conversion for blastn
    query_records = list(SeqIO.parse(args.query_fasta, "fasta"))
    if not query_records:
        raise ValueError("No sequences found in query-fasta")

    with tempfile.TemporaryDirectory() as td:
        qfa = Path(td) / "queries_for_blast.fa"
        with qfa.open("w") as fh:
            for rec in query_records:
                seq = str(rec.seq).upper().replace("-", "").replace("U", "T")
                fh.write(f">{rec.id}\n{seq}\n")

        per_q_hits = blast_all_queries(
            qfa,
            args.blast_db,
            args.blastn,
            evalue=args.evalue,
            threads=args.threads,
            task=args.task,
            max_target_seqs=args.max_target_seqs,
        )

    subj_map_re = re.compile(args.subject_id_regex) if args.subject_id_regex else None

    results = []
    random_rows = []

    subj_len_cache = {}

    for rec in query_records:
        qid = rec.id
        q_genome, q_start, q_end = parse_fasta_header(qid)
        hits = per_q_hits.get(qid, [])
        best = pick_best_hit(hits)

        if best is None:
            results.append({
                "qseqid": qid,
                "query_genome": q_genome,
                "query_start": q_start,
                "query_end": q_end,
                "blast_sseqid": None,
                "blast_start": None,
                "blast_end": None,
                "blast_len": None,
                "predicted_total_bp_in_subject": None,
                "overlap_bp": 0,
                "percent_overlap": float("nan"),
                "random_mean_percent": float("nan"),
                "random_sd_percent": float("nan"),
                "note": "No BLAST hits passing filters (or BLAST produced no output lines)."
            })
            continue

        sseqid = best.sseqid
        s0 = min(best.sstart, best.send)
        s1 = max(best.sstart, best.send)
        blast_len = s1 - s0 + 1

        # Map sseqid -> genome key used for predicted terminators
        if subj_map_re:
            m = subj_map_re.search(sseqid)
            if not m:
                raise ValueError(f"--subject-id-regex did not match sseqid={sseqid!r}")
            subject_candidates = [m.group(1)]
        else:
            subject_candidates = subject_key_candidates(sseqid)

        tree = None
        subject_key_used = None
        for k in subject_candidates:
            if k in trees:
                tree = trees[k]
                subject_key_used = k
                break

        if tree is None:
            tree = IntervalTree()
            subject_key_used = subject_candidates[0]  # for reporting only
        pred_total = total_bp(tree)

        ov = overlap_bp(tree, s0, s1)

        if args.overlap_denom == "predicted":
            denom = pred_total
        elif args.overlap_denom == "blast":
            denom = blast_len
        else:  # union (Jaccard-like denominator)
            denom = pred_total + blast_len - ov

        percent = compute_percent(ov, denom)

        # Random control: same subject sequence, random interval of same length as BLAST interval
        if sseqid not in subj_len_cache:
            subj_len_cache[sseqid] = blastdb_seq_length(sseqid, args.blast_db, args.blastdbcmd)
        subj_len = subj_len_cache[sseqid]

        rand_percents = []
        if args.random_n > 0 and subj_len >= blast_len:
            for i in range(args.random_n):
                r_start = random.randint(1, subj_len - blast_len + 1)
                r_end = r_start + blast_len - 1
                rov = overlap_bp(tree, r_start, r_end)
                if args.overlap_denom == "predicted":
                    rden = pred_total
                elif args.overlap_denom == "blast":
                    rden = blast_len
                else:
                    rden = pred_total + blast_len - rov
                rpercent = compute_percent(rov, rden)
                rand_percents.append(rpercent)
                random_rows.append({
                    "qseqid": qid,
                    "blast_sseqid": sseqid,
                    "rand_start": r_start,
                    "rand_end": r_end,
                    "rand_overlap_bp": rov,
                    "rand_percent_overlap": rpercent,
                })

        if rand_percents:
            mean_r = sum(rand_percents) / len(rand_percents)
            sd_r = (sum((x - mean_r) ** 2 for x in rand_percents) / max(1, (len(rand_percents) - 1))) ** 0.5
        else:
            mean_r = float("nan")
            sd_r = float("nan")

        results.append({
            "qseqid": qid,
            "query_genome": q_genome,
            "query_start": q_start,
            "query_end": q_end,
            "blast_sseqid": sseqid,
            "blast_start": s0,
            "blast_end": s1,
            "blast_len": blast_len,
            "predicted_total_bp_in_subject": pred_total,
            "overlap_bp": ov,
            "percent_overlap": percent,
            "random_mean_percent": mean_r,
            "random_sd_percent": sd_r,
            "note": "",
            "subject_key_used": subject_key_used
        })

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    res_df = pd.DataFrame(results)
    res_path = out_prefix.with_suffix(".overlap.tsv")
    res_df.to_csv(res_path, sep="\t", index=False)

    rnd_path = out_prefix.with_suffix(".random.tsv")
    if random_rows:
        pd.DataFrame(random_rows).to_csv(rnd_path, sep="\t", index=False)
    else:
        # still write an empty file with header for reproducibility
        pd.DataFrame(columns=["qseqid","blast_sseqid","rand_start","rand_end","rand_overlap_bp","rand_percent_overlap"]).to_csv(
            rnd_path, sep="\t", index=False
        )

    if args.plot:
        import matplotlib.pyplot as plt

        dfp = res_df.dropna(subset=["percent_overlap"]).copy()
        dfp["label"] = dfp["qseqid"]

        fig_h = max(3.0, 0.35 * len(dfp) + 1.5)
        fig, ax = plt.subplots(figsize=(11, fig_h))

        y = list(range(len(dfp)))
        ax.barh(y, dfp["percent_overlap"].values, label="Observed", alpha=0.85)

        # plot random mean as points with ±1 SD
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
        ax.set_yticklabels(dfp["label"].values)
        ax.set_xlabel(f"Percent overlap (denom={args.overlap_denom})")
        ax.set_title("Predicted terminators vs BLAST-hit interval overlap")
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
