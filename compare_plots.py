#!/usr/bin/env python3
"""
violin_stats.py
───────────────
Generate annotated seaborn violin plots from two or more CSV files and run
Prism 11-style statistical tests (Welch's t / paired-t / Mann–Whitney /
Wilcoxon / one-way ANOVA / RM-ANOVA / Kruskal–Wallis / Friedman) with
Tukey, Bonferroni, Sidak, Dunnett, or Dunn post-hoc corrections.

Dependencies
────────────
    pip install pandas numpy matplotlib seaborn scipy statsmodels pingouin scikit-posthocs

Usage examples
──────────────
    # Auto test selection with interactive prompts for distribution / pairing
    python violin_stats.py group1.csv group2.csv group3.csv \
        --column "expression" --interactive

    # Fully argparse-driven
    python violin_stats.py ctrl.csv treated.csv \
        --column "fold_change" --labels Control Treated \
        --distribution lognormal --posthoc tukey \
        --title "Fold Change by Treatment" --output results.png

    # Paired, nonparametric, three groups
    python violin_stats.py pre.csv mid.csv post.csv \
        --column "score" --paired --subject-column "patient_id" \
        --distribution nonparametric --posthoc dunn

    # Skip statistics, just plot
    python violin_stats.py a.csv b.csv --column "value" --no-stats
"""

from __future__ import annotations

import argparse
import sys
import warnings
from itertools import combinations
from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="seaborn")

# ── Optional heavy dependencies ───────────────────────────────────────────────
try:
    import pingouin as pg
    HAS_PINGOUIN = True
except ImportError:
    HAS_PINGOUIN = False

try:
    import scikit_posthocs as sp
    HAS_SP = True
except ImportError:
    HAS_SP = False

try:
    import statsmodels.stats.multicomp as smc
    HAS_SM = True
except ImportError:
    HAS_SM = False

# ── Constants ─────────────────────────────────────────────────────────────────
DISTRIBUTION_CHOICES = ["normal", "lognormal", "nonparametric"]
POSTHOC_CHOICES      = ["tukey", "bonferroni", "sidak", "dunnett", "dunn"]
TEST_CHOICES         = [
    "auto", "welch", "ttest_paired",
    "mannwhitney", "wilcoxon",
    "anova", "rm_anova", "kruskal", "friedman",
]
_STAT_LABEL = {
    "welch": "t",  "ttest_paired": "t",
    "mannwhitney": "U", "wilcoxon": "W",
    "anova": "F",  "rm_anova": "F",
    "kruskal": "H", "friedman": "χ²",
}


# ══════════════════════════════════════════════════════════════════════════════
# ARGUMENT PARSER
# ══════════════════════════════════════════════════════════════════════════════

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="violin_stats",
        description=(
            "Annotated violin plots (box + whiskers, mean + median) "
            "from ≥2 CSV files, with Prism 11-style statistical tests."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "Prism 11 test mapping:\n"
            "  2 groups  │ normal      │ unpaired → Welch's t          │ paired → Paired t\n"
            "  2 groups  │ lognormal   │ unpaired → Welch's t (log)    │ paired → Ratio paired t\n"
            "  2 groups  │ nonpar.     │ unpaired → Mann–Whitney U     │ paired → Wilcoxon SR\n"
            "  ≥3 groups │ normal      │ unpaired → One-way ANOVA      │ paired → RM-ANOVA\n"
            "  ≥3 groups │ lognormal   │ unpaired → ANOVA (log)        │ paired → RM-ANOVA (log)\n"
            "  ≥3 groups │ nonpar.     │ unpaired → Kruskal–Wallis     │ paired → Friedman\n"
        ),
    )

    # ── Data ──────────────────────────────────────────────────────────────────
    io = p.add_argument_group("Data I/O")
    io.add_argument(
        "csvfiles", nargs="+", metavar="FILE",
        help=(
            "Two or more delimited data files (.csv, .tsv, .txt). "
            "Delimiter is inferred from the extension; use --sep to override."
        ),
    )
    io.add_argument(
        "-c", "--column", required=True, metavar="HEADER",
        help="Column header containing the values to plot and analyse.",
    )
    io.add_argument(
        "-l", "--labels", nargs="+", metavar="LABEL",
        help="Group labels (must match CSV count). Defaults to file stems.",
    )
    io.add_argument(
        "--subject-column", metavar="HEADER", default=None,
        help="Column identifying subjects for paired / RM tests.",
    )
    io.add_argument(
        "--sep", "--delimiter", dest="sep", default=None, metavar="CHAR",
        help=(
            "Field delimiter. Inferred automatically from the file extension "
            "(.csv → comma, .tsv/.txt → tab) when not set. "
            "Use --sep $'\\t' in bash for an explicit tab."
        ),
    )


    # ── Plot cosmetics ────────────────────────────────────────────────────────
    vis = p.add_argument_group("Plot appearance")
    vis.add_argument("--title",   default="", help="Plot title.")
    vis.add_argument("--ylabel",  default="", help="Y-axis label (defaults to --column).")
    vis.add_argument("--palette", default="Set2",  help="Seaborn colour palette.")
    vis.add_argument("--width",   type=float, default=8.0,  help="Figure width (inches).")
    vis.add_argument("--height",  type=float, default=6.0,  help="Figure height (inches).")
    vis.add_argument("--dpi",     type=int,   default=300,  help="Output DPI.")
    vis.add_argument(
        "-o", "--output", default="violin_plot.png",
        help="Output image path (.png / .pdf / .svg).",
    )
    # ── New aesthetic controls ────────────────────────────────────────────────
    vis.add_argument(
        "--style",
        choices=["paper", "talk", "poster", "dark"],
        default="paper",
        help=(
            "Overall visual theme. "
            "'paper' = minimal, publication-ready; "
            "'talk' = larger text, bolder colours; "
            "'poster' = very large text; "
            "'dark' = dark background."
        ),
    )
    vis.add_argument(
        "--show-points", action="store_true", default=False,
        help="Overlay jittered individual data points on each violin (default: off).",
    )
    vis.add_argument(
        "--no-points", dest="show_points", action="store_false",
        help="Suppress individual data-point overlay.",
    )
    vis.add_argument(
        "--point-size",  type=float, default=3.5,
        help="Diameter of individual data points in points².",
    )
    vis.add_argument(
        "--point-alpha", type=float, default=0.55,
        help="Opacity of individual data points (0–1).",
    )
    vis.add_argument(
        "--violin-alpha", type=float, default=0.72,
        help="Opacity of the violin body fill (0–1).",
    )
    vis.add_argument(
        "--show-n", action="store_true", default=True,
        help="Print sample size (n=…) under each group label (default: on).",
    )
    vis.add_argument(
        "--no-n", dest="show_n", action="store_false",
        help="Suppress sample-size labels.",
    )


    # ── Statistics ────────────────────────────────────────────────────────────
    st = p.add_argument_group("Statistical tests")
    st.add_argument(
        "--distribution", "-d",
        choices=DISTRIBUTION_CHOICES, default=None,
        help=(
            "Distribution assumption used by Prism 11 to pick the test. "
            "normal → parametric; lognormal → parametric on log-transformed data; "
            "nonparametric → rank-based."
        ),
    )
    st.add_argument(
        "--paired", action="store_true",
        help="Use paired / repeated-measures tests.",
    )
    st.add_argument(
        "--test",
        choices=TEST_CHOICES, default="auto",
        help=(
            "Override automatic test selection. "
            "'auto' uses the Prism 11 decision tree (requires --distribution)."
        ),
    )
    st.add_argument(
        "--posthoc",
        choices=POSTHOC_CHOICES, default="tukey",
        help="Post-hoc correction for ≥3 group comparisons.",
    )
    st.add_argument(
        "--alpha", type=float, default=0.05,
        help="Significance threshold α.",
    )
    st.add_argument(
        "--no-stats", action="store_true",
        help="Skip all statistical analyses.",
    )
    st.add_argument(
        "--save-stats", metavar="CSV_PATH",
        help="Save pairwise post-hoc results to a CSV file.",
    )
    st.add_argument(
        "--interactive", action="store_true",
        help="Interactively prompt for distribution / pairing / post-hoc settings.",
    )

    return p


# ══════════════════════════════════════════════════════════════════════════════
# INTERACTIVE PROMPTS
# ══════════════════════════════════════════════════════════════════════════════

def interactive_settings(args: argparse.Namespace) -> argparse.Namespace:
    """Prompt the user to fill in any unset statistical options."""
    SEP = "─" * 62
    print(f"\n{SEP}")
    print("  Statistical test configuration")
    print(SEP)

    if args.distribution is None:
        mapping = {"1": "normal", "2": "lognormal", "3": "nonparametric"}
        while True:
            ans = input(
                "\n  Distribution assumption?\n"
                "    [1] Normal        – Gaussian; use parametric tests\n"
                "    [2] Log-normal    – right-skewed; parametric on log data\n"
                "    [3] Nonparametric – no distribution assumed; rank-based\n"
                "  Choice [1/2/3]: "
            ).strip()
            if ans in mapping:
                args.distribution = mapping[ans]
                break
            print("  Please enter 1, 2, or 3.")

    if not args.paired:
        ans = input(
            "\n  Are observations paired / repeated measures? [y/N]: "
        ).strip().lower()
        args.paired = ans in ("y", "yes")

    if args.test == "auto":
        ph_map = {
            "1": "tukey", "2": "bonferroni",
            "3": "sidak",  "4": "dunnett", "5": "dunn",
        }
        ans = input(
            "\n  Post-hoc correction (for ≥3 groups)?\n"
            "    [1] Tukey HSD          (recommended; controls familywise error)\n"
            "    [2] Bonferroni         (conservative; all pairwise)\n"
            "    [3] Sidak              (slightly less conservative than Bonferroni)\n"
            "    [4] Dunnett            (comparisons vs. first group / control only)\n"
            "    [5] Dunn               (nonparametric; follows Kruskal–Wallis / Friedman)\n"
            "  Choice [1–5, default=1]: "
        ).strip()
        args.posthoc = ph_map.get(ans, "tukey")

    print(f"\n  → distribution={args.distribution}, paired={args.paired}, "
          f"posthoc={args.posthoc}")
    print(SEP + "\n")
    return args


# ── Delimiter resolution ───────────────────────────────────────────────────────
_EXT_SEP: dict[str, str] = {
    ".csv":  ",",
    ".tsv":  "\t",
    ".txt":  "\t",   # plain-text tabular exports are almost always tab-delimited
}

def _infer_sep(path: Path, explicit_sep: str | None) -> str:
    """
    Resolve the delimiter for *path* using this priority order:
      1. Explicit --sep value from the user
      2. File extension mapping (.csv → ',' / .tsv|.txt → '\t')
      3. csv.Sniffer on the first 8 KB of the file (fallback)
    """
    if explicit_sep is not None:
        return explicit_sep

    ext = path.suffix.lower()
    if ext in _EXT_SEP:
        return _EXT_SEP[ext]

    # Unknown extension — sniff the file
    import csv as _csv
    try:
        with path.open(newline="", encoding="utf-8-sig") as fh:
            sample = fh.read(8192)
        dialect = _csv.Sniffer().sniff(sample, delimiters=",\t|;")
        sep = dialect.delimiter
        print(f"[INFO] Auto-detected delimiter for '{path.name}': {sep!r}")
        return sep
    except _csv.Error:
        print(
            f"[WARN] Could not detect delimiter for '{path.name}'; "
            "defaulting to comma. Use --sep to override."
        )
        return ","


# ══════════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════

def load_data(
    paths: list[str],
    column: str,
    labels: list[str] | None,
    subject_col: str | None,
    explicit_sep: str | None = None,
) -> tuple[pd.DataFrame, list[str]]:
    """
    Read each CSV, extract *column* (and optional subject column),
    and return a long-form DataFrame with '_group' column plus
    an ordered list of group labels.
    """
    if labels and len(labels) != len(paths):
        sys.exit(
            f"[ERROR] --labels count ({len(labels)}) ≠ file count ({len(paths)})."
        )

    frames: list[pd.DataFrame] = []
    for i, fp in enumerate(paths):
        path = Path(fp)
        if not path.exists():
            sys.exit(f"[ERROR] File not found: {fp}")
        # ── Resolve delimiter ────────────────────────────────────────────────
        sep = _infer_sep(path, explicit_sep)

        try:
            df = pd.read_csv(path, sep=sep, encoding="utf-8-sig")
        except Exception as e:
            sys.exit(f"[ERROR] Could not read '{fp}' (sep={sep!r}): {e}")

        if column not in df.columns:
            sys.exit(
                f"[ERROR] Column '{column}' not found in '{fp}'.\n"
                f"        Available columns: {list(df.columns)}"
            )

        keep = [column]
        if subject_col:
            if subject_col not in df.columns:
                sys.exit(
                    f"[ERROR] Subject column '{subject_col}' not found in '{fp}'.\n"
                    f"        Available columns: {list(df.columns)}"
                )
            keep.append(subject_col)

        sub = df[keep].dropna(subset=[column]).copy()
        lbl = (labels[i] if labels else None) or path.stem
        sub["_group"] = lbl
        frames.append(sub)

    long = pd.concat(frames, ignore_index=True)
    group_order = labels if labels else [Path(fp).stem for fp in paths]
    return long, group_order


# ══════════════════════════════════════════════════════════════════════════════
# PRISM 11 TEST SELECTION TREE
# ══════════════════════════════════════════════════════════════════════════════

def select_test(
    n_groups: int,
    distribution: str,
    paired: bool,
    forced_test: str,
) -> str:
    """
    Return the canonical test key following the Prism 11 decision tree.

    Prism 11 logic
    ──────────────
    Two groups
      Normal / Lognormal + Unpaired → Welch's t-test
      Normal / Lognormal + Paired   → Paired t-test (ratio paired for lognormal)
      Nonparametric  + Unpaired     → Mann–Whitney U
      Nonparametric  + Paired       → Wilcoxon matched-pairs signed-rank

    Three or more groups
      Normal / Lognormal + Unpaired → One-way ANOVA
      Normal / Lognormal + Paired   → Repeated-measures one-way ANOVA
      Nonparametric  + Unpaired     → Kruskal–Wallis H
      Nonparametric  + Paired       → Friedman
    """
    if forced_test != "auto":
        return forced_test

    if n_groups == 2:
        if distribution == "nonparametric":
            return "wilcoxon" if paired else "mannwhitney"
        return "ttest_paired" if paired else "welch"
    else:
        if distribution == "nonparametric":
            return "friedman" if paired else "kruskal"
        return "rm_anova" if paired else "anova"


# ══════════════════════════════════════════════════════════════════════════════
# STATISTICAL TESTS
# ══════════════════════════════════════════════════════════════════════════════

def _group_vals(long: pd.DataFrame, g: str, col: str) -> np.ndarray:
    return long.loc[long["_group"] == g, col].to_numpy(dtype=float)


def run_statistics(
    long: pd.DataFrame,
    column: str,
    group_order: list[str],
    test_key: str,
    distribution: str,
    posthoc_key: str,
    alpha: float,
    paired: bool,
    subject_col: str | None,
) -> dict:
    """
    Execute the selected test and return a results dict with keys:
        test_name, statistic, p_value, posthoc_df (optional), summary
    """
    result: dict = {}
    n_groups = len(group_order)

    # ── Log-transform for lognormal ───────────────────────────────────────────
    work = long.copy()
    if distribution == "lognormal":
        neg_mask = work[column] <= 0
        if neg_mask.any():
            n_neg = int(neg_mask.sum())
            warnings.warn(
                f"[WARN] {n_neg} non-positive values excluded before log-transform."
            )
            work = work[~neg_mask]
        work["_value"] = np.log(work[column])
    else:
        work["_value"] = work[column]

    groups_data = [_group_vals(work, g, "_value") for g in group_order]

    # ── Two-group tests ───────────────────────────────────────────────────────
    if test_key == "welch":
        stat, p = stats.ttest_ind(*groups_data, equal_var=False)
        result["test_name"] = (
            "Welch's unpaired t-test"
            if distribution != "lognormal"
            else "Welch's unpaired t-test (log-transformed, compares geometric means)"
        )

    elif test_key == "ttest_paired":
        stat, p = stats.ttest_rel(*groups_data)
        result["test_name"] = (
            "Paired t-test"
            if distribution != "lognormal"
            else "Ratio paired t-test (log-transformed, compares geometric means)"
        )

    elif test_key == "mannwhitney":
        stat, p = stats.mannwhitneyu(*groups_data, alternative="two-sided")
        result["test_name"] = "Mann–Whitney U test"

    elif test_key == "wilcoxon":
        if len(groups_data[0]) != len(groups_data[1]):
            sys.exit(
                "[ERROR] Wilcoxon signed-rank requires equal-length paired groups. "
                "Check your data or use --subject-column to align pairs."
            )
        stat, p = stats.wilcoxon(*groups_data)
        result["test_name"] = "Wilcoxon matched-pairs signed-rank test"

    # ── Multi-group omnibus ───────────────────────────────────────────────────
    elif test_key == "anova":
        stat, p = stats.f_oneway(*groups_data)
        result["test_name"] = (
            "One-way ANOVA"
            if distribution != "lognormal"
            else "One-way ANOVA (log-transformed)"
        )

    elif test_key == "rm_anova":
        if not HAS_PINGOUIN:
            print(
                "[WARN] pingouin not installed — falling back to Friedman test "
                "for repeated-measures ANOVA.\n"
                "       Install pingouin: pip install pingouin"
            )
            test_key = "friedman"
        else:
            if subject_col is None:
                sys.exit(
                    "[ERROR] --subject-column is required for repeated-measures ANOVA. "
                    "Provide the column that identifies individual subjects."
                )
            rm_df = work[["_group", subject_col, "_value"]].copy()
            try:
                aov = pg.rm_anova(
                    data=rm_df, dv="_value",
                    within="_group", subject=subject_col,
                    detailed=True,
                )
                stat = float(aov.loc[aov["Source"] == "_group", "F"].iloc[0])
                p    = float(aov.loc[aov["Source"] == "_group", "p-unc"].iloc[0])
            except Exception as e:
                sys.exit(f"[ERROR] RM-ANOVA failed: {e}")
            result["test_name"] = (
                "Repeated-measures one-way ANOVA"
                if distribution != "lognormal"
                else "Repeated-measures one-way ANOVA (log-transformed)"
            )

    if test_key == "kruskal":
        stat, p = stats.kruskal(*groups_data)
        result["test_name"] = "Kruskal–Wallis H test"

    elif test_key == "friedman":
        sizes = [len(d) for d in groups_data]
        if len(set(sizes)) > 1:
            sys.exit(
                "[ERROR] Friedman test requires equal group sizes (paired design). "
                f"Group sizes found: {dict(zip(group_order, sizes))}.\n"
                "        Use --subject-column to align observations."
            )
        stat, p = stats.friedmanchisquare(*groups_data)
        result["test_name"] = "Friedman test"

    result["statistic"] = stat
    result["p_value"]   = p

    # ── Post-hoc (only when omnibus is significant, ≥3 groups) ───────────────
    if n_groups >= 3 and p < alpha:
        result["posthoc_df"] = _run_posthoc(
            work, "_value", group_order,
            posthoc_key, distribution,
            subject_col, paired,
        )

    # ── Human-readable summary ────────────────────────────────────────────────
    sig_str = "significant" if p < alpha else "not significant"
    stat_sym = _STAT_LABEL.get(test_key, "stat")
    result["summary"] = (
        f"{result['test_name']}:  "
        f"{stat_sym} = {stat:.4f},  p = {p:.4e}  "
        f"({'*' * _p_stars_count(p, alpha) or 'ns'}  {sig_str} at α={alpha})"
    )
    return result


# ══════════════════════════════════════════════════════════════════════════════
# POST-HOC TESTS
# ══════════════════════════════════════════════════════════════════════════════

def _run_posthoc(
    work: pd.DataFrame,
    val_col: str,
    group_order: list[str],
    posthoc_key: str,
    distribution: str,
    subject_col: str | None,
    paired: bool,
) -> pd.DataFrame:
    """Dispatch to the appropriate post-hoc test; return symmetric p-value matrix."""
    n = len(group_order)
    mat = pd.DataFrame(
        np.ones((n, n)), index=group_order, columns=group_order
    )
    pair_list = list(combinations(group_order, 2))

    # ── Dunn's test (nonparametric post-hoc) ─────────────────────────────────
    if posthoc_key == "dunn":
        if HAS_SP:
            ph = sp.posthoc_dunn(
                work, val_col=val_col, group_col="_group",
                p_adjust="bonferroni",
            )
            return ph.loc[group_order, group_order]
        else:
            print(
                "[WARN] scikit-posthocs not installed — "
                "falling back to Bonferroni-corrected Mann–Whitney.\n"
                "       pip install scikit-posthocs"
            )
            return _pairwise_mw_bonferroni(work, val_col, group_order)

    # ── Tukey HSD (requires statsmodels) ─────────────────────────────────────
    if posthoc_key == "tukey":
        if not HAS_SM:
            print("[WARN] statsmodels not installed — falling back to Bonferroni t-tests.")
            return _pairwise_ttest(work, val_col, group_order, "bonferroni", paired)
        tukey = smc.pairwise_tukeyhsd(
            work[val_col].values, work["_group"].values, alpha=0.05
        )
        return _tukey_result_to_df(tukey, group_order)

    # ── Dunnett (comparisons vs. first group = control) ───────────────────────
    if posthoc_key == "dunnett":
        return _dunnett_test(work, val_col, group_order, paired)

    # ── Bonferroni / Sidak (pairwise t-tests with correction) ────────────────
    return _pairwise_ttest(work, val_col, group_order, posthoc_key, paired)


def _tukey_result_to_df(tukey_result, group_order: list[str]) -> pd.DataFrame:
    n = len(group_order)
    mat = pd.DataFrame(np.ones((n, n)), index=group_order, columns=group_order)
    for row in tukey_result.summary().data[1:]:
        g1, g2 = str(row[0]), str(row[1])
        p = float(row[3])
        if g1 in group_order and g2 in group_order:
            mat.loc[g1, g2] = p
            mat.loc[g2, g1] = p
    return mat


def _pairwise_ttest(
    work: pd.DataFrame,
    col: str,
    group_order: list[str],
    correction: str,
    paired: bool,
) -> pd.DataFrame:
    pairs = list(combinations(group_order, 2))
    raw_p = []
    for g1, g2 in pairs:
        a = _group_vals(work, g1, col)
        b = _group_vals(work, g2, col)
        if paired:
            _, p = stats.ttest_rel(a, b)
        else:
            _, p = stats.ttest_ind(a, b, equal_var=False)
        raw_p.append(p)

    k = len(pairs)
    if correction == "bonferroni":
        adj = [min(p * k, 1.0) for p in raw_p]
    else:  # sidak
        adj = [1.0 - (1.0 - p) ** k for p in raw_p]

    n = len(group_order)
    mat = pd.DataFrame(np.ones((n, n)), index=group_order, columns=group_order)
    for (g1, g2), p in zip(pairs, adj):
        mat.loc[g1, g2] = p
        mat.loc[g2, g1] = p
    return mat


def _pairwise_mw_bonferroni(
    work: pd.DataFrame, col: str, group_order: list[str]
) -> pd.DataFrame:
    pairs = list(combinations(group_order, 2))
    raw_p = [
        stats.mannwhitneyu(
            _group_vals(work, g1, col),
            _group_vals(work, g2, col),
            alternative="two-sided",
        )[1]
        for g1, g2 in pairs
    ]
    k   = len(pairs)
    adj = [min(p * k, 1.0) for p in raw_p]
    n   = len(group_order)
    mat = pd.DataFrame(np.ones((n, n)), index=group_order, columns=group_order)
    for (g1, g2), p in zip(pairs, adj):
        mat.loc[g1, g2] = p
        mat.loc[g2, g1] = p
    return mat


def _dunnett_test(
    work: pd.DataFrame,
    col: str,
    group_order: list[str],
    paired: bool,
) -> pd.DataFrame:
    """
    Compare each non-control group to the first group (control),
    Bonferroni-corrected for the number of comparisons.
    Uses scipy.stats.dunnett when available (scipy ≥ 1.10);
    falls back to Bonferroni-corrected Welch's t-test otherwise.
    """
    n = len(group_order)
    mat = pd.DataFrame(np.ones((n, n)), index=group_order, columns=group_order)
    ctrl = group_order[0]
    ctrl_data = _group_vals(work, ctrl, col)
    treatment_groups = group_order[1:]
    n_comp = len(treatment_groups)

    try:
        from scipy.stats import dunnett as scipy_dunnett
        treatments_data = [_group_vals(work, g, col) for g in treatment_groups]
        res = scipy_dunnett(*treatments_data, control=ctrl_data)
        for g, p in zip(treatment_groups, res.pvalue):
            mat.loc[ctrl, g] = p
            mat.loc[g, ctrl] = p
    except ImportError:
        # Fallback: Bonferroni-corrected t-test (approximate Dunnett)
        for g in treatment_groups:
            d = _group_vals(work, g, col)
            if paired:
                _, p_raw = stats.ttest_rel(ctrl_data, d)
            else:
                _, p_raw = stats.ttest_ind(ctrl_data, d, equal_var=False)
            p_adj = min(p_raw * n_comp, 1.0)
            mat.loc[ctrl, g] = p_adj
            mat.loc[g, ctrl] = p_adj
    return mat


# ══════════════════════════════════════════════════════════════════════════════
# P-VALUE ANNOTATIONS
# ══════════════════════════════════════════════════════════════════════════════

def _p_stars_count(p: float, alpha: float) -> int:
    if p < 0.0001: return 4
    if p < 0.001:  return 3
    if p < 0.01:   return 2
    if p < alpha:  return 1
    return 0


def p_to_stars(p: float, alpha: float) -> str:
    return "*" * _p_stars_count(p, alpha) or "ns"


def annotate_pairwise(
    ax: plt.Axes,
    group_order: list[str],
    posthoc_df: pd.DataFrame,
    alpha: float,
    y_top: float,
    y_range: float,
) -> None:
    """
    Draw bracket-and-star significance annotations above the violins.
    Only draws comparisons where Dunnett was used (vs. control) or
    all pairwise for Tukey / Bonferroni / Sidak / Dunn.
    """
    step   = y_range * 0.09
    level  = 0
    for (i, j) in combinations(range(len(group_order)), 2):
        g1, g2 = group_order[i], group_order[j]
        p      = posthoc_df.loc[g1, g2]
        label  = p_to_stars(p, alpha)
        top    = y_top + step * (level + 1.4)
        x1, x2 = float(i), float(j)
        bar_h  = top - step * 0.15
        ax.plot(
            [x1, x1, x2, x2],
            [bar_h - step * 0.15, bar_h, bar_h, bar_h - step * 0.15],
            lw=1.1, color="black", clip_on=False,
        )
        colour = "black" if label != "ns" else "dimgray"
        ax.text(
            (x1 + x2) / 2, bar_h + step * 0.05,
            label,
            ha="center", va="bottom", fontsize=9,
            color=colour, fontweight="bold" if label != "ns" else "normal",
        )
        level += 1

    ax.set_ylim(top=y_top + step * (level + 2.5))


# ── Per-style theme tokens ─────────────────────────────────────────────────────
_STYLE_CONFIG: dict[str, dict] = {
    "paper": dict(
        context="paper",
        font_scale=1.15,
        rc={
            "axes.facecolor":   "#FAFAFA",
            "figure.facecolor": "white",
            "axes.edgecolor":   "#BBBBBB",
            "axes.linewidth":   0.9,
            "grid.color":       "#E4E4E4",
            "grid.linewidth":   0.7,
            "xtick.color":      "#333333",
            "ytick.color":      "#333333",
            "axes.spines.top":  False,
            "axes.spines.right":False,
            "font.family":      "sans-serif",
        },
        violin_lw=1.1,
        box_lw=1.0,
        cap_lw=1.4,
        grid_axis="y",
        grid_alpha=0.55,
        bg_stripe=True,
    ),
    "talk": dict(
        context="talk",
        font_scale=1.0,
        rc={
            "axes.facecolor":   "white",
            "figure.facecolor": "white",
            "axes.edgecolor":   "#999999",
            "axes.linewidth":   1.2,
            "grid.color":       "#DDDDDD",
            "grid.linewidth":   0.8,
            "xtick.color":      "#222222",
            "ytick.color":      "#222222",
            "axes.spines.top":  False,
            "axes.spines.right":False,
            "font.family":      "sans-serif",
        },
        violin_lw=1.5,
        box_lw=1.3,
        cap_lw=1.8,
        grid_axis="y",
        grid_alpha=0.5,
        bg_stripe=False,
    ),
    "poster": dict(
        context="poster",
        font_scale=1.0,
        rc={
            "axes.facecolor":   "white",
            "figure.facecolor": "white",
            "axes.edgecolor":   "#888888",
            "axes.linewidth":   1.5,
            "grid.color":       "#DDDDDD",
            "grid.linewidth":   1.0,
            "xtick.color":      "#111111",
            "ytick.color":      "#111111",
            "axes.spines.top":  False,
            "axes.spines.right":False,
            "font.family":      "sans-serif",
        },
        violin_lw=1.8,
        box_lw=1.6,
        cap_lw=2.2,
        grid_axis="y",
        grid_alpha=0.45,
        bg_stripe=False,
    ),
    "dark": dict(
        context="paper",
        font_scale=1.15,
        rc={
            "axes.facecolor":   "#1C1C1E",
            "figure.facecolor": "#111111",
            "axes.edgecolor":   "#444444",
            "axes.linewidth":   0.8,
            "grid.color":       "#2E2E2E",
            "grid.linewidth":   0.7,
            "xtick.color":      "#CCCCCC",
            "ytick.color":      "#CCCCCC",
            "text.color":       "#EEEEEE",
            "axes.labelcolor":  "#EEEEEE",
            "axes.spines.top":  False,
            "axes.spines.right":False,
            "font.family":      "sans-serif",
        },
        violin_lw=1.1,
        box_lw=1.0,
        cap_lw=1.4,
        grid_axis="y",
        grid_alpha=0.6,
        bg_stripe=False,
    ),
}


def _apply_style(style_name: str) -> dict:
    """Apply the seaborn theme and return the style token dict."""
    cfg = _STYLE_CONFIG[style_name]
    sns.set_theme(
        context=cfg["context"],
        style="ticks",
        font_scale=cfg["font_scale"],
        rc=cfg["rc"],
    )
    return cfg



# ══════════════════════════════════════════════════════════════════════════════
# VIOLIN PLOT
# ══════════════════════════════════════════════════════════════════════════════

def make_violin_plot(
    long: pd.DataFrame,
    column: str,
    group_order: list[str],
    args: argparse.Namespace,
    stat_result: dict | None,
) -> plt.Figure:
    """
    Publication-quality violin plot featuring:
      • Semi-transparent KDE violin body
      • Alternating light column stripes (paper style)
      • Jittered individual data points (optional)
      • Elegant narrow box-and-whisker overlay
      • Mean (◆ filled) and median (◇ open) markers with numeric labels
      • n= sample-size tick sub-labels
      • Pairwise significance brackets
      • Omnibus test footer
    """
    cfg     = _apply_style(args.style)
    n_groups = len(group_order)
    palette  = sns.color_palette(args.palette, n_colors=n_groups)
    is_dark  = args.style == "dark"

    # ── Per-style colours ──────────────────────────────────────────────────────
    box_edge   = "#DDDDDD" if is_dark else "#2A2A2A"
    cap_col    = "#CCCCCC" if is_dark else "#2A2A2A"
    whisker_col= "#CCCCCC" if is_dark else "#444444"
    flier_col  = "#888888" if is_dark else "#888888"
    marker_text= "#BBBBBB" if is_dark else "#555555"
    stat_text  = "#AAAAAA" if is_dark else "#444444"
    bracket_col= "#CCCCCC" if is_dark else "#222222"

    fig, ax = plt.subplots(figsize=(args.width, args.height))

    # ── Subtle alternating column stripes (paper / talk only) ─────────────────
    if cfg.get("bg_stripe"):
        for xi in range(n_groups):
            if xi % 2 == 0:
                ax.axvspan(xi - 0.5, xi + 0.5, color="#F0F0F5", alpha=0.45, zorder=0)

    # ── Horizontal grid ───────────────────────────────────────────────────────
    ax.grid(
        True,
        axis="y",
        linestyle="--",
        linewidth=0.7,
        color=cfg["rc"]["grid.color"],
        alpha=cfg["grid_alpha"],
        zorder=0,
    )
    ax.set_axisbelow(True)

    # ── Violin body ───────────────────────────────────────────────────────────
    vp = sns.violinplot(
        data=long,
        x="_group", y=column,
        order=group_order,
        palette=args.palette,
        inner=None,
        cut=0,
        density_norm="width",
        linewidth=cfg["violin_lw"],
        saturation=0.82,
        ax=ax,
    )
    # Apply alpha to violin fills and soften the edge colour
    for i, poly in enumerate(ax.collections):
        if hasattr(poly, "get_paths") and len(poly.get_paths()) > 0:
            poly.set_alpha(args.violin_alpha)
            r, g_c, b, _ = poly.get_facecolor()[0]
            poly.set_edgecolor((r * 0.65, g_c * 0.65, b * 0.65, 0.9))

    # ── Jittered individual data points ───────────────────────────────────────
    if args.show_points:
        for xi, grp in enumerate(group_order):
            vals = _group_vals(long, grp, column)
            rng  = np.random.default_rng(seed=xi + 42)
            jitter_x = xi + rng.uniform(-0.12, 0.12, size=len(vals))
            r, g_c, b = palette[xi][:3]
            ax.scatter(
                jitter_x, vals,
                s=args.point_size ** 2 * 0.55,
                color=(r, g_c, b, args.point_alpha),
                edgecolors=(r * 0.55, g_c * 0.55, b * 0.55, 0.7),
                linewidths=0.35,
                zorder=2,
            )

    # ── Box-and-whisker overlay ───────────────────────────────────────────────
    sns.boxplot(
        data=long,
        x="_group", y=column,
        order=group_order,
        width=0.14,
        showfliers=False,          # outliers already shown by stripplot
        showcaps=True,
        boxprops    =dict(facecolor="white", edgecolor=box_edge,
                          linewidth=cfg["box_lw"], zorder=4),
        whiskerprops=dict(color=whisker_col,  linewidth=cfg["box_lw"],
                          linestyle="-", zorder=4),
        capprops    =dict(color=cap_col,      linewidth=cfg["cap_lw"], zorder=4),
        medianprops =dict(color=box_edge,     linewidth=0, zorder=5),
        ax=ax,
    )

    # ── Mean & median markers + numeric labels ────────────────────────────────
    y_all   = long[column].dropna()
    y_range = float(y_all.max() - y_all.min()) or 1.0
    lbl_x_off = 0.48

    for xi, grp in enumerate(group_order):
        vals     = _group_vals(long, grp, column)
        mean_v   = float(np.mean(vals))
        median_v = float(np.median(vals))
        r, g_c, b = palette[xi][:3]

        # Median — white/dark diamond
        ax.scatter(
            xi, median_v,
            marker="D", s=38,
            color="white" if not is_dark else "#2A2A2A",
            edgecolors=box_edge, linewidths=1.05,
            zorder=6,
        )
        # Mean — colour-filled circle
        ax.scatter(
            xi, mean_v,
            marker="o", s=38,
            color=(r, g_c, b, 1.0),
            edgecolors=box_edge, linewidths=1.05,
            zorder=6,
        )
        # Numeric labels
        for y_val, prefix in ((median_v, "Mdn"), (mean_v, "M")):
            ax.text(
                xi + lbl_x_off, y_val,
                f"{prefix} = {y_val:.3g}",
                fontsize=6.8, va="center", ha="left",
                color=marker_text, clip_on=False,
            )

    # ── n= sub-labels ─────────────────────────────────────────────────────────
    if args.show_n:
        current_xlabels = [grp for grp in group_order]
        new_labels = []
        for grp in group_order:
            n = int((_group_vals(long, grp, column)).shape[0])
            new_labels.append(f"{grp}\n$n={n}$")
        ax.set_xticks(range(n_groups))
        ax.set_xticklabels(new_labels, fontsize=9)

    # ── Significance brackets ─────────────────────────────────────────────────
    if stat_result and "posthoc_df" in stat_result:
        annotate_pairwise(
            ax, group_order, stat_result["posthoc_df"],
            args.alpha,
            y_top   = float(y_all.max()),
            y_range = y_range,
        )
        # Re-colour brackets for dark mode
        for line in ax.lines[-50:]:
            line.set_color(bracket_col)

    # ── Axis cosmetics ────────────────────────────────────────────────────────
    ax.set_xlabel("Group", fontsize=11, labelpad=8,
                  color=cfg["rc"].get("axes.labelcolor", "#333333"))
    ax.set_ylabel(
        args.ylabel or column, fontsize=11, labelpad=8,
        color=cfg["rc"].get("axes.labelcolor", "#333333"),
    )
    ax.set_title(
        args.title or f"Distribution of {column}",
        fontsize=13, fontweight="semibold", pad=12,
        color=cfg["rc"].get("text.color", "#111111"),
    )
    ax.tick_params(axis="both", which="both", length=3, width=0.8)
    sns.despine(ax=ax, trim=False)

    # ── Legend ────────────────────────────────────────────────────────────────
    legend_handles = [
        mlines.Line2D([], [], marker="D", linestyle="None",
                      markerfacecolor="white" if not is_dark else "#2A2A2A",
                      markeredgecolor=box_edge, markersize=7, label="Median"),
        mlines.Line2D([], [], marker="o", linestyle="None",
                      markerfacecolor="gray", markeredgecolor=box_edge,
                      markersize=7, label="Mean"),
    ]
    if args.show_points:
        legend_handles.append(
            mlines.Line2D([], [], marker="o", linestyle="None",
                          markerfacecolor=(0.5, 0.5, 0.5, 0.5),
                          markeredgecolor="none", markersize=5,
                          label="Individual observations")
        )
    if stat_result and "posthoc_df" in stat_result:
        legend_handles.append(
            mpatches.Patch(color="none",
                           label="ns  p≥0.05 | *  p<0.05 | **  p<0.01 | "
                                 "***  p<0.001 | ****  p<0.0001")
        )
    leg = ax.legend(
        handles=legend_handles,
        fontsize=7.2, framealpha=0.85,
        edgecolor="#CCCCCC" if not is_dark else "#444444",
        facecolor="white"  if not is_dark else "#1C1C1E",
        labelcolor=cfg["rc"].get("text.color", "#333333"),
        loc="upper right",
    )

    # ── Omnibus test footer ───────────────────────────────────────────────────
    if stat_result:
        fig.text(
            0.5, 0.005,
            stat_result["summary"],
            ha="center", va="bottom",
            fontsize=7.2, color=stat_text, style="italic",
            bbox=dict(
                boxstyle="round,pad=0.4",
                facecolor="#1A1A1A" if is_dark else "white",
                edgecolor="#444444" if is_dark else "#CCCCCC",
                alpha=0.88,
            ),
        )
        plt.subplots_adjust(bottom=0.11)

    plt.tight_layout(rect=[0, 0.06, 1, 1])
    return fig



# ══════════════════════════════════════════════════════════════════════════════
# CONSOLE OUTPUT & CSV EXPORT
# ══════════════════════════════════════════════════════════════════════════════

def print_descriptives(long: pd.DataFrame, column: str, group_order: list[str]) -> None:
    SEP = "─" * 70
    print(f"\n{SEP}")
    print(f"  Descriptive statistics  ({column})")
    print(SEP)
    header = f"{'Group':>22}  {'N':>5}  {'Mean':>10}  {'Median':>10}  {'SD':>10}  {'SEM':>9}"
    print(header)
    print("  " + "─" * 66)
    for g in group_order:
        v = _group_vals(long, g, column)
        n = len(v)
        print(
            f"  {g:>20}  {n:>5}  {np.mean(v):>10.4f}  "
            f"{np.median(v):>10.4f}  {np.std(v, ddof=1):>10.4f}  "
            f"{stats.sem(v):>9.4f}"
        )
    print(SEP)


def print_stat_results(
    stat_result: dict,
    group_order: list[str],
    alpha: float,
) -> None:
    SEP = "─" * 70
    print(f"\n{SEP}")
    print("  Statistical results")
    print(SEP)
    print(f"  Test    : {stat_result['test_name']}")
    print(f"  Result  : {stat_result['summary']}")

    if "posthoc_df" in stat_result:
        ph = stat_result["posthoc_df"]
        print(f"\n  Post-hoc pairwise adjusted p-values:")
        print(f"  {'Group A':>22}  {'Group B':<22}  {'p (adj)':>12}  Sig.")
        print("  " + "─" * 62)
        for g1, g2 in combinations(group_order, 2):
            p  = ph.loc[g1, g2]
            st = p_to_stars(p, alpha)
            print(f"  {g1:>22}  {g2:<22}  {p:>12.4e}  {st}")
    print(SEP + "\n")


def save_stats_csv(
    stat_result: dict,
    group_order: list[str],
    alpha: float,
    path: str,
) -> None:
    rows = []
    if "posthoc_df" in stat_result:
        ph = stat_result["posthoc_df"]
        for g1, g2 in combinations(group_order, 2):
            p = ph.loc[g1, g2]
            rows.append({
                "group_A":    g1,
                "group_B":    g2,
                "p_adjusted": p,
                "significance": p_to_stars(p, alpha),
                "test": stat_result["test_name"],
            })
    else:
        rows.append({
            "group_A": group_order[0] if len(group_order) == 2 else "all",
            "group_B": group_order[1] if len(group_order) == 2 else "all",
            "p_adjusted": stat_result["p_value"],
            "significance": p_to_stars(stat_result["p_value"], alpha),
            "test": stat_result["test_name"],
        })

    df_out = pd.DataFrame(rows)
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(out_path, index=False)
    print(f"[INFO] Statistics saved → {out_path}")


# ══════════════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    parser = build_parser()
    args   = parser.parse_args()

    if len(args.csvfiles) < 2:
        parser.error("At least two CSV files are required.")

    # ── Resolve statistical settings ──────────────────────────────────────────
    if args.interactive and not args.no_stats:
        args = interactive_settings(args)

    if args.distribution is None and not args.no_stats and args.test == "auto":
        print(
            "[INFO] No --distribution specified; defaulting to 'normal'.\n"
            "       Use --distribution {normal|lognormal|nonparametric} "
            "or --interactive to override."
        )
        args.distribution = "normal"

    # ── Load data ─────────────────────────────────────────────────────────────
    long, group_order = load_data(
        args.csvfiles, args.column, args.labels,
        args.subject_column, args.sep,
    )
    print(
        f"[INFO] Loaded {len(long):,} observations across "
        f"{len(group_order)} groups: {group_order}"
    )
    print_descriptives(long, args.column, group_order)

    # ── Statistics ────────────────────────────────────────────────────────────
    stat_result = None
    if not args.no_stats:
        test_key = select_test(
            len(group_order), args.distribution,
            args.paired, args.test,
        )
        print(
            f"[INFO] Test selected: {test_key}  "
            f"(distribution={args.distribution}, paired={args.paired})"
        )
        stat_result = run_statistics(
            long, args.column, group_order,
            test_key, args.distribution, args.posthoc,
            args.alpha, args.paired, args.subject_column,
        )
        print_stat_results(stat_result, group_order, args.alpha)

        if args.save_stats:
            save_stats_csv(stat_result, group_order, args.alpha, args.save_stats)

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig = make_violin_plot(long, args.column, group_order, args, stat_result)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    print(f"[INFO] Plot saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
