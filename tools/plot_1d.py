#!/usr/bin/env python3
"""Manhattan-style plot of -log10(p) vs genomic position.

Colors points by sign(effect): blue = positive, red = negative.
Chromosome boundaries are shown as alternating background bands with labels.
Boundaries are derived from --fasta (reference genome) or --boundaries FILE.
"""

import argparse
import glob
import math
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# ---------------------------------------------------------------------------
# Publication-quality defaults
# ---------------------------------------------------------------------------
BAND_COLORS = ["#eaeaea", "#ffffff"]
POINT_COLORS = {"pos": "#2166ac", "neg": "#b2182b", "zero": "#888888"}
THRESHOLD_COLOR = "#222222"


def parse_args():
    ap = argparse.ArgumentParser(
        description="Manhattan-style plot of -log10(p) vs genomic position."
    )
    ap.add_argument("hits",
                    help="TSV/whitespace file: pos effect p_value (no header)")
    ap.add_argument("--start", type=int, default=None,
                    help="Start position (inclusive)")
    ap.add_argument("--end", type=int, default=None,
                    help="End position (inclusive)")
    ap.add_argument("--out", default=None,
                    help="Output image path (png/svg/pdf). Omit to show window.")
    ap.add_argument("--title", default=None, help="Plot title")
    ap.add_argument("--point-size", type=float, default=4.0,
                    help="Marker size (default: 4)")
    ap.add_argument("--alpha", type=float, default=0.7,
                    help="Marker transparency (default: 0.7)")
    ap.add_argument("--threshold", type=float, default=None,
                    help="Draw horizontal line at this -log10(p) value")
    ap.add_argument("--fasta", default=None,
                    help="Reference FASTA to derive chromosome boundaries")
    ap.add_argument("--boundaries", default=None,
                    help="TSV: label<TAB>start<TAB>end (one row per chromosome)")
    ap.add_argument("--min-logp", type=float, default=None,
                    help="Filter cutoff -log10(p) value to mark on the y-axis")

    bonf = ap.add_mutually_exclusive_group()
    bonf.add_argument("--num-tests", "-m", type=int, default=None,
                      help="Number of tests for Bonferroni threshold")
    bonf.add_argument("--counts", type=str, nargs="+", default=None,
                      help=".count file(s) from kgwas")
    ap.add_argument("--significance", "-s", type=float, default=0.05,
                    help="Significance level for Bonferroni (default: 0.05)")
    return ap.parse_args()


def iter_hits(path, start=None, end=None):
    """Yield (pos, neglogp, sign) tuples with optional range filtering."""
    p_eps = 1e-300
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s[0] == "#":
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            try:
                pos = int(parts[0])
                eff = float(parts[1])
                p = float(parts[2])
            except ValueError:
                continue
            if start is not None and pos < start:
                continue
            if end is not None and pos > end:
                continue
            p = max(p, p_eps)
            neglogp = -math.log10(p)
            sign = 1 if eff > 0 else (-1 if eff < 0 else 0)
            yield pos, neglogp, sign


def load_boundaries(path):
    """Load chromosome boundaries from a TSV file.

    Expected format: label<TAB>start<TAB>end (no header).
    Returns list of (label, start, end) sorted by start.
    """
    boundaries = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            try:
                label = parts[0].strip()
                start = int(parts[1])
                end = int(parts[2])
                boundaries.append((label, start, end))
            except ValueError:
                continue
    boundaries.sort(key=lambda x: x[1])
    return boundaries


def boundaries_from_fasta(fasta_path):
    """Compute cumulative chromosome boundaries from a FASTA file.

    Sequences whose names match common chromosome patterns (Chr, Chrom,
    chromosome) are labelled by their number.  All remaining small scaffolds
    are collapsed into a single "Other" band at the end.

    Returns list of (label, start, end) in FASTA order.
    """
    chrom_pat = re.compile(
        r'^(?:chr(?:om(?:osome)?)?)[_\-]?0*(\d+)', re.IGNORECASE
    )
    sequences = []
    name = None
    length = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    sequences.append((name, length))
                name = line[1:].strip().split()[0]
                length = 0
            else:
                length += len(line.strip())
        if name is not None:
            sequences.append((name, length))

    boundaries = []
    offset = 0
    scaffold_start = None
    for seq_name, seq_len in sequences:
        m = chrom_pat.match(seq_name)
        if m:
            label = m.group(1)
            boundaries.append((label, offset, offset + seq_len))
        else:
            if scaffold_start is None:
                scaffold_start = offset
        offset += seq_len

    if scaffold_start is not None:
        boundaries.append(("Other", scaffold_start, offset))

    return boundaries


def main():
    args = parse_args()

    xs_pos, ys_pos = [], []
    xs_neg, ys_neg = [], []
    xs_zero, ys_zero = [], []

    for pos, nlp, sign in iter_hits(args.hits, args.start, args.end):
        if sign > 0:
            xs_pos.append(pos)
            ys_pos.append(nlp)
        elif sign < 0:
            xs_neg.append(pos)
            ys_neg.append(nlp)
        else:
            xs_zero.append(pos)
            ys_zero.append(nlp)

    if not (xs_pos or xs_neg or xs_zero):
        print("No points to plot (check range or input).", file=sys.stderr)
        sys.exit(0)

    # --- Chromosome boundaries ---
    boundaries = []
    if args.fasta:
        boundaries = boundaries_from_fasta(args.fasta)
        chrom_names = [b[0] for b in boundaries]
        print(f"Loaded {len(boundaries)} regions from FASTA: {', '.join(chrom_names)}",
              file=sys.stderr)
    elif args.boundaries:
        boundaries = load_boundaries(args.boundaries)
        print(f"Loaded {len(boundaries)} chromosome boundaries from file",
              file=sys.stderr)

    # --- Bonferroni threshold ---
    threshold = args.threshold
    threshold_label = "Significance threshold"
    num_tests = args.num_tests
    if args.counts:
        num_tests = 0
        for pattern in args.counts:
            for fpath in glob.glob(pattern):
                with open(fpath) as f:
                    for line in f:
                        parts = line.strip().split("\t")
                        if len(parts) >= 2:
                            num_tests += int(parts[1])
        print(f"Total tests from .count files: {num_tests:,}", file=sys.stderr)
    if num_tests is not None and num_tests > 0:
        p_thresh = args.significance / num_tests
        threshold = -math.log10(p_thresh)
        pct = args.significance * 100
        if pct == int(pct):
            threshold_label = f"Bonferroni {int(pct)}%"
        else:
            threshold_label = f"Bonferroni {pct}%"

    # --- Figure setup ---
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.linewidth": 0.8,
        "axes.labelsize": 13,
        "axes.titlesize": 16,
        "axes.titleweight": "semibold",
        "axes.titlepad": 14,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
    })

    fig, ax = plt.subplots(figsize=(14, 5))

    # --- Alternating chromosome background bands ---
    if boundaries:
        for i, (label, bstart, bend) in enumerate(boundaries):
            color = BAND_COLORS[i % 2]
            ax.axvspan(bstart, bend, facecolor=color, edgecolor="none",
                       zorder=0)

    # --- Scatter points ---
    scatter_kw = dict(s=args.point_size, alpha=args.alpha, marker="o",
                      edgecolors="none", linewidths=0, rasterized=True)
    if xs_pos:
        ax.scatter(xs_pos, ys_pos, c=POINT_COLORS["pos"],
                   label="Effect > 0", zorder=2, **scatter_kw)
    if xs_neg:
        ax.scatter(xs_neg, ys_neg, c=POINT_COLORS["neg"],
                   label="Effect < 0", zorder=2, **scatter_kw)
    if xs_zero:
        ax.scatter(xs_zero, ys_zero, c=POINT_COLORS["zero"],
                   label="Effect = 0", zorder=2, **scatter_kw)

    # --- Threshold line ---
    if threshold is not None:
        ax.axhline(y=threshold, color=THRESHOLD_COLOR, linestyle="--",
                   linewidth=1.0, zorder=3,
                   label=f"{threshold_label} ($-\\log_{{10}}$={threshold:.1f})")

    # --- Axis labels and title ---
    ax.set_ylabel("$-\\log_{10}(p)$")
    ttl = args.title or os.path.basename(args.hits)
    if args.start is not None or args.end is not None:
        s = args.start if args.start is not None else ""
        e = args.end if args.end is not None else ""
        ttl += f"  [{s}:{e}]"
    ax.set_title(ttl, fontfamily="serif")

    # --- X-axis: chromosome labels at midpoints ---
    if boundaries:
        midpoints = [(bstart + bend) / 2 for _, bstart, bend in boundaries]
        labels = [f"Chr {lbl}" if lbl != "Other" else lbl
                  for lbl, _, _ in boundaries]
        ax.set_xticks(midpoints)
        ax.set_xticklabels(labels)
        ax.set_xlim(boundaries[0][1], boundaries[-1][2])
    else:
        ax.set_xlabel("Position (bp)")
        ax.xaxis.set_major_formatter(
            ticker.FuncFormatter(lambda x, _: f"{x / 1e6:.1f} Mb"))

    # --- Y-axis: start at filter cutoff or zero ---
    if args.min_logp is not None:
        ax.set_ylim(bottom=args.min_logp)
        ax.axhline(y=args.min_logp, color="#999999", linestyle="-",
                   linewidth=0.8, zorder=1,
                   label=f"Filter cutoff ($-\\log_{{10}}$={args.min_logp:g})")
    else:
        ax.set_ylim(bottom=0)

    # --- Legend ---
    leg = ax.legend(loc="upper right", fontsize=9, frameon=True,
                    fancybox=False, edgecolor="#cccccc", framealpha=0.9)
    leg.get_frame().set_linewidth(0.6)

    # --- Remove top/right spines ---
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()

    if args.out:
        dpi = 300 if args.out.endswith(".png") else 150
        fig.savefig(args.out, dpi=dpi, bbox_inches="tight", facecolor="white")
        print(f"Plot saved to: {args.out}", file=sys.stderr)
    else:
        plt.show()


if __name__ == "__main__":
    main()
