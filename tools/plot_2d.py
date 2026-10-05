#!/usr/bin/env python3
"""3D Manhattan plot for pairwise epistasis data on a genome-vs-genome grid.

Binning strategy:
  1. High-resolution grid (auto-sized or user-specified)
  2. Count pairs per cell → drop cells below density percentile (noise)
  3. Keep max -log10(p) per surviving cell

Optional: --marginal embeds the 1D Manhattan as floor projections along both axes.
Optional: --pos1 / --pos2 for zoomed regional views.
"""

import argparse
import os
import re

import numpy as np
import pandas as pd
import plotly.graph_objects as go


def parse_args():
    p = argparse.ArgumentParser(description="3D epistasis Manhattan plot")
    p.add_argument("--fasta", required=True,
                    help="Assembly FASTA (chromosome lengths parsed from headers)")
    p.add_argument("--tsv", required=True,
                    help="Epistasis TSV: pos_a  pos_b  effect  pvalue")
    p.add_argument("--title", default="2D epistasis",
                    help="Plot subtitle")
    p.add_argument("-o", "--output", default="manhattan3d.html",
                    help="Output HTML file")
    p.add_argument("--bin-size", type=float, default=0,
                    help="Grid bin size in bp (0 = auto: genome / 200000)")
    p.add_argument("--density-percentile", type=float, default=90,
                    help="Percentile of density distribution used as threshold (default: 90)")
    p.add_argument("--marginal", default=None,
                    help="1D filtered TSV (pos, effect, pval) for floor projections")
    p.add_argument("--effect-percentile", type=float, default=90,
                    help="Keep only pairs with |effect| above this percentile (default: 90=top 10%%)")
    p.add_argument("--pos1", default=None,
                    help="Pos A range: START:END (bp) for zoomed view")
    p.add_argument("--pos2", default=None,
                    help="Pos B range: START:END (bp) for zoomed view")
    return p.parse_args()


def parse_chromosomes(fasta_path):
    names, lengths = [], []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line.strip().lstrip(">")
            if not header.startswith("Chrom"):
                continue
            m = re.search(r"length_(\d+)", header)
            if m:
                names.append(header.split("_")[0])
                lengths.append(int(m.group(1)))
    return names, np.array(lengths)


def parse_range(s):
    a, b = s.split(":")
    return int(float(a)), int(float(b))


def density_binning(df, bin_size, density_pct, x_range, y_range):
    """High-res grid → drop sparse cells by percentile → max -log10(p) per cell."""
    df = df.copy()
    df["bx"] = ((df["x"] - x_range[0]) / bin_size).astype(np.int64)
    df["by"] = ((df["y"] - y_range[0]) / bin_size).astype(np.int64)

    cell_counts = df.groupby(["bx", "by"]).size()
    n_occupied = len(cell_counts)

    densities = cell_counts.values
    threshold = np.percentile(densities, density_pct)
    dense_cells = cell_counts[cell_counts >= threshold].index
    pct_kept = 100 * len(dense_cells) / n_occupied if n_occupied > 0 else 0
    print(f"  Density distribution: min={densities.min()}, median={int(np.median(densities))}, "
          f"p{density_pct:.0f}={threshold:.0f}, max={densities.max()}")
    print(f"  Cells: {n_occupied:,} occupied, {len(dense_cells):,} kept ({pct_kept:.1f}%)")

    df = df.set_index(["bx", "by"])
    df = df.loc[df.index.isin(dense_cells)].reset_index()

    # Max -log10(p) per cell
    idx = df.groupby(["bx", "by"])["nlp"].idxmax()
    binned = df.loc[idx].copy()

    # Cell centers as positions
    binned["x"] = x_range[0] + (binned["bx"] + 0.5) * bin_size
    binned["y"] = y_range[0] + (binned["by"] + 0.5) * bin_size
    binned = binned.drop(columns=["bx", "by"]).reset_index(drop=True)

    return binned


def main():
    args = parse_args()

    # --- 1. Chromosome geometry ---
    chrom_names, chrom_lengths = parse_chromosomes(args.fasta)
    n_chr = len(chrom_names)
    chrom_short = [n.replace("Chrom", "") for n in chrom_names]
    cum_starts = np.concatenate([[0], np.cumsum(chrom_lengths[:-1])])
    cum_ends = np.cumsum(chrom_lengths)
    chrom_mids = (cum_starts + cum_ends) / 2
    genome_length = float(cum_ends[-1])

    print(f"Genome: {n_chr} chromosomes, {int(genome_length):,} bp")

    # --- 2. Load data ---
    import sys, time
    file_mb = os.path.getsize(args.tsv) / 1e6
    print(f"Loading {args.tsv} ({file_mb:.0f} MB)...", flush=True)
    t0 = time.time()
    df = pd.read_csv(args.tsv, sep="\t", header=None,
                      names=["x", "y", "effect", "pvalue"])
    df["nlp"] = -np.log10(df["pvalue"])
    print(f"  {len(df):,} raw pairs in {time.time()-t0:.1f}s, -log10(p): {df['nlp'].min():.1f} – {df['nlp'].max():.1f}", flush=True)

    # --- 2b. Effect percentile filter ---
    if args.effect_percentile > 0:
        abs_eff = df["effect"].abs()
        threshold = np.percentile(abs_eff, args.effect_percentile)
        before = len(df)
        df = df[abs_eff >= threshold]
        print(f"  Effect filter: |effect| >= {threshold:.6f} (p{args.effect_percentile:.0f}), "
              f"{before:,} → {len(df):,}")

    # --- 3. Filter to region ---
    x_range = (0.0, genome_length)
    y_range = (0.0, genome_length)

    if args.pos1:
        x_range = parse_range(args.pos1)
        df = df[(df["x"] >= x_range[0]) & (df["x"] < x_range[1])]
    if args.pos2:
        y_range = parse_range(args.pos2)
        df = df[(df["y"] >= y_range[0]) & (df["y"] < y_range[1])]

    if len(df) == 0:
        print("No pairs in specified region.")
        return

    print(f"  {len(df):,} pairs in region")

    # --- 4. Density binning ---
    span = max(x_range[1] - x_range[0], y_range[1] - y_range[0])
    bin_size = args.bin_size if args.bin_size > 0 else span / 5000000
    print(f"Binning at {bin_size/1e6:.2f} Mb resolution, dropping <p{args.density_percentile:.0f} density...")

    df = density_binning(df, bin_size, args.density_percentile, x_range, y_range)
    print(f"  {len(df):,} binned points")

    z_min, z_max = float(df["nlp"].min()), float(df["nlp"].max())
    z_floor = z_min - 0.3

    # --- 5. Load marginal if provided ---
    marg_df = None
    if args.marginal and os.path.isfile(args.marginal):
        print(f"Loading marginal: {args.marginal}")
        marg_df = pd.read_csv(args.marginal, sep="\t", header=None,
                               names=["pos", "effect", "pval"])
        marg_df["nlp"] = -np.log10(marg_df["pval"])
        print(f"  {len(marg_df):,} marginal points")

    # --- 6. Chromosome mapping for hover ---
    def pos_to_chrom(positions):
        chroms = np.empty(len(positions), dtype=object)
        local = np.empty(len(positions), dtype=np.int64)
        for i, p in enumerate(positions):
            idx = np.searchsorted(cum_ends, p, side="right")
            if idx < n_chr:
                chroms[i] = chrom_short[idx]
                local[i] = p - cum_starts[idx]
            else:
                chroms[i] = "scf"
                local[i] = p
        return chroms, local

    print("Building plot...")
    x_chr, x_loc = pos_to_chrom(df["x"].values)
    y_chr, y_loc = pos_to_chrom(df["y"].values)

    hover_text = [
        f"<b>{xc}:{xl:,}</b> × <b>{yc}:{yl:,}</b><br>"
        f"Effect: {eff:+.7f}<br>"
        f"P: {pv:.2e}<br>"
        f"-log₁₀(p): {nlp:.2f}"
        for xc, xl, yc, yl, eff, pv, nlp in zip(
            x_chr, x_loc, y_chr, y_loc,
            df["effect"], df["pvalue"], df["nlp"]
        )
    ]

    # --- 7. Derived columns ---
    size_min, size_max = 2.5, 8.0
    df["msize"] = size_min + (size_max - size_min) * (
        (df["nlp"] - z_min) / max(z_max - z_min, 1e-9)
    )

    # --- 8. Build figure ---
    fig = go.Figure()
    gl = genome_length

    # ---- 8a. Checkered floor tiles (stronger contrast) ----
    color_light = "rgba(255,255,255,1)"
    color_dark = "rgba(225,225,235,1)"
    color_diag_light = "rgba(248,242,228,1)"
    color_diag_dark = "rgba(228,220,200,1)"

    for ix in range(n_chr):
        for iy in range(n_chr):
            is_diag = ix == iy
            c = (color_diag_light if is_diag else color_light) if (ix+iy)%2==0 \
                else (color_diag_dark if is_diag else color_dark)
            x0, x1 = float(cum_starts[ix]), float(cum_ends[ix])
            y0, y1 = float(cum_starts[iy]), float(cum_ends[iy])
            fig.add_trace(go.Mesh3d(
                x=[x0,x1,x1,x0], y=[y0,y0,y1,y1], z=[z_floor]*4,
                i=[0,0], j=[1,2], k=[2,3],
                color=c, opacity=1.0, showlegend=False, hoverinfo="skip", flatshading=True,
            ))

    # ---- 8b. Floor sub-grid lines ----
    z_grid = z_floor + 0.05
    n_subdiv = 5
    grid_x, grid_y, grid_z = [], [], []
    sub_positions = []
    for i in range(n_chr):
        s, e = float(cum_starts[i]), float(cum_ends[i])
        step = (e - s) / n_subdiv
        for k in range(1, n_subdiv):
            sub_positions.append(s + k * step)
    for xp in sub_positions:
        grid_x.extend([xp, xp, None]); grid_y.extend([0, gl, None]); grid_z.extend([z_grid]*2+[None])
    for yp in sub_positions:
        grid_x.extend([0, gl, None]); grid_y.extend([yp, yp, None]); grid_z.extend([z_grid]*2+[None])
    fig.add_trace(go.Scatter3d(
        x=grid_x, y=grid_y, z=grid_z, mode="lines",
        line=dict(color="rgba(120,120,145,0.55)", width=1.5),
        showlegend=False, hoverinfo="skip",
    ))

    # ---- 8c. Outer wireframe box ----
    z_top = z_max + 0.3
    box_x, box_y, box_z = [], [], []
    for bx, by in [(0,0),(gl,0),(gl,gl),(0,gl)]:
        box_x.extend([bx,bx,None]); box_y.extend([by,by,None]); box_z.extend([z_floor,z_top,None])
    box_x.extend([0,gl,gl,0,0,None]); box_y.extend([0,0,gl,gl,0,None]); box_z.extend([z_top]*5+[None])
    fig.add_trace(go.Scatter3d(
        x=box_x, y=box_y, z=box_z, mode="lines",
        line=dict(color="rgba(180,180,200,0.25)", width=1.5),
        showlegend=False, hoverinfo="skip",
    ))

    # ---- 8d. Floor border ----
    fig.add_trace(go.Scatter3d(
        x=[0,gl,gl,0,0], y=[0,0,gl,gl,0], z=[z_floor]*5,
        mode="lines", line=dict(color="rgba(70,70,90,0.6)", width=3),
        showlegend=False, hoverinfo="skip",
    ))

    # ---- 8e. 1D marginal floor projections ----
    if marg_df is not None:
        MAX_MARGINAL_PTS = 50000
        if len(marg_df) > MAX_MARGINAL_PTS:
            print(f"  Subsampling marginal: {len(marg_df):,} → {MAX_MARGINAL_PTS:,}")
            marg_df = marg_df.sample(n=MAX_MARGINAL_PTS, random_state=42)
        m_pos = marg_df["pos"].values
        m_nlp = marg_df["nlp"].values
        m_eff = marg_df["effect"].values
        m_nlp_max = float(m_nlp.max())
        marg_spread = genome_length * 0.10
        m_disp = (m_nlp / m_nlp_max) * marg_spread

        m_mask_pos = m_eff > 0
        m_mask_neg = ~m_mask_pos

        # Build hover text for marginal points
        m_chr, m_loc = pos_to_chrom(m_pos)
        m_hover = [
            f"<b>Chr {c}:{l:,}</b><br>"
            f"Effect: {e:+.5f}<br>"
            f"-log₁₀(p): {n:.1f}"
            for c, l, e, n in zip(m_chr, m_loc, m_eff, m_nlp)
        ]

        print(f"  Adding {len(marg_df):,} marginal points to 3D floor")

        for axis_label, mx, my in [
            ("locus A", m_pos, -m_disp),
            ("locus B", -m_disp, m_pos),
        ]:
            # Positive effect (red, matching 3D scatter)
            fig.add_trace(go.Scatter3d(
                x=mx[m_mask_pos], y=my[m_mask_pos],
                z=np.full(m_mask_pos.sum(), z_floor),
                mode="markers",
                marker=dict(size=0.8, color="#b2182b", opacity=0.7,
                            symbol="circle", line=dict(width=0)),
                text=[m_hover[i] for i in np.where(m_mask_pos)[0]],
                hoverinfo="text",
                showlegend=False,
            ))
            # Negative effect (blue, matching 3D scatter)
            fig.add_trace(go.Scatter3d(
                x=mx[m_mask_neg], y=my[m_mask_neg],
                z=np.full(m_mask_neg.sum(), z_floor),
                mode="markers",
                marker=dict(size=0.8, color="#2166ac", opacity=0.7,
                            symbol="circle", line=dict(width=0)),
                text=[m_hover[i] for i in np.where(m_mask_neg)[0]],
                hoverinfo="text",
                showlegend=False,
            ))

    # ---- 8f. Main scatter — positive effect ----
    mask_pos = df["effect"] > 0
    fig.add_trace(go.Scatter3d(
        x=df.loc[mask_pos,"x"], y=df.loc[mask_pos,"y"], z=df.loc[mask_pos,"nlp"],
        mode="markers",
        marker=dict(
            size=df.loc[mask_pos,"msize"],
            color=df.loc[mask_pos,"nlp"],
            colorscale=[[0,"rgb(255,170,150)"],[0.5,"rgb(220,60,40)"],[1,"rgb(140,10,10)"]],
            cmin=z_min, cmax=z_max, opacity=0.85,
            line=dict(width=0.4, color="rgba(60,0,0,0.25)"), showscale=False,
        ),
        text=[hover_text[i] for i in df.index[mask_pos]],
        hoverinfo="text", name="<b>Positive</b> effect",
    ))

    # ---- 8h. Main scatter — negative effect ----
    mask_neg = df["effect"] <= 0
    fig.add_trace(go.Scatter3d(
        x=df.loc[mask_neg,"x"], y=df.loc[mask_neg,"y"], z=df.loc[mask_neg,"nlp"],
        mode="markers",
        marker=dict(
            size=df.loc[mask_neg,"msize"],
            color=df.loc[mask_neg,"nlp"],
            colorscale=[[0,"rgb(160,180,255)"],[0.5,"rgb(50,80,210)"],[1,"rgb(15,20,120)"]],
            cmin=z_min, cmax=z_max, opacity=0.85,
            line=dict(width=0.4, color="rgba(0,0,60,0.25)"), showscale=False,
        ),
        text=[hover_text[i] for i in df.index[mask_neg]],
        hoverinfo="text", name="<b>Negative</b> effect",
    ))

    # --- 9. Layout ---
    fig.update_layout(
        title=dict(
            text=f"<b>3D Epistasis Manhattan Plot</b> — {args.title}",
            font=dict(size=22, family="Arial", color="rgb(40,40,60)"),
            x=0.02, xanchor="left",
        ),
        scene=dict(
            xaxis=dict(
                title=dict(text="Locus A (chromosome)", font=dict(size=13)),
                tickvals=chrom_mids.tolist(), ticktext=chrom_short, tickfont=dict(size=11),
                showgrid=True, gridcolor="rgba(180,180,200,0.18)", zeroline=False,
                backgroundcolor="rgba(245,245,250,0.5)", linecolor="rgb(100,100,120)", showline=True,
            ),
            yaxis=dict(
                title=dict(text="Locus B (chromosome)", font=dict(size=13)),
                tickvals=chrom_mids.tolist(), ticktext=chrom_short, tickfont=dict(size=11),
                showgrid=True, gridcolor="rgba(180,180,200,0.18)", zeroline=False,
                backgroundcolor="rgba(245,250,245,0.5)", linecolor="rgb(100,100,120)", showline=True,
            ),
            zaxis=dict(
                title=dict(text="-log₁₀(p)", font=dict(size=13)), tickfont=dict(size=11),
                showgrid=True, gridcolor="rgba(160,160,180,0.25)", zeroline=False,
                backgroundcolor="rgba(248,248,255,0.5)",
            ),
            aspectmode="manual", aspectratio=dict(x=1, y=1, z=0.6),
            camera=dict(eye=dict(x=1.1, y=-1.4, z=1.1), up=dict(x=0, y=0, z=1)),
        ),
        legend=dict(
            yanchor="top", y=0.98, xanchor="left", x=0.01, font=dict(size=13),
            bgcolor="rgba(255,255,255,0.85)", bordercolor="rgba(100,100,120,0.3)",
            borderwidth=1, itemsizing="constant",
        ),
        margin=dict(l=0, r=0, t=60, b=0),
        paper_bgcolor="rgb(252,252,255)",
        font=dict(family="Arial"),
        uirevision="constant",
    )

    # --- 10. Save ---
    print(f"Writing {args.output}...")
    fig.write_html(args.output, include_plotlyjs=True)
    size_mb = os.path.getsize(args.output) / 1e6
    print(f"Done! {size_mb:.1f} MB — open {args.output} in a browser.")


if __name__ == "__main__":
    main()
