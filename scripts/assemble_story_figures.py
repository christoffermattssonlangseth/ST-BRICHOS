"""Assemble the three manuscript story figures from analysis outputs.

Reads TSVs from results/tables/attenuation/ and existing panel PNGs
from results/figures/manuscript/. Produces three Cell-styled SVG/PNG
composites:

    figure_story1_design_disease.svg
    figure_story2_regional_rescue.svg
    figure_story3_mechanistic_axis.svg

Style: Helvetica 7-8 pt, double-column 174 mm width, vector text.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
TBL = ROOT / "results" / "tables" / "attenuation"
FIG = ROOT / "results" / "figures" / "manuscript"
OUT = FIG  # write story figures alongside

# Cell column widths (inches)
W_DOUBLE = 6.85

C_WT, C_PBS, C_BRI = "#4a4a4a", "#1f77b4", "#d62728"
C_AX = "#222222"

plt.rcParams.update({
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "Helvetica",
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.titleweight": "regular",
    "axes.labelsize": 7,
    "axes.linewidth": 0.6,
    "axes.edgecolor": C_AX,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "legend.frameon": False,
    "lines.linewidth": 0.8,
    "figure.dpi": 150,
})


def _label(ax, txt: str):
    ax.text(-0.05, 1.05, txt, transform=ax.transAxes,
            fontsize=10, fontweight="bold", family="Helvetica",
            ha="right", va="bottom")


def _strip(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def _imshow(ax, path, label=None, title=None):
    if Path(path).exists():
        ax.imshow(mpimg.imread(path))
    else:
        ax.text(0.5, 0.5, f"missing\n{Path(path).name}",
                ha="center", va="center", transform=ax.transAxes,
                color="crimson", fontsize=7)
    _strip(ax)
    if label:
        _label(ax, label)
    if title:
        ax.set_title(title, loc="left", pad=3, fontsize=7.5)


def save(fig, name: str):
    out = OUT / f"{name}.svg"
    fig.savefig(out, format="svg", bbox_inches="tight", pad_inches=0.05,
                facecolor="white")
    fig.savefig(OUT / f"{name}.png", format="png", bbox_inches="tight",
                pad_inches=0.05, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


# ---------- Story 1: design + disease signature -----------------------------

def figure_story1():
    fig = plt.figure(figsize=(W_DOUBLE, 6.8))
    gs = GridSpec(2, 3, figure=fig,
                  height_ratios=[1.05, 0.95],
                  hspace=0.50, wspace=0.55)

    # A: schematic (re-use figure0)
    ax_a = fig.add_subplot(gs[0, :])
    _imshow(ax_a, FIG / "figure0_schematic.png", label="A")

    # B: disease signature top up-genes
    disease = pd.read_csv(TBL / "disease_signature.tsv",
                          sep="\t", index_col=0)
    top_up = disease.sort_values("lfc", ascending=False).head(15)
    top_dn = disease.sort_values("lfc").head(15)

    ax_b = fig.add_subplot(gs[1, 0])
    y = np.arange(len(top_up))
    ax_b.barh(y, top_up["lfc"], color="#d97c7c", edgecolor="#222",
              linewidth=0.4)
    ax_b.set_yticks(y); ax_b.set_yticklabels(top_up.index, fontsize=6.5)
    ax_b.invert_yaxis()
    ax_b.set_xlabel("log$_2$FC (PBS vs WT)")
    ax_b.set_title("Top 15 disease-up genes", loc="left",
                   fontsize=7.5, pad=2)
    _label(ax_b, "B")
    for s in ("top", "right"):
        ax_b.spines[s].set_visible(False)

    # C: disease signature top down-genes
    ax_c = fig.add_subplot(gs[1, 1])
    y = np.arange(len(top_dn))
    ax_c.barh(y, top_dn["lfc"], color="#7c9bd9", edgecolor="#222",
              linewidth=0.4)
    ax_c.set_yticks(y); ax_c.set_yticklabels(top_dn.index, fontsize=6.5)
    ax_c.invert_yaxis()
    ax_c.set_xlabel("log$_2$FC (PBS vs WT)")
    ax_c.set_title("Top 15 disease-down genes", loc="left",
                   fontsize=7.5, pad=2)
    _label(ax_c, "C")
    for s in ("top", "right"):
        ax_c.spines[s].set_visible(False)

    # D: PIG signature membership in disease set (overlap stat)
    from utils.gene_sets import PIG_CHEN2020
    pig_set = pd.Index(PIG_CHEN2020)
    in_disease = pig_set.intersection(disease.index)
    pig_lfc = disease.loc[in_disease, "lfc"].sort_values(ascending=False)

    ax_d = fig.add_subplot(gs[1, 2])
    colors = ["#d62728" if v > 0 else "#1f77b4" for v in pig_lfc]
    ax_d.barh(np.arange(len(pig_lfc)), pig_lfc.values,
              color=colors, edgecolor="#222", linewidth=0.4)
    ax_d.set_yticks(np.arange(len(pig_lfc)))
    ax_d.set_yticklabels(pig_lfc.index, fontsize=5.5)
    ax_d.invert_yaxis()
    ax_d.set_xlabel("log$_2$FC (PBS vs WT)")
    ax_d.set_title(f"PIG signature genes (n={len(pig_lfc)})",
                   loc="left", fontsize=7.5, pad=2)
    ax_d.axvline(0, color="#888", linewidth=0.4)
    _label(ax_d, "D")
    for s in ("top", "right"):
        ax_d.spines[s].set_visible(False)

    save(fig, "figure_story1_design_disease")


# ---------- Story 2: three orthogonal rescue tests --------------------------

def figure_story2():
    rescue = pd.read_csv(TBL / "regional_rescue_stats.tsv",
                         sep="\t", index_col=0).sort_values("slope")
    responder = pd.read_csv(TBL / "responder_fraction_by_region.tsv",
                            sep="\t", index_col=0)
    perm = pd.read_csv(TBL / "rescue_permutation_null.tsv",
                       sep="\t", index_col=0)

    # Align rows
    common = [r for r in rescue.index if r in responder.index]
    rescue = rescue.loc[common]
    responder = responder.loc[common]
    perm = perm.reindex(common)

    fig = plt.figure(figsize=(W_DOUBLE, 4.6))
    gs = GridSpec(1, 3, figure=fig, wspace=0.20,
                  width_ratios=[1.4, 1.0, 1.0])

    # A: sign concordance
    ax_a = fig.add_subplot(gs[0, 0])
    y = np.arange(len(rescue))
    ax_a.barh(y, rescue["frac_concordant"] - 0.5, left=0.5,
              color="#4d8db5", alpha=0.85, edgecolor="#222",
              linewidth=0.4)
    for i, p in enumerate(rescue["binom_p"].values):
        if p < 1e-3: star = "***"
        elif p < 1e-2: star = "**"
        elif p < 5e-2: star = "*"
        else: star = ""
        ax_a.text(rescue["frac_concordant"].iloc[i] + 0.005, i,
                  star, fontsize=7, va="center")
    ax_a.axvline(0.5, color="#bbb", linewidth=0.6, linestyle=":")
    ax_a.set_yticks(y); ax_a.set_yticklabels(rescue.index, fontsize=7)
    ax_a.set_xlim(0.3, 1.0)
    ax_a.set_xlabel("fraction of disease genes\nmoved by BRICHOS in rescue dir.")
    ax_a.set_title("Sign concordance", loc="left", fontsize=7.5, pad=2)
    _label(ax_a, "A")
    for s in ("top", "right"):
        ax_a.spines[s].set_visible(False)

    # B: rescue slope with bootstrap CI + permutation p-values overlaid
    ax_b = fig.add_subplot(gs[0, 1])
    y = np.arange(len(rescue))
    ax_b.errorbar(rescue["slope"], y,
                  xerr=[rescue["slope"] - rescue["slope_ci_low"],
                        rescue["slope_ci_high"] - rescue["slope"]],
                  fmt="o", color="#222", ecolor="#888",
                  capsize=2, markersize=4)
    ax_b.axvline(0, color="#bbb", linewidth=0.6, linestyle=":")
    ax_b.axvline(-1, color="#d62728", linewidth=0.6, linestyle="--",
                 alpha=0.4)
    # annotate permutation p
    for i, region in enumerate(rescue.index):
        if region in perm.index and pd.notna(perm.loc[region, "p_one_sided"]):
            p = perm.loc[region, "p_one_sided"]
            txt = f"p={p:.3f}" if p > 0.001 else "p≤.012"
            color = "#222" if p < 0.05 else "#999"
            ax_b.text(rescue["slope_ci_high"].iloc[i] + 0.01, i,
                      txt, fontsize=6, va="center", color=color)
    ax_b.set_yticks(y); ax_b.set_yticklabels([""] * len(y))
    ax_b.set_xlabel("rescue slope ($LFC_{BRI-PBS}$ on $LFC_{disease}$)")
    ax_b.set_title("Rescue slope + perm. null", loc="left",
                   fontsize=7.5, pad=2)
    _label(ax_b, "B")
    for s in ("top", "right"):
        ax_b.spines[s].set_visible(False)

    # C: spot responder fraction per mouse — bar with mwu p-values
    ax_c = fig.add_subplot(gs[0, 2])
    rf = responder
    y = np.arange(len(rf))
    ax_c.barh(y - 0.18, rf["mean_pbs"], height=0.36, color=C_PBS,
              alpha=0.85, edgecolor="#222", linewidth=0.4, label="PBS")
    ax_c.barh(y + 0.18, rf["mean_bri"], height=0.36, color=C_BRI,
              alpha=0.85, edgecolor="#222", linewidth=0.4, label="BRICHOS")
    for i, p in enumerate(rf["mwu_p"].values):
        if pd.isna(p): continue
        if p < 5e-2: star = "*"
        else: star = ""
        if star:
            xpos = max(rf["mean_pbs"].iloc[i], rf["mean_bri"].iloc[i]) + 0.005
            ax_c.text(xpos, i, star, fontsize=8, va="center")
    ax_c.set_yticks(y); ax_c.set_yticklabels([""] * len(y))
    ax_c.set_xlabel("PIG-rescued spot fraction\n(per-mouse mean)")
    ax_c.set_title("Spot responder fraction (MWU)", loc="left",
                   fontsize=7.5, pad=2)
    ax_c.legend(loc="lower right", fontsize=6, frameon=False)
    _label(ax_c, "C")
    for s in ("top", "right"):
        ax_c.spines[s].set_visible(False)

    save(fig, "figure_story2_regional_rescue")


# ---------- Story 3: mechanistic axis --------------------------------------

def _draw_cascade(ax):
    """Schematic cascade IFN -> IRM -> A1/Oligo + neuronal restoration.

    Layout: BRICHOS (left) -> IFN -> IRM, then branches to A1 (top
    branch) and Oligo stress (bottom branch). Neuronal IEG sits at
    far right as the downstream functional consequence of dampened
    inflammation.
    """
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    _strip(ax)
    DN = r"$\downarrow$"
    UP = r"$\uparrow$"

    bw, bh = 0.115, 0.22  # box width / height
    nodes = [
        (0.07, 0.50, "BRICHOS",                            "#666"),
        (0.27, 0.50, f"{DN} IFN-$\\gamma$ /\nantiviral",   "#aa6cd1"),
        (0.47, 0.50, f"{DN} IRM\nmicroglia",               "#d97c7c"),
        (0.69, 0.78, f"{DN} A1\nastrocytes",               "#d49a55"),
        (0.69, 0.22, f"{DN} Oligo\nstress",                "#7caa6c"),
        (0.91, 0.50, f"{UP} Neuronal\nIEG / synapse",      "#5b8def"),
    ]
    for x, y, text, color in nodes:
        ax.add_patch(FancyBboxPatch((x - bw / 2, y - bh / 2),
                                    bw, bh,
                                    boxstyle="round,pad=0.005,rounding_size=0.012",
                                    facecolor="white",
                                    edgecolor=color, linewidth=1.0))
        ax.text(x, y, text, ha="center", va="center",
                fontsize=6.5, color=color, fontweight="bold")

    def arrow(a, b, color="#666", lw=0.8, style="solid"):
        x0, y0 = nodes[a][:2]
        x1, y1 = nodes[b][:2]
        ls = "-" if style == "solid" else (0, (4, 2))
        ax.add_patch(FancyArrowPatch((x0 + bw / 2, y0), (x1 - bw / 2, y1),
                                     arrowstyle="-|>", mutation_scale=9,
                                     linewidth=lw, color=color,
                                     linestyle=ls,
                                     shrinkA=0, shrinkB=0))

    # Linear cascade
    arrow(0, 1)
    arrow(1, 2)
    # Branch to A1 / Oligo
    arrow(2, 3, color="#d49a55")
    arrow(2, 4, color="#7caa6c")
    # IRM downregulation feeds neuronal restoration (dashed = indirect)
    ax.add_patch(FancyArrowPatch((nodes[2][0] + bw / 2, 0.50),
                                 (nodes[5][0] - bw / 2, 0.50),
                                 arrowstyle="-|>", mutation_scale=9,
                                 linewidth=0.8, color="#5b8def",
                                 linestyle=(0, (3, 2)),
                                 shrinkA=0, shrinkB=0))


def _gsea_panel(ax, region, top_n=6, color_pos="#d97c7c", color_neg="#7caa6c"):
    f = TBL / f"gsea_{region.replace(' ', '_')}.tsv"
    if not f.exists():
        ax.text(0.5, 0.5, f"missing GSEA: {region}",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=6, color="crimson")
        _strip(ax)
        return
    df = pd.read_csv(f, sep="\t")
    if "NES" not in df.columns:
        # try other column names
        nes_col = next((c for c in df.columns if c.lower() == "nes"), None)
        if nes_col is None:
            _strip(ax); return
        df = df.rename(columns={nes_col: "NES"})
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    df = df.dropna(subset=["NES"])
    df["abs_nes"] = df["NES"].abs()
    df = df.sort_values("abs_nes", ascending=False).head(top_n)
    df = df.sort_values("NES")

    y = np.arange(len(df))
    colors = [color_pos if v > 0 else color_neg for v in df["NES"]]
    ax.barh(y, df["NES"], color=colors, edgecolor="#222", linewidth=0.4)

    # short term name — strip GO id, aggressive truncate to fit panel
    def _short(t):
        t = t.split(" (")[0]
        # collapse common phrases
        t = (t.replace("positive regulation of ", "+reg ")
              .replace("negative regulation of ", "-reg ")
              .replace("regulation of ", "reg ")
              .replace("response to ", "resp. ")
              .replace("signaling pathway", "sig.")
              .replace("biosynthetic process", "biosynth.")
              .replace("differentiation", "diff."))
        return t[:28] + ("..." if len(t) > 28 else "")
    terms = [_short(t) for t in df["Term"].astype(str)]
    ax.set_yticks(y)
    ax.set_yticklabels(terms, fontsize=5.5)
    ax.axvline(0, color="#888", linewidth=0.4)
    ax.set_xlabel("NES")
    ax.set_title(region, loc="left", fontsize=7, pad=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def _heatmap_from_data(ax):
    from utils.gene_sets import SIGNATURE_DIRECTION, SIGNATURE_REGISTRY
    df = pd.read_csv(TBL / "broader_signature_rescue.tsv", sep="\t")
    pivot = (df.pivot(index="signature", columns="region",
                      values="median_signed"))
    order_up = [s for s, d in SIGNATURE_DIRECTION.items()
                if d == "up" and s in pivot.index]
    order_dn = [s for s, d in SIGNATURE_DIRECTION.items()
                if d == "down" and s in pivot.index]
    ordered = order_up + order_dn
    pivot = pivot.reindex(ordered).fillna(0)

    ylabels = [
        (r"$\uparrow$ " if SIGNATURE_DIRECTION.get(s) == "up" else r"$\downarrow$ ") + s
        for s in pivot.index
    ]
    vmax = max(0.3, np.nanpercentile(np.abs(pivot.values), 95))
    im = ax.imshow(pivot.values, cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=40, ha="right",
                       fontsize=6.5)
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(ylabels, fontsize=6)
    ax.axhline(len(order_up) - 0.5, color="#222", linewidth=0.5)
    for s in ax.spines.values():
        s.set_visible(False)
    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("median signed LFC\n(rescue+)", fontsize=6.5)
    cbar.ax.tick_params(labelsize=5.5)


def figure_story3():
    # Order: rescue regions first (by sign-concordance strength),
    # then corpus callosum as the discordant comparator.
    panels = [
        ("Supragranualar layers",  "Supragranular layers"),
        ("Infragranular layers",   "Infragranular layers"),
        ("Caudoputamen",           "Caudoputamen"),
        ("Olfactory areas",        "Olfactory areas"),
        ("Meninges",               "Meninges"),
        ("Corpus_callosum",        "Corpus callosum (discordant)"),
    ]

    fig = plt.figure(figsize=(W_DOUBLE, 11.0))
    gs = GridSpec(4, 3, figure=fig,
                  height_ratios=[0.7, 1.6, 1.0, 1.0],
                  hspace=0.85, wspace=1.30)

    # A: cascade diagram (banner)
    ax_a = fig.add_subplot(gs[0, :])
    _draw_cascade(ax_a)
    _label(ax_a, "A")
    ax_a.set_title("Mechanistic cascade reconstructed from independent tests",
                   loc="left", fontsize=8, pad=2)

    # B: broader-signatures heatmap
    ax_b = fig.add_subplot(gs[1, :])
    _heatmap_from_data(ax_b)
    _label(ax_b, "B")
    ax_b.set_title("Cell-type / pathway signature rescue across regions",
                   loc="left", fontsize=8, pad=2)

    # C-H: GSEA per region (2 rows x 3 columns)
    panel_letters = list("CDEFGH")
    for i, ((region_key, title), letter) in enumerate(zip(panels, panel_letters)):
        row = 2 + i // 3
        col = i % 3
        ax = fig.add_subplot(gs[row, col])
        _gsea_panel(ax, region_key, top_n=6)
        ax.set_title(title, loc="left", fontsize=7.5, pad=2)
        _label(ax, letter)

    save(fig, "figure_story3_mechanistic_axis")


def main():
    figure_story1()
    figure_story2()
    figure_story3()


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(ROOT))
    main()
