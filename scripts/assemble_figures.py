"""Assemble manuscript figures in Cell-journal style.

Outputs SVGs (vector text, raster panels) plus PNG previews to
results/figures/manuscript/. Includes a schematic for the experimental
setup (figure0_schematic) and five data figures pulling existing panels
from results/html/nb_figures/.

Style follows Cell guidelines: Helvetica 7-9 pt, double-column 174 mm
(6.85 in) width, no figure suptitle (caption goes in manuscript), bold
uppercase panel labels in the upper-left outside the axes.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
NB = ROOT / "results" / "html" / "nb_figures"
FIG = ROOT / "results" / "figures"
OUT = FIG / "manuscript"
OUT.mkdir(parents=True, exist_ok=True)

# Cell column widths (inches)
W_DOUBLE = 6.85   # 174 mm
W_1P5 = 4.49      # 114 mm
W_SINGLE = 3.35   #  85 mm

# Palette — colourblind-safe, muted
C_WT = "#4a4a4a"
C_PBS = "#1f77b4"
C_BRI = "#d62728"
C_NEU = "#9aa0a6"
C_BG = "#f5f5f5"
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


# ---------- helpers ----------------------------------------------------------

def _label(ax, txt: str):
    """Bold uppercase panel label, upper-left, outside axes (Cell style)."""
    ax.text(-0.04, 1.04, txt, transform=ax.transAxes,
            fontsize=10, fontweight="bold", family="Helvetica",
            ha="right", va="bottom")


def _strip(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def _panel(ax, img_path: Path, label: str | None = None, title: str | None = None):
    if img_path.exists():
        ax.imshow(mpimg.imread(img_path))
    else:
        ax.add_patch(mpatches.Rectangle((0, 0), 1, 1, transform=ax.transAxes,
                                        facecolor=C_BG, edgecolor="none"))
        ax.text(0.5, 0.5, f"missing\n{img_path.name}",
                ha="center", va="center", transform=ax.transAxes,
                color="#999", fontsize=6, family="Helvetica")
    _strip(ax)
    if title:
        ax.set_title(title, loc="left", pad=3, fontsize=7.5)
    if label:
        _label(ax, label)


def save(fig, name: str):
    out_svg = OUT / f"{name}.svg"
    fig.savefig(out_svg, format="svg", bbox_inches="tight",
                pad_inches=0.05, dpi=300, facecolor="white")
    fig.savefig(OUT / f"{name}.png", format="png", bbox_inches="tight",
                pad_inches=0.05, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"wrote {out_svg.relative_to(ROOT)}")


# ---------- schematic --------------------------------------------------------

def _mouse(ax, cx, cy, scale=1.0, color="#dcdcdc", edge=C_AX):
    """Tiny mouse silhouette — total bbox roughly (0.20 x 0.10) * scale.

    `(cx, cy)` is the centre of the body. Drawn facing right.
    """
    s = scale * 0.10  # base size unit
    body = mpatches.Ellipse((cx, cy), 2.0 * s, 1.1 * s,
                            facecolor=color, edgecolor=edge, linewidth=0.7)
    head = mpatches.Ellipse((cx + 1.0 * s, cy + 0.05 * s),
                            0.85 * s, 0.75 * s,
                            facecolor=color, edgecolor=edge, linewidth=0.7)
    ear = mpatches.Circle((cx + 1.20 * s, cy + 0.45 * s),
                          0.20 * s, facecolor=color,
                          edgecolor=edge, linewidth=0.7)
    eye = mpatches.Circle((cx + 1.30 * s, cy + 0.05 * s),
                          0.05 * s, facecolor=edge, edgecolor=edge)
    for p in (body, head, ear, eye):
        ax.add_patch(p)
    ax.plot([cx - 1.0 * s, cx - 1.6 * s, cx - 1.7 * s],
            [cy, cy + 0.2 * s, cy - 0.1 * s],
            color=edge, linewidth=0.7, solid_capstyle="round")


def _syringe(ax, cx, cy, scale=1.0, fill="#e6e6e6", edge=C_AX):
    """Mini syringe pointing right; (cx, cy) is centre of barrel."""
    s = scale * 0.10
    barrel = FancyBboxPatch((cx - 0.9 * s, cy - 0.20 * s),
                            1.8 * s, 0.40 * s,
                            boxstyle="round,pad=0.002,rounding_size=0.002",
                            facecolor=fill, edgecolor=edge, linewidth=0.6)
    ax.add_patch(barrel)
    ax.plot([cx + 0.9 * s, cx + 1.4 * s],
            [cy, cy], color=edge, linewidth=0.8)
    plunger = mpatches.Rectangle((cx - 1.20 * s, cy - 0.30 * s),
                                 0.30 * s, 0.60 * s,
                                 facecolor=edge, edgecolor=edge)
    ax.add_patch(plunger)


def _brain(ax, x, y, scale=1.0, fill="#fde9d9", edge=C_AX):
    """Coronal-ish brain icon."""
    outer = mpatches.Ellipse((x, y), 0.36 * scale, 0.30 * scale,
                             facecolor=fill, edgecolor=edge, linewidth=0.7)
    ax.add_patch(outer)
    # midline
    ax.plot([x, x], [y - 0.13 * scale, y + 0.13 * scale],
            color=edge, linewidth=0.5)
    # ventricles
    for sx in (-0.07, 0.07):
        v = mpatches.Ellipse((x + sx * scale, y + 0.02 * scale),
                             0.04 * scale, 0.08 * scale,
                             facecolor="white", edgecolor=edge, linewidth=0.4)
        ax.add_patch(v)


def _arrow(ax, x0, y0, x1, y1, color=C_AX, lw=0.8, mut=8):
    a = FancyArrowPatch((x0, y0), (x1, y1),
                        arrowstyle="-|>", mutation_scale=mut,
                        linewidth=lw, color=color,
                        shrinkA=0, shrinkB=0)
    ax.add_patch(a)


def _capture_chip(ax, x, y, w=0.5, h=0.32, edge=C_AX):
    """Visium-style capture area: outer chip + dotted spot grid."""
    chip = FancyBboxPatch((x, y), w, h,
                          boxstyle="round,pad=0.005,rounding_size=0.01",
                          facecolor="#ffffff", edgecolor=edge, linewidth=0.7)
    ax.add_patch(chip)
    # spot grid
    n_cols, n_rows = 12, 8
    pad = 0.04
    xs = [x + pad + i * (w - 2 * pad) / (n_cols - 1) for i in range(n_cols)]
    ys = [y + pad + j * (h - 2 * pad) / (n_rows - 1) for j in range(n_rows)]
    for j, yy in enumerate(ys):
        for i, xx in enumerate(xs):
            offset = 0 if j % 2 == 0 else (w - 2 * pad) / (2 * (n_cols - 1))
            ax.add_patch(mpatches.Circle((xx + offset, yy), 0.006,
                                         facecolor="#5b8def",
                                         edgecolor="none"))


def figure0_schematic():
    """Experimental design + analysis workflow."""
    fig = plt.figure(figsize=(W_DOUBLE, 4.4))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 0.8],
                  width_ratios=[1.0, 1.05], wspace=0.05, hspace=0.30)

    # --- A: cohort + dosing ------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.set_xlim(0, 1); ax_a.set_ylim(0, 1)
    _strip(ax_a); _label(ax_a, "A")
    ax_a.set_title("Cohort & dosing", loc="left", pad=2, fontsize=7.5)

    # column anchors
    x_label = 0.02
    x_mouse = 0.34
    x_syr = 0.58
    x_dose = 0.78
    rows = [
        (0.80, "WT",                       "n = 1", C_WT,  False, "#e1e1e1"),
        (0.55, "APP$_{NLGF}$ + PBS",       "n = 3", C_PBS, True,  "#cfe1f2"),
        (0.30, "APP$_{NLGF}$ + BRICHOS",   "n = 3", C_BRI, True,  "#f4cccc"),
    ]
    for y, label, n, col, dosed, body in rows:
        # row separator
        ax_a.add_patch(mpatches.Rectangle((0.0, y - 0.10), 1.0, 0.20,
                                          facecolor="white",
                                          edgecolor="#eee", linewidth=0.4))
        # color tag
        ax_a.add_patch(mpatches.Rectangle((0.0, y - 0.08), 0.008, 0.16,
                                          facecolor=col, edgecolor="none"))
        ax_a.text(x_label, y + 0.045, label, fontsize=7.4, fontweight="bold",
                  color=C_AX, va="center")
        ax_a.text(x_label, y - 0.035, n, fontsize=6.3, color="#666",
                  va="center")
        _mouse(ax_a, x_mouse, y, scale=0.85, color=body)
        if dosed:
            _syringe(ax_a, x_syr, y, scale=0.85, fill=col)
            ax_a.text(x_dose, y, "i.p., 2x/wk", fontsize=6.2,
                      color="#555", va="center")
        else:
            ax_a.text(x_syr, y, "untreated", fontsize=6.3,
                      color="#999", va="center", style="italic")

    # timeline along bottom
    ax_a.plot([0.06, 0.94], [0.06, 0.06], color=C_AX, linewidth=0.7)
    for frac, lab in [(0.06, "wk 0"), (0.50, "wk 6"), (0.94, "wk 12 — sacrifice")]:
        ax_a.plot([frac, frac], [0.05, 0.07], color=C_AX, linewidth=0.7)
        ha = "center" if frac < 0.9 else "right"
        ax_a.text(frac, 0.015, lab, fontsize=6, ha=ha, color="#444")
    ax_a.text(0.5, 0.115, "12 weeks treatment",
              ha="center", fontsize=6.2, color="#555", style="italic")

    # --- B: ST workflow ----------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.set_xlim(0, 1); ax_b.set_ylim(0, 1)
    _strip(ax_b); _label(ax_b, "B")
    ax_b.set_title("Spatial transcriptomics workflow", loc="left", pad=2, fontsize=7.5)

    yy = 0.62
    # 4 stages, evenly spaced
    stage_x = [0.11, 0.36, 0.62, 0.88]

    # 1) coronal section (brain)
    _brain(ax_b, stage_x[0], yy, scale=0.6)
    ax_b.text(stage_x[0], yy - 0.27, "coronal\nsection", ha="center",
              fontsize=6.2, color="#444")

    _arrow(ax_b, stage_x[0] + 0.10, yy, stage_x[1] - 0.13, yy, mut=7)

    # 2) Visium capture chip
    cw, ch = 0.20, 0.26
    _capture_chip(ax_b, stage_x[1] - cw / 2, yy - ch / 2, w=cw, h=ch)
    ax_b.text(stage_x[1], yy - 0.27, "Visium capture\n(10x Genomics)",
              ha="center", fontsize=6.2, color="#444")

    _arrow(ax_b, stage_x[1] + cw / 2 + 0.01, yy, stage_x[2] - 0.10, yy, mut=7)

    # 3) spots x genes matrix glyph
    mw, mh = 0.16, 0.24
    mat_x, mat_y = stage_x[2] - mw / 2, yy - mh / 2
    ax_b.add_patch(mpatches.Rectangle((mat_x, mat_y), mw, mh,
                                      facecolor="white",
                                      edgecolor=C_AX, linewidth=0.7))
    cols, rows = 6, 8
    for i in range(cols):
        for j in range(rows):
            shade = 0.95 - 0.13 * ((i * 2 + j) % 5)
            ax_b.add_patch(mpatches.Rectangle(
                (mat_x + 0.005 + i * (mw - 0.01) / cols,
                 mat_y + 0.005 + j * (mh - 0.01) / rows),
                (mw - 0.01) / cols * 0.92, (mh - 0.01) / rows * 0.92,
                facecolor=str(shade), edgecolor="none"))
    ax_b.text(stage_x[2], yy - 0.27, "spots x genes",
              ha="center", fontsize=6.2, color="#444")

    _arrow(ax_b, stage_x[2] + mw / 2 + 0.01, yy, stage_x[3] - 0.10, yy, mut=7)

    # 4) Allen atlas-coloured brain
    _brain(ax_b, stage_x[3], yy, scale=0.55, fill="#eef3e8")
    region_cols = ["#9bbed6", "#f6b26b", "#a5c4a3", "#d8a4d4"]
    for i, c in enumerate(region_cols):
        wedge = mpatches.Wedge((stage_x[3], yy), 0.085,
                               90 + i * 90, 90 + (i + 1) * 90,
                               width=0.035, facecolor=c,
                               edgecolor="none", alpha=0.9)
        ax_b.add_patch(wedge)
    ax_b.text(stage_x[3], yy - 0.27, "Allen atlas\nregistration",
              ha="center", fontsize=6.2, color="#444")

    # --- C: analysis pipeline ---------------------------------------------
    ax_c = fig.add_subplot(gs[1, :])
    ax_c.set_xlim(0, 1); ax_c.set_ylim(0, 1)
    _strip(ax_c); _label(ax_c, "C")
    ax_c.set_title("Analysis modules", loc="left", pad=2, fontsize=7.5)

    modules = [
        ("Disease\nbaseline", "PBS vs WT\nDE & PIG signature", "#e6eef7", C_PBS),
        ("Global\nrescue", "BRI vs WT, BRI vs PBS\nattenuation", "#f7e6e6", C_BRI),
        ("Regional\nheterogeneity", "DE per region\nrescue ranking", "#ecf3ec", "#5a8c5a"),
        ("Spatial\ngradients", "meninges to parenchyma\nPIG flattening", "#f5efe1", "#b8862c"),
        ("Cell-type\nbasis", "tangram deconv.\nmicroglia, endo.", "#efe7f3", "#7e57a8"),
    ]
    n = len(modules)
    bw = 0.165
    gap = (1.0 - n * bw) / (n + 1)
    for i, (title, sub, fill, col) in enumerate(modules):
        x0 = gap + i * (bw + gap)
        box = FancyBboxPatch((x0, 0.18), bw, 0.66,
                             boxstyle="round,pad=0.006,rounding_size=0.012",
                             facecolor=fill, edgecolor=col, linewidth=0.9)
        ax_c.add_patch(box)
        ax_c.text(x0 + bw / 2, 0.66, title, ha="center", va="center",
                  fontsize=7.4, fontweight="bold", color=col)
        ax_c.text(x0 + bw / 2, 0.40, sub, ha="center", va="center",
                  fontsize=6.3, color="#444")
        # tag with figure number
        ax_c.text(x0 + bw / 2, 0.10, f"Fig {i + 1}", ha="center",
                  fontsize=6.3, color="#666", style="italic")
        if i < n - 1:
            _arrow(ax_c, x0 + bw + 0.005,
                   0.51, x0 + bw + gap - 0.005, 0.51,
                   color="#888", lw=0.7, mut=6)

    save(fig, "figure0_schematic")


# ---------- data figures -----------------------------------------------------

def figure1_disease_baseline():
    fig = plt.figure(figsize=(W_DOUBLE, 5.0))
    gs = GridSpec(2, 3, figure=fig, hspace=0.30, wspace=0.20,
                  height_ratios=[1, 1])
    _panel(fig.add_subplot(gs[0, 0]),
           NB / "nb00_volcano_pbs_vs_wt.png",
           "A", "PBS vs WT — global DE")
    _panel(fig.add_subplot(gs[0, 1]),
           NB / "nb04_spatial_regions.png",
           "B", "Region annotation")
    _panel(fig.add_subplot(gs[0, 2]),
           NB / "nb04_umap_regions.png",
           "C", "UMAP, by region")
    _panel(fig.add_subplot(gs[1, 0]),
           NB / "nb04_spatial_pig.png",
           "D", "PIG signature in space")
    _panel(fig.add_subplot(gs[1, 1]),
           NB / "nb04_pig_by_region.png",
           "E", "PIG score per region")
    _panel(fig.add_subplot(gs[1, 2]),
           NB / "nb04_region_composition.png",
           "F", "Regional composition")
    save(fig, "figure1_disease_baseline")


def figure2_global_rescue():
    fig = plt.figure(figsize=(W_DOUBLE, 4.6))
    gs = GridSpec(2, 3, figure=fig, hspace=0.30, wspace=0.22)
    _panel(fig.add_subplot(gs[0, 0]),
           NB / "nb01_volcano_bri_vs_wt.png",
           "A", "BRICHOS vs WT")
    _panel(fig.add_subplot(gs[0, 1]),
           NB / "nb02_volcano_bri_vs_pbs.png",
           "B", "BRICHOS vs PBS")
    _panel(fig.add_subplot(gs[0, 2]),
           NB / "nb03_lfc_correlation.png",
           "C", "LFC concordance")
    _panel(fig.add_subplot(gs[1, 0]),
           NB / "nb05_lfc_comparison.png",
           "D", "LFC: BRI vs PBS comparison")
    _panel(fig.add_subplot(gs[1, 1]),
           NB / "nb05_attenuation_scatter.png",
           "E", "Per-gene attenuation")
    _panel(fig.add_subplot(gs[1, 2]),
           NB / "nb03_stacked_summary.png",
           "F", "DEG class summary")
    save(fig, "figure2_global_rescue")


def figure3_regional_heterogeneity():
    fig = plt.figure(figsize=(W_DOUBLE, 6.2))
    gs = GridSpec(3, 3, figure=fig, hspace=0.32, wspace=0.20,
                  height_ratios=[1.05, 1, 1])
    # top: panoramic regional summary
    _panel(fig.add_subplot(gs[0, :]),
           FIG / "regional_BRICHOS_vs_PBS_summary.png",
           "A", "BRICHOS vs PBS — regional summary")
    _panel(fig.add_subplot(gs[1, 0]),
           NB / "nb04_de_summary_barplot.png",
           "B", "DEG count per region")
    _panel(fig.add_subplot(gs[1, 1]),
           NB / "nb04_sig_genes_region.png",
           "C", "Significant gene ranks")
    _panel(fig.add_subplot(gs[1, 2]),
           NB / "nb05_pig_barplot.png",
           "D", "PIG response per region")
    _panel(fig.add_subplot(gs[2, 0]),
           NB / "nb04_heatmap_regions.png",
           "E", "Region × gene heatmap")
    _panel(fig.add_subplot(gs[2, 1]),
           NB / "nb04_regional_summary.png",
           "F", "Regional rescue overview")
    _panel(fig.add_subplot(gs[2, 2]),
           NB / "nb03_heatmap.png",
           "G", "Gene-set heatmap")
    save(fig, "figure3_regional_heterogeneity")


def figure4_spatial_gradient():
    fig = plt.figure(figsize=(W_DOUBLE, 6.0))
    gs = GridSpec(3, 3, figure=fig, hspace=0.28, wspace=0.18)
    _panel(fig.add_subplot(gs[0, 0]),
           NB / "nb07_tissue_rim.png",
           "A", "Tissue rim / meninges")
    _panel(fig.add_subplot(gs[0, 1]),
           NB / "nb07_spatial_pig.png",
           "B", "PIG signature in space")
    _panel(fig.add_subplot(gs[0, 2]),
           NB / "nb07_umap_regions.png",
           "C", "UMAP, regional embedding")
    _panel(fig.add_subplot(gs[1, 0]),
           NB / "nb07_pig_gradient_pbs.png",
           "D", "PIG gradient — PBS")
    _panel(fig.add_subplot(gs[1, 1]),
           NB / "nb07_pig_gradient_bri.png",
           "E", "PIG gradient — BRICHOS")
    _panel(fig.add_subplot(gs[1, 2]),
           NB / "nb07_spatial_dist_meninges.png",
           "F", "Distance from meninges")
    _panel(fig.add_subplot(gs[2, 0]),
           NB / "nb07_pig_compartment_pbs.png",
           "G", "Compartment quant — PBS")
    _panel(fig.add_subplot(gs[2, 1]),
           NB / "nb07_pig_compartment_bri.png",
           "H", "Compartment quant — BRICHOS")
    save(fig, "figure4_spatial_gradient")


def figure5_celltype():
    """Cell-type basis — placeholder until tangram outputs are exported."""
    fig = plt.figure(figsize=(W_DOUBLE, 4.6))
    gs = GridSpec(2, 3, figure=fig, hspace=0.28, wspace=0.20)
    candidates = [
        (NB / "nb06_celltype_proportions.png",   "A", "Cell-type proportions"),
        (NB / "nb06_celltype_region_heatmap.png", "B", "Cell-type × region"),
        (NB / "nb06_microglia_rescue.png",        "C", "Microglia rescue"),
        (NB / "nb06_endocytosis_module.png",      "D", "Endocytosis module"),
        (NB / "nb06_tangram_overlay.png",         "E", "Tangram overlay"),
        (NB / "nb06_celltype_dge.png",            "F", "Cell-type DGE"),
    ]
    for ax_pos, (p, lab, title) in zip(
        [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)], candidates
    ):
        _panel(fig.add_subplot(gs[ax_pos]), p, lab, title)
    save(fig, "figure5_celltype")


def main():
    figure0_schematic()
    figure1_disease_baseline()
    figure2_global_rescue()
    figure3_regional_heterogeneity()
    figure4_spatial_gradient()
    figure5_celltype()


if __name__ == "__main__":
    main()
