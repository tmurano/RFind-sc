#!/usr/bin/env python3
"""Supplementary Fig. S2 — choice of the per-cell statistic (Delta): panels b-e.

  a  from UMI counts to a per-cell ranked list, in the design of Fig. 1c, with no formula of any candidate (drawn
     in the slide, not here)
  b  the six candidates: a table of definitions, advantages and drawbacks, every cell on one line, 162 mm wide, 6 pt
     (FigS2b_statistics_table.png). P2 is not marked: which candidate is adopted follows from the comparison in d, e.
  c  each candidate over the genes of one example cell, labelled "One example cell" (FigS2c_one_cell.png, 151.5 x 39.3 mm):
     the genes split first by UMI = 0 or >= 1 (two rows), then by above or below the expected count; one key block on the left
     (heading 7 pt bold, the number of genes in grey, the two groups as 7-pt bold sub-headings, each entry on one line at
     6 pt); the six plots with the same width and the same gaps, one y axis title ("Number of genes", on one line) for both rows and tick labels on the
     first plot of each row only; simpleFC on a log axis (0, 1, 100, with unlabelled ticks at 0.1 and 10) with its zeros as a
     separate bar at the left; the z-score and P1-P3 on one axis (PLO..PHI), the genes outside it drawn as bars at the two
     ends (no gene is dropped), with small break marks on the axis
  d  the schematic of the check (FigS2d_check_schematic.png, 68 x 56.2 mm): 1 the clusters of Hammond et al. are the answer,
     2 UP / DOWN of each cell, 3 check against the markers of its cluster
  e  the result (FigS2e_benchmark.png, 105 x 56.2 mm = the height of d, with a 3-mm white margin on the left): 4 the two shares averaged over the
     13 clusters and 5 their sum for each cluster (points of one cluster joined across the candidates); the candidates
     ordered by the composite, P2 in black; legend 7 pt inside the plot, clear of the bars
Every list length is 500, as in the pipeline. The share of single-molecule genes among the top 500 of a cell is written to
output/data/figS2_table_check.csv (not drawn).

Style: Arial; tables as in Fig. 4c (header bar, one rule, zebra rows, no grid); red = UP / top of the list, blue = DOWN / bottom
of the list (paper-wide rule); candidates themselves are neutral grey with P2 in black in the results (e), so red and blue keep a
single meaning. No panel letters (they are added in the final figure). Each panel is drawn at the size at which it is pasted,
600 dpi, axes placed in mm, no bbox_inches="tight", with the text sizes of all figures: tick labels 6 pt, axis titles, legends
and labels 7 pt, dense text (table, schematic notes, histogram keys) 6 pt; axes and ticks 0.5 pt.

Input : figs/Supp_S2_Delta/output/data/{figS2_delta_bench.csv, figS2_delta_bench_percluster.csv,
        figS2b_one_cell.csv (all three by figS2_delta_benchmark.R)}
Output: figs/Supp_S2_Delta/output/FigS2{b_statistics_table, c_one_cell, d_check_schematic, e_benchmark}.png (PNG only)
Run   : python figs/Supp_S2_Delta/analysis/figS2_panels.py   (project root; after the two R scripts)
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.patches import Rectangle

if any(f.name == "Arial" for f in font_manager.fontManager.ttflist):
    matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["mathtext.fontset"] = "custom"
matplotlib.rcParams["mathtext.it"] = "Arial:italic"
matplotlib.rcParams["mathtext.rm"] = "Arial"
## sizes at the pasted size: text in pt, lengths in mm
FT, FL, FN = 6, 7, 6                                                    # tick labels / axis titles, legends, labels / dense notes
matplotlib.rcParams.update({"font.size": FL, "axes.labelsize": FL, "xtick.labelsize": FT, "ytick.labelsize": FT,
                            "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                            "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
                            "xtick.major.size": 2.0, "ytick.major.size": 2.0, "xtick.major.pad": 1.2, "ytick.major.pad": 1.2,
                            "axes.labelpad": 1.5, "axes.titlepad": 2.0, "legend.fontsize": FL})
MM, DPI = 1 / 25.4, 600

OUT = "figs/Supp_S2_Delta/output"; DAT = os.path.join(OUT, "data")
RED, BLUE = "#C0392B", "#1F5FA8"
HEAD, ZEBRA, INK, MUTED, P2C, OTHER = "#222222", "#F2F2F2", "#222222", "#8A8A8A", "#111111", "#B8B8B8"
LABEL = {"P2": "P2", "P1": "P1", "log2diff": "log2diff", "P3_NB": "P3 (NB)", "simpleFC": "simpleFC", "zscore": "z-score"}

bench = pd.read_csv(os.path.join(DAT, "figS2_delta_bench.csv")).set_index("model")
ORDER = bench.sort_values("composite", ascending=False).index.tolist()   # e: by the composite of the benchmark
perc = pd.read_csv(os.path.join(DAT, "figS2_delta_bench_percluster.csv"))
one = pd.read_csv(os.path.join(DAT, "figS2b_one_cell.csv"))


def mm_fig(w, h):
    """A figure of w x h mm (rounded to whole pixels; + 0.01 px so that the pixel count is not truncated)."""
    fig = plt.figure(figsize=((round(w * MM * DPI) + 0.01) / DPI, (round(h * MM * DPI) + 0.01) / DPI), dpi=DPI)
    fig.W, fig.H = w, h
    return fig


def mm_axes(fig, x, y, w, h, **kw):
    """Axes placed in mm from the lower left corner of the figure."""
    return fig.add_axes([x / fig.W, y / fig.H, w / fig.W, h / fig.H], **kw)


def save(fig, name):
    fig.savefig(os.path.join(OUT, f"{name}.png"), dpi=DPI, facecolor="white")   # PNG only
    plt.close(fig)
    print("wrote", os.path.join(OUT, f"{name}.png"))


_mfig = plt.figure(dpi=DPI)                                             # for measuring texts


def text_wh(s, **kw):
    """Width and height of a text in mm."""
    t = _mfig.text(0, 0, s, **kw); b = t.get_window_extent(_mfig.canvas.get_renderer()); t.remove()
    return b.width / DPI * 25.4, b.height / DPI * 25.4


def badge(ax, x, y, n):
    """A step number in a dark disc (2.7 mm across; a marker, so it stays round when the axes are not to scale)."""
    ax.scatter([x], [y], s=(2.7 * MM * 72) ** 2, color="#333333", zorder=5, clip_on=False, linewidths=0)
    ax.text(x, y, str(n), ha="center", va="center", color="white", fontsize=FN, fontweight="bold", zorder=6)


# ------------------------------------------------------------------ b: the statistics (table, Fig. 4c style), every cell on one line
rows = [
    ("simpleFC", r"CPM / mean CPM", "Same as the bulk fold change",
     "All zeros tie at 0; one molecule of a rare gene gives a huge value"),
    ("log2diff", r"$\log_2$(CPM + 1) − $\log_2$(mean CPM + 1)", "Widely used; no −∞",
     "Zeros fall far below 0, below detected genes"),
    ("z-score", r"$(y_{ij} - \bar{y}_i)\,/\,s_i$", "Standard in Seurat; zeros of rare genes near 0",
     "One molecule of a rare gene gives a large value (up to 10)"),
    ("P1", r"$(x_{ij} - e_{ij})\,/\,\sqrt{e_{ij}}$", "Zeros ordered by $e$; near 0 when $e$ is small",
     "One molecule at a small $e$ gives a huge value"),
    ("P2", r"$(x_{ij} - e_{ij})\,/\,\sqrt{e_{ij} + 1}$", "As P1; one molecule stays below 1",
     "Noise of highly expressed genes can rank high"),
    ("P3 (NB)", r"$(x_{ij} - e_{ij})\,/\,\sqrt{e_{ij} + e_{ij}^2/\theta_i}$", "As P1; noise of highly expressed genes is damped",
     "One molecule at a small $e$ gives a huge value; θ for each gene"),
]
cols = ["Statistic", "Definition", "Advantages", "Drawbacks"]
TFS, PAD, VPAD = FN, 1.0, 0.15                                          # 6 pt; padding left/right of a cell, above/below its text (mm)
HROW, MINROW = 3.3, 3.8                                                 # header height, least row height (mm)
wb = lambda k: "bold" if k == 0 else "normal"
need = [max(max(text_wh(r[k], fontsize=TFS, fontweight=wb(k))[0] for r in rows),
            text_wh(cols[k], fontsize=TFS, fontweight="bold")[0]) + 2 * PAD for k in range(4)]
TW = max(160.0, math.ceil(sum(need) * 2) / 2)                           # about 160 mm, wider if one line needs it
colw = [w + (TW - sum(need)) / 4 for w in need]
rh = [max(MINROW, max(text_wh(v, fontsize=TFS, fontweight=wb(k))[1] for k, v in enumerate(r)) + 2 * VPAD) for r in rows]
base = max(h for h in rh if h <= 1.2 * MINROW)                          # one height for the rows of plain text
rh = [max(base, h) for h in rh]                                         # only the row of P3 (its tall formula) is taller
TH = HROW + sum(rh)
print(f"table: {TW:.1f} x {TH:.1f} mm; columns {', '.join(f'{c:.1f}' for c in colw)}; rows {', '.join(f'{h:.2f}' for h in rh)}")
fig = mm_fig(TW, TH); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, TW); ax.set_ylim(0, TH); ax.axis("off")
xs = np.concatenate([[0], np.cumsum(colw)])
ax.add_patch(Rectangle((0, TH - HROW), TW, HROW, facecolor=HEAD, edgecolor="none"))
for k, c in enumerate(cols):
    ax.text(xs[k] + PAD, TH - HROW / 2, c, ha="left", va="center", color="white", fontsize=TFS, fontweight="bold")
ax.plot([0, TW], [TH - HROW, TH - HROW], color="#333333", lw=0.5)
y1 = TH - HROW
for i, r in enumerate(rows):
    y0 = y1 - rh[i]
    if i % 2 == 0:
        ax.add_patch(Rectangle((0, y0), TW, rh[i], facecolor=ZEBRA, edgecolor="none"))
    for k, v in enumerate(r):
        ax.text(xs[k] + PAD, (y0 + y1) / 2, v, ha="left", va="center", fontsize=TFS, fontweight=wb(k), color=INK)
    y1 = y0
save(fig, "FigS2b_statistics_table")

# ------------------------------------------------------------------ c: each statistic over the genes of one example cell
## Genes split first by UMI = 0 or >= 1 (rows), then (colours): UMI >= 1 = share at or above the mean share (UP) or below it
## (DOWN-1); UMI = 0 = expected count < ln 2 (a zero is more likely than not, DOWN-2) or >= ln 2 (DOWN-3) (figS2_delta_benchmark.R).
## Two rows. Counts on a square-root axis. simpleFC on a log axis, with
## its zeros as a separate bar at the left (a ratio of zero has no logarithm). P1-P3 share the axis PLO..PHI (so that extreme
## outliers do not stretch it); the genes outside it are drawn as the bars at the two ends, so no gene is dropped.
KCOL = {"UP": RED, "DOWN-1": "#9ECAE1", "DOWN-2": "#BDBDBD", "DOWN-3": BLUE}
GROUPS = [("UMI ≥ 1", ["DOWN-1", "UP"], {"UP": "above expected", "DOWN-1": "below expected"}),
          ("UMI = 0", ["DOWN-3", "DOWN-2"], {"DOWN-2": "expected count < ln 2", "DOWN-3": "expected count ≥ ln 2"})]
TORDER = ["simpleFC", "log2diff", "zscore", "P1", "P2", "P3_NB"]       # the order of the table
PLO, PHI, PBW, POUT = -4.0, 8.0, 0.1, 1.1
YMAX = [200, 8000]; YT = [[0, 50, 100, 200], [0, 1000, 4000, 8000]]    # fewer labels in the lower row
## layout in mm (size in the slide 151.5 x 39.3): key block | y axis | six plots
HW, HH = 151.5, 39.3
XP0, XP1, GAPc = 35.8, 150.5, 3.6                                       # left edge of the first plot, right edge of the last, gap
PWc = (XP1 - XP0 - 5 * GAPc) / 6                                        # one width for all six candidates
RB, RG, RT = 5.9, 2.4, HH - 3.4                                         # bottom of the lower row, gap between rows, top of the upper row
RHc = (RT - RB - RG) / 2
fig = mm_fig(HW, HH)
ov = fig.add_axes([0, 0, 1, 1], zorder=10); ov.set_xlim(0, HW); ov.set_ylim(0, HH); ov.axis("off"); ov.patch.set_visible(False)


def break_mark(xm, ym):
    """Axis-break mark at (xm, ym) mm: a short gap in the axis line and two small slanted strokes."""
    ov.add_patch(Rectangle((xm - 0.3, ym - 0.25), 0.6, 0.5, facecolor="white", edgecolor="none", zorder=1))
    for dx in (-0.3, 0.3):
        ov.plot([xm + dx - 0.25, xm + dx + 0.25], [ym - 0.55, ym + 0.55], color="#555555", lw=0.5, solid_capstyle="butt", zorder=2)


nout = {}
for c, m in enumerate(TORDER):
    xl = XP0 + c * (PWc + GAPc)
    v_all = one[m].values; k_all = one["kind"].values
    for r, (glab, kinds, klab) in enumerate(GROUPS):
        yb = RT - RHc if r == 0 else RB
        ax = mm_axes(fig, xl, yb, PWc, RHc)
        sel = np.isin(k_all, kinds); v = v_all[sel]; k = k_all[sel]
        if m == "simpleFC":
            lvall = np.log10(v_all[v_all > 0]); edges = np.linspace(lvall.min(), lvall.max(), 61)
            x0 = edges[0] - 0.75; w0 = (edges[1] - edges[0]) * 2.2          # the bar of zeros, left of a gap
            bottom = 0.0
            for kk in kinds:
                n = ((v == 0) & (k == kk)).sum(); ax.bar(x0, n, width=w0, bottom=bottom, color=KCOL[kk], linewidth=0); bottom += n
            pos = v > 0; bottom = np.zeros(len(edges) - 1)
            for kk in kinds:
                n, _ = np.histogram(np.log10(v[pos & (k == kk)]), bins=edges)
                ax.bar(edges[:-1], n, width=np.diff(edges), bottom=bottom, align="edge", color=KCOL[kk], linewidth=0); bottom += n
            ax.axvline(0, color="#777777", lw=0.4, ls=(0, (3, 2)), zorder=0)
            ax.set_xticks([x0, 0, 2]); ax.set_xticklabels(["0", "1", "100"])
            ax.set_xticks([-1, 1], minor=True)                         # 0.1 and 10, unlabelled
            ax.set_xlim(x0 - w0, edges[-1] + 0.1)
            brk = [(x0 + w0 / 2 + edges[0]) / 2]
        elif m == "log2diff":
            edges = np.arange(-11.5, 7.51, 0.25); bottom = np.zeros(len(edges) - 1)
            for kk in kinds:
                n, _ = np.histogram(v[k == kk], bins=edges)
                ax.bar(edges[:-1], n, width=np.diff(edges), bottom=bottom, align="edge", color=KCOL[kk], linewidth=0); bottom += n
            ax.axvline(0, color="#777777", lw=0.4, ls=(0, (3, 2)), zorder=0)
            ax.set_xlim(edges[0] - 0.3, edges[-1] + 0.3); ax.set_xticks([-10, -5, 0, 5]); brk = []
        else:
            edges = np.arange(PLO, PHI + PBW / 2, PBW); bottom = np.zeros(len(edges) - 1)
            for kk in kinds:
                n, _ = np.histogram(v[k == kk], bins=edges)
                ax.bar(edges[:-1], n, width=np.diff(edges), bottom=bottom, align="edge", color=KCOL[kk], linewidth=0); bottom += n
            for side, xo, test in (("lo", PLO - POUT, v < PLO), ("hi", PHI + POUT, v > PHI)):   # genes outside the axis
                b0 = 0.0
                for kk in kinds:
                    n = (test & (k == kk)).sum(); ax.bar(xo, n, width=0.45, bottom=b0, color=KCOL[kk], linewidth=0); b0 += n
                nout[(m, side)] = nout.get((m, side), 0) + int(b0)
            ax.axvline(0, color="#777777", lw=0.4, ls=(0, (3, 2)), zorder=0)
            ax.set_xlim(PLO - POUT - 0.5, PHI + POUT + 0.5)
            ax.set_xticks([PLO - POUT, 0, 4, PHI + POUT]); ax.set_xticklabels([f"<{PLO:g}".replace("-", "−"), "0", "4", f">{PHI:g}"])
            brk = [PLO - POUT / 2, PHI + POUT / 2]
        lo_, hi_ = ax.get_xlim()
        for bx in brk:
            break_mark(xl + (bx - lo_) / (hi_ - lo_) * PWc, yb)
        ax.set_yscale("function", functions=(np.sqrt, np.square)); ax.set_ylim(0, YMAX[r])
        ax.set_yticks(YT[r]); ax.set_yticklabels([f"{t:,}" for t in YT[r]] if c == 0 else [])
        if r == 0:
            ax.set_title(LABEL[m], fontsize=FL, fontweight="bold", pad=2.0); ax.tick_params(axis="x", which="both", labelbottom=False)
        ax.tick_params(which="major", length=1.5, pad=1.0); ax.tick_params(which="minor", length=1.0)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
# key block on the left: the heading, then each row's group next to its row
axl = mm_axes(fig, 0, 0, XP0 - 13.0, HH); axl.set_xlim(0, XP0 - 13.0); axl.set_ylim(0, HH); axl.axis("off")
axl.text(0.0, HH - 0.3, "One example cell", ha="left", va="top", fontsize=FL, fontweight="bold")
axl.text(0.0, HH - 3.0, "(8,960 genes)", ha="left", va="top", fontsize=FN, color="#555555")
for r, (glab, kinds, klab) in enumerate(GROUPS):
    yc = (RT - RHc / 2) if r == 0 else (RB + RHc / 2)                   # centre of the row
    axl.text(0.0, yc + 2.8, glab, ha="left", va="center", fontsize=FL, fontweight="bold")
    for j, kk in enumerate(kinds[::-1]):
        yy_ = yc - j * 2.8
        axl.add_patch(Rectangle((0.0, yy_ - 0.8), 1.6, 1.6, facecolor=KCOL[kk], edgecolor="none"))
        axl.text(2.2, yy_, klab[kk], ha="left", va="center", fontsize=FN)
fig.text((XP0 + XP1) / 2 / HW, 0.3 / HH, "Value of each gene", ha="center", va="bottom", fontsize=FL)
fig.text((XP0 - 8.2) / HW, (RB + RT) / 2 / HH, "Number of genes", rotation=90, ha="center", va="center", fontsize=FL)   # one y title for both rows, on one line
print("genes outside the P axis:", nout)
save(fig, "FigS2c_one_cell")

# ------------------------------------------------------------------ d: the schematic of the check (steps 1-3)
## list length 500 everywhere (cell lists, cluster markers, references), as in the pipeline
## The drawing is given in its own coordinates (a 55 x 50 box), spread to
## the size of the panel by the axes; text and lines keep their size in pt.
DW, DT = 68.0, 43.2                                                     # width in the slide (mm); top of the drawing (units)
DH = DT * 65.0 / 50.0                                                   # height in the slide (mm): 65 mm per 50 units of the box
fig = mm_fig(DW, DH)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 55.0); ax.set_ylim(0, DT); ax.axis("off")
X0, XU, XD, X1, BH = 26.5, 35.5, 45.5, 54.5, 3.6                        # bar: UP / markers [X0, XU], DOWN / negative markers [XD, X1]
rng = np.random.default_rng(3)


def bar(y, left, right, lab, ticks=None):
    ax.add_patch(Rectangle((X0, y - BH / 2), X1 - X0, BH, facecolor="white", edgecolor="#555555", lw=0.5))
    ax.add_patch(Rectangle((X0, y - BH / 2), XU - X0, BH, **left))
    ax.add_patch(Rectangle((XD, y - BH / 2), X1 - XD, BH, **right))
    ax.text(X0 - 0.8, y, lab, ha="right", va="center", fontsize=FN)
    if ticks:                                                        # genes of the cell's list that are markers of its cluster
        for (lo, hi, frac, col) in ticks:
            for xx in lo + (hi - lo) * np.sort(rng.random(int(round(frac * 14)))):
                ax.plot([xx, xx], [y - BH / 2, y + BH / 2], color=col, lw=0.6, solid_capstyle="butt")


XL = 3.6                                                                # left edge of the step names
# 1 the answer: the markers of the cell's cluster
ya = 38.0
badge(ax, 1.4, ya, 1); ax.text(XL, ya, "Answer", ha="left", va="center", fontsize=FL, fontweight="bold")
ax.text(X1, ya + BH / 2 + 1.0, "Hammond et al.: cluster k vs. other cells", ha="right", va="bottom", fontsize=FN, color="#555555")
bar(ya, dict(facecolor=RED, edgecolor="none"), dict(facecolor=BLUE, edgecolor="none"), "cluster k")
ax.text((X0 + XU) / 2, ya, "markers", ha="center", va="center", color="white", fontsize=FN, fontweight="bold")
ax.text((XD + X1) / 2, ya, "negative\nmarkers", ha="center", va="center", color="white", fontsize=FN, fontweight="bold",
        linespacing=0.95)
# 2 the lists of each cell of cluster k (by each candidate)
yc = [27.0, 22.0, 17.0]
badge(ax, 1.4, yc[1], 2); ax.text(XL, yc[1] + 0.3, "Each cell", ha="left", va="bottom", fontsize=FL, fontweight="bold")
ax.text(XL, yc[1] - 0.2, "each candidate", ha="left", va="top", fontsize=FN, color="#555555")
ax.text((X0 + XU) / 2, yc[0] + BH / 2 + 0.5, "UP (500)", ha="center", va="bottom", fontsize=FN)
ax.text(X1, yc[0] + BH / 2 + 0.5, "DOWN (500)", ha="right", va="bottom", fontsize=FN)
for k, y in enumerate(yc):
    bar(y, dict(facecolor=RED, alpha=0.22, edgecolor="none"), dict(facecolor=BLUE, alpha=0.22, edgecolor="none"), f"cell {k + 1}",
        ticks=[(X0 + 0.3, XU - 0.3, 0.45 + 0.1 * (k == 1), RED), (XD + 0.3, X1 - 0.3, 0.8, BLUE)])
ax.scatter([X0 - 4.5] * 3, [yc[-1] - 3.4, yc[-1] - 4.4, yc[-1] - 5.4], s=0.8, color="#333333", linewidths=0)   # more cells
# 3 check: dotted guides from the answer to the lists, and what is counted
for xx in ((X0 + XU) / 2, (XD + X1) / 2):
    ax.plot([xx, xx], [ya - BH / 2 - 0.4, yc[0] + BH / 2 + 3.2], color="#777777", lw=0.5, ls=(0, (2, 2)))
yk = 8.2
badge(ax, 1.4, yk, 3); ax.text(XL, yk, "Check", ha="left", va="center", fontsize=FL, fontweight="bold")
ax.text((X0 + XU) / 2 - 1.2, yk, "share of UP\nin markers", ha="center", va="center", fontsize=FN, color=RED, linespacing=1.05)
ax.text(X1, yk, "share of DOWN in\nnegative markers", ha="right", va="center", fontsize=FN, color=BLUE, linespacing=1.05)
save(fig, "FigS2d_check_schematic")

# ------------------------------------------------------------------ e: the shares, mean of the 13 clusters (4) and each cluster (5)
EW, EH = 105.0, 56.2                                                    # size in the slide (mm); the height of d
ML = 3.0                                                                # white margin on the left (space between d and e)
fig = mm_fig(EW, EH)
YB, YT = 6.6, EH - 7.0                                                  # bottom and top of the two plots (mm)
YLIM = (-1.0, 5.6)                                                      # room for the legend under the z-score bars
AX1, AW1 = ML + 10.4, 42.0                                              # left edge and width of the two plots (mm); 10.4 = the names
AX2 = AX1 + AW1 + 13.2
AW2 = EW - 1.6 - AX2
yy = np.arange(len(ORDER))[::-1]


def names(ax):
    ax.set_yticks(yy); ax.set_yticklabels([LABEL[m] for m in ORDER])
    for t, m in zip(ax.get_yticklabels(), ORDER):
        t.set_fontweight("bold" if m == "P2" else "normal"); t.set_color(P2C if m == "P2" else "#555555")
    ax.set_ylim(*YLIM)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def title(x, n, text):
    """Badge and title above a plot; x = left edge of the badge (mm)."""
    tax = mm_axes(fig, x, YT + 2.2, 40.0, 3.4); tax.set_xlim(0, 40.0); tax.set_ylim(-1.7, 1.7); tax.axis("off")
    badge(tax, 1.35, 0, n)
    tax.text(3.4, 0, text, ha="left", va="center", fontsize=FL, fontweight="bold")


ax = mm_axes(fig, AX1, YB, AW1, YT - YB)
h = 0.36
ax.barh(yy + h / 2, [bench.loc[m, "marker_top"] for m in ORDER], height=h, color=RED, label="markers in UP")
ax.barh(yy - h / 2, [bench.loc[m, "negmarker_bot"] for m in ORDER], height=h, color=BLUE, label="negative markers\nin DOWN")
names(ax)
ax.set_xlabel("Share"); ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 0.02), handlelength=1.2, handletextpad=0.5,
          labelspacing=0.35, borderaxespad=0, borderpad=0)               # right of the simpleFC bars, under the z-score bars
title(ML, 4, "Mean of the 13 clusters")
# the points of one cluster are joined across the candidates, so that the candidates can be compared within each
# cluster (Supplementary Note 2)
ax = mm_axes(fig, AX2, YB, AW2, YT - YB)
Wc = perc.pivot(index="cluster", columns="model", values="composite")[ORDER]
for _, r in Wc.iterrows():
    ax.plot(r.values, yy, color="#C4C4C4", lw=0.45, zorder=1)
for k, m in enumerate(ORDER):
    v = Wc[m].values
    y = yy[k]
    col = P2C if m == "P2" else "#9A9A9A"
    ax.plot([np.percentile(v, 25), np.percentile(v, 75)], [y, y], color=col, lw=3.6, solid_capstyle="butt", alpha=0.35, zorder=2)
    ax.scatter(v, np.full(len(v), y), s=8, color=col, zorder=3, edgecolor="white", linewidth=0.3)
    ax.plot([np.median(v)] * 2, [y - 0.3, y + 0.3], color=col, lw=1.0, zorder=4)
names(ax)
ax.set_xticks([0.2, 0.4, 0.6])
ax.set_xlabel("Share in UP + share in DOWN")
title(AX2 - 10.4, 5, "Each of the 13 clusters")
save(fig, "FigS2e_benchmark")
plt.close(_mfig)
