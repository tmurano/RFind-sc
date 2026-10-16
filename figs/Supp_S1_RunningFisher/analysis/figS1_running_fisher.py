#!/usr/bin/env python3
"""Supplementary Fig. S1 — the Running Fisher test (Kupershmidt et al., 2010), the basis of RFind-sc.

The panels follow the order of Supplementary Note 1:

  a  the two ranked gene lists (titled "Gene list b1" / "Gene list b2"), each drawn as a table of genes ordered by the change,
     with a row of dots for the genes in between; b1_up / b1_down are marked on the table with arrows that run from each end
     inwards (the strongest change first); the four directional pairs between the lists
  b  one pair (b1_up x b2_up) at one pair of cut-offs: the top G and top K genes (G and K marked as movable by double arrows),
     their overlap M and the 2 x 2 table, drawn at about the size of the lists, each cell with its symbol and the number of
     this example (G = 4, K = 5, M = 3, U = 26), and the P of this test
  c  every pair of cut-offs tried for that pair (one square = one test; colour key "-log10 P of the test at each (G, K)";
     arrows along both axes; a thin frame around the map): the dashed square is the test of b, and the smallest P is kept as
     the representative P value of the pair
  d  the four scores (-log10 of the representative P values; small bars) and their signed sum, S_total, written in words and
     drawn as in Fig. 1d (a line from opposite to same direction, 0 in the middle, the value of this example marked; no numbers)

It is a SCHEMATIC. The gene lists are toy lists: each ranks all 26 genes of the universe (GeneA-GeneZ), 13 up and 13 down, so
the map of c is 13 x 13 (b draws the top 8 and dots). Every number behind the drawing (M, P, the heat map, the scores, S_total)
is computed from these toy lists with the one-sided Fisher's exact test and a gene universe of U = 26, so the panels agree
with each other. P values are a capital P, as in the paper. Panel titles say what each panel shows.
Style follows the main figures: Arial, header bars with white bold text and zebra rows (Fig. 4c), no outer frames,
red #C0392B = same direction, blue #1F5FA8 = opposite direction, no green; mechanical markers are neutral grey.
Run   : python figs/Supp_S1_RunningFisher/analysis/figS1_running_fisher.py   (from the project root)
Output: figs/Supp_S1_RunningFisher/output/FigS1_running_fisher.png   (pasted size, 180 x 125.6 mm, 600 dpi) and each panel
        on its own, FigS1{a_gene_lists, b_one_cutoff, c_every_cutoff, d_signed_sum}.png (see the end of the script)
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy.stats import hypergeom

if any(f.name == "Arial" for f in font_manager.fontManager.ttflist):
    matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["mathtext.fontset"] = "custom"
matplotlib.rcParams["mathtext.it"] = "Arial:italic"
matplotlib.rcParams["mathtext.rm"] = "Arial"

OUT = "figs/Supp_S1_RunningFisher/output"
## drawn at the size it is pasted into the figure (180 x 125.6 mm, 600 dpi). The layout is given in axes fractions; line
## widths, arrow heads and markers are given for a 12.6 x 8.8 in drawing and scaled by SC; text is set in points: panel letters 8 pt bold, panel headings 7 pt bold, labels 7 pt, tick-like numbers 6 pt.
W_MM, H_MM, DPI = 180.0, 125.6, 600
MM = 1 / 25.4
SC = W_MM / (12.6 * 25.4)                 # scale of line widths, arrow heads and markers to the pasted size
FS_LETTER, FS_HEAD, FS_TXT, FS_SMALL = 8, 7, 7, 6
matplotlib.rcParams.update({"axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                            "xtick.major.pad": 2, "ytick.major.pad": 2})
RED, BLUE = "#C0392B", "#1F5FA8"          # same direction / opposite direction (paper-wide, Fig. 1)
HEAD, ZEBRA, INK, MUTED = "#222222", "#F2F2F2", "#222222", "#8A8A8A"
P_UNIVERSE = 26                           # gene universe of the toy example = the 26 genes GeneA-GeneZ

# ---------------------------------------------------------------- toy lists (rank 1 first)
## Each list ranks all 26 genes of the universe (GeneA-GeneZ): the 13 most increased form the UP list and the 13 most decreased the
## DOWN list (so the map of c is 13 x 13). b2 changes in the same direction as b1 for most genes.
LISTS = {
    "b1_up":   ["GeneA", "GeneB", "GeneJ", "GeneC", "GeneK", "GeneL", "GeneD", "GeneE", "GeneM", "GeneN", "GeneO", "GeneT", "GeneU"],
    "b1_down": ["GeneF", "GeneP", "GeneG", "GeneQ", "GeneR", "GeneH", "GeneS", "GeneI", "GeneV", "GeneW", "GeneX", "GeneY", "GeneZ"],
    "b2_up":   ["GeneB", "GeneA", "GeneC", "GeneM", "GeneN", "GeneO", "GeneD", "GeneK", "GeneJ", "GeneL", "GeneE", "GeneS", "GeneW"],
    "b2_down": ["GeneG", "GeneT", "GeneF", "GeneU", "GeneH", "GeneV", "GeneP", "GeneQ", "GeneR", "GeneI", "GeneX", "GeneY", "GeneZ"],
}
assert all(sorted(LISTS[f"{b}_up"] + LISTS[f"{b}_down"]) == [f"Gene{c}" for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"] for b in ("b1", "b2"))
EX_G, EX_K = 4, 5                         # the pair of cut-offs drawn in b


def fisher_p(l1, l2, G, K, U=P_UNIVERSE):
    """One-sided Fisher's exact test for the overlap of the top G of l1 and the top K of l2."""
    M = len(set(l1[:G]) & set(l2[:K]))
    return M, hypergeom.sf(M - 1, U, K, G)


def scan(l1, l2):
    """-log10 p for every (G, K); rows G = 1.., columns K = 1.."""
    return np.array([[-np.log10(fisher_p(l1, l2, g, k)[1]) for k in range(1, len(l2) + 1)]
                     for g in range(1, len(l1) + 1)])


PAIRS = [("b1_up", "b2_up", "up × up", RED), ("b1_down", "b2_down", "down × down", RED),
         ("b1_up", "b2_down", "up × down", BLUE), ("b1_down", "b2_up", "down × up", BLUE)]
SCORES = [scan(LISTS[a], LISTS[b]).max() for a, b, _, _ in PAIRS]
S_TOTAL = SCORES[0] + SCORES[1] - SCORES[2] - SCORES[3]
H_UPUP = scan(LISTS["b1_up"], LISTS["b2_up"])
BEST_G, BEST_K = [int(i) + 1 for i in np.unravel_index(np.argmax(H_UPUP), H_UPUP.shape)]
EX_M, EX_P = fisher_p(LISTS["b1_up"], LISTS["b2_up"], EX_G, EX_K)
print(f"example (G={EX_G}, K={EX_K}): M={EX_M}, p={EX_P:.2e}; best (G={BEST_G}, K={BEST_K}): "
      f"p={10 ** -H_UPUP.max():.2e}; scores={[round(s, 2) for s in SCORES]}; S_total={S_TOTAL:.2f}")


def sci(p):
    """p as 'a × 10^b' in mathtext."""
    e = int(np.floor(np.log10(p))); m = p / 10 ** e
    return rf"{m:.1f} \times 10^{{{e}}}"


# ---------------------------------------------------------------- drawing helpers (axes coordinates 0-1)
def letter(ax, s, title):
    # the panel letter (s) is not drawn: all panel letters are added in the final figure; its space is kept
    ax.text(0.07, 1.045, title, transform=ax.transAxes, fontsize=FS_HEAD, fontweight="bold", va="bottom", ha="left", color=INK)


def gene_column(ax, x, ytop, genes, header, w=0.22, rh=0.082, window=None, fs=FS_TXT):
    """Mini table in the style of Fig. 4c: header bar, zebra rows, rank and gene; optional shaded top-n window."""
    ax.add_patch(Rectangle((x, ytop), w, rh, facecolor=HEAD, edgecolor="none"))
    ax.text(x + w / 2, ytop + rh / 2, header, ha="center", va="center", color="white", fontsize=fs, fontweight="bold")
    centers = {}
    for i, g in enumerate(genes):
        y = ytop - (i + 1) * rh
        if window and i < window:
            ax.add_patch(Rectangle((x, y), w, rh, facecolor="#DADADA", edgecolor="none"))
        elif i % 2 == 0:
            ax.add_patch(Rectangle((x, y), w, rh, facecolor=ZEBRA, edgecolor="none"))
        ax.text(x + 0.18 * w, y + rh / 2, str(i + 1), ha="center", va="center", fontsize=fs, color=INK)
        ax.text(x + 0.62 * w, y + rh / 2, g, ha="center", va="center", fontsize=fs, color=INK)
        centers[g] = y + rh / 2
    return centers


def bracket(ax, x, ytop, ybot, label, side):
    d = 0.018 if side == "right" else -0.018
    ax.plot([x, x + d, x + d, x], [ytop, ytop, ybot, ybot], color=INK, lw=1.8 * SC, solid_capstyle="butt")
    ax.text(x + 2.6 * d, (ytop + ybot) / 2, label, ha="left" if side == "right" else "right", va="center",
            fontsize=FS_TXT, color=INK, fontstyle="italic")


fig = plt.figure(figsize=(round(W_MM * MM * DPI) / DPI, round(H_MM * MM * DPI) / DPI), dpi=DPI)   # whole pixels
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.28], height_ratios=[1.0, 0.92], wspace=0.16, hspace=0.34,
                      left=0.029, right=0.982, top=0.957, bottom=0.06)   # outer margins of ~1 mm


def vdots(ax, x, y, gap=0.028, s=7 * SC ** 2):
    """Three vertical dots (the vertical-ellipsis glyph is missing in Arial)."""
    ax.scatter([x] * 3, [y + gap, y, y - gap], s=s, color=INK, zorder=4)


# ---------------------------------------------------------------- a: the two gene lists and the four directional pairs
ax = fig.add_subplot(gs[0, 0]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
letter(ax, "a", "Two ranked gene lists and four pairs")
## Each list as a spreadsheet: genes ordered by the change (most increased at the top, most decreased at the bottom), hundreds
## of genes in between. The values are illustrative; the gene names match b.
TABLE = {"b1": [("GeneA", "+3.1"), ("GeneB", "+2.8"), ("GeneJ", "+2.5"), None, ("GeneG", "−2.4"), ("GeneP", "−2.8"), ("GeneF", "−3.3")],
         "b2": [("GeneB", "+3.0"), ("GeneA", "+2.7"), ("GeneC", "+2.3"), None, ("GeneF", "−2.5"), ("GeneT", "−2.9"), ("GeneG", "−3.4")]}
colx = {"b1": 0.16, "b2": 0.595}; tw, rh_a, ytab = 0.245, 0.085, 0.84
for lst, x in colx.items():
    ax.text(x + tw / 2, ytab + rh_a + 0.03, f"Gene list {lst}", ha="center", va="bottom", fontsize=FS_TXT, fontweight="bold")
    ax.add_patch(Rectangle((x, ytab), tw, rh_a, facecolor=HEAD, edgecolor="none"))
    ax.text(x + 0.33 * tw, ytab + rh_a / 2, "gene", ha="center", va="center", color="white", fontsize=FS_SMALL, fontweight="bold")
    ax.text(x + 0.77 * tw, ytab + rh_a / 2, "log2FC", ha="center", va="center", color="white", fontsize=FS_SMALL, fontweight="bold")
    for i, row in enumerate(TABLE[lst]):
        y = ytab - (i + 1) * rh_a
        if row is None:                                           # the hundreds of genes in between
            vdots(ax, x + 0.33 * tw, y + rh_a / 2, gap=0.022); vdots(ax, x + 0.77 * tw, y + rh_a / 2, gap=0.022)
            continue
        up = i < 3; strength = (3 - i) / 3 if up else (i - 3) / 3    # deeper colour towards the ends
        col = matplotlib.colors.to_rgba(RED if up else BLUE, 0.12 + 0.33 * strength)
        ax.add_patch(Rectangle((x, y), tw, rh_a, facecolor=col, edgecolor="white", lw=1.0 * SC))
        ax.text(x + 0.33 * tw, y + rh_a / 2, row[0], ha="center", va="center", fontsize=FS_TXT, color=INK)
        ax.text(x + 0.77 * tw, y + rh_a / 2, row[1], ha="center", va="center", fontsize=FS_TXT, color=INK)
    # the two sub-lists, each read from its end inwards (strongest change first)
    side = -1 if lst == "b1" else 1
    xa = x - 0.028 if lst == "b1" else x + tw + 0.028
    top_y0, top_y1 = ytab - 0.005, ytab - 3 * rh_a + 0.005
    bot_y0, bot_y1 = ytab - 7 * rh_a + 0.005, ytab - 4 * rh_a - 0.005
    ax.add_patch(FancyArrowPatch((xa, top_y0), (xa, top_y1), arrowstyle="-|>", mutation_scale=13 * SC, color=RED, lw=2.0 * SC))
    ax.add_patch(FancyArrowPatch((xa, bot_y0), (xa, bot_y1), arrowstyle="-|>", mutation_scale=13 * SC, color=BLUE, lw=2.0 * SC))
    ha = "right" if lst == "b1" else "left"
    ax.text(xa + side * 0.015, (top_y0 + top_y1) / 2, f"{lst}_up", ha=ha, va="center", fontsize=FS_TXT, color=RED, fontweight="bold")
    ax.text(xa + side * 0.015, (bot_y0 + bot_y1) / 2, f"{lst}_down", ha=ha, va="center", fontsize=FS_TXT, color=BLUE, fontweight="bold")
# the four pairs between the lists
yu, yd = ytab - 1.5 * rh_a, ytab - 5.5 * rh_a
xl, xr = colx["b1"] + tw + 0.01, colx["b2"] - 0.01
arr = dict(arrowstyle="<|-|>", mutation_scale=13 * SC, shrinkA=0, shrinkB=0)
ax.add_patch(FancyArrowPatch((xl, yu + 0.04), (xr, yu + 0.04), color=RED, lw=2.2 * SC, **arr))
ax.add_patch(FancyArrowPatch((xl, yd - 0.04), (xr, yd - 0.04), color=RED, lw=2.2 * SC, **arr))
ax.add_patch(FancyArrowPatch((xl, yu - 0.04), (xr, yd + 0.04), color=BLUE, lw=1.8 * SC, linestyle=(0, (4, 2.5)), **arr))
ax.add_patch(FancyArrowPatch((xl, yd + 0.04), (xr, yu - 0.04), color=BLUE, lw=1.8 * SC, linestyle=(0, (4, 2.5)), **arr))
ax.text(0.5, yu + 0.07, "up × up", ha="center", va="bottom", fontsize=FS_TXT, color=RED, fontweight="bold")
ax.text(0.5, yd - 0.07, "down ×\ndown", ha="center", va="top", fontsize=FS_TXT, color=RED, fontweight="bold", linespacing=1.05)
ax.text(0.5, (yu + yd) / 2, "up × down\ndown × up", ha="center", va="center", fontsize=FS_TXT, color=BLUE, fontweight="bold",
        linespacing=1.2, bbox=dict(facecolor="white", edgecolor="none", pad=1.2))

# ---------------------------------------------------------------- b: one pair, one pair of cut-offs (the 2 x 2 table as large as the lists)
ax = fig.add_subplot(gs[0, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
letter(ax, "b", "Overlap at one pair of cut-offs")
x1, x2, yt, rh, cw_b = 0.07, 0.33, 0.86, 0.082, 0.20
NB = 8                                                          # rows drawn; the lists go on to 13 (dots)
c1 = gene_column(ax, x1, yt, LISTS["b1_up"][:NB], "b1_up", w=cw_b, window=EX_G)
c2 = gene_column(ax, x2, yt, LISTS["b2_up"][:NB], "b2_up", w=cw_b, window=EX_K)
for xx in (x1, x2):
    vdots(ax, xx + 0.62 * cw_b, yt - (NB + 0.55) * rh, gap=0.02)
bracket(ax, x1 - 0.012, yt, yt - EX_G * rh, rf"$G$ = {EX_G}", "left")
bracket(ax, x2 + cw_b + 0.012, yt, yt - EX_K * rh, rf"$K$ = {EX_K}", "right")
for g in set(LISTS["b1_up"][:NB]) & set(LISTS["b2_up"][:NB]):
    inside = LISTS["b1_up"].index(g) < EX_G and LISTS["b2_up"].index(g) < EX_K
    ax.plot([x1 + cw_b + 0.002, x2 - 0.002], [c1[g], c2[g]], color=INK if inside else "#B0B0B0",
            lw=1.9 * SC if inside else 1.2, linestyle="-" if inside else (0, (3, 2)), zorder=3)
tx, ty, cwt, rht = 0.715, 0.62, 0.135, 0.15
ax.text(tx + cwt, 0.93, rf"$M$ = {EX_M}", ha="center", va="center", fontsize=FS_TXT, color=INK)
ax.text(tx + cwt, 0.80, "one-sided Fisher's\nexact test " + r"$\rightarrow$ " + rf"$P$ = {EX_P:.2g}", ha="center", va="center",
        fontsize=FS_TXT, color=INK, linespacing=1.25)   # the one P value of the figure (G = 4, K = 5, M = 3)
hdr = [rf"in top $G$", rf"not"]; rows = [rf"in top $K$", "not"]
cells = [[r"$M$", r"$K-M$"], [r"$G-M$", r"$U-G-K+M$"]]
VALS = [[EX_M, EX_K - EX_M], [EX_G - EX_M, P_UNIVERSE - EX_G - EX_K + EX_M]]   # 3, 2 / 1, 20 for G = 4, K = 5, M = 3, U = 26
for j, h in enumerate(hdr):
    ax.text(tx + (j + 0.5) * cwt, ty + 0.03, h, ha="center", va="bottom", fontsize=FS_TXT, color=INK)
for i, r in enumerate(rows):
    ax.text(tx - 0.012, ty - (i + 0.5) * rht, r, ha="right", va="center", fontsize=FS_TXT, color=INK)
    for j in range(2):
        ax.add_patch(Rectangle((tx + j * cwt, ty - (i + 1) * rht), cwt, rht,
                               facecolor="#DADADA" if (i, j) == (0, 0) else "white", edgecolor="#9A9A9A", lw=1.0 * SC))
        ax.text(tx + (j + 0.5) * cwt, ty - (i + 0.5) * rht + 0.022, cells[i][j], ha="center", va="center",
                fontsize=FS_TXT if (i, j) != (1, 1) else 5.5, color=INK)   # U-G-K+M: 5.5 pt to fit the cell
        # the numbers of this example (M = 3, ...), small, under the symbols
        ax.text(tx + (j + 0.5) * cwt, ty - (i + 0.5) * rht - 0.042, f"= {VALS[i][j]}", ha="center", va="center", fontsize=FS_SMALL, color=MUTED)
ax.text(tx + cwt, ty - 2 * rht - 0.035, rf"$U$ = all genes ({P_UNIVERSE} here)", ha="center", va="top", fontsize=FS_TXT, color=MUTED)   # "all genes" rather than the term "gene universe"

# ---------------------------------------------------------------- c: every pair of cut-offs; the smallest P represents the pair
ax = fig.add_subplot(gs[1, 0])
letter(ax, "c", r"Every pair of cut-offs ($G$, $K$) for up × up")
n = H_UPUP.shape[0]
im = ax.imshow(H_UPUP, cmap="Greys", vmin=0, vmax=np.ceil(H_UPUP.max()), origin="upper",
               extent=(0.5, n + 0.5, n + 0.5, 0.5), aspect="equal")
ax.set_xticks(range(1, n + 1)); ax.set_yticks(range(1, n + 1)); ax.tick_params(labelsize=FS_SMALL)
ax.set_xlabel(""); ax.set_ylabel("")
ax.set_xlim(0.5, n + 5.2); ax.set_ylim(n + 0.5, 0.5)
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
# every K and every G is tried: arrows along the two sides of the map
ax.add_patch(Rectangle((0.5, 0.5), n, n, fill=False, edgecolor="#6A6A6A", lw=1.0 * SC, zorder=4))   # outer frame (marks the edge where the map is pale)
# the words sit on the arrows, which run along the axes next to the tick labels
ax.annotate("", xy=(n + 0.3, n + 1.45), xytext=(0.7, n + 1.45), annotation_clip=False,
            arrowprops=dict(arrowstyle="-|>", mutation_scale=14 * SC, color=INK, lw=1.6 * SC))
ax.text((n + 1) / 2, n + 1.75, r"every $K$ (top genes of b2_up)", ha="center", va="top", fontsize=FS_TXT, color=INK, clip_on=False)
ax.annotate("", xy=(-0.75, n + 0.3), xytext=(-0.75, 0.7), annotation_clip=False,
            arrowprops=dict(arrowstyle="-|>", mutation_scale=14 * SC, color=INK, lw=1.6 * SC))
ax.text(-1.05, (n + 1) / 2, r"every $G$ (top genes of b1_up)", ha="right", va="center", fontsize=FS_TXT, color=INK, rotation=90, clip_on=False)
# the test of b and the smallest P
ax.add_patch(Rectangle((EX_K - 0.5, EX_G - 0.5), 1, 1, fill=False, edgecolor="white", lw=1.8 * SC, linestyle=(0, (2.5, 1.5))))
yb, ys = (EX_G + 0.6, BEST_G) if abs(EX_G - BEST_G) >= 2 else (EX_G - 1.4, BEST_G + 1.6)   # keep the labels apart (and below the colour key)
ax.annotate(r"the test of b" + "\n" + r"($G$ = 4, $K$ = 5)", xy=(EX_K + 0.5, EX_G), xytext=(n + 0.95, yb), fontsize=FS_TXT, ha="left", linespacing=1.25,
            va="center", color=INK, arrowprops=dict(arrowstyle="-", color=INK, lw=1.0 * SC, linestyle=(0, (2, 1.5))))
ax.add_patch(Rectangle((BEST_K - 0.5, BEST_G - 0.5), 1, 1, fill=False, edgecolor=RED, lw=3.0 * SC, zorder=5))
ax.annotate(r"smallest $P$" + "\n" + r"= the $P$ of up × up", xy=(BEST_K + 0.5, BEST_G + 0.35), xytext=(n + 0.95, ys), fontsize=FS_TXT,
            ha="left", va="center", color=RED, linespacing=1.3, arrowprops=dict(arrowstyle="-", color=RED, lw=1.4 * SC))
cax = ax.inset_axes([0.81, 0.955, 0.17, 0.035])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_ticks([]); cb.outline.set_visible(False)
cax.text(0.56, -0.9, r"$-\log_{10} P$ of the test" + "\n" + r"at each ($G$, $K$)", transform=cax.transAxes, ha="center", va="top", fontsize=FS_TXT, color=INK, linespacing=1.25)

# ---------------------------------------------------------------- d: the signed sum gives the direction (as in Fig. 1d; no numbers)
def draw_d(fig, ax, bar_box, lx0, lx1, lzero, ly, formula_y, tri_x, key_right=False):
    """Panel d in the axes ax (axes coordinates): the four scores (bar_box = inset position), the direction line from lx0 to
    lx1 with 0 at lzero, and the S_total formula centred at formula_y. key_right: 'S_total of this example' ends at lx1."""
    letter(ax, "d", "Signed sum gives the direction")
    bx = ax.inset_axes(bar_box)                                   # the four scores, small
    xs = np.arange(4)
    bx.bar(xs, SCORES, color=[c for *_, c in PAIRS], width=0.62)
    bx.set_xticks(xs); bx.set_xticklabels([lab.replace(" × ", "×\n") for _, _, lab, _ in PAIRS], fontsize=FS_SMALL)
    for t, (*_, c) in zip(bx.get_xticklabels(), PAIRS):
        t.set_color(c)
    bx.set_yticks([]); bx.set_ylabel(r"Score = $-\log_{10}$" + "\n" + r"(smallest $P$)", fontsize=FS_TXT)
    for sp in ("top", "right"):
        bx.spines[sp].set_visible(False)
    bx.spines["left"].set_color("#8A8A8A"); bx.spines["bottom"].set_color("#8A8A8A")
    ## "S_total" and " = " as two pieces (the width of a mathtext piece does not include a space after "=")
    parts = [(r"$S_{\mathrm{total}}$", INK), (" = ", INK), ("up × up", RED), (" + ", INK), ("down × down", RED),
             (" − ", INK), ("up × down", BLUE), (" − ", INK), ("down × up", BLUE)]
    rend = fig.canvas.get_renderer()
    pieces = [ax.text(0, formula_y, p_, transform=ax.transAxes, fontsize=FS_TXT, color=c, va="center", ha="left") for p_, c in parts]
    fig.canvas.draw()
    widths = [t.get_window_extent(renderer=rend).width for t in pieces]
    axw = ax.get_window_extent(renderer=rend).width
    x0 = 0.5 - sum(widths) / axw / 2                          # centred in the panel
    for t, w in zip(pieces, widths):
        t.set_x(x0); x0 += w / axw
    # the direction line (as Fig. 1d): opposite to the left, same to the right, 0 in between; S_total of this example is positive
    ax.add_patch(FancyArrowPatch((lzero, ly), (lx1, ly), arrowstyle="-|>", mutation_scale=15 * SC, color=RED, lw=2.0 * SC))
    ax.add_patch(FancyArrowPatch((lzero, ly), (lx0, ly), arrowstyle="-|>", mutation_scale=15 * SC, color=BLUE, lw=2.0 * SC))
    ax.plot([lzero, lzero], [ly - 0.03, ly + 0.03], color=INK, lw=1.4 * SC)
    ax.text(lzero, ly - 0.06, "0", ha="center", va="top", fontsize=FS_SMALL, color=INK)
    ax.text(lx1, ly - 0.06, "+  same\ndirection", ha="right", va="top", fontsize=FS_TXT, color=RED, linespacing=1.15)
    ax.text(lx0, ly - 0.06, "−  opposite\ndirection", ha="left", va="top", fontsize=FS_TXT, color=BLUE, linespacing=1.15)
    ax.scatter([tri_x], [ly + 0.075], marker="v", s=110 * SC ** 2, color=INK, zorder=5)
    if key_right:
        ax.text(lx1, ly + 0.11, r"$S_{\mathrm{total}}$ of this example", ha="right", va="bottom", fontsize=FS_TXT, color=INK)
    else:
        ax.text(tri_x, ly + 0.11, r"$S_{\mathrm{total}}$ of this example", ha="center", va="bottom", fontsize=FS_TXT, color=INK)


ax = fig.add_subplot(gs[1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
draw_d(fig, ax, [0.05, 0.40, 0.36, 0.48], 0.50, 0.97, 0.735, 0.62, 0.18, 0.90)

os.makedirs(OUT, exist_ok=True)
OUT_PPTX = OUT                                  # the pasted-size figure; no bbox_inches, the figure is the pasted size
os.makedirs(OUT_PPTX, exist_ok=True)
fig.savefig(os.path.join(OUT_PPTX, "FigS1_running_fisher.png"), dpi=DPI, facecolor="white")
print("wrote", os.path.join(OUT_PPTX, "FigS1_running_fisher.png"))

## ---- each panel as its own PNG, so that the panels can be arranged in the final figure ----
## Cut out of the whole figure above at the same scale (nothing is resized): a panel = what is drawn by the axes whose centre
## lies in its grid cell (texts and the heading included), plus 0.5 mm of white around. Output: FigS1{a,b,c,d}_*.png, 600 dpi.
from matplotlib.transforms import Bbox
import matplotlib.collections
fig.canvas.draw(); rend = fig.canvas.get_renderer()
PANELS = {"FigS1a_gene_lists": gs[0, 0], "FigS1b_one_cutoff": gs[0, 1], "FigS1c_every_cutoff": gs[1, 0]}   # d: drawn on its own below
PAD_IN = 0.5 * MM                                 # 0.5 mm in inches
fig_in = Bbox.from_bounds(0, 0, fig.get_figwidth(), fig.get_figheight())
for name, cell in PANELS.items():
    box = cell.get_position(fig)
    members = [a for a in fig.axes if box.x0 <= (a.get_position().x0 + a.get_position().x1) / 2 <= box.x1
               and box.y0 <= (a.get_position().y0 + a.get_position().y1) / 2 <= box.y1]
    ## the extent of what is drawn (texts, patches, lines, images, dots), not of the whole grid cell (the axes are switched off)
    ext = []
    for a in members:
        for art in a.get_children():
            if not art.get_visible() or art is a.patch or isinstance(art, (matplotlib.spines.Spine, matplotlib.axis.Axis)):
                continue
            try:
                if isinstance(art, matplotlib.collections.Collection) and len(art.get_offsets()):
                    ## dots (scatter): their positions plus the marker radius (get_window_extent misses the offsets)
                    pts = art.get_offset_transform().transform(art.get_offsets())
                    r = np.sqrt(np.max(art.get_sizes()) if len(art.get_sizes()) else 1.0) / 2 * fig.dpi / 72
                    e = Bbox.from_extents(pts[:, 0].min() - r, pts[:, 1].min() - r, pts[:, 0].max() + r, pts[:, 1].max() + r)
                else:
                    e = art.get_window_extent(rend)
            except Exception:
                continue
            if e.width > 0 and e.height > 0 and np.all(np.isfinite(e.get_points())):
                ext.append(e)
    bb = Bbox.union(ext).transformed(fig.dpi_scale_trans.inverted())
    bb = Bbox.from_extents(max(bb.x0 - PAD_IN, 0), max(bb.y0 - PAD_IN, 0), min(bb.x1 + PAD_IN, fig_in.x1), min(bb.y1 + PAD_IN, fig_in.y1))
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=DPI, facecolor="white", bbox_inches=bb)
    print(f"wrote {name}.png ({bb.width / MM:.1f} x {bb.height / MM:.1f} mm, {len(members)} axes)")

## d on its own, narrower than in the whole figure: 74 mm wide,
## the bars and the direction line closer together; the formula (about 62 mm at 7 pt) sets the width. Same text sizes.
D_W, D_H = 74.0, 40.0
figd = plt.figure(figsize=(round(D_W * MM * DPI) / DPI, round(D_H * MM * DPI) / DPI), dpi=DPI)
axd = figd.add_axes([0.0, 0.0, 1.0, 0.885]); axd.set_xlim(0, 1); axd.set_ylim(0, 1); axd.axis("off")
draw_d(figd, axd, [0.18, 0.40, 0.34, 0.55], 0.575, 0.985, 0.775, 0.62, 0.09, 0.91, key_right=True)
figd.savefig(os.path.join(OUT, "FigS1d_signed_sum.png"), dpi=DPI, facecolor="white")
print(f"wrote FigS1d_signed_sum.png ({D_W:.1f} x {D_H:.1f} mm, drawn on its own)")
