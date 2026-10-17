# -*- coding: utf-8 -*-
"""Fig4a (schematic): what in vivo Perturb-seq makes possible.
  1 a knock-out library in the brain of Cas9 mice (one guide per neuron); 2 single-nucleus RNA-seq
  reads the knocked-out gene, the cell type and the expression change of each nucleus; 3 a reference panel: every cell type x every
  knocked-out gene (m x n; no concrete numbers, they depend on the experiment); 4 the RFind score of one target (a ranked list, e.g.
  CA1 Glut, Stress vs No Stress) against the reference panels (one card = one knocked-out gene in one cell type).
Schematic only: no data values. Design facts of the data (AAV library, Cas9 mice, sorted neuronal
nuclei, probe-based single-nucleus RNA-seq) are in the legend Notes.
  <python> figs/Fig4_Perturbseq/analysis/fig4a_schematic.py      (from the project root)
Output: figs/Fig4_Perturbseq/output/Fig4a_schematic.png
Drawn at the size the panel has in the final figure (W_MM x H_MM below, 600 dpi, Arial). The drawing is laid out in mm of a
  180 x 45 layout (data unit = 1 mm) and scaled to the final size (SC); the shapes are sized in units of U and the line widths
  and marker sizes in units of LW. The header boxes span the whole width, each at least as wide as its text, 'each nucleus' goes
  on two lines under the nucleus, and the columns of 3 are spaced so that their labels stay apart. Bold only in the headers.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Ellipse, Rectangle
from matplotlib.transforms import Bbox

OUT = os.path.join('figs', 'Fig4_Perturbseq', 'output'); os.makedirs(OUT, exist_ok=True)
## The layout below (mm of the 180 x 45 drawing = drawing units) is drawn at the final size, i.e. at SC = 168.7 / 180 = 0.937,
## lines and points included; text 7 pt, the four headers 7 pt bold (at 7.5 pt the four do not fit in one row of 168.7 mm).
W_MM, H_MM, DPI = 168.7, 42.0, 600  # size in the final figure
DW = 180.0; SC = W_MM / DW; DH = H_MM / SC   # drawing units: mm of the 180 x 45 layout; same scale on x and y
FS, FSH = 7, 7                      # pt: every text (Nature: text <= 7 pt); the headers (bold)
LW = 0.5 * SC                       # unit of the line widths and marker sizes (pt), scaled with the drawing (SC)
U = 4.05                             # size unit of the shapes (mm of the 180 x 45 layout)
NAVY, GREY, DARK = '#1F3A5F', '#9AA3AD', '#2B2B2B'
KO = {'A': '#D1495B', 'B': '#EDAE49', 'C': '#7B5EA7', 'D': '#00798C', 'ctrl': '#9AA3AD'}
RED, BLUE, TEAL = '#C0392B', '#1F5FA8', '#028090'
matplotlib.rcParams['font.family'] = 'Arial'

## whole pixels: round(mm / 25.4 * 600) (+0.1 px so that matplotlib does not floor to one pixel less); data unit = 1 mm
fig = plt.figure(figsize=[(round(v / 25.4 * DPI) + 0.1) / DPI for v in (W_MM, H_MM)])
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, DW); ax.set_ylim(0, DH); ax.axis('off')

## layout (mm). The header boxes sit at fixed positions (XS below); _tw = the width of each 7-pt bold header text, used to check
## that each box holds its text with PADH on each side. YA = height of the three arrows between the boxes.
HEADERS = ['1. Knock-out library in the brain', '2. Single-nucleus RNA-seq', '3. A comprehensive reference panel',
           '4. RFind score against a target']
PADH = 0.9 / SC
_r = fig.canvas.get_renderer()
_tw = [fig.text(0, 0, h, fontsize=FSH, fontweight='bold').get_window_extent(_r).width / fig.dpi * 25.4 / SC for h in HEADERS]
for t in list(fig.texts): t.remove()
## the row spans the whole width and the spare width goes to panel 2 and 3, not to the margins; each box still at least as
## wide as its text
XS = [(0.4, 42.3), (43.4, 83.2), (84.3, 137.6), (138.7, 179.6)]
for (a, b), w in zip(XS, _tw):
    assert b - a >= w + 2 * PADH - 0.05, (a, b, w + 2 * PADH)
X1, X2, X3, X4 = XS
print('header boxes (drawing units):', [(round(a, 2), round(b, 2)) for a, b in XS])
YH, HH = 42.75, 3.7   # header row: centre, height
YA = 19.5

def header(x0, x1, text):
    ax.add_patch(FancyBboxPatch((x0, YH - HH / 2), x1 - x0, HH, boxstyle=f'round,pad=0,rounding_size={HH / 2}', fc=NAVY, ec='none'))
    ax.text((x0 + x1) / 2, YH, text, ha='center', va='center', color='white', fontsize=FSH, fontweight='bold')

def guide(x, y, col, s=1.0):   # sgRNA hairpin: stem + loop (s = size in units of U)
    S = s * U
    ax.plot([x, x], [y, y + 0.9 * S], color=col, lw=2.6 * s * LW, solid_capstyle='round')
    ax.add_patch(Circle((x, y + 1.08 * S), 0.2 * S, fc='none', ec=col, lw=2.2 * s * LW))
    ax.plot([x - 0.25 * S, x + 0.25 * S], [y, y], color=col, lw=2.6 * s * LW, solid_capstyle='round')

def arrow(x0, x1, y=YA):   # no shrink at the ends (2 pt each by default): at this size it would leave only the head
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle='-|>', mutation_scale=20 * LW, lw=2.2 * LW, color=DARK, shrinkA=0, shrinkB=0))

def ladder(xc, ytop, w, h, n, ec=DARK, lw=0.8):   # a ranked gene list: n cells from red (up, top) to blue (down, bottom)
    cm = LinearSegmentedColormap.from_list('rank', [RED, '#FFFFFF', BLUE])
    for i in range(n):
        ax.add_patch(Rectangle((xc - w / 2, ytop - (i + 1) * h / n), w, h / n, fc=cm(i / (n - 1)), ec=ec, lw=lw * LW))

def bars(xc, yc, w, h, n, base=True):   # one expression change drawn as sorted bars (a waterfall): genes from the most
    # up (red, above the line) on the left to the most down (blue, below the line) on the right
    t = np.linspace(1, -1, n); v = np.sign(t) * np.abs(t) ** 0.7
    bw = w / n
    for i, a in enumerate(v):
        ax.add_patch(Rectangle((xc - w / 2 + i * bw + 0.1 * bw, yc), 0.8 * bw, a * h / 2, fc=RED if a > 0 else BLUE, ec='none'))
    if base: ax.plot([xc - w / 2, xc + w / 2], [yc, yc], color=DARK, lw=0.6 * LW)

rng = np.random.default_rng(3)
keys = list(KO)

# ---- 1 pooled knock-out library in the brain ----
header(*X1, HEADERS[0])
for i, (k, c) in enumerate(KO.items()):
    guide(5.8 + i * 0.74 * U, 32.6, c, 0.92)   # centred over the text below
ax.text(11.8, 31.4, 'AAV library\n(4 guides per gene)', ha='center', va='top', fontsize=FS, color=DARK)   # 4 guides per gene (Methods)
# brain with neurons, each carrying one guide
BX, BY, UB = 28.6, 18.2, 1.3 * U   # centre of the brain; the brain drawn in units of UB = 1.3 U, to fill panel 1
ax.add_patch(FancyArrowPatch((21.0, 37.6), (27.6, 26.6), arrowstyle='-|>', mutation_scale=20 * LW, lw=2.2 * LW, color=DARK,
                             connectionstyle='arc3,rad=-0.35'))   # guides into the brain
## a round, soft sagittal mouse brain: soft brown outline, pale fill (no highlight), a round cerebellum with folia, a small
## olfactory bulb; neurons as white discs with a coloured dot (the colour of each neuron's guide)
BRAIN, BRAIN2, BEDGE = '#FBEFE3', '#F5E1CF', '#8A6A55'
from matplotlib.path import Path
from scipy.interpolate import splprep, splev
def smooth(P, n=300):
    P = np.asarray(P, float)
    tck, _ = splprep([np.r_[P[:, 0], P[0, 0]], np.r_[P[:, 1], P[0, 1]]], s=0, per=True)
    xx, yy = splev(np.linspace(0, 1, n), tck)
    return np.c_[BX + xx * UB, BY + yy * UB]
STEM = smooth([(0.7, -0.4), (1.15, -0.55), (1.38, -0.9), (1.22, -1.12), (0.95, -0.95), (0.72, -0.7)])
ax.add_patch(matplotlib.patches.Polygon(STEM, closed=True, fc=BRAIN2, ec=BEDGE, lw=1.5 * LW, joinstyle='round'))
CBX, CBY = 1.42, -0.1
ax.add_patch(Circle((BX + CBX * UB, BY + CBY * UB), 0.56 * UB, fc=BRAIN2, ec=BEDGE, lw=1.5 * LW))
for r in (0.4, 0.27):
    th = np.linspace(np.radians(-75), np.radians(75), 40)
    ax.plot(BX + (CBX + 0.04 + r * np.cos(th)) * UB, BY + (CBY + r * np.sin(th)) * UB, color=BEDGE, lw=0.8 * LW, solid_capstyle='round')
OUTL = smooth([(-2.15, -0.15), (-2.05, 0.18), (-1.75, 0.3), (-1.45, 0.6), (-0.9, 1.08), (-0.1, 1.3), (0.65, 1.22),
               (1.15, 0.9), (1.32, 0.42), (1.18, -0.05), (0.8, -0.45), (0.05, -0.72), (-0.8, -0.72), (-1.45, -0.52),
               (-1.8, -0.48), (-2.08, -0.4)])
ax.add_patch(matplotlib.patches.Polygon(OUTL, closed=True, fc=BRAIN, ec=BEDGE, lw=1.9 * LW, joinstyle='round'))
yy = np.linspace(0.42, -0.44, 20)
ax.plot(BX + (-1.68 + 0.06 * np.sin((yy + 0.44) / 0.86 * np.pi)) * UB, BY + yy * UB, color=BEDGE, lw=1.0 * LW, solid_capstyle='round')   # olfactory bulb
path = Path(OUTL)
pts, tries = [], 0
while len(pts) < 14 and tries < 8000:
    tries += 1
    q = (rng.uniform(-1.4, 1.0), rng.uniform(-0.5, 0.95))
    xy = (BX + q[0] * UB, BY + q[1] * UB)
    ring = [(xy[0] + 0.42 * UB * np.cos(a), xy[1] + 0.42 * UB * np.sin(a)) for a in np.linspace(0, 2 * np.pi, 12, endpoint=False)]
    if not (path.contains_point(xy) and all(path.contains_point(r) for r in ring)) or q[0] < -1.45: continue
    if q[0] < -0.4 and q[1] > 0.55: continue   # no neurons in the upper left of the brain
    if any((xy[0] - a) ** 2 + (xy[1] - b) ** 2 < (0.48 * UB) ** 2 for a, b in pts): continue
    pts.append(xy)
for i, (x, y) in enumerate(pts):
    ax.add_patch(Circle((x, y), 0.19 * UB, fc='white', ec=DARK, lw=0.8 * LW))   # neurons: white disc with a coloured dot
    ax.add_patch(Circle((x, y), 0.09 * UB, fc=list(KO.values())[i % len(KO)], ec='none'))
ax.text(BX - 0.55 * UB, BY - 1.25 * UB - 1.0, 'Cas9 mouse brain:\none knock-out per neuron', ha='center', va='top', fontsize=FS, color=DARK, linespacing=1.15)
## (the source, Shi et al., 2026, is cited in the legend; no DOI in the figure)

arrow(39.6, 43.6)
# ---- 2 single-nucleus RNA-seq ----
header(*X2, HEADERS[1])
CX, R = 48.6, 1.0 * U   # the nucleus (inside panel 2)
ax.add_patch(Circle((CX, YA), R, fc='white', ec=DARK, lw=1.6 * LW))
guide(CX - 0.4 * U, YA - 0.15 * U, KO['B'], 0.55)
for dy in (-0.38, 0.0, 0.38):
    xs = np.linspace(0.02, 0.67, 30); ax.plot(CX + xs * U, YA + (dy + 0.07 * np.sin(xs * 20)) * U, color=DARK, lw=1.3 * LW)
ax.text(CX, YA - R - 0.7, 'each\nnucleus', ha='center', va='top', fontsize=FS, color=DARK, linespacing=1.1)
## three read-outs in one column, joined to the nucleus by lines; icons centred on IC, labels from LX
IC, LX = 61.0, 65.6   # close together so that the labels end inside panel 2
ROWS = [(YA + 13.0, 'knocked-out\ngene', 'guide'), (YA + 1.5, 'cell\ntype', 'type'), (YA - 9.5, 'expression\nchange', 'ladder')]
for y, txt, kind in ROWS:
    ax.plot([CX + R * 0.98, IC - {'type': 3.4, 'guide': 1.8}.get(kind, 3.6)], [YA + (y - YA - 1.5) * 0.12, y], color=GREY, lw=1.0 * LW, solid_capstyle='round')
    if kind == 'guide':
        guide(IC, y - 0.6 * U, KO['B'], 0.7)
    elif kind == 'type':
        ax.add_patch(FancyBboxPatch((IC - 3.0, y - 2.3), 6.0, 4.6, boxstyle='round,pad=0.08', fc='#E8EEF5', ec=NAVY, lw=1.1 * LW))
        ax.text(IC, y, 'CA1\nGlut', ha='center', va='center', fontsize=6, color=NAVY, linespacing=1.0)   # two lines, to fit panel 2
    else:
        bars(IC, y, 6.0, 5.0, 9)
    ax.text(LX, y, txt, ha='left', va='center', fontsize=FS, color=DARK, linespacing=1.1)

arrow(81.0, 84.6)   # between rows KO 3 and KO n of 3
# ---- 3 a comprehensive reference panel: every cell type x every knocked-out gene ----
# no concrete numbers: the 1,947 knock-outs are what this experiment happened to give and can change with the conditions,
# and the reference panel has as many cell types as Perturb-seq identifies, not only these 4; 4a shows what is possible
N_CT, N_KO = 'm', 'n'
header(*X3, HEADERS[2])
CW, CH = 1.05 * U, 0.75 * U   # card size
def list_card(xc, yc, w=CW, h=CH, ec='#8A94A0', lw=1.0):   # a card = the expression change of one knock-out (sorted bars)
    ax.add_patch(FancyBboxPatch((xc - w / 2, yc - h / 2), w, h, boxstyle=f'round,pad={0.02 * U},rounding_size={0.07 * U}', fc='white', ec=ec, lw=lw * LW))
    bars(xc, yc, 0.78 * w, 0.78 * h, 7)
def vdots(x, y, d=0.22 * U):
    ax.plot([x] * 3, [y + d, y, y - d], 'o', color=DARK, markersize=2.4 * LW)
def hdots(x, y, d=0.22 * U):
    ax.plot([x - d, x, x + d], [y] * 3, 'o', color=DARK, markersize=2.4 * LW)
def grid(cx, ry, clab, rlab, cell):
    xl = cx[0] - CW / 2 - 0.9   # right end of the row labels
    for c, lab in zip(cx, clab):
        if lab: ax.text(c, ry[0] + 5.0, lab, ha='center', va='center', fontsize=FS, color=NAVY, linespacing=1.0)
        else: hdots(c, ry[0] + 5.0)
    for r, (y, lab) in enumerate(zip(ry, rlab)):
        if lab: ax.text(xl, y, lab, ha='right', va='center', fontsize=FS, color=DARK)
        else: vdots(xl - 2.85, y)
        for k, (c, cl) in enumerate(zip(cx, clab)):
            if lab and cl: cell(c, y, r, k)
            elif lab: hdots(c, y)
            elif cl: vdots(c, y)
rlab = ['KO 1', 'KO 2', 'KO 3', None, f'KO {N_KO}']
clab = ['Cell\ntype 1', 'Cell\ntype 2', 'Cell\ntype 3', None, f'Cell\ntype {N_CT}']
ry = [29.5, 25.5, 21.5, 17.5, 13.5]
cx3 = [93.8, 103.6, 113.4, 119.8, 126.8]   # labelled columns 9.8 apart, so that the labels stay apart
grid(cx3, ry, clab, rlab, lambda x, y, r, k: list_card(x, y))
xr = cx3[-1] + CW / 2 + 0.6   # braces: m cell types (right), n knocked-out genes (below)
ax.plot([xr - 0.8, xr, xr, xr - 0.8], [ry[0] + 1.6, ry[0] + 1.6, ry[-1] - 1.6, ry[-1] - 1.6], color=DARK, lw=1.3 * LW)
ax.text(xr + 0.7, 26.0, f'{N_KO} KOs', ha='left', va='center', fontsize=FS, color=DARK)   # one line, above the arrow
yb = ry[-1] - CH / 2 - 1.0
ax.plot([cx3[0] - CW / 2, cx3[0] - CW / 2, cx3[-1] + CW / 2, cx3[-1] + CW / 2], [yb + 0.7, yb, yb, yb + 0.7], color=DARK, lw=1.3 * LW)
ax.text((cx3[0] + cx3[-1]) / 2, yb - 0.8, f'{N_CT} cell types', ha='center', va='top', fontsize=FS, color=DARK)
## key of the cards
xk, yk = cx3[0] - CW / 2 + 0.6, 4.6
bars(xk + 0.6, yk, 2.4, 2.4, 7)
ax.text(xk + 2.4, yk, 'each card = genes sorted, up (red) to down (blue)', ha='left', va='center', fontsize=6, color=DARK)

arrow(141.0, 145.0, y=YA + 4.0)   # after 'n KOs'
# ---- 4 one target matched against every reference -> an RFind score for each ----
# (with a comprehensive panel, one target, e.g. dentate gyrus Stress vs No Stress, can be matched against every
#  reference, a CA3 panel as well as the dentate-gyrus one; each match gives an RFind score)
header(*X4, HEADERS[3])
# the reference panels (left, next to panel 3; one card = one knocked-out gene in one cell type) vs. the target (right; a ranked
# list, red = up, blue = down)
## the reference and the target are both ranked gene lists (the same ladder), matched by RFind-sc
nT, hT, wT, yT0 = 10, 0.5 * U, 1.1 * U, 35.4   # sized and placed to fill panel 4
for xc, top, bottom in [(153.0, 'Reference', 'KO i in\ncell type j'), (171.0, 'Target', 'e.g. CA1 Glut,\nStress vs\nNo Stress')]:
    ladder(xc, yT0, wT, nT * hT, nT)
    ax.text(xc, yT0 + 1.0, top, ha='center', va='bottom', fontsize=FS, color=DARK)
    ax.text(xc, yT0 - nT * hT - 1.0, bottom, ha='center', va='top', fontsize=FS, color=DARK, linespacing=1.15)
ax.text(162.0, yT0 - nT * hT / 2, 'vs.', ha='center', va='center', fontsize=FS, color=DARK)
ax.text(153.0 - wT / 2 - 0.8, yT0 - hT / 2, 'up', ha='right', va='center', fontsize=6, color=RED)
ax.text(153.0 - wT / 2 - 0.8, yT0 - nT * hT + hT / 2, 'down', ha='right', va='center', fontsize=6, color=BLUE)

## check: everything inside the figure (mm from each edge)
fig.canvas.draw(); rd = fig.canvas.get_renderer()
bb = Bbox.union([a.get_tightbbox(rd) for a in ax.get_children() if a is not ax.patch and a.get_visible()
                 and not isinstance(a, (matplotlib.spines.Spine, matplotlib.axis.Axis))
                 and not (isinstance(a, matplotlib.text.Text) and not a.get_text())])
(x0, y0), (x1, y1) = ax.transData.inverted().transform(bb)
print(f'drawn: x {x0:.2f}-{x1:.2f}, y {y0:.2f}-{y1:.2f} drawing units (figure {DW} x {DH:.2f}; = {W_MM} x {H_MM} mm)')
## check: texts that overlap each other (drawing units)
_tb = [(t.get_text().replace(chr(10), ' ')[:22], ax.transData.inverted().transform(t.get_window_extent(rd))) for t in ax.texts if t.get_text()]
for i in range(len(_tb)):
    for j in range(i + 1, len(_tb)):
        (a0, b0), (a1, b1) = _tb[i][1]; (c0, d0), (c1, d1) = _tb[j][1]
        if min(a1, c1) - max(a0, c0) > 0 and min(b1, d1) - max(b0, d0) > 0:
            print(f'  TEXT OVERLAP: {_tb[i][0]!r} / {_tb[j][0]!r}')
fig.savefig(os.path.join(OUT, 'Fig4a_schematic.png'), dpi=DPI, facecolor='white')   # no bbox_inches='tight': the size is the final size
print('wrote Fig4a_schematic')
