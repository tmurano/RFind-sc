# -*- coding: utf-8 -*-
"""Fig4b (schematic): the whole panel b: the target (CA1 Glut, Stress vs No Stress; Fig. 3) is matched by RFind-sc
against every knock-out of the same cell type (Fig. 4a), and the knock-outs come out ranked (Knock-out 1 = most opposed, as in
Fig4c, Fig4d and Supplementary Table 1).
The cloud on the left is the UMAP of the CA1 glutamatergic neurons of the control line (fixed coordinates of
data/DeDonno2026/output/cell_meta_qc.csv), coloured by condition; everything else is schematic (no data values or gene names).
  <python> figs/Fig4_Perturbseq/analysis/fig4b_schematic.py      (from the project root)
Output: figs/Fig4_Perturbseq/output/Fig4b_schematic.png
Drawn at the size the panel has in the final figure (600 dpi, Arial): the layout is drawn at 80.0 x 40.0 mm with 7-pt text and
  scaled to fill that area (fit_to_mm), then scaled as a whole to the final width (OUT_W_MM below).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import Bbox

OUT = os.path.join('figs', 'Fig4_Perturbseq', 'output'); os.makedirs(OUT, exist_ok=True)
W_MM, H_MM, DPI, PAD_MM = 80.0, 40.0, 600, 0.4   # size the layout is drawn at (7-pt text); outer margin
## b is 74.1 x 37.1 mm in the final figure (92.7%). Drawn at that size with 7-pt text, the two titles and 'Knock-out n-2' and
## 'aligned' would collide, so the 80 x 40 drawing is scaled as a whole to the final size: text 6.5 pt, lines and points 92.7%.
## Written at 600 dpi (pixels = OUT_W_MM at 600 dpi).
OUT_W_MM = 74.1                                    # size in the final figure (width; the height follows, 37.0 mm)
FS = 7                                             # pt, every text
RED, BLUE, GREY = '#C0392B', '#1F5FA8', '#333333'          # Fig1 red / blue (red = up, aligned; blue = down, opposed)
KO = ['#D1495B', '#EDAE49', '#7B5EA7', '#2E9E44', '#9AA3AD']  # guide colours as in Fig4a
matplotlib.rcParams['font.family'] = 'Arial'

## whole pixels: round(mm / 25.4 * 600) (+0.1 px so that matplotlib does not floor to one pixel less)
fig = plt.figure(figsize=[(round(v / 25.4 * DPI) + 0.1) / DPI for v in (W_MM, H_MM)])
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 19.2); ax.set_ylim(0, 8.8); ax.axis('off')

def fit_to_mm(fig, ax, w_mm, h_mm, pad_mm):
    """Set the data limits so that everything drawn fills the figure (w x h mm, less pad_mm on each side) at one scale on x
    and y (aspect equal). Text keeps its size in pt, so its extent in data units depends on the scale: repeated until stable."""
    for a in ax.get_children():
        a.set_clip_on(False)
    for _ in range(6):
        fig.canvas.draw(); r = fig.canvas.get_renderer()
        bbs = [a.get_tightbbox(r) for a in ax.get_children()
               if a is not ax.patch and a.get_visible() and not isinstance(a, (matplotlib.spines.Spine, matplotlib.axis.Axis))
               and not (isinstance(a, matplotlib.text.Text) and not a.get_text())]
        (x0, y0), (x1, y1) = ax.transData.inverted().transform(Bbox.union([b for b in bbs if b is not None and b.width > 0]))
        k = min((w_mm - 2 * pad_mm) / (x1 - x0), (h_mm - 2 * pad_mm) / (y1 - y0))   # mm per data unit
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        ax.set_xlim(cx - w_mm / 2 / k, cx + w_mm / 2 / k); ax.set_ylim(cy - h_mm / 2 / k, cy + h_mm / 2 / k)
    return k

# ---- target: CA1 Glut of Fig. 3, Stress vs No Stress (the example of Fig. 4c) ----
## the title on three lines, centred over the UMAP and level with the title of the knock-outs
ax.text(2.7, 7.95, 'Target: CA1 Glut,\nStress vs. No Stress\n(Fig. 3)', ha='center', va='center', fontsize=FS, color=GREY, linespacing=1.2)
m = pd.read_csv(os.path.join('data', 'DeDonno2026', 'output', 'cell_meta_qc.csv'), usecols=['cell_type', 'condition', 'line', 'umap_x', 'umap_y'])
m = m[(m.line == 'Ctrl') & (m.cell_type == 'CA1 Glut')].sample(frac=1, random_state=1)   # shuffled so that neither condition hides the other
au = ax.inset_axes([0.6, 2.35, 4.2, 4.5], transform=ax.transData); au.axis('off')
## the view is the main cloud of the CA1 Glut cells: a small separate group of cells (about 2%) far from it is left out, as in
## the final figure (largest group of cells linked within 0.6 UMAP units; single linkage)
from scipy.cluster.hierarchy import linkage, fcluster
lab = fcluster(linkage(m[['umap_x', 'umap_y']].to_numpy(), 'single'), t=0.6, criterion='distance')
main = lab == np.bincount(lab).argmax(); print(f'main cloud: {main.sum()} of {len(m)} cells')
mm = m[main]
au.scatter(mm.umap_x, mm.umap_y, s=1.6, c=np.where(mm.condition == 'Stress', RED, BLUE), alpha=0.75, linewidths=0, rasterized=True)   # point size and alpha set so that Stress and No Stress are easy to tell apart
au.set_aspect('equal', adjustable='datalim')
for y, lab, col in [(1.75, 'No Stress', BLUE), (1.1, 'Stress', RED)]:
    ax.plot([1.35], [y], 'o', color=col, markersize=3.5)
    ax.text(1.75, y, lab, ha='left', va='center', fontsize=FS, color=GREY)

ax.text(5.72, 4.45, 'vs.', ha='center', va='center', fontsize=FS, color=GREY)

# ---- every knock-out of the same cell type (Fig. 4a) ----
ax.text(8.0, 7.95, 'Every knock-out\nin the same cell type\n(Fig. 4a)', ha='center', va='center', fontsize=FS, color=GREY, linespacing=1.2)
ax.add_patch(FancyBboxPatch((6.65, 0.95), 2.7, 5.55, boxstyle='round,pad=0.02,rounding_size=0.2', fc='#F5F6F8', ec='none'))
def guide(x, y, col, s=0.55):   # sgRNA hairpin, as in Fig4a
    ax.plot([x, x], [y, y + 0.9 * s], color=col, lw=1.2, solid_capstyle='round')
    ax.add_patch(Circle((x, y + 1.08 * s), 0.2 * s, fc='none', ec=col, lw=1.0))
    ax.plot([x - 0.25 * s, x + 0.25 * s], [y, y], color=col, lw=1.2, solid_capstyle='round')
rng = np.random.default_rng(4)
for i, (y, col) in enumerate(zip([5.75, 4.85, 3.95, 3.05, 2.15], KO)):
    guide(7.35, y - 0.35, col)
    ax.add_patch(Circle((8.65, y), 0.36, fc='white', ec=col, lw=0.8, ls=(0, (2, 1.5))))   # a nucleus carrying that knock-out
    ax.add_patch(Circle((8.65, y), 0.13, fc=col, ec='none'))
    for k in range(6):
        a = rng.uniform(0, 2 * np.pi); r = rng.uniform(0.18, 0.27)
        ax.add_patch(Circle((8.65 + r * np.cos(a), y + r * np.sin(a)), 0.03, fc=col, ec='none', alpha=0.7))
ax.plot([8.0] * 3, [1.55, 1.3, 1.05], 'o', color=GREY, markersize=1.5)   # more knock-outs

# ---- RFind-sc -> knock-outs ranked from opposed to aligned (Knock-out 1 = most opposed) ----
ax.add_patch(FancyArrowPatch((9.75, 4.45), (11.45, 4.45), arrowstyle='-|>', mutation_scale=11, lw=1.1, color=GREY))
ax.text(10.6, 4.8, 'RFind-sc', ha='center', va='bottom', fontsize=FS, color=GREY)
cmap = LinearSegmentedColormap.from_list('rank', [RED, '#F2F2F2', BLUE])   # origin='lower': top = opposed (blue)
ax.imshow(np.linspace(0, 1, 256).reshape(-1, 1), extent=(11.85, 12.25, 1.0, 7.1), origin='lower', cmap=cmap, aspect='auto', zorder=2)
ax.text(11.85, 7.55, 'Ranked knock-outs', ha='left', va='center', fontsize=FS, color=GREY)
rows = [(6.6, 'Knock-out 1', BLUE), (6.0, 'Knock-out 2', BLUE), (5.4, 'Knock-out 3', BLUE),
        (2.7, 'Knock-out n−2', RED), (2.1, 'Knock-out n−1', RED), (1.5, 'Knock-out n', RED)]
for y, lab, col in rows:
    ax.text(12.6, y, lab, ha='left', va='center', fontsize=FS, color=col)
ax.plot([13.45] * 3, [4.55, 4.05, 3.55], 'o', color=GREY, markersize=1.5)
ax.text(16.45, 6.0, 'opposed', ha='center', va='center', fontsize=FS, color=BLUE, rotation=90)   # clear of 'Knock-out n−2' at 7 pt
ax.text(16.45, 2.1, 'aligned', ha='center', va='center', fontsize=FS, color=RED, rotation=90)

k = fit_to_mm(fig, ax, W_MM, H_MM, PAD_MM); print(f'scale: {k:.2f} mm per unit')
png = os.path.join(OUT, 'Fig4b_schematic.png')
fig.savefig(png, dpi=DPI * OUT_W_MM / W_MM, facecolor='white')   # the whole drawing at 92.7%; no bbox_inches='tight'
from PIL import Image
im = Image.open(png); im.load(); im.save(png, dpi=(DPI, DPI))   # mark it 600 dpi -> pasted at 100% it is OUT_W_MM wide
print(f'{im.size[0]} x {im.size[1]} px = {im.size[0] / DPI * 25.4:.1f} x {im.size[1] / DPI * 25.4:.1f} mm at {DPI} dpi')
print('wrote Fig4b_schematic')
