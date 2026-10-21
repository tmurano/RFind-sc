# -*- coding: utf-8 -*-
"""Fig4c — knock-out ranking tables: the 20 KOs most aligned with, and the 20 most opposed to, the
stress response of one cell type (De Donno, Stress vs No Stress; CA1 Glut in Fig. 4c) among the Perturb-seq KO effect
vectors of the same cell type.

  <python> figs/Fig4_Perturbseq/analysis/fig4c_rank_tables.py [axis] [context]     (default: CA1 CA1 = the example of Fig. 4c,
  the neuronal type with the highest GR expression [Fig. 4e])

Input : analysis/perturbseq/results/axis_<axis>_vs_<context>.csv   (analysis/perturbseq/03_axis_vs_ko.py)
Output: figs/Fig4_Perturbseq/output/Fig4c_<axis>_<context>_{aligned,opposed}.png
  Drawn at the size the panels have in the final figure (FIG_W_MM x FIG_H_MM below, 600 dpi, Arial): the title bar and the
  column heads 5.0 pt bold (the title has to stay on one line), the rows 5.5 pt.
Layout: two separate tables, 20 rows each, columns #, Knock-out, RFind score. The purpose is the list itself;
no significance marks.
"""
import os, sys, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import font_manager

if any(f.name == 'Arial' for f in font_manager.fontManager.ttflist):
    matplotlib.rcParams['font.family'] = 'Arial'
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
AX = sys.argv[1] if len(sys.argv) > 1 else 'CA1'
CTX = sys.argv[2] if len(sys.argv) > 2 else 'CA1'
OUT = os.path.join(BASE, 'figs', 'Fig4_Perturbseq', 'output'); os.makedirs(OUT, exist_ok=True)
RED, TEAL = '#C0392B', '#1F5FA8'   # the red and blue of Fig1
N_SHOW = 20
## MARK = knock-outs drawn with a diamond. The panel of the paper has none: the two regulators of glucocorticoid-receptor
## signalling in the DG lists (loss in MICE alters GR signalling or HPA-axis function: Ncoa1, Winnay 2006, Lachize 2009; Crebrf,
## Martyn 2012) are named in the Results. Kmt2d (aligned 1st) stabilises GR only in cell lines (Wu 2024).
MARK = set()
CT_NAME = {'DG': 'DG Glut', 'CA1': 'CA1 Glut', 'CA3': 'CA3 Glut', 'GABA': 'GABA'}   # the names of De Donno et al.

rows = [r for r in csv.DictReader(open(os.path.join(BASE, 'analysis', 'perturbseq', 'results', f'axis_{AX}_vs_{CTX}.csv')))
        if r['kind'] == 'KO']
rows.sort(key=lambda r: int(r['rank_among_KO']))                    # rank_among_KO 1 = most aligned (pipeline)
## the tables number the knock-outs from the OPPOSED end (opposed at the top):
## # 1 = the knock-out whose effect is most opposed to the stress response; the opposed table comes first
N_KO = len(rows)
for r in rows:
    r['rank_opp'] = str(N_KO - int(r['rank_among_KO']) + 1)
print(f'{len(rows):,} KOs; aligned: ' + ', '.join(r['target'] for r in rows[:N_SHOW]))
print('opposed: ' + ', '.join(r['target'] for r in rows[::-1][:N_SHOW]))

## geometry in mm, at the size of the panel in the final figure (row : bar : header heights 1 : 1.42 : 1.26)
MM = 1 / 25.4
FIG_W_MM, FIG_H_MM, DPI = 43.7, 59.5, 600   # size in the final figure
FS_BAR, FS_HEAD, FS_ROW = 5.0, 5.0, 5.5      # pt (dense table); at 5.0 pt bold the title bar fits the panel width (FIG_W_MM) on one line
ROW_H = 2.5 * MM * FIG_H_MM / 57.4           # row height grows with the panel (2.5 mm at 57.4 mm); text stays in pt
BAR_H, HDR_H = ROW_H * 1.42, ROW_H * 1.26
PAD_X, PAD_Y = 0.3 * MM, 0.35 * MM           # outer margin, left/right and top/bottom
RULE_LW = 0.5                                # pt
COLW = [.13, .49, .38]

def draw(ax, data, title, colour, AX_H):
    ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fbar, fhdr, frow = BAR_H / AX_H, HDR_H / AX_H, ROW_H / AX_H
    ax.add_patch(Rectangle((0, 1 - fbar), 1, fbar, facecolor=colour, edgecolor='none', transform=ax.transAxes, clip_on=False))
    ax.text(.5, 1 - fbar / 2, title, ha='center', va='center', color='white', fontsize=FS_BAR, fontweight='bold',
            transform=ax.transAxes)
    x = 0.0
    for w, h in zip(COLW, ['#', 'Knock-out', 'RFind score']):
        ax.text(x + w / 2, 1 - fbar - fhdr * .40, h, ha='center', va='center', fontsize=FS_HEAD, fontweight='bold', color='#333333')
        x += w
    rule = 1 - fbar - fhdr * .78
    ax.plot([0, 1], [rule, rule], color='#333333', lw=RULE_LW, clip_on=False)
    y = 1 - fbar - fhdr - frow / 2
    for i, r in enumerate(data):
        if i % 2 == 0:
            ax.add_patch(Rectangle((0, y - frow / 2), 1, frow, facecolor='#F2F2F2', edgecolor='none', zorder=0))
        vals = [r['rank_opp'], r['target'], f"{float(r['rfind']):+.2f}"]
        x = 0.0
        for j, (w, v) in enumerate(zip(COLW, vals)):
            ax.text(x + w / 2, y, v, ha='center', va='center', fontsize=FS_ROW, fontweight='bold' if j == 1 else 'normal',
                    color=colour if j == 1 else '#222222', zorder=2)
            if j == 1 and r['target'] in MARK:   # diamond left of the name, as in Fig4d
                ax.plot([x + w / 2 - 0.035 - 0.0115 * len(v)], [y], marker='D', markersize=FS_ROW * 0.55, color='#222222', transform=ax.transAxes, zorder=3, clip_on=False)
            x += w
        y -= frow

def one(data, title, colour, path):
    AX_H = BAR_H + HDR_H + len(data) * ROW_H
    ## whole pixels: round(mm / 25.4 * 600) (+0.1 px so that matplotlib does not floor to one pixel less)
    W, H = [(round(v * MM * DPI) + 0.1) / DPI for v in (FIG_W_MM, FIG_H_MM)]
    fig = plt.figure(figsize=(W, H)); fig.patch.set_facecolor('white')
    ax = fig.add_axes([PAD_X / W, (H - PAD_Y - AX_H) / H, 1 - 2 * PAD_X / W, AX_H / H])
    draw(ax, data, title, colour, AX_H)
    fig.savefig(path, dpi=DPI, facecolor='white'); plt.close(fig)   # no bbox_inches='tight': the size is the final size
    print('wrote', os.path.relpath(path, BASE))

one(rows[::-1][:N_SHOW], f'MOST OPPOSED to the {CT_NAME.get(AX, AX)} stress response', TEAL, os.path.join(OUT, f'Fig4c_{AX}_{CTX}_opposed.png'))
one(rows[:N_SHOW], f'MOST ALIGNED with the {CT_NAME.get(AX, AX)} stress response', RED, os.path.join(OUT, f'Fig4c_{AX}_{CTX}_aligned.png'))
