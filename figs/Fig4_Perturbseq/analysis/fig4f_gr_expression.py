# -*- coding: utf-8 -*-
"""Fig4f — glucocorticoid-receptor (Nr3c1) expression and the RFind score of the Nr3c1 knock-out, in the four cell
types of Fig. 4 (CA1, CA3, DG, GABA).

  <python> figs/Fig4_Perturbseq/analysis/fig4f_gr_expression.py

Input : analysis/perturbseq/results/axis_<cell type>_vs_<context>.csv   (analysis/perturbseq/03_axis_vs_ko.py)
        x = target_cpm_in_non (CPM of Nr3c1 in the non-targeting control column of that cell type)
        y = rfind (RFind score of the Nr3c1 knock-out against the De Donno stress response of the same cell type;
            negative = opposed to the stress response)
Output: figs/Fig4_Perturbseq/output/Fig4f_GR_expression_vs_score.png   (41 x 40 mm, 600 dpi, Arial)
        figs/Fig4_Perturbseq/output/fig4f_GR_expression_vs_score.csv        (the four cell types)
        figs/Fig4_Perturbseq/output/fig4f_GR_expression_vs_score_stats.csv  (Spearman rho, descriptive; no test with four
        cell types, as in Methods)
"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

if any(f.name == 'Arial' for f in font_manager.fontManager.ttflist):
    matplotlib.rcParams['font.family'] = 'Arial'
    matplotlib.rcParams['mathtext.fontset'] = 'custom'
    matplotlib.rcParams['mathtext.rm'] = 'Arial'
    matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
RES = os.path.join(BASE, 'analysis', 'perturbseq', 'results')
OUT = os.path.join(BASE, 'figs', 'Fig4_Perturbseq', 'output'); os.makedirs(OUT, exist_ok=True)
GENE = 'Nr3c1'
## cell type of Fig. 4 -> Perturb-seq context of the axis file (as in fig4d_heatmap.R)
CTX = {'CA1': 'CA1', 'CA3': 'g017_CA3', 'DG': 'g017_DG', 'GABA': 'GABA3'}
NAME = {'CA1': 'CA1 Glut', 'CA3': 'CA3 Glut', 'DG': 'DG Glut', 'GABA': 'GABA'}   # labels = the names of De Donno et al.
FIG_W_MM, FIG_H_MM = 41.0, 40.0   # size in the final figure
MM = 1 / 25.4

rows = []
for ct, ctx in CTX.items():
    r = [r for r in csv.DictReader(open(os.path.join(RES, f'axis_{ct}_vs_{ctx}.csv')))
         if r['kind'] == 'KO' and r['target'] == GENE]
    assert len(r) == 1, (ct, len(r))
    r = r[0]
    rows.append(dict(cell_type=ct, perturbseq_context=ctx, nr3c1_cpm_non_targeting=float(r['target_cpm_in_non']),
                     rfind_score=float(r['rfind']), ko_nuclei=int(r['n_cells'])))
x = np.array([r['nr3c1_cpm_non_targeting'] for r in rows])
y = np.array([r['rfind_score'] for r in rows])

def spearman(a, b):
    ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])

stats = [dict(x='Nr3c1 CPM in the non-targeting controls', y='RFind score of the Nr3c1 knock-out', statistic='Spearman rho',
              value=round(spearman(x, y), 4), n_cell_types=len(x))]

with open(os.path.join(OUT, 'fig4f_GR_expression_vs_score.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
with open(os.path.join(OUT, 'fig4f_GR_expression_vs_score_stats.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(stats[0])); w.writeheader(); w.writerows(stats)

## ---- the panel ----
plt.rcParams.update({'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
                     'xtick.major.size': 2, 'ytick.major.size': 2, 'xtick.labelsize': 6, 'ytick.labelsize': 6})
fig = plt.figure(figsize=(FIG_W_MM * MM, FIG_H_MM * MM))
ax = fig.add_axes([0.25, 0.24, 0.70, 0.71])
ax.axhline(0, color='0.6', lw=0.5, ls=(0, (3, 2)), zorder=1)
ax.scatter(x, y, s=9, color='black', lw=0, zorder=3)
## labels beside the points (CA3 to the right, the others to the left)
for r in rows:
    right = r['cell_type'] == 'CA3'
    ax.annotate(NAME[r['cell_type']], (r['nr3c1_cpm_non_targeting'], r['rfind_score']), xytext=(4 if right else -4, 0),
                textcoords='offset points', ha='left' if right else 'right', va='center', fontsize=7)
ax.set_xlim(0, 130); ax.set_xticks([0, 40, 80, 120])
ax.set_ylim(-30, 5); ax.set_yticks([-30, -20, -10, 0])
ax.set_xlabel(r'$\it{Nr3c1}$ expression in' + '\nnon-targeting controls (CPM)', fontsize=7, labelpad=1.5)
ax.set_ylabel(r'RFind score of' + '\n' + r'$\it{Nr3c1}$ knock-out', fontsize=7, labelpad=1.5)
rho = stats[0]['value']
## the Spearman rho is not drawn; it is written only to the stats CSV (Source Data)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
ax.tick_params(pad=1.5)
png = os.path.join(OUT, 'Fig4f_GR_expression_vs_score.png')
fig.savefig(png, dpi=600)
print('wrote', png)
for r in rows:
    print(r)
for s in stats:
    print(s)
