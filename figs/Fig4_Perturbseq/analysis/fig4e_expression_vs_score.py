# -*- coding: utf-8 -*-
"""Fig4e — expression of Nr3c1 and Ncoa1 in the No Stress cells of De Donno et al. and the RFind score of the knock-out of the same
gene, in the four cell types of Fig. 4 (DG Glut, CA1 Glut, CA3 Glut, GABA).

  <python> figs/Fig4_Perturbseq/analysis/fig4e_expression_vs_score.py

Input : data/Perturbseq/axis_dedonno_<cell type>_stress.csv   (analysis/perturbseq/02_stress_axis.py)
          x = CPM_NoStress (CPM of the gene in the summed counts of the No Stress cells of that cell type)
        analysis/perturbseq/results/axis_<cell type>_vs_<context>.csv   (analysis/perturbseq/03_axis_vs_ko.py)
          y = rfind (RFind score of the knock-out against the De Donno stress response of the same cell type;
              negative = opposed to the stress response)
Output: figs/Fig4_Perturbseq/output/Fig4e_expression_vs_score.png   (82 x 40 mm, two plots, 600 dpi, Arial)
        figs/Fig4_Perturbseq/output/fig4e_expression_vs_score.csv        (the two genes x the four cell types)
        figs/Fig4_Perturbseq/output/fig4e_expression_vs_score_stats.csv  (Spearman rho per gene, descriptive; no test with four
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
    matplotlib.rcParams['mathtext.rm'] = 'Arial:italic'      # only the gene names are mathtext: digits italic too
    matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
RES = os.path.join(BASE, 'analysis', 'perturbseq', 'results')
DD = os.path.join(BASE, 'data', 'Perturbseq')
OUT = os.path.join(BASE, 'figs', 'Fig4_Perturbseq', 'output'); os.makedirs(OUT, exist_ok=True)
GENES = ['Nr3c1', 'Ncoa1']
## cell type of Fig. 4 -> Perturb-seq context of the axis file (as in fig4d_go_set.R)
CTX = {'DG': 'g017_DG', 'CA1': 'CA1', 'CA3': 'g017_CA3', 'GABA': 'GABA3'}
NAME = {'DG': 'DG Glut', 'CA1': 'CA1 Glut', 'CA3': 'CA3 Glut', 'GABA': 'GABA'}   # labels = the names of De Donno et al.
COL = {'DG': '#0072B2', 'CA1': '#E69F00', 'CA3': '#009E73', 'GABA': '#CC79A7'}    # Okabe-Ito colours, distinct for colour-blind readers
FIG_W_MM, FIG_H_MM = 82.0, 40.0   # size in the final figure
MM = 1 / 25.4
RIGHT = {('Ncoa1', 'GABA'): True}   # label to the right of the point, away from the CA3 Glut label

rows = []
for ct, ctx in CTX.items():
    ko = {r['target']: r for r in csv.DictReader(open(os.path.join(RES, f'axis_{ct}_vs_{ctx}.csv'))) if r['kind'] == 'KO'}
    dd = {r['Symbol']: r for r in csv.DictReader(open(os.path.join(DD, f'axis_dedonno_{ct}_stress.csv')))}
    for g in GENES:
        rows.append(dict(gene=g, cell_type=ct, cpm_no_stress=float(dd[g]['CPM_NoStress']),
                         rfind_score=float(ko[g]['rfind']), ko_nuclei=int(ko[g]['n_cells'])))


def spearman(a, b):
    ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


stats = []
for g in GENES:
    sub = [r for r in rows if r['gene'] == g]
    stats.append(dict(gene=g, x=f'{g} CPM in the No Stress cells of De Donno et al.', y=f'RFind score of the {g} knock-out',
                      statistic='Spearman rho', value=round(spearman([r['cpm_no_stress'] for r in sub], [r['rfind_score'] for r in sub]), 4),
                      n_cell_types=len(sub)))
with open(os.path.join(OUT, 'fig4e_expression_vs_score.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
with open(os.path.join(OUT, 'fig4e_expression_vs_score_stats.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(stats[0])); w.writeheader(); w.writerows(stats)

## ---- the panel: two plots side by side (left Nr3c1, right Ncoa1) ----
plt.rcParams.update({'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
                     'xtick.major.size': 2, 'ytick.major.size': 2, 'xtick.labelsize': 6, 'ytick.labelsize': 6})
fig = plt.figure(figsize=(FIG_W_MM * MM, FIG_H_MM * MM))
LIM = {'Nr3c1': dict(x=(0, 140), xt=[0, 40, 80, 120], y=(-30, 5), yt=[-30, -20, -10, 0]),
       'Ncoa1': dict(x=(70, 140), xt=[80, 100, 120, 140], y=(-35, 5), yt=[-30, -20, -10, 0])}
for k, g in enumerate(GENES):
    ax = fig.add_axes([0.125 + 0.5 * k, 0.24, 0.35, 0.71])
    ax.axhline(0, color='0.6', lw=0.5, ls=(0, (3, 2)), zorder=1)
    for r in (q for q in rows if q['gene'] == g):
        ax.scatter(r['cpm_no_stress'], r['rfind_score'], s=12, color=COL[r['cell_type']], lw=0, zorder=3)
        right = r['cpm_no_stress'] < np.mean([q['cpm_no_stress'] for q in rows if q['gene'] == g])
        near0 = -3 < r['rfind_score'] < 0                       # a label just below the dashed line of 0 is moved down
        right = RIGHT.get((g, r['cell_type']), right)
        ax.annotate(NAME[r['cell_type']], (r['cpm_no_stress'], r['rfind_score']), xytext=(4 if right else -4, -1 if near0 else 0),
                    textcoords='offset points', ha='left' if right else 'right', va='top' if near0 else 'center', fontsize=7)
    L = LIM[g]
    ax.set_xlim(*L['x']); ax.set_xticks(L['xt']); ax.set_ylim(*L['y']); ax.set_yticks(L['yt'])
    ax.set_xlabel(r'$\mathrm{' + g + '}$ expression in' + '\nNo Stress cells (CPM)', fontsize=7, labelpad=1.5)
    ax.set_ylabel(r'RFind score of' + '\n' + r'$\mathrm{' + g + '}$ knock-out', fontsize=7, labelpad=1.5)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.tick_params(pad=1.5)
## the Spearman rho is not drawn; it is written only to the stats CSV (Source Data)
png = os.path.join(OUT, 'Fig4e_expression_vs_score.png')
fig.savefig(png, dpi=600)
print('wrote', png)
for r in rows:
    print(r)
for s in stats:
    print(s)
