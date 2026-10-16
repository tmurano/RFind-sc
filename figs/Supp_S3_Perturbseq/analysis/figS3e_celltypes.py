# -*- coding: utf-8 -*-
"""Supplementary Fig. S3e — how often the pairs that Shi et al. report as similar come back, in three cell types.
Bars: the thalamus, the hippocampus (the Delta of the dentate gyrus, CA1 and CA3 summed) and the GABAergic neurons; their values
are also written to figS3e_celltypes.txt.
Measure: of the pairs of knock-outs that Shi et al. name together as similar in the main text (figS3_shi_pairs.py; a pair counts when both
knock-outs are in the data), the share that falls in the top 1% of all knock-out pairs of that cell type (rank by RFind score, as in
figS3_rfind_panels.py); by chance 1%. Cell types: the thalamus (b-d), the GABAergic neurons of Fig. 4 and the hippocampus = dentate gyrus,
CA1 and CA3 with the Delta of each cell type summed (figS3_rfind_allpairs.py HIP).
Input : figs/Supp_S3_Perturbseq/output/data/rfind_allpairs_{g151,GABA,HIP}.npz
Output: figs/Supp_S3_Perturbseq/output/FigS3e_celltypes.png and figS3e_celltypes.txt
Size: drawn directly at its size in the whole figure (180 mm wide in the slide), 53.3 x 44.7 mm, 600 dpi, axes
placed in mm; Arial, tick labels 6 pt, axis title 7 pt, the numbers over the bars and the key of the chance line 6 pt; axes and ticks
0.5 pt.
Run   : python3 figs/Supp_S3_Perturbseq/analysis/figS3e_celltypes.py   (project root)
"""
import os, sys
from itertools import combinations
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
OUT = os.path.join(BASE, 'figs', 'Supp_S3_Perturbseq', 'output')
sys.path.insert(0, HERE)
from figS3_shi_pairs import SHI
RED = '#C0392B'
## bars: thalamus, hippocampus, GABAergic neurons; the numbers of nuclei are given in the figure legend only
SHOWN = [('g151', 'TH Prkcd\nGrin2c Glut\n(thalamus)'), ('HIP', 'Hippocampus'), ('GABA', 'GABA')]   # names of Shi et al. / De Donno et al.
FT, FL, FN = 6, 7, 6                                 # tick labels / axis title / notes (pt)
MM, DPI = 1 / 25.4, 600
plt.rcParams.update({'font.family': 'Arial', 'font.size': FL, 'axes.labelsize': FL, 'xtick.labelsize': FT, 'ytick.labelsize': FT,
                     'savefig.dpi': DPI, 'figure.dpi': DPI, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5,
                     'ytick.major.width': 0.5, 'xtick.major.size': 2.0, 'ytick.major.size': 2.0, 'xtick.major.pad': 1.2,
                     'ytick.major.pad': 1.2, 'axes.labelpad': 1.5})
PPTX = OUT                                           # panel at its size in the whole figure, in output/ (put together
                                                     # in the slide)
FW, FH = 53.3, 44.7                                  # mm; the size of e in the whole figure, so that the text is 6 / 7 pt there

def shi_top1(key):
    d = np.load(os.path.join(OUT, 'data', f'rfind_allpairs_{key}.npz'), allow_pickle=True)
    S, t, k, n = d['S'].astype(float), d['targets'].astype(str), d['kind'].astype(str), d['nuclei']
    ko = np.where(k == 'KO')[0]; tk = t[ko]
    iu = np.triu_indices(len(ko), 1); Sk = S[np.ix_(ko, ko)][iu]
    order = np.argsort(-Sk, kind='stable'); rpos = np.empty_like(order); rpos[order] = np.arange(1, len(Sk) + 1)
    ix = {g: i for i, g in enumerate(tk)}; N = len(tk)
    seen, top = set(), []
    for _, genes in SHI:
        for a, b in combinations([g for g in genes if g in ix], 2):
            i, j = sorted((ix[a], ix[b]))
            if (i, j) in seen:
                continue
            seen.add((i, j)); q = i * N - i * (i + 1) // 2 + (j - i - 1)      # index of (i, j) in the upper triangle (row-major)
            top.append(100 * rpos[q] / len(Sk) <= 1)
    return int(np.sum(top)), len(top), int(np.median(n[ko]))

res = {key: shi_top1(key) for key, _ in SHOWN}
lines = ['Pairs named as similar by Shi et al. that fall in the top 1% of all knock-out pairs of the cell type (by chance 1%).']
for key, name in SHOWN:
    a, n, m = res[key]
    lines.append(f'{name.replace(chr(10), " "):14s} {a:3d}/{n:3d} = {100 * a / n:5.1f}%   median nuclei per knock-out {m}')

fig = plt.figure(figsize=(round(FW * MM * DPI) / DPI, round(FH * MM * DPI) / DPI))
ax = fig.add_axes([12.5 / FW, 7.5 / FH, (FW - 13.0) / FW, (FH - 7.5 - 1.5) / FH])
xs = np.arange(len(SHOWN))
ax.bar(xs, [100 * res[k][0] / res[k][1] for k, _ in SHOWN], width=0.62, color=RED)
for x, (k, _) in zip(xs, SHOWN):
    a, n, _m = res[k]
    ax.text(x, 100 * a / n + 1.5, f'{100 * a / n:.0f}%\n({a}/{n})', ha='center', va='bottom', fontsize=FN)
ax.axhline(1, color='black', lw=0.5, ls='--', label='chance (1%)')
## the chance line is named in a key at the top right (at 53 mm wide there is no room right of the bars for the text,
## and the x labels need the whole width)
ax.legend(loc='upper right', frameon=False, fontsize=FN, handlelength=1.8, borderaxespad=0.2, labelcolor='#4D4D4D')
ax.set_xticks(xs); ax.set_xticklabels([name for k, name in SHOWN], fontsize=FT)
ax.set_xlim(-0.5, len(SHOWN) - 0.5); ax.set_ylim(0, 90); ax.set_yticks(range(0, 91, 20))   # every 20
ax.set_ylabel('Pairs reported by Shi et al.\nin the top 1% of all pairs (%)')
ax.spines[['top', 'right']].set_visible(False)
os.makedirs(PPTX, exist_ok=True)
for ext in ('png',):  # PNG only
    fig.savefig(os.path.join(PPTX, f'FigS3e_celltypes.{ext}'), facecolor='white')
open(os.path.join(OUT, 'figS3e_celltypes.txt'), 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines)); print('wrote FigS3e_celltypes.png')
