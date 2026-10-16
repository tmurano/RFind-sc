# -*- coding: utf-8 -*-
"""Supplementary Fig. S3 b, c, d — RFind score (4-term Running Fisher) between all pairs of columns (Delta), thalamic glutamatergic neurons.

The all-pairs comparison of the columns is a plain Running Fisher comparison (the step is labelled "Running Fisher" in a); the
score keeps the name RFind score.
The default context is g151 = '151 TH Prkcd Grin2c Glut', the one cell type that Shi et al. name in their main-text comparison
of knock-outs (Fig. 2H, 2I, 3C), so that the knock-out pairs they describe as similar are compared in the same cell type.
Fold change is not shown (Delta only).
  b  heat map of the RFind score among the knock-outs that have at least one knock-out partner above the highest
     safe-targeting x safe-targeting pair (the dashed line of c), clustered (average linkage on 1 - Pearson correlation of
     their rows); brackets = the known complexes that appear in d (proteasome: ticks on the left); flat colour bar, values
     beyond +-40 in the end colours; the title gives the number of knock-outs shown
  c  distributions by type of pair (KO x KO, KO x safe, safe x safe; three colours), and the right end (above +30) enlarged in
     its own panel below, linked to a box on the right end, with the safe x safe pairs as one tick per pair; dashed line =
     the highest safe x safe pair, labelled with its value; the text gives how many knock-out pairs lie above it and how
     many would if they were distributed like the control pairs (n_KO pairs / (n_control pairs + 1); a description of the
     tail, not a test of single pairs)
  d  the 20 most similar knock-out pairs (no statistical test), after each bar the complex or pathway that both knock-outs
     belong to; red = both knock-outs named together as similar in the main text of Shi et al.
Input : figs/Supp_S3_Perturbseq/output/data/rfind_allpairs_<context>.npz (figS3_rfind_allpairs.py <n> delta <context>)
Output: figs/Supp_S3_Perturbseq/output/FigS3{b_delta_heatmap,c_delta_distribution,d_top_pairs}.png (g151 only),
        figS3_rfind_pairs_<context>.txt (numbers for the text, in output/)
Each panel is drawn at its size in the whole figure (180 mm wide): b 97 x 62.5 mm, c 75 x 62.5 mm,
d 98 x 70 mm, 600 dpi, axes placed in mm; Arial, tick labels 6 pt, axis titles, titles and legends 7 pt, notes and labels of
complexes 6 pt; axes and ticks 0.5 pt.
Run   : python3 figs/Supp_S3_Perturbseq/analysis/figS3_rfind_panels.py [g151|DG|CA1|CA3|GABA|HIP]   (from the project root)
"""
import os, sys
from itertools import combinations
import numpy as np
from scipy.stats import gaussian_kde
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
OUT = os.path.join(BASE, 'figs', 'Supp_S3_Perturbseq', 'output')
KOC, SAFEC = '#8C8C8C', '#F0A030'                    # knock-out / safe-targeting
RED, BLUE = '#C0392B', '#1F5FA8'                     # the paper's red / blue
CMAP = LinearSegmentedColormap.from_list('rb', [BLUE, '#AFC9E8', 'white', '#F2B8B0', RED])
VL = 20.0                                            # colour scale of the heat map (RFind score)
CTX = sys.argv[1] if len(sys.argv) > 1 else 'g151'
sys.path.insert(0, HERE)
from figS3_shi_pairs import SHI                      # pairs named as similar by Shi et al. (main text)
## sizes at the pasted size: text in pt, lengths in mm
FT, FL, FN = 6, 7, 6                                 # tick labels / axis titles, titles, legends / notes
MM, DPI = 1 / 25.4, 600
plt.rcParams.update({'font.family': 'Arial', 'font.size': FL, 'axes.labelsize': FL, 'xtick.labelsize': FT,
                     'ytick.labelsize': FT, 'legend.fontsize': FL, 'savefig.dpi': DPI, 'figure.dpi': DPI,
                     'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2.0,
                     'ytick.major.size': 2.0, 'xtick.major.pad': 1.2, 'ytick.major.pad': 1.2, 'axes.labelpad': 1.5,
                     'axes.titlepad': 2.0})
PPTX = OUT                                           # panels at their size in the whole figure, in output/ (put together
                                                     # in the slide)
SIZE = {'b': (97.0, 62.5), 'c': (75.0, 62.5), 'd': (98.0, 70.0)}   # mm (the slide places them 1:1)


def mm_fig(w, h):
    """A figure of w x h mm (rounded to whole pixels)."""
    fig = plt.figure(figsize=(round(w * MM * DPI) / DPI, round(h * MM * DPI) / DPI))
    fig.W, fig.H = w, h
    return fig


def mm_axes(fig, x, y, w, h):
    """Axes placed in mm from the lower left corner of the figure."""
    return fig.add_axes([x / fig.W, y / fig.H, w / fig.W, h / fig.H])

d = np.load(os.path.join(OUT, 'data', f'rfind_allpairs_{CTX}.npz'), allow_pickle=True)
S, tgt, kind, nuc = d['S'].astype(np.float64), d['targets'].astype(str), d['kind'].astype(str), d['nuclei']
ko, sf = kind == 'KO', kind == 'safe'
tri = lambda m: S[np.ix_(m, m)][np.triu_indices(int(m.sum()), 1)]
cls = {'Knock-out × knock-out': tri(ko), 'Knock-out × safe-targeting': S[np.ix_(ko, sf)].ravel(),
       'Safe-targeting × safe-targeting': tri(sf)}
kk, ss = cls['Knock-out × knock-out'], cls['Safe-targeting × safe-targeting']
smax = ss.max()
n_above, n_exp = int((kk > smax).sum()), len(kk) / (len(ss) + 1)

def shi_of(a, b):
    return [f for f, g in SHI if a in g and b in g]

iu = np.triu_indices(int(ko.sum()), 1); Sk = S[np.ix_(ko, ko)][iu]; tk = tgt[ko]
top = np.argsort(-Sk, kind='stable')[:20]
TOP = [(tk[iu[0][k]], tk[iu[1][k]], Sk[k]) for k in top]
lines = [f'RFind score (4-term Running Fisher) between all pairs of columns, {CTX} ({int(ko.sum()):,} knock-outs + {int(sf.sum())} safe-targeting groups; '
         f'median nuclei knock-out {int(np.median(nuc[ko]))}, safe-targeting {int(np.median(nuc[sf]))}).']
for k, v in cls.items():
    lines.append(f'{k:32s} mean {v.mean():+7.3f}  median {np.median(v):+7.3f}  sd {v.std():6.3f}  positive {(v > 0).mean():.3f}  '
                 f'min {v.min():+8.1f}  max {v.max():+8.1f}  n {len(v):,}')
lines.append(f'highest safe x safe pair {smax:+.1f}; knock-out pairs above it {n_above:,} of {len(kk):,} '
             f'({100 * n_above / len(kk):.2f}%); expected if distributed like the control pairs {n_exp:.0f} -> {n_above / n_exp:.1f} times')
lines.append('20 most similar knock-out pairs (rank, pair, score, statement of Shi et al.):')
for r, (a, b, v) in enumerate(TOP, 1):
    lines.append(f'  {r:2d} {a} x {b}  {v:+.1f}  {", ".join(shi_of(a, b)) or "-"}')
## every pair named by Shi et al.: score and rank among all knock-out pairs
order_all = np.argsort(-Sk, kind='stable'); rpos = np.empty_like(order_all); rpos[order_all] = np.arange(1, len(Sk) + 1)
ix = {t: i for i, t in enumerate(tk)}
key = lambda a, b: (min(ix[a], ix[b]), max(ix[a], ix[b]))
lines.append('pairs named by Shi et al. (score, rank among all knock-out pairs, top %):')
seen = set()
for f, genes in SHI:
    have = [g for g in genes if g in ix]
    miss = [g for g in genes if g not in ix]
    for a, b in combinations(have, 2):
        i, j = key(a, b)
        if (i, j) in seen:
            continue
        seen.add((i, j))
        n = len(tk); k = i * n - i * (i + 1) // 2 + (j - i - 1)          # index of (i, j) in the upper triangle (row-major)
        lines.append(f'  {f:9s} {a:>8} x {b:<8} {Sk[k]:+8.1f}  rank {rpos[k]:>9,}  top {100 * rpos[k] / len(Sk):6.3f}%')
    if miss:
        lines.append(f'  {f:9s} not in the data: {", ".join(miss)}')
open(os.path.join(OUT, f'figS3_rfind_pairs_{CTX}.txt'), 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))

TAG = '' if CTX == 'g151' else f'_{CTX}'                   # other contexts do not overwrite the panels of the figure
def save(fig, fname):
    fname += TAG
    os.makedirs(PPTX, exist_ok=True)
    for ext in ('png',):  # PNG only
        fig.savefig(os.path.join(PPTX, f'{fname}.{ext}'), facecolor='white')
    plt.close(fig); print(f'wrote {fname}.png')

if CTX != 'g151':                     # other cell types: the numbers only, no panels
    open(os.path.join(OUT, f'figS3_rfind_pairs_{CTX}.txt'), 'w').write('\n'.join(lines) + '\n')
    sys.exit(0)

## ---- b: knock-outs with at least one partner above the highest control pair, clustered ----
## only the knock-outs with a partner above the highest safe-targeting pair (the dashed line of c), average linkage on
## 1 - Pearson correlation of their rows; brackets = known complexes; flat colour bar (values beyond +-VLB in the end colours).
VLB = 40.0
Skk = S[np.ix_(ko, ko)].copy(); np.fill_diagonal(Skk, -np.inf)
strong = np.flatnonzero((Skk > smax).any(1))
tks = tgt[ko][strong]
Xs = S[np.ix_(np.flatnonzero(ko)[strong], np.flatnonzero(ko)[strong])].copy(); np.fill_diagonal(Xs, 0.0)
Rs = np.corrcoef(Xs); Ds = np.clip(1 - Rs, 0, None); np.fill_diagonal(Ds, 0.0)
ordr = hierarchy.leaves_list(hierarchy.linkage(squareform(Ds, checks=False), method='average'))
M = Xs[np.ix_(ordr, ordr)].copy(); np.fill_diagonal(M, np.nan)
pos = {g: i for i, g in enumerate(tks[ordr])}
## only the complexes that appear in d, named as in d
MODULES = [('UFMylation (E1–E2)', ['Uba5', 'Ufc1']), ('AP-2 adaptor', ['Ap2s1', 'Ap2m1']),
           ('TSC complex', ['Tsc1', 'Tsc2']), ('SWI/SNF', ['Smarcb1', 'Smarcc2', 'Smarce1']), ('U5 snRNP', ['Ddx23', 'Prpf6']),
           ('Elongator', ['Elp1', 'Elp2'])]
lines.append(f'b: {len(strong)} knock-outs with at least one knock-out partner above the highest safe x safe pair ({smax:+.1f})')
FW, FH = SIZE['b']; LM, HM, TM = 14.0, 58.0, 4.0                    # left margin, heat map, top margin (mm)
fig = mm_fig(FW, FH)
axm = mm_axes(fig, LM, FH - TM - HM, HM, HM)
im = axm.imshow(M, cmap=CMAP, vmin=-VLB, vmax=VLB, interpolation='nearest', aspect='auto', rasterized=True, extent=(0, len(M), len(M), 0))
axm.set_xticks([]); axm.set_yticks([])
## the title gives the number of knock-outs shown; the figure legend says which knock-outs these are
axm.set_title(f'{len(strong)} of {int(np.sum(ko)):,} knock-outs', fontsize=FL, pad=2.5)
for name, genes in MODULES:
    p_ = sorted(pos[g] for g in genes if g in pos)
    if len(p_) < 2 or p_[-1] - p_[0] > 0.06 * len(M):
        lines.append(f'   module {name}: positions {p_} (not bracketed)'); continue
    y0, y1 = p_[0], p_[-1] + 1
    axm.plot([len(M) * 1.015] * 2, [y0, y1], color='black', lw=0.7, clip_on=False)
    axm.text(len(M) * 1.03, (y0 + y1) / 2, name, va='center', ha='left', fontsize=FN, clip_on=False)
    lines.append(f'   module {name}: rows {y0}-{y1 - 1}')
## the proteasome knock-outs fall in two places (Psmb4 and Pomp next to the V-ATPase, Psmc5 and Psmc1 next to the U5 snRNP):
## ticks and a label on the left, where there is room
for gs in (['Psmb4', 'Pomp'], ['Psmc5', 'Psmc1']):
    p_ = [pos[g] for g in gs if g in pos]
    for q in p_:
        axm.plot([-len(M) * 0.03, -len(M) * 0.005], [q + 0.5] * 2, color='black', lw=0.6, clip_on=False)
    if p_:
        axm.text(-len(M) * 0.04, np.mean(p_) + 0.5, 'proteasome', va='center', ha='right', fontsize=FN, clip_on=False)
    lines.append(f'   proteasome ticks: {dict((g, pos.get(g)) for g in gs)}')
cax = mm_axes(fig, LM + HM + 14.5, FH - TM - HM + 0.17 * HM, 2.2, 0.66 * HM)   # between the labels of the complexes
cb = fig.colorbar(im, cax=cax); cb.set_label('RFind score', labelpad=1.5); cb.set_ticks([-40, -20, 0, 20, 40])
cb.outline.set_linewidth(0.5); cb.ax.tick_params(width=0.5, length=2.0)
axm.set_xlim(0, len(M)); axm.set_ylim(len(M), 0)            # the labels outside must not widen the axes
save(fig, 'FigS3b_delta_heatmap')

## ---- c: distributions by type of pair (linear); the right end enlarged in a separate panel below ----
## the right end is enlarged in its own panel under the density (same colours and lines, linked to a box on the right end);
## the three types of pair have three colours (KO x safe = blue dashed).
from matplotlib.patches import Rectangle, ConnectionPatch
style = {'Knock-out × knock-out': ('#4D4D4D', 0.9, '-'), 'Knock-out × safe-targeting': ('#3A8FD0', 0.75, '--'),
         'Safe-targeting × safe-targeting': (SAFEC, 0.95, '-')}
SHORT = {'Knock-out × knock-out': 'KO × KO', 'Knock-out × safe-targeting': 'KO × safe', 'Safe-targeting × safe-targeting': 'safe × safe'}
## bins of 1 (main) and 2.5 (right end), drawn as lines through the bin centres; the dashed line is labelled as the highest
## safe x safe pair with its value (smax, computed from the data), and the axis runs past it.
FW, FH = SIZE['c']
fig = mm_fig(FW, FH)
ax = mm_axes(fig, 10.0, 33.5, FW - 11.5, 27.5); axt = mm_axes(fig, 10.0, 8.2, FW - 11.5, 13.0)   # density | right end (mm)
def dens(v, edges):
    h, _ = np.histogram(v, bins=edges); return (edges[:-1] + edges[1:]) / 2, h / (len(v) * np.diff(edges))
bins = np.arange(-40, 47, 1.0)
for lab, v in cls.items():
    c, lw, ls = style[lab]; x_, y_ = dens(v, bins)
    ax.plot(x_, y_, color=c, lw=lw, ls=ls, label=SHORT[lab])
ax.axvline(smax, color='black', lw=0.5, ls='--')
ax.text(smax - 1.2, ax.get_ylim()[1] * 0.62, f'highest\nsafe × safe\n({smax:.2f})', ha='right', va='center', fontsize=FN, color='#333333')
ax.set_xlim(-42, 47); ax.set_ylim(0, None); ax.set_xlabel('RFind score'); ax.set_ylabel('Density')
ax.set_xticks([-40, -20, 0, 20, 40])                        # fixed ticks (the narrow panel would otherwise get more)
ax.spines[['top', 'right']].set_visible(False)
ax.legend(loc='upper left', frameon=False, fontsize=FL, handlelength=1.8, borderaxespad=0.1, labelspacing=0.3)
zb = np.arange(30, 152.5, 2.5)
for lab, v in cls.items():                                  # safe x safe has only a few pairs here: one tick per pair instead of a density
    c, lw, ls = style[lab]
    if lab == 'Safe-targeting × safe-targeting':
        vv = v[v >= 30]; axt.plot(vv, np.zeros_like(vv), '|', color=c, ms=7, mew=0.9, clip_on=False, zorder=4)
        continue
    x_, y_ = dens(v, zb); axt.plot(x_, y_, color=c, lw=lw, ls=ls)
axt.axvline(smax, color='black', lw=0.5, ls='--')
axt.set_xlim(30, 150); axt.set_ylim(0, None); axt.set_yticks([]); axt.set_ylabel('Density')
axt.set_xticks([40, 60, 80, 100, 120, 140]); axt.set_xlabel('RFind score (right end)')
axt.spines[['top', 'right']].set_visible(False)
axt.text(0.98, 0.93, f'{n_above:,} KO × KO pairs above the highest\nsafe × safe pair (dashed; expected {n_exp:.0f})\nticks: the {int((ss >= 30).sum())} safe × safe pairs above 30',
         transform=axt.transAxes, fontsize=FN, ha='right', va='top')
bh = ax.get_ylim()[1] * 0.06
ax.add_patch(Rectangle((30, 0), 17, bh, fill=False, ec='black', lw=0.5, clip_on=False, zorder=5))
for xa, xb in ((30, 0), (47, 1)):
    fig.add_artist(ConnectionPatch(xyA=(xa, 0), coordsA=ax.transData, xyB=(xb, 1), coordsB=axt.transAxes, color='#4D4D4D', lw=0.4, ls=':'))
save(fig, 'FigS3c_delta_distribution')

## ---- d: the 20 most similar knock-out pairs ----
## after each bar: the complex both knock-outs belong to, when both are subunits of one known complex or pathway;
## the title says that these are the top KO x KO pairs; red = 'Reported by Shi et al.'
COMPLEX = {('Psmb4', 'Psmc5'): 'proteasome', ('Psmb4', 'Psmc1'): 'proteasome', ('Psmc1', 'Psmc5'): 'proteasome',
           ('Pomp', 'Psmb4'): 'proteasome assembly', ('Pomp', 'Psmc1'): 'proteasome assembly', ('Pomp', 'Psmc5'): 'proteasome assembly',
           ('Uba5', 'Ufc1'): 'UFMylation (E1–E2)', ('Ap2m1', 'Ap2s1'): 'AP-2 adaptor', ('Tsc1', 'Tsc2'): 'TSC complex',
           ('Smarcb1', 'Smarcc2'): 'SWI/SNF', ('Ddx23', 'Prpf6'): 'U5 snRNP', ('Atxn1', 'Cic'): 'ATXN1–CIC',
           ('Elp1', 'Elp2'): 'Elongator', ('Hspa5', 'Hspa9'): 'HSP70 family'}
FW, FH = SIZE['d']
fig = mm_fig(FW, FH); ax = mm_axes(fig, 21.5, 7.6, FW - 23.0, FH - 7.6 - 4.4)
labs = [f'{a} × {b}' for a, b, _ in TOP]
vals = [v for _, _, v in TOP]
cols = [RED if shi_of(a, b) else KOC for a, b, _ in TOP]
y = np.arange(len(TOP))
ax.barh(y, vals, color=cols, height=0.72)
for yi, (a, b, v) in zip(y, TOP):
    lab = COMPLEX.get(tuple(sorted((a, b))), '')
    if lab:
        ax.text(v + 8, yi, lab, va='center', ha='left', fontsize=FN, color='#4D4D4D')
ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=FT); ax.invert_yaxis()
for t_ in ax.get_yticklabels():
    t_.set_fontstyle('italic')
ax.set_xlabel('RFind score'); ax.set_xlim(0, max(vals) * 1.42)
ax.set_title('20 most similar knock-out pairs', fontsize=FL, pad=2.5)
ax.spines[['top', 'right']].set_visible(False)
ax.legend(handles=[Patch(facecolor=RED, label='Reported by Shi et al.'), Patch(facecolor=KOC, label='Other')],
          loc='lower right', frameon=False, fontsize=FL, handlelength=1.6, labelspacing=0.3)
save(fig, 'FigS3d_top_pairs')
open(os.path.join(OUT, f'figS3_rfind_pairs_{CTX}.txt'), 'w').write('\n'.join(lines) + '\n')
