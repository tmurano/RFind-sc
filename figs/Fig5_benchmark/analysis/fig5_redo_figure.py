"""Fig. 5b — the DAM example of the benchmark (panels c and d of this script), drawn at the printed size.
Its two outputs, Fig5c_embedding.png and Fig5d_auc.png, are panel b of Fig. 5 (see "Use in the paper" below).
The figure is 180 mm wide as the slides of Figs. 1-4 and uses their fonts: Arial, text 6-7 pt.
  c  the DAM task on the Hammond t-SNE of Fig. 2: the cells in / not in the state, and their DAM-index (RFind-sc; 0 = white);
  d  the distributions of the DAM-index with the ROC curve as an inset (the AUC strip of d, drawn right of them, is not saved).
No panel title in d: what the panel shows is said in the legend.
draw_e, draw_f and draw_key (one small panel per task, the mean AUC by the length of the reference, their legend) are not used.
Input : figs/Fig5_benchmark/output/results/fig5_redo_{auc,percell}_L<L>.*
Output: figs/Fig5_benchmark/output/data/Fig5c_embedding.png and Fig5d_auc.png
Run   : python figs/Fig5_benchmark/analysis/fig5_redo_figure.py [--test] [--palette=tab|type]   (from the project root; any Python
        with numpy, pandas, matplotlib; figs/run_all_figs.sh runs it with $PYTHON)
Use in the paper: Fig. 5 of the paper (whole figure and panels) is drawn by fig5_rebuild_figure_draft.py; this script
        is run before it because panel b of that figure is cut from Fig5c_embedding.png and Fig5d_auc.png made here
        (Fig5d_auc.png is saved without the AUC strip of d).
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.transforms import Bbox
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

TEST = '--test' in sys.argv
TAG = next((a.split('=', 1)[1] for a in sys.argv if a.startswith('--tag=')), '')
PALETTE = next((a.split('=', 1)[1] for a in sys.argv if a.startswith('--palette=')), 'tab')
O = os.path.join(os.environ.get('FIG5_TEST_DIR', 'logs/runs/fig5_test'), '') if TEST else 'figs/Fig5_benchmark/output/'   # test: as fig5_bench_redo.R
OUT = O                                                            # panels at their pasted size (as in all figures)
os.makedirs(OUT, exist_ok=True)
FW, DPI = 180.0, 600                                              # width of the figure (mm; the slides of Figs. 1-4) and resolution
T6, T65, T7 = 6.0, 7.0, 7.0                                       # text (pt): ticks, axis labels and keys, titles (all figures: 6 / 7 pt)
plt.rcParams.update({'font.family': 'Arial', 'font.size': T65, 'axes.labelsize': T65, 'xtick.labelsize': T6, 'ytick.labelsize': T6,
                     'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2, 'ytick.major.size': 2,
                     'xtick.major.pad': 1.2, 'ytick.major.pad': 1.2, 'axes.labelpad': 1.5, 'axes.spines.top': False,
                     'axes.spines.right': False, 'pdf.fonttype': 42})

MEMB = ['AddModuleScore', 'UCell (up genes)', 'AUCell']           # gene-set membership: on the signed tasks = up genes only, outside their intended use
SETS = ['singscore', 'ssGSEA', 'VISION', 'Beyondcell']            # up and down sets
RANK = ['scFOCAL-type']                                           # ranked reference needed
FIRST = MEMB + SETS + RANK
ORDER = ['RFind-sc'] + FIRST                                      # one order of the methods: d, e and the key
LABEL = {'UCell (up genes)': 'UCell'}
RED, BLUE = '#C0392B', '#1F5FA8'                                  # the red and blue of Figs. 1 and 4
PAL = {'tab': {'AddModuleScore': '#1F77B4', 'UCell (up genes)': '#17BECF', 'AUCell': '#9467BD', 'singscore': '#2CA02C',   # one saturated colour per method
               'ssGSEA': '#BCBD22', 'VISION': '#8C564B', 'Beyondcell': '#E377C2', 'scFOCAL-type': '#FF7F0E'},
       'type': {'AddModuleScore': '#003F6B', 'UCell (up genes)': '#0072B2', 'AUCell': '#56B4E9', 'singscore': '#00563F',     # shades of one colour per type
                'ssGSEA': '#009E73', 'VISION': '#5CC1A0', 'Beyondcell': '#A6DBA0', 'scFOCAL-type': '#B36B00'}}[PALETTE]
COL = {'RFind-sc': RED, **PAL}
MARK = {'RFind-sc': 'o', **{m: 'o' for m in MEMB}, **{m: 's' for m in SETS}, **{m: 'D' for m in RANK}}   # marker = type of method
GENESET = ['H2', 'H3', 'H5']                                      # gene sets (type B in the paper; labelled type A in this layout)
SIGNED = ['H1', 'H4', 'D1', 'D2', 'D3']                           # signed lists (type A in the paper; labelled type B in this layout)
NUM = {'H2': 1, 'H3': 2, 'H5': 3, 'H1': 4, 'H4': 5, 'D1': 6, 'D2': 7, 'D3': 8}   # the task numbers (e and f)
TNAME = {'H2': ('Interferon', 'cluster 9'), 'H3': ('Interferon', 'P540 vs P100'), 'H5': ('G2M', 'cluster 2b'), 'H1': ('DAM', 'LPC vs saline'),
         'H4': ('Maturity', 'P100 vs E14–P5'), 'D1': ('Glucocorticoid', 'Oligodendrocytes'), 'D2': ('Glucocorticoid', 'Astro'),
         'D3': ('Glucocorticoid', 'Endothelial')}
EX = 'H1'                                                         # the example of c and d: DAM
LS = ['100', '250', '500', '1000', '2000', 'all']


def f(stem, L, ext='csv'):
    return f'{O}results/{stem}_L{L}{TAG}.{ext}'                  # the benchmark tables in output/results/


def lab(m):
    return LABEL.get(m, m)


def axmm(fig, H, x, y, w, h):                                     # axes at (x, y) mm from the top left of a figure H mm high
    return fig.add_axes([x / FW, 1 - (y + h) / H, w / FW, h / H])


def txmm(fig, H, x, y, s, **kw):                                  # text at (x, y) mm from the top left
    return fig.text(x / FW, 1 - y / H, s, **kw)


auc = pd.read_csv(f('fig5_redo_auc', '500'))
auc = auc[auc.id.isin(GENESET + SIGNED)]
meths = [m for m in ORDER if m in set(auc.method)]
pc = pd.read_csv(f('fig5_redo_percell', '500', 'csv.gz'))
ALL = pd.concat([pd.read_csv(f('fig5_redo_auc', L), dtype={'L': str}) for L in LS if os.path.exists(f('fig5_redo_auc', L))], ignore_index=True)
Lh = [L for L in LS if L in set(ALL.L)]


def dot(ax, x, y, m, off=False, zb=3):                            # one method: RFind-sc large with a black edge; open = outside the intended use
    big = m == 'RFind-sc'
    ax.scatter([x], [y], s=20 if big else 9, marker=MARK[m], facecolors='white' if off else COL[m], edgecolors='black' if big else COL[m],
               linewidths=0.5 if big else 0.8, zorder=zb + 2 if big else zb, clip_on=False)


## ---- c: the DAM task on the Hammond t-SNE (the same as Fig. 2; said in the legend) ----
def draw_c(fig, H, x0, y0):
    xy = pd.read_csv('data/Hammond2019/raw/cell_meta_hammond.csv', usecols=['cell_id', 'tsne_x', 'tsne_y']).rename(columns={'tsne_x': 'x', 'tsne_y': 'y'})
    e = pc[pc.id == EX][['cell_id', 'positive', 'RFind-sc']].merge(xy, on='cell_id', how='left')
    assert e[['x', 'y']].notna().all().all(), EX
    v = e['RFind-sc'].to_numpy(dtype=float); pos = e['positive'].to_numpy(dtype=bool)
    axs = [axmm(fig, H, x0 + 1 + k * 35, y0 + 6, 34, 34) for k in range(2)]
    for ax in axs:
        ax.scatter(xy.x, xy.y, s=0.03, c='#E3E3E3', linewidths=0, rasterized=True, zorder=1)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect('equal', adjustable='datalim')
        for sp in ax.spines.values():
            sp.set_visible(False)
    o = np.random.default_rng(1).permutation(len(e))                  # random order: no group hides the other
    axs[0].scatter(e.x.values[o], e.y.values[o], s=0.3, c=np.where(pos[o], RED, BLUE), alpha=0.8, linewidths=0, rasterized=True, zorder=2)
    axs[0].scatter([], [], s=6, c=RED, label='LPC (positive)'); axs[0].scatter([], [], s=6, c=BLUE, label='saline (negative)')
    axs[0].legend(loc='upper left', bbox_to_anchor=(0.12, 1.12), frameon=False, fontsize=T7, handletextpad=0.1, borderaxespad=0, labelspacing=0.3)
    DIV = LinearSegmentedColormap.from_list('div', [BLUE, '#F2F2F2', RED])   # the score itself, 0 = white (the values are shown as they are)
    vl = np.percentile(np.abs(v), 99); nrm = Normalize(vmin=-vl, vmax=vl); o2 = np.argsort(np.abs(v), kind='mergesort')   # symmetric about 0
    im = axs[1].scatter(e.x.values[o2], e.y.values[o2], s=0.3, c=np.clip(v[o2], -vl, vl), cmap=DIV, norm=nrm, linewidths=0, rasterized=True, zorder=2)
    cb = fig.colorbar(im, cax=axmm(fig, H, x0 + 72, y0 + 11, 1.6, 24), extend='both'); cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=T6, width=0.5, length=1.5, pad=1); cb.set_label('DAM-index (RFind-sc)', fontsize=T6, labelpad=2)


## ---- d: score distributions with the ROC curve as an inset; the AUC of each method (95% CI) ----
def draw_d(fig, H, x0, y0):
    d = pc[pc.id == EX]; v = d['RFind-sc'].to_numpy(dtype=float); pos = d['positive'].to_numpy(dtype=bool)
    r = auc[auc.id == EX].set_index('method')
    h1 = axmm(fig, H, x0 + 6, y0 + 6, 44, 30)
    lo, hi = np.nanpercentile(v, [0.5, 99.5]); bins = np.linspace(lo, hi, 40)
    n1, _, _ = h1.hist(np.clip(v[~pos], lo, hi), bins=bins, color=BLUE, alpha=0.55, density=True, label='saline (negative)')
    n2, _, _ = h1.hist(np.clip(v[pos], lo, hi), bins=bins, color=RED, alpha=0.55, density=True, label='LPC (positive)')
    h1.set_ylim(0, 1.4 * max(n1.max(), n2.max()))                  # room above the peaks for the key and the ROC curve
    h1.axvline(0, ymax=0.7, color='k', lw=0.5, ls=':'); h1.set_yticks([]); h1.set_xlabel('DAM-index (RFind-sc)'); h1.set_ylabel('Density')
    h1.legend(loc='upper left', frameon=False, fontsize=T7, handlelength=0.9, handletextpad=0.4, borderaxespad=0.1, labelspacing=0.3)
    o = np.argsort(-v, kind='mergesort'); pp = pos[o]
    tpr = np.concatenate([[0], np.cumsum(pp) / pp.sum()]); fpr = np.concatenate([[0], np.cumsum(~pp) / (~pp).sum()])
    ins = inset_axes(h1, width='30%', height='44%', loc='upper right', borderpad=0.15)
    ins.plot(fpr, tpr, color=RED, lw=0.8); ins.fill_between(fpr, tpr, color=RED, alpha=0.15, lw=0); ins.plot([0, 1], [0, 1], color='k', lw=0.4, ls=':')
    ins.set_xlim(0, 1); ins.set_ylim(0, 1); ins.set_xticks([0, 1]); ins.set_yticks([0, 1]); ins.tick_params(labelsize=5.5, length=1.5, pad=0.8)
    ins.set_xlabel('FPR', fontsize=5.5, labelpad=-3); ins.set_ylabel('TPR', fontsize=5.5, labelpad=-3)
    ins.text(0.96, 0.08, f'AUC {r.loc["RFind-sc", "AUC"]:.3f}', ha='right', fontsize=5.5, color=RED); ins.set_title('ROC', fontsize=T6, pad=1.5)
    h2 = axmm(fig, H, x0 + 62, y0 + 6, 31, 30)
    ms = [m for m in ORDER if m in r.index]
    h2.axhline(r.loc['RFind-sc', 'AUC'], color=RED, lw=0.6, ls=(0, (3, 2)), zorder=1)
    for xi, m in enumerate(ms):
        h2.plot([xi, xi], [r.loc[m, 'lo'], r.loc[m, 'hi']], color=COL[m], lw=0.8, zorder=2)
        dot(h2, xi, r.loc[m, 'AUC'], m, off=m in MEMB and EX in SIGNED)
    h2.set_xticks(range(len(ms))); h2.set_xticklabels([lab(m) for m in ms], rotation=55, ha='right', rotation_mode='anchor')
    for tl, m in zip(h2.get_xticklabels(), ms):
        if m == 'RFind-sc':
            tl.set_color(RED); tl.set_fontweight('bold')
    h2.set_xlim(-0.7, len(ms) - 0.3); h2.set_ylabel('AUC (95% CI)')


## ---- e: one small panel per task; AUC axis fitted to each task; RFind-sc = large red point + dashed line ----
PW, PH, PSTEP = 19.0, 17.0, 24.5                                  # plot area of a small panel and the step between panels (mm)
SPAN = 0.2                                                        # the AUC range of every small panel (the same scale in all)


def draw_e(fig, H, x0, y0, A, key=True):
    rows = [('Reference type A: gene set (which genes only)', GENESET, ['RFind-sc'] + MEMB + SETS),
            ('Reference type B: signed list (direction and size of change)', SIGNED, ORDER)]
    for rw, (head, tasks, ms) in enumerate(rows):
        yr = y0 + rw * 31
        txmm(fig, H, x0 + 6, yr + 3.4, head, fontsize=T65, fontweight='bold', va='baseline')
        mm = [m for m in ms if m in meths]
        for k, t in enumerate(tasks):
            xa = x0 + 9 + k * PSTEP; ax = axmm(fig, H, xa, yr + 11, PW, PH)
            txmm(fig, H, xa + PW / 2, yr + 6.8, f'{NUM[t]}  {TNAME[t][0]}', fontsize=T6, fontweight='bold', ha='center', va='baseline')
            txmm(fig, H, xa + PW / 2, yr + 9.3, TNAME[t][1], fontsize=T6, ha='center', va='baseline', color='0.3')
            v = A[A.id == t].set_index('method'); ci = {'lo', 'hi'} <= set(v.columns)
            ok = [m for m in mm if m in v.index and np.isfinite(v.loc[m, 'AUC'])]
            yhi = min(1.012, max(v.loc[m, 'hi'] if ci and np.isfinite(v.loc[m, 'hi']) else v.loc[m, 'AUC'] for m in ok) + 0.012); ylo = yhi - SPAN
            ax.axhline(v.loc['RFind-sc', 'AUC'], color=RED, lw=0.6, ls=(0, (3, 2)), zorder=1)
            for m in ok:
                i = mm.index(m); y = v.loc[m, 'AUC']; off = m in MEMB and t in SIGNED
                if ci and np.isfinite(v.loc[m, 'lo']) and v.loc[m, 'hi'] > ylo:
                    ax.plot([i, i], [max(v.loc[m, 'lo'], ylo), v.loc[m, 'hi']], color=COL[m], lw=0.7, zorder=2, solid_capstyle='butt')
                if y < ylo:                                                   # below the range: a triangle at the bottom edge
                    ax.scatter([i], [ylo], marker='v', s=10, facecolors='white' if off else COL[m], edgecolors=COL[m], linewidths=0.8, clip_on=False, zorder=3)
                else:
                    dot(ax, i, y, m, off=off)
            ax.set_ylim(ylo, yhi); ax.set_xlim(-0.8, len(mm) - 0.2); ax.set_xticks([])
            ax.yaxis.set_major_locator(MultipleLocator(0.05)); ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            if k == 0:
                ax.set_ylabel('AUC')
    if key:
        draw_key(fig, H, x0 + 9 + 3 * PSTEP - 2, y0 + 5)   # 2 mm to the left, so that the 7-pt key clears the y label of f


def style(m):                                                     # the line and marker of a method (f and the legend)
    off = m in MEMB; rf = m == 'RFind-sc'
    return dict(marker=MARK[m], ms=3.8 if rf else 2.6, lw=1.4 if rf else 0.8, color=COL[m], ls=(0, (1.2, 1.2)) if off else '-',
                markerfacecolor='white' if off else COL[m], markeredgecolor='black' if rf else COL[m], markeredgewidth=0.5 if rf else 0.7)


def draw_key(fig, H, x, y):                                       # one legend for e and f
    ax = axmm(fig, H, x, y, 44, 26); ax.axis('off')
    lg = ax.legend([Line2D([], [], **style(m)) for m in meths], [lab(m) for m in meths], ncol=2, loc='upper left', frameon=False, fontsize=T7,
                   handlelength=1.6, columnspacing=0.5, handletextpad=0.3, labelspacing=0.35, borderaxespad=0, title='Methods',
                   title_fontproperties={'weight': 'bold', 'size': T65})
    lg._legend_box.align = 'left'


## ---- f: mean AUC over the tasks 4-8 by the length of the reference (the key is shared with e) ----
def draw_f(fig, H, x0, y0):
    ax = axmm(fig, H, x0 + 8, y0 + 7, 40, 43)
    for m in meths:
        sel = [(ALL.method == m) & (ALL.L == L) & ALL.id.isin(SIGNED) for L in Lh]
        y = [ALL[s].AUC.mean() if s.sum() == len(SIGNED) else np.nan for s in sel]
        ax.plot(range(len(Lh)), y, zorder=5 if m == 'RFind-sc' else 3, **style(m))
    ax.set_xticks(range(len(Lh))); ax.set_xticklabels([{'1000': '1,000', '2000': '2,000'}.get(L, L) for L in Lh])   # x label in two lines: at 7 pt it does not fit in 40 mm
    ax.axvline(Lh.index('500'), color='0.5', lw=0.4, ls=':', zorder=0)
    ax.set_xlabel('Reference length\n(genes per direction)'); ax.set_ylabel(f'Mean AUC (tasks {NUM[SIGNED[0]]}–{NUM[SIGNED[-1]]})')


## ---- the two panels used by the paper: c and d (Fig. 5 of the paper = fig5_rebuild_figure_draft.py) ----
DATA = O + 'data/'; os.makedirs(DATA, exist_ok=True)               # intermediate images read by fig5_rebuild_figure_draft.py
XD, YC = 86.0, 0.5                                                 # left of d; top of c and d (mm)
H = YC + 53.0 + 0.5
fig = plt.figure(figsize=(FW / 25.4, H / 25.4))
draw_c(fig, H, 0, YC); draw_d(fig, H, XD, YC)
for name, (x, y, w, h) in {'Fig5c_embedding': (0, YC, XD, 53), 'Fig5d_auc': (XD, YC, 52, 53)}.items():   # d: without its AUC strip
    fig.savefig(f'{DATA}{name}.png', dpi=DPI, facecolor='white', bbox_inches=Bbox.from_bounds(x / 25.4, (H - y - h) / 25.4, w / 25.4, h / 25.4))
plt.close(fig)
print('Fig5c_embedding.png and Fig5d_auc.png in', DATA)
