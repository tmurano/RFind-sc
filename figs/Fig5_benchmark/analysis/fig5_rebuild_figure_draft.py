"""Fig. 5 — the capability benchmark, the whole figure of the paper (output/Fig5.png) and its panels as separate images for the slide
(output/, next to it), drawn by the same functions.
Style of fig5_redo_figure.py: 180 mm wide, Arial 6-7 pt, one colour per method, marker = the form of reference the method is built for,
RFind-sc = the large red point and the red dashed line, small panels with an AUC range of 0.2, a triangle at the bottom edge for a value
below the range, 95% interval as a bar.
  a  what each method can use from a reference (table, with the form of the reference each method is given in the type A and the
     type B tasks; the key to the symbols is in the legend, not under the table)
  b  the DAM example (t-SNE, the DAM-index distributions and the ROC curve of RFind-sc, and the AUC of every method in this task)
  c  type A (2 tasks: DAM and Glucocorticoid in Astro; the Glucocorticoid tasks in Oligodendrocytes and Endothelial are in
     Supplementary Fig. S4b): every method once, in the most complete form it can take, side by side; the method key right of them
  d  type B (3 tasks: MSigDB Hallmark sets): every method that takes a gene set
  e  the mean over the tasks of each type: the mean AUC of RFind-sc minus that of each method, with its 95% interval, over the eight
     type A tasks (500 genes per direction) and the three type B tasks (results/fig5_task_means.csv; fig5_task_means.py)
  f  reference length (type A): mean AUC over the eight type A tasks (the four DAM and Glucocorticoid tasks and the four PS2APP and
     Forced swim tasks: DA1, the cells of the DAM task with the PS2APP reference; AS1-AS3, the cells of the three Glucocorticoid tasks
     with the Forced swim reference; fig5_extra_score.R) by the number of reference genes per direction
All symbols are filled: the figure shows no comparison with RFind-sc (the comparisons are in Supplementary Table 3). The panels have
no header line (the legend names the type of reference of c and d). Layout: c (with the key) above d on the left, three small panels
wide; e above f on the right; the whole figure fits the 180 mm height of the slide.
Intervals: bootstrap over the mice of the task (all cells of each drawn mouse; 1,000 resamples, the same draws for every method), as in
the benchmark.
The panels are blocks placed in rows by render(). For the slide, render_panels() draws each block alone (no panel letter: the letters
belong to the slide) on a canvas of the figure's width and cuts it to its drawn content with a 0.5 mm margin, so that the images keep
the physical size of the figure at 600 dpi (paste at 100%).
Inputs: figs/Fig5_benchmark/output/results/fig5_redo_auc_L<L>.csv (all tasks of the figure), extra_refs/fig5_extra_auc_L<L>.csv (the added
        type A tasks, f only), results/fig5_task_means.csv (e), Fig5c_embedding.png and Fig5d_auc.png (b; fig5_redo_figure.py)
Output: figs/Fig5_benchmark/output/Fig5.png (the figure of the paper) and, in the same folder, the
        panels Fig5a_methods_table.png, Fig5b_DAM_example.png, Fig5c_typeA.png, Fig5d_typeB.png, Fig5e_task_means.png,
        Fig5f_reference_length.png and Fig5_methods_key.png (the method key; the note on filled / open symbols is in the legend, KEY_NOTE): the panels for the
        slide, with their sizes in mm printed
Run from the project root: python3 figs/Fig5_benchmark/analysis/fig5_rebuild_figure_draft.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.lines import Line2D
from matplotlib.font_manager import findfont, get_font
from matplotlib.ft2font import Kerning, LoadFlags

O = 'figs/Fig5_benchmark/output/'; OUT = O; os.makedirs(OUT, exist_ok=True)   # the figure and its panels straight into output/
## AUC tables: the benchmark run of fig5_bench_redo.R.
S = O
FW, DPI = 180.0, 600
T6, T7 = 6.0, 7.0
plt.rcParams.update({'font.family': 'Arial', 'font.size': T7, 'axes.labelsize': T7, 'xtick.labelsize': T6, 'ytick.labelsize': T6,
                     'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2, 'ytick.major.size': 2,
                     'xtick.major.pad': 1.2, 'ytick.major.pad': 1.2, 'axes.labelpad': 1.5, 'axes.spines.top': False,
                     'axes.spines.right': False, 'pdf.fonttype': 42})
RED, NAVY = '#C0392B', '#34495E'
RF = 'RFind-sc'
## colour = method (tab10 without red; the decoupleR methods get darker tab20b colours), marker = the form the method is built for
COL = {RF: RED, 'AddModuleScore': '#1F77B4', 'UCell (up genes)': '#17BECF', 'UCell': '#17BECF', 'AUCell': '#9467BD',
       'ORA': '#7B4173', 'singscore': '#2CA02C', 'ssGSEA': '#BCBD22', 'VISION': '#8C564B', 'Beyondcell': '#E377C2', 'ULM': '#393B79',
       'WMEAN': '#637939', 'scFOCAL-type': '#FF7F0E', 'ULM (log2FC)': '#393B79', 'WMEAN (log2FC)': '#637939'}
MARK = {RF: 'D', 'AddModuleScore': 'o', 'UCell (up genes)': 'o', 'AUCell': 'o', 'ORA': 'o', 'UCell': 's', 'singscore': 's',
        'ssGSEA': 's', 'VISION': 's', 'Beyondcell': 's', 'ULM': 'D', 'WMEAN': 's', 'scFOCAL-type': 'D', 'ULM (log2FC)': 'D', 'WMEAN (log2FC)': 'D'}
LABEL = {'UCell (up genes)': 'UCell (up genes)', 'UCell': 'UCell', 'ULM': 'ULM (decoupleR)', 'ULM (log2FC)': 'ULM (decoupleR)'}
## c: every method once, in the most complete form it can take (RFind-sc, scFOCAL-type and decoupleR with log2FC weights: the ranked
## signed list; the up/down-set methods, UCell included (its signed signature, up and down genes): up and down names; the membership
## methods: up names only). UCell with the up genes only is reported in Supplementary Table 3.
## The forms are given in the table (a) and the text; c shows the methods side by side in the order of the table, as d does.
## Reported in Supplementary Table 3 only: decoupleR's WMEAN and ORA (decoupleR is represented by one method, ULM). Beyondcell is shown
## although it gives no score to the cells expressing fewer than 10% of the reference genes (counted as 0, as in the benchmark; the
## numbers of these cells are in Supplementary Table 3, sheet Beyondcell missing).
CMETH = [RF, 'AddModuleScore', 'AUCell', 'UCell', 'singscore', 'ssGSEA', 'VISION', 'Beyondcell', 'ULM (log2FC)', 'scFOCAL-type']
SETM = [RF, 'AddModuleScore', 'AUCell', 'UCell', 'singscore', 'ssGSEA', 'VISION', 'Beyondcell', 'ULM']
TA = ['H1', 'D1', 'D2', 'D3']                                      # the DAM and Glucocorticoid tasks (averaged in e and f with EXTRA)
TA_C = ['H1', 'D2']                                                # the type A tasks drawn in c (and in b: H1); D1 and D3 in Supplementary Fig. S4b
EXTRA = ['DA1', 'AS1', 'AS2', 'AS3']                               # the added type A tasks (Supplementary Fig. S4)
PAIR = {'DA1': 'H1', 'AS1': 'D1', 'AS2': 'D2', 'AS3': 'D3'}        # added task -> the DAM or Glucocorticoid task on the same cells
TA_EVAL = TA + EXTRA                                               # type A tasks averaged in e (the eight)
TB = ['H2', 'H3', 'H5']                                            # Maturity (H4) is left out of the main figure: at the ceiling
                                                                      # (median AUC of the methods >= 0.99 at L = 500); shown only in Supplementary Table 3
TNAME = {'H2': ('Interferon', 'cluster 9'), 'H3': ('Interferon', 'P540 vs P100'), 'H5': ('G2M', 'cluster 2b'), 'H1': ('DAM', 'LPC vs saline'),
         'H4': ('Maturity', 'P100 vs E14–P5'), 'D1': ('Glucocorticoid', 'Oligodendrocytes'), 'D2': ('Glucocorticoid', 'Astro'),
         'D3': ('Glucocorticoid', 'Endothelial'), 'DA1': ('PS2APP', 'LPC vs saline'), 'AS1': ('Forced swim', 'Oligodendrocytes'),
         'AS2': ('Forced swim', 'Astro'), 'AS3': ('Forced swim', 'Endothelial')}
LS = ['100', '250', '500', '1000', '2000', 'all']


def lab(m):
    return LABEL.get(m, m)


def axmm(fig, H, x, y, w, h):
    return fig.add_axes([x / FW, 1 - (y + h) / H, w / FW, h / H])


def txmm(fig, H, x, y, s, **kw):
    return fig.text(x / FW, 1 - y / H, s, **kw)


## G2M with a subscript 2: mathtext draws a subscript at 0.7 x the label (4.2 pt for 6 pt), below Nature's minimum of 5 pt for
## any text, and its size cannot be set. The label is therefore drawn as two texts in the same Arial (bold in a bold label): G and M at
## the size of the label, as one text with an en space (0.5 em) between them (one baseline for both), and the 2 at 5 pt in that space,
## 0.06 pt after G (as mathtext placed it) and its baseline lowered by 0.98 pt (as mathtext). Display only: TNAME (also read by
## figS4_capability_map.py) keeps 'G2M'.
SUB_PT, SUB_DROP, SUB_KERN, EN = 5.0, 0.98, 0.06, '\u2002'        # size of the 2, lowering of its baseline, gap after G (pt), en space
PT = 25.4 / 72                                                     # mm per point


def ink_pt(s, prop):                                               # advance, left and right ink edges of s in the font prop (pt; unhinted)
    f = get_font(findfont(prop)); f.set_size(prop.get_size_in_points(), 72)
    pen, left, right, prev = 0.0, None, None, None
    for c in s:
        gi = f.get_char_index(ord(c))
        pen += f.get_kerning(prev, gi, Kerning.UNFITTED) / 64 if prev is not None else 0.0
        g = f.load_glyph(gi, flags=LoadFlags.NO_HINTING)
        if g.bbox[2] > g.bbox[0]:                                  # a glyph with ink (not a space)
            left = pen + g.bbox[0] / 64 if left is None else left
            right = pen + g.bbox[2] / 64
        pen, prev = pen + g.linearHoriAdvance / 65536, gi
    return pen, left, right


def txmm_g2m(fig, H, x, y, s, fontsize, **kw):                     # a label centred at x on the baseline y (mm); G2M with a subscript 2
    if 'G2M' not in s:
        return txmm(fig, H, x, y, s, fontsize=fontsize, ha='center', va='baseline', **kw)
    i = s.index('G2M'); pre = s[:i] + 'G'
    t = txmm(fig, H, x, y, pre + EN + 'M' + s[i + 3:], fontsize=fontsize, ha='left', va='baseline', **kw)
    t2 = txmm(fig, H, x, y + SUB_DROP * PT, '2', fontsize=SUB_PT, ha='left', va='baseline', **kw)
    p, p2 = t.get_fontproperties(), t2.get_fontproperties()
    assert findfont(p2) == findfont(p) and SUB_KERN + ink_pt('2', p2)[0] < ink_pt(EN, p)[0], s   # same font; the 2 fits in the space
    _, l, r = ink_pt(t.get_text(), p)
    x0 = x - (l + r) / 2 * PT                                      # pen of the label (mm): its ink centred at x
    t.set_x(x0 / FW); t2.set_x((x0 + (ink_pt(pre, p)[0] + SUB_KERN) * PT) / FW)
    return t, t2


## ---- data ----
A = pd.read_csv(S + 'results/fig5_redo_auc_L500.csv')[['id', 'method', 'AUC', 'lo', 'hi', 'diff_vs_RFind', 'd_lo', 'd_hi']]   # 95% intervals: bootstrap over the mice
CI = 'bootstrap over mice (1,000 resamples)'
ALL = []                                                           # every length: the tasks of the figure, then the added type A tasks (e)
for L in LS:
    if not os.path.exists(f'{S}results/fig5_redo_auc_L{L}.csv'):
        continue
    m0 = pd.read_csv(f'{S}results/fig5_redo_auc_L{L}.csv', dtype={'L': str}).assign(L=L)
    x0 = pd.read_csv(f'{S}extra_refs/fig5_extra_auc_L{L}.csv', dtype={'L': str}).assign(L=L)
    for t in EXTRA:                                                # each added task once with every method of c, the cells of its pair
        n = x0[(x0.id == t) & x0.method.isin(CMETH)].method.value_counts()
        assert len(n) == len(CMETH) and n.max() == 1, (t, L)
        nn = x0[x0.id == t][['n_pos', 'n_neg']].drop_duplicates()
        assert nn.shape[0] == 1 and (nn.to_numpy() == m0[m0.id == PAIR[t]][['n_pos', 'n_neg']].drop_duplicates().to_numpy()).all(), (t, L)
    ALL += [m0, x0[x0.id.isin(EXTRA)]]
ALL = pd.concat(ALL, ignore_index=True)
TM = pd.read_csv(S + 'results/fig5_task_means.csv')              # e: the mean over the tasks of each type (fig5_task_means.py)
HAVE = set(A.method)
assert {'ULM (log2FC)'} <= HAVE                                   # decoupleR of c and e: ULM with the log2FC weights
print('intervals:', CI, '| methods with results at L = 500:', sorted(HAVE))


def differs(r):                                                   # the comparison rule of the benchmark (one row of results/fig5_redo_auc_L500.csv)
    return bool((r.d_lo > 0 and r.diff_vs_RFind > 0.01) or (r.d_hi < 0 and r.diff_vs_RFind < -0.01))


def dot(ax, x, y, m, zb=3, hollow=False):
    rf = m.startswith('RFind-sc')
    ax.scatter([x], [y], s=(20 if m == RF else 14) if rf else 9, marker=MARK[m], facecolors='white' if hollow else COL[m],
               edgecolors='black' if rf else COL[m], linewidths=0.5 if rf else 0.8, zorder=zb + 2 if rf else zb, clip_on=False)


SPAN = 0.2


def small(ax, v, ms, xs, rfline=True):                             # one small panel
    ok = [(m, x) for m, x in zip(ms, xs) if m in v.index and np.isfinite(v.loc[m, 'AUC'])]
    yhi = min(1.012, max(v.loc[m, 'hi'] if np.isfinite(v.loc[m, 'hi']) else v.loc[m, 'AUC'] for m, _ in ok) + 0.012); ylo = yhi - SPAN
    if rfline and RF in v.index:
        ax.axhline(v.loc[RF, 'AUC'], color=RED, lw=0.6, ls=(0, (3, 2)), zorder=1)
    for m, x in ok:
        y = v.loc[m, 'AUC']; hollow = False                          # no open symbols: no win/lose marks in Fig. 5
        if np.isfinite(v.loc[m, 'lo']) and v.loc[m, 'hi'] > ylo:
            ax.plot([x, x], [max(v.loc[m, 'lo'], ylo), v.loc[m, 'hi']], color=COL[m], lw=0.7, zorder=2, solid_capstyle='butt')
        if y < ylo:
            ax.scatter([x], [ylo], marker='v', s=10, facecolors='white' if hollow else COL[m], edgecolors=COL[m], linewidths=0.8,
                       clip_on=False, zorder=3)
        else:
            dot(ax, x, y, m, hollow=hollow)
    ax.set_ylim(ylo, yhi); ax.set_xticks([])
    ax.yaxis.set_major_locator(MultipleLocator(0.05)); ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))


## ---- a: what each method can use from a reference ----
ROWS = [  # method; as designed: built for, direction, magnitude, gene sets; given in this study: type A tasks, type B tasks; measured
    ('RFind-sc', 'External ranked DE lists; gene sets', '✓', 'Rank', '✓', 'Ranked lists', 'Gene set', 'Best overlap of cell and reference lists'),
    ('Names only (membership)', None),
    ('AddModuleScore', 'Gene sets', '–', '–', '✓', 'Up names only', 'Gene set', 'Set expression above control genes'),
    ('AUCell', 'Gene sets', '–', '–', '✓', 'Up names only', 'Gene set', "Set genes among the cell's top 5%"),
    ('Names and direction (up and down sets)', None),
    ('UCell', 'Gene sets (up/down allowed)', '✓', '–', '✓', 'Up/down names only', 'Gene set', 'Set-gene ranks; up − down'),
    ('singscore', 'Up/down gene sets', '✓', '–', '✓', 'Up/down names only', 'Gene set', 'Mean up / down set ranks'),
    ('ssGSEA', 'Gene sets', '△', '–', '✓', 'Up/down names only', 'Gene set', 'Rank enrichment (up − down)'),
    ('VISION', 'Signatures (±1)', '✓', '–', '✓', 'Up/down names only', 'Gene set', 'Signed mean, standardized'),
    ('Beyondcell', 'Drug signatures (up/down)', '✓', '–', '✓', 'Up/down names only', 'Gene set', 'Up − down expression, scaled'),
    ('Names, direction and magnitude (ranks or weights)', None),
    ('ULM (decoupleR)', 'Weighted networks', '✓', 'Weight', '✓', 'log2FC weights', 'Gene set', 't value of a linear model'),
    ('scFOCAL-type', 'Drug signatures (log2FC)', '✓', 'Rank', '–', 'Ranked lists', 'Not run', 'Rank correlation with log2FC')]
## The header, in two rows, separates what each method uses as designed (built for, direction, magnitude, gene sets) from the
## reference given in this study (the form in the type A tasks, in words; the gene set or 'Not run' in the type B tasks); rank and
## numerical weight in one column, Magnitude; the key to the symbols is in the legend, not under the table
CW = [26.5, 34, 12.5, 12.5, 11, 21, 15.5, 40]                  # column widths (mm), sum 173
A_HH, A_RH = 8.6, 2.95                                             # heights (mm) of the two header rows together and of a row
HEAD = ['Method', 'Built for', 'Direction', 'Magnitude', 'Gene\nsets', 'Type A\ntasks', 'Type B\ntasks', 'What is measured']
CENT = [False, False, True, True, True, True, True, False]
GROUPS = [(1, 4, 'What the method uses, as designed'), (5, 6, 'Reference given in this study')]   # first and last column, label
TOP = {0, 7}                                                       # columns whose header spans both header rows


def draw_a(fig, H, y0):
    x0 = 0.5; g1 = 3.4; hh = A_HH; rh = A_RH                       # g1 = height of the group header row
    ax = axmm(fig, H, x0, y0, sum(CW), hh + rh * len(ROWS) + 1); ax.set_xlim(0, sum(CW)); ax.set_ylim(hh + rh * len(ROWS) + 1, 0); ax.axis('off')
    ax.add_patch(plt.Rectangle((0, 0), sum(CW), hh, color=NAVY, lw=0))
    xs = np.cumsum([0] + CW[:-1])
    for a, b, lab in GROUPS:                                       # group labels, with a white rule under each
        xa, xb = xs[a], xs[b] + CW[b]
        ax.text((xa + xb) / 2, g1 / 2 + 0.2, lab, color='white', fontsize=T6, fontweight='bold', ha='center', va='center')
        ax.plot([xa + 0.8, xb - 0.8], [g1 + 0.1, g1 + 0.1], color='white', lw=0.5)
    for k, (x, w, h, c) in enumerate(zip(xs, CW, HEAD, CENT)):
        yc = hh / 2 if k in TOP else (g1 + hh) / 2 + 0.1
        ax.text(x + (w / 2 if c else 1.2), yc, h, color='white', fontsize=T6, fontweight='bold', ha='center' if c else 'left', va='center', linespacing=0.95)
    y = hh
    for r in ROWS:
        if r[1] is None:
            ax.text(1.2, y + rh * 0.62, r[0], fontsize=T6, color='0.35', va='center'); ax.plot([0, sum(CW)], [y, y], color='0.6', lw=0.4)
            y += rh; continue
        isrf = r[0] == 'RFind-sc'
        for k, (x, w, c) in enumerate(zip(xs, CW, CENT)):
            s = r[k]; sym = s in ('✓', '–', '△')
            ax.text(x + (w / 2 if c else 1.2), y + rh / 2, s, fontsize=T6, ha='center' if c else 'left', va='center',
                    color=RED if isrf else 'black', fontweight='bold' if (k == 0) else 'normal', fontfamily='DejaVu Sans' if sym else 'Arial')
        y += rh
    ax.plot([0, sum(CW)], [y, y], color='0.6', lw=0.4)
    return hh + rh * len(ROWS) + 1.0


## ---- c: type A (2 tasks): every method once, side by side ----
LEFT_W = 104.0                                                     # width of the left column (c above d), about 6/10 of the figure
PW_S, GAP_S = (LEFT_W - 9.3 - (3 - 1) * 5.8) / 3, 5.8              # three small panels fill it (c: two tasks and the key)


def draw_b(fig, H, y0, X1=9.3, GAP=GAP_S):                         # X1 = left edge of the first small panel (room for its axis label)
    PW, PH = PW_S, 24.0                                             # c: the two tasks and the key fill the left column (c and d left, e and f right)
    STEP = PW + GAP
    for k, t in enumerate(TA_C):
        xa = X1 + k * STEP; ax = axmm(fig, H, xa, y0 + 8, PW, PH)
        txmm(fig, H, xa + PW / 2, y0 + 3.6, TNAME[t][0], fontsize=T6, fontweight='bold', ha='center', va='baseline')
        txmm(fig, H, xa + PW / 2, y0 + 6.1, TNAME[t][1], fontsize=T6, ha='center', va='baseline', color='0.3')
        small(ax, A[A.id == t].set_index('method'), CMETH, list(range(len(CMETH))))
        ax.set_xlim(-0.8, len(CMETH) - 0.2)
        if k == 0:
            ax.set_ylabel('AUC (95% CI)')
    return 8 + PH + 2


## ---- d: type B arena (Hallmark sets) ----
def draw_c(fig, H, x0, y0):
    PW, PH, STEP = PW_S, 24.0, PW_S + GAP_S                           # the small panels of d as those of c
    ms = [m for m in SETM]
    for k, t in enumerate(TB):
        xa = x0 + 9.3 + k * STEP; ax = axmm(fig, H, xa, y0 + 8, PW, PH)
        txmm_g2m(fig, H, xa + PW / 2, y0 + 3.6, TNAME[t][0], fontsize=T6, fontweight='bold')   # G2M: subscript 2 at 5 pt
        txmm(fig, H, xa + PW / 2, y0 + 6.1, TNAME[t][1], fontsize=T6, ha='center', va='baseline', color='0.3')
        v = A[A.id == t].set_index('method')
        v = v.rename(index={RF: RF})
        small(ax, v, ms, list(range(len(ms))))
        ax.set_xlim(-0.8, len(ms) - 0.2)
        if k == 0:
            ax.set_ylabel('AUC (95% CI)')


## ---- e: the mean over the tasks of each type, RFind-sc minus each method (the methods of c and d at their places) ----
MEAN_TITLE = {'A': ('Type A', '8 tasks'), 'B': ('Type B', '3 tasks')}


E_PW, E_GAP = 26.5, 6.0                                            # e: width of its two panels and the gap between them (mm)


def draw_means(fig, H, x0, y0, X1=9.3):
    PW, PH, STEP = E_PW, 24.0, E_PW + E_GAP                          # the height of the small panels of c and d
    for k, (kind, ms) in enumerate((('A', CMETH), ('B', SETM))):
        xa = x0 + X1 + k * STEP; ax = axmm(fig, H, xa, y0 + 8, PW, PH)
        txmm(fig, H, xa + PW / 2, y0 + 3.6, MEAN_TITLE[kind][0], fontsize=T6, fontweight='bold', ha='center', va='baseline')
        txmm(fig, H, xa + PW / 2, y0 + 6.1, MEAN_TITLE[kind][1], fontsize=T6, ha='center', va='baseline', color='0.3')
        v = TM[TM.type == kind].set_index('method')
        ax.axhline(0, color='0.45', lw=0.5, ls=(0, (2, 2)), zorder=1)
        top = 0.0
        for x, m in enumerate(ms):
            if m == RF or m not in v.index:
                continue
            d, lo, hi = v.loc[m, ['RFind-sc minus method', 'difference 95% CI low', 'difference 95% CI high']]
            ax.plot([x, x], [lo, hi], color=COL[m], lw=0.7, zorder=2, solid_capstyle='butt')
            dot(ax, x, d, m); top = max(top, hi)
        ax.set_xlim(-0.8, len(ms) - 0.2); ax.set_xticks([])
        step = 0.05 if top > 0.1 else 0.01
        ax.set_ylim(-step / 2, (np.floor(top / step) + 1) * step)
        ax.yaxis.set_major_locator(MultipleLocator(step)); ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        if k == 0:
            ax.set_ylabel('RFind-sc minus method\n(mean AUC; 95% CI)')


## ---- f: reference length (type A) ----
def style(m):
    rf = m == RF
    return dict(marker=MARK[m], ms=3.6 if rf else 2.6, lw=1.4 if rf else 0.8, color=COL[m], ls='-', markerfacecolor=COL[m],
                markeredgecolor='black' if rf else COL[m], markeredgewidth=0.5 if rf else 0.7)


E_MAIN = 28.0                                                      # height (mm) of the plot of f (the whole figure within the slide's 180 mm)
E_W = 32.0                                                         # width (mm) of the plot of f


def draw_d(fig, H, x0, y0):
    ax = axmm(fig, H, x0 + 10, y0 + 2, E_W, E_MAIN)                  # f: right of d (below e), to the right edge of the figure
    Lh = [L for L in LS if L in set(ALL.L)]
    for m in CMETH:
        y = []
        for L in Lh:
            s = ALL[(ALL.method == m) & (ALL.L == L) & ALL.id.isin(TA_EVAL)]
            y.append(s.AUC.mean() if len(s) == len(TA_EVAL) else np.nan)
        if np.all(np.isnan(y)):
            continue
        ax.plot(range(len(Lh)), y, zorder=5 if m == RF else 3, **style(m))
    ax.set_xticks(range(len(Lh))); ax.set_xticklabels([{'1000': '1,000', '2000': '2,000'}.get(L, L) for L in Lh])
    ax.axvline(Lh.index('500'), color='0.5', lw=0.4, ls=':', zorder=0)
    ax.set_xlabel('Reference length\n(genes per direction)'); ax.set_ylabel(f'Mean AUC\n({len(TA_EVAL)} type A tasks)')


MARKS_H = 7.0                                                      # height (mm) of the two lines on filled / open symbols under the key


def key_h(ncol, marks=False):                                      # height (mm) of the key in ncol columns
    return 5.0 + 3.0 * -(-len(CMETH) // ncol) + (MARKS_H if marks else 0.0)


def draw_key(fig, H, x, y, w=175, ncol=5, marks=False):
    ax = axmm(fig, H, x, y, w, key_h(ncol)); ax.axis('off')
    hs = [Line2D([], [], **style(m)) for m in CMETH]
    ls = [lab(m) for m in CMETH]
    lg = ax.legend(hs, ls, ncol=ncol, loc='upper left', frameon=False, fontsize=T6, handlelength=1.6, columnspacing=0.6, handletextpad=0.3,
                   labelspacing=0.3, borderaxespad=0, title='Methods', title_fontproperties={'weight': 'bold', 'size': T7})
    lg._legend_box.align = 'left'
    if marks:                                                      # b-d: filled = different from RFind-sc, open = not distinguished
        y1 = y + key_h(ncol) + 1.8
        a2 = axmm(fig, H, x, y1 - 1.5, 4, 3); a2.axis('off'); a2.set_xlim(0, 4); a2.set_ylim(0, 3)
        a2.scatter([1.0], [1.5], s=9, marker='o', facecolors='0.3', edgecolors='0.3', linewidths=0.8, clip_on=False)
        a2.scatter([3.0], [1.5], s=9, marker='o', facecolors='white', edgecolors='0.3', linewidths=0.8, clip_on=False)
        txmm(fig, H, x + 5.0, y1 + 0.9, 'Filled, AUC higher or lower than that of RFind-sc (95% CI of the', fontsize=T6, va='baseline')
        txmm(fig, H, x + 5.0, y1 + 3.6, 'difference excludes 0 and the difference > 0.01); open, not distinguished', fontsize=T6, va='baseline')


## ---- b: the DAM example (t-SNE and the DAM-index distributions with the ROC curve: the panels of fig5_redo_figure.py as they
## are, Fig5c_embedding.png and Fig5d_auc.png, blank margins trimmed), then the AUC strip of the DAM task drawn here (the methods
## of c, with their names) ----
IE = plt.imread(O + 'data/Fig5c_embedding.png'); IDD = plt.imread(O + 'data/Fig5d_auc.png')   # Fig5d_auc.png is saved without its AUC strip
MM = DPI / 25.4                                                     # pixels per mm


def ink(im):                                                       # pixels with anything drawn (neither white nor transparent)
    k = im[..., :3].min(axis=2) < 0.97
    if im.shape[2] == 4:
        k &= im[..., 3] > 0.05
    return k


rows = np.concatenate([np.flatnonzero(ink(im).any(axis=1)) for im in (IE, IDD)])
r0, r1 = max(rows.min() - int(MM), 0), rows.max() + int(MM) + 1     # 1 mm kept above and below
XAX = int(np.argmax((IDD[..., :3].max(axis=2) < 0.2).sum(axis=1))) # the x axis of the distribution plot (its longest dark row)
cols = [np.flatnonzero(ink(im).any(axis=0)) for im in (IE, IDD)]
IE, IDD = [im[r0:r1, max(c.min() - int(MM / 2), 0):c.max() + int(MM / 2) + 1] for im, c in zip((IE, IDD), cols)]   # 0.5 mm kept left and right
EX_AX = ((XAX - r0) / MM - 30.0, (XAX - r0) / MM)                  # top and bottom of the strip (mm below the top of b): bottom on that x axis
EX_H = max(IE.shape[0] / MM, EX_AX[1] + 13.0)                       # room for the slanted method names


def draw_example(fig, H, y0):
    x = 2.5
    for im in (IE, IDD):
        w, h = im.shape[1] / MM, im.shape[0] / MM; a_ = axmm(fig, H, x, y0, w, h); a_.imshow(im, interpolation='antialiased'); a_.axis('off')
        x += w + 1.0
    xs = x + 8.0; ax = axmm(fig, H, xs, y0 + EX_AX[0], FW - 0.8 - xs, EX_AX[1] - EX_AX[0])
    small(ax, A[A.id == 'H1'].set_index('method'), CMETH, list(range(len(CMETH))))
    ax.set_xlim(-0.7, len(CMETH) - 0.3)
    ax.set_xticks(range(len(CMETH))); ax.set_xticklabels([RF if m == RF else lab(m) for m in CMETH], rotation=45, ha='right', rotation_mode='anchor')
    ax.get_xticklabels()[0].set_color(RED); ax.get_xticklabels()[0].set_fontweight('bold')
    ax.set_ylabel('AUC (95% CI)')


## ---- the whole figure: blocks in rows (a row = blocks side by side, each with its left edge x0 in mm); fixed panel letters (LETTER) ----
TABLE_H = A_HH + A_RH * len(ROWS) + 1.0                           # the height draw_a returns (two header rows, no key below)
BLOCKS = {  # name: (height in mm, drawing function (fig, H, x0, y0 = top of the block), baseline of the panel letter below the top)
    'table': (3.2 + TABLE_H, lambda fig, H, x0, y0: draw_a(fig, H, y0 + 3.2), 2.4),
    'example': (EX_H, lambda fig, H, x0, y0: draw_example(fig, H, y0), 3.0),
    'typeA': (34.0, lambda fig, H, x0, y0: draw_b(fig, H, y0), 3.0),
    'typeB': (34.0, lambda fig, H, x0, y0: draw_c(fig, H, x0, y0), 3.0),
    'means': (34.0, lambda fig, H, x0, y0: draw_means(fig, H, x0, y0), 3.0),
    'length': (2 + E_MAIN + 10.0, lambda fig, H, x0, y0: draw_d(fig, H, x0, y0), 3.0)}
KEY_UNDER_D = (3, 100, 4)                                          # unused alternative: the key under d (x0, width, columns)
KEY_RIGHT_OF_D = (9.3 + 2 * (PW_S + GAP_S) + 0.5, 1)               # the key right of the two panels of c (the third slot): x0 (mm), one column
RIGHT_X = LEFT_W + 7.0                                             # left edge of the right column (e above f), 7 mm right of the left column
KEY_NOTE = False                                                   # the note on filled / open symbols is in the legend, not the figure
CURRENT = [[('table', 0)], [('example', 0)], [('typeA', 0), ('means', RIGHT_X)], [('typeB', 0), ('length', RIGHT_X)]]   # e right of c, f right of d
ROWSPAN = {'length'}                                               # blocks that may be taller than their row (f reaches below d)
LETTER = {'table': 'a', 'example': 'b', 'typeA': 'c', 'typeB': 'd', 'means': 'e', 'length': 'f'}   # as in the text and legend


def render(rows, path, note=None):
    def bh(b):                                                     # height of a block in its row (a key block: the height of the key)
        if b[0] == 'key':
            return key_h(b[3])
        return BLOCKS[b[0]][0]
    hs = [max([bh(b) for b in row if b[0] not in ROWSPAN] or [0.0]) for row in rows]   # f may reach below its row: it does not set its row's height
    bottom = 3.0 if note else 1.0                                  # the draft note needs room below the panels
    H = 0.8 + sum(hs) + 1.5 * (len(rows) - 1) + bottom
    ys = np.cumsum([0.8] + [h + 1.5 for h in hs])
    for row, y in zip(rows, ys):                                   # the figure is at least as tall as f
        for b in row:
            if b[0] in ROWSPAN:
                H = max(H, y + BLOCKS[b[0]][0] + bottom)
    fig = plt.figure(figsize=(FW / 25.4, H / 25.4))
    y, k = 0.8, 0
    for row, h in zip(rows, hs):
        for b in row:
            if b[0] == 'key':
                draw_key(fig, H, b[1], y, w=b[2], ncol=b[3]); continue
            BLOCKS[b[0]][1](fig, H, b[1], y)
            if b[0] == 'typeA':                                    # the key right of the two panels of c, one column
                draw_key(fig, H, b[1] + KEY_RIGHT_OF_D[0], y + 2.0, w=RIGHT_X - 1.0 - KEY_RIGHT_OF_D[0], ncol=KEY_RIGHT_OF_D[1], marks=KEY_NOTE)
            txmm(fig, H, b[1], y + BLOCKS[b[0]][2], LETTER[b[0]], fontsize=8, fontweight='bold', va='baseline'); k += 1   # fixed letters
        y += h + 1.5
    if note:
        txmm(fig, H, FW - 1, H - 1.0, f'DRAFT {note} — intervals: {CI}; length: mean over the {len(TA_EVAL)} type A tasks', fontsize=5, color='0.5', ha='right')
    fig.savefig(path, dpi=DPI, facecolor='white'); plt.close(fig)
    print('written', path, f'{FW:.0f} x {H:.0f} mm')


## ---- the panels for the slide: each block alone, without its panel letter (the letters belong to the slide), then the key ----
PANELS = [('a', 'methods_table', 'table'), ('b', 'DAM_example', 'example'), ('c', 'typeA', 'typeA'), ('d', 'typeB', 'typeB'),
          ('e', 'task_means', 'means'), ('f', 'reference_length', 'length')]   # panel letter, name of the file, block


def save_trimmed(fig, path):                                       # saved, then cut to its drawn content with a 0.5 mm margin
    fig.savefig(path, dpi=DPI, facecolor='white'); plt.close(fig)
    im = plt.imread(path)
    ys, xs = np.nonzero(ink(im))
    m = int(MM / 2)
    im = im[max(ys.min() - m, 0):ys.max() + m + 1, max(xs.min() - m, 0):xs.max() + m + 1]
    plt.imsave(path, im, dpi=DPI)
    print(f'{os.path.basename(path)}: {im.shape[1] / MM:.1f} x {im.shape[0] / MM:.1f} mm')


def render_panels(out):                                            # every block on a canvas of the figure's width, then the key
    os.makedirs(out, exist_ok=True)
    for letter, name, block in PANELS:
        h, draw, _ = BLOCKS[block]
        H = h + 2.0
        fig = plt.figure(figsize=(FW / 25.4, H / 25.4))
        draw(fig, H, 5.0, 1.0)                                     # x0 = 5 mm: room for labels at the left edge (trimmed below)
        save_trimmed(fig, f'{out}Fig5{letter}_{name}.png')
    H = key_h(1, marks=KEY_NOTE) + 2.0
    fig = plt.figure(figsize=(FW / 25.4, H / 25.4))
    draw_key(fig, H, 3, 1.0, w=RIGHT_X - 1.0 - KEY_RIGHT_OF_D[0], ncol=1, marks=KEY_NOTE)   # one column, as right of c
    save_trimmed(fig, f'{out}Fig5_methods_key.png')            # the key of the methods


if __name__ == '__main__':
    render(CURRENT, OUT + 'Fig5.png')                              # the whole figure (Word; the check of the slide)
    render_panels(OUT)                                             # the panels pasted in the slide
