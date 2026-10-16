"""Supplementary Fig. S4 — the tasks of the Fig. 5 benchmark and the four added type A tasks.
Added tasks (fig5_extra_score.R): DA1 = the cells of the DAM task (H1; LPC- vs saline-injected white matter, microglia) with the
PS2APP reference; AS1, AS2, AS3 = the cells of the three Glucocorticoid tasks (D1-D3; stress vs no stress in Oligodendrocytes,
Astro and Endothelial) with the Forced swim reference. Fig. 5c shows the four tasks of the paper; Fig. 5e averages all eight.
Panels a (table of the tasks) and b (added tasks), drawn without panel letters (the slide adds them), no figure title, no
notes in the image (the notes are in the legend; the numbers they need are printed):
  part 0  every task of the benchmark (a), in the style of the table of Fig. 5a (cell numbers from the sheet Tasks of
          Supplementary Table 3)
  part 1  the four added tasks: AUC with 95% CI of every method of Fig. 5c at 500 genes per direction, drawn with the functions and
          the method colours of fig5_rebuild_figure_draft.py (no key; the colours are those of Fig. 5). A panel's lower edge goes
          below the AUC range of 0.2 only to show a whole CI; points below the edge are triangles, as in Fig. 5.
The name of the script refers to the capability map of the benchmark (type A with each method at its best reference length, and
type B); these comparisons are in Supplementary Table 3 (sheets Best length, Comparisons and Mean AUC differences), not in this figure.
All text >= 5 pt (the subscript 2 of G2M too: drawn at 5 pt by SubText), width 180 mm, 600 dpi.
Inputs: figs/Fig5_benchmark/output/extra_refs/fig5_extra_auc_L500.csv, SuppTable3.xlsx (sheet Tasks); the module
        fig5_rebuild_figure_draft.py (imported: data and drawing functions of Fig. 5; on import it also reads
        Fig5c_embedding.png and Fig5d_auc.png, the outputs of fig5_redo_figure.py, which therefore runs first)
Outputs (figs/Supp_S4_Benchmark/output/; the folder of Supplementary Fig. S4): FigS4.png (the whole figure, a and b, without panel
        letters; Word and the check of the slide) and the panels for the slide, directly in output/: FigS4a_tasks.png and
        FigS4b_added_tasks.png, each alone and without its letter, cut to its drawn content with a 0.5 mm margin
Run from the project root: python figs/Supp_S4_Benchmark/analysis/figS4_capability_map.py
"""
import os
import sys

sys.dont_write_bytecode = True                                     # no __pycache__ written into figs/
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'Fig5_benchmark', 'analysis'))   # fig5_rebuild_figure_draft
import numpy as np                                                 # noqa: E402
import pandas as pd                                                # noqa: E402
import matplotlib.pyplot as plt                                    # noqa: E402
from matplotlib.artist import Artist, allow_rasterization         # noqa: E402
from matplotlib.backends.backend_agg import get_hinting_flag      # noqa: E402
from matplotlib.font_manager import FontProperties, findfont, get_font   # noqa: E402
from matplotlib.ft2font import Kerning                             # noqa: E402
from matplotlib.text import Text                                 # noqa: E402
from matplotlib.transforms import IdentityTransform               # noqa: E402
from PIL import Image                                              # noqa: E402
import fig5_rebuild_figure_draft as F                              # noqa: E402  data and drawing functions of Fig. 5

O = F.O                                                            # figs/Fig5_benchmark/output/
RF, CM = F.RF, F.CMETH
PAPER = list(F.TA)                                                 # H1, D1, D2, D3 (Fig. 5c)
ADDED = list(F.EXTRA)                                              # DA1, AS1, AS2, AS3
TYPEA = PAPER + ADDED
assert F.TA_EVAL == TYPEA
## names of the added tasks as in the paper (PS2APP, Forced swim); the other task names are those of Fig. 5
TNAME = {**F.TNAME, 'DA1': ('PS2APP', F.TNAME['DA1'][1]), **{t: ('Forced swim', F.TNAME[t][1]) for t in ('AS1', 'AS2', 'AS3')}}

## ---- AUC of the added tasks at 500 genes (each task once with each method of Fig. 5c) ----
X5 = pd.read_csv(f'{O}extra_refs/fig5_extra_auc_L500.csv')
for t in ADDED:
    n = X5[(X5.id == t) & X5.method.isin(CM)].method.value_counts()
    assert len(n) == len(CM) and n.max() == 1, t
XA = X5[X5.method.isin(CM)][F.A.columns]

## ---- part 1: the added tasks, drawn as the panels of Fig. 5c ----
SPAN0 = F.SPAN                                                     # 0.2: the AUC range of the small panels
FLOOR = {}


def small_c(ax, v, t):                                             # F.small; the lower edge goes down just enough (to 0.005 below
    v = v[v.index.isin(CM) & np.isfinite(v.AUC)]                   # the CI) when the 95% CI of a point inside the range would be
    yhi = min(1.012, np.where(np.isfinite(v.hi), v.hi, v.AUC).max() + 0.012)   # cut there; points below the edge stay triangles
    ylo = yhi - SPAN0
    while ((v.AUC >= ylo) & (v.lo < ylo)).any():
        ylo = v[(v.AUC >= ylo) & (v.lo < ylo)].lo.min() - 0.005
    if yhi - ylo > SPAN0 + 1e-12:
        FLOOR[t] = (ylo, yhi)
    F.SPAN = yhi - ylo
    try:
        F.small(ax, v, CM, list(range(len(CM))))
    finally:
        F.SPAN = SPAN0


P1_H = 34.0                                                        # height of part 1 (mm): task labels 8, panels 24, 2 below


def draw_part1(fig, H, y0, X1=9.3, GAP=5.8, letter=True):          # the geometry of Fig. 5c (4 small panels of the same width)
    PW, PH = F.PW_S, 24.0                                           # the small panels of Fig. 5c (b about 118 mm wide)
    for k, t in enumerate(ADDED):
        xa = X1 + k * (PW + GAP)
        ax = F.axmm(fig, H, xa, y0 + 8, PW, PH)
        if k == 0 and letter:
            F.txmm(fig, H, 0.2, y0 + 3.6, 'b', fontsize=12, fontweight='bold', va='baseline')   # 12 pt as on the slides
        F.txmm(fig, H, xa + PW / 2, y0 + 3.6, TNAME[t][0], fontsize=F.T6, fontweight='bold', ha='center', va='baseline')
        F.txmm(fig, H, xa + PW / 2, y0 + 6.1, TNAME[t][1], fontsize=F.T6, ha='center', va='baseline', color='0.3')
        small_c(ax, XA[XA.id == t].set_index('method'), t)
        ax.set_xlim(-0.8, len(CM) - 0.2)
        if k == 0:
            ax.set_ylabel('AUC (95% CI)')


## ---- part 0: every task of the benchmark (a, the summary of all tasks used) ----
TK = pd.read_excel(f'{O}SuppTable3.xlsx', sheet_name='Tasks')


def ncell(prefix, which):                                          # the cell numbers of a task of the Tasks sheet
    r = TK[TK.task.str.startswith(prefix)]
    assert len(r) == 1, prefix
    return f"{int(r.iloc[0]['positive cells' if which == 'pos' else 'negative cells']):,}"


## ---- G2M with a subscript 2: the subscript at 5 pt, the smallest text of the figures (mathtext would draw it at
## 0.7 x 6 = 4.2 pt) ----
G2M = 'G2M'                                                        # in a label, drawn by SubText (G, a 5-pt '2', then M); keys untouched
SUB_PT = 5.0                                                       # size of the subscript 2 (pt)


class SubText(Artist):
    """One line of text in which 'G2M' is drawn as G, a subscript 2 of SUB_PT and M; the rest at `fontsize`, all in the same font and
    weight. The 2 sits on a baseline lowered by the font's own subscript offset (OS/2 ySubscriptYOffset; 283/2048 em for Arial and
    Arial Bold) and takes its own advance at SUB_PT. The line is placed where ax.text would place it in plain text: ha 'left' or
    'center' (centred on the advance width), va 'baseline' or 'center' (line `line` of a block of `nlines` lines with `linespacing`,
    laid out as matplotlib lays out lines no taller than 'lp'). Positions are taken from the renderer when drawn (pen advances at the
    output dpi), so the pieces join as one word."""

    def __init__(self, x, y, s, fontsize, ha='left', va='center', line=0, nlines=1, linespacing=1.2, transform=None, **kw):
        super().__init__()
        i = s.index(G2M)
        self._xy, self._ha, self._va = (x, y), ha, va
        self._line, self._nlines, self._ls = line, nlines, linespacing
        self._parts = [Text(0, 0, p, fontsize=size, ha='left', va='baseline', transform=IdentityTransform(), clip_on=False, **kw)
                       for p, size in ((s[:i + 1], fontsize), ('2', SUB_PT), (s[i + 2:], fontsize))]
        if transform is not None:
            self.set_transform(transform)
        self.set_zorder(3)                                         # as Text

    def set_figure(self, fig):
        super().set_figure(fig)
        for t in self._parts:
            t.set_figure(fig)

    @staticmethod
    def _advance(t, renderer):                                     # pen advance of a piece (px), as RendererAgg.draw_text lays it out
        s, prop = t.get_text(), t.get_fontproperties()
        font = get_font(findfont(prop))
        font.set_size(prop.get_size_in_points(), renderer.dpi)
        xys = font.set_text(s + s[-1], 0.0, flags=get_hinting_flag())   # the last glyph repeated: its pen position ends s
        return xys[-1][0] / 64.0                                   # (no kerning of a glyph with itself in Arial: checked below)

    @allow_rasterization
    def draw(self, renderer):
        if not self.get_visible():
            return
        main, sub, rest = self._parts
        x, y = self.get_transform().transform(self._xy)
        prop = main.get_fontproperties()
        _, lp_h, lp_d = renderer.get_text_width_height_descent('lp', prop, ismath=False)
        if self._va == 'center':
            step = lp_d + self._ls * (lp_h - lp_d)                 # baseline to baseline
            y = y + (lp_h + (self._nlines - 1) * step) / 2 - (lp_h - lp_d) - self._line * step
        font = get_font(findfont(prop))
        drop = font.get_sfnt_table('OS/2')['ySubscriptYOffset'] / font.units_per_EM * renderer.points_to_pixels(prop.get_size_in_points())
        if renderer.flipy():                                       # Agg puts each piece on whole pixels (round(y + descent)): the
            hgt = renderer.get_canvas_width_height()[1]            # baseline 1/4 px past a whole pixel and the drop in whole pixels,
            y = hgt - (np.floor(hgt - y) + 0.25)                   # so that G and M share one baseline, as in one string
            drop = float(round(drop))
        a1, a2 = self._advance(main, renderer), self._advance(sub, renderer)
        if self._ha == 'center':
            x -= (a1 + a2 + self._advance(rest, renderer)) / 2
        for t, (dx, dy) in zip(self._parts, ((0.0, 0.0), (a1, -drop), (a1 + a2, 0.0))):
            t.set_position((x + dx, y + dy))
            t.draw(renderer)
        self.stale = False


def text_g2m(ax, x, y, s, **kw):                                   # ax.text, with G2M drawn by SubText (data coordinates of ax)
    if G2M not in s:
        return ax.text(x, y, s, **kw)
    return ax.add_artist(SubText(x, y, s, transform=ax.transData, **kw))


for _w in ('normal', 'bold'):                                      # the font of the labels: no kerning of a glyph with itself
    _f = get_font(findfont(FontProperties(weight=_w)))             # (SubText._advance relies on it)
    assert all(_f.get_kerning(_f.get_char_index(ord(c)), _f.get_char_index(ord(c)), Kerning.DEFAULT) == 0 for c in 'G2Mt'), _w


TROWS = [  # task, cells, positive (n), negative (n), reference, source, shown in; None = a group row
    ('Type A: ranked lists of up- and down-regulated genes', None),
    ('DAM', 'Microglia, Fig. 2', f"LPC-injected ({ncell('DAM:', 'pos')})", f"Saline-injected ({ncell('DAM:', 'neg')})",
     'DAM vs homeostatic microglia', 'Keren-Shaul et al., 2017', 'Fig. 5c'),
    ('PS2APP', 'Microglia, Fig. 2', f"LPC-injected ({ncell('PS2APP:', 'pos')})", f"Saline-injected ({ncell('PS2APP:', 'neg')})",
     'PS2APP vs non-transgenic microglia', 'Srinivasan et al., 2016', 'S4b'),
    ('Glucocorticoid', 'Oligodendrocytes, Fig. 3', f"Stress ({ncell('Glucocorticoid: stress vs no stress, oligodendrocytes', 'pos')})",
     f"No stress ({ncell('Glucocorticoid: stress vs no stress, oligodendrocytes', 'neg')})", 'Corticosterone 3 h, ventral hippocampus', 'Viho et al., 2022', 'Fig. 5c'),
    ('Glucocorticoid', 'Astro, Fig. 3', f"Stress ({ncell('Glucocorticoid: stress vs no stress, astrocytes', 'pos')})",
     f"No stress ({ncell('Glucocorticoid: stress vs no stress, astrocytes', 'neg')})", 'Corticosterone 3 h, ventral hippocampus', 'Viho et al., 2022', 'Fig. 5c'),
    ('Glucocorticoid', 'Endothelial, Fig. 3', f"Stress ({ncell('Glucocorticoid: stress vs no stress, endothelial', 'pos')})",
     f"No stress ({ncell('Glucocorticoid: stress vs no stress, endothelial', 'neg')})", 'Corticosterone 3 h, ventral hippocampus', 'Viho et al., 2022', 'Fig. 5c'),
    ('Forced swim', 'Oligodendrocytes, Fig. 3', f"Stress ({ncell('Forced swim: stress vs no stress, oligodendrocytes', 'pos')})",
     f"No stress ({ncell('Forced swim: stress vs no stress, oligodendrocytes', 'neg')})", 'Forced swim 3 h, ventral hippocampus', 'von Ziegler et al., 2022', 'S4b'),
    ('Forced swim', 'Astro, Fig. 3', f"Stress ({ncell('Forced swim: stress vs no stress, astrocytes', 'pos')})",
     f"No stress ({ncell('Forced swim: stress vs no stress, astrocytes', 'neg')})", 'Forced swim 3 h, ventral hippocampus', 'von Ziegler et al., 2022', 'S4b'),
    ('Forced swim', 'Endothelial, Fig. 3', f"Stress ({ncell('Forced swim: stress vs no stress, endothelial', 'pos')})",
     f"No stress ({ncell('Forced swim: stress vs no stress, endothelial', 'neg')})", 'Forced swim 3 h, ventral hippocampus', 'von Ziegler et al., 2022', 'S4b'),
    ('Type B: gene sets', None),
    ('Interferon', 'Microglia, Fig. 2', f"Cluster 9 ({ncell('Interferon: cluster 9', 'pos')})", f"Other microglia ({ncell('Interferon: cluster 9', 'neg')})",
     'Hallmark interferon alpha response', 'Liberzon et al., 2015', 'Fig. 5d'),
    ('Interferon', 'Microglia, Fig. 2', f"P540 ({ncell('Interferon: P540', 'pos')})", f"P100 ({ncell('Interferon: P540', 'neg')})",
     'Hallmark interferon alpha response', 'Liberzon et al., 2015', 'Fig. 5d'),
    (G2M, 'Microglia, Fig. 2', f"Cluster 2b ({ncell('G2M', 'pos')})", f"Other microglia ({ncell('G2M', 'neg')})",
     f'Hallmark {G2M} checkpoint', 'Liberzon et al., 2015', 'Fig. 5d')]
TCW = [22, 30, 25, 31, 44, 31, 15]                                 # column widths (mm), sum 198 -> scaled to 179 below
TCW = [w * 179.0 / sum(TCW) for w in TCW]
THEAD = ['Task', 'Cells', 'Positive (n)', 'Negative (n)', 'Reference', 'Source', 'Shown in']
TRH, THH = 2.9, 4.6                                                # row and header heights (mm)
P0_H = THH + TRH * len(TROWS) + 1.0                                # height of part 0 (mm)


def draw_part0(fig, H, y0, letter=True):                           # the style of Fig. 5a (header bar, group rows, rules)
    ax = F.axmm(fig, H, 0.5, y0 + 4.2, sum(TCW), P0_H); ax.set_xlim(0, sum(TCW)); ax.set_ylim(P0_H, 0); ax.axis('off')
    if letter:
        F.txmm(fig, H, 0.2, y0 + 3.6, 'a', fontsize=12, fontweight='bold', va='baseline')
    ax.add_patch(plt.Rectangle((0, 0), sum(TCW), THH, color=F.NAVY, lw=0))
    xs = np.cumsum([0] + TCW[:-1])
    for x, h in zip(xs, THEAD):
        ax.text(x + 1.2, THH / 2, h, color='white', fontsize=F.T6, fontweight='bold', ha='left', va='center')
    y = THH
    for r in TROWS:
        if r[1] is None:
            ax.text(1.2, y + TRH * 0.62, r[0], fontsize=F.T6, color='0.35', va='center'); ax.plot([0, sum(TCW)], [y, y], color='0.6', lw=0.4)
            y += TRH; continue
        for k, x in enumerate(xs):
            text_g2m(ax, x + 1.2, y + TRH / 2, r[k], fontsize=F.T6, ha='left', va='center', fontweight='bold' if k == 0 else 'normal')
        y += TRH
    ax.plot([0, sum(TCW)], [y, y], color='0.6', lw=0.4)


## ---- outputs: the whole S4 (180 mm wide; a and b) and its panels a and b, all directly in output/ ----
O4 = 'figs/Supp_S4_Benchmark/output/'; os.makedirs(O4, exist_ok=True)   # the outputs of S4 (the inputs are in figs/Fig5_benchmark/output/)
P0_T = P0_H + 5.2                                                  # part 0 (a) with room above the table
Hp = P1_H + 2.0
GAPP = 4.0                                                         # between the parts (mm)
HS = 0.5 + P0_T + GAPP + P1_H + 0.5
fig = plt.figure(figsize=(F.FW / 25.4, HS / 25.4))
draw_part0(fig, HS, 0.5, letter=False)                             # no panel letters in the figure (the slide adds them)
draw_part1(fig, HS, 0.5 + P0_T + GAPP, letter=False)
OUTS = [O4 + 'FigS4.png']
fig.savefig(OUTS[0], dpi=F.DPI, facecolor='white'); plt.close(fig)
## the panels for the slide, each alone and without its letter
OUTS.append(O4 + 'FigS4a_tasks.png')
fig = plt.figure(figsize=(F.FW / 25.4, (P0_T + 1.0) / 25.4)); draw_part0(fig, P0_T + 1.0, 0.0, letter=False); F.save_trimmed(fig, OUTS[-1])
OUTS.append(O4 + 'FigS4b_added_tasks.png')
fig = plt.figure(figsize=(F.FW / 25.4, Hp / 25.4)); draw_part1(fig, Hp, 1.0, letter=False); F.save_trimmed(fig, OUTS[-1])

## ---- report: the numbers that the legend needs ----
for t, (lo, hi) in FLOOR.items():
    print(f'part 1, {t}: lower edge {lo:.3f} instead of {hi - SPAN0:.3f} (range {hi - lo:.3f}), so that no 95% CI is cut')
for t in ADDED:
    v = X5[(X5.id == t) & X5.method.isin(CM)].set_index('method'); top = v.drop(index=RF).AUC.idxmax()
    print(f'part 1, {t}: RFind-sc AUC {v.loc[RF, "AUC"]:.3f}; highest other {top} {v.loc[top, "AUC"]:.3f}; cells '
          f'{int(v.loc[RF, "n_pos"]):,} positive / {int(v.loc[RF, "n_neg"]):,} negative')
for p in OUTS:  # the whole figure, then the panels a and b
    with Image.open(p) as im:
        print(f'{os.path.basename(p)}: {im.size[0]} x {im.size[1]} px, {im.size[0] / F.MM:.1f} x {im.size[1] / F.MM:.1f} mm, '
              f'dpi {tuple(round(v) for v in im.info.get("dpi", (0, 0)))}')
