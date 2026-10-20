"""Fig. 5f / Supplementary Table 3 — the comparisons of fig5_best_length.py over the eight type A tasks: the four
DAM and Glucocorticoid tasks (H1, D1, D2, D3) and the four PS2APP and Forced swim tasks (DA1 = the cells of H1 with the PS2APP
reference; AS1, AS2, AS3 = the cells of D1, D2, D3 with the Forced swim reference; fig5_extra_score.R).
The computation is that of fig5_best_length.py (replicated here, not imported, so that its four-task outputs are never rewritten):
  - per-cell scores of the benchmark; Beyondcell without a score = 0; AUC = Mann-Whitney with ties averaged
  - bootstrap over mice (fig5_mice.py: all cells of a drawn mouse): default_rng(1), one generator, 1,000 resamples, the same
    resamples for every method and length; the tasks ask for their draws in the order of fig5_best_length.py (H1, D1, D2, D3: the
    same draws), then DA1, AS1, AS2, AS3; the tasks with the same mice share the draws (Hammond: H1, DA1; De Donno: D1-D3, AS1-AS3).
  - mean over tasks: point = mean of the point AUCs; resample b = mean of the tasks' resampled AUCs; CI = 2.5% and 97.5% quantiles of
    RFind-sc minus the method; RFind-sc higher / lower when the CI excludes 0 and |difference| > 0.01, otherwise not distinguished;
    'comparable' when the whole CI lies within -0.01..+0.01
  - comparisons: fixed lengths 100, 250, 500, 1,000, 2,000, all; best length per task (each method, RFind-sc included, at the length
    with its highest AUC in that task); leave one task out (for task t, the length with the highest mean AUC over the other seven
    tasks); per task, also with the length re-selected in every resample
  - methods: the 10 main-figure methods and the table-only methods UCell (up genes), WMEAN (log2FC) and ORA
Check first: the four DAM and Glucocorticoid tasks must reproduce results/fig5_mean_auc_diff.csv and results/fig5_best_length.csv (every row and column, to
1e-9; the script stops otherwise, before the eight tasks). These two files are read, never written, here.
Inputs: figs/Fig5_benchmark/output/results/fig5_redo_percell_L<L>.csv.gz and results/fig5_redo_auc_L<L>.csv (H1, D1-D3),
        extra_refs/fig5_extra_percell_L<L>.csv.gz and extra_refs/fig5_extra_auc_L<L>.csv (DA1, AS1-AS3),
        results/fig5_mean_auc_diff.csv and results/fig5_best_length.csv (the check)
Outputs (figs/Fig5_benchmark/output/):
  results/fig5_best_length_8tasks.csv    the eight tasks, each method at its best length vs RFind-sc at its best, per task and leave one task
                                 out (the other seven tasks); columns of results/fig5_best_length.csv ('task set': main = the DAM and Glucocorticoid
                                 tasks, added = the PS2APP and Forced swim tasks)
  results/fig5_mean_auc_diff_8tasks.csv  mean AUC over the eight tasks, RFind-sc minus each method with 95% CI (columns of results/fig5_mean_auc_diff.csv)
Run from the project root (about 5 min, 0.3 GB of memory): python figs/Fig5_benchmark/analysis/fig5_best_length_8tasks.py
"""
import sys

import os
import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig5_mice import MouseDraws, mouse_of  # noqa: E402

O = 'figs/Fig5_benchmark/output/'
os.makedirs(O + 'results', exist_ok=True)                        # the tables in output/results/
X = O + 'extra_refs/'
LS = ['100', '250', '500', '1000', '2000', 'all']
METH = ['RFind-sc', 'AddModuleScore', 'AUCell', 'UCell', 'singscore', 'ssGSEA', 'VISION', 'Beyondcell', 'ULM (log2FC)', 'scFOCAL-type']
TONLY = ['UCell (up genes)', 'WMEAN (log2FC)', 'ORA']
ALLM = METH + TONLY
B, SEED = 1000, 1
TA4 = ['H1', 'D1', 'D2', 'D3']                                     # DAM and Glucocorticoid (drawn first, in the order of fig5_best_length.py)
ADD = {'DA1': 'H1', 'AS1': 'D1', 'AS2': 'D2', 'AS3': 'D3'}         # PS2APP and Forced swim: added task -> task on the same cells
TA8 = TA4 + list(ADD)                                              # drawing order of the added tasks: DA1, AS1, AS2, AS3
TSET = {**{t: 'main' for t in TA4}, **{t: 'added' for t in ADD}}


def aucs(S, pos):
    """AUC of every column of S (cells x methods) for the labels pos (Mann-Whitney, ties averaged)"""
    r = rankdata(S, axis=0); n1 = pos.sum(); n0 = len(pos) - n1
    return (r[pos].sum(axis=0) - n1 * (n1 + 1) / 2) / (n1 * n0)


def verdict(diff, lo, hi):
    return 'RFind-sc higher' if (lo > 0 and diff > 0.01) else ('RFind-sc lower' if (hi < 0 and diff < -0.01) else 'not distinguished')


def comparable(lo, hi):
    return bool(lo >= -0.01 and hi <= 0.01)


## ---- per-cell scores at every length (only the needed columns and tasks) and the stored point AUCs ----
UC = ['id', 'cell_id', 'positive'] + ALLM
pc, ref = {}, {}
for L in LS:
    parts = [ch[ch.id.isin(TA4)] for ch in pd.read_csv(f'{O}results/fig5_redo_percell_L{L}.csv.gz', usecols=UC, chunksize=100000)]
    e = pd.read_csv(f'{X}fig5_extra_percell_L{L}.csv.gz', usecols=UC)
    parts.append(e[e.id.isin(list(ADD))])
    d = pd.concat(parts, ignore_index=True)
    assert d.positive.dtype == bool and set(d.id) == set(TA8), L
    pc[L] = d.assign(Beyondcell=lambda x: x.Beyondcell.fillna(0))
    assert not pc[L][ALLM].isna().any().any(), L
    ref[L] = pd.concat([pd.read_csv(f'{O}results/fig5_redo_auc_L{L}.csv'), pd.read_csv(f'{X}fig5_extra_auc_L{L}.csv')],
                       ignore_index=True).set_index(['id', 'method']).AUC

## the added tasks score the same cells as the DAM and Glucocorticoid tasks, in the same order
for t, p in ADD.items():
    a, b = pc['500'][pc['500'].id == t], pc['500'][pc['500'].id == p]
    assert a.cell_id.tolist() == b.cell_id.tolist() and (a.positive.to_numpy() == b.positive.to_numpy()).all(), (t, p)

## ---- point AUCs and bootstrap AUCs (B x columns) per task; columns = (method, length) ----
DRAWS = MouseDraws(np.random.default_rng(SEED), B)
point, boot, col = {}, {}, {}
maxdiff = 0.0


def prepare_and_draw(t):                                           # scores (cells x (method, length)), point AUCs, then B resamples
    global maxdiff
    base = pc['500'][pc['500'].id == t][['cell_id', 'positive']].reset_index(drop=True)
    cols, names = [], []
    for L in LS:
        x = pc[L][pc[L].id == t].set_index('cell_id').loc[base.cell_id]
        assert (x.positive.to_numpy() == base.positive.to_numpy()).all()
        for m in ALLM:
            cols.append(x[m].to_numpy(dtype=float)); names.append((m, L))
    S = np.column_stack(cols); pos = base.positive.to_numpy(dtype=bool)
    for (m, L), v in zip(names, aucs(S, pos)):
        point[(t, m, L)] = v
        assert abs(v - ref[L][(t, m)]) < 1e-5, (t, m, L, v, ref[L][(t, m)])   # the stored point AUCs are reproduced (the per-cell
        maxdiff = max(maxdiff, abs(v - ref[L][(t, m)]))                         # files keep ~15 digits, so near-ties can move ~1e-7)
    col[t] = {nm: k for k, nm in enumerate(names)}
    bt = np.empty((B, S.shape[1]))
    for b, i in enumerate(DRAWS.cells(mouse_of(base.cell_id), pos)):
        bt[b] = aucs(S[i], pos[i])
    boot[t] = bt


def analyse(tasks):                                                # fig5_best_length.py for the task set 'tasks'
    def bmat(t, m):
        return boot[t][:, [col[t][(m, L)] for L in LS]]

    best_task = {(t, m): max(LS, key=lambda L: point[(t, m, L)]) for t in tasks for m in ALLM}
    best_loto = {(t, m): max(LS, key=lambda L: np.mean([point[(u, m, L)] for u in tasks if u != t])) for t in tasks for m in ALLM}
    sel_task = {(t, m): bmat(t, m).argmax(axis=1) for t in tasks for m in ALLM}
    sel_loto = {(t, m): np.mean([bmat(u, m) for u in tasks if u != t], axis=0).argmax(axis=1) for t in tasks for m in ALLM}

    def pick(t, m, sel):
        return bmat(t, m)[np.arange(B), sel]

    rows = []
    for version in ['per task', 'leave one task out']:
        for t in tasks:
            if version == 'per task':
                bl = {m: best_task[(t, m)] for m in ALLM}; sl = {m: sel_task[(t, m)] for m in ALLM}
            else:
                bl = {m: best_loto[(t, m)] for m in ALLM}; sl = {m: sel_loto[(t, m)] for m in ALLM}
            rf = boot[t][:, col[t][('RFind-sc', bl['RFind-sc'])]]
            rf_sel = pick(t, 'RFind-sc', sl['RFind-sc'])
            share = float(np.mean(sl['RFind-sc'] == LS.index(bl['RFind-sc'])))
            for m in ALLM[1:]:
                mb = boot[t][:, col[t][(m, bl[m])]]
                diff = point[(t, 'RFind-sc', bl['RFind-sc'])] - point[(t, m, bl[m])]
                lo, hi = np.quantile(rf - mb, [0.025, 0.975])
                slo, shi = np.quantile(rf_sel - pick(t, m, sl[m]), [0.025, 0.975])
                rows.append({'version': version, 'task id': t, 'method': m, 'RFind-sc length': bl['RFind-sc'],
                             'RFind-sc AUC': point[(t, 'RFind-sc', bl['RFind-sc'])], 'method length': bl[m],
                             'method AUC': point[(t, m, bl[m])], 'RFind-sc minus method': diff, 'difference 95% CI low': lo,
                             'difference 95% CI high': hi, 'result': verdict(diff, lo, hi), 'comparable': comparable(lo, hi),
                             'task set': TSET[t], 'method set': 'main' if m in METH else 'table only',
                             'reselected difference 95% CI low': slo, 'reselected difference 95% CI high': shi,
                             'reselected result': verdict(diff, slo, shi), 'reselected comparable': comparable(slo, shi),
                             'RFind-sc share of resamples at its point-estimate length': share})
    out = pd.DataFrame(rows)

    def mean_point_boot(m, kind):
        if kind in LS:
            p = np.mean([point[(t, m, kind)] for t in tasks]); bb = np.mean([boot[t][:, col[t][(m, kind)]] for t in tasks], axis=0)
        elif kind == 'best':
            p = np.mean([point[(t, m, best_task[(t, m)])] for t in tasks])
            bb = np.mean([boot[t][:, col[t][(m, best_task[(t, m)])]] for t in tasks], axis=0)
        else:
            p = np.mean([point[(t, m, best_loto[(t, m)])] for t in tasks])
            bb = np.mean([boot[t][:, col[t][(m, best_loto[(t, m)])]] for t in tasks], axis=0)
        return p, bb

    mrows = []
    for cname, kind in [(f'fixed length {L}', L) for L in LS] + [('best length per task', 'best'), ('leave one task out', 'loto')]:
        pr, br = mean_point_boot('RFind-sc', kind)
        for m in ALLM[1:]:
            p, bb = mean_point_boot(m, kind)
            diff = pr - p
            lo, hi = np.quantile(br - bb, [0.025, 0.975])
            mrows.append({'comparison': cname, 'method': m, 'method set': 'main' if m in METH else 'table only', 'mean AUC RFind-sc': pr,
                          'mean AUC method': p, 'difference': diff, 'CI low': lo, 'CI high': hi, 'result': verdict(diff, lo, hi),
                          'comparable': comparable(lo, hi)})
    return out, pd.DataFrame(mrows), best_task


## ---- 1. the four DAM and Glucocorticoid tasks with the draws of fig5_best_length.py; check against its outputs, stop if not reproduced ----
for t in TA4:
    prepare_and_draw(t)
out4, md4, best4 = analyse(TA4)
pub_md = pd.read_csv(f'{O}results/fig5_mean_auc_diff.csv'); pub_bl = pd.read_csv(f'{O}results/fig5_best_length.csv')


def same(mine, pub):                                               # (max |numeric difference|, list of differing text columns)
    assert list(mine.columns) == list(pub.columns) and len(mine) == len(pub), 'columns or rows differ'
    num = [c for c in pub.columns if pd.api.types.is_float_dtype(pub[c])]
    dmax = max(float(np.abs(mine[c].to_numpy(float) - pub[c].to_numpy(float)).max()) for c in num)
    txt = [c for c in pub.columns if c not in num and not (mine[c].astype(str).to_numpy() == pub[c].astype(str).to_numpy()).all()]
    return dmax, txt


d_md, t_md = same(md4, pub_md); d_bl, t_bl = same(out4, pub_bl)
print(f'point AUCs vs the stored AUC tables (4 tasks): max |diff| {maxdiff:.1e}')
print(f'check, 4 tasks vs results/fig5_mean_auc_diff.csv: max |diff| {d_md:.1e}, differing text columns {t_md or "none"}')
print(f'check, 4 tasks vs results/fig5_best_length.csv:   max |diff| {d_bl:.1e}, differing text columns {t_bl or "none"}')
if d_md > 1e-9 or d_bl > 1e-9 or t_md or t_bl:
    sys.exit('the four DAM and Glucocorticoid tasks are NOT reproduced: stopped before the eight tasks')
print('  reproduced exactly (to 1e-9): continuing with the eight tasks')

## ---- 2. the added tasks (continuing the generator), then the eight tasks ----
for t in ADD:
    prepare_and_draw(t)
print(f'point AUCs vs the stored AUC tables (8 tasks): max |diff| {maxdiff:.1e}')
out8, md8, best8 = analyse(TA8)
out8.to_csv(O + 'results/fig5_best_length_8tasks.csv', index=False)
md8.to_csv(O + 'results/fig5_mean_auc_diff_8tasks.csv', index=False)
print('written', O + 'results/fig5_best_length_8tasks.csv', out8.shape, '|', O + 'results/fig5_mean_auc_diff_8tasks.csv', md8.shape)


## ---- 3. summary (the numbers of Results 5, Supplementary Note 4 and Supplementary Table 3) ----
def counts(g, rc='result'):
    c = g[rc].value_counts()
    return f'{c.get("RFind-sc higher", 0)}/{c.get("not distinguished", 0)}/{c.get("RFind-sc lower", 0)}'


def f3(x):
    return f'{x:.3f}'


main9 = METH[1:]
for lab, md in [('4 DAM and Glucocorticoid tasks', md4), ('8 type A tasks', md8)]:
    print(f'\n== mean AUC over the {lab} (RFind-sc minus ULM (log2FC), 95% CI; vs the 9 main methods higher/not distinguished/lower)')
    for cname in md.comparison.unique():
        g = md[(md.comparison == cname) & md.method.isin(main9)]
        u = g[g.method == 'ULM (log2FC)'].iloc[0]
        top = g.loc[g['mean AUC method'].idxmax()]
        print(f'  {cname:22s} RFind-sc {u["mean AUC RFind-sc"]:.4f} | ULM {u["mean AUC method"]:.4f} diff {u.difference:+.4f} '
              f'[{u["CI low"]:+.4f}, {u["CI high"]:+.4f}] {u.result}{" (comparable)" if u.comparable else ""} | highest other '
              f'{top.method} {f3(top["mean AUC method"])} | 9 methods {counts(g)}')
for lab, out, tasks, best in [('4 DAM and Glucocorticoid tasks', out4, TA4, best4), ('8 type A tasks', out8, TA8, best8)]:
    g = out[(out.version == 'per task') & out.method.isin(main9)]
    print(f'\n== best length per task, {lab}: RFind-sc vs the 8 methods other than ULM {counts(g[g.method != "ULM (log2FC)"])} '
          f'(re-selected {counts(g[g.method != "ULM (log2FC)"], "reselected result")}); vs ULM {counts(g[g.method == "ULM (log2FC)"])} '
          f'(re-selected {counts(g[g.method == "ULM (log2FC)"], "reselected result")})')
    for t in tasks:
        r = g[(g['task id'] == t) & (g.method == 'ULM (log2FC)')].iloc[0]
        others = sorted({best[(t, m)] for m in main9 if m != 'ULM (log2FC)'}, key=LS.index)
        print(f'  {t}: RFind-sc L {r["RFind-sc length"]} AUC {f3(r["RFind-sc AUC"])} | ULM L {r["method length"]} {f3(r["method AUC"])} '
              f'diff {r["RFind-sc minus method"]:+.3f} [{r["difference 95% CI low"]:+.3f}, {r["difference 95% CI high"]:+.3f}] '
              f'{r.result} | 8 others: {counts(g[(g["task id"] == t) & (g.method != "ULM (log2FC)")])}, their lengths {others}')
    gl =out[(out.version == 'leave one task out') & out.method.isin(main9)]
    print(f'  leave one task out: vs the 8 methods other than ULM {counts(gl[gl.method != "ULM (log2FC)"])}; vs ULM '
          f'{counts(gl[gl.method == "ULM (log2FC)"])} (re-selected {counts(gl[gl.method == "ULM (log2FC)"], "reselected result")})')
    low = g[g.result == 'RFind-sc lower']
    print('  lower (per task):', [(r['task id'], r.method, r['method length'], round(r['RFind-sc minus method'], 3)) for _, r in low.iterrows()])
    ch = g[g.result != g['reselected result']]
    print('  outcomes changed by re-selection (per task):', [(r['task id'], r.method, r.result, r['reselected result']) for _, r in ch.iterrows()] or 'none')
