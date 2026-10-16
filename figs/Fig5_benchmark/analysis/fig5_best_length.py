"""Fig. 5 — each method at its best reference length (type A tasks). For every method (RFind-sc included) and task, the reference length
among 100, 250, 500, 1,000, 2,000 and all genes per direction with the highest AUC is taken ('per task'); as a second version, for
task t the length with the highest mean AUC over the other three main tasks ('leave one task out'). RFind-sc at its best length is then
compared with each method at its best length by the 95% CI of the difference in AUC from a bootstrap over cells (1,000 resamples, the
same resamples for every method and length; RFind-sc higher / lower when the CI excludes 0 and |difference| > 0.01, otherwise not
distinguished; 'comparable' when the whole CI lies within -0.01..+0.01).
Re-selection columns: the best length is chosen again inside every bootstrap resample (for each method, RFind-sc included, the length
with the highest resampled AUC in that resample; for 'leave one task out', the highest mean of the resampled AUCs over the other three
main tasks), the resampled difference uses those lengths; the point difference stays the point-estimate one.
The share of resamples in which RFind-sc chose its point-estimate length is reported.
Methods: the 10 main-figure methods and the table-only methods UCell (up genes), WMEAN (log2FC) and ORA.
Tasks: the four type A tasks H1, D1, D2, D3. The same comparisons over the eight type A tasks (these four and the four added tasks of
Supplementary Fig. S4) are made by fig5_best_length_8tasks.py, which reproduces the two outputs of this script first.
Random draws: default_rng(1), one generator for all tasks, drawn task by task in the order H1, D1, D2, D3.
Scores: the per-cell scores of the benchmark (results/fig5_redo_percell_L<L>.csv.gz); Beyondcell without a score = 0, as in the benchmark.
The point AUCs are checked against results/fig5_redo_auc_L<L>.csv.
Outputs (figs/Fig5_benchmark/output/):
  results/fig5_best_length.csv   per task / leave one task out versions, with re-selection columns and 'task set'
  results/fig5_mean_auc_diff.csv mean AUC over the 4 main tasks, RFind-sc minus each method with 95% CI (resample b: mean over tasks of the
                         resampled AUCs), per fixed length, best length per task, leave one task out
Run from the project root: python figs/Fig5_benchmark/analysis/fig5_best_length.py
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import rankdata

O = 'figs/Fig5_benchmark/output/'
os.makedirs(O + 'results', exist_ok=True)                        # the tables in output/results/
LS = ['100', '250', '500', '1000', '2000', 'all']
TA = ['H1', 'D1', 'D2', 'D3']                       # main tasks (bootstrap draws first, in this order)
TN = []                                             # no other tasks (an empty list; the code below runs over TA + TN)
METH = ['RFind-sc', 'AddModuleScore', 'AUCell', 'UCell', 'singscore', 'ssGSEA', 'VISION', 'Beyondcell', 'ULM (log2FC)', 'scFOCAL-type']
TONLY = ['UCell (up genes)', 'WMEAN (log2FC)', 'ORA']
ALLM = METH + TONLY
B, SEED = 1000, 1
TSET = {**{t: 'main' for t in TA}, **{t: 'neuron (not in the main figure)' for t in TN}}


def aucs(S, pos):
    """AUC of every column of S (cells x methods) for the labels pos (Mann-Whitney, ties averaged)"""
    r = rankdata(S, axis=0); n1 = pos.sum(); n0 = len(pos) - n1
    return (r[pos].sum(axis=0) - n1 * (n1 + 1) / 2) / (n1 * n0)


def verdict(diff, lo, hi):
    return 'RFind-sc higher' if (lo > 0 and diff > 0.01) else ('RFind-sc lower' if (hi < 0 and diff < -0.01) else 'not distinguished')


def comparable(lo, hi):
    return bool(lo >= -0.01 and hi <= 0.01)


## per-cell scores at every length, aligned by cell
pc = {}
for L in LS:
    d = pd.read_csv(f'{O}results/fig5_redo_percell_L{L}.csv.gz', usecols=['id', 'cell_id', 'positive'] + ALLM)
    pc[L] = d[d.id.isin(TA + TN)].assign(Beyondcell=lambda x: x.Beyondcell.fillna(0))
    assert not pc[L][ALLM].isna().any().any(), L
ref = {L: pd.read_csv(f'{O}results/fig5_redo_auc_L{L}.csv').set_index(['id', 'method']).AUC for L in LS}

## point AUCs and bootstrap AUCs (B x columns) per task; columns = (method, length)
rng = np.random.default_rng(SEED)
point, boot, col = {}, {}, {}
maxdiff = 0.0
for t in TA + TN:
    base = pc['500'][pc['500'].id == t][['cell_id', 'positive']].reset_index(drop=True)
    cols, names = [], []
    for L in LS:
        x = pc[L][pc[L].id == t].set_index('cell_id').loc[base.cell_id]
        assert (x.positive.to_numpy() == base.positive.to_numpy()).all()
        for m in ALLM:
            cols.append(x[m].to_numpy(dtype=float)); names.append((m, L))
    S = np.column_stack(cols); pos = base.positive.to_numpy(dtype=bool)
    a = aucs(S, pos)
    for (m, L), v in zip(names, a):
        point[(t, m, L)] = v
        assert abs(v - ref[L][(t, m)]) < 1e-5, (t, m, L, v, ref[L][(t, m)])     # the stored point AUCs are reproduced (the per-cell
        maxdiff = max(maxdiff, abs(v - ref[L][(t, m)]))                         # files keep ~15 digits, so near-ties can move ~1e-7)
    bt = np.empty((B, S.shape[1]))
    n = len(pos)
    for b in range(B):
        i = rng.integers(0, n, n)
        bt[b] = aucs(S[i], pos[i])
    boot[t] = bt
    col[t] = {nm: k for k, nm in enumerate(names)}


def bmat(t, m):
    """resampled AUCs of method m on task t, B x lengths (order LS)"""
    return boot[t][:, [col[t][(m, L)] for L in LS]]


## point-estimate lengths per version
best_task = {(t, m): max(LS, key=lambda L: point[(t, m, L)]) for t in TA + TN for m in ALLM}
best_loto = {(t, m): max(LS, key=lambda L: np.mean([point[(u, m, L)] for u in TA if u != t])) for t in TA for m in ALLM}

## re-selected lengths inside every resample (index into LS, B per task and method)
sel_task = {(t, m): bmat(t, m).argmax(axis=1) for t in TA + TN for m in ALLM}
sel_loto = {(t, m): np.mean([bmat(u, m) for u in TA if u != t], axis=0).argmax(axis=1) for t in TA for m in ALLM}


def pick(t, m, sel):
    """resampled AUC of m on t, at the length index sel[b] in resample b"""
    return bmat(t, m)[np.arange(B), sel]


rows = []
VERSIONS = [('per task', TA + TN), ('leave one task out', TA)]
for version, tasks in VERSIONS:
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
                         'RFind-sc AUC': point[(t, 'RFind-sc', bl['RFind-sc'])], 'method length': bl[m], 'method AUC': point[(t, m, bl[m])],
                         'RFind-sc minus method': diff, 'difference 95% CI low': lo, 'difference 95% CI high': hi,
                         'result': verdict(diff, lo, hi), 'comparable': comparable(lo, hi),
                         'task set': TSET[t], 'method set': 'main' if m in METH else 'table only',
                         'reselected difference 95% CI low': slo, 'reselected difference 95% CI high': shi,
                         'reselected result': verdict(diff, slo, shi), 'reselected comparable': comparable(slo, shi),
                         'RFind-sc share of resamples at its point-estimate length': share})
out = pd.DataFrame(rows)
out.to_csv(f'{O}results/fig5_best_length.csv', index=False)

## mean AUC over the 4 main tasks: point, resampled (B), per comparison and method
COMPS = [(f'fixed length {L}', L) for L in LS] + [('best length per task', 'best'), ('leave one task out', 'loto')]


def mean_point_boot(m, kind):
    if kind in LS:
        p = np.mean([point[(t, m, kind)] for t in TA]); bb = np.mean([boot[t][:, col[t][(m, kind)]] for t in TA], axis=0)
    elif kind == 'best':
        p = np.mean([point[(t, m, best_task[(t, m)])] for t in TA])
        bb = np.mean([boot[t][:, col[t][(m, best_task[(t, m)])]] for t in TA], axis=0)
    else:
        p = np.mean([point[(t, m, best_loto[(t, m)])] for t in TA]); bb = np.mean([boot[t][:, col[t][(m, best_loto[(t, m)])]] for t in TA], axis=0)
    return p, bb


mrows = []
for cname, kind in COMPS:
    pr, br = mean_point_boot('RFind-sc', kind)
    for m in ALLM[1:]:
        p, bb = mean_point_boot(m, kind)
        diff = pr - p
        lo, hi = np.quantile(br - bb, [0.025, 0.975])
        mrows.append({'comparison': cname, 'method': m, 'method set': 'main' if m in METH else 'table only', 'mean AUC RFind-sc': pr,
                      'mean AUC method': p, 'difference': diff, 'CI low': lo, 'CI high': hi, 'result': verdict(diff, lo, hi),
                      'comparable': comparable(lo, hi)})
md = pd.DataFrame(mrows); md.to_csv(f'{O}results/fig5_mean_auc_diff.csv', index=False)

## summary
print(f'point AUCs vs results/fig5_redo_auc_L<L>.csv: max |diff| {maxdiff:.1e}')


def counts(g, rc='result'):
    c = g[rc].value_counts()
    return (f'higher {c.get("RFind-sc higher", 0)}, not distinguished {c.get("not distinguished", 0)}, lower {c.get("RFind-sc lower", 0)} '
            f'of {len(g)}; comparable {int(g[rc.replace("result", "comparable")].sum())}')


main = out[out['task id'].isin(TA) & out.method.isin(METH)]
for v in ['per task', 'leave one task out']:
    g = main[main.version == v]
    for lab, gg in [('main 9 methods', g), ('8 non-decoupleR', g[g.method != 'ULM (log2FC)']), ('ULM (log2FC)', g[g.method == 'ULM (log2FC)'])]:
        print(f'{v} | {lab}: {counts(gg)} || re-selected: {counts(gg, "reselected result")}')
    print(f'   RFind-sc lengths: {dict(g.groupby("task id", sort=False)["RFind-sc length"].first())}; share at point length: '
          f'{dict(g.groupby("task id", sort=False)["RFind-sc share of resamples at its point-estimate length"].first())}')
    print('   lower:', [(r['task id'], r.method, r['method length'], round(r['RFind-sc minus method'], 3)) for _, r in g[g.result == 'RFind-sc lower'].iterrows()])
print('leave one task out lengths (task: RFind-sc):', {t: best_loto[(t, 'RFind-sc')] for t in TA})
def mean_aucs(cname):
    g = md[md.comparison == cname]
    return {'RFind-sc': round(float(g['mean AUC RFind-sc'].iloc[0]), 4), **{m: round(float(v), 4) for m, v in zip(g.method, g['mean AUC method'])}}


print('mean AUC over the 4 main tasks, leave one task out:', mean_aucs('leave one task out'))
print('mean AUC over the 4 main tasks, best length per task:', mean_aucs('best length per task'))
for m in ['ULM (log2FC)', 'WMEAN (log2FC)']:
    print(f'mean AUC difference, RFind-sc minus {m}:')
    for _, r in md[md.method == m].iterrows():
        print(f'   {r.comparison}: {r["mean AUC RFind-sc"]:.4f} vs {r["mean AUC method"]:.4f}, diff {r.difference:+.4f} '
              f'[{r["CI low"]:+.4f}, {r["CI high"]:+.4f}] {r.result}{" (comparable)" if r.comparable else ""}')
print('table-only methods, per task (main tasks):')
for _, r in out[(out.version == 'per task') & out['task id'].isin(TA) & out.method.isin(TONLY)].iterrows():
    print(f'   {r["task id"]} {r.method} (L {r["method length"]}): diff {r["RFind-sc minus method"]:+.4f} [{r["difference 95% CI low"]:+.4f}, '
          f'{r["difference 95% CI high"]:+.4f}] {r.result} | re-selected [{r["reselected difference 95% CI low"]:+.4f}, '
          f'{r["reselected difference 95% CI high"]:+.4f}] {r["reselected result"]}')
