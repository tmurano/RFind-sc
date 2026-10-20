"""Fig. 5e / Supplementary Table 3 — the mean AUC over the tasks of each type, RFind-sc minus each method, with 95% CIs:
the eight type A tasks with references of 500 genes per direction (the DAM and Glucocorticoid tasks H1, D1, D2, D3 and the PS2APP and Forced
swim tasks DA1, AS1, AS2, AS3) and the three type B tasks (H2, H3, H5 of Fig. 5d).
  - per-cell scores of the benchmark; Beyondcell without a score = 0; AUC = Mann-Whitney with ties averaged
  - bootstrap over mice (fig5_mice.py: all cells of a drawn mouse): default_rng(1), one generator, 1,000 resamples, the same resamples
    for every method; the tasks ask for their draws in the order H1, D1, D2, D3, DA1, AS1, AS2, AS3, H2, H3, H5, and the tasks with the
    same mice share the draws (Hammond: H1, DA1; De Donno: D1-D3, AS1-AS3; H2, H5); the type A draws are therefore those of
    fig5_best_length_8tasks.py, whose mean AUC differences at 500 genes per direction are reproduced here (the script stops otherwise)
  - mean over tasks: point = mean of the point AUCs; resample b = mean of the tasks' resampled AUCs; CI = 2.5% and 97.5% quantiles;
    RFind-sc higher / lower when the CI of the difference excludes 0 and |difference| > 0.01, otherwise not distinguished
  - methods: the methods of Fig. 5c (type A) and 5d (type B), and the table-only methods UCell (up genes), WMEAN and ORA
Inputs: figs/Fig5_benchmark/output/results/fig5_redo_percell_L500.csv.gz (H1, D1-D3, H2, H3, H5),
        extra_refs/fig5_extra_percell_L500.csv.gz (DA1, AS1-AS3), results/fig5_mean_auc_diff_8tasks.csv (the check)
Output: figs/Fig5_benchmark/output/results/fig5_task_means.csv  one row per type and method: tasks, mean AUC with its 95% CI, RFind-sc
        minus the method with its 95% CI, result
Run from the project root (about 6 min): python figs/Fig5_benchmark/analysis/fig5_task_means.py
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig5_mice import MouseDraws, mouse_of  # noqa: E402

O = 'figs/Fig5_benchmark/output/'
B, SEED = 1000, 1
TYPE_A = ['H1', 'D1', 'D2', 'D3', 'DA1', 'AS1', 'AS2', 'AS3']      # the order of fig5_best_length_8tasks.py
TYPE_B = ['H2', 'H3', 'H5']
METH_A = ['RFind-sc', 'AddModuleScore', 'AUCell', 'UCell', 'singscore', 'ssGSEA', 'VISION', 'Beyondcell', 'ULM (log2FC)', 'scFOCAL-type',
          'UCell (up genes)', 'WMEAN (log2FC)', 'ORA']               # RFind-sc first: column 0 is the reference of the differences
METH_B = ['RFind-sc', 'AddModuleScore', 'AUCell', 'UCell', 'singscore', 'ssGSEA', 'VISION', 'Beyondcell', 'ULM', 'WMEAN', 'ORA']


def aucs(S, pos):
    r = rankdata(S, axis=0); n1 = pos.sum(); n0 = len(pos) - n1
    return (r[pos].sum(axis=0) - n1 * (n1 + 1) / 2) / (n1 * n0)


def verdict(diff, lo, hi):
    return 'RFind-sc higher' if (lo > 0 and diff > 0.01) else ('RFind-sc lower' if (hi < 0 and diff < -0.01) else 'not distinguished')


pc = pd.concat([pd.read_csv(O + 'results/fig5_redo_percell_L500.csv.gz', low_memory=False),
                pd.read_csv(O + 'extra_refs/fig5_extra_percell_L500.csv.gz', low_memory=False)], ignore_index=True, sort=False)
assert pc.positive.dtype == bool and set(TYPE_A + TYPE_B) <= set(pc.id), sorted(set(pc.id))
DRAWS = MouseDraws(np.random.default_rng(SEED), B)
rows = []
for kind, tasks, ms in (('A', TYPE_A, METH_A), ('B', TYPE_B, METH_B)):
    point, boot = [], []
    for t in tasks:
        g = pc[pc.id == t]
        pos = g.positive.to_numpy(dtype=bool)
        S = g[ms].to_numpy(float)
        S[:, ms.index('Beyondcell')] = np.nan_to_num(S[:, ms.index('Beyondcell')], nan=0.0)   # cells without a score = 0
        assert not np.isnan(S).any(), (t, [m for m, k in zip(ms, np.isnan(S).any(0)) if k])
        point.append(aucs(S, pos))
        boot.append(np.array([aucs(S[i], pos[i]) for i in DRAWS.cells(mouse_of(g.cell_id), pos)]))
    p, bm = np.mean(point, axis=0), np.mean(boot, axis=0)              # methods; B x methods (mean over tasks in every resample)
    for k, m in enumerate(ms):
        lo, hi = np.percentile(bm[:, k], [2.5, 97.5])
        d_lo, d_hi = np.percentile(bm[:, 0] - bm[:, k], [2.5, 97.5])
        diff = p[0] - p[k]
        rows.append({'type': kind, 'tasks': ', '.join(tasks), 'method': m, 'mean AUC': p[k], 'mean AUC 95% CI low': lo,
                     'mean AUC 95% CI high': hi, 'RFind-sc minus method': np.nan if k == 0 else diff,
                     'difference 95% CI low': np.nan if k == 0 else d_lo, 'difference 95% CI high': np.nan if k == 0 else d_hi,
                     'result': '' if k == 0 else verdict(diff, d_lo, d_hi)})
res = pd.DataFrame(rows)

## check: the type A differences are those of fig5_best_length_8tasks.py at 500 genes per direction (the same draws)
ref = pd.read_csv(O + 'results/fig5_mean_auc_diff_8tasks.csv')
ref = ref[ref.comparison == 'fixed length 500'].set_index('method')
a = res[(res.type == 'A') & (res.method != 'RFind-sc')].set_index('method')
d = max(np.abs(a['RFind-sc minus method'] - ref.loc[a.index, 'difference']).max(),
        np.abs(a['difference 95% CI low'] - ref.loc[a.index, 'CI low']).max(), np.abs(a['difference 95% CI high'] - ref.loc[a.index, 'CI high']).max())
assert d < 1e-9, f'type A differences differ from results/fig5_mean_auc_diff_8tasks.csv (max |diff| {d:.2g})'
print(f'check, type A vs results/fig5_mean_auc_diff_8tasks.csv (500 genes per direction): max |diff| {d:.1e}')

os.makedirs(O + 'results', exist_ok=True)
out = O + 'results/fig5_task_means.csv'
res.to_csv(out, index=False)
main = METH_A[:10] + ['ULM']
for kind in ('A', 'B'):
    x = res[(res.type == kind) & res.method.isin(main)]
    print(f'\ntype {kind} ({x.tasks.iloc[0]}):')
    print(x[['method', 'mean AUC', 'RFind-sc minus method', 'difference 95% CI low', 'difference 95% CI high', 'result']].round(3).to_string(index=False))
print('wrote', out, f'({len(res)} rows)')
